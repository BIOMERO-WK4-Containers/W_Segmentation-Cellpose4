# -*- coding: utf-8 -*-

# * Copyright (c) 2009-2018. Authors: see NOTICE file.
# *
# * Licensed under the Apache License, Version 2.0 (the "License");
# * you may not use this file except in compliance with the License.
# * You may obtain a copy of the License at
# *
# *      http://www.apache.org/licenses/LICENSE-2.0
# *
# * Unless required by applicable law or agreed to in writing, software
# * distributed under the License is distributed on an "AS IS" BASIS,
# * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# * See the License for the specific language governing permissions and
# * limitations under the License.

import sys
import os
import shutil
import subprocess
import time
import numpy as np
from tifffile import imwrite, TiffFile
import skimage
import skimage.color
from cytomine.models import Job
from biaflows import CLASS_OBJSEG
from biaflows.helpers import BiaflowsJob, prepare_data, upload_data, upload_metrics, get_discipline

def convert_to_5d_from_tifffile(volume, axes, target="XYZCT"):
    """
    Convert a numpy array from TiffFile to 5D dimensions suitable for OMERO
    
    Parameters
    ----------
    volume : numpy.ndarray
        Image data from tifffile's asarray()
    axes : str
        Axes string from tifffile (e.g., 'TZCYX', 'YX', etc.)
    target : str, optional
        String specifying the desired dimension order, default is "XYZCT"
        
    Returns
    -------
    img_5d : numpy.ndarray or tuple
        5D numpy array with dimensions ordered according to target
        When unpacked as a tuple, returns (img_5d, target)
    """
    # Validate input volume is a numpy array
    if not isinstance(volume, np.ndarray):
        raise TypeError("Input volume must be a numpy.ndarray")
    
    # Standardize to uppercase
    axes = axes.upper()
    target = target.upper()
    
    # Validate axes dimensions match array dimensions
    if len(axes) != volume.ndim:
        raise ValueError(f"Axes string '{axes}' does not match array dimensions {volume.ndim}")
    
    # Some TIFF files use 'S' for samples/channels, convert to 'C' for consistency
    axes = axes.replace('S', 'C')
        
    # Validate target dimensions
    if len(target) != 5:
        raise ValueError(f"Target dimensions must have exactly 5 dimensions, got '{target}'")
    
    if set(target) != set("XYZCT"):
        raise ValueError("Target dimensions must contain letters X, Y, Z, C, and T exactly once")
    
    # Create a 5D array by adding missing dimensions
    img_5d = volume
    current_order = axes
    
    # Add missing dimensions
    for dim in "XYZCT":
        if dim not in current_order:
            img_5d = np.expand_dims(img_5d, axis=-1)
            current_order += dim
    
    # Reorder dimensions if needed
    if current_order != target:
        # Create list of current positions for each dimension
        current_positions = []
        for dim in target:
            current_positions.append(current_order.index(dim))
        
        # Rearrange dimensions
        img_5d = np.moveaxis(img_5d, current_positions, range(len(target)))
    
    # Return both the array and target, allowing for flexible unpacking
    class ReturnValue(tuple):
        """Custom return class to allow both direct access and unpacking"""
        def __new__(cls, img, axes):
            return tuple.__new__(cls, (img, axes))
            
        def __repr__(self):
            return repr(self[0])
            
        # Make the first element (the image) accessible directly
        def __array__(self, dtype=None):
            return np.asarray(self[0], dtype=dtype)
    
    return ReturnValue(img_5d, target)

def process_image(img_path, out_dir, bj, params):
    """
    Process a 5D image with Cellpose, handling individual dimensions appropriately
    
    Args:
        img_path: Path to the input image
        out_dir: Directory to save processed slices 
        bj: BiaFlows job object
        params: Parameter object with processing settings
        
    Returns:
        tuple: (output_path, dimension_mapping)
            - output_path: Path to the output image
            - dimension_mapping: Dictionary of processed dimensions
    """
    # Create unique subdirectory for this image to avoid conflicts
    img_name = os.path.basename(img_path)
    img_tmp_dir = os.path.join(out_dir, f"tmp_{os.path.splitext(img_name)[0]}_{int(time.time())}")
    os.makedirs(img_tmp_dir, exist_ok=True)
    
    with TiffFile(img_path) as tif:
        # Load image data and get axes information
        volume = tif.asarray()
        original_axes = tif.series[0].axes
        bj.job.update(status=Job.RUNNING, progress=30,
                     statusComment=f"Processing image with original axes: {original_axes}")
        
        # Convert to 5D with TZCYX order for processing
        bj.job.update(status=Job.RUNNING, progress=35,
                     statusComment="Converting to 5D TZCYX format...")
        img_5d, axes_5d = convert_to_5d_from_tifffile(volume, original_axes, target="TZCYX")
        
        # Track dimension sizes for clarity
        dims = {
            'T': img_5d.shape[0],  # Time
            'Z': img_5d.shape[1],  # Depth/slices
            'C': img_5d.shape[2],  # Channels
            'Y': img_5d.shape[3],  # Height
            'X': img_5d.shape[4]   # Width
        }
        bj.job.update(status=Job.RUNNING, progress=40,
                     statusComment=f"Dimensions (TZCYX): {dims}")
        
        # Get channel, time, and z-slice parameters with defensive handling
        # Get the nuclear channel with fallback to 0
        nuc_channel = getattr(params, 'nuc_channel', 0)
        
        # Get time point parameter with fallback to -1 (all frames)
        time_point = -1  # Default to all frames
        for param_name in ['time_series', 'time_point', 't']:
            if hasattr(params, param_name):
                time_point = getattr(params, param_name)
                break
        
        # Get z-slice parameter with fallback to -1 (all slices)
        z_slice = -1  # Default to all slices
        for param_name in ['z_slices', 'z_slice', 'z']:
            if hasattr(params, param_name):
                z_slice = getattr(params, param_name)
                break
                
        bj.job.update(status=Job.RUNNING, progress=45,
                      statusComment=f"Using parameters: channel={nuc_channel}, time={time_point}, z={z_slice}")
        
        # Prepare channel list to process
        if isinstance(nuc_channel, (int, np.integer)):
            if nuc_channel == -1:
                # Process all channels
                channels_to_process = list(range(dims['C']))
            else:
                channels_to_process = [nuc_channel]
        elif isinstance(nuc_channel, (list, tuple, np.ndarray)):
            channels_to_process = list(nuc_channel)
        else:
            # Default to first channel if type is unexpected
            channels_to_process = [0]
        
        # Prepare time points to process
        if isinstance(time_point, (int, np.integer)):
            if time_point == -1:
                # Process all time points
                time_points_to_process = list(range(dims['T']))
            else:
                time_points_to_process = [time_point]
        elif isinstance(time_point, (list, tuple, np.ndarray)):
            time_points_to_process = list(time_point)
        else:
            # Default to all time points if type is unexpected
            time_points_to_process = list(range(dims['T']))
        
        # Prepare z-slices to process
        if isinstance(z_slice, (int, np.integer)):
            if z_slice == -1:
                # Process all z-slices
                z_slices_to_process = list(range(dims['Z']))
            else:
                z_slices_to_process = [z_slice]
        elif isinstance(z_slice, (list, tuple, np.ndarray)):
            z_slices_to_process = list(z_slice)
        else:
            # Default to all z-slices if type is unexpected
            z_slices_to_process = list(range(dims['Z']))
            
        # Validate indices
        for channel in channels_to_process:
            if channel >= dims['C']:
                raise ValueError(f"Invalid channel index {channel}. Image has {dims['C']} channels (0-{dims['C']-1})")
        
        for t in time_points_to_process:
            if t >= dims['T']:
                raise ValueError(f"Invalid time point {t}. Image has {dims['T']} time points (0-{dims['T']-1})")
        
        for z in z_slices_to_process:
            if z >= dims['Z']:
                raise ValueError(f"Invalid z-slice {z}. Image has {dims['Z']} z-slices (0-{dims['Z']-1})")
        
        # Create output array - shape matches the selected dimensions
        output_shape = (
            len(time_points_to_process),
            len(z_slices_to_process),
            len(channels_to_process),
            dims['Y'],
            dims['X']
        )
        output = np.zeros(output_shape, dtype=np.uint16)
        
        # Calculate total slices for progress tracking
        total_slices = len(channels_to_process) * len(z_slices_to_process) * len(time_points_to_process)
        current_slice = 0
        
        # Process each dimension
        processed_files = []
        slice_metadata = []
        
        for ch_idx, channel in enumerate(channels_to_process):
            bj.job.update(status=Job.RUNNING,
                         statusComment=f"Processing channel {channel} ({ch_idx+1}/{len(channels_to_process)})")
            
            for t_idx, t in enumerate(time_points_to_process):
                for z_idx, z in enumerate(z_slices_to_process):
                    # Extract YX slice from current TZC position
                    img_slice = img_5d[t, z, channel, :, :]  # Note TZCYX order
                    
                    # Create a unique filename for this slice
                    slice_filename = f"t{t:03d}_z{z:03d}_c{channel:03d}.tif"
                    slice_path = os.path.join(img_tmp_dir, slice_filename)
                    
                    # Set minimum dimension handling - pad if needed
                    minshape = min(img_slice.shape[:2])
                    maxshape = max(img_slice.shape[:2])
                    original_shape = img_slice.shape
                    
                    # Ensure dimensions are at least 224x224 (Cellpose minimum)
                    if minshape != maxshape or minshape < 224:
                        # Get image dimensions
                        height, width = img_slice.shape[:2]  # Get first two dimensions
                        
                        # Calculate padding needed for height and width
                        pad_height = max(0, max(224, maxshape) - height)
                        pad_width = max(0, max(224, maxshape) - width)
                        
                        # Create appropriate padding configuration based on image dimensions
                        if len(img_slice.shape) == 2:  # Grayscale
                            padshape = ((0, pad_height), (0, pad_width))
                        elif len(img_slice.shape) == 3:  # With channels
                            padshape = ((0, pad_height), (0, pad_width), (0, 0))
                        else:
                            # Add padding for each dimension
                            padshape = [(0, 0)] * len(img_slice.shape)
                            padshape[0] = (0, pad_height)
                            padshape[1] = (0, pad_width)
                        
                        # Apply padding
                        img_slice = np.pad(img_slice, padshape, 'constant', constant_values=0)
                    
                    # Save the slice
                    imwrite(slice_path, img_slice)
                    processed_files.append(slice_path)
                    
                    # Store metadata for reconstruction
                    slice_metadata.append({
                        'path': slice_path,
                        'filename': slice_filename,
                        't': t,
                        'z': z,
                        'c': channel,
                        't_idx': t_idx,
                        'z_idx': z_idx,
                        'c_idx': ch_idx,
                        'original_shape': original_shape,
                        'padded': original_shape != img_slice.shape,
                        'padded_shape': img_slice.shape if original_shape != img_slice.shape else None
                    })
                    
                    # Update progress
                    current_slice += 1
                    progress = 40 + (20 * current_slice / total_slices)
                    bj.job.update(status=Job.RUNNING, progress=int(progress),
                                statusComment=f"Prepared slice {current_slice}/{total_slices}: t{t}, z{z}, c{channel}")
        
        # Run Cellpose on all prepared slices
        bj.job.update(status=Job.RUNNING, progress=60,
                    statusComment="Running Cellpose on prepared slices...")
        
        # Cellpose parameters with defensive handling
        # Get model type with fallback
        model_type = getattr(params, 'cp_model', 'cyto')
        if not model_type:
            model_type = 'cyto'  # Default to cyto model
            
        # Get diameter with fallback
        diameter = getattr(params, 'diameter', 30)
        
        # Get probability threshold with fallbacks
        cellprob_threshold = 0.0  # Default
        for param_name in ['cellprob_threshold', 'prob_threshold', 'cell_probability']:
            if hasattr(params, param_name):
                cellprob_threshold = getattr(params, param_name)
                break
                
        # Get flow threshold with fallback
        flow_threshold = getattr(params, 'flow_threshold', 0.4)
        
        # Get minimum size with fallback
        min_size = getattr(params, 'min_size', 15)
        
        # Get GPU flag with fallback
        use_gpu = getattr(params, 'use_gpu', False)

        batch_size = getattr(params, 'batch_size', 8)
        
        # Building Cellpose command line
        cmd = [
            "conda", "run", "-n", "cellpose", "python", "-m", "cellpose", 
            "--dir", img_tmp_dir,
            "--pretrained_model", model_type,
            "--save_tif",
            "--no_npy",
            "--diameter", f"{diameter}",
            "--cellprob_threshold", f"{cellprob_threshold}",
            "--flow_threshold", f"{flow_threshold}",
            "--min_size", f"{min_size}",
            "--batch_size", f"{batch_size}",
            '--verbose'
        ]
        
        # Add optional parameters
        if use_gpu:
            cmd.append("--use_gpu")
        
        # Add tiling parameters with defensive handling
        #cmd.append("--tile")  # Enable tiling
        
        # Get tile overlap with fallback
        tile_overlap = getattr(params, 'tile_overlap', 0.1)
        #cmd.extend(["--tile_overlap", f"{tile_overlap}"])
        
        # Get block size with fallback
        bsize = getattr(params, 'bsize', 224)
        #cmd.extend(["--bsize", f"{bsize}"])
        
        # Run Cellpose
        bj.job.update(status=Job.RUNNING, progress=65,
                    statusComment=f"Running Cellpose with command: {' '.join(cmd)}")
        
        status = subprocess.run(cmd)
        
        if status.returncode != 0:
            error_msg = f"Cellpose failed with exit code {status.returncode}"
            bj.job.update(status=Job.FAILED, progress=65, statusComment=error_msg)
            raise RuntimeError(error_msg)
        
        # Process results
        bj.job.update(status=Job.RUNNING, progress=80,
                    statusComment="Processing Cellpose results...")
        
        # Reconstruct 5D array from processed slices
        for metadata in slice_metadata:
            t_idx = metadata['t_idx']
            z_idx = metadata['z_idx']
            c_idx = metadata['c_idx']
            
            # Get the mask file (Cellpose adds _cp_masks.tif to the filename)
            base_filename = os.path.splitext(metadata['filename'])[0]
            mask_filename = f"{base_filename}_cp_masks.tif"
            mask_path = os.path.join(img_tmp_dir, mask_filename)
            
            if os.path.exists(mask_path):
                # Load mask
                mask = skimage.io.imread(mask_path)
                
                # Crop back to original dimensions if padded
                if metadata['padded']:
                    orig_h, orig_w = metadata['original_shape']
                    mask = mask[:orig_h, :orig_w]
                
                # Store in output array
                output[t_idx, z_idx, c_idx, :mask.shape[0], :mask.shape[1]] = mask
            else:
                bj.job.update(status=Job.RUNNING, 
                             statusComment=f"Warning: Mask not found for {metadata['filename']}")
        
        # Create final output path
        output_name = f"{os.path.splitext(img_name)[0]}_cellpose.tif"
        output_path = os.path.join(out_dir, output_name)
        
        # Add dimension mapping information to output metadata
        dimension_mapping = {
            'T': time_points_to_process,
            'Z': z_slices_to_process,
            'C': channels_to_process
        }
        
        # Save final 5D result
        metadata = {
            'axes': 'TZCYX',
            'dimension_mapping': str(dimension_mapping),
            'cellpose_params': {
                'model': model_type,
                'diameter': diameter,
                'cellprob_threshold': cellprob_threshold,
                'flow_threshold': flow_threshold,
                'min_size': min_size                
            }
        }
        
        imwrite(output_path, output,
               metadata=metadata,
               photometric='minisblack',
               ome=True,
               description='Processed with Cellpose, standardized to TZCYX format')
        
        # Clean up temporary directory
        try:
            shutil.rmtree(img_tmp_dir)
        except Exception as e:
            bj.job.update(status=Job.RUNNING, statusComment=f"Warning: Could not remove temp directory: {str(e)}")
        
        return output_path, dimension_mapping

def main(argv):
    base_path = "{}".format(os.getenv("HOME"))  # Mandatory for Singularity
    problem_cls = CLASS_OBJSEG

    with BiaflowsJob.from_cli(argv) as bj:
        bj.job.update(status=Job.RUNNING, progress=0, statusComment="Initialization...")

        # 1. Prepare data for workflow
        in_imgs, gt_imgs, in_path, gt_path, out_path, tmp_path = prepare_data(problem_cls, bj, is_2d=False, **bj.flags)
        
        # MAKE SURE TMP PATH IS UNIQUE
        try:
            tmp_path = os.path.join(tmp_path, f"{int(time.time() * 1000)}")  # timestamp in ms
            os.makedirs(tmp_path, exist_ok=True)  # setup tmp
        except Exception as e:
            # If error, try with a different timestamp
            tmp_path = os.path.join(tmp_path, f"{int(time.time() * 10000)}")  # timestamp in ms
            os.makedirs(tmp_path, exist_ok=True)  # setup tmp
        
        # Get list of image files
        list_imgs = [image.filepath for image in in_imgs]
        
        # Process each image
        bj.job.update(progress=10, statusComment=f"Number of images to process: {len(list_imgs)}")
        
        for img_index, img_path in enumerate(list_imgs):
            img_name = os.path.basename(img_path)
            bj.job.update(progress=int(10 + (70 * img_index / len(list_imgs))),
                        statusComment=f"Processing image {img_index+1}/{len(list_imgs)}: {img_name}")
            
            try:
                # Process image with Cellpose
                output_path, dimension_mapping = process_image(img_path, out_path, bj, bj.parameters)
                
                # Log information about the processed image
                bj.job.update(progress=int(80 + (10 * (img_index + 1) / len(list_imgs))),
                            statusComment=f"Completed {img_name}: Saved to {output_path}")
            
            except Exception as e:
                error_msg = f"Error processing {img_name}: {str(e)}"
                bj.job.update(status=Job.RUNNING, progress=80, statusComment=error_msg)
                import traceback
                traceback.print_exc()
                # Continue with next image
                continue

        # 3. Upload data to BIAFLOWS
        bj.job.update(progress=90, statusComment="Uploading results...")
        upload_data(problem_cls, bj, in_imgs, out_path, **bj.flags, monitor_params={
            "start": 90, "end": 95, "period": 0.1,
            "prefix": "Extracting and uploading polygons from masks"})
        
        # 4. Compute metrics
        bj.job.update(progress=95, statusComment="Computing and uploading metrics...")
        upload_metrics(problem_cls, bj, in_imgs, gt_path, out_path, tmp_path, **bj.flags)
        
        # 5. Clean up
        try:
            shutil.rmtree(tmp_path)  # cleanup tmp
        except Exception as e:
            bj.job.update(status=Job.RUNNING, 
                         statusComment=f"Warning: Could not remove temp directory: {str(e)}")
        
        # 6. Finish
        bj.job.update(progress=100, status=Job.TERMINATED, status_comment="Finished.")

if __name__ == "__main__":
    main(sys.argv[1:])