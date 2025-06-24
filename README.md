# W_Segmentation-Cellpose4

2D/3D cell segmentation using Cellpose 4 with Cellpose-SAM pretrained models. This container processes multi-dimensional images (TZCYX) and handles channel selection, time points, and z-slices.

## References
- Cellpose-SAM: superhuman generalization for cellular segmentation  
  Marius Pachitariu, Michael Rariden, Carsen Stringer  
  bioRxiv 2025.04.28.651001; https://doi.org/10.1101/2025.04.28.651001

## Build Container
```bash
docker build -t cellpose_biomero .
```

## Key Parameters
- `diameter`: Cell diameter in pixels (0 for auto-estimation)
- `prob_threshold`: Probability threshold (default: 0.5)
- `cp_model`: Model type (default: 'cpsam', options: 'cyto', 'nuclei', or custom path)
- `nuc_channel`: Channel to segment (0-based, -1 for all)
- `use_gpu`: Enable GPU acceleration
- `auto_tiling`: Automatic tiling for large images (default: false)
- `bsize`: Tile size for processing (default: 512)

## Test Locally
Create test directories and run:
```bash
mkdir -p E:\tmp\cellpose\{in,out,gt}
docker run --rm --gpus=all -v E:\tmp\cellpose:/data cellpose_biomero --local --infolder /data/in --outfolder /data/out --gtfolder /data/gt --diameter 30
```

## Saving Cellpose models

By default Cellpose is downloading models to ```/tmp/models/cellpose/```, by mounting this folder to the container it is possible to store the model files externally and reload the models at the next run of the container.   
For biomero it is possible to define in the 'slurm_data_bind_path' a custom binding to the singularity container, by adding e.g. ```/data1/models:/tmp/models``` models will be saved at ```/data1/models/cellpose``` on the HPC.

## Complete Parameter Documentation

All parameters below are actively used in the segmentation workflow. No obsolete parameters exist in the current version.

### Core Segmentation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `diameter` | Number | 30 | Cell diameter in pixels (0 for auto-estimation) |
| `prob_threshold` | Number | 0.5 | Probability threshold, centered at 0.0 |
| `cp_model` | String | cpsam | Cellpose model ('cpsam', 'cyto', 'nuclei', or custom path) |
| `nuc_channel` | Number | 0 | Channel to segment (0-based, -1 for all channels) |
| `use_gpu` | Boolean | true | Enable GPU acceleration |

### Image Processing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `time_series` | Number | -1 | Process specific time point (0-based, -1 for all) |
| `z_slices` | Number | -1 | Process specific z-slice (0-based, -1 for all) |
| `flow_threshold` | Number | 0.4 | Flow error threshold (0 disables QC step) |
| `min_size` | Number | 15 | Minimum pixels per mask (-1 to disable) |
| `exclude_on_edges` | Boolean | false | Discard masks touching image edges |

### Performance Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `batch_size` | Number | 8 | Number of images to process in parallel |
| `bsize` | Number | 512 | Tile size when tiling is enabled |
| `auto_tiling` | Boolean | false | Auto-enable tiling for large images (>2048x2048) |

### 3D Processing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `do_3D` | Boolean | false | Process images as 3D stacks |
| `anisotropy` | Number | 1.0 | Z-axis anisotropy factor for 3D processing |
