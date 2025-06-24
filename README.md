# W_NucleiSegmentation-Cellpose

2D/3D nuclei segmentation using Cellpose with pretrained models. This container processes multi-dimensional images (TZCYX) and handles channel selection, time points, and z-slices.

## References
- Cellpose-SAM: superhuman generalization for cellular segmentation  
  Marius Pachitariu, Michael Rariden, Carsen Stringer  
  bioRxiv 2025.04.28.651001; https://doi.org/10.1101/2025.04.28.651001

## Build Container
```bash
docker build -t cellpose_biomero .
```

## Key Parameters
- `diameter`: Nuclei diameter in pixels (0 for auto-estimation)
- `prob_threshold`: Probability threshold (default: 0.5)
- `cp_model`: Model type ('nuclei' or 'cyto')
- `nuc_channel`: Channel to segment (0-based, -1 for all)
- `use_gpu`: Enable GPU acceleration
- `auto_tiling`: Automatic tiling for large images
- `bsize`: Tile size for processing (default: 224)

## Test Locally
Create test directories and run:
```bash
mkdir -p E:\tmp\cellpose\{in,out,gt}
docker run --rm --gpus=all -v E:\tmp\cellpose:/data cellpose_biomero --local --infolder /data/in --outfolder /data/out --gtfolder /data/gt --diameter 30
```

## Saving Cellpose models

By default Cellpose is downloading models to ```/tmp/models/cellpose/```, by mounting this folder to the container it is possible to store the model files externally and reload the models at the next run of the container.   
For biomero it is possible to define in the 'slurm_data_bind_path' a custom binding to the singularity container, by adding e.g. ```/data1/models:/tmp/models``` models will be saved at ```/data1/models/cellpose``` on the HPC.
