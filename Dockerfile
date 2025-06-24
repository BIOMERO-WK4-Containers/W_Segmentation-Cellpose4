# Choose a base image with CUDA support.
# Select a CUDA version compatible with the PyTorch needed for cellpose and your host drivers.
# Using an Ubuntu base often works well with Conda. Using a -devel image includes compilers if needed.
# Example: CUDA 11.8 on Ubuntu 20.04
FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu20.04
# Example: CUDA 12.1 on Ubuntu 22.04
# FROM nvidia/cuda:12.1.1-cudnn8-devel-ubuntu22.04

# Set environment variables to prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH

# >>> ADDED: Define default directories for models <<<
ENV CELLPOSE_LOCAL_MODELS_PATH=/tmp/models/cellpose/

# Install base dependencies and cleanup
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        git \
        bzip2 \
        ca-certificates \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        libgeos-dev \
        libgl1-mesa-dev \
        build-essential \
        libxcb-cursor0 \
        libxcb-xinerama0 \
        libxcb-icccm4 \
        libxcb-image0 \
        libxcb-keysyms1 \
        libxcb-randr0 \
        libxcb-render-util0 \
        libxcb-xkb1 \
        libxkbcommon-x11-0 \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*
# >>> ADDED: Create the default cache directories <<<
# Ensure the default cache directories exist within the container
RUN mkdir -p ${CELLPOSE_LOCAL_MODELS_PATH} && chmod 777 ${CELLPOSE_LOCAL_MODELS_PATH}

# ------------------------------------------------------------------------------
# Install Miniconda
# ------------------------------------------------------------------------------
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -a -y && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    conda init bash # Initialize conda for the SHELL

# Make conda available in RUN instructions using the initialized shell
SHELL ["/bin/bash", "--login", "-c"]

# Update conda
RUN conda update -n base -c defaults conda --yes

# ------------------------------------------------------------------------------
# Create Cytomine environment (Python 3.7)
# ------------------------------------------------------------------------------
ENV CYTOMINE_ENV_NAME=cytomine_py37
RUN conda create -n $CYTOMINE_ENV_NAME python=3.7 -y

RUN conda run -n $CYTOMINE_ENV_NAME pip install --no-cache-dir \
        git+https://github.com/cytomine-uliege/Cytomine-python-client.git@v2.7.3

RUN conda run -n $CYTOMINE_ENV_NAME pip install --no-cache-dir \
        git+https://github.com/Neubias-WG5/biaflows-utilities.git@v0.9.2


# ------------------------------------------------------------------------------
# Create a dedicated Cellpose environment with CUDA support
# ------------------------------------------------------------------------------
ENV CELLPOSE_ENV_NAME=cellpose_env

# Create environment with Python 3.10 for better compatibility with PyTorch/CUDA
RUN conda create -n $CELLPOSE_ENV_NAME -c conda-forge python=3.10 -y

# Install PyTorch with CUDA 11.8 support
RUN conda run -n $CELLPOSE_ENV_NAME conda install -c pytorch -c nvidia \
    pytorch torchvision torchaudio pytorch-cuda==11.8 -y

# Install Cellpose with GPU support
RUN conda run -n $CELLPOSE_ENV_NAME pip install --no-cache-dir cellpose[distributed]==4.0.4

# Clean up conda cache
RUN conda clean -a -y

# ------------------------------------------------------------------------------
# Application Code & Entrypoint
# ------------------------------------------------------------------------------
WORKDIR /app
# >>> Ensure you are copying the correct wrapper script <<<
COPY run.py /app/run.py
COPY descriptor.json /app/descriptor.json

# >>> Ensure Entrypoint points to the correct script <<<
ENTRYPOINT ["/opt/conda/envs/cytomine_py37/bin/python", "/app/run.py"]
