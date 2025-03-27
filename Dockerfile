FROM condaforge/miniforge3:24.11.3-2 

# System dependencies
RUN apt-get update && apt-get install -y \
    git \
    libgeos-dev \
    vim \
    gcc \
    g++ \
    zlib1g-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libtiff-dev \
    libwebp-dev \
    libgl1-mesa-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /app

RUN mkdir -p /root/.cellpose/models && \
    cd /root/.cellpose/models && \
    wget http://www.cellpose.org/models/nuclei_0 && \
    wget http://www.cellpose.org/models/nuclei_1 && \
    wget http://www.cellpose.org/models/nuclei_2 && \
    wget http://www.cellpose.org/models/nuclei_3 && \
    wget http://www.cellpose.org/models/size_nuclei_0.npy && \
    wget --no-check-certificate https://www.cellpose.org/models/cyto_0 && \
    wget --no-check-certificate https://www.cellpose.org/models/cyto_1 && \
    wget --no-check-certificate https://www.cellpose.org/models/cyto_2 && \
    wget --no-check-certificate https://www.cellpose.org/models/cyto_3 && \
    wget --no-check-certificate https://www.cellpose.org/models/size_cyto_0.npy && \
    wget --no-check-certificate https://www.cellpose.org/models/cytotorch_0 && \
    wget --no-check-certificate https://www.cellpose.org/models/cytotorch_1 && \
    wget --no-check-certificate https://www.cellpose.org/models/cytotorch_2 && \
    wget --no-check-certificate https://www.cellpose.org/models/cytotorch_3 && \
    wget --no-check-certificate https://www.cellpose.org/models/size_cytotorch_0.npy && \
    wget --no-check-certificate https://www.cellpose.org/models/nucleitorch_0 && \
    wget --no-check-certificate https://www.cellpose.org/models/nucleitorch_1 && \
    wget --no-check-certificate https://www.cellpose.org/models/nucleitorch_2 && \
    wget --no-check-certificate https://www.cellpose.org/models/nucleitorch_3 && \
    wget --no-check-certificate https://www.cellpose.org/models/size_nucleitorch_0.npy

# Install BIAFLOWS requirements in base environment
RUN conda install -y python=3.7

# ------------------------------------------------------------------------------
# Install Cytomine python client
RUN git clone https://github.com/cytomine-uliege/Cytomine-python-client.git && \
    cd Cytomine-python-client && git checkout tags/v2.7.3 && pip install . && \
    cd .. && rm -rf Cytomine-python-client

# ------------------------------------------------------------------------------
# Install BIAFLOWS utilities (annotation exporter, compute metrics, helpers,...)
RUN apt-get update && apt-get install libgeos-dev -y && apt-get clean
RUN git clone https://github.com/Neubias-WG5/biaflows-utilities.git && \
    cd biaflows-utilities/ && git checkout tags/v0.9.2 && pip install . && \
    cd ..

# install utilities binaries
RUN chmod +x /app/biaflows-utilities/bin/*
RUN cp /app/biaflows-utilities/bin/* /usr/bin/ && \
    rm -rf /app/biaflows-utilities



# Create separate conda environment for micro-sam
RUN conda create -n cellpose -y python=3.10
RUN conda run -n cellpose pip install cellpose[distributed]
RUN conda run -n cellpose pip install zarr

RUN pip install numpy==1.19.4
RUN pip install numba==0.50.1

# Add wrapper and descriptor
ADD wrapper.py /app/wrapper.py
ADD descriptor.json /app/descriptor.json

ENTRYPOINT ["python3.7","/app/wrapper.py"]