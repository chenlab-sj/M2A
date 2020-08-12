FROM ubuntu:19.10 as builder
  
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install wget -y && \
    apt-get install parallel -y && \
    apt-get install g++ -y && \ 
    rm -r /var/lib/apt/lists/*

RUN wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O miniconda.sh && \
    /bin/bash miniconda.sh -b -p /opt/conda/ && \
    rm miniconda.sh

RUN conda update -n base -c defaults conda -y && \
    conda install \
    -c conda-forge \
    -c bioconda \
    python=3.6.5 \
    coreutils=8.31 \
    numpy=1.17.1 \
    pandas=0.25.1 \
    h5py=2.9.0 \ 
    keras=2.2.4 \ 
    pyBigWig \
    scikit-learn=0.20.0 \
    mkl-service \
    ipywidgets \
    psutil=5.6.1 \
    keras-preprocessing=1.0.5 \
    keras-applications=1.0.6 \
    scipy=1.3.1 \
    blas \
    tensorflow=1.10.0 \
    matplotlib \
    -y && \
    conda clean --all -y

RUN  pip install pandarallel

COPY 0_PromoterDefinitions /opt/M2A/0_PromoterDefinitions
COPY 1_ResponseVariable /opt/M2A/1_ResponseVariable
COPY 2_MethylationFeatures /opt/M2A/2_MethylationFeatures
COPY 3_Combine /opt/M2A/3_Combine
COPY 4_RunModel /opt/M2A/4_RunModel
COPY 5_TransferLearning /opt/M2A/5_TransferLearning
COPY cwl /opt/M2A/cwl
    
ENV PATH /opt/M2A/0_PromoterDefinitions:/opt/M2A/1_ResponseVariable:/opt/M2A/2_MethylationFeatures:/opt/M2A/3_Combine:/opt/M2A/4_RunModel:/opt/M2A/5_TransferLearning:/opt/M2A/cwl:${PATH}
