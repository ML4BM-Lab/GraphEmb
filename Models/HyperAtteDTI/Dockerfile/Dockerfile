FROM python:3.6.8

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils \
    pkg-config build-essential wget software-properties-common libssl-dev curl \
    libcurl4-openssl-dev liblapack-dev libopenblas-base libopenblas-dev \ 
    gfortran libpng-dev libcairo2-dev libxml2-dev libfontconfig1-dev \ 
    libmagick++-dev libpq-dev libhdf5-dev git nano locate


RUN git clone https://github.com/zhaoqichang/HpyerAttentionDTI.git
WORKDIR /HpyerAttentionDTI

RUN pip install --upgrade pip


RUN pip install torch
RUN pip install numpy
RUN pip install sklearn
RUN pip install tqdm
RUN pip install tensorboardX
RUN pip install prefetch_generator
RUN pip install six
