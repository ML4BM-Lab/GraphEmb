FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
	apt-utils pkg-config build-essential python3.10 python3-pip \ 
	python3-setuptools python3-dev wget software-properties-common \ 
	libssl-dev curl libcurl4-openssl-dev liblapack-dev libopenblas-base \
	libopenblas-dev gfortran libpng-dev libcairo2-dev libxml2-dev \
	libfontconfig1-dev libmagick++-dev libpq-dev libhdf5-dev git nano mlocate


WORKDIR /
RUN git clone https://github.com/MahaThafar/DTi2Vec.git
WORKDIR /DTi2Vec
RUN python3 -m pip install -r requirements.txt

RUN cp /DTi2Vec/DTi2vec_main.py /DTi2Vec/DTi2vec_main_bck.py
COPY DTi2vec_main.py /DTi2Vec/DTi2vec_main.py
COPY DTi2vec_newDrug_seting_generatedEMBED.py /DTi2Vec/DTi2vec_newDrug_seting_generatedEMBED.py