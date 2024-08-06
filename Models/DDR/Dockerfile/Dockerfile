FROM python:2.7.18


ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
	apt-utils pkg-config build-essential\ 
	wget software-properties-common \ 
	libssl-dev curl libcurl4-openssl-dev liblapack-dev libopenblas-base \
	libopenblas-dev gfortran libpng-dev libcairo2-dev libxml2-dev \
	libfontconfig1-dev libmagick++-dev libpq-dev libhdf5-dev git nano mlocate


WORKDIR /
RUN git clone https://bitbucket.org/RSO24/ddr.git
WORKDIR /ddr
RUN pip install scikit-learn numpy pandas
