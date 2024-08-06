FROM tensorflow/tensorflow:1.15.5-gpu-py3
# FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
	apt-utils pkg-config build-essential git nano mlocate wget python3-pip

RUN apt-get install apt-transport-https ca-certificates -y \
	&& update-ca-certificates

RUN pip3 install --upgrade pip


RUN git clone https://github.com/kexinhuang12345/MolTrans.git
WORKDIR /MolTrans
RUN pip3 install -r requirements.txt


