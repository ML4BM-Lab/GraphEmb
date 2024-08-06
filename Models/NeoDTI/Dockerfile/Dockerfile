FROM tensorflow/tensorflow:1.2.0
# FROM tensorflow/tensorflow:1.2.0-gpu
# FROM tensorflow/tensorflow:1.15.4

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils \ 
pkg-config git nano mlocate libprotobuf-dev

WORKDIR /
RUN git clone https://github.com/FangpingWan/NeoDTI.git
WORKDIR /NeoDTI 

RUN pip install numpy==1.14.0 scikit-learn==0.19.0 tflearn==0.3.2

#Test run 
WORKDIR /NeoDTI/data
RUN unzip data.zip
# WORKDIR /NeoDTI/src
# RUN python NeoDTI_cv.py


# Tensorflow (tested on version 1.0.1 and version 1.2.0)
# tflearn
# numpy (tested on version 1.13.3 and version 1.14.0)
# sklearn (tested on version 0.18.1 and version 0.19.0)


