# FROM mathworks/matlab:latest
FROM mathworks/matlab-deep-learning:latest


ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils \ 
pkg-config git nano mlocate libprotobuf-dev wget build-essential



WORKDIR /
RUN git clone https://github.com/luoyunan/DTINet.git
WORKDIR /DTINet/src

