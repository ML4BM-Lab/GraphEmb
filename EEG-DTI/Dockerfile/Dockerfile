FROM tensorflow/tensorflow:1.15.0-gpu-py3

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils pkg-config git nano mlocate libprotobuf-dev wget 


RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install bidict  scikit-learn numpy tqdm matplotlib


WORKDIR /
RUN git clone https://github.com/MedicineBiology-AI/EEG-DTI.git
WORKDIR /EEG-DTI 

# ADD ./DTI_data/ /EEG-DTI/DTI_data/
# add here the codes that change 
# for three cases: default, splits, rmsd

# remove for sanity !
RUN rm -r ./DTI_data 
RUN rm /EEG-DTI/decagon/deep/minibatch.py 


# add new files

# For Default run
ADD ./main_modified_eegdti.py /EEG-DTI
ADD ./minibatch_mod.py /EEG-DTI/decagon/deep/

# For split run !
ADD ./main_modified_eegdti_splits.py /EEG-DTI
ADD ./minibatch_splits.py /EEG-DTI/decagon/deep/
ADD ./optimizer_splits.py /EEG-DTI/decagon/deep/

