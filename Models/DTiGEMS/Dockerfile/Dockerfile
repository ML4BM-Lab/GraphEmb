FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils pkg-config build-essential python3.10 python3-pip python3-setuptools python3-dev wget software-properties-common libssl-dev curl libcurl4-openssl-dev liblapack-dev libopenblas-base libopenblas-dev gfortran libpng-dev libcairo2-dev libxml2-dev libfontconfig1-dev libmagick++-dev libpq-dev libhdf5-dev git nano mlocate


WORKDIR /
RUN git clone https://github.com/MahaThafar/Drug-Target-Interaction-Prediciton-Method.git
WORKDIR /Drug-Target-Interaction-Prediciton-Method/DTiGEMS+
RUN pip3 install scikit-learn==0.24
RUN pip3 install networkx
RUN pip3 install gensim
RUN pip3 install -r requirements.txt


RUN mkdir BCK
#make a backUp
RUN cp *.py ./BCK/
RUN sed -i 's/fit_sample/fit_resample/g' *.py

# make a backup and replace the map behavoiur to adpat to python3
RUN cp /usr/local/lib/python3.8/dist-packages/gensim/models/word2vec.py /usr/local/lib/python3.8/dist-packages/gensim/models/word2vec_bck.py
RUN sed -i 's/total_words += len(sentence)/total_words += len(list(sentence))/g' /usr/local/lib/python3.8/dist-packages/gensim/models/word2vec.py
RUN sed -i 's/return sum(len(sentence) for sentence in job)/return sum(len(list(sentence)) for sentence in job)/g' /usr/local/lib/python3.8/dist-packages/gensim/models/word2vec.py

# make a backup and replace some arguments to adapt to new naming
RUN cp /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/n2v_mainFunctions.py /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/n2v_mainFunctions_bck.py
RUN sed -i 's/model = Word2Vec(walks, size=args.dimensions, window=args.window_size, min_count=0, sg=1, workers=args.workers, iter=args.iter)/model = Word2Vec(walks, vector_size=args.dimensions, window=args.window_size, min_count=0, sg=1, workers=args.workers, epochs=args.iter)/g' /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/n2v_mainFunctions.py

# make a backup and replace some arguments to adapt the print function to python3
RUN cp /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/node2vec.py /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/node2vec_bck.py
RUN sed -i "s/print 'Walk iteration:'/print('Walk iteration:')/g" /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/node2vec.py
RUN sed -i "s/print str(walk_iter+1), '\/', str(num_walks)/print(str(walk_iter+1), '\/', str(num_walks))/g" /Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/node2vec.py
