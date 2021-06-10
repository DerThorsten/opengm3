FROM continuumio/anaconda3:5.0.1

# install sys requirements
RUN apt-get update
RUN apt-get install -y build-essential git

# clone in repo
RUN git clone https://github.com/DerThorsten/opengm

# update conda
RUN conda update --yes conda

# go into repo
RUN cd opengm && conda env create -f opengm-dev-requirements.yml

# activate
RUN echo "source activate opengm-dev-requirements" > ~/.bash
ENV PATH /opt/conda/envs/opengm-dev-requirements/bin:$PATH

# activate
#RUN conda activate opengm-dev-requirements

# make build dir
RUN cd opengm && mkdir build

# run cmake
RUN  cd opengm/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=$(conda info --base)/envs/opengm-dev-requirements

# build and install project
RUN cd opengm/build && \
    make -j2 install