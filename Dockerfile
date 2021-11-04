FROM nvidia/cuda:11.1.1-devel-centos8

LABEL maintainer="k.masaru.k@gmail.com"

# Install miniconda
WORKDIR /
RUN \
  yum install -y wget git && \
  wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh -O Miniconda3-py37_4.10.3-Linux-x86_64.sh && \
  bash Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p miniconda3 && \
  rm -f Miniconda3-py37_4.10.3-Linux-x86_64.sh
RUN \
  source /miniconda3/bin/activate && \
  conda update conda && \
  conda update --all

# Install MENTR
WORKDIR /
# TODO:
# RUN git clone https://github.com/koido/MENTR.git
COPY ./mentr_env.txt /MENTR/
COPY ./resources/example/ /MENTR/resources/example/
COPY ./resources/mutgen_cmn/ /MENTR/resources/mutgen_cmn/
COPY ./resources/trained/ /MENTR/resources/trained/
COPY ./resources/F5.cage_cluster.hg19.info.tsv.gz /MENTR/resources/
COPY ./resources/supp_table_10.sample_ontology_information.tsv.gz /MENTR/resources/
COPY ./src/ /MENTR/src/
COPY ./bin/ /MENTR/bin/
COPY ./LICENSE /MENTR/
COPY ./README.md /MENTR/
WORKDIR /MENTR

RUN \
   /miniconda3/bin/conda create \
     --name mentr \
     -c conda-forge -c bioconda \
     --file mentr_env.txt python=3.6.4

ENV PATH /miniconda3/bin:${PATH}
ENV PATH /miniconda3/envs/mentr/bin/:${PATH}

CMD ["/bin/bash"]
