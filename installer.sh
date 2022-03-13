#!/bin/bash

conda install -c r rpy2 r-cluster
conda install -c conda-forge r-devtools \
	r-ggplot2 \
	r-fitdistrplus \
	r-lda \
	r-plyr \
	r-text2vec \
	r-rann \
	r-tidyr \
	r-mclust \
	r-igraph

conda install -c bioconda bioconductor-biocinstaller \
	bioconductor-singler \
	bioconductor-flowcore \
	bioconductor-rcistarget \
	bioconductor-rtracklayer \
	bioconductor-aucell

conda install -c conda-forge dash dash-bootstrap-components dash-bio
conda install -c conda-forge scikit-learn
conda install -c conda-forge scikit-learn-extra
conda install -c conda-forge scikit-image
conda install -c conda-forge scanpy
conda install -c conda-forge matplotlib
conda install -c plotly plotly
conda install -c conda-forge python-kaleido
conda install -c conda-forge tifffile
conda install -c conda-forge leidenalg
conda install -c conda-forge faiss-cpu
conda install -c conda-forge zodb
conda install -c bioconda anndata2ri gseapy pyensembl

pip install pydiffmap diffxpy bintogene
