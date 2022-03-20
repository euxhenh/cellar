#!/bin/bash

# Use exact versions for r-base and rpy2
conda install -c conda-forge -c bioconda r-base=4.1.2 \
	rpy2=3.4.5 \
	r-devtools \
	r-ggplot2 \
	r-fitdistrplus \
	r-lda \
	r-plyr \
	r-text2vec \
	r-rann \
	r-tidyr \
	r-mclust \
	r-igraph \
	r-cluster \
	bioconductor-biocinstaller \
	bioconductor-singler \
	bioconductor-flowcore \
	bioconductor-rcistarget \
	bioconductor-rtracklayer \
	bioconductor-aucell

conda install -c conda-forge -c plotly -c bioconda dash=2.2.0 \
	dash-bootstrap-components=1.0.3 \
	dash-bio \
	scikit-learn \
	scikit-learn-extra \
	scikit-image \
	scanpy \
	matplotlib \
	python-kaleido \
	tifffile \
	leidenalg \
	faiss-cpu \
	zodb \
	plotly \
	anndata2ri \
	gseapy \
	pyensembl \
	gunicorn

pip install pydiffmap diffxpy bintogene
