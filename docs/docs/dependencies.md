---
title: Dependencies
layout: default
nav_order: 2
---

# Dependencies
{: .no_toc}

Cellar depends on several Python and R libraries. These can be both installed
in a conda environment.
{: .fs-6 .fw-300}

---

View Cellar's [conda environment](https://github.com/ferrocactus/CellarV/blob/master/env.yml) for a full list of libraries and their
versions. Below we describe only the libraries that Cellar depends on
explicitly.

## Python Dependencies

- [dash](https://plotly.com/dash/) (1.21.0), [dash-bootstrap-components](https://dash-bootstrap-components.opensource.faculty.ai/) (0.13.0), [dash-bio](https://dash.plotly.com/dash-bio) (0.7.1)
- [numpy](https://numpy.org/) (1.21.1), [scipy](https://scipy.org/) (1.6.3)
- [scikit-learn](https://scikit-learn.org/) (0.24.2), [scikit-learn-extra](https://scikit-learn-extra.readthedocs.io/en/stable/) (0.2.0), [scikit-image](https://scikit-image.org/) (0.18.1)
- [pandas](https://pandas.pydata.org/) (1.3.1), [anndata](https://anndata.readthedocs.io/en/latest/) (0.7.6), [scanpy](https://scanpy.readthedocs.io/en/stable/) (1.8.1)
- [matplotlib](https://matplotlib.org/) (3.4.2), [plotly](https://plotly.com/) (5.1.0), [python-kaleido](https://pypi.org/project/kaleido/) (0.2.1), [tifffile](https://pypi.org/project/tifffile/) (2019.7.26.2)
- [igraph](https://igraph.org/) (0.9.4), [leidenalg](https://leidenalg.readthedocs.io/en/latest/) (0.8.7)
- [umap-learn](https://umap-learn.readthedocs.io/en/latest/) (0.5.1), [pydiffmap](https://pydiffmap.readthedocs.io/en/master/) (0.2.0.1)
- [faiss-cpu](https://github.com/facebookresearch/faiss) (1.7.1), [joblib](https://joblib.readthedocs.io/en/latest/) (1.0.1)
- [zodb](https://zodb.org/en/latest/) (5.6.0)
- [r-base](https://www.r-project.org/) (4.1.0), [r-stringi](https://cran.r-project.org/web/packages/stringi/index.html) (1.7.3), [rpy2](https://rpy2.github.io/) (3.4.5), [anndata2ri](https://github.com/theislab/anndata2ri) (1.0.6)
- [diffxpy](https://diffxpy.readthedocs.io/en/latest/) (0.7.4), [gseapy](https://gseapy.readthedocs.io/en/latest/introduction.html) (0.10.5), [pyensembl](https://readthedocs.org/projects/pyensembl/downloads/pdf/latest/) (1.9.4), [bintogene](https://github.com/ferrocactus/BinToGene) (1.23)

## R Dependencies

- [SingleR](https://github.com/dviraran/SingleR) (1.6.1)
- [cisTopic](https://github.com/aertslab/cisTopic) (0.3.0)
- [STvEA](https://github.com/CamaraLab/STvEA) (0.2.0)

NOTE: We prefer to install R (r-base) via conda as seen above.
The few remaining R libraries can be installed
by using the [install_Rdeps.py](https://github.com/ferrocactus/CellarV/blob/master/install_Rdeps.py) script while inside the conda environment.

---

### Description
We use *Dash* open-source for the interface and all the UI components. *NumPy*
is divine, and *Scikit-learn* is at the core of many of the tools
that Cellar uses, including several dimensionality reduction and
clustering methods. *AnnData* is the data structure that we use to store
annotations and is also used as a session file.
*Scanpy* is at the core of the preprocessing tools as well as the Ingest
integration algorithm. *Leidenalg* is used for Leiden clustering, which is the
default clustering algorithm that Cellar uses, and *faiss-cpu* is used for fast
neighbors computation when constructing the neighbors graph for Leiden.
*UMAP-learn* and *pydiffmap* are used for dimensionality reduction. *Rpy2* serves
as an interface for running R functions from within Python. We use *diffxpy*
for differential gene testing and *GSEAPY* for enrichment analysis.
*SingleR* and *STvEA* are two other data integration algorithms,
and *cisTopic* is used to obtain cell-by-topic matrices from scATAC-seq data.
*BinToGene* is used to obtain a cell-by-gene matrix from scATAC-seq data.