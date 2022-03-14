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

View Cellar's [conda environment](https://github.com/euxhenh/cellar/blob/main/env.yml) for a full list of libraries and their
versions. Below we describe only the libraries that Cellar depends on
explicitly.

## Python Dependencies

- [dash](https://plotly.com/dash/), [dash-bootstrap-components](https://dash-bootstrap-components.opensource.faculty.ai/), [dash-bio](https://dash.plotly.com/dash-bio)
- [numpy](https://numpy.org/), [scipy](https://scipy.org/)
- [scikit-learn](https://scikit-learn.org/), [scikit-learn-extra](https://scikit-learn-extra.readthedocs.io/en/stable/), [scikit-image](https://scikit-image.org/)
- [pandas](https://pandas.pydata.org/), [anndata](https://anndata.readthedocs.io/en/latest/), [scanpy](https://scanpy.readthedocs.io/en/stable/)
- [matplotlib](https://matplotlib.org/), [plotly](https://plotly.com/), [python-kaleido](https://pypi.org/project/kaleido/), [tifffile](https://pypi.org/project/tifffile/)
- [igraph](https://igraph.org/), [leidenalg](https://leidenalg.readthedocs.io/en/latest/)
- [umap-learn](https://umap-learn.readthedocs.io/en/latest/), [pydiffmap](https://pydiffmap.readthedocs.io/en/master/)
- [faiss-cpu](https://github.com/facebookresearch/faiss), [joblib](https://joblib.readthedocs.io/en/latest/)
- [zodb](https://zodb.org/en/latest/)
- [r-base](https://www.r-project.org/), [r-stringi](https://cran.r-project.org/web/packages/stringi/index.html), [rpy2](https://rpy2.github.io/), [anndata2ri](https://github.com/theislab/anndata2ri)
- [diffxpy](https://diffxpy.readthedocs.io/en/latest/), [gseapy](https://gseapy.readthedocs.io/en/latest/introduction.html), [pyensembl](https://readthedocs.org/projects/pyensembl/downloads/pdf/latest/), [bintogene](https://github.com/euxhenh/BinToGene)

## R Dependencies

- [SingleR](https://github.com/dviraran/SingleR)
- [cisTopic](https://github.com/aertslab/cisTopic)
- [STvEA](https://github.com/CamaraLab/STvEA)

NOTE: We prefer to install R (r-base) via conda as seen above.
The few remaining R libraries can be installed
by using the [install_Rdeps.py](https://github.com/euxhenh/cellar/blob/main/install_Rdeps.py) script while inside the conda environment.

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