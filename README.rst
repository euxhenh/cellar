.. -*- mode: rst -*-

|PythonVersion|_ |MITLicense|_ |Website|_

.. |PythonVersion| image:: https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue
.. _PythonVersion: https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue
.. |MITLicense| image:: https://img.shields.io/badge/License-MIT-blue
.. _MITLicense: https://raw.githubusercontent.com/euxhenh/CellarV/master/LICENSE.txt
.. |Website| image:: https://img.shields.io/website-up-down-green-red/http/shields.io
.. _Website: https://data.test.hubmapconsortium.org/app/cellar

.. |PythonMinVersion| replace:: 3.7

.. image:: https://raw.githubusercontent.com/euxhenh/CellarV/master/assets/cellar-logo.png
  :width: 400
  :target: https://data.test.hubmapconsortium.org/app/cellar

**Cellar** is an interactive tool for analyzing single-cell omics data. Cellar
is built in Python using the `Dash <https://plotly.com/dash/>`__ framework
and relies on several open-source packages.

The app is developed and actively maintained by the
`Systems Biology Group <http://www.sb.cs.cmu.edu/>`__ at
`Carnegie Mellon University <https://www.cmu.edu/>`__. Our web-server
running Cellar can be accessed
`here <https://data.test.hubmapconsortium.org/app/cellar>`__. See below
for a local installation.

An accompanying preprint and supplementary files can be accessed through
`bioRxiv <https://www.biorxiv.org/content/10.1101/2021.03.19.436162v1?rss=1>`__.

The `documentation <https://euxhenh.github.io/CellarV/>`__
includes details on how to use Cellar and the data types
it supports. These include but are not limited to scRNA-seq, scATAC-seq,
CODEX, SNARE-seq, sciRNA-seq, Visium. Cellar supports preprocessing,
dimensionality reduction, clustering, DE gene testing, enrichment analysis,
cluster and gene visualization modules, projection to spatial tiles,
label transfer, and semi-supervised clustering among others. The documentation
also contains several written tutorials.
`Video tutorials
<https://www.youtube.com/playlist?list=PL5sLSLkTYpWgfBQ0M8ObfBIqDMAzx0-D2>`__
are also available.

Links
_____

- Cellar Web Server: https://data.test.hubmapconsortium.org/app/cellar
- Official Source Code Repository: https://github.com/euxhenh/CellarV
- Documentation & Tutorials: https://euxhenh.github.io/CellarV/
- Video Tutorials: https://www.youtube.com/playlist?list=PL5sLSLkTYpWgfBQ0M8ObfBIqDMAzx0-D2
- Issue Tracker: https://github.com/euxhenh/CellarV/issues
- Preprint: https://www.biorxiv.org/content/10.1101/2021.03.19.436162v1?rss=1

Local Installation
__________________

Docker Installation
~~~~~~~~~~~~~~~~~~~

Probably the easiest way to install Cellar locally is using ``Docker``.
The image name is ``euxhen/cellarv`` and can be pulled with::

    docker pull euxhen/cellarv

After the pull is complete, running Cellar is as simple as::

    docker run --rm -p 8050:8050 euxhen/cellarv

and visiting ``localhost:8050`` on your web browser.

Manual Installation
~~~~~~~~~~~~~~~~~~~

A manual installation involves cloning the Cellar repository and installing
the necessary Python and R packages. To run Cellar you will need at least
**Python 3.7** and **R 4.0**. We recommend using a ``Conda`` environment
for installing the dependencies.

Cellar's source code can be downloaded using::

    git clone https://github.com/euxhenh/CellarV

Dependencies
++++++++++++

Cellar requires the following Python packages:

- dash (== 1.21.0), dash-bootstrap-components (== 0.13.0), dash-bio (>= 0.7.1)
- numpy (>= 1.21.1), scipy (>= 1.6.3)
- scikit-learn (>= 0.24.2), scikit-learn-extra (>= 0.2.0), scikit-image (>= 0.18.1)
- pandas (>= 1.3.1), anndata (>= 0.7.6), scanpy (>= 1.8.1)
- matplotlib (>= 3.4.2), plotly (>= 5.1.0), python-kaleido (>= 0.2.1)
- tifffile (>= 2019.7.26.2)
- igraph (>= 0.9.4), leidenalg (>= 0.8.7)
- umap-learn (>= 0.5.1)
- faiss-cpu (>= 1.7.1)
- joblib (>= 1.0.1)
- zodb (>= 5.6.0)
- r-base (>= 4.1.0), r-stringi (>= 1.7.3), rpy2 (>= 3.4.5), anndata2ri (>= 1.0.6)
- pydiffmap (>= 0.2.0.1)
- diffxpy (>= 0.7.4), gseapy (>= 0.10.5), pyensembl (>=1.9.4), bintogene (>= 1.23)
- kneed (>= 0.7.0)

and the following R packages

- SingleR (>= 1.6.1)
- aertslab/cisTopic (>= 0.3.0)
- STvEA (>= 0.2.0)

These R packages can be installed using the provided python script
``install_Rdeps.py``.

The entry point of the application is ``main.py`` and can be run with::

    python main.py

Contributing
____________

We welcome code contributions as well as feature requests. To request
new features please raise an issue in the links provided above or directly
`send us an email <mailto:ehasanaj@cs.cmu.edu>`__.