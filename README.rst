.. -*- mode: rst -*-

|PythonVersion|_ |MITLicense|_ |Website|_ |DOI|_

.. |PythonVersion| image:: https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue
.. _PythonVersion: https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue
.. |MITLicense| image:: https://img.shields.io/badge/License-MIT-blue
.. _MITLicense: https://raw.githubusercontent.com/euxhenh/cellar/main/LICENSE.txt
.. |Website| image:: https://img.shields.io/website-up-down-green-red/http/shields.io
.. _Website: https://cellar.cmu.hubmapconsortium.org/app/cellar
.. |DOI| image:: https://zenodo.org/badge/372980254.svg
.. _DOI: https://zenodo.org/badge/latestdoi/372980254

.. |PythonMinVersion| replace:: 3.7

.. image:: https://raw.githubusercontent.com/euxhenh/cellar/main/assets/cellar-logo.png
  :width: 400
  :target: https://cellar.cmu.hubmapconsortium.org/app/cellar

**Cellar** is an interactive tool for analyzing single-cell omics data. Cellar
is built in Python using the `Dash <https://plotly.com/dash/>`__ framework
and relies on several open-source packages.

The app is developed and actively maintained by the
`Systems Biology Group <http://www.sb.cs.cmu.edu/>`__ at
`Carnegie Mellon University <https://www.cmu.edu/>`__. Our web-server
running Cellar can be accessed
`here <https://cellar.cmu.hubmapconsortium.org/app/cellar>`__. See below
for a local installation.

An accompanying preprint and supplementary files can be accessed through
`bioRxiv <https://www.biorxiv.org/content/10.1101/2021.03.19.436162v1?rss=1>`__.

The `documentation <https://euxhenh.github.io/cellar/>`__
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

- Cellar Web Server: https://cellar.cmu.hubmapconsortium.org/app/cellar
- Official Source Code Repository: https://github.com/euxhenh/cellar
- Documentation & Tutorials: https://euxhenh.github.io/cellar/
- Video Tutorials: https://www.youtube.com/playlist?list=PL5sLSLkTYpWgfBQ0M8ObfBIqDMAzx0-D2
- Issue Tracker: https://github.com/euxhenh/cellar/issues
- Preprint: https://www.biorxiv.org/content/10.1101/2021.03.19.436162v1?rss=1
- DOI: https://zenodo.org/badge/latestdoi/372980254

Local Installation
__________________

Docker Installation
~~~~~~~~~~~~~~~~~~~

Probably the easiest way to install Cellar locally is using ``Docker``.
The image name is ``euxhen/cellar`` and can be pulled with::

    docker pull euxhen/cellar

After the pull is complete, running Cellar is as simple as::

    docker run --rm -p 8050:8050 euxhen/cellar

and visiting ``localhost:8050`` on your web browser.

Manual Installation
~~~~~~~~~~~~~~~~~~~

A manual installation involves cloning the Cellar repository and installing
the necessary Python and R packages. To run Cellar you will need at least
**Python 3.7** and **R 4.0**. We recommend using a ``Conda`` environment
for installing the dependencies. For a full list of dependencies and
installation instructions please refer to the
`documentation <https://euxhenh.github.io/cellar/>`__.

Contributing
____________

We welcome code contributions as well as feature requests. To request
new features please raise an issue in the links provided above or directly
`send us an email <mailto:ehasanaj@cs.cmu.edu>`__.