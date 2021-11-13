---
title: File Formats
layout: default
nav_order: 4
---

# Accepted File Formats
{: .no_toc}

Cellar accepts `.h5ad`, `.csv`, and `.gz` files.
{: .fs-6 .fw-300}

---

Internally, Cellar uses the [AnnData](https://anndata.readthedocs.io/en/latest/)
data structure. These objects can be
uploaded to Cellar as `.h5ad` files. We also accept
`.csv` files or feature barcode sequences in
[cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) or
[spaceranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices) format as `.gz` files.
In the case of `.h5ad` or `.csv` files,
please make sure the data follows the *cell* x *feature*
format, i.e., rows of the matrix correspond to cells and columns
to features (genes, proteins).
Additionally, row and columns names must be provided with the file.

- Example cell barcode: **TTACCATTCCAGCACG**. <br>
  Example gene name (ensembl format or HGNC symbol): **ENSG00000211699** or **TRGV3**.

We recommend uploading the data as an `.h5ad` file.
An AnnData object stores expression matrices and annotations
efficiently while also providing convenient methods for accessing
and storing any additional unstructured objects. Furthermore,
an AnnData object can be highly compressed, thus saving precious
upload (download) time.

Here we provide a Python snippet that can
convert a `.csv` file into a maximally compressed AnnData object that
can be uploaded to Cellar. All you need is the `anndata` Python package
that can be installed via `pip install anndata`. The following code
reads a csv file `path/to/your/csv/file.csv` and writes it to another
file `new/file/name.h5ad`.

```python
import anndata
adata = anndata.read_csv('path/to/your/csv/file.csv')
adata.write('new/file/name.h5ad', compression=9)
```

The file can now be uploaded to Cellar. This procedure may
reduce a dataset's size from several GBs to just a few MBs.