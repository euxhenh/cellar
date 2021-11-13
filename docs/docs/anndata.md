---
title: AnnData
layout: default
nav_order: 6
---

# AnnData Structure
{: .no_toc}

[AnnData](https://anndata.readthedocs.io/en/latest/) is the internal object
that Cellar uses for storing annotations and exporting session files.
Here we describe the keys that Cellar uses to store its results.
{: .fs-6 .fw-300}

---

For an overview of the AnnData parameters check the AnnData
[API](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html).

* `adata.X`: The main data matrix

* `adata.obs['labels']`: Cluster IDs, integer
* `adata.obs['annotations']`: Cluster annotations (cell types), strings

* `adata.var['gene_symbols']`: HGNC gene symbols, strings

* `adata.uns['labels']`: Clustering settings, dictionary
* `adata.uns['x_emb']`: Dimensionality reduction settings, dictionary
* `adata.uns['x_emb_2d']`: 2D dimensionality reduction settings, dictionary
* `adata.uns['neighs']`: Neighbors graph settings, dictionary

* `adata.obsm['x_emb']`: Reduced representation of the data, array
* `adata.obsm['x_emb_2d']`: 2D embeddings of the data, array
* `adata.obsm['proteins']`: Protein data for CITE-seq, array
* `adata.obsm['genes']`: Mapped gene expression data from CITE-seq for CODEX

* `adata.obsp['neighs']`: Neighbors adjacency matrix, sparse matrix