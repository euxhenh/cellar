---
title: 'DIMENSIONALITY REDUCTION'

category: User Interface
layout: null
type: ui
---
This menu consists of two components. The data is first reduced by using
the method chosen under **Obtain Embeddings**. These embeddings are then
used for both the **Obtain 2D Embeddings** (further data reduction) and **CLUSTERING** steps.

* **Dimensionality Reduction**
    * The default method is PCA. Truncated SVD and UMAP
    can also work with sparse matrices.

    * After a method is chosen from the dropdown menu, the user can expand a tab with the settings and documentation of the method by clicking the button on the right.

    * Most methods require you to manually set the number of components
    (default: 40). For a select few (namely, PCA and its variants), we provide
    an <span class='keyword'>Automatic</span> option which will
    determine the best number of components
    for you. This runs the chosen method with a large number of components (200)
    and uses the [Kneedle algorithm](https://ieeexplore.ieee.org/document/5961514)
    to detect the "knee" in the graph of the explained variance.
    <span class="good">GP<span class="tooltip">Good Practice</span></span>
    We recommend running PCA in <span class='keyword'>Automatic</span>
    mode once, take note of the computed number of components, and then
    manually set this value on the next run to save computational time.

    [//]:# "* [cisTopic](https://www.nature.com/articles/s41592-019-0367-1) is a"
    [//]:# "probabilistic framework that uses Latent Dirichlet Allocation to model"
    [//]:# "cis-regulatory topics. We use cisTopic to mainly analyze scATAC-seq data."
    [//]:# "Instead of running the scATAC-seq preprocessing step, you may run "
    [//]:# "cisTopic"
    [//]:# "on the raw cell $$\small\times$$ peak matrix to obtain a cell"
    [//]:# "$$\small\times$$ cisTopic matrix."

    * <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    This step populates the
    `adata.obsm['x_emb']` and `adata.uns['dim_reduction_info']` keys.

* **Obtain 2D Embeddings**
    * This list of algorithms is used to further reduce the embeddings
    obtained via the "Dimensionality Reduction" step to 2D for visualization
    purposes. [UMAP](https://umap-learn.readthedocs.io/en/latest/) generally
    gives good results, but methods like TSNE or PCA should also be considered.

    * Again the button on the right of the dropdown manu allows the user to configure detailed setting of the method.

    * <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    This step populates the
    `adata.obsm['x_emb_2d']` and `adata.uns['visualization_info_2d']` keys.

* **<span class='mbutton'>RUN</span>**
    * Pressing this button will reduce the data and generate a scatter plot.
    This may take a while depending on the size of the data and the
    complexity of the chosen algorithm. Once data reduction is done, you
    should see a scatter plot of your data. Initially all cells belong to
    a single cluster (cluster 0).
    <span class="warn">!!<span class="tooltip">Attention</span></span>
    You might experience a short delay after processing has finished due to
    the time it takes for the plot to render.

