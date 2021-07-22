---
title: 'CLUSTERING'

category: User Interface
layout: null
type: ui
---
This menu consists of three tabs for obtaining labels: **UNSUPERVISED**, **SEMI-SUPERVISED**, and **LABEL TRANSFER**. In short, **UNSUPERVISED** allows the user to cluster the dimensionality reduced data and generate labels for cells. **SEMI-SEPERVISED** uses already labeled data. In this tab, the user can fix cell labels with high certainty and re-cluster the uncertain cells. **LABEL TRANSFER** should be used under dual mode and also uses dimensionality reduced data. It labels the activated dataset based on the deactivated one.



* **UNSUPERVISED**
    * The clustering dropdown menu includes
    [Leiden](https://www.nature.com/articles/s41598-019-41695-z) (default)
    and several other classical algorithms such as KMeans and Spectral
    Clustering.
      * Leiden does not require the user to manually set the number
        of clusters, however, this number can be controlled via the
        <span class='mbox'>Resolution</span>
        parameter (default: 1). Higher resolution means more clusters and lower
        resolution means less clusters. You may want to experiment with this
        value until you get a desirable cluster configuration. Leiden is a
        community detection algorithm which takes as input a connectivity graph.
        We construct this graph via a nearest neighbors approach. The
        <span class='mbox'>Number of neighbors for the Graph</span> hyperparameter can also
        be tuned by the user.
      * The remaining algorithms share a similar setup. Unlike Leiden, these
        algorithms require the user to specify the number of clusters in advance.
        To simplify tuning this hyperparameter, we allow the users to specify
        a range (or list) of values. Doing this will run the selected
        algorithm once for each value and will automatically select the number
        of clusters that achieves the highest score (default:
        [Silhouette Score](https://en.wikipedia.org/wiki/Silhouette_(clustering))).
        To specify a range, use the notation $$(\text{start}, \text{end}+1,
        \text{step})$$.
        For example, $$(3, 7, 1)$$ will run four instances of the selected
        algorithm, one for each number of clusters in $$[3, 4, 5, 6]$$.
        Alternatively, you can use a list of values. (Note the distinction
        betweeen parentheses ( ) and square brackets [ ]) E.g., $$[4, 8, 16]$$ will
        spawn three instances of the algorithm. Finally, you could use a single
        integer value.
      * There is another special algorithm called Uncertainty clustering. You can only run it after cells are labeled by other algorithms. Uncertainty clustering generates a new cluster (cluster ID -1) with cells having high uncertainty score. The uncertainty score is computed according to cells' distance to the cluster centers in the dimensionality reduced space. After this, re-clustering the new generated uncertain cluster using Constrained Leiden could often improve the clustering results (getting closer to the ground truth) if the original labels are obtained using algorithms based on local neighbors (e.g. Leiden). Intuitively, the improvements come from incorporating the information about the cluster centers, which is not considered in Leiden. 
        
      [//]:# "* An [Ensemble](https://github.com/GGiecold/Cluster_Ensembles) algorithm"
      [//]:# "  based on Hypergraph Partitioning has been added which allows the"
      [//]:# "  selection and integration of several clustering algorithms into an"
      [//]:# "  ensemble."

    * <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    This step populates the
    `adata.obs['labels']`, `adata.uns['cluster_info']` and
    `adata.uns['cluster_names']` keys.



* **SEMI-SEPERVISED**
    * Cellar incorporates an additional semi-supervised clustering algorithm.
    This can be useful when the user intervenes in the data by manually
    splitting/merging clusters or after running label transfer (see below).

    * This includes an algorithms [Constrained Leiden] for refining clustering results.
        * Constrained Leiden allows the user to "preserve"
        a user-defined set of clusters. These clusters will be left as is, while
        the remaining clusters will be subject to change. To do this, click the setting button on the right of the dropdown menu and select clusters to preserve using the checkboxes for each cluster in the setting menu. 

    * <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    This step updates the `adata.obs['labels']`, `adata.uns['cluster_info']` and
    `adata.uns['cluster_names']` keys.

* **LABEL TRANSFER**
  * Label transfer should be used under the dual mode. It can be used to transfer
    the labels from a deactivated reference dataset 
    to the current activated one. Currently Cellar supports label transfer with
    [Scanpy Ingest](https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html)
    and [SingleR](https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html).
    <span class="warn">!!<span class="tooltip">Attention</span></span>
    When running Scanpy Ingest, we use only
    the genes which are common in both the active and the reference dataset.
    If you get bad results, make sure the two datasets use the same format
    for the gene names. Label transfer typically takes a few minutes to finish.

  * <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    Updates the `adata.obs['labels']`, `adata.uns['cluster_info']` and
    `adata.uns['cluster_names']` keys.


* **<span class='mbutton'>RUN</span>**
    * Pressing this button will run the selected clustering algorithm. After
    a short delay, each cell in the scatter plot will be colored based on its
    assigned cluster.
    <span class="warn">!!<span class="tooltip">Attention</span></span>
    It is important to
    distinguish between cluster IDs and cluster names. Cluster IDs are simple
    integers, beginning with 0, that refer to the ID of that cluster.
    Cluster names on the other hand, can be
    modified by the user (typically you want to set it equal to the cell
    type comprising that cluster). Initially, every cluster name is
    equal to the cluster ID.