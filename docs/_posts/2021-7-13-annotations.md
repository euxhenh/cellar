---
title: 'ANNOTATIONS'

category: User Interface
layout: null
type: ui
---
This menu is where you can add new subsets of cells or annotate existing ones.
<br>
![Annotations](images/annotations.png)
<br>


* **Annotation**
    * This is where you select a subset/cluster of cells and annotate them (typically with cell types).
    <span class="warn">!!</span> If the subset contains cells from multiple
    clusters, the labels of those cells will be lost and replaced by
    the new type. There will also be a table under the input field for showing the cluster ID and the corresponding annotation.

* **Store Subset**
    * This allows you to select any subset of cells from the plot and give
    them a unique identifier. To select a subset of cells, first make sure the
    <span class='mbox'>Lasso Select</span> tool is selected (see image below).
    Then, by clicking anywhere on the plot and dragging your mouse around,
    you should be able to select a group of cells
    (double-click on the plot to reset the selection).<br>
    ![Lasso Select](images/lasso-select.png)<br>
    Once you are happy with your selection, click <span class="mbutton">Store</span>. These identifiers
    can be used to find DE genes for the corresponding subset under
    the <span class='mbox'>Analysis</span> menu, and you can also assign a cell
    type to the subset.
    Subsets that start with `Cluster` are reserved for clusters and no new
    subsets can start with this prefix.


* **Merge**
    * This allows the users to select a subset and store it as a new cluster. The color in the main plot will also change accordingly.
    * <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    Updates the `adata.obs['labels']`, `adata.uns['cluster_info']` and
    `adata.uns['cluster_names']` keys.