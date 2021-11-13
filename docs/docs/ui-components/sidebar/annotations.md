---
layout: default
title: Annotations
nav_order: 3
parent: Sidebar
usemathjax: true
grand_parent: UI Components
---

# Annotations
{: .no_toc }

The end goal of the Cellar workflow is to assign a cell type to each sample
point. This is achieved via the clustering step above which finds clusters
of "similar" cells, which are then analyzed as a group. To assign a cell
type to a cluster, first select that cluster from the dropdown menu and
type the annotation you wish to assign. The annotations table will only
list those clusters which have been assigned some cell type.

<img src="../../../images/annotations.png" width="400" class="center"/>

Hovering at any point in the scatter plot will show the annotation for
that point or an empty string if there is none.

In the AnnData file, these annotations are stored under
`adata.obs['annotations']`. Annotations are limited to a 200 character length.