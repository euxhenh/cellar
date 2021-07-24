---
title: 'Tutorial: CODEX and 10X Spatial Transcriptomics'

category: Tutorials
layout: null
type: null
---
This is a step-by-step tutorial for analyzing CODEX data and 10x Genomics's Spatial Transcriptomics data. We will take a server dataset `CODEX_Florida_19-003-kymph-node-R2` as example. The analysis of 10X Spatial Transcriptomics follows almost exactly the same steps. We will focus on how to generate spatial tiles for the dataset. Details about other basic analysis based on gene expressions are covered in the **Tutorial: Basic Analysis Pipeline** section.

To begin the analysis, we need to load a spatial dataset. At the top of the page, click <span class='mbutton'>LOAD DATA</span> button in the navigation bar. 
This will expand a panel for loading datasets. Then we can select a server dataset that we are interested in from the dropdown menu. Alternatively, you can upload your own CODEX or 10X dataset and select it in the dropdown menu. See details in the **LOAD DATA** section. After selecting the dataset, click the <span class='mbutton'>LOAD</span> button on the right. 

<br>
<img src="images/codex-load-dataset.png" alt="drawing" width="300"/>
<br>

This will load the dataset and visualize any annotations in it. In our case, the server dataset that we loaded has 2D embedding, cluster labels, and cell type annotations stored. So once the dataset is loaded, the main plot is populated with the 2D visualization, where each point is a cell and the color indicates the clusters. 
<br>
![CODEX Main Plot](images/codex-mainplot.png)
<br>
Number of samples (cells) and number of features (genes) will also be shown on the top of the side bar.
<br>
![CODEX No. Cells and No. Genes](images/codex-cells-genes.png)
<br>
And the cell type annotations are also shown in the **ANNOTATIONS** panel in the side bar.
<br>
![CODEX Cell Types](images/codex-cell-types.png)
<br>
And this is how these results look in Cellar:
<br>
![CODEX Results](images/codex-results.png)
<br>
Please refer to the section **Tutorial: Basic Analysis Pipeline** for details about how to do these analysis and how to do more gene expression visualizations.

Next, we show how to visualize cell type annotations with spatial information. To do so, click the **SPATIAL DATA** tab uncer the main plot and select `CODEX` (or `10X Genomics`, depending on your data type) in the dropdown menu. Then, if you are working with your own dataset, you also need to upload a file including the spatial image and spatial information about the cells by clicking the <span class='mbutton'>UPLOAD</span> button (See the **SPATIAL DATA** section for more details). But we can skip this step since these files are already stored on the server for server datasets. Finally, click the <span class='mbutton'>GENERATE TILE</span> button. After a few seconds, the spatial tile will be generated.
<br>
![CODEX Results](images/codex-spatial-image.png)
<br>
Note that the colors are consistent with the one in the main plot. This means if you change the labels in the main plot (e.g. by re-clustering, merging clusters), and click <span class='mbutton'>GENERATE TILE</span> again, the color in the spatial image will be changed accordingly.

That is all for our functions spatial data. In the end, we would also like to show the results for a 10X Genomics Spatial Transcriptomics dataset [Human Breast Cancer Dataset](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.3.0/Visium_FFPE_Human_Breast_Cancer). You can reproduce the result by going through exactly the same steps metioned above. 
<br>
![10x clustering](images/10x-clustering.png)
<br>
<br>
![10x spatial image](images/10x-spatial-image.png)
<br>


