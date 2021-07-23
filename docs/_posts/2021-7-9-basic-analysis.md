---
title: 'Tutorial: Basic Gene Expression Analysis Pipeline'

category: Tutorials
layout: null
type: null
---
This is a step-by-step tutorial for the basic single cell gene expression data analysis using Cellar. THe functions included in this tutorial work for all types of datasets. In this tutorial, we walk through the pipeline using default settings. However, feel free to check out detailed information about other available options of each step in the corresponding sections, so that you can experiment with different settings and see which one works the best for you. 

**1. Load Dataset**
<br>
    We begin the analysis by loading a dataset. To load a dataset, click the <span class='mbutton'>LOAD DATA</span> button in the navigation bar, which is at the top of the page. This will expand a panel for loading datasets. We can either select a server dataset from the dropdown menu, or upload a dataset and select it from the menu. See the **LOAD DATA** section for specific instructions regarding uploading a dataset. Here we select take server dataset `HBMP1_lymph_node_1` as an example. After select a dataset, we click the <span class='mbutton'>LOAD</span> to load the dataset.

<br>
![Basic Load Dataset](images/basic-load-dataset.png)
<br>
    After loading the dataset, information stored in the dataset (e.g. 2D plot, celltype annotations, etc.) will be visualized. However, this dataset only has the original expression matrix. So only the number of cells and number of genes are shown. 

<br>
![Basic Cells Genes](images/basic-cells-genes.png)
<br>

**2. Preprocessing (optional)**
<br>
    If your dataset is not preprocessed, which is not the case for server datasets, you can use click the  <span class='mbutton'>PREPROCESSING</span> button in the navigation bar and use our preprocessing tools. See the **PREPROCESSING** section for more details.

<br>
![Preprocessing](images/preprocessing1.png)
<br>

**3. Dimensionality Reduction**
<br>
    Go to the **DIMENSIONALITY REDUCTION** panel in the sidebar and run dimensionality reduction by: selecting embedding method (e.g. PCA) -> selecting 2D embedding method (e.g. UMAP, TSNE) -> click the <span class='mbutton'>RUN</span> button. You can also configure the parameters by clicking the setting buttons on the right of the dropdown menus. Here, "embedding" is a low dimensional (default 40) embedding used for subsequent analysis (e.g. clustering and DE analysis), and 2D embedding is for the visalization in the main plot. After this step, a scatter plot (2D embedding) will be shown in the main plot where each point is a cell. 

<br>
![Dim Settings](images/dim-settings.png)
<br>

<br>
![2D Embedding](images/2d.png)
<br>

**4. Clustering**
<br>
    Go to the **CLUSTERING** panel in the sidebar and stay in the **UNSUPERVISED** tab. Choose a clustering algorithm in the dropdown menu. The default option is Leiden. Click the button on the right to expand a setting panel and if you want to change the parameters of the algorithm. Finally, click the <span class='mbutton'>PREPROCESSING</span> button to start the clustering. When it's finished, the colors of the points will be changed to indicate clusters.

<br>
![Settings](images/cluster-settings.png)
<br>

<br>
![Leiden](images/leiden.png)
<br>

**5. Annotations**
<br>
    Then in the **ANNOTATIONS** panel, you can: 1. annotate a cluster with a cell type; 2. store a group of cells as a subset for subsequent analysis (DE analysis); 3. merge a subset into a new cluster.
<br>
<br>
    To annotate a cluster, select a cluster using the dropdown menu in the first row, enter the cell type name (e.g. "cell type 1"), and click <span class='mbutton'>STORE</span>.
<br>
![Cell type annotation](images/basic-cell-type.png)
<br>
<br>
    To store a subset, we use second row of this panel. You need to use the lasso tool to select a group of cells in the main plot, enter the new subset name, and click <span class='mbutton'>STORE</span>.
<br>
![Store Subset](images/basic-subset.png)
<br>
<br>
    Then, we can use the last row to merge the cells in this new subset into a new cluster. To do so, simply select this new subset, and click <span class='mbutton'>MERGE</span>.
<br>
![Merge](images/basic-merge.png)
<br>

**6. Differentially Expressed Genes Analysis**
<br>
To do DE analysis, scroll down a bit and go to the **DE Analysis** panel in the **ANALYSIS** tab under the main plot. Select two subsets that you want to find DE genes for. These subsets include the subsets the you defined in previous steps, and all the clusters. Then, after clicking <span class='mbutton'>FIND DE GENES</span>, a table will be shown to present the DE genes. For example, here we find the DE genes of the new subset that we just defined against all the rest cells.

<br>
![DE](images/basic-de.png)
<br>

**7. Feature Visualizations**
<br>
You can also visualize one or more genes in the **Feature Visualization** panel, which is on the right of the **DE Analysis** panel. For example, if we are interested in the first DE gene "CD74", we select this gene by either typing or finding it in the first dropdown. Then we select the type of visualization in the second dropdown, which has 3 options `PLOT EXPRESSION`, `HEATMAP`, and `VIOLIN PLOT`. For example let's first take a look at the violin plot. 
<br>

![violin](images/basic-violin.png)
<br>

We can see this DE gene is indeed expressed higher in the new cluster that we just defined (Cluster 9). 
<br>

If we select `PLOT EXPRESSION`, then the color of the main plot will be changed to indicate expression value. Brighter means higher expression. In this case this DE gene is indeed expressed higher in the new subset.

<br>

![expression](images/basic-plot-expression.png)
<br>
Finally, we plot a heatmap for top 3 DE genes:

<br>
![exprheatmapession](images/basic-heatmap.png)
<br>

**8. Enrichment Analysis**
<br>
In the **Enrichment Analysis** panel at the bottom, we provide several marker lists for enrichment analysis. You can use them to, for example, identify cell types. You only need to choose a marker list, enter the number of genes to use, and click <span class='mbutton'>RUN</span>. Then the results will be presented in a table. Here we found top cell type candidates for the cluster that we ran DE genes analysis for.
<br>
![enrichment](images/basic-enrichment.png)
<br>

**9. Export Session**
<br>
When you are done with the analysis, you can export your results by using the **SESSION** panel in the side bar. **EXPORT SESSION** will download the dataset with all the results as an  <span class='extension'>h5ad</span> file. **EXPORT ANNOTATIONS ONLY** will download the cell annotations as a  <span class='extension'>csv</span> file. 
<br>
![session](images/session.png)
<br>

