---
title: 'ANALYSIS'

category: User Interface
layout: null
type: ui
---
These panels consist of tools for differentially
expressed genes analysis, enrichment analysis and visualizations that
can be used to get a better understanding of the data.

<br>
![Analysis](images/analysis.png)
<br>


* **DE Analysis**
    * Here you can perform differentially expressed gene analysis.
    First, specify the subset that you wish to find DE genes for.
    A subset can be any cluster or any user-defined subset of cells. To
    assign a group of cells into a subset, see
    <span class='mbox'>Store Subset</span> under the
    <span class='mbox'>ANNOTATIONs</span> menu. If
    the second subset is set to
    <span class='keyword'>rest</span>,
    then Cellar will find DE genes for Subset 1 against all the remaining
    cells. Otherwise, Cellar finds DE genes of Subset 1 vs Subset 2.
    After specifying the subsets, click <span class="mbox">FIND DE GENES</span>.
    A table will be generated for showing the DE genes and relevant information.
<br>
![DE](images/de.png)
<br>

* **Feature Visualization**
    * Three types of feature visualization are implemented: **PLOT EXPRESSION**,
    **HEATMAP**, and **VIOLIN PLOT**. To perform the visualization, first select 
    genes you wish to visualize in the first dropdown menu, then select the visualization 
    method in the second dropdown menu **PLOTTING**. 

    * When multiple genes are selected for plotting the expression value, the color will be used for showing the co-expression value. If we let the expression matrix be
    $$M\in\mathbb{R}^{m\small\times n}$$, then given a cell
    $$i\in[m]$$, we compute its co-expression value
    for genes $$p\in[n]$$ and $$q\in[n]$$ as
    $$c_i = \min \Big\{ \frac{M_{i, p} - v_p}{V_p - v_p},
    \frac{M_{i, q} - v_q}{V_q - v_q} \Big\}$$
    where $$v_t = \min_j M_{j, t}$$ and $$V_t = \max_j M_{j, t}$$, i.e., the
    minimum and maximum expression values for gene $$t$$. Notice that
    $$0\leq c_i\leq 1$$. We do this normalization so as not to bias the
    co-expression towards any one gene. We choose to display the minimum
    between the two values in the expression above as the minimum satisfies
    the nice property that $$\min\{0, p\}=\min\{p, 0\}=\min\{0, 0\}=0$$
    for $$p\geq 0$$. In all these three cases we would want the co-expression
    to be $$0$$.

    Plot Gene Expression:
    <br>
    ![Plot Expression](images/expression.png)
    <br>

    Heatmap:
    <br>
    ![Heatmap](images/heatmap.png)
    <br>

    Violin Plot:
    <br>
    ![Violin Plot](images/violin.png)
    <br>

* **Enrichment Anlysis**
    * After DE analysis, this can be used for performing different enrichment analysis via a hypergeometric test. To do this analysis, simply choose a marker
    list in the first dropdown menu, enter the number of DE genes to use for the analysis, and click <span class="mbox">RUN</span>. Then a table will be generated for showing the results. 
    The entries are sorted by p-value.
    An <span class="mbox">EXPORT</span> button above the table
    will download the dataframe as a <span class="extension">csv</span> file.<br>

<br>
![Enrichment](images/enrichment.png)
<br>