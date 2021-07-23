---
title: 'PREPROCESSING'

category: User Interface
layout: null
type: ui
---
The preprocessing panel has two tabs, the first tab includes preprocessing tools that can be applied to any data type, and the other contains a special scATAC-seq preprocessing pipeline.<br>
**Note**: To save the user some time, all Server Datasets have
been preprocessed with the default options. Therefore, there is no need to
run preprocessing again and you can directly run dimensionality
reduction/clustering analysis.


* **PREPROCESSING**
    * Normally, before running any clustering analysis, the data undergoes
    quality control and basic filtering steps. After that the data is
    normalized and **log1p** transformed. We use the excellent preprocessing
    tools implemented in the [scanpy](https://scanpy.readthedocs.io/en/stable/)
    Python package. Our choice of preprocessing options is guided by this
    [tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html).
    For a full description of the implemented preprocessing options see
    [here](https://scanpy.readthedocs.io/en/stable/api/index.html#basic-preprocessing).

    * The checkboxes can be used for selecting the filters to be applied to the dataset.

    <br>
    ![Preprocessing1](images/preprocessing1.png)
    <br>


* **SCATAC-SEQ**
    * To allow the analysis of scATAC-seq data, we provide an option to convert
    a cell $$\small\times$$ topic matrix to a cell $$\small\times$$ gene one.
    This allows scATAC-seq data to be analyzed via the existing tools
    build in Cellar, just like you normally would with scRNA-seq data.
    The conversion is based on the open chromatin accessibility associated
    with the nearby region of all genes (as listed in
    [GENCODE v35](https://www.gencodegenes.org/human/release_35.html)).

    * Formally, denote by $$P\in\mathbb{R}^{c\times p}$$ the
    cell $$\small\times$$ peak matrix and by $$G\in\mathbb{R}^{c\times g}$$ the
    cell $$\small\times$$ gene matrix to be created. Let the total number of
    protein-coding genes in GENCODE v35 be $$n$$. Fix $$j\in[n]$$.
    <br><br>
    The first step is to extend the length of gene $$j$$. Let the
    $$j^{\text{th}}$$ gene location be the interval
    $$I_j=[\text{start}, \text{end}]$$ where $$\text{start}$$ and $$\text{end}$$
    are taken from GENCODE v35. We extend this interval by multiples of the
    gene length $$(\text{end} - \text{start})$$ or by a fixed number of base
    pairs. This value can be specified by the user in the
    <span class="mbox">Extend</span> box. Any value that ends in `x` will
    be automatically considered a multiple of the gene length, otherwise
    it will be considered a fixed number of base pairs. Additionally, a
    <span class="mbox">Max Extend</span> input box is provided which will limit
    the value of <span class="mbox">Extend</span>. This is only useful when
    using a combination of gene length and a fixed number of base pairs.
    <br><br>
    Assume for convenience that <span class="mbox">Extend $$=5x$$</span> and
    <span class="mbox">Max Extend $$=5000$$</span>. Then
    $$\text{start} \leftarrow \text{start} - \min(5\cdot(\text{end} -
    \text{start}), 5000)$$. The same value is added to $$\text{end}$$.
    If <span class="mbox">Consider strand orientation</span> is set to `Yes`,
    you can specify different extension values to the end of the interval that
    corresponds to the downstream direction of the gene. This direction is
    determined by the `strand` key in the GENCODE v35 file.
    <br><br>
    Denote the updated interval by $$I_j^u$$. The next step is to find peaks
    (bins) whose range intersects with $$I_j^u$$. The values of these peaks
    are added or averaged (depending on the value of
    <span class="mbox">Operation to apply to bins</span>) Formally, if we let
    $$K$$ denote peak ranges, then<br><br>
    $$G_{i, j} = \text{op}(\{P_{i, k}: \text{for all } k \text{ such that} K_k
    \cap I_j^u \neq \emptyset\})$$.
    The newly obtained cell $$\small\times$$ gene matrix undergoes basic
    preprocessing and filtering just like in the Preprocess tab.

    <br>
    ![Preprocessing2](images/preprocessing2.png)
    <br>

* <span class="warn">!!<span class="tooltip">Attention</span></span>
    Running any preprocessing
    will replace the raw data matrix with the preprocessed one.
    If you wish to run preprocessing with different options, you need
    to reload the dataset from the Dataset menu.
    <span class="pn">PN<span class="tooltip">Note to Programmers</span></span>
    In case you export the session as an <span class="extension">h5ad</span>
    file from the Import/Export menu, you can still find the raw data under
    `adata.raw` where `adata` is the AnnData object.