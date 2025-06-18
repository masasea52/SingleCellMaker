# SingleCellMaker: <br>Interactive Single-cell RNA-seq Analysis Web Application<br>

## Overview

'[SingleCellMaker](https://many-happy-cat.shinyapps.io/SingleCellMaker/)' is an interactive web application developed to facilitate exploratory data analysis of single-cell RNA sequencing (scRNA-seq) data. 
It leverages leading R packages such as Seurat and decoupleR as its backend, providing a user-friendly graphical interface (GUI). This enables researchers, 
even those without prior R coding experience, to intuitively execute standard analysis workflows including quality control (QC), 
clustering, dimensionality reduction (UMAP), marker gene identification, and transcription factor (TF) activity inference. Users can also easily visualize and download their results.

## Key Features

* **Create & Load Demo Data**: <p>* Experience SingleCellMaker: Download Sample Data to dive in! </p>

* **10x Genomics Data Upload**: <p>* Directly upload 10x Genomics data folders in ZIP format (containing barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz) to initiate analysis. </p>

* **Interactive Quality Control (QC)**: <p> Adjust QC filtering parameters based on cell feature counts, total RNA counts, and mitochondrial gene percentages, observing their effects in real-time. </p>
    * **`Min Features (nFeature_RNA)`**: Low-quality cells, often characterized by low unique molecular identifier (UMI) counts or few detected genes, can be indicative of damaged or dying cells. This parameter helps to remove such low-quality cells from the analysis.<br>
    * **`Max Features (nFeature_RNA)`**: Doublet cells (two or more cells captured together as one) often exhibit unusually high unique identifier (UMI) counts and large number of detected genes. This parameter helps to exclude suc doublet cells from the analysis.<br>
    * **`Max Mitochondrial gene %`**:  Cells with a high proportion of mitochondrial gene reads often indicate damaged, stressed, or dying cells, representing poor cellular quality. This parameter helps to exclude such low-quality cells from the analysis.<br><br>
   
* **Seurat-Based Analysis**:
     * Adjust key parameters such as the number of principal component analysis (PCA) dimensions and clustering resolution to identify cell populations.<br>
     * Visualize cell populations using UMAP (Uniform Manifold Approximation and Projection).<br><br>
    
* **Marker Gene Identification**: <p> Identify characteristic marker genes for each cluster. The calculation is performed using the `FindAllMarkers` function with the following parameters: `only.pos = TRUE`, `min.pct = 0.3`, `logfc.threshold = 0.4`, and `test.use = "MAST"`. Results are displayed in a data table and can be downloaded in CSV format.</p>
    * **`only.pos = TRUE`**: This parameter ensures that only **positive marker genes** are identified for each cluster. A positive marker gene is defined as a gene that is expressed at a higher level in a specific cluster compared to all other cells (or clusters). This setting is commonly used to find genes that uniquely characterize a cell type. <br>
    * **`min.pct = 0.3`**: This sets a **minimum percentage threshold for detection**. A gene must be detected (i.e., expressed in at least one cell) in at least 30% of the cells in either the cluster of interest OR in all other cells, to be considered for differential expression analysis. This helps to filter out genes that are only sporadically expressed.<br>
    * **`logfc.threshold = 0.4`**: This establishes a **minimum log-fold change (logFC) threshold**. Only genes that show at least a 0.4-fold difference (on a log2 scale) in expression between the cluster of interest and all other cells are considered significant marker genes. A logFC of 0.4 (log2(1.32)) means the gene is expressed approximately 1.32 times higher. This parameter helps to identify genes with substantial expression differences.
    * **`test.use = "MAST"`**: This specifies the **statistical test to be used** for differential expression analysis. MAST (Model-based Analysis of Single-cell Transcriptomics) is a popular and robust statistical framework specifically designed for single-cell RNA-seq data. It accounts for the zero-inflated and multimodal nature of scRNA-seq data, making it suitable for identifying true differentially expressed genes.
    Other common tests include `Wilcox` (Wilcoxon Rank Sum test) or `bimod` (Likelihood-ratio test for zero-inflated data), but MAST is often recommended for its performance with single-cell data.<br><br>

* **Gene/TF Expression Visualization**:<p></p>
    * **`FeaturePlot`** to visualize the expression patterns of selected genes or TFs on the UMAP.<br>
    * **`DotPlot`** and **`VlnPlot`** to compare expression levels across clusters.<br>
    * **`RNA and TF Activity Switching`**: Flexibly switch between "RNA" (gene expression) and "tfsulm" (transcription factor activity) assays to display plots.<br><br>
    
* **Heatmap Generation**: <p>
    * Generate heatmaps of identified marker genes or inferred transcription factor activities to visually capture differences between cell populations. </p>

* **Results Download**: <p>
   * Download all generated plots (UMAP, FeaturePlot, DotPlot, VlnPlot, Heatmap) and marker gene lists in PNG/CSV formats.</p>

* **Estimation of Transcription Factor Activity**:<p></p>
    * **`Species`**: Please select species (Human or mouse) of your data sets.
    * **`Top N most active TFs for Heatmap`**: Please fill in number you want to show transcription factor activity in heatmap.

* **Calculation Pathway Activity**:<p></p>
    * **'Species'**: Please select species (Human or Mouse) of your data sets.
    * **'Run Pathway Activity Analysis'**: Click the "Run Pathway Activity Analysis" button to initiate the estimation of major signaling pathway activities based on gene expression.
    * **'Analysis Progress'**: During Analysis, a progress message, "Estimating pathway activities with Progeny..." will be displayed.
    * **'Visualize Results'**: Upon completion, the pathway activity results will be displayed as a heatmap or table.
    * **'Download Results'**: You can download tha analysis results both as a PNG-formatted heatmap and as a CSV file table.
