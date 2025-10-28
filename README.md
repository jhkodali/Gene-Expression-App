# Gene-Expression-App

## Overview
App that allows user to explore mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue from Huntington's Disease and neurologically normal individuals. Data can be found on the GEO database [(GSE64810)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810). 

[Go to App](https://gene-expression-explorer.shinyapps.io/gene_exp_explorer/)

### Tour of Functionalities
- Sample Information Exploration
  - Summary of metadata columns
  - Table of metadata with sorting and filtering functionality
  - Violin plot of user-chosen continuous column grouped by diagnosis
    <img width="1909" height="707" alt="Screenshot 2025-09-18 214627" src="https://github.com/user-attachments/assets/c05ac92c-ecda-4972-972b-9c4b7de71b56" />
    <img width="1915" height="924" alt="Screenshot 2025-09-18 214717" src="https://github.com/user-attachments/assets/a6279c33-a685-43a9-8fe2-458ccab8e444" />
- Counts Matrix Exploration
  - For this tab, the user is able to set filters for gene percentile of variance and non-zero samples.
  - Based on this filtering, the user is able to explore:
    - Summary of the effect of user-defined filtering on genes
    - Scatterplot visualizing the effect of user-defined filtering
    - Clustered heatmap visualizing counts after filtering
  - User is also able to choose PCs for X and Y axes and plot a PCA plot
    <img width="1919" height="927" alt="Screenshot 2025-09-18 214751" src="https://github.com/user-attachments/assets/a30ac716-d0a0-4174-a6da-e55c0729451a" />
    <img width="1916" height="1001" alt="Screenshot 2025-09-18 220500" src="https://github.com/user-attachments/assets/91b1760f-10be-4e5d-86c8-d07ea9c81436" />
    <img width="1915" height="950" alt="Screenshot 2025-09-18 220527" src="https://github.com/user-attachments/assets/a3fd9d41-7b00-4ed0-b0f5-0ef38332a0c9" />
- Differential Expression Exploration
  - User is able to sort and filter through a table of Deseq2 results
  - User is able to plot a volcano plot of selected X and Y axes and define padj threshold for coloring
    <img width="1918" height="994" alt="Screenshot 2025-09-18 220607" src="https://github.com/user-attachments/assets/4add327a-23ad-4bbb-84ec-454466bfe4da" />
- Individual Gene Expression Exploration
  - User is able to select a categorical variable and a gene and choose a type of visualization (bar plot, boxplot, violin plot, or beeswarm) to plot.
    <img width="1918" height="802" alt="Screenshot 2025-09-18 220702" src="https://github.com/user-attachments/assets/8a821f2d-c569-4dbc-94a5-d8e456fc38ec" />
    <img width="1919" height="785" alt="Screenshot 2025-09-18 220722" src="https://github.com/user-attachments/assets/d2466c15-8f96-44d3-aaf8-67828ebb0a17" />

### Explanation of Data Files
- GSE64810_series_matrix.csv: sample info acquired by GEOquery
- **GSE64810_samples_input.csv**: sample info formatted for app input on Samples tab
- GSE64810_mlhd_DESeq2_norm_counts_adjust.txt: deseq normalized counts matrix downloaded from GEO
- **GSE64810_counts_input.csv**: deseq normalized counts matrix converted to CSV for app input on Counts tab
- GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt: deseq analysis results downloaded from GEO
- **GSE64810_DESeq2_input.csv**: deseq analysis results converted to CSV for app input on Differential Expression tab

*Note: Files in bold are the files that should be uploaded onto the app for its respective tab. The "Individual Gene Expression" tab requires upload of the files from the previous three tabs.*

### Additional R Script
- data_converter.R: converts and formats available data files into input files for upload onto app (files in bold)
