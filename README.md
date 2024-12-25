# Gene-Expression-App

## Overview
App that allows user to explore mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue from Huntington's Disease and neurologically normal individuals. Data can be found on the GEO database [(GSE64810)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810). 

### Functionalities
- Sample Information Exploration
  - Summary of metadata columns
  - Table of metadata with sorting and filtering functionality
  - Violin plot of user-chosen continuous column grouped by diagnosis
- Counts Matrix Exploration
  - For this tab, the user is able to set filters for gene percentile of variance and non-zero samples.
  - Based on this filtering, the user is able to explore:
    - Summary of the effect of user-defined filtering on genes
    - Scatterplot visualizing the effect of user-defined filtering
    - Clustered heatmap visualizing counts after filtering
  - User is also able to choose PCs for X and Y axes and plot a PCA plot 
- Differential Expression Exploration
  - User is able to sort and filter through a table of Deseq2 results
  - User is able to plot a volcano plot of selected X and Y axes and define padj threshold for coloring
- Individual Gene Expression Exploration
  - User is able to select a categorical variable and a gene and choose a type of visualization (bar plot, boxplot, violin plot, or beeswarm) to plot. 

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
