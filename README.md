# project-1-dachsund
project-1-dachsund created by GitHub Classroom

# Microarray Based Tumor Classification Project 1
Microarray technologies allow scientists to quickly and accurately profile the expression of tens of thousands of genes at a time at a fraction of the cost of modern day sequencing protocols. Repositories and databases such as the Gene Expression Omnibus (GEO) are home to tens of thousands of published and unpublished microarray experiments, and millions of individual samples. Microarray gene expression data are useful in a broad range of research tasks, including disease diagnosis, drug discovery, and toxicological research. Several companies have been born out of the application of this technology to disease prediction such as the Oncotype Dx colon, breast, and prostate cancer prediction assays or Agendia’s MammaPrint breast cancer test. With such a tremendous amount of publicly available data and the translational impact of modern microarray analyses, it is essential that every bioinformatics researcher possess the knowledge and skill set to understand, analyze, and interpret these data. This project will give you first-hand experience in acquiring and analyzing a public microarray dataset.

This analysis will focus only on reproducing the results from the comparison of the C3 and C4 tumor subtypes. The study was conducted in a two-phase design, where an initial set of “discovery” samples was used to identify patterns among the samples, and a separate set of “validation” samples was used to test if the results from the discovery set were robust. For this analysis, we have combined the discover and validation set samples into a single dataset that you will use. There are 135 samples in total. (see class website for more details on project, https://bf528.readthedocs.io/en/latest/).


## Contributors

Allison Nau  
Mae Rose  
Sheila Yee   
Abhishek Thakar 

## Repository Contents

Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

**data_normalization_QC.R** -- R script to pre-process and normalize all Affymetrix chip data from Marisa et al. through the Robust Multi-array Average (RMA) method. Quality controls such as relative log expression (RLE) and normalized unscaled standard errors (NUSE) were performed. Batch effects were corrected using the ComBat algorithm. Resulting data was plotted on principal component analysis (PCA).

**20210217_proj1_biologist.Rmd** -- R markdown script to assign gene names to the probes used in the dataset, identify which genes were most differentially expressed, and to identify which GO, KEGG, and Hallmark gene sets were most enriched in the differentially expressed genes. Can be executed from within Rstudio, does not require SCC. Requires packages "hgu133plus2.db", GSEABase, tidyverse, dplyr, knitr, and kableExtra. Requires input file "project_1_step_5_6_2_fromanalyst_20210222v3_mod.csv". Also requires .gmt files for KEGG, GO, and Hallmark genesets ("c2.cp.kegg.v7.2.symbols_kegg_sym.gmt", "c5.go.v7.2.symbols_go_sym.gmt", and "h.all.v7.2.symbols_hallmark_sym.gmt"). 

