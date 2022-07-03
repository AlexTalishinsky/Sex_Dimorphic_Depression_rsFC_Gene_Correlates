---
output:
  html_document: default
  pdf_document: default
---
# Sexually_Dimorphic_Depression_rsFC_Gene_Correlates
 Scripts and instructions for uncovering brain gene expression correlates of sexually dimorphic resting state functional connectivity changes associated with depression. 

# README

This software package contains three functions:

* ‘map_gene2GlasserROIs.m’ 
    * Maps Allen Human Brain Atlas normalized microarray data to ROIs in the 360-region. ‘Glasser’ functional brain parcellation. 
*	‘run_PLS.m’ 
    * Uses ‘map_gene2GlasserROIs.m’ to generate a normalized brain regional gene expression matrix and use it in a partial least squares (PLS) regression model to predict a neuroimaging-derived response variable. 
* ‘run_ENGLM.m’ 
    * Uses the MATLAB functions ‘lassoglm’ and ‘perfcurve’ to generate elastic-net regularized general linear models predicting diagnostic status using rsFC predictor variables. 

The functions ‘map_gene2GlasserROIs.m’ and ‘run_PLS.m’ are intended to reproduce the gene-neuroimaging PLS results in Figures 4-5, Supplementary Figures 5-7, and Supplementary Tables 2-4 of the submitted manuscript. The ‘run_ENGLM.m’ function is intended to reproduce the results in Figure 6c-e of the submitted manuscript. 


# 1. System requirements

### 1.1 Hardware Dependencies:

This package requires a standard computer with enough RAM to support in-memory operations.

### 1.2 Software Dependencies: 

This package has been tested on Mac OS Catalina 10.15.6, Mac OS Mojave 10.14.6, and Windows 10. 

The custom code in this package is run in MATLAB, which must be downloaded and installed from <https://www.mathworks.com/downloads/>. Version 2018a or above is recommended. We have tested our software package on versions 2018a, 2019a, and 2020a 

In addition to MATLAB software, the ‘map_gene2GlasserROIs.m’ and ‘run_PLS.m’ functions in this package require the following software/files to be downloaded, installed, and added to the MATLAB search path: 

* CONN toolbox
    * available for download at <https://web.conn-toolbox.org/>
* Midnight Scan Club Codebase
    * available for download at <https://github.com/MidnightScanClub/MSCcodebase>
* ‘spls' package by Joao M. Monteiro
    * available for download at <https://github.com/jmmonteiro/spls>
* ‘rotate-parcellation' package by Frantisek Vasa
    * available for download at <https://github.com/frantisekvasa/rotate_parcellation>
* ‘Richiardi_Data_File_S2.csv’ file 
    * available for download in the supplement of the article, <Richiardi et al, Correlated gene expression supports synchronous activity in brain networks', Science, 2015. DOI: 10.1126/science.1255905> 
* ‘Glasser_SubCort.dtseries.nii’ 
    * included in this folder
* ‘lh.midthickness.32k_fs_LR.surf.gii’ 
    * included in this folder
* ‘rh.midthickness.32k_fs_LR.surf.gii’ 
    * included in this folder
 

# 2. Installation guide

After installing the software dependencies listed above, the three custom functions in this package may be downloaded and run in MATLAB. Download time should be <1 minute on a normal computer. 


# 3. Demo

### 3.1 Instructions to run on data: 

Instructions for each function are included in the function script header. The ‘map_gene2GlasserROIs.m’ and ‘run_PLS.m’ functions require Allen Human Brain Atlas normalized microarray expression data to be downloaded from <http://human.brain-map.org/static/download> and added to the MATLAB search path. This data can be downloaded directly using the following links: 

* http://human.brain-map.org/api/v2/well_known_file_download/178238387
* http://human.brain-map.org/api/v2/well_known_file_download/178238373
* http://human.brain-map.org/api/v2/well_known_file_download/178238359
* http://human.brain-map.org/api/v2/well_known_file_download/178238316
* http://human.brain-map.org/api/v2/well_known_file_download/178238266
* http://human.brain-map.org/api/v2/well_known_file_download/178236545

Downloading these files may take up to 30 minutes. 

Once downloaded, these files, and all the software dependencies listed above, should be added to the MATLAB search path. Once they are added, the ‘map_gene2GlasserROIs.m’ function can run with the file paths to any of these folders as input, and the ‘run_PLS.m’ function can run with only the attached ‘rsFC_vector’ and ‘sphereHCP’ files as input. 

Expected outputs are described in the function script headers. If run correctly using the provided ‘rsFC_vector’, the output of the ‘run_PLS.m’ function should be identical to the results shown in Supplementary Figure 5 for the ‘Female BA25’ PLS regression model. 

Expected run time is 1-3 minutes for the ‘map_gene2GlasserROIs.m’ function, and 30-40 minutes for the ‘run_PLS.m’  function, depending on the available RAM. Expected run time for the ‘run_ENGLM.m’ function varies greatly depending on the size of inputs. For a 533 x 71631 predictor matrix and 100 bootstrap iterations, 2-4 hours run time is expected depending on the available RAM. 


# 4. Instructions for use

### 4.1 How to run the software on your data: 

See above. Instructions for each function are provided in the function script header. 

### 4.2 Reproduction instructions: 

Quantitative results for the ‘female BA25’ PLS regression model in Supplementary Figure 5 and Supplementary Tables 2-4 may be reproduced by running the ‘run_PLS.m’ function with only the attached ‘rsFC_vector’ and ‘sphereHCP’ files as input, so long as all the function dependencies described above are added to the MATLAB search path. 

Quantitative gene set enrichment analysis results for the ‘female BA25’ PLS regression model in Figures 4 and 5 and Supplementary Figure 6 may be reproduced by inputting the ‘gene_LW_ranklist’ output of the ‘run_PLS.m’ function (run as described above, with attached ‘rsFC_vector’ and ‘sphereHCP’ files as input) into the ‘fGSEA’ function of the R ‘fGSEA’ package as the ‘stats’ argument, and using each of the gene sets in Supplementary Table 5 as input for the ‘pathways’ argument. 
The ‘fGSEA package can be downloaded by running the following code in R:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fgsea")
```

Quantitative Gene Ontology (GO) analysis results for the ‘female BA25’ PLS regression model in Figure 5 and Supplementary Figure 7 may be reproduced by inputting the HGNC gene symbols associated with the Entrez-gene IDs in the first column of ‘gene_LW_ranklist’ output of the ‘run_PLS.m’ function (run as described above, with attached ‘rsFC_vector’ and ‘sphereHCP’ files as input) into the online GOrilla analysis tool, which can be found at <http://cbl-gorilla.cs.technion.ac.il/> .
