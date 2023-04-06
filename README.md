# Dvoretskova_et_al_2023
This repository includes the necessary code to reproduce the analysis presented in the paper titled 'Spatial enhancer activation determines inhibitory neuron identity.'

## Data and code associated with the study

Spatial enhancer activation determines inhibitory neuron identity

doi: https://doi.org/10.1101/2023.01.30.525356

For feedback or suggestions, please contact <mailto:christian.mayer@bi.mpg.de>.

## Authors

Elena Dvoretskova
May C. Ho
Volker Kittke
Ilaria Vitali
Daniel D. Lam
Irene Delgado
Chao Feng
Miguel Torres
Juliane Winkelmann
Christian Mayer

Correspondence: <mailto:christian.mayer@bi.mpg.de>

## Sripts
There are individual scripts for the different parts of the analysis. 

### tCropSeq

This directory contains scripts used to process and analyze tCROP-seq data.

#### 1. Clustering (**clustering.Rmd**)
    - reads in output from the 10X cellranger pipeline, removes cells that fail to pass quallity check, normalizes data, and clusters cells using the Louvain algorithm
#### 2. Integration (**integration.R**)
    - integrates data based on the pipeline described in **Stuart, Butler et al, 2019**. Detailed description of the protocol can be found [here](https://satijalab.org/seurat/articles/integration_introduction.html).
#### 3. Differential Expression analysis (**DE_analysis.Rmd**)
    - reads in seurat objects and performs pseudobullk differential expression testing using [Libra](https://github.com/neurorestore/Libra).
#### 4. Proportion Change (**proportion_change.Rmd**)
    - outputs proportion chanage analysis (effect size dot plot and bar plot in Figure 1.
#### 5. Lineage comparison (**lineage_analysis.Rmd**)
    - outputs UpSet plots shown in Figure 2. 
#### 6. Constructing modules (**module.ipynb**)
    - uses the Hotspot package to identify modules. Scripts follows the steps outlined at their [Github repository](https://github.com/neurorestore/Libra). 
#### 7. Modules (**modules.Rmd**)
    - Calculates the effect size of perturbation on module score and outputs heatmap of top genes in modules. 


### ChipSeq_Meis2

#### preprocessing_peak_calling_MEIS2_mouse_GE.sh

    - MEIS1/2 ChIP-seq analysis starting from .fastq files to peak calling

#### peak overlaps.R

    - Overlap of ChIP-seq peak regions for MEIS1/2, DLX5 and LHX6 datasets.

    - Overlap of peaks with promoters from EPD, enhancers from Gorkin et al. and with Vista enhancers

plot Fig S3e

#### peak-TSS_distance.R

    - Calculation of distance from MEIS1/2 peaks to the nearest TSS; plot for Fig. 3A

#### plot_enhancer-promoter_vista.R 

    - Upset plots for overlap of MEIS1/2 and DLX5 peaks with promoters from EPD, enhancers from Gorkin et al. and with Vista enhancers

    - Fig. 3e, S3d

#### HOMER_motif_analysis_preprocessing.R

    - Extraction of peak sequences (coordinates) of MEIS1/2 binding sites

#### HOMER_motif_anaylsis1.sh

    - Motif discovery and enrichment in MEIS1/2 peaks

#### HOMER_motif_anaylsis2.R
    - Plots for motif density and motif enrichment, Figs. 3c&d, S3b

#### SpaMo_input_sequences.R

    - Extraction of peak sequences for motif spacing analysis of MEIS1/2-DLX5 overlapping peaks

#### SpaMo_motif_analysis.txt

    - Motif spacing analysis parameters, 

    - Fig S3c
