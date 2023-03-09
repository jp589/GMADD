# Gut microbiome alpha diversity decreases in relation to body weight, antibiotic exposure, and infection with multidrug-resistant organisms
### Jonathan Panzer
### 3/9/23

## Introduction {#Introduction}

This GMADD github repository contains all code used to process and analyze data and metadata associated with the brief report published in the journal *Open Forum Infectious Diseases* (DOI:___). The following abstract summarizes the results of the brief report:

Abstract:
16S rRNA gene sequencing of 45 fecal samples revealed that subjects with multiple multidrug-resistant organisms (MDROs), subjects weighing greater than 80kg infected with MDRO E. coli, and subjects weighing less than 80kg with exposure to vancomycin and carbapenem antibiotics during hospitalization had significantly decreased gut microbiome richness.  

`.fastq.gz` files associated with this study have been uploaded to the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) under the BioProject number PRJNA942348 and can be accessed with accession numbers SRR23752615-SRR23752663.

The recommended way to replicate the code is to download `Gut_microbiota_alpha_diversity_decrease_analyses.rmd` and "knit" it in R Studio after [installing](#Installation) all required packages. While knitting, the `.rmd` will create a new `GMADD` directory in the current working directory. Then it will download and unzip the `GMADD` repository into this new directory. From there all analyses will be performed and figures/tables in their raw form will be generated.

For those who wish to explore the R environment after analyses have been run, `GMADD.rda`, which is a snapshot of all R objects and functions after analyses were run, can be loaded into R studio.

Please see the `.Rmd` for the complete analysis code.
