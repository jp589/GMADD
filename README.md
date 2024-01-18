# Gut microbiome alpha diversity decreases in relation to body weight, antibiotic exposure, and infection with multidrug-resistant organisms

### Jonathan Panzer

### 1/24

## Introduction

This GMADD github repository contains all code used to process and analyze data and metadata associated with the major article published in *American Journal of Infection Control* (<https://doi.org/10.1016/j.ajic.2023.12.017>). The following abstract summarizes the results:

### Abstract:

#### Background

The human gastrointestinal tract is home to a dense and diverse microbiome, predominated by bacteria. Despite the conservation of critical functionality across most individuals, the composition of the gut microbiome is highly individualized, leading to differential responses to perturbations such as oral antibiotics or multidrug-resistant organism (MDRO) infection. Herein, subject responses to these perturbations based on their body weight were evaluated.

#### Methods

Fecal samples were collected from 45 subjects at the Detroit Medical Center to evaluate the effects of perturbations on subjects' gut microbiome composition. Bacterial profiling was completed using 16S rRNA gene sequencing.

#### Results

Subjects with multiple MDROs, subjects weighing greater than 80 kg infected with MDRO *E coli*, and subjects weighing less than 80 kg with exposure to vancomycin and carbapenem antibiotics during hospitalization had significantly decreased gut microbiome richness.

#### Conclusions

Both administration of oral antibiotics and MDRO infections decreased gut microbiome alpha diversity, but the magnitude of these gut microbiome perturbations was body weight dependent.

### Data Availability and Replication

`.fastq.gz` files associated with this study have been uploaded to the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) under the BioProject number PRJNA942348 and can be accessed with accession numbers SRR23752615-SRR23752663.

The recommended way to replicate the code is to download `Gut_microbiota_alpha_diversity_decrease_analyses.rmd` and "knit" it in R Studio after [installing](#Installation) all required packages. While knitting, the `.rmd` will create a new `GMADD` directory in the current working directory. Then it will download and unzip the `GMADD` repository into this new directory. From there all analyses will be performed and figures/tables in their raw form will be generated.

For those who wish to explore the R environment after analyses have been run, `GMADD.rda`, which is a snapshot of all R objects and functions after analyses were run, can be loaded into R studio.

Please see the `.Rmd` for the complete analysis code.
