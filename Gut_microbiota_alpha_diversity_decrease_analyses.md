---
title: "Gut microbiota alpha diversity decreases in relation to body weight, antibiotic exposure, and infection with multidrug-resistant organisms"
author: "Jonathan Panzer"
date: "2/16/23"
output:
  html_document:
    keep_md: true
    toc: true
    toc_float: true
---



# Table of Contents

***

1. [Introduction](#Introduction)
2. [DADA2 Processing](#DADA2_processing)
3. [Setup](#Setup)
    * [Installation](#Installation)
    * [Package Loading](#Package_loading)
4. [Function Definitions](#Function_definitions)
5. [Data Processing](#Data_processing)
    * [DECONTAM Processing](#DECONTAM_processing)
    * [Normalization](#Normalization)
    * [MDRO Infection Assessment](#MDRO_setup)
    * [Cohort Body Weight Distribution](#Weight_distribution)
6. [Alpha diversity](#Alpha_diversity)
    * [Exposure to Antibiotics](#Antibiotic_exposure)
    * [Body Weight](#Body_weight)
    * [Number of MDROs](#Number_of_MDROs)
    * [MDRO Status](#MDRO_status)
    * [MDRO *E. coli* infection](#MDRO_E_coli)
    * [Antibiotic Classes](#Antibiotic_classes)
7. [LEfSe](#LEfSe)
8. [Cohort Demographics](#Cohort_demographics)

***

## Introduction {#Introduction}

The GMADD github repository contains all code used to process and analyze data and metadata associated with the short report published at (Journal) (DOI:___).

Abstract:
16S rRNA gene sequencing of 45 fecal samples revealed that subjects with multiple multidrug-resistant organisms, subjects weighing greater than 80kg infected with MDRO E. coli, and subjects weighing less than 80kg with exposure to vancomycin and carbapenem antibiotics during hospitalization had significantly decreased gut microbiome richness. 

`.fastq` files associated with this study have been uploaded to the (Database_name) public database. 

The recommended way to replicate the code is to download the `.rmd` and knit it after installing all required packages.


## DADA2 Processing {#DADA2_processing}

DADA2 processing of fastq files was performed on the HPC grid at Wayne State University similar to the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html). Processing was completed in two steps utilizing two R scripts.


```r
library(dada2)
library(pdftools)
library(RPushbullet)

options(error = recover)

path <- "./"

Study_Name <- "GMADD"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230, 215),
              maxN=0, maxEE=c(2,7), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)

filtered_reads <- out

write.csv(filtered_reads , file = "filtered_reads.csv")

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

save.image(file= paste0("Dada2_Grid_Analysis_pre_ERR", Study_Name, ".rda"))
savehistory(file= paste0("Dada2_Grid_Analysis_pre_ERR", Study_Name, ".rhistory"))

pbPost("note", "Srt Err Lrn", "Now")

errF <- learnErrors(derepFs, multithread=TRUE)

errR <- learnErrors(derepRs, multithread=TRUE)

Errors_plot <- plotErrors(errF, nominalQ=TRUE)

pdf('Errors_plot.pdf')
plot(Errors_plot)
dev.off()

save.image(file= paste0("Dada2_Grid_Analysis_pre_D2", Study_Name, ".rda"))
savehistory(file= paste0("Dada2_Grid_Analysis_pre_D2", Study_Name, ".rhistory"))

pbPost("note", "Srt d2", "Now")

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

save.image(file= paste0("Dada2_Grid_Analysis_post_D2", Study_Name, ".rda"))
savehistory(file= paste0("Dada2_Grid_Analysis_post_D2", Study_Name, ".rhistory"))

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)

pbPost("note", "Input", "Needed")

table <- table(nchar(getSequences(seqtab))) 
write.csv(table, file = "merged_seq_length.csv")

save.image(file= paste0("Dada2_Grid_Analysis_post_merge_", Study_Name, ".rda"))
savehistory(file= paste0("Dada2_Grid_Analysis_post_merge_", Study_Name, ".rhistory"))
```

Based on the distribution of merged sequence lengths, minimum and maximum length cutoffs were chosen.


```r
library(dada2)
library(RPushbullet)

load("Dada2_Grid_Analysis_post_merge_GMADD.rda")

#user input for minimum and maximum sequence lengths based on expected read length
min_len <- 252
max_len <- 253

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% min_len:max_len]

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="pooled", multithread=FALSE, verbose=TRUE)

#transpose seqtab to make easier to copy over
seqtab.nochim.t <- t(seqtab.nochim)
write.csv(seqtab.nochim.t, file = "seqtab.nochim.csv")

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

pbPost("note", "Tax", "Start")
taxa <- assignTaxonomy(seqtab.nochim, "~/D2tax/silva_nr_v138_train_set.fa.gz", minBoot=80, multithread=TRUE)

write.csv(taxa , file = "taxa.csv")

taxa.species <- addSpecies(taxa, "~/D2tax/silva_species_assignment_v138.fa.gz")

write.csv(taxa.species , file = "taxa_with_species.csv")

merged <- cbind(taxa.species, seqtab.nochim.t)
write.csv(merged, file = "merged.csv")

save.image(file= paste0("Dada2_Grid_Analysis_RObjects_", Study_Name, ".rda"))
pbPost("note", "Finished", "Saved")
```

## Setup {#Setup}

This `.Rmd` document was successfully run with RStudio (Version 2022.07.1+554 "Spotted Wakerobin") and R (Version 4.0.3) to produce the figures and tables in Panzer et al. 2023 DOI:___. The version of RStudio can be checked under `Help` -> `About RStudio`. The version of R can be checked and changed under `Tools` -> `Global Options`. In the `General` side tab, there should be an R version listed near the top of the window, which can be changed, but this requires RStudio to restart before the version change takes effect. RStudio version 2022.07.1+554 is stable and can be downloaded [here](https://dailies.rstudio.com/version/2022.07.1+554/). R version 4.0.3 for Windows can be downloaded [here](https://cran.r-project.org/bin/windows/base/old/4.0.3/). Distributions for Linux and macOS can be downloaded [here](https://cran.r-project.org/).

In addition to downloading `RStudio` and `R`, `Rtools` is also necessary since it is required for the installation of the R package `devtools`, which is used here to install older versions of packages and packages residing on github. `Rtools` version 4.0 for windows, which is compatible with `R/4.0.3`, can be downloaded [here](https://cran.r-project.org/bin/windows/Rtools/history.html). Be sure to leave the installation directory as the default `C:/Rtools/`. 

### Installation {#Installation}

The installation chunks are not evaluated when the entire document is knit, but these packages need to be installed before the document will successfully knit. Packages can be installed with several commands. The following chunk lists the code for the specific packages, which were used to successfully knit the document. If a package has never been installed before,  it may be necessary to first install the package with `BiocManager::install("package_name")`, and then use `devtools::install_version("package_name", version = "1.2.3")` for the exact version. If a package has been installed but was built under a different R 4.0 version, it may be necessary to re-install packages using `install.packages("package_name", force = TRUE)`, or to select a specific version with `devtools::install_version("package_name", version = "1.2.3")`. `install.packages("Package_Name")` is also another option.

For exact versions please **skip updates to other packages when asked**.


```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

#allow R to find Rtools via environment variable
Sys.setenv(PATH=paste("C:/Rtools/bin",Sys.getenv("PATH"),sep=";"))
Sys.setenv(BINPREF="C:/Rtools/mingw_$(WIN)/bin/")
#double checking that the path is correct
Sys.which("make")

#devtools required to install package `dada2tools`
if (!requireNamespace("devtools", quietly = TRUE)){
  BiocManager::install("devtools")
}
#installs version 2.1.6 of the usethis package if not currently installed
if(!packageVersion("usethis") == "2.1.6"){
  devtools::install_version("usethis", version = "2.1.6")
}

#installs version 2.4.4 of the devtools package if not currently installed
if(!packageVersion("devtools") == "2.4.4"){
  devtools::install_version("devtools", version = "2.4.4")
}

library(devtools)
```

Note that the R package `reticulate` is only required for saving `plotly` dynamic plots as static images. This process also requires the installation of miniconda. Please refer to [this link](https://search.r-project.org/CRAN/refmans/plotly/html/save_image.html) for installation instructions within R. Otherwise, omit all lines with `save_image()`.


```r
#installs dada2tools version 1.6 from github.
if (!requireNamespace("dada2tools", quietly = TRUE)){
  devtools::install_github("jp589/dada2tools", build_vignettes = TRUE)
}
devtools::install_version("vegan", version = "2.5.6")
devtools::install_version("phyloseq", version = "1.34.0")
devtools::install_version("openxlsx", version = "4.2.4")
devtools::install_version("dplyr", version = "1.0.7")
devtools::install_version("ggplot2", version = "3.3.6")
devtools::install_version("plotly", version = "4.10.0.9001")
devtools::install_version("reticulate", version = "1.24")
devtools::install_version("SummarizedExperiment", version = "1.20.0")
devtools::install_version("finalfit", version = "1.0.6")
devtools::install_version("lefser", version = "1.0.0")
devtools::install_version("ggpubr", version = "0.4.0")
devtools::install_version("curl", version = "4.3.2")
```

### Package Loading {#Package_Loading}


```r
#loading packages
library(dada2tools)
library(phyloseq)
library(openxlsx)
library(dplyr)
library(vegan)
library(ggplot2)
library(plotly)
library(reticulate)
library(SummarizedExperiment)
library(finalfit)
library(lefser)
library(ggpubr)
library(curl)

sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19044)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] curl_4.3.2                  ggpubr_0.4.0               
##  [3] lefser_1.0.0                finalfit_1.0.6             
##  [5] SummarizedExperiment_1.20.0 Biobase_2.50.0             
##  [7] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
##  [9] IRanges_2.24.1              S4Vectors_0.28.1           
## [11] BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
## [13] matrixStats_0.60.1          reticulate_1.24            
## [15] plotly_4.10.0.9001          ggplot2_3.3.6              
## [17] vegan_2.6-2                 lattice_0.20-41            
## [19] permute_0.9-7               dplyr_1.0.7                
## [21] openxlsx_4.2.4              phyloseq_1.34.0            
## [23] dada2tools_1.6             
## 
## loaded via a namespace (and not attached):
##  [1] TH.data_1.1-0          colorspace_2.0-3       ggsignif_0.6.3        
##  [4] ellipsis_0.3.2         modeltools_0.2-23      XVector_0.30.0        
##  [7] rstudioapi_0.14        mice_3.14.0            fansi_1.0.3           
## [10] mvtnorm_1.1-3          coin_1.4-2             codetools_0.2-16      
## [13] splines_4.0.3          cachem_1.0.6           libcoin_1.0-9         
## [16] knitr_1.40             ade4_1.7-19            jsonlite_1.8.0        
## [19] broom_0.7.9            cluster_2.1.0          png_0.1-7             
## [22] compiler_4.0.3         httr_1.4.4             backports_1.2.1       
## [25] assertthat_0.2.1       Matrix_1.5-1           fastmap_1.1.0         
## [28] lazyeval_0.2.2         cli_3.4.0              htmltools_0.5.2       
## [31] tools_4.0.3            igraph_1.3.0           gtable_0.3.1          
## [34] glue_1.6.2             GenomeInfoDbData_1.2.4 reshape2_1.4.4        
## [37] Rcpp_1.0.8.3           carData_3.0-5          jquerylib_0.1.4       
## [40] vctrs_0.4.1            Biostrings_2.58.0      rhdf5filters_1.2.1    
## [43] multtest_2.46.0        ape_5.5                nlme_3.1-149          
## [46] iterators_1.0.14       xfun_0.30              stringr_1.4.1         
## [49] lifecycle_1.0.2        rstatix_0.7.0          zlibbioc_1.36.0       
## [52] MASS_7.3-53            zoo_1.8-9              scales_1.2.1          
## [55] biomformat_1.18.0      sandwich_3.0-1         rhdf5_2.34.0          
## [58] yaml_2.3.5             sass_0.4.0             stringi_1.7.6         
## [61] foreach_1.5.2          boot_1.3-25            zip_2.2.0             
## [64] rlang_1.0.5            pkgconfig_2.0.3        bitops_1.0-7          
## [67] evaluate_0.16          purrr_0.3.4            Rhdf5lib_1.12.1       
## [70] htmlwidgets_1.5.4      tidyselect_1.1.2       plyr_1.8.7            
## [73] magrittr_2.0.3         R6_2.5.1               generics_0.1.2        
## [76] multcomp_1.4-18        DelayedArray_0.16.3    DBI_1.1.1             
## [79] pillar_1.7.0           withr_2.5.0            mgcv_1.8-33           
## [82] abind_1.4-5            survival_3.2-7         RCurl_1.98-1.4        
## [85] tibble_3.1.7           car_3.0-12             crayon_1.4.1          
## [88] utf8_1.2.2             rmarkdown_2.10         grid_4.0.3            
## [91] data.table_1.14.2      forcats_0.5.1          digest_0.6.29         
## [94] tidyr_1.1.3            munsell_0.5.0          viridisLite_0.4.1     
## [97] bslib_0.4.0
```

## Function definitions {#Function_definitions}

These first two functions are wrappers for `phyloseq` and `decontam` functions, which can be applied to samples with normal levels of biomass to remove ASVs more prevalent in negative control samples.


```r
#both use isContaminant
#physeq is a phyloseq object
#thresh is a numeric value between 0 and 1 denoting the cutoff between true and contaminant taxa
#study_name is a character vector which will be listed in the title of the generated plots
decontam_histo_prev_plots2 <- function(physeq, thresh = 0.5, study_name){

  p <- prev <- physeq.neg <- physeq.pos <- contaminant <- NULL
  #creates new column of type logical to determine which samples are controls (TRUE) or samples (FALSE) based on data in Type column.
  phyloseq::sample_data(physeq)$is.neg <- phyloseq::sample_data(physeq)$data == "control"

  #runs isNotContaminant() which splits taxa up into contaminants and true taxa.
  #isNotContaminant() is chosen here over isContaminant() since samples are presumed low biomass and the majority of taxa are assumed to be contaminants.
  #Keep in mind that taxa with true are not contaminants.
  contamdf.prev <- decontam::isContaminant(physeq, neg = "is.neg", detailed = TRUE, threshold = thresh)
  #counts how many true taxa and contaminants there are and prints it to the console.
  print("True represents contaminants:")
  print(table(contamdf.prev$contaminant))

  #code for histogram plotting to evaluate for an appropriate threshold
  #Green tones used for the prevalence palette
  prevalencePalette <- c("2" = "#edf8e9", "3-5" = "#bae4b3", "6-10" = "#74c476", "11+" = "#238b45")
  #getPalette local function calls colorRamp Palette on the prevalence palette
  getPalette = grDevices::colorRampPalette(prevalencePalette)
  #determines the number of steps that the histogram should have (how many unique counts were there at each decontam score level)
  steps <- length(table(contamdf.prev$prev))
  #splits the palette to number of steps determined.
  colr <- getPalette(steps)
  #creates the title for histogram
  histo_title <- paste("Decontam Histgram for", study_name)
  #fill has to be in factor form and then use scale fill manual() to select palatte to use.
  decontam_histogram <- ggplot2::ggplot(contamdf.prev, ggplot2::aes(x = p, fill = factor(prev))) + ggplot2::geom_histogram(bins = 100) + ggplot2::labs(x = 'Decontam score', y = 'Number of ASVs') + ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.05)) + ggplot2::scale_fill_manual(values = colr) + ggplot2::theme_dark() + ggplot2::ggtitle(histo_title)

  #code for prevalence plotting
  #converts sample counts into prevalence counts by turning values greater than 0 into 1's.
  physeq.pa <- phyloseq::transform_sample_counts(physeq, function(abund) 1*(abund>0))

  #subsets samples which are controls into phyloseq object physeq.pa.neg
  physeq.pa.neg <- phyloseq::prune_samples(phyloseq::sample_data(physeq.pa)$data == "control", physeq.pa)

  #subsets samples which are actual samples into phyloseq object physeq.pa.pos
  physeq.pa.pos <- phyloseq::prune_samples(phyloseq::sample_data(physeq.pa)$data == "sample", physeq.pa)

  #creates a data frame where taxa prevalence is summed up across control and true samples. Contaminants are indicated in !contamdf.prev$not.contaminant
  df.pa <- data.frame(physeq.pos= phyloseq::taxa_sums(physeq.pa.pos), physeq.neg = phyloseq::taxa_sums(physeq.pa.neg),
                      contaminant=contamdf.prev$contaminant)

  #creates title for prevalence plot
  prev_title <- paste("Contaminant Prevalence at Threshold:", thresh, "for", study_name, "Study")

  #plots prevalence with taxa technical control prevalence on the x-axis and taxa sample prevalence on the y-axis.
  prev_plot <- ggplot2::ggplot(data=df.pa, ggplot2::aes(x=physeq.neg, y=physeq.pos, color=contaminant)) + ggplot2::geom_point() +
    ggplot2::xlab("Prevalence in Technical Controls") + ggplot2::ylab("Prevalence in Samples") + ggplot2::ggtitle(prev_title)

  #returns plot list
  list(Histogram = decontam_histogram, Prevalence_Plot = prev_plot)
}

#df is a dataframe with ASVs on rows and samples on columns
#physeq is a phyloseq object generated from df
#thresh is a numeric value between 0 and 1 denoting the cutoff between true and contaminant taxa
decontaminate2 <- function(df, physeq, thresh){

  #creates new column of type logical to determine which samples are controls (TRUE) or samples (FALSE) based on data in Type column.
  phyloseq::sample_data(physeq)$is.neg <- phyloseq::sample_data(physeq)$data == "control"

  #isNotContaminant() is chosen here over isContaminant() since samples are presumed low biomass and the majority of taxa are assumed to be contaminants.
  #Keep in mind that Trues are not contaminants
  contamdf.prev <- decontam::isContaminant(physeq, neg = "is.neg", detailed = TRUE, threshold = thresh)

  #subsetting contaminants and 'true taxa'
  true_taxa <- subset(df, !contamdf.prev$contaminant)
  contaminants <- subset(df, contamdf.prev$contaminant)

  #returns list of dataframes true taxa and contaminants along with the threshold used.
  list(TrueTaxa = true_taxa, Contaminants = contaminants, Threshold = thresh)
}
```

These are slightly modified versions of functions found in the R package `lefser`. The main difference is that the grouping column is set by `lefse()` and not by `lefse_plot` so that the grouping column in the returned LEfSe dataframe can be descriptive rather than only `0` or `1`. This impacts the legend plotted by `lefse_plot`. Helper functions `.numeric01`, `filterKruskal`, and `createUniqueValues` are unchanged, but made globally available to `lefse` and `lefse_plot`.

For an overview of LEfSe function parameters and usage, please read the overview found [here](https://waldronlab.io/lefser/articles/lefser.html).


```r
lefse_plot <- function (df, colors = c("red", "forestgreen"), trim.names = TRUE) 
{
  df <- lefser:::trunc(df, trim.names)
  plt <- ggplot2::ggplot(df, aes(reorder(Names, scores), scores)) + 
    ylab("LDA SCORE (log 10)") + theme(axis.title.y = element_blank(), 
                                       axis.title.x = element_text(size = 11, face = "bold"), 
                                       axis.text.y = element_text(vjust = 0.7, size = 9, face = "bold"), 
                                       axis.text.x = element_text(vjust = 0.7, size = 9, face = "bold"), 
                                       plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) + 
    geom_bar(stat = "identity", aes(fill = Group)) + scale_fill_manual(values = colors) + 
    coord_flip()
  return(plt)
}

.numeric01 <- function (x) 
{
  x <- as.factor(x)
  uvals <- levels(x)
  ifelse(x == uvals[1L], 0L, 1L)
}

filterKruskal <- function (expr, group, p.value) 
{
  kw.res <- apply(expr, 1L, function(x) {
    kruskal.test(x ~ group)[["p.value"]]
  })
  kw.sub <- kw.res <= p.value
  kw.sub[is.na(kw.sub)] <- FALSE
  expr[kw.sub, ]
}

createUniqueValues <- function (df, group) 
{
  orderedrows <- rownames(df)
  splitdf <- split(df, group)
  maxim <- vapply(table(group), function(x) max(x * 0.5, 4), 
                  numeric(1L))
  for (i in seq_along(splitdf)) {
    sdat <- splitdf[[i]]
    splitdf[[i]][] <- lapply(sdat, function(cols) {
      if (length(unique(cols)) > maxim[i]) 
        cols
      else abs(cols + rnorm(length(cols), mean = 0, sd = max(cols * 
                                                               0.05, 0.01)))
    })
  }
  df <- do.call(rbind, unname(splitdf))
  df[match(orderedrows, rownames(df)), , drop = FALSE]
}

lefse <- function (expr, kruskal.threshold = 0.05, wilcox.threshold = 0.05, 
          lda.threshold = 2, groupCol = "GROUP", blockCol = NULL, 
          assay = 1L, trim.names = FALSE) 
{
  groupf <- colData(expr)[[groupCol]]
  if (is.null(groupf)) 
    stop("A valid group assignment 'groupCol' must be provided")
  groupf <- as.factor(groupf)
  groupsf <- levels(groupf)
  if (length(groupsf) != 2L) 
    stop("Group classification is not dichotomous:\n", "Found (", 
         paste(groupsf, collapse = ", "), ")")
  group <- lefser:::.numeric01(groupf)
  groups <- levels(groupsf)
  expr_data <- assay(expr, i = assay)
  expr_sub <- lefser:::filterKruskal(expr_data, group, kruskal.threshold)
  if (!is.null(blockCol)) {
    block <- as.factor(colData(expr)[[blockCol]])
    expr_sub <- fillPmatZmat(groupf, block, expr_sub, wilcox.threshold)
  }
  expr_sub_t <- t(expr_sub)
  expr_sub_t_df <- as.data.frame(expr_sub_t)
  expr_sub_t_df <- lefser:::createUniqueValues(expr_sub_t_df, groupf)
  expr_sub_t_df <- cbind(expr_sub_t_df, class = group)
  lfk <- nrow(expr_sub_t_df)
  rfk <- floor(lfk * 2/3)
  ncl <- length(groups)
  min_cl <- as.integer(min(table(expr_sub_t_df$class)) * 2/3 * 
                         2/3 * 0.5)
  min_cl <- max(min_cl, 1)
  eff_size_mat <- replicate(30, suppressWarnings(lefser:::ldaFunction(expr_sub_t_df, 
                                                             lfk, rfk, min_cl, ncl, groups)), simplify = TRUE)
  raw_lda_scores <- rowMeans(eff_size_mat)
  processed_scores <- sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 
                                                 10)
  processed_sorted_scores <- sort(processed_scores)
  scores_df <- data.frame(Names = names(processed_sorted_scores), 
                          scores = as.vector(processed_sorted_scores), stringsAsFactors = FALSE)
  scores_df <- lefser:::trunc(scores_df, trim.names)
  threshold_scores <- abs(scores_df$scores) >= lda.threshold
  scores_df <- scores_df[threshold_scores, ]
  scores_df$Group <- dplyr::case_when(
    scores_df$scores < 0 ~ groupsf[1],
    scores_df$scores > 0 ~ groupsf[2]
  )
  return(scores_df)
}
```

`split_violin_plotly()` produces split violin plots for both Shannon and Chao1 alpha diversity metrics using the R package `plotly`. A subplot is returned which displays the Shannon and Chao1 split violin plots together.


```r
#diversity_df data frame with metadata and alpha diversity data
#shared_axis group column to allow for data to be grouped into a single or multiple split violin plots
#split_axis group column in diversity_df on which to split data into positive and negative sides of the violin
#META2 optional metadata specifically listing MDRO information by subject
#flipped is a boolean value which determines if groups are flipped to the positive or negative side of the split violins
split_violin_plotly <- function(diversity_df, shared_axis, split_axis, META2, flipped = FALSE){
  neg_num <- 1
  pos_num <- 2
  if(!missing(flipped)){
    if(flipped==TRUE){
      neg_num <- 2
      pos_num <- 1
    }
  }
  if(!missing(META2)){
    #adds negative MDRO data
    META2_neg <- META2 %>% filter(Patient_number %in% diversity_df$Patient_number) %>% filter(diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num])
    textly_neg <- ~paste("Subject: ", diversity_df$Patient_number[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num]], "<br>MDRO", META2_neg$MRDO1, "<br>MDRO", META2_neg$MDRO2, "<br>MDRO", META2_neg$MDRO3, "<br>MDRO", META2_neg$MDRO4)
    #adds positive MDRO data
    META2_pos <- META2 %>% filter(Patient_number %in% diversity_df$Patient_number) %>% filter(diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num])
    textly_pos <- ~paste("Subject: ", diversity_df$Patient_number[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num]], "<br>MDRO", META2_pos$MDRO1, "<br>MDRO", META2_pos$MDRO2, "<br>MDRO", META2_pos$MDRO3, "<br>MDRO", META2_pos$MDRO4)
  } else {
    textly_neg <- ~paste("Subject:", diversity_df$Patient_number[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num]])
    textly_pos <- ~paste("Subject:", diversity_df$Patient_number[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num]])
  }
  
  A_plotly_Shannon <- diversity_df %>% plot_ly(type = 'violin')
  #adds negative side to Shannon violin plot
  A_plotly_Shannon <- A_plotly_Shannon %>% add_trace(
  x = diversity_df[[shared_axis]][diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num]],
  y = ~Shannon[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num]],
  legendgroup = unique(sort(diversity_df[[split_axis]]))[neg_num],
  scalegroup = unique(sort(diversity_df[[split_axis]]))[neg_num],
  name = unique(sort(diversity_df[[split_axis]]))[neg_num],
  side = "negative",
  box = list(
    visible = TRUE,
    fillcolor = "white",
    line = list(
      color = "black"
    )
  ),
  meanline = list(
    visible = TRUE
  ),
  color = I("#1B9E77"),
  points = "all",
  pointpos = -.5,
  jitter = 0.5,
  hovertext = textly_neg
)
#adds positive side to Shannon violin plot

A_plotly_Shannon <- A_plotly_Shannon %>% add_trace(
  x = diversity_df[[shared_axis]][diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num]],
  y = ~Shannon[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num]],
  legendgroup = unique(sort(diversity_df[[split_axis]]))[pos_num],
  scalegroup = unique(sort(diversity_df[[split_axis]]))[pos_num],
  name = unique(sort(diversity_df[[split_axis]]))[pos_num],
  side = "positive",
  box = list(
    visible = TRUE,
    fillcolor = "white",
    line = list(
      color = "black"
    )
  ),
  meanline = list(
    visible = TRUE
  ),
  color = I("#D95F02"),
  points = "all",
  pointpos = .5,
  jitter = 0.5,
  hovertext = textly_pos
)

A_plotly_Shannon <- A_plotly_Shannon %>%
  layout(
    xaxis = list(
      title = ""
    ),
    yaxis = list(
      title = "Shannon",
      zeroline = F)
);A_plotly_Shannon

#adds Chao split violin plot
A_plotly_Chao <- diversity_df %>% plot_ly(type = 'violin')

A_plotly_Chao <- A_plotly_Chao %>% add_trace(
  x = diversity_df[[shared_axis]][diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num]],
  y = ~Chao[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[neg_num]],
  legendgroup = unique(sort(diversity_df[[split_axis]]))[neg_num],
  name = unique(sort(diversity_df[[split_axis]]))[neg_num],
  side = "negative",
  scalegroup = "Chao",
  box = list(
    visible = TRUE,
    fillcolor = "white",
    line = list(
      color = "black"
    )
  ),
  meanline = list(
    visible = TRUE
  ),
  color = I("#1B9E77"),
  points = "all",
  pointpos = -.5,
  jitter = 0.5,
  showlegend = F,
  hovertext = textly_neg
)

#adds positive side to Chao violin plot
A_plotly_Chao <- A_plotly_Chao %>% add_trace(
  x = diversity_df[[shared_axis]][diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num]],
  y = ~Chao[diversity_df[[split_axis]] == unique(sort(diversity_df[[split_axis]]))[pos_num]],
  legendgroup = unique(sort(diversity_df[[split_axis]]))[pos_num],
  name = unique(sort(diversity_df[[split_axis]]))[pos_num],
  side = "positive",
  scalegroup = "Chao",
  box = list(
    visible = TRUE,
    fillcolor = "white",
    line = list(
      color = "black"
    )
  ),
  meanline = list(
    visible = TRUE
  ),
  color = I("#D95F02"),
  points = "all",
  pointpos = .5,
  jitter = 0.5,
  showlegend = F,
  hovertext = textly_pos
)
A_plotly_Chao <- A_plotly_Chao %>%
  layout(
    xaxis = list(
      title = ""
    ),
    yaxis = list(
      title = "Chao1",
      zeroline = F)
);A_plotly_Chao

A_plotly <- subplot(nrows = 1, A_plotly_Shannon, A_plotly_Chao, shareX = F, shareY = F, titleX = F, titleY = T, margin = 0.05);A_plotly
}
```

## Data Processing {#Data_processing}

First the data and metadata are downloaded from the [github repository](https://github.com/jp589/GMADD/). Then the data are cleaned and sorted for ease of analysis.


```r
#if directory doesn't exist, create new directory and make current working directory
if(!dir.exists("./GMADD_analysis/")){
  dir.create("./GMADD_analysis/")
};setwd(dir = "./GMADD_analysis/")

url <- "https://github.com/jp589/GMADD/zipball/master"
tmp <- tempfile()
curl::curl_download(url, tmp, mode = "wb")
unzip(tmp);unlink(tmp)

#set working directory to location of unzipped files - test if this works
setwd(list.dirs(full.names = TRUE)[2])

#loading ASV and sample data
merged <- read.csv(file = "GMADD_merged.csv", header = TRUE)
#cleans sequence data and renames rows with numbered ASVs
merged_clean <- dada2tools::dadaset_clean(df = merged, read_thresh = 100);rownames(merged_clean) <- paste0("ASV", 1:nrow(merged_clean))
```

```
## Originally there were 2169 taxa.
##  216 were not bacterial or classified at phylum level.
##  0 were mitochondrial. 5 were chloroplast.
##  1948 taxa remain. 0 samples eliminated with less than 100 reads:
```

```r
#splits merged_clean into taxonomy and ASV table objects
TAX <- merged_clean[,1:8];ASV <- merged_clean[,9:ncol(merged_clean)]

#loads metadata
meta <- read.xlsx(xlsxFile = "GMADD_study_metadata.xlsx", sheet = 1, startRow = 1)
#loads metadata key
meta_key <- read.xlsx(xlsxFile = "GMADD_study_metadata.xlsx", sheet = 2, startRow = 1)


#substitutes values based on key
meta$Place_of_admission <- case_when(
    meta$Place_of_admission == "1" ~ "Home",
    meta$Place_of_admission == "2" ~ "Nursing Home",
    meta$Place_of_admission == "3" ~ "Subacute Rehab",
    meta$Place_of_admission == "4" ~ "Other Hospital",
    meta$Place_of_admission == "5" ~ "Homeless",
    meta$Place_of_admission == "NA" ~ "NA"
)
#substitutes values based on key
meta$MDRO_status <- case_when(
  meta$MDRO_organism == "12" ~ "Negative",
  meta$MDRO_organism == "13" ~ "Negative",
  meta$MDRO_organism == "NA" ~ "NA",
  TRUE ~ "Positive"
)
#substitutes values based on key
meta$C_diff <- case_when(
  meta$C_diff == "1" ~ "C. diff positive",
  meta$C_diff == "0" ~ "C. diff negative",
  TRUE ~ "NA"
)
#substitutes values based on key
meta$C_diff_hist <- case_when(
  meta$C_diff_prior_history == "1" ~ "Prior C.diff",
  meta$C_diff_prior_history == "0" ~ "No Prior C.diff",
  TRUE ~ "NA"
)
#substitutes values based on key
meta$COVID19_Pos <- case_when(
  meta$COVID19_Pos == "0" ~ "COVID19 Negative",
  meta$COVID19_Pos == "1" ~ "COVID19 Positive",
  TRUE ~ "NA"
)
#substitutes values based on key
meta$Exposure_to_antibiotics_during_visit <- case_when(
  meta$Exposure_to_antibiotics_during_visit == "0" ~ "No_antibiotics",
  meta$Exposure_to_antibiotics_during_visit == "1" ~ "Given_antibiotics"
)
#creates `Patient_Status` column to summarize whether patients had COVID19, MDRO, or C. difficile infections in combination or in isolation.
meta$Patient_Status <- case_when(
  meta$COVID19_Pos == "COVID19 Positive" & meta$MDRO_status == "Positive" & meta$C_diff == "C. diff positive" ~ "COVID19/MDRO/Cdiff",
  meta$COVID19_Pos == "COVID19 Positive" & meta$MDRO_status == "Positive" & meta$C_diff == "C. diff negative" ~ "COVID19/MDRO",
  meta$COVID19_Pos == "COVID19 Positive" & meta$MDRO_status == "Negative" & meta$C_diff == "C. diff positive" ~ "COVID19/Cdiff",
  meta$COVID19_Pos == "COVID19 Positive" & meta$MDRO_status == "Negative" & meta$C_diff == "C. diff negative" ~ "COVID19",
  meta$COVID19_Pos == "COVID19 Negative" & meta$MDRO_status == "Positive" & meta$C_diff == "C. diff positive" ~ "MDRO/Cdiff",
  meta$COVID19_Pos == "COVID19 Negative" & meta$MDRO_status == "Positive" & meta$C_diff == "C. diff negative" ~ "MDRO",
  meta$COVID19_Pos == "COVID19 Negative" & meta$MDRO_status == "Negative" & meta$C_diff == "C. diff positive" ~ "Cdiff",
  meta$COVID19_Pos == "COVID19 Negative" & meta$MDRO_status == "Negative" & meta$C_diff == "C. diff negative" ~ "Control", TRUE ~ "Blank")
#divides subjects into weight groups with the split at 80kg
meta$Weight_kg <- as.numeric(meta$Weight_kg);split_weight <- 80
meta$Weight_group <- case_when(
  meta$Weight_kg > split_weight ~ "High Weight",
  meta$Weight_kg < split_weight ~ "Low Weight"
)

#substitutes values based on key
meta$Sex <- case_when(grepl("0", meta$Sex) ~ "Male",
                      grepl("1", meta$Sex) ~ "Female")

#move Case_Control to end
meta <- meta[,c(1:63,65:68,64)]
#adds extraction blanks to metadata
meta_blanks <- cbind(c("E1", "E2", "E3", "E4"), matrix(data = rep("NA", 66*4), nrow = 4, ncol = 66), c(rep("Blank",4)))
colnames(meta_blanks) <- colnames(meta); meta_blanks <- as.data.frame(meta_blanks)
#combines sample metadata with extraction blank metadata
META <- rbind(meta, meta_blanks)

#orders samples in ASV table by metadata
ASV <- ASV[,META$Patient_number]

#adding a blank column for shared axis in split violin plots
META$Blank <- ""
#converting columns to numeric
META$Age <- as.numeric(META$Age);META$Weight_kg <- as.numeric(META$Weight_kg);META$BMI <- as.numeric(META$BMI);META$Length_of_stay <- as.numeric(META$Length_of_stay)
```

```
## Warning: NAs introduced by coercion

## Warning: NAs introduced by coercion

## Warning: NAs introduced by coercion

## Warning: NAs introduced by coercion
```

### DECONTAM Processing {#DECONTAM_processing}

Removes ASVs which were more prevalent in extraction blanks than in fecal samples.


```r
META$DECONTAM_type <- case_when(
  META$Case_Control == "Case" ~ "sample",
  META$Case_Control == "Control" ~ "sample",
  META$Case_Control == "Blank" ~ "control"
)

merged_order <- cbind(TAX, ASV)
DECO_prep <- dada2tools::decontam_prep(
  df = merged_order,
  meta = META,
  type = "DECONTAM_type", 
  sample_col = "Patient_number"
)
decontam_plots <- decontam_histo_prev_plots2(physeq = DECO_prep, thresh = 0.5, study_name = "GMADD")
```

```
## [1] "True represents contaminants:"
## 
## FALSE  TRUE 
##  1908    40
```

```r
decontam_plots$Prevalence_Plot
```

![](Gut_microbiota_alpha_diversity_decrease_analyses_files/figure-html/DECONTAM-1.png)<!-- -->

```r
decontaminated <- decontaminate2(df = merged_order, physeq = DECO_prep, thresh = 0.5)
true_taxa <- decontaminated$TrueTaxa
#removal of blanks
true_taxa <- subset(true_taxa, select = -c(E1, E2, E3, E4))
ASV_Deco <- true_taxa[,9:ncol(true_taxa)]
TAX_Deco <- true_taxa[,1:8]
```


```r
decontam_plots$Histogram
```

```
## Warning: Removed 1349 rows containing non-finite values (stat_bin).
```

![](Gut_microbiota_alpha_diversity_decrease_analyses_files/figure-html/DECONTAM_histogram-1.png)<!-- -->

### Normalization {#Normalization}

Allows for diversity comparisons between samples by equalizing sampling depth across all samples.


```r
#normalizing data
phy <- dada2tools::phyloseqize(merged_df = ASV_Deco, tax_df = TAX_Deco)
#sampling depth after normalization
min(sample_sums(phy))
```

```
## [1] 18797
```

```r
phy_rare <- phyloseq::rarefy_even_depth(phy, replace = FALSE, rngseed = 1)
```

```
## `set.seed(1)` was used to initialize repeatable random subsampling.
```

```
## Please record this for your records so others can reproduce.
```

```
## Try `set.seed(1); .Random.seed` for the full vector
```

```
## ...
```

```
## 388OTUs were removed because they are no longer 
## present in any sample after random subsampling
```

```
## ...
```

```r
rarefied <- dada2tools::dephy(phy_rare)
TAX_rare <- rarefied[,1:8]
ASV_rare <- rarefied[,9:ncol(rarefied)]

META <- META %>% filter(Patient_number %in% colnames(ASV_rare))
```

### MDRO infection assessment {#MRDO_setup}

Some subjects had multiple MDRO infections at the same time. Here we break down the MDRO infection data by subject.


```r
#code to get META set up.
META_MDRO <- data.frame(Patient_number = META$Patient_number[1:45], MDRO = META$MDRO_organism[1:45])
META_MDRO <- META_MDRO %>% mutate(MDRO = strsplit(as.character(MDRO), ",")) %>% tidyr::unnest(MDRO)
META_MDRO$MDRO <- as.numeric(META_MDRO$MDRO)
META_MDRO_wide <- data.frame(Patient_number = unique(META_MDRO$Patient_number), MDRO1 = rep(12, 45), MDRO2 = rep(12, 45), MDRO3 = rep(12, 45), MDRO4 = rep(12, 45))

#need to make 4 columns representing four directions of arrows to be plotted with each row as a different subject
for (i in 1:length(unique(META_MDRO$Patient_number))){
  which_ones <- META_MDRO$Patient_number == unique(META_MDRO$Patient_number)[i]
  how_many <- sum(META_MDRO$Patient_number == unique(META_MDRO$Patient_number)[i])
  if (how_many == 1){
    subset <- META_MDRO[which_ones,]
    META_MDRO_wide[i,2] <- subset$MDRO
  }else{
      for(k in 1:how_many){
        subset <- META_MDRO[which_ones,]
        META_MDRO_wide[i,k+1] <- subset$MDRO[k]
      }
  }
}
#recodes 13 to 12 since both indicate that there were no MDROs
META_MDRO_wide$MDRO1 <- case_when(
  META_MDRO_wide$MDRO1 == 13 ~ 12,
  TRUE ~ META_MDRO_wide$MDRO1
)

#takes key value pairs from MDRO organism key and splits by ` = `
MDRO_key <- as.data.frame(t(as.data.frame(strsplit(meta_key$MDRO_organism, split = " = "))));rownames(MDRO_key) <- 1:13
MDRO_key$V1 <- as.numeric(MDRO_key$V1)
#extracts character vector of MDROs in key
MDRO_key_vec <- as.character(MDRO_key$V2)
#names vector by number and changes value for 12 to NA
names(MDRO_key_vec) <- MDRO_key$V1; MDRO_key_vec[12] <- "NA"
#instantiates META_MDRO_wide_recoded
META_MDRO_wide_recoded <- META_MDRO_wide
#Recodes by key vector
META_MDRO_wide_recoded$MDRO1 <- recode(META_MDRO_wide$MDRO1, !!!MDRO_key_vec)
META_MDRO_wide_recoded$MDRO2 <- recode(META_MDRO_wide$MDRO2, !!!MDRO_key_vec)
META_MDRO_wide_recoded$MDRO3 <- recode(META_MDRO_wide$MDRO3, !!!MDRO_key_vec)
META_MDRO_wide_recoded$MDRO4 <- recode(META_MDRO_wide$MDRO4, !!!MDRO_key_vec)

META_MDRO$MDRO_name <- recode(META_MDRO$MDRO, !!!MDRO_key_vec)

#converts `META_MDRO_wide` to number of MDRO per subject
META_MDRO_wide$Num_MDROs <- as.numeric(table(META_MDRO$Patient_number)[unique(META_MDRO$Patient_number)]) - (META_MDRO_wide$MDRO1 == 12)*1

#most prevalent MDROs include Escherichia coli and Pseudomonas aeruginosa
sort(table(META_MDRO$MDRO_name), decreasing = TRUE)
```

```
## 
##                                  NA                    Escherichia coli 
##                                  14                                   8 
##                 No screen for MDROs              Pseudomonas aeruginosa 
##                                   8                                   8 
##                Enterococcus faecium               Klebsiella pneumoniae 
##                                   7                                   6 
##          Staphylococcus epidermidis              Acinetobacter baumanii 
##                                   4                                   3 
##        Enterobacter colacae complex                   Proteus mirabilis 
##                                   2                                   2 
## Klebsiella (Enterobacter) aerogenes       Pseudomonas fluorescens group 
##                                   1                                   1 
##                 Rothia mucilaginosa 
##                                   1
```

```r
#append to META
META <- cbind(META, META_MDRO_wide[,2:6])
```

### Cohort body weight distribution {#Weight_distribution}

Weight distribution across the cohort.


```r
fig <- plot_ly(
    data = META,
    type = "scatter",
    mode = "markers",
    x = ~reorder(Patient_number, Weight_kg),
    y = ~Weight_kg,
    color = ~Sex,
    colors = c("pink", "blue")
) %>% add_lines(y = 80,line = list(
  color = "grey"
)) %>% layout(
    yaxis = list(
        zeroline = F,
        title = "Weight (kg)"
    ),
    xaxis = list(
      title = "Subjects"
    ),
    legend = list(
      title = "Sex"
    )
);fig
```

```{=html}
<div id="htmlwidget-f7eea6fd9a425d09fcf6" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-f7eea6fd9a425d09fcf6">{"x":{"visdat":{"300419e55995":["function () ","plotlyVisDat"]},"cur_data":"300419e55995","attrs":{"300419e55995":{"mode":"markers","x":{},"y":{},"color":{},"colors":["pink","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"},"300419e55995.1":{"mode":"lines","x":{},"y":80,"color":{},"colors":["pink","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","line":{"color":"grey"},"inherit":true}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"yaxis":{"domain":[0,1],"automargin":true,"zeroline":false,"title":"Weight (kg)"},"xaxis":{"domain":[0,1],"automargin":true,"title":"Subjects","type":"category","categoryorder":"array","categoryarray":["F7","F9","C10","C15","F22","C11","F26","F2","F30","C2","F24","F44","C9","F27","F21","F32","F43","C16","F45","F48","F13","F35","C1","F46","F11","F31","C7","F18","F34","F16","F33","F10","F5","F28","F38","F37","F3","F6","F12","F4","F23","C5","F47","F41","C8"]},"legend":{"title":"Sex"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"mode":"markers","x":["C9","C15","C16","F3","F7","F12","F13","F22","F23","F24","F31","F34","F37","F43","F44","F45","F46","F47","F48"],"y":[66,51,70,105,45,112,76.7,53.2,114.7,63,82,87.5,102.9,69.1,63.9,70,78.2,132,71.4],"type":"scatter","name":"Female","marker":{"color":"rgba(255,192,203,1)","line":{"color":"rgba(255,192,203,1)"}},"textfont":{"color":"rgba(255,192,203,1)"},"error_y":{"color":"rgba(255,192,203,1)"},"error_x":{"color":"rgba(255,192,203,1)"},"line":{"color":"rgba(255,192,203,1)"},"xaxis":"x","yaxis":"y","frame":null},{"mode":"markers","x":["C1","C2","C5","C7","C8","C10","C11","F2","F4","F5","F6","F9","F10","F11","F16","F18","F21","F26","F27","F28","F30","F32","F33","F35","F38","F41"],"y":[78,63,124,85,166,49,54,59,114,93.6,107.2,48.3,92.3,81.8,88.6,85.9,68,57.6,66.2,93.8,60.5,68,90,77.5,98.6,139.5],"type":"scatter","name":"Male","marker":{"color":"rgba(0,0,255,1)","line":{"color":"rgba(0,0,255,1)"}},"textfont":{"color":"rgba(0,0,255,1)"},"error_y":{"color":"rgba(0,0,255,1)"},"error_x":{"color":"rgba(0,0,255,1)"},"line":{"color":"rgba(0,0,255,1)"},"xaxis":"x","yaxis":"y","frame":null},{"mode":"lines","x":["F7","C15","F22","F24","F44","C9","F43","C16","F45","F48","F13","F46","F31","F34","F37","F3","F12","F23","F47"],"y":[80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80],"type":"scatter","line":{"color":"grey"},"name":"Female","marker":{"color":"rgba(255,192,203,1)","line":{"color":"rgba(255,192,203,1)"}},"textfont":{"color":"rgba(255,192,203,1)"},"error_y":{"color":"rgba(255,192,203,1)"},"error_x":{"color":"rgba(255,192,203,1)"},"xaxis":"x","yaxis":"y","frame":null},{"mode":"lines","x":["F9","C10","C11","F26","F2","F30","C2","F27","F21","F32","F35","C1","F11","C7","F18","F16","F33","F10","F5","F28","F38","F6","F4","C5","F41","C8"],"y":[80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80],"type":"scatter","line":{"color":"grey"},"name":"Male","marker":{"color":"rgba(0,0,255,1)","line":{"color":"rgba(0,0,255,1)"}},"textfont":{"color":"rgba(0,0,255,1)"},"error_y":{"color":"rgba(0,0,255,1)"},"error_x":{"color":"rgba(0,0,255,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(fig, file = "Fig_S1_weight_distribution_of_cohort_by_sex.png")
```

## Alpha diversity calculations {#Alpha_diversity}

Chao1 and Shannon alpha diversity metrics were calculated with `vegan`.


```r
#alpha diversity metrics
chao <- c()
#for each sample on columns calculate Chao1
for(i in 1:dim(ASV_rare)[2]){
  chao[i] <- vegan::estimateR(ASV_rare[,i])[[2]]
}
#combines metadata with Shannon and Chao1 alpha diversity metric calculations
diversity_Alpha_with_M <- cbind(
  data.frame(
    Shannon = vegan::diversity(x = ASV_rare, index = "shannon", MARGIN = 2), 
    Chao = chao,
    META)
  )
#convert antibiotic columns to numeric
diversity_Alpha_with_M[,55:64] <- lapply(diversity_Alpha_with_M[,55:64], as.numeric)

div_Alpha_Low <- diversity_Alpha_with_M %>% filter(Weight_group == "Low Weight")
div_Alpha_High <- diversity_Alpha_with_M %>% filter(Weight_group == "High Weight")
```

### Exposure to antibiotics {#Antibiotic_exposure}

All subjects exposed to antibiotics during their visit had significantly decreased alpha diversity compared to those not exposed.


```r
p <- split_violin_plotly(diversity_df = diversity_Alpha_with_M, shared_axis = "Blank", split_axis = "Exposure_to_antibiotics_during_visit", flipped = TRUE);p
```

```
## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis
```

```{=html}
<div id="htmlwidget-3392c9696adb2ad02f8c" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-3392c9696adb2ad02f8c">{"x":{"data":[{"fillcolor":"rgba(31,119,180,0.498)","type":"violin","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","",""],"y":[2.60603121385823,3.85919539902935,3.32875845211335,3.87838317971592],"legendgroup":"No_antibiotics","scalegroup":"No_antibiotics","name":"No_antibiotics","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":["Subject: C8","Subject: C11","Subject: F27","Subject: F33"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"y":[2.10921401105009,2.21852492591197,2.80298990780815,2.53480454429879,3.60289719753111,2.85908547271524,2.89381232528356,0.988094338169116,3.12169633121412,0.15068687281452,3.44773511150555,1.32541682553659,2.57070426636493,2.18516084263117,0.231413345233765,3.38233487480608,3.73572332009627,2.73213856677722,2.42801897711492,2.67967017681985,1.59432822416783,0.531679872114147,2.02696133421695,3.47952755027271,3.30223570624823,2.64973451797741,3.02229914155243,2.21378777648883,2.73246914323887,3.12554469046479,3.05706139869675,0.555137945783053,3.70300730614575,2.52617676898279,3.87687505038611,2.01579394773494,0.737327733215049,1.72865040130223,1.43237240902259,2.44790233863895,1.57470353630208],"legendgroup":"Given_antibiotics","scalegroup":"Given_antibiotics","name":"Given_antibiotics","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":["Subject: C1","Subject: C2","Subject: C5","Subject: C7","Subject: C9","Subject: C10","Subject: C15","Subject: C16","Subject: F2","Subject: F3","Subject: F4","Subject: F5","Subject: F6","Subject: F7","Subject: F9","Subject: F10","Subject: F11","Subject: F12","Subject: F13","Subject: F16","Subject: F18","Subject: F21","Subject: F22","Subject: F23","Subject: F24","Subject: F26","Subject: F28","Subject: F30","Subject: F31","Subject: F32","Subject: F34","Subject: F35","Subject: F37","Subject: F38","Subject: F41","Subject: F43","Subject: F44","Subject: F45","Subject: F46","Subject: F47","Subject: F48"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(214,39,40,0.498)","type":"violin","marker":{"color":"rgba(214,39,40,1)","line":{"color":"rgba(214,39,40,1)"}},"line":{"color":"rgba(214,39,40,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","",""],"y":[189,224.052631578947,213,207.1],"legendgroup":"No_antibiotics","name":"No_antibiotics","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject: C8","Subject: C11","Subject: F27","Subject: F33"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"y":[44,48,100,102.142857142857,204.235294117647,116.625,193.368421052632,22,182.076923076923,11,133,13,74.4285714285714,84,10,84,199.111111111111,109.666666666667,38,62,53.2,18,26,234.176470588235,144.428571428571,44,53,63.75,64,104.875,110.142857142857,26,178.5,121,188,106.230769230769,15.3333333333333,47.4285714285714,155.6,34,25],"legendgroup":"Given_antibiotics","name":"Given_antibiotics","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject: C1","Subject: C2","Subject: C5","Subject: C7","Subject: C9","Subject: C10","Subject: C15","Subject: C16","Subject: F2","Subject: F3","Subject: F4","Subject: F5","Subject: F6","Subject: F7","Subject: F9","Subject: F10","Subject: F11","Subject: F12","Subject: F13","Subject: F16","Subject: F18","Subject: F21","Subject: F22","Subject: F23","Subject: F24","Subject: F26","Subject: F28","Subject: F30","Subject: F31","Subject: F32","Subject: F34","Subject: F35","Subject: F37","Subject: F38","Subject: F41","Subject: F43","Subject: F44","Subject: F45","Subject: F46","Subject: F47","Subject: F48"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x2","yaxis":"y2","frame":null}],"layout":{"xaxis":{"domain":[0,0.45],"automargin":true,"anchor":"y"},"xaxis2":{"domain":[0.55,1],"automargin":true,"anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"title":"Chao1","zeroline":false,"anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"title":"Shannon","zeroline":false,"anchor":"x"},"annotations":[],"shapes":[],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"hovermode":"closest","showlegend":true},"attrs":{"300468a4744e":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"300468a4744e.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","",""],"y":{},"legendgroup":"No_antibiotics","scalegroup":"No_antibiotics","name":"No_antibiotics","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":{},"inherit":true},"300468a4744e.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"y":{},"legendgroup":"Given_antibiotics","scalegroup":"Given_antibiotics","name":"Given_antibiotics","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":{},"inherit":true},"300473884424":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"300473884424.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","",""],"y":{},"legendgroup":"No_antibiotics","name":"No_antibiotics","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true},"300473884424.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"y":{},"legendgroup":"Given_antibiotics","name":"Given_antibiotics","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(p=p, file = "Fig_1A_split_violin_exposure_to_antibiotics_alpha.png", width = 8, height = 8)
```

```
## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis
```

```r
anova_results <- aov(Shannon ~ Exposure_to_antibiotics_during_visit, data = diversity_Alpha_with_M); summary(anova_results) #p = 0.0422
```

```
##                                      Df Sum Sq Mean Sq F value Pr(>F)  
## Exposure_to_antibiotics_during_visit  1   4.16   4.161   4.385 0.0422 *
## Residuals                            43  40.80   0.949                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Chao ~ Exposure_to_antibiotics_during_visit, data = diversity_Alpha_with_M); summary(anova_results) #p = 0.000552
```

```
##                                      Df Sum Sq Mean Sq F value   Pr(>F)    
## Exposure_to_antibiotics_during_visit  1  51980   51980   13.93 0.000552 ***
## Residuals                            43 160447    3731                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


### Body weight {#Body_weight}

Alpha diversity after antibiotic administration body weight dependent, but only Shannon.


```r
#significant shannon but not Chao1. Increased alpha diversity in high weight group signifying more but not necessarily different ASVs in high weight group
anova_results <- aov(Shannon ~ Exposure_to_antibiotics_during_visit+Weight_group, data = diversity_Alpha_with_M); summary(anova_results) #p = 0.0351
```

```
##                                      Df Sum Sq Mean Sq F value Pr(>F)  
## Exposure_to_antibiotics_during_visit  1   4.16   4.161   4.766 0.0347 *
## Weight_group                          1   4.14   4.140   4.742 0.0351 *
## Residuals                            42  36.66   0.873                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Chao ~ Exposure_to_antibiotics_during_visit+Weight_group, data = diversity_Alpha_with_M); summary(anova_results) #p = 0.297
```

```
##                                      Df Sum Sq Mean Sq F value   Pr(>F)    
## Exposure_to_antibiotics_during_visit  1  51980   51980  13.968 0.000556 ***
## Weight_group                          1   4149    4149   1.115 0.297042    
## Residuals                            42 156298    3721                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#Stats for figure 1B
anova_results <- aov(Shannon ~ Weight_group, data = diversity_Alpha_with_M %>% filter(Exposure_to_antibiotics_during_visit == "Given_antibiotics")); summary(anova_results) #p = 0.0225
```

```
##              Df Sum Sq Mean Sq F value Pr(>F)  
## Weight_group  1   5.03   5.028    5.65 0.0225 *
## Residuals    39  34.70   0.890                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Chao ~ Weight_group, data = diversity_Alpha_with_M %>% filter(Exposure_to_antibiotics_during_visit == "Given_antibiotics")); summary(anova_results) #p = 0.247
```

```
##              Df Sum Sq Mean Sq F value Pr(>F)
## Weight_group  1   5463    5463    1.38  0.247
## Residuals    39 154340    3957
```

```r
p <- split_violin_plotly(diversity_df = diversity_Alpha_with_M %>% filter(Exposure_to_antibiotics_during_visit == "Given_antibiotics"), shared_axis = "Blank", split_axis = "Weight_group", flipped = TRUE);p
```

```
## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis
```

```{=html}
<div id="htmlwidget-990531a1ae14759a89b4" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-990531a1ae14759a89b4">{"x":{"data":[{"fillcolor":"rgba(31,119,180,0.498)","type":"violin","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","","","","","","","","","","","","","","","","","","","",""],"y":[2.10921401105009,2.21852492591197,3.60289719753111,2.85908547271524,2.89381232528356,0.988094338169116,3.12169633121412,2.18516084263117,0.231413345233765,2.42801897711492,0.531679872114147,2.02696133421695,3.30223570624823,2.64973451797741,2.21378777648883,3.12554469046479,0.555137945783053,2.01579394773494,0.737327733215049,1.72865040130223,1.43237240902259,1.57470353630208],"legendgroup":"Low Weight","scalegroup":"Low Weight","name":"Low Weight","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":["Subject: C1","Subject: C2","Subject: C9","Subject: C10","Subject: C15","Subject: C16","Subject: F2","Subject: F7","Subject: F9","Subject: F13","Subject: F21","Subject: F22","Subject: F24","Subject: F26","Subject: F30","Subject: F32","Subject: F35","Subject: F43","Subject: F44","Subject: F45","Subject: F46","Subject: F48"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","","","","","","","","","","","","","","","","","",""],"y":[2.80298990780815,2.53480454429879,0.15068687281452,3.44773511150555,1.32541682553659,2.57070426636493,3.38233487480608,3.73572332009627,2.73213856677722,2.67967017681985,1.59432822416783,3.47952755027271,3.02229914155243,2.73246914323887,3.05706139869675,3.70300730614575,2.52617676898279,3.87687505038611,2.44790233863895],"legendgroup":"High Weight","scalegroup":"High Weight","name":"High Weight","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":["Subject: C5","Subject: C7","Subject: F3","Subject: F4","Subject: F5","Subject: F6","Subject: F10","Subject: F11","Subject: F12","Subject: F16","Subject: F18","Subject: F23","Subject: F28","Subject: F31","Subject: F34","Subject: F37","Subject: F38","Subject: F41","Subject: F47"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(214,39,40,0.498)","type":"violin","marker":{"color":"rgba(214,39,40,1)","line":{"color":"rgba(214,39,40,1)"}},"line":{"color":"rgba(214,39,40,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","","","","","","","","","","","","","","","","","","","",""],"y":[44,48,204.235294117647,116.625,193.368421052632,22,182.076923076923,84,10,38,18,26,144.428571428571,44,63.75,104.875,26,106.230769230769,15.3333333333333,47.4285714285714,155.6,25],"legendgroup":"Low Weight","name":"Low Weight","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject: C1","Subject: C2","Subject: C9","Subject: C10","Subject: C15","Subject: C16","Subject: F2","Subject: F7","Subject: F9","Subject: F13","Subject: F21","Subject: F22","Subject: F24","Subject: F26","Subject: F30","Subject: F32","Subject: F35","Subject: F43","Subject: F44","Subject: F45","Subject: F46","Subject: F48"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","","","","","","","","","","","","","","","","","",""],"y":[100,102.142857142857,11,133,13,74.4285714285714,84,199.111111111111,109.666666666667,62,53.2,234.176470588235,53,64,110.142857142857,178.5,121,188,34],"legendgroup":"High Weight","name":"High Weight","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject: C5","Subject: C7","Subject: F3","Subject: F4","Subject: F5","Subject: F6","Subject: F10","Subject: F11","Subject: F12","Subject: F16","Subject: F18","Subject: F23","Subject: F28","Subject: F31","Subject: F34","Subject: F37","Subject: F38","Subject: F41","Subject: F47"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x2","yaxis":"y2","frame":null}],"layout":{"xaxis":{"domain":[0,0.45],"automargin":true,"anchor":"y"},"xaxis2":{"domain":[0.55,1],"automargin":true,"anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"title":"Chao1","zeroline":false,"anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"title":"Shannon","zeroline":false,"anchor":"x"},"annotations":[],"shapes":[],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"hovermode":"closest","showlegend":true},"attrs":{"30042379f55":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"30042379f55.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","","","","","","","","","",""],"y":{},"legendgroup":"Low Weight","scalegroup":"Low Weight","name":"Low Weight","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":{},"inherit":true},"30042379f55.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","","","","","","",""],"y":{},"legendgroup":"High Weight","scalegroup":"High Weight","name":"High Weight","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":{},"inherit":true},"30044d2a2a5f":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"30044d2a2a5f.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","","","","","","","","","",""],"y":{},"legendgroup":"Low Weight","name":"Low Weight","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true},"30044d2a2a5f.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","","","","","","",""],"y":{},"legendgroup":"High Weight","name":"High Weight","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(p = p, file = "Fig_1B_alpha_diversity_Weight_group_Shannon_significant_only_exposed.png", width = 8, height=8)
```

```
## Warning: Can't display both discrete & non-discrete data on same axis

## Warning: Can't display both discrete & non-discrete data on same axis
```

### Number of MDROs {#Number_of_MDROs}

There is a significant inverse relationship between the number of MDROs detected per subject and Chao1 richness (This significance is eliminated when accounting for antibiotic exposure, but number of MDROs is closely tied to antibiotic exposure. These results may imply that a polymicrobial MDRO infection either leads to a decrease in unique bacteria or a decrease in unique bacteria makes one more susceptible to polymicrobial MDRO infection.


```r
#number of MDROs significantly decreases estimated richness only when exposure to antibiotics not considered. Likely because the two are correlated
lm_result <- lm(Chao~ Num_MDROs, data = diversity_Alpha_with_M); summary(lm_result) #p = 0.0255 with adjusted R=0.09002
```

```
## 
## Call:
## lm(formula = Chao ~ Num_MDROs, data = diversity_Alpha_with_M)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -95.81 -53.81 -11.58  68.59 135.55 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  117.807     12.665   9.302  7.4e-12 ***
## Num_MDROs    -19.182      8.291  -2.314   0.0255 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 66.28 on 43 degrees of freedom
## Multiple R-squared:  0.1107,	Adjusted R-squared:  0.09002 
## F-statistic: 5.353 on 1 and 43 DF,  p-value: 0.02553
```

```r
lm_result <- lm(Chao~ Exposure_to_antibiotics_during_visit + Num_MDROs, data = diversity_Alpha_with_M); summary(lm_result) #p=0.10426 with R=0.2575
```

```
## 
## Call:
## lm(formula = Chao ~ Exposure_to_antibiotics_during_visit + Num_MDROs, 
##     data = diversity_Alpha_with_M)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -80.334 -51.488  -1.188  18.666 144.688 
## 
## Coefficients:
##                                                    Estimate Std. Error t value
## (Intercept)                                         102.334     12.380   8.266
## Exposure_to_antibiotics_during_visitNo_antibiotics  105.954     32.395   3.271
## Num_MDROs                                           -12.846      7.736  -1.661
##                                                    Pr(>|t|)    
## (Intercept)                                        2.39e-10 ***
## Exposure_to_antibiotics_during_visitNo_antibiotics  0.00215 ** 
## Num_MDROs                                           0.10426    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 59.87 on 42 degrees of freedom
## Multiple R-squared:  0.2912,	Adjusted R-squared:  0.2575 
## F-statistic: 8.629 on 2 and 42 DF,  p-value: 0.0007255
```

```r
#wilcox test between number of MDROs and exposure to antibiotics
wilcox.test(x = diversity_Alpha_with_M$Num_MDROs, y = as.numeric(as.factor(diversity_Alpha_with_M$Exposure_to_antibiotics_during_visit)))# W =733 with p = 0.01081
```

```
## Warning in wilcox.test.default(x = diversity_Alpha_with_M$Num_MDROs, y =
## as.numeric(as.factor(diversity_Alpha_with_M$Exposure_to_antibiotics_during_visit))):
## cannot compute exact p-value with ties
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  diversity_Alpha_with_M$Num_MDROs and as.numeric(as.factor(diversity_Alpha_with_M$Exposure_to_antibiotics_during_visit))
## W = 733, p-value = 0.01081
## alternative hypothesis: true location shift is not equal to 0
```

```r
#only subjects who had MDROs were given antibiotics
ggboxplot(diversity_Alpha_with_M, 
          x = "Exposure_to_antibiotics_during_visit", 
          y = "Num_MDROs", 
          color = "Exposure_to_antibiotics_during_visit", 
          palette = c("#00AFBB", "#E7B800"),
          ylab = "Number of MDROs", 
          xlab = "Exposure to antibiotics during visit")
```

![](Gut_microbiota_alpha_diversity_decrease_analyses_files/figure-html/Alpha_diversity_Number_of_MDROs-1.png)<!-- -->

```r
#Shannon diversity not impacted by number of MDROs
lm_result <- lm(Shannon~ Num_MDROs, data = diversity_Alpha_with_M); summary(lm_result) #p=0.294 with R=0.002923
```

```
## 
## Call:
## lm(formula = Shannon ~ Num_MDROs, data = diversity_Alpha_with_M)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.2072 -0.4636  0.1597  0.7178  1.3056 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   2.5728     0.1929  13.339   <2e-16 ***
## Num_MDROs    -0.1342     0.1263  -1.063    0.294    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.009 on 43 degrees of freedom
## Multiple R-squared:  0.02558,	Adjusted R-squared:  0.002923 
## F-statistic: 1.129 on 1 and 43 DF,  p-value: 0.2939
```

```r
#Shannon diversity as a function of number of detected MDROs
fig <- plot_ly(
  data = diversity_Alpha_with_M,
  type = "scatter",
  mode = "markers",
  x = ~Num_MDROs,
  y = ~Shannon,
  color = "black",colors = "black",
  hovertext = ~paste("Subject:", diversity_Alpha_with_M$Patient_number)
) %>% add_lines(
    x = ~Num_MDROs,
    y = fitted(lm_result),
    line=
      list(color = "black")
    ) %>% layout(
    yaxis = list(
        zeroline = F  
    )
)

lm_result <- lm(Chao~ Num_MDROs, data = diversity_Alpha_with_M); summary(lm_result) #p = 0.0255 with adjusted R=0.09002
```

```
## 
## Call:
## lm(formula = Chao ~ Num_MDROs, data = diversity_Alpha_with_M)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -95.81 -53.81 -11.58  68.59 135.55 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  117.807     12.665   9.302  7.4e-12 ***
## Num_MDROs    -19.182      8.291  -2.314   0.0255 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 66.28 on 43 degrees of freedom
## Multiple R-squared:  0.1107,	Adjusted R-squared:  0.09002 
## F-statistic: 5.353 on 1 and 43 DF,  p-value: 0.02553
```

```r
#Chao1 richness as a function of number of detected MDROs
fig2 <- plot_ly(
  data = diversity_Alpha_with_M,
  type = "scatter",
  mode = "markers",
  x = ~Num_MDROs,
  y = ~Chao,
  color = "black",colors = "black",
  hovertext = ~paste("Subject:", diversity_Alpha_with_M$Patient_number)) %>% add_lines(
    x = ~Num_MDROs,
    y = fitted(lm_result),
    line=
      list(color = "black")
    ) %>% layout(
    yaxis = list(
        zeroline = F  
    )
  )
p <- subplot(fig, fig2, shareY = FALSE, titleX = TRUE, titleY = TRUE, margin = 0.05);p
```

```{=html}
<div id="htmlwidget-d9c392ead024041f7a33" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-d9c392ead024041f7a33">{"x":{"data":[{"mode":"markers","x":[0,0,0,0,0,0,0,0,0,0,1,3,1,3,3,1,1,2,1,4,1,3,4,1,3,1,2,2,0,2,1,0,0,0,0,0,1,0,0,0,1,0,0,1,0],"y":[2.10921401105009,2.21852492591197,2.80298990780815,2.53480454429879,2.60603121385823,3.60289719753111,2.85908547271524,3.85919539902935,2.89381232528356,0.988094338169116,3.12169633121412,0.15068687281452,3.44773511150555,1.32541682553659,2.57070426636493,2.18516084263117,0.231413345233765,3.38233487480608,3.73572332009627,2.73213856677722,2.42801897711492,2.67967017681985,1.59432822416783,0.531679872114147,2.02696133421695,3.47952755027271,3.30223570624823,2.64973451797741,3.32875845211335,3.02229914155243,2.21378777648883,2.73246914323887,3.12554469046479,3.87838317971592,3.05706139869675,0.555137945783053,3.70300730614575,2.52617676898279,3.87687505038611,2.01579394773494,0.737327733215049,1.72865040130223,1.43237240902259,2.44790233863895,1.57470353630208],"hovertext":["Subject: C1","Subject: C2","Subject: C5","Subject: C7","Subject: C8","Subject: C9","Subject: C10","Subject: C11","Subject: C15","Subject: C16","Subject: F2","Subject: F3","Subject: F4","Subject: F5","Subject: F6","Subject: F7","Subject: F9","Subject: F10","Subject: F11","Subject: F12","Subject: F13","Subject: F16","Subject: F18","Subject: F21","Subject: F22","Subject: F23","Subject: F24","Subject: F26","Subject: F27","Subject: F28","Subject: F30","Subject: F31","Subject: F32","Subject: F33","Subject: F34","Subject: F35","Subject: F37","Subject: F38","Subject: F41","Subject: F43","Subject: F44","Subject: F45","Subject: F46","Subject: F47","Subject: F48"],"type":"scatter","name":"black","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"mode":"lines","x":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,4,4],"y":[2.57277988758213,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.43861645339972,2.30445301921731,2.30445301921731,2.30445301921731,2.30445301921731,2.17028958503489,2.17028958503489,2.17028958503489,2.17028958503489,2.17028958503489,2.03612615085247,2.03612615085247],"hovertext":["Subject: C1","Subject: C2","Subject: C5","Subject: C7","Subject: C8","Subject: C9","Subject: C10","Subject: C11","Subject: C15","Subject: C16","Subject: F27","Subject: F31","Subject: F32","Subject: F33","Subject: F34","Subject: F35","Subject: F38","Subject: F41","Subject: F43","Subject: F45","Subject: F46","Subject: F48","Subject: F2","Subject: F4","Subject: F7","Subject: F9","Subject: F11","Subject: F13","Subject: F21","Subject: F23","Subject: F30","Subject: F37","Subject: F44","Subject: F47","Subject: F10","Subject: F24","Subject: F26","Subject: F28","Subject: F3","Subject: F5","Subject: F6","Subject: F16","Subject: F22","Subject: F12","Subject: F18"],"type":"scatter","line":{"color":"black"},"name":"black","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"mode":"markers","x":[0,0,0,0,0,0,0,0,0,0,1,3,1,3,3,1,1,2,1,4,1,3,4,1,3,1,2,2,0,2,1,0,0,0,0,0,1,0,0,0,1,0,0,1,0],"y":[44,48,100,102.142857142857,189,204.235294117647,116.625,224.052631578947,193.368421052632,22,182.076923076923,11,133,13,74.4285714285714,84,10,84,199.111111111111,109.666666666667,38,62,53.2,18,26,234.176470588235,144.428571428571,44,213,53,63.75,64,104.875,207.1,110.142857142857,26,178.5,121,188,106.230769230769,15.3333333333333,47.4285714285714,155.6,34,25],"hovertext":["Subject: C1","Subject: C2","Subject: C5","Subject: C7","Subject: C8","Subject: C9","Subject: C10","Subject: C11","Subject: C15","Subject: C16","Subject: F2","Subject: F3","Subject: F4","Subject: F5","Subject: F6","Subject: F7","Subject: F9","Subject: F10","Subject: F11","Subject: F12","Subject: F13","Subject: F16","Subject: F18","Subject: F21","Subject: F22","Subject: F23","Subject: F24","Subject: F26","Subject: F27","Subject: F28","Subject: F30","Subject: F31","Subject: F32","Subject: F33","Subject: F34","Subject: F35","Subject: F37","Subject: F38","Subject: F41","Subject: F43","Subject: F44","Subject: F45","Subject: F46","Subject: F47","Subject: F48"],"type":"scatter","name":"black","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"mode":"lines","x":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,4,4],"y":[117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,98.6246327476201,98.62463274762,98.62463274762,98.6246327476201,98.6246327476201,98.62463274762,98.6246327476201,98.6246327476201,98.62463274762,98.6246327476201,98.6246327476201,98.6246327476201,79.4423449052243,79.4423449052243,79.4423449052243,79.4423449052243,60.2600570628286,60.2600570628286,60.2600570628286,60.2600570628286,60.2600570628286,41.0777692204328,41.0777692204328],"hovertext":["Subject: C1","Subject: C2","Subject: C5","Subject: C7","Subject: C8","Subject: C9","Subject: C10","Subject: C11","Subject: C15","Subject: C16","Subject: F27","Subject: F31","Subject: F32","Subject: F33","Subject: F34","Subject: F35","Subject: F38","Subject: F41","Subject: F43","Subject: F45","Subject: F46","Subject: F48","Subject: F2","Subject: F4","Subject: F7","Subject: F9","Subject: F11","Subject: F13","Subject: F21","Subject: F23","Subject: F30","Subject: F37","Subject: F44","Subject: F47","Subject: F10","Subject: F24","Subject: F26","Subject: F28","Subject: F3","Subject: F5","Subject: F6","Subject: F16","Subject: F22","Subject: F12","Subject: F18"],"type":"scatter","line":{"color":"black"},"name":"black","marker":{"color":"rgba(0,0,0,1)","line":{"color":"rgba(0,0,0,1)"}},"textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"xaxis":"x2","yaxis":"y2","frame":null}],"layout":{"xaxis":{"domain":[0,0.45],"automargin":true,"title":"Num_MDROs","anchor":"y"},"xaxis2":{"domain":[0.55,1],"automargin":true,"title":"Num_MDROs","anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"zeroline":false,"title":"Chao","anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"zeroline":false,"title":"Shannon","anchor":"x"},"annotations":[],"shapes":[],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"hovermode":"closest","showlegend":true},"attrs":{"300451f554fc":{"mode":"markers","x":{},"y":{},"hovertext":{},"color":"black","colors":"black","alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"},"300451f554fc.1":{"mode":"lines","x":{},"y":[2.57277988758213,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.43861645339972,2.17028958503489,2.43861645339972,2.17028958503489,2.17028958503489,2.43861645339972,2.43861645339972,2.30445301921731,2.43861645339972,2.03612615085247,2.43861645339972,2.17028958503489,2.03612615085247,2.43861645339972,2.17028958503489,2.43861645339972,2.30445301921731,2.30445301921731,2.57277988758214,2.30445301921731,2.43861645339972,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.57277988758214,2.43861645339972,2.57277988758214,2.57277988758214,2.57277988758214,2.43861645339972,2.57277988758214,2.57277988758214,2.43861645339972,2.57277988758214],"hovertext":{},"color":"black","colors":"black","alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","line":{"color":"black"},"inherit":true},"300446f16b19":{"mode":"markers","x":{},"y":{},"hovertext":{},"color":"black","colors":"black","alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"},"300446f16b19.1":{"mode":"lines","x":{},"y":[117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,98.6246327476201,60.2600570628286,98.62463274762,60.2600570628286,60.2600570628286,98.62463274762,98.6246327476201,79.4423449052243,98.6246327476201,41.0777692204328,98.62463274762,60.2600570628286,41.0777692204328,98.6246327476201,60.2600570628286,98.6246327476201,79.4423449052243,79.4423449052243,117.806920590016,79.4423449052243,98.62463274762,117.806920590016,117.806920590016,117.806920590016,117.806920590016,117.806920590016,98.6246327476201,117.806920590016,117.806920590016,117.806920590016,98.6246327476201,117.806920590016,117.806920590016,98.6246327476201,117.806920590016],"hovertext":{},"color":"black","colors":"black","alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter","line":{"color":"black"},"inherit":true}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(p, file = "Fig_2A_number_of_MDROs_alpha.png", width = 8, height = 8)
```

### MDRO status {#MDRO_status}

Taking into account exposure to antibiotics during visit, MDRO status as a whole does not affect alpha diversity. However, MDROs need to be looked at in isolation rather than as a whole.


```r
diversity_Alpha_MDRO <- diversity_Alpha_with_M %>% filter(MDRO_status == "Positive" | MDRO_status == "Negative")

#Shannon not affected by MDRO status regardless of whether exposure to antibiotics is taken into account
anova_results <- aov(Shannon ~ Exposure_to_antibiotics_during_visit + MDRO_status, data = diversity_Alpha_MDRO); summary(anova_results) #p = 0.9137
```

```
##                                      Df Sum Sq Mean Sq F value Pr(>F)  
## Exposure_to_antibiotics_during_visit  1   4.16   4.161   4.284 0.0447 *
## MDRO_status                           1   0.01   0.012   0.012 0.9137  
## Residuals                            42  40.79   0.971                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Shannon ~ MDRO_status, data = diversity_Alpha_MDRO); summary(anova_results)#p = 0.463
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)
## MDRO_status  1   0.57  0.5674    0.55  0.463
## Residuals   43  44.40  1.0325
```

```r
#Chao1 almost significantly affected by MDRO status when exposure to antibiotics not taken into account
anova_results <- aov(Chao ~ Exposure_to_antibiotics_during_visit + MDRO_status, data = diversity_Alpha_MDRO); summary(anova_results) #p = 0.362
```

```
##                                      Df Sum Sq Mean Sq F value   Pr(>F)    
## Exposure_to_antibiotics_during_visit  1  51980   51980  13.882 0.000575 ***
## MDRO_status                           1   3178    3178   0.849 0.362164    
## Residuals                            42 157269    3744                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Chao ~ MDRO_status, data = diversity_Alpha_MDRO); summary(anova_results) #p = 0.0687
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## MDRO_status  1  15935   15935   3.487 0.0687 .
## Residuals   43 196491    4570                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### MDRO *E. coli* infection {#MDRO_E_coli}

One of the most prevalent MDROs in the cohort was *E. coli*. Intriguingly,  MDRO *E.coli* significantly affected Chao1 richness in the high weight group while accounting for exposure to antibiotics and Shannon was almost significant.


```r
diversity_Alpha_MDRO$Ecoli <- diversity_Alpha_MDRO$Patient_number %in% META_MDRO$Patient_number[META_MDRO$MDRO_name == "Escherichia coli"]

p <- split_violin_plotly(diversity_df = diversity_Alpha_MDRO %>% filter(Weight_group == "High Weight")%>% filter(MDRO_status == "Negative" | MDRO_status == "Positive" & Ecoli == TRUE), shared_axis = "Blank", split_axis = "MDRO_status", META2 = META_MDRO_wide_recoded);p
```

```{=html}
<div id="htmlwidget-77c42872179c5d88cf98" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-77c42872179c5d88cf98">{"x":{"data":[{"fillcolor":"rgba(31,119,180,0.498)","type":"violin","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","","","","","",""],"y":[2.80298990780815,2.53480454429879,2.60603121385823,2.73246914323887,3.87838317971592,3.05706139869675,2.52617676898279,3.87687505038611],"legendgroup":"Negative","scalegroup":"Negative","name":"Negative","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":["Subject:  C5 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C7 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C8 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F31 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F33 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F34 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F38 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F41 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","","","",""],"y":[1.32541682553659,2.57070426636493,2.73213856677722,2.67967017681985,1.59432822416783],"legendgroup":"Positive","scalegroup":"Positive","name":"Positive","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":["Subject:  F5 <br>MDRO Enterococcus faecium <br>MDRO Klebsiella pneumoniae <br>MDRO Escherichia coli <br>MDRO NA","Subject:  F6 <br>MDRO Escherichia coli <br>MDRO Klebsiella (Enterobacter) aerogenes <br>MDRO Staphylococcus epidermidis <br>MDRO NA","Subject:  F12 <br>MDRO Escherichia coli <br>MDRO Pseudomonas aeruginosa <br>MDRO Enterococcus faecium <br>MDRO Proteus mirabilis","Subject:  F16 <br>MDRO Escherichia coli <br>MDRO Pseudomonas aeruginosa <br>MDRO Enterococcus faecium <br>MDRO NA","Subject:  F18 <br>MDRO Escherichia coli <br>MDRO Pseudomonas aeruginosa <br>MDRO Acinetobacter baumanii <br>MDRO Enterococcus faecium"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(214,39,40,0.498)","type":"violin","marker":{"color":"rgba(214,39,40,1)","line":{"color":"rgba(214,39,40,1)"}},"line":{"color":"rgba(214,39,40,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","","","","","",""],"y":[100,102.142857142857,189,64,207.1,110.142857142857,121,188],"legendgroup":"Negative","name":"Negative","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject:  C5 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C7 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C8 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F31 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F33 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F34 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F38 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F41 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","","","",""],"y":[13,74.4285714285714,109.666666666667,62,53.2],"legendgroup":"Positive","name":"Positive","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject:  F5 <br>MDRO Enterococcus faecium <br>MDRO Klebsiella pneumoniae <br>MDRO Escherichia coli <br>MDRO NA","Subject:  F6 <br>MDRO Escherichia coli <br>MDRO Klebsiella (Enterobacter) aerogenes <br>MDRO Staphylococcus epidermidis <br>MDRO NA","Subject:  F12 <br>MDRO Escherichia coli <br>MDRO Pseudomonas aeruginosa <br>MDRO Enterococcus faecium <br>MDRO Proteus mirabilis","Subject:  F16 <br>MDRO Escherichia coli <br>MDRO Pseudomonas aeruginosa <br>MDRO Enterococcus faecium <br>MDRO NA","Subject:  F18 <br>MDRO Escherichia coli <br>MDRO Pseudomonas aeruginosa <br>MDRO Acinetobacter baumanii <br>MDRO Enterococcus faecium"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x2","yaxis":"y2","frame":null}],"layout":{"xaxis":{"domain":[0,0.45],"automargin":true,"anchor":"y"},"xaxis2":{"domain":[0.55,1],"automargin":true,"anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"title":"Chao1","zeroline":false,"anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"title":"Shannon","zeroline":false,"anchor":"x"},"annotations":[],"shapes":[],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"hovermode":"closest","showlegend":true},"attrs":{"30042b6b5ba4":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"30042b6b5ba4.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","",""],"y":{},"legendgroup":"Negative","scalegroup":"Negative","name":"Negative","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":{},"inherit":true},"30042b6b5ba4.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","",""],"y":{},"legendgroup":"Positive","scalegroup":"Positive","name":"Positive","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":{},"inherit":true},"300460ff5617":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"300460ff5617.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","",""],"y":{},"legendgroup":"Negative","name":"Negative","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true},"300460ff5617.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","",""],"y":{},"legendgroup":"Positive","name":"Positive","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(p = p, file = "Fig_2B_Ecoli_high_weight_alpha.png", width = 8, height = 8)

anova_results <- aov(Chao ~ Exposure_to_antibiotics_during_visit+Ecoli, data = diversity_Alpha_MDRO%>% filter(Weight_group == "High Weight") %>% filter(MDRO_status == "Negative" | MDRO_status == "Positive" & Ecoli == TRUE)); summary(anova_results) #p = 0.04197
```

```
##                                      Df Sum Sq Mean Sq F value  Pr(>F)   
## Exposure_to_antibiotics_during_visit  1  19506   19506  14.509 0.00343 **
## Ecoli                                 1   7305    7305   5.434 0.04197 * 
## Residuals                            10  13444    1344                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Shannon ~ Exposure_to_antibiotics_during_visit+Ecoli, data = diversity_Alpha_MDRO%>% filter(Weight_group == "High Weight") %>% filter(MDRO_status == "Negative" | MDRO_status == "Positive" & Ecoli == TRUE)); summary(anova_results) #p = 0.0776
```

```
##                                      Df Sum Sq Mean Sq F value Pr(>F)  
## Exposure_to_antibiotics_during_visit  1  0.731  0.7314   1.887 0.1996  
## Ecoli                                 1  1.499  1.4986   3.866 0.0776 .
## Residuals                            10  3.877  0.3877                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

However, significant effects did not appear in the low weight group.


```r
p <- split_violin_plotly(diversity_df = diversity_Alpha_MDRO %>% filter(Weight_group == "Low Weight")%>% filter(MDRO_status == "Negative" | MDRO_status == "Positive" & Ecoli == TRUE), shared_axis = "Blank", split_axis = "MDRO_status", META2 = META_MDRO_wide_recoded);p
```

```{=html}
<div id="htmlwidget-31910e77114e3b6dbc2b" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-31910e77114e3b6dbc2b">{"x":{"data":[{"fillcolor":"rgba(31,119,180,0.498)","type":"violin","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","","","","","","","","","","","",""],"y":[2.10921401105009,2.21852492591197,3.60289719753111,2.85908547271524,3.85919539902935,2.89381232528356,0.988094338169116,3.32875845211335,3.12554469046479,0.555137945783053,2.01579394773494,1.72865040130223,1.43237240902259,1.57470353630208],"legendgroup":"Negative","scalegroup":"Negative","name":"Negative","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":["Subject:  C1 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C2 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C9 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C10 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C11 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C15 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C16 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F27 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F32 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F35 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F43 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F45 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F46 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F48 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","",""],"y":[3.12169633121412,2.02696133421695,0.737327733215049],"legendgroup":"Positive","scalegroup":"Positive","name":"Positive","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":["Subject:  F2 <br>MDRO Escherichia coli <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F22 <br>MDRO Klebsiella pneumoniae <br>MDRO Escherichia coli <br>MDRO Proteus mirabilis <br>MDRO NA","Subject:  F44 <br>MDRO Escherichia coli <br>MDRO NA <br>MDRO NA <br>MDRO NA"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x","yaxis":"y","frame":null},{"fillcolor":"rgba(214,39,40,0.498)","type":"violin","marker":{"color":"rgba(214,39,40,1)","line":{"color":"rgba(214,39,40,1)"}},"line":{"color":"rgba(214,39,40,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(27,158,119,0.5)","type":"violin","x":["","","","","","","","","","","","","",""],"y":[44,48,204.235294117647,116.625,224.052631578947,193.368421052632,22,213,104.875,26,106.230769230769,47.4285714285714,155.6,25],"legendgroup":"Negative","name":"Negative","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject:  C1 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C2 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C9 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C10 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C11 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C15 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  C16 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F27 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F32 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F35 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F43 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F45 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F46 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F48 <br>MDRO  <br>MDRO NA <br>MDRO NA <br>MDRO NA"],"marker":{"color":"rgba(27,158,119,1)","line":{"color":"rgba(27,158,119,1)"}},"line":{"color":"rgba(27,158,119,1)"},"xaxis":"x2","yaxis":"y2","frame":null},{"fillcolor":"rgba(217,95,2,0.5)","type":"violin","x":["","",""],"y":[182.076923076923,26,15.3333333333333],"legendgroup":"Positive","name":"Positive","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":["Subject:  F2 <br>MDRO Escherichia coli <br>MDRO NA <br>MDRO NA <br>MDRO NA","Subject:  F22 <br>MDRO Klebsiella pneumoniae <br>MDRO Escherichia coli <br>MDRO Proteus mirabilis <br>MDRO NA","Subject:  F44 <br>MDRO Escherichia coli <br>MDRO NA <br>MDRO NA <br>MDRO NA"],"marker":{"color":"rgba(217,95,2,1)","line":{"color":"rgba(217,95,2,1)"}},"line":{"color":"rgba(217,95,2,1)"},"xaxis":"x2","yaxis":"y2","frame":null}],"layout":{"xaxis":{"domain":[0,0.45],"automargin":true,"anchor":"y"},"xaxis2":{"domain":[0.55,1],"automargin":true,"anchor":"y2"},"yaxis2":{"domain":[0,1],"automargin":true,"title":"Chao1","zeroline":false,"anchor":"x2"},"yaxis":{"domain":[0,1],"automargin":true,"title":"Shannon","zeroline":false,"anchor":"x"},"annotations":[],"shapes":[],"images":[],"margin":{"b":40,"l":60,"t":25,"r":10},"hovermode":"closest","showlegend":true},"attrs":{"3004131a6753":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"3004131a6753.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","",""],"y":{},"legendgroup":"Negative","scalegroup":"Negative","name":"Negative","side":"negative","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"hovertext":{},"inherit":true},"3004131a6753.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","",""],"y":{},"legendgroup":"Positive","scalegroup":"Positive","name":"Positive","side":"positive","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"hovertext":{},"inherit":true},"300415013d4f":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin"},"300415013d4f.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","","","","","","","","","","","","",""],"y":{},"legendgroup":"Negative","name":"Negative","side":"negative","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#1B9E77"],"points":"all","pointpos":-0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true},"300415013d4f.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"violin","x":["","",""],"y":{},"legendgroup":"Positive","name":"Positive","side":"positive","scalegroup":"Chao","box":{"visible":true,"fillcolor":"white","line":{"color":"black"}},"meanline":{"visible":true},"color":["#D95F02"],"points":"all","pointpos":0.5,"jitter":0.5,"showlegend":false,"hovertext":{},"inherit":true}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"subplot":true,"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(p = p, file = "Fig_S2_Ecoli_low_weight_alpha.png", width = 8, height = 8)

anova_results <- aov(Chao ~ Exposure_to_antibiotics_during_visit+ Ecoli, data = diversity_Alpha_MDRO%>% filter(Weight_group == "Low Weight") %>% filter(MDRO_status == "Negative" | MDRO_status == "Positive" & Ecoli == TRUE)); summary(anova_results)#p = 0.7115
```

```
##                                      Df Sum Sq Mean Sq F value Pr(>F)  
## Exposure_to_antibiotics_during_visit  1  30165   30165   6.463 0.0235 *
## Ecoli                                 1    665     665   0.142 0.7115  
## Residuals                            14  65341    4667                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova_results <- aov(Shannon ~ Exposure_to_antibiotics_during_visit+ Ecoli, data = diversity_Alpha_MDRO%>% filter(Weight_group == "Low Weight") %>% filter(MDRO_status == "Negative" | MDRO_status == "Positive" & Ecoli == TRUE)); summary(anova_results)#p = 0.8310
```

```
##                                      Df Sum Sq Mean Sq F value Pr(>F)  
## Exposure_to_antibiotics_during_visit  1  4.120   4.120   4.803 0.0458 *
## Ecoli                                 1  0.041   0.041   0.047 0.8310  
## Residuals                            14 12.009   0.858                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Antibiotic classes {#Antibiotic_classes}


```r
#Only Carbapenem in high weight group significantly impacts Shannon alpha diversity
LW_abx_model <- glm(Shannon ~ Vancomycin+Cephalosporin+Floroquinolone+Carbapenem+Tetracycline+Penicillin+Aminoglycoside+Monobactam+Macrolide, data = diversity_Alpha_with_M %>% filter(Weight_group == "Low Weight"));summary(LW_abx_model)
```

```
## 
## Call:
## glm(formula = Shannon ~ Vancomycin + Cephalosporin + Floroquinolone + 
##     Carbapenem + Tetracycline + Penicillin + Aminoglycoside + 
##     Monobactam + Macrolide, data = diversity_Alpha_with_M %>% 
##     filter(Weight_group == "Low Weight"))
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -1.49548  -0.55779  -0.00516   0.70550   1.38508  
## 
## Coefficients: (1 not defined because of singularities)
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      3.2211     0.4923   6.543  9.3e-06 ***
## Vancomycin      -1.0033     0.5160  -1.944   0.0708 .  
## Cephalosporin   -0.3150     0.3133  -1.006   0.3306    
## Floroquinolone   1.2179     1.2486   0.975   0.3448    
## Carbapenem      -0.5919     0.3935  -1.504   0.1533    
## Tetracycline    -0.2745     0.8452  -0.325   0.7499    
## Penicillin           NA         NA      NA       NA    
## Aminoglycoside   0.3049     0.6506   0.469   0.6461    
## Monobactam       0.2925     1.2728   0.230   0.8214    
## Macrolide       -0.1742     1.1335  -0.154   0.8799    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 1.139814)
## 
##     Null deviance: 24.285  on 23  degrees of freedom
## Residual deviance: 17.097  on 15  degrees of freedom
## AIC: 79.97
## 
## Number of Fisher Scoring iterations: 2
```

```r
HW_abx_model <- glm(Shannon ~ Vancomycin+ Cephalosporin+Floroquinolone+Carbapenem+Tetracycline+Penicillin+Aminoglycoside+Monobactam+Macrolide, data = diversity_Alpha_with_M %>% filter(Weight_group == "High Weight"));summary(HW_abx_model)
```

```
## 
## Call:
## glm(formula = Shannon ~ Vancomycin + Cephalosporin + Floroquinolone + 
##     Carbapenem + Tetracycline + Penicillin + Aminoglycoside + 
##     Monobactam + Macrolide, data = diversity_Alpha_with_M %>% 
##     filter(Weight_group == "High Weight"))
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7498  -0.3664   0.0000   0.6498   1.1945  
## 
## Coefficients: (2 not defined because of singularities)
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      3.0859     0.4058   7.605 3.87e-06 ***
## Vancomycin       0.4430     0.5250   0.844   0.4141    
## Cephalosporin   -0.1651     0.2172  -0.760   0.4608    
## Floroquinolone       NA         NA      NA       NA    
## Carbapenem      -0.8552     0.3625  -2.359   0.0346 *  
## Tetracycline    -0.6016     1.4017  -0.429   0.6748    
## Penicillin       0.1208     0.6895   0.175   0.8637    
## Aminoglycoside   0.2248     0.3952   0.569   0.5793    
## Monobactam      -0.5219     1.1965  -0.436   0.6699    
## Macrolide            NA         NA      NA       NA    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 0.8273133)
## 
##     Null deviance: 16.367  on 20  degrees of freedom
## Residual deviance: 10.755  on 13  degrees of freedom
## AIC: 63.543
## 
## Number of Fisher Scoring iterations: 2
```

```r
#Vancomycin and Carbapenems significantly decrease Chao1 diversity in low weight group but not high weight group.
LW_abx_model <- glm(Chao ~ Vancomycin+Cephalosporin+Floroquinolone+Carbapenem+Tetracycline+Penicillin+Aminoglycoside+Monobactam+Macrolide, data = diversity_Alpha_with_M %>% filter(Weight_group == "Low Weight"));summary(LW_abx_model)
```

```
## 
## Call:
## glm(formula = Chao ~ Vancomycin + Cephalosporin + Floroquinolone + 
##     Carbapenem + Tetracycline + Penicillin + Aminoglycoside + 
##     Monobactam + Macrolide, data = diversity_Alpha_with_M %>% 
##     filter(Weight_group == "Low Weight"))
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -95.56  -43.13    0.00   42.45   93.34  
## 
## Coefficients: (1 not defined because of singularities)
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     183.306     29.467   6.221 1.64e-05 ***
## Vancomycin      -72.409     30.886  -2.344   0.0332 *  
## Cephalosporin   -34.731     18.751  -1.852   0.0838 .  
## Floroquinolone   79.453     74.742   1.063   0.3046    
## Carbapenem      -57.139     23.553  -2.426   0.0283 *  
## Tetracycline     29.602     50.593   0.585   0.5672    
## Penicillin           NA         NA      NA       NA    
## Aminoglycoside    6.443     38.946   0.165   0.8708    
## Monobactam       36.122     76.188   0.474   0.6422    
## Macrolide       -28.738     67.851  -0.424   0.6779    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 4083.995)
## 
##     Null deviance: 119208  on 23  degrees of freedom
## Residual deviance:  61260  on 15  degrees of freedom
## AIC: 276.38
## 
## Number of Fisher Scoring iterations: 2
```

```r
HW_abx_model <- glm(Chao ~ Vancomycin+ Cephalosporin+Floroquinolone+Carbapenem+Tetracycline+Penicillin+Aminoglycoside+Monobactam+Macrolide, data = diversity_Alpha_with_M %>% filter(Weight_group == "High Weight"));summary(HW_abx_model)
```

```
## 
## Call:
## glm(formula = Chao ~ Vancomycin + Cephalosporin + Floroquinolone + 
##     Carbapenem + Tetracycline + Penicillin + Aminoglycoside + 
##     Monobactam + Macrolide, data = diversity_Alpha_with_M %>% 
##     filter(Weight_group == "High Weight"))
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -113.74   -25.52    13.08    30.66    89.21  
## 
## Coefficients: (2 not defined because of singularities)
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     163.389     27.672   5.904  5.2e-05 ***
## Vancomycin      -31.798     35.806  -0.888   0.3906    
## Cephalosporin     3.648     14.813   0.246   0.8093    
## Floroquinolone       NA         NA      NA       NA    
## Carbapenem      -45.948     24.722  -1.859   0.0859 .  
## Tetracycline     -8.930     95.590  -0.093   0.9270    
## Penicillin      -43.051     47.022  -0.916   0.3766    
## Aminoglycoside   -3.226     26.955  -0.120   0.9066    
## Monobactam      -39.610     81.597  -0.485   0.6354    
## Macrolide            NA         NA      NA       NA    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 3847.678)
## 
##     Null deviance: 88435  on 20  degrees of freedom
## Residual deviance: 50020  on 13  degrees of freedom
## AIC: 240.88
## 
## Number of Fisher Scoring iterations: 2
```

```r
#2d versions
plot_ly(data = diversity_Alpha_with_M, x=~Carbapenem, y=~Chao, color = ~Weight_group, type = "scatter", colors = c("red", "blue"))
```

```
## No scatter mode specifed:
##   Setting the mode to markers
##   Read more about this attribute -> https://plotly.com/r/reference/#scatter-mode
```

```{=html}
<div id="htmlwidget-fcad3680bff6bbcb60b8" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-fcad3680bff6bbcb60b8">{"x":{"visdat":{"300470f63600":["function () ","plotlyVisDat"]},"cur_data":"300470f63600","attrs":{"300470f63600":{"x":{},"y":{},"color":{},"colors":["red","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"automargin":true,"title":"Carbapenem"},"yaxis":{"domain":[0,1],"automargin":true,"title":"Chao"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0,0,0,1,0,1,2,1,0,0,0,2,0,0,0,0,0,1,1,0,1],"y":[100,102.142857142857,189,11,133,13,74.4285714285714,84,199.111111111111,109.666666666667,62,53.2,234.176470588235,53,64,207.1,110.142857142857,178.5,121,188,34],"type":"scatter","mode":"markers","name":"High Weight","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,1,0,0,0,0,0,1,0,0,1,2,2,0,1,0,1,0,1,0,0,0,0,0],"y":[44,48,204.235294117647,116.625,224.052631578947,193.368421052632,22,182.076923076923,84,10,38,18,26,144.428571428571,44,213,63.75,104.875,26,106.230769230769,15.3333333333333,47.4285714285714,155.6,25],"type":"scatter","mode":"markers","name":"Low Weight","marker":{"color":"rgba(0,0,255,1)","line":{"color":"rgba(0,0,255,1)"}},"textfont":{"color":"rgba(0,0,255,1)"},"error_y":{"color":"rgba(0,0,255,1)"},"error_x":{"color":"rgba(0,0,255,1)"},"line":{"color":"rgba(0,0,255,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
plot_ly(data = diversity_Alpha_with_M, x=~Vancomycin, y=~Chao, color = ~Weight_group, type = "scatter", colors = c("red", "blue"))
```

```
## No scatter mode specifed:
##   Setting the mode to markers
##   Read more about this attribute -> https://plotly.com/r/reference/#scatter-mode
```

```{=html}
<div id="htmlwidget-b73b0c39ebba9edbdaeb" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-b73b0c39ebba9edbdaeb">{"x":{"visdat":{"300451131498":["function () ","plotlyVisDat"]},"cur_data":"300451131498","attrs":{"300451131498":{"x":{},"y":{},"color":{},"colors":["red","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"automargin":true,"title":"Vancomycin"},"yaxis":{"domain":[0,1],"automargin":true,"title":"Chao"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0,1,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1],"y":[100,102.142857142857,189,11,133,13,74.4285714285714,84,199.111111111111,109.666666666667,62,53.2,234.176470588235,53,64,207.1,110.142857142857,178.5,121,188,34],"type":"scatter","mode":"markers","name":"High Weight","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1],"y":[44,48,204.235294117647,116.625,224.052631578947,193.368421052632,22,182.076923076923,84,10,38,18,26,144.428571428571,44,213,63.75,104.875,26,106.230769230769,15.3333333333333,47.4285714285714,155.6,25],"type":"scatter","mode":"markers","name":"Low Weight","marker":{"color":"rgba(0,0,255,1)","line":{"color":"rgba(0,0,255,1)"}},"textfont":{"color":"rgba(0,0,255,1)"},"error_y":{"color":"rgba(0,0,255,1)"},"error_x":{"color":"rgba(0,0,255,1)"},"line":{"color":"rgba(0,0,255,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
#only significant predictors
abx_model <- glm(Chao ~ Vancomycin+Carbapenem, data = diversity_Alpha_with_M); summary(abx_model)
```

```
## 
## Call:
## glm(formula = Chao ~ Vancomycin + Carbapenem, data = diversity_Alpha_with_M)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -107.905   -44.607     2.826    47.206   120.893  
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   151.90      16.10   9.434 6.21e-12 ***
## Vancomycin    -52.59      18.62  -2.824  0.00722 ** 
## Carbapenem    -41.71      13.68  -3.049  0.00396 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 3575.627)
## 
##     Null deviance: 212427  on 44  degrees of freedom
## Residual deviance: 150176  on 42  degrees of freedom
## AIC: 500.79
## 
## Number of Fisher Scoring iterations: 2
```

```r
HW_abx_model <- glm(Chao ~ Vancomycin+Carbapenem, data = diversity_Alpha_with_M %>% filter(Weight_group == "High Weight")); summary(HW_abx_model)
```

```
## 
## Call:
## glm(formula = Chao ~ Vancomycin + Carbapenem, data = diversity_Alpha_with_M %>% 
##     filter(Weight_group == "High Weight"))
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -114.913   -45.648     0.495    41.199   100.851  
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   157.91      22.18   7.120 1.24e-06 ***
## Vancomycin    -48.26      28.81  -1.676    0.111    
## Carbapenem    -32.00      20.47  -1.563    0.135    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 3383.23)
## 
##     Null deviance: 88435  on 20  degrees of freedom
## Residual deviance: 60898  on 18  degrees of freedom
## AIC: 235.02
## 
## Number of Fisher Scoring iterations: 2
```

```r
LW_abx_model <- glm(Chao ~ Vancomycin+Carbapenem, data = diversity_Alpha_with_M %>% filter(Weight_group == "Low Weight")) ;summary(LW_abx_model)
```

```
## 
## Call:
## glm(formula = Chao ~ Vancomycin + Carbapenem, data = diversity_Alpha_with_M %>% 
##     filter(Weight_group == "Low Weight"))
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -112.493   -43.235     1.626    41.768   116.251  
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   156.49      24.48   6.392 2.46e-06 ***
## Vancomycin    -68.51      27.00  -2.537   0.0192 *  
## Carbapenem    -57.22      20.42  -2.802   0.0107 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 3726.328)
## 
##     Null deviance: 119208  on 23  degrees of freedom
## Residual deviance:  78253  on 21  degrees of freedom
## AIC: 270.26
## 
## Number of Fisher Scoring iterations: 2
```

```r
#coefficients of the model
HW_cf.mod <- coef(HW_abx_model)
LW_cf.mod <- coef(LW_abx_model)
### Calculate z on a grid of x-y values
xseq <- seq(0,2,length.out=25)
yseq <- seq(0,1,length.out=25)
Hz <- t(outer(xseq, yseq, function(x,y) HW_cf.mod[1]+HW_cf.mod[2]*x+HW_cf.mod[3]*y))
Lz <- t(outer(xseq, yseq, function(x,y) LW_cf.mod[1]+LW_cf.mod[2]*x+LW_cf.mod[3]*y))

Hz_col <- Hz*0
Lz_col <- Lz*0+1

#saved as Chao1 richness Vancomycin Carbapenem planes
p <- plot_ly(colors = c("red","blue")) %>% add_surface(
  x =~xseq,
  y=~yseq, 
  z=~Hz,
  opacity = 0.7,
  surfacecolor =Hz_col,
  cauto = F,
  cmax=1,
  cmin=0,
  showscale = FALSE,
  legendgroup= "High Weight",
  name= "High Weight"
) %>% add_surface(
  x=~xseq, 
  y=~yseq, 
  z=~Lz,
  opacity = 0.7,
  surfacecolor = Lz_col,
  cauto = F,
  cmax=1,
  cmin=0,
  showscale = FALSE,
  legendgroup = "Low Weight",
  name= "Low Weight"
) %>% add_trace(
  data = diversity_Alpha_with_M, 
  x=~Carbapenem, 
  y=~Vancomycin, 
  z=~Chao, 
  color = ~Weight_group, 
  type = "scatter3d", 
  mode = "markers",
  legendgroup = ~Weight_group
  ) %>% layout(
    scene = list(
      xaxis = list(
        title = "Carbapenem",
        tickvals = c(0,1,2)
      ),
      yaxis = list(
        title = "Vancomycin",
        tickvals = c(0,1)
      ),
      zaxis = list(
        title = "Chao1"
      )
    )
  );p
```

```{=html}
<div id="htmlwidget-05f4b6912f53f51ab5fe" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-05f4b6912f53f51ab5fe">{"x":{"visdat":{"3004413c76aa":["function () ","plotlyVisDat"],"3004653ca33":["function () ","data"]},"cur_data":"3004653ca33","attrs":{"3004413c76aa":{"colors":["red","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"z":{},"type":"surface","x":{},"y":{},"opacity":0.7,"surfacecolor":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],"cauto":false,"cmax":1,"cmin":0,"showscale":false,"legendgroup":"High Weight","name":"High Weight","inherit":true},"3004413c76aa.1":{"colors":["red","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"z":{},"type":"surface","x":{},"y":{},"opacity":0.7,"surfacecolor":[[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]],"cauto":false,"cmax":1,"cmin":0,"showscale":false,"legendgroup":"Low Weight","name":"Low Weight","inherit":true},"3004653ca33":{"colors":["red","blue"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":{},"y":{},"z":{},"color":{},"type":"scatter3d","mode":"markers","legendgroup":{},"inherit":true}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"Carbapenem","tickvals":[0,1,2]},"yaxis":{"title":"Vancomycin","tickvals":[0,1]},"zaxis":{"title":"Chao1"}},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"colorbar":{"title":"Hz<br />Lz<br />Weight_group","ticklen":2},"colorscale":[["0","rgba(255,0,0,1)"],["0.0416666666666667","rgba(251,0,23,1)"],["0.0833333333333333","rgba(248,0,37,1)"],["0.125","rgba(244,0,49,1)"],["0.166666666666667","rgba(240,0,60,1)"],["0.208333333333333","rgba(236,0,70,1)"],["0.25","rgba(232,0,80,1)"],["0.291666666666667","rgba(228,0,89,1)"],["0.333333333333333","rgba(223,0,99,1)"],["0.375","rgba(218,0,108,1)"],["0.416666666666667","rgba(213,0,118,1)"],["0.458333333333333","rgba(207,0,127,1)"],["0.5","rgba(202,0,136,1)"],["0.541666666666667","rgba(195,0,146,1)"],["0.583333333333333","rgba(188,0,156,1)"],["0.625","rgba(181,0,165,1)"],["0.666666666666667","rgba(173,0,175,1)"],["0.708333333333333","rgba(164,0,185,1)"],["0.75","rgba(154,0,195,1)"],["0.791666666666667","rgba(143,0,204,1)"],["0.833333333333333","rgba(130,0,214,1)"],["0.875","rgba(115,0,224,1)"],["0.916666666666667","rgba(96,0,235,1)"],["0.958333333333333","rgba(69,0,245,1)"],["1","rgba(0,0,255,1)"]],"showscale":false,"z":[[157.912389970738,153.890365123434,149.86834027613,145.846315428826,141.824290581523,137.802265734219,133.780240886915,129.758216039612,125.736191192308,121.714166345004,117.6921414977,113.670116650397,109.648091803093,105.626066955789,101.604042108486,97.5820172611818,93.559992413878,89.5379675665743,85.5159427192706,81.4939178719669,77.4718930246632,73.4498681773594,69.4278433300557,65.405818482752,61.3837936354483],[156.579092133412,152.557067286108,148.535042438804,144.513017591501,140.490992744197,136.468967896893,132.44694304959,128.424918202286,124.402893354982,120.380868507678,116.358843660375,112.336818813071,108.314793965767,104.292769118464,100.27074427116,96.2487194238561,92.2266945765524,88.2046697292486,84.1826448819449,80.1606200346412,76.1385951873375,72.1165703400337,68.09454549273,64.0725206454263,60.0504957981226],[155.245794296086,151.223769448783,147.201744601479,143.179719754175,139.157694906871,135.135670059568,131.113645212264,127.09162036496,123.069595517656,119.047570670353,115.025545823049,111.003520975745,106.981496128442,102.959471281138,98.9374464338341,94.9154215865304,90.8933967392267,86.8713718919229,82.8493470446192,78.8273221973155,74.8052973500118,70.7832725027081,66.7612476554043,62.7392228081006,58.7171979607969],[153.912496458761,149.890471611457,145.868446764153,141.846421916849,137.824397069546,133.802372222242,129.780347374938,125.758322527635,121.736297680331,117.714272833027,113.692247985723,109.67022313842,105.648198291116,101.626173443812,97.6041485965084,93.5821237492047,89.560098901901,85.5380740545973,81.5160492072935,77.4940243599898,73.4719995126861,69.4499746653824,65.4279498180786,61.4059249707749,57.3839001234712],[152.579198621435,148.557173774131,144.535148926827,140.513124079524,136.49109923222,132.469074384916,128.447049537613,124.425024690309,120.402999843005,116.380974995701,112.358950148398,108.336925301094,104.31490045379,100.292875606486,96.2708507591827,92.248825911879,88.2268010645753,84.2047762172716,80.1827513699678,76.1607265226641,72.1387016753604,68.1166768280567,64.094651980753,60.0726271334492,56.0506022861455],[151.245900784109,147.223875936805,143.201851089502,139.179826242198,135.157801394894,131.135776547591,127.113751700287,123.091726852983,119.069702005679,115.047677158376,111.025652311072,107.003627463768,102.981602616464,98.9595777691608,94.9375529218571,90.9155280745533,86.8935032272496,82.8714783799459,78.8494535326422,74.8274286853384,70.8054038380347,66.783378990731,62.7613541434273,58.7393292961235,54.7173044488198],[149.912602946783,145.89057809948,141.868553252176,137.846528404872,133.824503557569,129.802478710265,125.780453862961,121.758429015657,117.736404168354,113.71437932105,109.692354473746,105.670329626443,101.648304779139,97.6262799318351,93.6042550845314,89.5822302372276,85.5602053899239,81.5381805426202,77.5161556953165,73.4941308480127,69.472106000709,65.4500811534053,61.4280563061016,57.4060314587979,53.3840066114941],[148.579305109458,144.557280262154,140.53525541485,136.513230567547,132.491205720243,128.469180872939,124.447156025635,120.425131178332,116.403106331028,112.381081483724,108.359056636421,104.337031789117,100.315006941813,96.2929820945094,92.2709572472057,88.2489323999019,84.2269075525982,80.2048827052945,76.1828578579908,72.1608330106871,68.1388081633833,64.1167833160796,60.0947584687759,56.0727336214722,52.0507087741684],[147.246007272132,143.223982424828,139.201957577525,135.179932730221,131.157907882917,127.135883035614,123.11385818831,119.091833341006,115.069808493702,111.047783646399,107.025758799095,103.003733951791,98.9817091044874,94.9596842571837,90.93765940988,86.9156345625763,82.8936097152725,78.8715848679688,74.8495600206651,70.8275351733614,66.8055103260577,62.7834854787539,58.7614606314502,54.7394357841465,50.7174109368428],[145.912709434806,141.890684587503,137.868659740199,133.846634892895,129.824610045592,125.802585198288,121.780560350984,117.75853550368,113.736510656377,109.714485809073,105.692460961769,101.670436114465,97.6484112671617,93.626386419858,89.6043615725543,85.5823367252506,81.5603118779469,77.5382870306431,73.5162621833394,69.4942373360357,65.472212488732,61.4501876414282,57.4281627941245,53.4061379468208,49.3841130995171],[144.579411597481,140.557386750177,136.535361902873,132.51333705557,128.491312208266,124.469287360962,120.447262513658,116.425237666355,112.403212819051,108.381187971747,104.359163124444,100.33713827714,96.3151134298361,92.2930885825323,88.2710637352286,84.2490388879249,80.2270140406212,76.2049891933174,72.1829643460137,68.16093949871,64.1389146514063,60.1168898041026,56.0948649567988,52.0728401094951,48.0508152621914],[143.246113760155,139.224088912851,135.202064065548,131.180039218244,127.15801437094,123.135989523636,119.113964676333,115.091939829029,111.069914981725,107.047890134422,103.025865287118,99.0038404398141,94.9818155925104,90.9597907452066,86.9377658979029,82.9157410505992,78.8937162032955,74.8716913559918,70.849666508688,66.8276416613843,62.8056168140806,58.7835919667769,54.7615671194731,50.7395422721694,46.7175174248657],[141.912815922829,137.890791075526,133.868766228222,129.846741380918,125.824716533614,121.802691686311,117.780666839007,113.758641991703,109.7366171444,105.714592297096,101.692567449792,97.6705426024884,93.6485177551847,89.626492907881,85.6044680605772,81.5824432132735,77.5604183659698,73.5383935186661,69.5163686713623,65.4943438240586,61.4723189767549,57.4502941294512,53.4282692821475,49.4062444348437,45.38421958754],[140.579518085504,136.5574932382,132.535468390896,128.513443543593,124.491418696289,120.469393848985,116.447369001681,112.425344154378,108.403319307074,104.38129445977,100.359269612466,96.3372447651627,92.315219917859,88.2931950705553,84.2711702232516,80.2491453759478,76.2271205286441,72.2050956813404,68.1830708340367,64.1610459867329,60.1390211394292,56.1169962921255,52.0949714448218,48.072946597518,44.0509217502143],[139.246220248178,135.224195400874,131.202170553571,127.180145706267,123.158120858963,119.136096011659,115.114071164356,111.092046317052,107.070021469748,103.047996622444,99.0259717751408,95.003946927837,90.9819220805333,86.9598972332296,82.9378723859259,78.9158475386221,74.8938226913184,70.8717978440147,66.849772996711,62.8277481494072,58.8057233021035,54.7836984547998,50.7616736074961,46.7396487601924,42.7176239128886],[137.912922410852,133.890897563549,129.868872716245,125.846847868941,121.824823021637,117.802798174334,113.78077332703,109.758748479726,105.736723632423,101.714698785119,97.6926739378151,93.6706490905113,89.6486242432076,85.6265993959039,81.6045745486002,77.5825497012964,73.5605248539927,69.538500006689,65.5164751593853,61.4944503120816,57.4724254647778,53.4504006174741,49.4283757701704,45.4063509228667,41.3843260755629],[136.579624573527,132.557599726223,128.535574878919,124.513550031615,120.491525184312,116.469500337008,112.447475489704,108.425450642401,104.403425795097,100.381400947793,96.3593761004894,92.3373512531857,88.3153264058819,84.2933015585782,80.2712767112745,76.2492518639708,72.2272270166671,68.2052021693633,64.1831773220596,60.1611524747559,56.1391276274522,52.1171027801484,48.0950779328447,44.073053085541,40.0510282382373],[135.246326736201,131.224301888897,127.202277041593,123.18025219429,119.158227346986,115.136202499682,111.114177652379,107.092152805075,103.070127957771,99.0481031104674,95.0260782631637,91.00405341586,86.9820285685563,82.9600037212525,78.9379788739488,74.9159540266451,70.8939291793414,66.8719043320376,62.8498794847339,58.8278546374302,54.8058297901265,50.7838049428228,46.761780095519,42.7397552482153,38.7177304009116],[133.913028898875,129.891004051572,125.868979204268,121.846954356964,117.82492950966,113.802904662357,109.780879815053,105.758854967749,101.736830120445,97.7148052731417,93.692780425838,89.6707555785343,85.6487307312306,81.6267058839268,77.6046810366231,73.5826561893194,69.5606313420157,65.538606494712,61.5165816474082,57.4945568001045,53.4725319528008,49.4505071054971,45.4284822581933,41.4064574108896,37.3844325635859],[132.57973106155,128.557706214246,124.535681366942,120.513656519638,116.491631672335,112.469606825031,108.447581977727,104.425557130424,100.40353228312,96.3815074358161,92.3594825885123,88.3374577412086,84.3154328939049,80.2934080466012,76.2713831992974,72.2493583519937,68.22733350469,64.2053086573863,60.1832838100825,56.1612589627788,52.1392341154751,48.1172092681714,44.0951844208677,40.0731595735639,36.0511347262602],[131.246433224224,127.22440837692,123.202383529616,119.180358682313,115.158333835009,111.136308987705,107.114284140402,103.092259293098,99.0702344457941,95.0482095984904,91.0261847511866,87.0041599038829,82.9821350565792,78.9601102092755,74.9380853619718,70.916060514668,66.8940356673643,62.8720108200606,58.8499859727569,54.8279611254531,50.8059362781494,46.7839114308457,42.761886583542,38.7398617362382,34.7178368889345],[129.913135386898,125.891110539594,121.869085692291,117.847060844987,113.825035997683,109.80301115038,105.780986303076,101.758961455772,97.7369366084684,93.7149117611647,89.692886913861,85.6708620665572,81.6488372192535,77.6268123719498,73.6047875246461,69.5827626773423,65.5607378300386,61.5387129827349,57.5166881354312,53.4946632881274,49.4726384408237,45.45061359352,41.4285887462163,37.4065638989126,33.3845390516088],[128.579837549572,124.557812702269,120.535787854965,116.513763007661,112.491738160358,108.469713313054,104.44768846575,100.425663618446,96.4036387711427,92.381613923839,88.3595890765353,84.3375642292315,80.3155393819278,76.2935145346241,72.2714896873204,68.2494648400166,64.2274399927129,60.2054151454092,56.1833902981055,52.1613654508018,48.139340603498,44.1173157561943,40.0952909088906,36.0732660615869,32.0512412142831],[127.246539712247,123.224514864943,119.202490017639,115.180465170336,111.158440323032,107.136415475728,103.114390628424,99.0923657811208,95.070340933817,91.0483160865133,87.0262912392096,83.0042663919059,78.9822415446021,74.9602166972984,70.9381918499947,66.916167002691,62.8941421553872,58.8721173080835,54.8500924607798,50.8280676134761,46.8060427661724,42.7840179188686,38.7619930715649,34.7399682242612,30.7179433769575],[125.913241874921,121.891217027617,117.869192180314,113.84716733301,109.825142485706,105.803117638403,101.781092791099,97.7590679437951,93.7370430964913,89.7150182491876,85.6929934018839,81.6709685545802,77.6489437072764,73.6269188599727,69.604894012669,65.5828691653653,61.5608443180616,57.5388194707578,53.5167946234541,49.4947697761504,45.4727449288467,41.4507200815429,37.4286952342392,33.4066703869355,29.3846455396318]],"type":"surface","x":[0,0.0833333333333333,0.166666666666667,0.25,0.333333333333333,0.416666666666667,0.5,0.583333333333333,0.666666666666667,0.75,0.833333333333333,0.916666666666667,1,1.08333333333333,1.16666666666667,1.25,1.33333333333333,1.41666666666667,1.5,1.58333333333333,1.66666666666667,1.75,1.83333333333333,1.91666666666667,2],"y":[0,0.0416666666666667,0.0833333333333333,0.125,0.166666666666667,0.208333333333333,0.25,0.291666666666667,0.333333333333333,0.375,0.416666666666667,0.458333333333333,0.5,0.541666666666667,0.583333333333333,0.625,0.666666666666667,0.708333333333333,0.75,0.791666666666667,0.833333333333333,0.875,0.916666666666667,0.958333333333333,1],"opacity":0.7,"surfacecolor":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],"cauto":false,"cmax":1,"cmin":0,"legendgroup":"High Weight","name":"High Weight","frame":null},{"colorbar":{"title":"Hz<br />Lz<br />Weight_group","ticklen":2},"colorscale":[["0","rgba(255,0,0,1)"],["0.0416666666666667","rgba(251,0,23,1)"],["0.0833333333333333","rgba(248,0,37,1)"],["0.125","rgba(244,0,49,1)"],["0.166666666666667","rgba(240,0,60,1)"],["0.208333333333333","rgba(236,0,70,1)"],["0.25","rgba(232,0,80,1)"],["0.291666666666667","rgba(228,0,89,1)"],["0.333333333333333","rgba(223,0,99,1)"],["0.375","rgba(218,0,108,1)"],["0.416666666666667","rgba(213,0,118,1)"],["0.458333333333333","rgba(207,0,127,1)"],["0.5","rgba(202,0,136,1)"],["0.541666666666667","rgba(195,0,146,1)"],["0.583333333333333","rgba(188,0,156,1)"],["0.625","rgba(181,0,165,1)"],["0.666666666666667","rgba(173,0,175,1)"],["0.708333333333333","rgba(164,0,185,1)"],["0.75","rgba(154,0,195,1)"],["0.791666666666667","rgba(143,0,204,1)"],["0.833333333333333","rgba(130,0,214,1)"],["0.875","rgba(115,0,224,1)"],["0.916666666666667","rgba(96,0,235,1)"],["0.958333333333333","rgba(69,0,245,1)"],["1","rgba(0,0,255,1)"]],"showscale":false,"z":[[156.493027495406,150.783929208043,145.074830920679,139.365732633316,133.656634345952,127.947536058589,122.238437771225,116.529339483862,110.820241196498,105.111142909135,99.4020446217711,93.6929463344076,87.9838480470441,82.2747497596806,76.5656514723171,70.8565531849536,65.1474548975901,59.4383566102266,53.7292583228631,48.0201600354996,42.3110617481361,36.6019634607726,30.8928651734091,25.1837668860456,19.4746685986821],[154.108900608252,148.399802320889,142.690704033525,136.981605746162,131.272507458798,125.563409171435,119.854310884071,114.145212596708,108.436114309344,102.727016021981,97.0179177346173,91.3088194472538,85.5997211598903,79.8906228725268,74.1815245851633,68.4724262977998,62.7633280104363,57.0542297230728,51.3451314357092,45.6360331483457,39.9269348609822,34.2178365736187,28.5087382862552,22.7996399988917,17.0905417115282],[151.724773721099,146.015675433735,140.306577146371,134.597478859008,128.888380571644,123.179282284281,117.470183996917,111.761085709554,106.05198742219,100.342889134827,94.6337908474634,88.9246925600999,83.2155942727364,77.5064959853729,71.7973976980094,66.0882994106459,60.3792011232824,54.6701028359189,48.9610045485554,43.2519062611919,37.5428079738284,31.8337096864649,26.1246113991014,20.4155131117379,14.7064148243744],[149.340646833945,143.631548546581,137.922450259218,132.213351971854,126.504253684491,120.795155397127,115.086057109764,109.3769588224,103.667860535037,97.9587622476731,92.2496639603096,86.5405656729461,80.8314673855826,75.1223690982191,69.4132708108556,63.7041725234921,57.9950742361286,52.2859759487651,46.5768776614016,40.8677793740381,35.1586810866746,29.449582799311,23.7404845119476,18.031386224584,12.3222879372205],[146.956519946791,141.247421659427,135.538323372064,129.8292250847,124.120126797337,118.411028509973,112.70193022261,106.992831935246,101.283733647883,95.5746353605193,89.8655370731558,84.1564387857923,78.4473404984288,72.7382422110653,67.0291439237018,61.3200456363383,55.6109473489747,49.9018490616112,44.1927507742477,38.4836524868842,32.7745541995207,27.0654559121572,21.3563576247937,15.6472593374302,9.9381610500667],[144.572393059637,138.863294772273,133.15419648491,127.445098197546,121.735999910183,116.026901622819,110.317803335456,104.608705048092,98.8996067607289,93.1905084733654,87.4814101860019,81.7723118986384,76.0632136112749,70.3541153239114,64.6450170365479,58.9359187491844,53.2268204618209,47.5177221744574,41.8086238870939,36.0995255997304,30.3904273123669,24.6813290250034,18.9722307376399,13.2631324502764,7.55403416291286],[142.188266172483,136.47916788512,130.770069597756,125.060971310393,119.351873023029,113.642774735666,107.933676448302,102.224578160939,96.5154798735751,90.8063815862116,85.0972832988481,79.3881850114846,73.6790867241211,67.9699884367576,62.2608901493941,56.5517918620306,50.8426935746671,45.1335952873036,39.42449699994,33.7153987125765,28.006300425213,22.2972021378495,16.588103850486,10.8790055631225,5.16990727575902],[139.804139285329,134.095040997966,128.385942710602,122.676844423239,116.967746135875,111.258647848512,105.549549561148,99.8404512737848,94.1313529864213,88.4222546990578,82.7131564116943,77.0040581243308,71.2949598369673,65.5858615496037,59.8767632622402,54.1676649748767,48.4585666875132,42.7494684001497,37.0403701127862,31.3312718254227,25.6221735380592,19.9130752506957,14.2039769633322,8.49487867596869,2.78578038860518],[137.420012398175,131.710914110812,126.001815823448,120.292717536085,114.583619248721,108.874520961358,103.165422673994,97.4563243866309,91.7472260992674,86.0381278119039,80.3290295245404,74.6199312371769,68.9108329498134,63.2017346624499,57.4926363750864,51.7835380877229,46.0744398003594,40.3653415129959,34.6562432256324,28.9471449382689,23.2380466509054,17.5289483635419,11.8198500761784,6.11075178881485,0.401653501451342],[135.035885511022,129.326787223658,123.617688936295,117.908590648931,112.199492361568,106.490394074204,100.781295786841,95.0721974994771,89.3630992121136,83.6540009247501,77.9449026373866,72.2358043500231,66.5267060626596,60.8176077752961,55.1085094879326,49.3994112005691,43.6903129132056,37.981214625842,32.2721163384785,26.563018051115,20.8539197637515,15.144821476388,9.43572318902452,3.726624901661,-1.9824733857025],[132.651758623868,126.942660336504,121.233562049141,115.524463761777,109.815365474414,104.10626718705,98.3971688996868,92.6880706123233,86.9789723249598,81.2698740375962,75.5607757502327,69.8516774628692,64.1425791755057,58.4334808881422,52.7243826007787,47.0152843134152,41.3061860260517,35.5970877386882,29.8879894513247,24.1788911639612,18.4697928765977,12.7606945892342,7.05159630187068,1.34249801450716,-4.36660027285634],[130.267631736714,124.55853344935,118.849435161987,113.140336874623,107.43123858726,101.722140299896,96.0130420125329,90.3039437251694,84.5948454378059,78.8857471504424,73.1766488630789,67.4675505757154,61.7584522883519,56.0493540009884,50.3402557136249,44.6311574262614,38.9220591388979,33.2129608515344,27.5038625641709,21.7947642768073,16.0856659894438,10.3765677020803,4.66746941471684,-1.04162887264668,-6.75072716001018],[127.88350484956,122.174406562197,116.465308274833,110.75620998747,105.047111700106,99.3380134127426,93.6289151253791,87.9198168380156,82.2107185506521,76.5016202632885,70.792521975925,65.0834236885615,59.374325401198,53.6652271138345,47.956128826471,42.2470305391075,36.537932251744,30.8288339643805,25.119735677017,19.4106373896535,13.70153910229,7.99244081492649,2.283342527563,-3.42575575980052,-9.13485404716402],[125.499377962406,119.790279675043,114.081181387679,108.372083100316,102.662984812952,96.9538865255887,91.2447882382252,85.5356899508617,79.8265916634982,74.1174933761347,68.4083950887712,62.6992968014077,56.9901985140442,51.2811002266807,45.5720019393172,39.8629036519537,34.1538053645902,28.4447070772267,22.7356087898632,17.0265105024997,11.3174122151362,5.60831392777265,-0.100784359590836,-5.80988264695435,-11.5189809343179],[123.115251075252,117.406152787889,111.697054500525,105.987956213162,100.278857925798,94.5697596384349,88.8606613510714,83.1515630637079,77.4424647763444,71.7333664889809,66.0242682016174,60.3151699142539,54.6060716268904,48.8969733395269,43.1878750521634,37.4787767647999,31.7696784774364,26.0605801900729,20.3514819027093,14.6423836153458,8.93328532798233,3.22418704061882,-2.48491124674467,-8.19400953410819,-13.9031078214717],[120.731124188099,115.022025900735,109.312927613372,103.603829326008,97.8947310386446,92.185632751281,86.4765344639175,80.767436176554,75.0583378891905,69.349239601827,63.6401413144635,57.9310430271,52.2219447397365,46.512846452373,40.8037481650095,35.094649877646,29.3855515902825,23.676453302919,17.9673550155555,12.258256728192,6.54915844082849,0.840060153464968,-4.86903813389852,-10.578136421262,-16.2872347086255],[118.346997300945,112.637899013581,106.928800726218,101.219702438854,95.5106041514907,89.8015058641272,84.0924075767637,78.3833092894002,72.6742110020367,66.9651127146732,61.2560144273097,55.5469161399462,49.8378178525827,44.1287195652192,38.4196212778557,32.7105229904922,27.0014247031287,21.2923264157652,15.5832281284017,9.87412984103815,4.16503155367465,-1.54406673368887,-7.25316502105235,-12.9622633084159,-18.6713615957794],[115.962870413791,110.253772126427,104.544673839064,98.8355755517004,93.1264772643369,87.4173789769734,81.7082806896099,75.9991824022464,70.2900841148829,64.5809858275194,58.8718875401559,53.1627892527924,47.4536909654288,41.7445926780653,36.0354943907018,30.3263961033383,24.6172978159748,18.9081995286113,13.1991012412478,7.49000295388431,1.78090466652081,-3.92819362084271,-9.63729190820619,-15.3463901955697,-21.0554884829332],[113.578743526637,107.869645239274,102.16054695191,96.4514486645466,90.7423503771831,85.0332520898195,79.324153802456,73.6150555150925,67.905957227729,62.1968589403655,56.487760653002,50.7786623656385,45.069564078275,39.3604657909115,33.651367503548,27.9422692161845,22.233170928821,16.5240726414575,10.814974354094,5.10587606673047,-0.603222220633029,-6.31232050799655,-12.02141879536,-17.7305170827236,-23.4396153700871],[111.194616639483,105.48551835212,99.7764200647562,94.0673217773927,88.3582234900292,82.6491252026657,76.9400269153022,71.2309286279387,65.5218303405752,59.8127320532117,54.1036337658482,48.3945354784847,42.6854371911212,36.9763389037577,31.2672406163942,25.5581423290307,19.8490440416672,14.1399457543037,8.43084746694014,2.72174917957663,-2.98734910778687,-8.69644739515039,-14.4055456825139,-20.1146439698774,-25.8237422572409],[108.810489752329,103.101391464966,97.3922931776024,91.6831948902389,85.9740966028754,80.2649983155119,74.5559000281484,68.8468017407849,63.1377034534214,57.4286051660578,51.7195068786943,46.0104085913308,40.3013103039673,34.5922120166038,28.8831137292403,23.1740154418768,17.4649171545133,11.7558188671498,6.0467205797863,0.337622292422793,-5.37147599494071,-11.0805742823042,-16.7896725696677,-22.4987708570312,-28.2078691443947],[106.426362865176,100.717264577812,95.0081662904485,89.299068003085,83.5899697157215,77.880871428358,72.1717731409945,66.462674853631,60.7535765662675,55.044478278904,49.3353799915405,43.626281704177,37.9171834168135,32.20808512945,26.4989868420865,20.789888554723,15.0807902673595,9.37169197999597,3.66259369263246,-2.04650459473105,-7.75560288209455,-13.4647011694581,-19.1737994568216,-24.8828977441851,-30.5919960315486],[104.042235978022,98.3331376906582,92.6240394032947,86.9149411159312,81.2058428285677,75.4967445412042,69.7876462538407,64.0785479664772,58.3694496791137,52.6603513917502,46.9512531043867,41.2421548170231,35.5330565296596,29.8239582422961,24.1148599549326,18.4057616675691,12.6966633802056,6.98756509284213,1.27846680547862,-4.43063148188489,-10.1397297692484,-15.8488280566119,-21.5579263439754,-27.2670246313389,-32.9761229187024],[101.658109090868,95.9490108035044,90.2399125161409,84.5308142287774,78.8217159414139,73.1126176540504,67.4035193666869,61.6944210793233,55.9853227919598,50.2762245045963,44.5671262172328,38.8580279298693,33.1489296425058,27.4398313551423,21.7307330677788,16.0216347804153,10.3125364930518,4.6034382056883,-1.10566008167522,-6.81475836903872,-12.5238566564022,-18.2329549437657,-23.9420532311292,-29.6511515184927,-35.3602498058562],[99.273982203714,93.5648839163505,87.855785628987,82.1466873416235,76.43758905426,70.7284907668965,65.019392479533,59.3102941921695,53.601195904806,47.8920976174425,42.182999330079,36.4739010427155,30.764802755352,25.0557044679885,19.346606180625,13.6375078932615,7.92840960589795,2.21931131853445,-3.48978696882907,-9.19888525619257,-14.9079835435561,-20.6170818309196,-26.3261801182831,-32.0352784056466,-37.7443766930101]],"type":"surface","x":[0,0.0833333333333333,0.166666666666667,0.25,0.333333333333333,0.416666666666667,0.5,0.583333333333333,0.666666666666667,0.75,0.833333333333333,0.916666666666667,1,1.08333333333333,1.16666666666667,1.25,1.33333333333333,1.41666666666667,1.5,1.58333333333333,1.66666666666667,1.75,1.83333333333333,1.91666666666667,2],"y":[0,0.0416666666666667,0.0833333333333333,0.125,0.166666666666667,0.208333333333333,0.25,0.291666666666667,0.333333333333333,0.375,0.416666666666667,0.458333333333333,0.5,0.541666666666667,0.583333333333333,0.625,0.666666666666667,0.708333333333333,0.75,0.791666666666667,0.833333333333333,0.875,0.916666666666667,0.958333333333333,1],"opacity":0.7,"surfacecolor":[[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]],"cauto":false,"cmax":1,"cmin":0,"legendgroup":"Low Weight","name":"Low Weight","frame":null},{"x":[0,0,0,1,0,1,2,1,0,0,0,2,0,0,0,0,0,1,1,0,1],"y":[0,1,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1],"z":[100,102.142857142857,189,11,133,13,74.4285714285714,84,199.111111111111,109.666666666667,62,53.2,234.176470588235,53,64,207.1,110.142857142857,178.5,121,188,34],"type":"scatter3d","mode":"markers","legendgroup":"High Weight","name":"High Weight","marker":{"color":"rgba(255,0,0,1)","line":{"color":"rgba(255,0,0,1)"}},"textfont":{"color":"rgba(255,0,0,1)"},"error_y":{"color":"rgba(255,0,0,1)"},"error_x":{"color":"rgba(255,0,0,1)"},"line":{"color":"rgba(255,0,0,1)"},"frame":null},{"x":[0,1,0,0,0,0,0,1,0,0,1,2,2,0,1,0,1,0,1,0,0,0,0,0],"y":[0,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1],"z":[44,48,204.235294117647,116.625,224.052631578947,193.368421052632,22,182.076923076923,84,10,38,18,26,144.428571428571,44,213,63.75,104.875,26,106.230769230769,15.3333333333333,47.4285714285714,155.6,25],"type":"scatter3d","mode":"markers","legendgroup":"Low Weight","name":"Low Weight","marker":{"color":"rgba(0,0,255,1)","line":{"color":"rgba(0,0,255,1)"}},"textfont":{"color":"rgba(0,0,255,1)"},"error_y":{"color":"rgba(0,0,255,1)"},"error_x":{"color":"rgba(0,0,255,1)"},"line":{"color":"rgba(0,0,255,1)"},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
save_image(p, file = "Fig_1C_Antibiotic_planes_alpha_diversity.png")

diversity_Alpha_with_M$Carbapenem_01 <- case_when(
  diversity_Alpha_with_M$Carbapenem > 0 ~ "Treated",
  TRUE ~ "Not treated"
)
diversity_Alpha_with_M$Vancomycin_01 <- case_when(
  diversity_Alpha_with_M$Vancomycin > 0 ~ "Treated",
  TRUE ~ "Not treated"
)
```

## LEfSe {#LEfSe}

ASVs enriched in vancomycin and carbapenem treated subjects.


```r
#LEfSe to see which ASVs were enriched after antibiotic treatment
Abx_SE <- SummarizedExperiment(assays = ASV_rare , rowData = TAX_rare, colData = diversity_Alpha_with_M)

#Fig_1D_Vancomycin_LEfSe.pdf (saved as pdf 5x2)
Abx_lefse <- lefse(expr = Abx_SE, groupCol = "Vancomycin_01")
tax <- TAX_rare %>% filter(rownames(.) %in% Abx_lefse$Names);tax <- tax[Abx_lefse$Names,];tax <- tax$Genus
Abx_lefse$Names <- paste0(Abx_lefse$Names, "-", tax)
lefse_plot(Abx_lefse, colors = c("forestgreen", "red"))
```

![](Gut_microbiota_alpha_diversity_decrease_analyses_files/figure-html/LEfSe-1.png)<!-- -->

```r
#Fig_1D_Carbapenem_LEfSe (saved as pdf 5x2)
Abx_lefse <- lefse(expr = Abx_SE, groupCol = "Carbapenem_01")
tax <- TAX_rare %>% filter(rownames(.) %in% Abx_lefse$Names);tax <- tax[Abx_lefse$Names,];tax <- tax$Genus
Abx_lefse$Names <- paste0(Abx_lefse$Names, "-", tax)
lefse_plot(Abx_lefse, colors = c("forestgreen", "red"))
```

![](Gut_microbiota_alpha_diversity_decrease_analyses_files/figure-html/LEfSe-2.png)<!-- -->

## Cohort demographics {#Cohort_demographics}


```r
META_demographics <- META

META_demographics$Race <- case_when(grepl("1", META_demographics$Race) ~ "White",
                                   grepl("2", META_demographics$Race) ~ "Black/African American",
                                   grepl("3", META_demographics$Race) ~ "Asian",
                                   grepl("4", META_demographics$Race) ~ "Other", TRUE ~ META_demographics$Race)

META_demographics$Smoking <- case_when(grepl("0", META_demographics$Smoking) ~ "No history",
                                       grepl("1", META_demographics$Smoking) ~ "Past history",
                                       grepl("2", META_demographics$Smoking) ~ "Current smoker",
                                       grepl("unk", META_demographics$Smoking) ~ "Unknown",
                                       TRUE ~ META_demographics$Smoking)
META_demographics$Smoking <- factor(META_demographics$Smoking, levels = unique(META_demographics$Smoking))

META_demographics$Alcohol <- case_when(grepl("0", META_demographics$Alcohol) ~ "No history",
                                       grepl("1", META_demographics$Alcohol) ~ "Past history",
                                       grepl("2", META_demographics$Alcohol) ~ "Current drinker",
                                       grepl("unk", META_demographics$Alcohol) ~ "Unknown",
                                       TRUE ~ META_demographics$Alcohol)
META_demographics$Alcohol <- factor(META_demographics$Alcohol, levels = unique(META_demographics$Alcohol))

META_demographics$MDRO <- case_when(grepl("12", META_demographics$MDRO_organism) ~ "Negative",
                                             grepl("13", META_demographics$MDRO_organism) ~ "No screen for MDRO",
                                             TRUE ~ "Positive")
META_demographics$MDRO <- factor(META_demographics$MDRO, levels = unique(META_demographics$MDRO))

META_demographics$Exposure_to_antibiotics_during_visit <- case_when(grepl("No", META_demographics$Exposure_to_antibiotics_during_visit) ~ "No",
                                       grepl("Given", META_demographics$Exposure_to_antibiotics_during_visit) ~ "Yes")

META_demographics$COVID <- case_when(grepl("Negative", META_demographics$COVID19_Pos) ~ "Negative",
                                       grepl("Positive", META_demographics$COVID19_Pos) ~ "Positive")
META_demographics$C_diff <- case_when(grepl("negative", META_demographics$C_diff) ~ "Negative",
                                       grepl("positive", META_demographics$C_diff) ~ "Positive")
META_demographics$Mortality_90_day <- case_when(grepl("0", META_demographics$Mortality_90_day) ~ "No",
                                       grepl("1", META_demographics$Mortality_90_day) ~ "Yes",
                                       TRUE ~ "NA")

# Specify explanatory variables of interest
explanatory <- c("Age", "Weight_kg", "BMI", "Length_of_stay","Sex",
                "Race", "Smoking", "Alcohol", "Exposure_to_antibiotics_during_visit", "COVID", "C_diff", "MDRO",  "Mortality_90_day")

demographics <- META_demographics %>% 
  finalfit::summary_factorlist("Case_Control", explanatory,
                     p=TRUE, na_include=TRUE)

knitr::kable(demographics, format = "simple")
```



label                                  levels                   Case          Control       p      
-------------------------------------  -----------------------  ------------  ------------  -------
Age                                    Mean (SD)                60.5 (13.9)   63.0 (16.9)   0.630  
Weight_kg                              Mean (SD)                83.3 (22.9)   80.6 (37.1)   0.775  
BMI                                    Mean (SD)                28.6 (7.9)    28.3 (13.3)   0.936  
Length_of_stay                         Mean (SD)                26.3 (20.6)   16.7 (9.6)    0.162  
Sex                                    Female                   16 (45.7)     3 (30.0)      0.600  
                                       Male                     19 (54.3)     7 (70.0)             
Race                                   Black/African American   25 (71.4)     6 (60.0)      0.528  
                                       Declined                 1 (2.9)                            
                                       Other                    2 (5.7)       2 (20.0)             
                                       White                    7 (20.0)      2 (20.0)             
Smoking                                No history               19 (54.3)     2 (20.0)      0.001  
                                       Current smoker           14 (40.0)     4 (40.0)             
                                       Unknown                  0 (0.0)       4 (40.0)             
                                       Past history             2 (5.7)       0 (0.0)              
Alcohol                                No history               23 (65.7)     4 (40.0)      0.001  
                                       Current drinker          10 (28.6)     2 (20.0)             
                                       Unknown                  0 (0.0)       4 (40.0)             
                                       Past history             2 (5.7)       0 (0.0)              
Exposure_to_antibiotics_during_visit   No                       2 (5.7)       2 (20.0)      0.441  
                                       Yes                      33 (94.3)     8 (80.0)             
COVID                                  Negative                 21 (61.8)     10 (100.0)    0.053  
                                       Positive                 13 (38.2)                          
C_diff                                 Negative                 24 (68.6)     10 (100.0)    0.105  
                                       Positive                 11 (31.4)                          
MDRO                                   Negative                 4 (11.4)      10 (100.0)    <0.001 
                                       Positive                 23 (65.7)     0 (0.0)              
                                       No screen for MDRO       8 (22.9)      0 (0.0)              
Mortality_90_day                       Missing                  5 (14.3)                    0.414  
                                       No                       20 (57.1)     6 (60.0)             
                                       Yes                      10 (28.6)     4 (40.0)             

```r
write.csv(demographics, file = "Table_1_Cohort_Demographics.csv")
```
