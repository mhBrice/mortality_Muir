# Mortality analyses in an old-growth forest

![DOI](https://zenodo.org/badge/191195897.svg)](https://zenodo.org/badge/latestdoi/191195897)

This repository includes the data and R scripts to reproduce analyses and figures found in the article *Long-term impact of a major ice storm on tree mortality in an old-growth forest* by Deschênes, Brice and Brisson published in Forest Ecology and Management.


## Installation

The analyses were carried out with [R version 3.5.1 (a free software environment for statistical computing and graphics)](https://www.r-project.org/) and require the installation of a recent version of it.

Analyses were reproduced in the MacOSX Mojave environment:

<details>
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin18.0.0 (64-bit)
Running under: macOS  10.14.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] units_0.6-1         sf_0.7-1            graphicsutils_1.2-1 scales_1.0.0       
 [5] colourlovers_0.2.2  RColorBrewer_1.1-2  rms_5.1-3           SparseM_1.77       
 [9] Hmisc_4.2-0         Formula_1.2-3       lattice_0.20-35     sjPlot_2.6.1       
[13] survminer_0.4.3     ggpubr_0.2          magrittr_1.5        ggplot2_3.1.0      
[17] survival_2.42-3     dplyr_0.7.8         gtools_3.8.1       

loaded via a namespace (and not attached):
 [1] TH.data_1.0-9       minqa_1.2.4         colorspace_1.4-0    class_7.3-14       
 [5] modeltools_0.2-22   ggridges_0.5.1      sjlabelled_1.0.14   estimability_1.3   
 [9] snakecase_0.9.2     htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.9.0   
[13] glmmTMB_0.2.2.0     MatrixModels_0.4-1  mvtnorm_1.0-8       coin_1.2-2         
[17] codetools_0.2-15    splines_3.5.1       mnormt_1.5-5        knitr_1.21         
[21] sjmisc_2.7.6        bayesplot_1.6.0     jsonlite_1.6        nloptr_1.2.1       
[25] ggeffects_0.7.0     broom_0.5.0         km.ci_0.5-2         cluster_2.0.7-1    
[29] png_0.1-7           compiler_3.5.1      sjstats_0.17.2      emmeans_1.3.0      
[33] backports_1.1.3     assertthat_0.2.0    Matrix_1.2-14       lazyeval_0.2.1     
[37] acepack_1.4.1       htmltools_0.3.6     quantreg_5.36       tools_3.5.1        
[41] bindrcpp_0.2.2      coda_0.19-2         gtable_0.2.0        glue_1.3.0         
[45] Rcpp_1.0.0          nlme_3.1-137        psych_1.8.10        xfun_0.4           
[49] stringr_1.3.1       lme4_1.1-19         XML_3.98-1.16       stringdist_0.9.5.1
[53] polspline_1.1.13    MASS_7.3-50         zoo_1.8-4           hms_0.4.2          
[57] parallel_3.5.1      sandwich_2.5-0      pwr_1.2-2           TMB_1.7.15         
[61] yaml_2.2.0          gridExtra_2.3       KMsurv_0.1-5        rpart_4.1-13       
[65] latticeExtra_0.6-28 stringi_1.2.4       e1071_1.7-0         checkmate_1.9.1    
[69] spData_0.2.9.4      rlang_0.3.1         pkgconfig_2.0.2     purrr_0.2.5        
[73] prediction_0.3.6    bindr_0.1.1         htmlwidgets_1.3     cmprsk_2.2-7       
[77] tidyselect_0.2.5    plyr_1.8.4          R6_2.3.0            multcomp_1.4-8     
[81] DBI_1.0.0           pillar_1.3.1        haven_1.1.2         foreign_0.8-70     
[85] withr_2.1.2         nnet_7.3-12         tibble_2.0.1        modelr_0.1.2       
[89] crayon_1.3.4        survMisc_0.5.5      grid_3.5.1          data.table_1.12.0  
[93] forcats_0.3.0       classInt_0.2-3      digest_0.6.18       xtable_1.8-3       
[97] tidyr_0.8.2         stats4_3.5.1        munsell_0.5.0     
</details>

## The following packages must be installed to run the scripts:

- survival
- survminer
- dplyr
- sf
- units
- gtools
- rms
- sjPlot
- RColorBrewer
- colourlovers
- scales
- [graphicsutils](https://github.com/inSileco/graphicsutils) (not on CRAN)

Below are the R commands to install them all:

```R
install.packages(
  c("survival", "survminer", "dplyr", "sf", "units", "gtools", "rms", "sjPlot",
  "RColorBrewer", "colourlovers", "scales", "remotes")
)
remotes::install_github("inSileco/graphicsutils")
```

## Details

To reproduce the entire analysis including data cleaning, analyses and figures, run:

```R
source("scripts/1_dataFormatting.R")
source("scripts/2_coxph.R")
```

All data used for the analyses can be found in the [data](https://github.com/mhBrice/mortality_Muir/tree/master/data) folder.

Script #1, `scripts/1_dataFormatting.R`, prepares and formats the data for the analyses. The cleaned data are also included in the [data](https://github.com/mhBrice/mortality_Muir/tree/master/data) folder, hence the first script (`scripts/1_dataFormatting.R`) can be skipped.

Script #2, `scripts/2_coxph.R`, performed all the analyses and produced the figures.

The figures and tables produced are saved in the [results](https://github.com/mhBrice/mortality_Muir/tree/master/results) folder.
