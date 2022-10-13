Scripts for re-running and recreating figures for the manuscript titled: "When traits vary across species but performance doesn't: One solution to the paradox"

All figures saved in \Output\Figures

To recreate figures and results in order they are presented in the manuscript:

1. Data_Wrangling.R - this will import the original data file used in Nolting et al. (2021) and trim and clean it to create a Protea_data object that will be used in subsequent analyses for the current manuscript. 

2. multiple_regressions.R - this script will run through the models to fit the structural trait-physiological performance trait multiple regressions, demonstrating that structural traits explain variation in performance well. *Figure 1* in the manuscript is generated from this script (legend added manually). All model fit files have been saved for reproducibility of results.

3. variance_partitioning.R - this script will run through the models to perform a variance partitioning on all structural traits and on all physiological performance traits. *Figure 2* in the manuscript is generated from this script. All model fit files have been saved for reproducibility of results.

4. two-trait-model_simulation.R - this script runs through the simple two-trait model to simulate trait data under scenario 1 and scenario 2, as presented in the manuscript. *Figure 3* in the manuscript is generated from this script.

5. pc-functional-axis-alignment.R - this script runs through fitting multiple regression models as in 2 above, but without the random effects (for ease of comparison with the generalized model presented in the text). Using this model output, the script then runs through calculating the alignment between each structural trait PC axis and the functional axis corresponding to each physiological performance trait. *Figure 4* in the manuscript is generated from this script.

6. pca-performance_plots.R - this script runs through plotting the Protea_data in structural trait space, and coloring points based on physiological performance. *Figure 5* in the manuscript is generated from this script.



## All session info reported below:

R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.4

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3   mvtnorm_1.1-3        tidybayes_3.0.2      cowplot_1.1.1        brms_2.18.0         
 [6] rstanarm_2.21.3      Rcpp_1.0.9           rstan_2.21.7         StanHeaders_2.21.0-7 forcats_0.5.2       
[11] stringr_1.4.1        purrr_0.3.5          readr_2.1.3          tidyr_1.2.1          tibble_3.1.8        
[16] ggplot2_3.3.6        tidyverse_1.3.2      dplyr_1.0.10        

loaded via a namespace (and not attached):
  [1] googledrive_2.0.0    minqa_1.2.4          colorspace_2.0-3     ellipsis_0.3.2       ggridges_0.5.4      
  [6] markdown_1.1         base64enc_0.1-3      fs_1.5.2             rstudioapi_0.14      farver_2.1.1        
 [11] svUnit_1.0.6         DT_0.25              fansi_1.0.3          lubridate_1.8.0      xml2_1.3.3          
 [16] bridgesampling_1.1-2 codetools_0.2-18     splines_4.2.1        shinythemes_1.2.0    bayesplot_1.9.0     
 [21] jsonlite_1.8.0       nloptr_2.0.3         broom_1.0.1          dbplyr_2.2.1         ggdist_3.2.0        
 [26] shiny_1.7.2          compiler_4.2.1       httr_1.4.4           backports_1.4.1      assertthat_0.2.1    
 [31] Matrix_1.4-1         fastmap_1.1.0        gargle_1.2.1         cli_3.4.1            later_1.3.0         
 [36] htmltools_0.5.3      prettyunits_1.1.1    tools_4.2.1          igraph_1.3.5         coda_0.19-4         
 [41] gtable_0.3.1         glue_1.6.2           posterior_1.3.1      reshape2_1.4.4       cellranger_1.1.0    
 [46] vctrs_0.4.1          nlme_3.1-157         crosstalk_1.2.0      tensorA_0.36.2       ps_1.7.1            
 [51] lme4_1.1-30          rvest_1.0.3          mime_0.12            miniUI_0.1.1.1       lifecycle_1.0.2     
 [56] gtools_3.9.3         googlesheets4_1.0.1  MASS_7.3-57          zoo_1.8-11           scales_1.2.1        
 [61] colourpicker_1.1.1   ragg_1.2.3           Brobdingnag_1.2-7    hms_1.1.2            promises_1.2.0.1    
 [66] parallel_4.2.1       inline_0.3.19        shinystan_2.6.0      gridExtra_2.3        loo_2.5.1           
 [71] stringi_1.7.8        dygraphs_1.1.1.6     checkmate_2.1.0      boot_1.3-28          pkgbuild_1.3.1      
 [76] systemfonts_1.0.4    rlang_1.0.6          pkgconfig_2.0.3      matrixStats_0.62.0   distributional_0.3.1
 [81] lattice_0.20-45      labeling_0.4.2       rstantools_2.2.0     htmlwidgets_1.5.4    processx_3.7.0      
 [86] tidyselect_1.1.2     plyr_1.8.7           magrittr_2.0.3       R6_2.5.1             generics_0.1.3      
 [91] DBI_1.1.3            mgcv_1.8-40          pillar_1.8.1         haven_2.5.1          withr_2.5.0         
 [96] xts_0.12.1           survival_3.3-1       abind_1.4-5          modelr_0.1.9         crayon_1.5.1        
[101] arrayhelpers_1.1-0   utf8_1.2.2           tzdb_0.3.0           grid_4.2.1           readxl_1.4.1        
[106] callr_3.7.2          threejs_0.3.3        reprex_2.0.2         digest_0.6.29        xtable_1.8-4        
[111] httpuv_1.6.6         textshaping_0.3.6    RcppParallel_5.1.5   stats4_4.2.1         munsell_0.5.0       
[116] viridisLite_0.4.1    shinyjs_2.1.0       