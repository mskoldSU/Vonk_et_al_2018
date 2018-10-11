Vonk et al temporal deconvolution
================

This repository contains sample code for reproducing results in *Temporal deconvolution of vascular plant-derived fatty acids exported from terrestial watersheds* by Vonk et.al., to appear in [*Geochimica et Cosmochimica Acta*](https://www.journals.elsevier.com/geochimica-et-cosmochimica-acta).

The model considered in the paper can be described by the equation (\[4\] in the paper)

<p align="center">
<img src="tex/d16301d632c377abbabd558b00415cf4.svg?invert_in_darkmode&sanitize=true" align=middle width=609.7429701pt height=38.242408049999995pt/>
</p>

for <img src="tex/e777e1ef6ee7e48f0276ad38390a560a.svg?invert_in_darkmode&sanitize=true" align=middle width=82.19635874999999pt height=21.68300969999999pt/>, where <img src="tex/b0a852ae71cd56ce00cb7c1a592a1edc.svg?invert_in_darkmode&sanitize=true" align=middle width=161.97330434999998pt height=24.65753399999998pt/> and <img src="tex/1cd32b0756da515bc59142b9318ff797.svg?invert_in_darkmode&sanitize=true" align=middle width=11.323291649999991pt height=14.15524440000002pt/> is a normally distributed (Gaussian) measurement error with known standard deviation. Prior to fitting, the integral is approximated by the sum

<p align="center">
<img src="tex/022e39b38d956874f3112ada2eaa29d2.svg?invert_in_darkmode&sanitize=true" align=middle width=592.0066052999999pt height=73.6915278pt/>
</p>

where <img src="tex/fb15704e32b633db97db1b14d37cc11e.svg?invert_in_darkmode&sanitize=true" align=middle width=125.00495204999999pt height=30.632847300000012pt/> and <img src="tex/3f6b797c2dfc37f90de17acc95eb614b.svg?invert_in_darkmode&sanitize=true" align=middle width=91.34822069999998pt height=31.36100879999999pt/>. Here, the second approximation uses that <img src="tex/06ec6f89507400dfedf2a279f5d04475.svg?invert_in_darkmode&sanitize=true" align=middle width=65.20135049999999pt height=26.76175259999998pt/> is essentially zero prior to the bomb-spike.

We set vague prior distributions as follows; <img src="tex/3aabf826ac9ea07e956a7a48ca601b94.svg?invert_in_darkmode&sanitize=true" align=middle width=29.46735659999999pt height=22.831056599999986pt/> is uniform on <img src="tex/acf5ce819219b95070be2dbeb8a671e9.svg?invert_in_darkmode&sanitize=true" align=middle width=32.87674994999999pt height=24.65753399999998pt/>, <img src="tex/1edede1e27fbb6f1b2b3ec665f22f03b.svg?invert_in_darkmode&sanitize=true" align=middle width=86.71228499999998pt height=30.632847300000012pt/> is <img src="tex/33c91016030fd67399f9d5d1920723e9.svg?invert_in_darkmode&sanitize=true" align=middle width=134.42240625pt height=24.65753399999998pt/> distributed, <img src="tex/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?invert_in_darkmode&sanitize=true" align=middle width=7.054796099999991pt height=22.831056599999986pt/> is <img src="tex/1117d0e319c13cf635b97a2ccf61f8b3.svg?invert_in_darkmode&sanitize=true" align=middle width=75.34258214999998pt height=26.76175259999998pt/> and finally <img src="tex/3697af5e84c4a5802cdf35c401c9590c.svg?invert_in_darkmode&sanitize=true" align=middle width=161.16873465pt height=27.91243950000002pt/> is apriori uniform on <img src="tex/0d7179a315f6448786cdd1dcd7e17b0d.svg?invert_in_darkmode&sanitize=true" align=middle width=70.31981219999999pt height=24.65753399999998pt/> (corresponding to an exponential prior for <img src="tex/933cc196ec5f448dc64a42cfdfc45064.svg?invert_in_darkmode&sanitize=true" align=middle width=28.605474149999992pt height=14.15524440000002pt/> with mean <img src="tex/1a871532f9a24c565a9cae2c3f30402f.svg?invert_in_darkmode&sanitize=true" align=middle width=26.02750259999999pt height=24.65753399999998pt/>).

We take a Bayesian approach and fit the model using [JAGS](http://mcmc-jags.sourceforge.net/) called from the [R](https://www.r-project.org/) computing environment. Model code can be found in file [`deconvolve_f.jags`](deconvolve_f.jags) and an MCMC sample from the posterior distribution of parameters is drawn using the helper function [`deconvolve_f`](functions.R). This is illustrated below for the Mackenzie data loaded together with the atmospheric data by [`load.R`](load.R):

``` r
source("load.R") # Load data and libraries (rjags requires installation of JAGS)
source("functions.R") # Load helper functions
mcmc_sample <- deconvolve_f(D14_Mackenzie, D14_atm) # Draw sample from posterior
```

The function [`draw_tau_hist`](functions.R) draws a sample from the posterior predictive distribution of residence times in the younger segment (0-50 years) and plots a histogram:

``` r
draw_tau_hist(mcmc_sample, 0:50) + ggtitle("Mackenzie delta, posterior predictive residence times (younger)")
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

In order to check model fit, [`plot_fit`](functions.R) plots the posterior mean curve together with a shaded 90% credibility interval and data with error bars (four standard deviations wide):

``` r
plot_fit(mcmc_sample, D14_Mackenzie, D14_atm) + ggtitle("Mackenzie delta, fitted curve (posterior mean) with data")
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

Finally [`table_pars`](functions.R) creates a table of posterior parameter summaries:

``` r
table_pars(mcmc_sample)
```

| parameter | mean (sd)          |
|:----------|:-------------------|
| b         | 1.33 (0.13)        |
| f\_Old    | 0.52 (0.02)        |
| tau\_Old  | 28386.23 (9424.77) |

``` r
sessionInfo()
```

    ## R version 3.4.4 (2018-03-15)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 16299)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Swedish_Sweden.1252  LC_CTYPE=Swedish_Sweden.1252   
    ## [3] LC_MONETARY=Swedish_Sweden.1252 LC_NUMERIC=C                   
    ## [5] LC_TIME=Swedish_Sweden.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2     rjags_4-6          coda_0.19-1       
    ##  [4] forcats_0.3.0      stringr_1.3.1      dplyr_0.7.6       
    ##  [7] purrr_0.2.5        readr_1.1.1        tidyr_0.8.1       
    ## [10] tibble_1.4.2       ggplot2_3.0.0.9000 tidyverse_1.2.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.4 haven_1.1.2      lattice_0.20-35  colorspace_1.3-2
    ##  [5] htmltools_0.3.6  yaml_2.2.0       rlang_0.2.2      pillar_1.3.0    
    ##  [9] glue_1.3.0       withr_2.1.2      modelr_0.1.2     readxl_1.1.0    
    ## [13] bindr_0.1.1      plyr_1.8.4       munsell_0.5.0    gtable_0.2.0    
    ## [17] cellranger_1.1.0 rvest_0.3.2      codetools_0.2-15 evaluate_0.11   
    ## [21] labeling_0.3     knitr_1.20       highr_0.7        broom_0.5.0     
    ## [25] Rcpp_0.12.19     scales_1.0.0     backports_1.1.2  jsonlite_1.5    
    ## [29] hms_0.4.2        digest_0.6.17    stringi_1.1.7    grid_3.4.4      
    ## [33] rprojroot_1.3-2  cli_1.0.1        tools_3.4.4      magrittr_1.5    
    ## [37] lazyeval_0.2.1   crayon_1.3.4     pkgconfig_2.0.2  xml2_1.2.0      
    ## [41] lubridate_1.7.4  assertthat_0.2.0 rmarkdown_1.10   httr_1.3.1      
    ## [45] rstudioapi_0.8   R6_2.2.2         nlme_3.1-131.1   compiler_3.4.4
