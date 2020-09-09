Bering Sea and Aleutian Islands Other Rockfish stock assessment
Contact info: jane.sullivan@noaa.gov

```
> session_info()
- Session info --------------------------------------------------------------------------------------------
 setting  value                       
 version  R version 4.0.2 (2020-06-22)
 os       Windows 10 x64              
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United States.1252  
 ctype    English_United States.1252  
 tz       America/Anchorage           
 date     2020-09-08                  

- Packages ------------------------------------------------------------------------------------------------
 package     * version  date       lib source        
 assertthat    0.2.1    2019-03-21 [1] CRAN (R 4.0.2)
 backports     1.1.7    2020-05-13 [1] CRAN (R 4.0.0)
 blob          1.2.1    2020-01-20 [1] CRAN (R 4.0.2)
 broom         0.7.0    2020-07-09 [1] CRAN (R 4.0.2)
 callr         3.4.3    2020-03-28 [1] CRAN (R 4.0.2)
 cellranger    1.1.0    2016-07-27 [1] CRAN (R 4.0.2)
 cli           2.0.2    2020-02-28 [1] CRAN (R 4.0.2)
 coda          0.19-3   2019-07-05 [1] CRAN (R 4.0.2)
 colorspace    1.4-1    2019-03-18 [1] CRAN (R 4.0.2)
 crayon        1.3.4    2017-09-16 [1] CRAN (R 4.0.2)
 DBI           1.1.0    2019-12-15 [1] CRAN (R 4.0.2)
 dbplyr        1.4.4    2020-05-27 [1] CRAN (R 4.0.2)
 desc          1.2.0    2018-05-01 [1] CRAN (R 4.0.2)
 devtools    * 2.3.1    2020-07-21 [1] CRAN (R 4.0.2)
 digest        0.6.25   2020-02-23 [1] CRAN (R 4.0.2)
 dplyr       * 1.0.1    2020-07-31 [1] CRAN (R 4.0.2)
 ellipsis      0.3.1    2020-05-15 [1] CRAN (R 4.0.2)
 evaluate      0.14     2019-05-28 [1] CRAN (R 4.0.2)
 fansi         0.4.1    2020-01-08 [1] CRAN (R 4.0.2)
 farver        2.0.3    2020-01-16 [1] CRAN (R 4.0.2)
 forcats     * 0.5.0    2020-03-01 [1] CRAN (R 4.0.2)
 fs            1.4.2    2020-06-30 [1] CRAN (R 4.0.2)
 generics      0.0.2    2018-11-29 [1] CRAN (R 4.0.2)
 ggplot2     * 3.3.2    2020-06-19 [1] CRAN (R 4.0.2)
 ggridges    * 0.5.2    2020-01-12 [1] CRAN (R 4.0.2)
 glue          1.4.1    2020-05-13 [1] CRAN (R 4.0.2)
 gtable        0.3.0    2019-03-25 [1] CRAN (R 4.0.2)
 haven         2.3.1    2020-06-01 [1] CRAN (R 4.0.2)
 hms           0.5.3    2020-01-08 [1] CRAN (R 4.0.2)
 htmltools     0.5.0    2020-06-16 [1] CRAN (R 4.0.2)
 httr          1.4.2    2020-07-20 [1] CRAN (R 4.0.2)
 jsonlite      1.7.0    2020-06-25 [1] CRAN (R 4.0.2)
 knitr         1.29     2020-06-23 [1] CRAN (R 4.0.2)
 labeling      0.3      2014-08-23 [1] CRAN (R 4.0.0)
 lattice       0.20-41  2020-04-02 [2] CRAN (R 4.0.2)
 lifecycle     0.2.0    2020-03-06 [1] CRAN (R 4.0.2)
 lubridate     1.7.9    2020-06-08 [1] CRAN (R 4.0.2)
 magrittr      1.5      2014-11-22 [1] CRAN (R 4.0.2)
 memoise       1.1.0    2017-04-21 [1] CRAN (R 4.0.2)
 modelr        0.1.8    2020-05-19 [1] CRAN (R 4.0.2)
 munsell       0.5.0    2018-06-12 [1] CRAN (R 4.0.2)
 pillar        1.4.6    2020-07-10 [1] CRAN (R 4.0.2)
 pkgbuild      1.1.0    2020-07-13 [1] CRAN (R 4.0.2)
 pkgconfig     2.0.3    2019-09-22 [1] CRAN (R 4.0.2)
 pkgload       1.1.0    2020-05-29 [1] CRAN (R 4.0.2)
 plyr          1.8.6    2020-03-03 [1] CRAN (R 4.0.2)
 prettyunits   1.1.1    2020-01-24 [1] CRAN (R 4.0.2)
 processx      3.4.3    2020-07-05 [1] CRAN (R 4.0.2)
 ps            1.3.3    2020-05-08 [1] CRAN (R 4.0.2)
 purrr       * 0.3.4    2020-04-17 [1] CRAN (R 4.0.2)
 R2admb      * 0.7.16.2 2020-03-09 [1] CRAN (R 4.0.2)
 R6            2.4.1    2019-11-12 [1] CRAN (R 4.0.2)
 Rcpp          1.0.5    2020-07-06 [1] CRAN (R 4.0.2)
 readr       * 1.3.1    2018-12-21 [1] CRAN (R 4.0.2)
 readxl        1.3.1    2019-03-13 [1] CRAN (R 4.0.2)
 remotes       2.2.0    2020-07-21 [1] CRAN (R 4.0.2)
 reprex        0.3.0    2019-05-16 [1] CRAN (R 4.0.2)
 rlang         0.4.7    2020-07-09 [1] CRAN (R 4.0.2)
 rmarkdown     2.3      2020-06-18 [1] CRAN (R 4.0.2)
 RODBC       * 1.3-17   2020-05-11 [1] CRAN (R 4.0.0)
 rprojroot     1.3-2    2018-01-03 [1] CRAN (R 4.0.2)
 rstudioapi    0.11     2020-02-07 [1] CRAN (R 4.0.2)
 rvest         0.3.6    2020-07-25 [1] CRAN (R 4.0.2)
 scales        1.1.1    2020-05-11 [1] CRAN (R 4.0.2)
 sessioninfo   1.1.1    2018-11-05 [1] CRAN (R 4.0.2)
 stringi       1.4.6    2020-02-17 [1] CRAN (R 4.0.0)
 stringr     * 1.4.0    2019-02-10 [1] CRAN (R 4.0.2)
 testthat      2.3.2    2020-03-02 [1] CRAN (R 4.0.2)
 tibble      * 3.0.3    2020-07-10 [1] CRAN (R 4.0.2)
 tidyr       * 1.1.1    2020-07-31 [1] CRAN (R 4.0.2)
 tidyselect    1.1.0    2020-05-11 [1] CRAN (R 4.0.2)
 tidyverse   * 1.3.0    2019-11-21 [1] CRAN (R 4.0.2)
 tinytex       0.25     2020-07-24 [1] CRAN (R 4.0.2)
 usethis     * 1.6.1    2020-04-29 [1] CRAN (R 4.0.2)
 utf8          1.1.4    2018-05-24 [1] CRAN (R 4.0.2)
 vctrs         0.3.2    2020-07-15 [1] CRAN (R 4.0.2)
 withr         2.2.0    2020-04-20 [1] CRAN (R 4.0.2)
 xfun          0.16     2020-07-24 [1] CRAN (R 4.0.2)
 xml2          1.3.2    2020-04-23 [1] CRAN (R 4.0.2)
 yaml          2.2.1    2020-02-01 [1] CRAN (R 4.0.2)

[1] C:/Users/Jane.Sullivan/Work/R/win-library/4.0
[2] C:/Program Files/R/R-4.0.2/library
```
