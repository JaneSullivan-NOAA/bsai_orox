Bering Sea and Aleutian Islands Other Rockfish stock assessment
Contact info: jane.sullivan@noaa.gov

```
> session_info()
─ Session info ───────────────────────────
 setting  value
 version  R version 4.2.0 (2022-04-22 ucrt)
 os       Windows 10 x64 (build 19044)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United States.utf8
 ctype    English_United States.utf8
 tz       America/Anchorage
 date     2022-12-14
 rstudio  2022.02.3+492 Prairie Trillium (desktop)
 pandoc   NA

─ Packages ───────────────────────────────
 package     * version date (UTC) lib source
 assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.2.0)
 bit           4.0.4   2020-08-04 [1] CRAN (R 4.2.0)
 bit64         4.0.5   2020-08-30 [1] CRAN (R 4.2.0)
 blob          1.2.3   2022-04-10 [1] CRAN (R 4.2.0)
 brio          1.1.3   2021-11-30 [1] CRAN (R 4.2.0)
 cachem        1.0.6   2021-08-19 [1] CRAN (R 4.2.0)
 callr         3.7.0   2021-04-20 [1] CRAN (R 4.2.0)
 cli           3.3.0   2022-04-25 [1] CRAN (R 4.2.0)
 codetools     0.2-18  2020-11-04 [1] CRAN (R 4.2.0)
 colorspace    2.0-3   2022-02-21 [1] CRAN (R 4.2.0)
 cowplot     * 1.1.1   2020-12-30 [1] CRAN (R 4.2.0)
 crayon        1.5.1   2022-03-26 [1] CRAN (R 4.2.0)
 DBI           1.1.3   2022-06-18 [1] CRAN (R 4.2.0)
 desc          1.4.1   2022-03-06 [1] CRAN (R 4.2.0)
 devtools    * 2.4.3   2021-11-30 [1] CRAN (R 4.2.0)
 dplyr       * 1.0.9   2022-04-28 [1] CRAN (R 4.2.0)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
 fansi         1.0.3   2022-03-24 [1] CRAN (R 4.2.0)
 fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
 fs            1.5.2   2021-12-08 [1] CRAN (R 4.2.0)
 generics      0.1.2   2022-01-31 [1] CRAN (R 4.2.0)
 ggplot2     * 3.3.6   2022-05-03 [1] CRAN (R 4.2.0)
 glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
 gtable        0.3.0   2019-03-25 [1] CRAN (R 4.2.0)
 hms           1.1.1   2021-09-26 [1] CRAN (R 4.2.0)
 lattice       0.20-45 2021-09-22 [1] CRAN (R 4.2.0)
 lifecycle     1.0.1   2021-09-24 [1] CRAN (R 4.2.0)
 lubridate     1.8.0   2021-10-07 [1] CRAN (R 4.2.0)
 magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
 memoise       2.0.1   2021-11-26 [1] CRAN (R 4.2.0)
 munsell       0.5.0   2018-06-12 [1] CRAN (R 4.2.0)
 odbc          1.3.3   2021-11-30 [1] CRAN (R 4.2.1)
 pillar        1.7.0   2022-02-01 [1] CRAN (R 4.2.0)
 pkgbuild      1.3.1   2021-12-20 [1] CRAN (R 4.2.0)
 pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
 pkgload       1.2.4   2021-11-30 [1] CRAN (R 4.2.0)
 plyr          1.8.7   2022-03-24 [1] CRAN (R 4.2.0)
 prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.2.0)
 processx      3.5.3   2022-03-25 [1] CRAN (R 4.2.0)
 ps            1.7.0   2022-04-23 [1] CRAN (R 4.2.0)
 purrr         0.3.4   2020-04-17 [1] CRAN (R 4.2.0)
 R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
 raster        3.5-15  2022-01-22 [1] CRAN (R 4.2.0)
 Rcpp          1.0.8.3 2022-03-17 [1] CRAN (R 4.2.0)
 readr       * 2.1.2   2022-01-30 [1] CRAN (R 4.2.0)
 rema        * 0.1.0   2022-09-16 [1] Github (afsc-assessments/rema@9b40438)
 remotes       2.4.2   2021-11-30 [1] CRAN (R 4.2.0)
 rlang         1.0.2   2022-03-04 [1] CRAN (R 4.2.0)
 rprojroot     2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
 rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.2.0)
 scales        1.2.0   2022-04-13 [1] CRAN (R 4.2.0)
 sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
 sp            1.5-0   2022-06-05 [1] CRAN (R 4.2.0)
 terra         1.5-34  2022-06-09 [1] CRAN (R 4.2.0)
 testthat    * 3.1.4   2022-04-26 [1] CRAN (R 4.2.0)
 tibble        3.1.7   2022-05-03 [1] CRAN (R 4.2.0)
 tidyr       * 1.2.0   2022-02-01 [1] CRAN (R 4.2.0)
 tidyselect    1.1.2   2022-02-21 [1] CRAN (R 4.2.0)
 tzdb          0.3.0   2022-03-28 [1] CRAN (R 4.2.0)
 usethis     * 2.1.6   2022-05-25 [1] CRAN (R 4.2.0)
 utf8          1.2.2   2021-07-24 [1] CRAN (R 4.2.0)
 vctrs         0.4.1   2022-04-13 [1] CRAN (R 4.2.0)
 withr         2.5.0   2022-03-03 [1] CRAN (R 4.2.0)
 zoo           1.8-10  2022-04-15 [1] CRAN (R 4.2.0)
```
