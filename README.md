
## asmcjr <img src="man/figures/logo.png" width="140" align="right" /> <br />

[Build
Status](https://travis-ci.com/yl17124/asmcjr.svg?branch=master)\](<https://travis-ci.com/yl17124/asmcjr>)
[![codecov](https://codecov.io/gh/davidycliao/asmcjr/branch/master/graph/badge.svg?token=OJKOF5SX9X)](https://codecov.io/gh/davidycliao/asmcjr)

This package supports the book [“2nd Edition Analyzing Spatial Models of
Choice and
Judgment”](https://www.routledge.com/Analyzing-Spatial-Models-of-Choice-and-Judgment/II-Bakker-Carroll-Hare-Poole-Rosenthal/p/book/9781138715332).
In its second edition, much of the R code has been streamlined. This
package contains all of the data and functions to replicate the analyses
in the book.

<br />
<img src="https://raw.githack.com/yl17124/asmcjr/master/vignettes/book_image.jpg" width="200" align="center" />  
 

## Installation

You will need lastest installation of
[*R*](https://cran.r-project.org/mirrors.html) (preferably version 3.6
or above) and
[RStudio](https://rstudio.com/products/rstudio/download/#download).
Visit [Installation](articles/installation.html) for further
instructions.

<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
install.packages("devtools", dependencies=TRUE)
#> 
#> The downloaded binary packages are in
#>  /var/folders/60/xlhp81vn4_79cp2dbjgy07tm0000gn/T//RtmpwFTSnk/downloaded_packages
library(devtools)
#> Loading required package: usethis
#> Warning: package 'usethis' was built under R version 4.3.3
devtools::install_github("uniofessex/asmcjr")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo uniofessex/asmcjr@HEAD
#> Skipping 1 packages not available: basicspace
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/60/xlhp81vn4_79cp2dbjgy07tm0000gn/T/RtmpwFTSnk/remotes916a750e2b86/uniofessex-asmcjr-8d7d9be/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/60/xlhp81vn4_79cp2dbjgy07tm0000gn/T/RtmpwFTSnk/remotes916a750e2b86/uniofessex-asmcjr-8d7d9be/DESCRIPTION’
#>   ─  preparing ‘asmcjr’:
#>      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  cleaning src
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>        NB: this package now depends on R (>= 3.5.0)
#>        WARNING: Added dependency on R >= 3.5.0 because serialized objects in
#>      serialize/load version 3 cannot be read in older versions of R.
#>      File(s) containing such objects:
#>        ‘asmcjr/data/ANES1968.rda’ ‘asmcjr/data/ANES2004.rda’
#>        ‘asmcjr/data/ANES2004_OOC.rda’ ‘asmcjr/data/CDS2000.rda’
#>        ‘asmcjr/data/SOTUcorpus.rda’ ‘asmcjr/data/candidatetherms2008.rda’
#>        ‘asmcjr/data/denmarkEES2009.rda’ ‘asmcjr/data/france4.rda’
#>        ‘asmcjr/data/hr108.rda’ ‘asmcjr/data/hr111.rda’
#>        ‘asmcjr/data/rc_ep.rda’ ‘asmcjr/data/rcv_ep1.rda’
#>        ‘asmcjr/data/rcx.rda’ ‘asmcjr/data/uk.rda’
#>        ‘asmcjr/vignettes/legis_7th_Taiwan.rda’
#>   ─  building ‘asmcjr_1.0.2.tar.gz’
#>      
#> 
```
