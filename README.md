
<!-- README.md is generated from README.Rmd. Please edit that file -->

# waterfowl

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The waterfowl package is aimed at simplifying functions and workflows
for waterfowl biologist and graduate students. Ultimately, its function
are simple wrappers around functions from other packages, but with more
automation and with a preferred output. Looking to do some habitat
selection analyses with tracking data? Simply use landcover_points() to
generate a table of land cover classes for each location, or use
landcover_percent() to generate a GLMM-ready table of land cover
percentages for a given home range or area. No need to fuss with all the
intermediate steps; the waterfowl package will do that for you!
Suggestions for additions or changes to the package are more than
welcome!

The waterfowl package is a work in progress!

## Installation

You can install the development version of waterfowl from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("weberbirding/waterfowl")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(waterfowl)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
