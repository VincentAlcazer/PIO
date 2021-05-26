
# Panel Informativity Optimizer (PIO)
 
PIO is an open-source software created to help researchers and biologists in designing their NGS panels.

Based on individual mutational profiles available on public cancer genomics databases or in private databases, PIO proposes a tool to either select the minimal set of genes or exons to achieve absolute informativity (100% of the patients with at least one mutation) or assess the informativity of a given panel in a cancer of interest. Several others features are proposed by PIO, such as the exploration of mutation data frequency and distribution. By combining metrics on informativity and panel length, PIO can help making accurate choices in the design of panels, and might also be used to easily benchmark available commercial solutions. 
 
 
# Getting started

## Online version

PIO has a ready-to-use online version available at
<https://vincentalcazer.shinyapps.io/Panel_informativity_optimizer/>.

## Local version

You can install the development version from
[GitHub](https://github.com/VincentAlcazer/PIO) else by cloning the
repository or directly by downloading the package in R:

``` r

# install.packages("remotes")
# remotes::install_github("VincentAlcazer/PIO")

# PIO::run_app()
```