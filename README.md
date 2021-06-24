
# Panel Informativity Optimizer (PIO)
 
PIO is an open-source software created to help researchers and biologists in designing their NGS panels.

Based on individual mutational profiles available on public cancer genomics databases or in private databases, PIO proposes a tool to either select the minimal set of genes or exons to achieve absolute informativity (100% of the patients with at least one mutation) or assess the informativity of a given panel in a cancer of interest. Several others features are proposed by PIO, such as the exploration of mutation data frequency and distribution. By combining metrics on informativity and panel length, PIO can help making accurate choices in the design of panels, and might also be used to easily benchmark available commercial solutions. 
 
 
# Getting started

## Online version

PIO has a ready-to-use online version available at
<https://vincentalcazer.shinyapps.io/Panel_informativity_optimizer/>.

## Local version

You can install the development version from
[GitHub](https://github.com/VincentAlcazer/PIO) either by cloning the
repository or directly downloading the package in R:

``` r

install.packages("remotes")
remotes::install_github("VincentAlcazer/PIO")

PIO::run_app()
```

# Quick tutorial

## Parameters

Parameters can be selected in the left panel.

### 1. Datasets

Preloaded datasets from 91 independent cohorts spanning 31 different cancer types are available in the base package. A custom mutation dataset can be uploaded by selecting browse and the correct file format. Custom datasets should contain at least the 10 following columns:  "patient_id", "gene_id", "exon_id", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "HGVSp", "cohort".  Column names can differ. Columns can be left empty in case the variable is not available (e.g. for a mutation dataset where exons have not been annotated, exon_id will be an empty column). Example of datasets from cbioportal are provided in the raw_mutation_data folder.

### 2. Analysis parameters

#### Analysis mode

For the optimal panel selection (optimal mode), the most informative gene on the overall cohort is first selected. Patients presenting this mutation are then removed from the cohort, and the most informative gene in the remaining cohort is then selected. The algorithm reiterates until all patients are removed or no patients are added. 

For the custom panel interrogation (custom mode), only genes/exons from the provided list of genomic intervals are considered. An additional list of genes/exons, established by running the algorithm on the patients without mutation on the proposed genomic intervals is provided to complete the panel if all patients are not diagnosed with the custom panel.

#### Group mutations by

Should mutations (and their respective length) be grouped by gene or exon/intron?

#### Informativity metric

The number of unique patients per kilobase (UPKB) or the number of unique patients (UP) can be used as informativity metric. 

#### Min patients/mutation

A minimal number of patients per mutation, calculated on the overall cohort, can be set in order to avoid overfitting to private mutations. 

#### Panel upload

For custom analysis: a custom panel should be uploaded. The uploaded file should contain only one column containing the list of gene or exon names.