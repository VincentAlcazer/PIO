
#usethis::create_package("PIO")

library(BiocManager)
options(repos = BiocManager::repositories())
library(shinydashboard)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(viridis)
library(data.table)
library(DT)
library(forcats)
library(ggrepel)
library(shinyjs)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(shinycssloaders)
library(shinyalert)


########## ========== Parameters

source("functions.R")

### Ggplot 2 default theme

custom_theme <- theme_bw() + theme(
  plot.title = element_text(size = 20, face = "bold"),
  axis.text = element_text(size = 16, color = "black"),
  axis.title = element_text(size = 18, face = "bold"),
  legend.title = element_text(size = 18, face = "bold"),
  legend.text = element_text(size = 16),
  legend.position = "none"
)



exon_length <- read_rds("references/2021_04_19_Refseq_GRCH37_stats.rds")
clin_impact <- read_rds("references/2021_12_07_CIVIC_clinical_evidences.rds")

#### Run app

run_app <- function(){shiny::runApp()}
