

usethis::create_package("./")

usethis::use_package_doc()

usethis::use_mit_license(copyright_holder = "Vincent Alcazer")

lapply(c("BiocManager","shinydashboard","dplyr","readr","ggplot2",
            "viridis","data.table","DT","forcats","ggrepel","shinyjs","limma",
            "AnnotationDbi","org.Hs.eg.db"), usethis::use_package)


roxygen2::roxygenise()

