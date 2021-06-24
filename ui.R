
ui <- dashboardPage(
  dashboardHeader(title = "PIO v1.0"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      # menuItem("Introduction/Guide", tabName = "introduction", icon = icon("play-circle")),
      # menuItem("Panel assistant", tabName = "panel_assistant"),
      h3("1. Datasets", style = "margin-bottom: 0;"),
      selectInput("disease",
        label = "Preloaded datasets",
        choices = c(
          "ACC", "ALL", "AML", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
          "DLBCL", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC",
          "LUAD", "LUSC", "OV", "PAAD", "PRAD", "SARC", "SKCM", "STAD", "TGCT",
          "THCA", "THYM", "UCEC", "UCS", "UVM"
        )
      ),
      fileInput("custom_data",
        label = "Custom dataset",
        # label = h4("3. Load your panel (for custom analysis only)"),
        accept = c(
          "text/tab-separated-values",
          "text/comma-separated-values",
          "text/plain",
          "text/csv",
          ".csv",
          ".tsv"
        )
      ),
      div(style = "margin-top: -50px"),
      radioButtons("sep_custom", "Dataset file type (delim)",
        choices = c(
          "comma-delim (.csv1)" = ",",
          "semi colon-delim (.csv2)" = ";",
          "tab-delim (.tsv/.txt)" = "\t"
        ),
        selected = "\t"
      ),
      # actionButton("load_data", "Load dataset"),

      h3("2. Analysis parameters", style = "margin-bottom: 0;"),
      radioButtons("analysis_type",
        label = ("Analysis mode"),
        choices = c(
          "PIO Optimal" = "optimal", "PIO Custom" = "custom",
          "Panel Test" = "test"
        ),
        selected = "optimal"
      ),
      div(style = "margin-top: -20px"),
      radioButtons("mutation_group", "Group mutations by",
        choices = c(
          "Exon/intron" = "exon_id",
          "Gene" = "gene_id"
        ), selected = "gene_id"
      ),
      div(style = "margin-top: -20px"),
      radioButtons("metric_type",
        label = "Informativity metric",
        choices = c(
          "Unique patients (UP)" = "unique_pts",
          "Unique patients / kb (UPKB)" = "pts_per_kb"
        ), selected = "unique_pts"
      ),
      div(style = "margin-top: -20px"),
      numericInput("min_pts",
        label = "Min. patients / mutation",
        value = 2, min = 1
      ),

      # div(style = "margin-top: -25px"),
      fileInput("panel_list",
        label = "Panel upload (custom analysis only)",
        # label = h4("3. Load your panel (for custom analysis only)"),
        accept = c(
          "text/tab-separated-values",
          "text/comma-separated-values",
          "text/plain",
          "text/csv",
          ".csv",
          ".tsv"
        )
      ),
      div(style = "margin-top: -40px"),
      radioButtons("sep", "Panel file type (delim)",
        choices = c(
          "comma-delim (.csv1)" = ",",
          "semi colon-delim (.csv2)" = ";",
          "tab-delim (.tsv/.txt)" = "\t"
        ),
        selected = "\t"
      ),
      actionButton("run_analysis", "Run analysis")
    ) # sidebar Menu
  ), # Side bar
  dashboardBody(
    fluidPage(
      useShinyjs(),
      h2("Panel Informativity Optimizer v1.0"),
      column(
        9,
        tabsetPanel(
          id = "panassis", type = "tabs",
          tabPanel(
            "Informativity",
            h3("Methods"),
            column(
              12,
              p("The optimal set of mutated genes for disease diagnosis is proposed as follow:
                                 The most frequent mutation in the given cohort is first selected.
                                 Patients presenting this mutation are removed from the cohort, and
                                 the most frequent mutation in the remaining patients is then selected.
                                 The process is reiterated untill no more patients are added or the total number
                                 of patients in the cohort is reached.
                                 The cumulative percentage of diagnosed patient is then represented
                                 on yellow/blue color scale, according to the number of independent cohorts
                                 the gene was found mutated in.
                                 The cumulative exonic length of the panel
                                 (total exonic length for genes analysis, exonic length for exons)
                                 is represented in red as a percent of the total panel length.


                                   ")
            ),
            h3("Merged cohort analysis"),
            column(
              12,
              downloadButton("download_panel", "Download merged cohort table (.tsv)"),
              shinycssloaders::withSpinner(plotOutput("graph_informative_merged"), type = 6)
            ),
            br(),
            h3("Suggested genes to add"),
            column(
              12,
              downloadButton("download_sug_genes", "Download suggested genes (.tsv)"),
              textOutput("sup_genes_message"),
              DT::DTOutput("sug_genes")
            ),
            br(),
            h3("Individual cohort analysis"),
            column(
              12,
              br(),
              downloadButton("download_panel_indiv", "Download individual cohorts table (.tsv)"),
              shinycssloaders::withSpinner(plotOutput("graph_informative"), type = 6),
            )
          ), # tabpanel
          tabPanel(
            "Mutations",
            tabsetPanel(
              id = "mutations", type = "tabs",
              tabPanel(
                "Statistics",
                column(
                  12,
                  h3("Mutations stats (overall cohort)"),
                  p("Overall mutations statistics are computed on the overall merged cohort.
                                                        Mutation are grouped either by genes
                                                        (length = total exonic length)
                                                        or exons (length = exon length).
                                                        NB: The Min. patients / mutation filter is not applied here.
                                                        "),
                  downloadButton("download_mut_stat", "Download table (.tsv)"),
                  shinycssloaders::withSpinner(DT::DTOutput("table_mut_stat"))
                ),
                column(
                  12, h3("Cumulated mutations"),
                  downloadButton("download_table_des", "Download table (.tsv)"),
                  p("Cumulated mutations statistics are reported for mutations
                                               selected in the panel (informativity tab).
                                               The number of patients with at least or exactly n mutations is reported.

                                                        "),
                  shinycssloaders::withSpinner(DT::DTOutput("table_des"))
                )
              ), # tabpanel
              tabPanel(
                "Types & frequencies",
                column(
                  12,
                  h3("Mutations types & frequencies"),
                  p("Overall mutations frequencies are computed on the overall merged cohort.
                                                        NB: The Min. patients / mutation filter is not applied here.
                                                        "),
                  downloadButton("download_mutfreq", "Download table (.tsv)"),
                  column(10), column(2),
                  shinycssloaders::withSpinner(plotOutput("graph_mut_freq"), type = 6)
                )
              ),
              tabPanel(
                "Distribution",
                column(
                  12,
                  h3("Mutations distribution"),
                  p("The mutation status of all the mutations
                                                      retained in the panel (informativity tab)
                                                      is represented across all patients from the
                                                        merged cohorts"),
                  br(),
                  shinycssloaders::withSpinner(plotOutput("graph_heat"), type = 6)
                ) # column
              )
            ) # tabsetpanel
          ), # tabpanel
          tabPanel(
            "Length analysis",
            column(
              12,
              br(),
              p("Mutations are grouped and counted by gene or exon/intron for the overall cohort (merged datasets).
                                   Four approaches are compared to select the most informative panel
                                   within the size limit: two classic approaches where genes or exons/introns are
                                   ranked according to the corresponding metric on the overall cohort and then
                                   linearly added to the panel untill reaching the size limit, and
                                   two approaches using PIO algorithm with the corresponding metrics.
                                    "),
              numericInput("max_length",
                label = "Max panel length (kb)",
                min = 1, value = 1000
              ), actionButton("run_length_analysis", "Run length analysis"),
              downloadButton("download_length", "Download full table (.tsv)"),
              shinycssloaders::withSpinner(plotOutput("graph_length"), type = 6)
            ) # column
          ) # tabpanel
        ) # tabset panel
      ), # column
      column(
        3,
        wellPanel(
          h4("Graph Parameters"),
          numericInput("max_freq", label = "Max frequency", value = 1, min = 0, max = 1),
          numericInput("max_genes", label = "Max genes/exons", value = 100, min = 1),
          sliderInput("x_size", label = "Genes/Exons fontsize", value = 12, min = 1, max = 30),
          # radioButtons("freq_type", label="Frequencies:",
          #              choices = c("Count"="cum_sum",
          #                          "Percent"="percent_tot"),
          #              selected = "cum_sum"),
          actionButton("apply_param", "Apply")
        )
      ) # column
    ) # fluid Page
  ) # body
) # Page
