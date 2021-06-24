
server <- function(input, output) {

  ########## ========== DATAFRAME: Loading

  ##### ===== Load raw datasets

  mutations_raw_data <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)

    isolate({
      if (is.null(input$custom_data) == F) {

        ### Custom dataset read

        df <- data.table::fread(input$custom_data$datapath, sep = input$sep_custom) %>%
          as.data.frame() %>%
          dplyr::select(1:10)
        colnames(df) <- c(
          "patient_id", "gene_id", "exon_id", "Chromosome", "Start_Position", "End_Position",
          "Strand", "Variant_Classification", "HGVSp", "cohort"
        )

        ## Update symbols names
        df$gene_id <- limma::alias2SymbolTable(df$gene_id)

        if (all(is.na(df$exon_id))) {
          df <- df %>%
            dplyr::select(-exon_id) %>%
            left_join(exon_length, by = "gene_id") %>%
            mutate(diff = Start_Position - start) %>%
            arrange(abs(diff)) %>%
            distinct(gene_id, patient_id, Variant_Classification, Start_Position, End_Position, .keep_all = T) %>%
            mutate(outlier = abs(diff) > 10000 & abs(diff) > tot_gene_length) %>%
            filter(outlier == F)
        } else {
          df <- df %>%
            left_join(exon_length, by = c("gene_id", "exon_id"))
        }

        return(df)
      } else {
        files <- list.files(paste0("raw_mutation_data/", input$disease))

        file_list <- list()

        for (f in files) {
          name <- gsub(".rds", "", f)

          file_list[[name]] <- read_rds(paste0("raw_mutation_data/", input$disease, "/", f)) %>%
            mutate(cohort = name)
        }

        df <- bind_rows(file_list)
      }



      return(df)
    })
  })

  ### Custom panel df & files

  panel_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)
    isolate({
      df <- data.table::fread(input$panel_list$datapath, sep = input$sep) %>% as.data.frame()
      colnames(df)[1] <- "mutation_id"

      return(df)
    })
  })

  ########## ========== Get informative mutations data

  ### Merged dataset
  informative_merged_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)

    isolate({
      merged_dataset <- list(merged = mutations_raw_data())

      if (is.null(input$panel_list) | input$analysis_type == "optimal") {
        custom_panel <- NULL
      } else {
        custom_panel <- unlist(panel_df() %>% dplyr::select(mutation_id)) %>% as.character()
      }

      if (input$analysis_type == "test") {
        df <- panel_tester(merged_dataset,
          group_by = input$mutation_group,
          info_metric = input$metric_type, min_pts = input$min_pts,
          custom_panel = custom_panel, keep_order = F
        )
      } else {
        df <- panel_finder(merged_dataset,
          group_by = input$mutation_group,
          info_metric = input$metric_type, min_pts = input$min_pts,
          custom_panel = custom_panel
        )
      }


      ## Add the number of independent datasets
      independent_sets <- mutations_raw_data() %>%
        dplyr::select(mutation_id = input$mutation_group, everything()) %>%
        group_by(mutation_id) %>%
        summarise(n_cohort = length(unique(cohort)))

      df <- df %>%
        left_join(independent_sets, by = "mutation_id") %>%
        droplevels()

      return(df)
    }) # isolate
  })

  ### Table: suggested genes to reach full informativity

  suggested_genes <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)

    isolate({
      if (input$analysis_type == "optimal" | max(informative_merged_df() %>% dplyr::select(percent_tot)) >= 99.9) {
        NULL
      } else {
        optimal_panel <- informative_merged_df() %>%
          distinct(mutation_id) %>%
          unlist() %>%
          as.character()

        diag_pts <- mutations_raw_data() %>%
          filter(gene_id %in% optimal_panel | exon_id %in% optimal_panel) %>%
          distinct(patient_id) %>%
          unlist() %>%
          as.character()

        tot_pts <- nrow(mutations_raw_data() %>% distinct(patient_id))

        new_data_list <- list(merged = (mutations_raw_data() %>%
          filter(!patient_id %in% diag_pts)))

        # if(is.null(input$panel_list) | input$analysis_type == "optimal"){
        #   custom_panel = NULL
        # } else {
        #   custom_panel <- unlist(panel_df() %>% dplyr::select(mutation_id)) %>% as.character()
        # }
        #


        df <- panel_finder(new_data_list,
          group_by = input$mutation_group,
          info_metric = input$metric_type, min_pts = input$min_pts,
          custom_panel = NULL
        ) %>%
          mutate(
            percent_remaining = round(percent_tot, 2),
            percent_tot = round(n_pts / tot_pts, 3)
          ) %>%
          dplyr::select(mutation_id, n_pts, cum_sum, percent_remaining, percent_tot, everything())


        return(df)
      }
    })
  })

  output$sug_genes <- DT::renderDT(

    suggested_genes(), # data
    class = "display nowrap compact", # style
    filter = "none", # location of column filters
    server = T,
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      lengthChange = TRUE,
      sDom = '<"top">lrt<"bottom">ip',
      columnDefs = list(list(className = "dt-left", targets = "_all"))
    )
  )

  # sup_genes_message <- reactive({
  #   validate(!is.null(suggested_genes()), "Full informativity reached: no genes to add")
  # })

  output$sup_genes_message <- renderText({
    optimal_panel <- informative_merged_df() %>%
      distinct(mutation_id) %>%
      unlist() %>%
      as.character()

    diag_pts <- mutations_raw_data() %>%
      filter(gene_id %in% optimal_panel | exon_id %in% optimal_panel) %>%
      distinct(patient_id) %>%
      unlist() %>%
      as.character()

    remaining_pts <- mutations_raw_data() %>%
      filter(!patient_id %in% diag_pts) %>%
      distinct(patient_id) %>%
      nrow() %>%
      as.numeric()


    if (is.null(suggested_genes())) {
      ("Full informativity reached or optimal mode selected: no genes to add")
    } else {
      paste0("Suggested genes based on the ", remaining_pts, " remaining patients:")
    }
  })

  ### Mutation stats
  table_mut_stat_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)
    isolate({
      if (input$mutation_group == "exon_id") {
        feature_length <- "exon_length"
      } else {
        feature_length <- "tot_exons_length"
      }

      if (is.null(input$panel_list) | input$analysis_type == "optimal") {
        custom_panel <- mutations_raw_data() %>%
          dplyr::select(mutation_id = input$mutation_group) %>%
          unlist() %>%
          as.character()
      } else {
        custom_panel <- unlist(panel_df() %>% dplyr::select(mutation_id)) %>% as.character()
      }


      df <- mutations_raw_data() %>%
        mutate(n_pts = length(unique(patient_id))) %>%
        dplyr::select(mutation_id = input$mutation_group, length = feature_length, everything()) %>%
        filter(mutation_id %in% custom_panel) %>%
        group_by(mutation_id, length) %>%
        summarise(
          unique_pts = length(unique(patient_id)),
          unique_mut = length(unique(Start_Position)),
          mut_types = length(unique(Variant_Classification)),
          cohort_pts = unique(n_pts),
          n_cohort = length(unique(cohort))
        ) %>%
        arrange(desc(unique_pts)) %>%
        mutate(
          pts_per_kb = round(unique_pts * 1000 / length, 4),
          percent_total = round(unique_pts / cohort_pts, 3)
        ) %>%
        ungroup() %>%
        dplyr::select(
          mutation_id, unique_pts, percent_total, unique_mut, mut_types,
          length, pts_per_kb, n_cohort, cohort_pts, everything()
        )

      if (input$mutation_group == "exon_id") {
        df <- df %>%
          left_join(exon_length, by = c("mutation_id" = "exon_id", "length" = "exon_length")) %>%
          dplyr::select(-gene_id, -tot_exons_length, -tot_gene_length)
      }

      return(df)
    })
  })

  output$table_mut_stat <- DT::renderDT(
    table_mut_stat_df(), # data
    class = "display nowrap compact", # style
    filter = "none", # location of column filters
    server = T,
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      lengthChange = TRUE,
      sDom = '<"top">lrt<"bottom">ip',
      columnDefs = list(list(className = "dt-left", targets = "_all"))
    )
  )

  ### Table des
  table_des_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)
    isolate({
      mutations <- informative_merged_df() %>%
        dplyr::select(mutation_id) %>%
        unlist() %>%
        as.character()

      df <- mutations_raw_data() %>%
        dplyr::select(mutation_id = input$mutation_group, everything())

      df$mutation_id[!df$mutation_id %in% mutations] <- NA

      df_stat <- df %>%
        filter(is.na(mutation_id) == F) %>%
        group_by(patient_id) %>%
        summarise(n_unique_mut = length(unique(mutation_id))) %>%
        group_by(n_mutations = n_unique_mut) %>%
        summarise(exact = n()) %>%
        arrange(desc(n_mutations)) %>%
        mutate(at_least = cumsum(exact)) %>%
        rbind(data.frame(
          n_mutations = 0,
          exact = length(setdiff(
            unique(df$patient_id),
            unique(df$patient_id[df$mutation_id %in% mutations])
          )),
          at_least = length(setdiff(
            unique(df$patient_id),
            unique(df$patient_id[df$mutation_id %in% mutations])
          ))
        )) %>%
        arrange((n_mutations))

      return(df_stat)
    })
  })


  output$table_des <- DT::renderDT(
    table_des_df(), # data
    class = "display nowrap compact", # style
    filter = "none", # location of column filters
    server = T,
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      lengthChange = TRUE,
      sDom = '<"top">lrt<"bottom">ip',
      columnDefs = list(list(className = "dt-left", targets = "_all"))
    )
  )


  ### Individuals datasets

  informative_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)

    isolate({
      dataset_list <- split(mutations_raw_data(), mutations_raw_data()[, "cohort"])

      if (is.null(input$panel_list) | input$analysis_type == "optimal") {
        custom_panel <- NULL
      } else {
        custom_panel <- unlist(panel_df() %>% dplyr::select(mutation_id)) %>% as.character()
      }



      if (length(names(dataset_list)) == 1) {
        df <- informative_merged_df() %>% mutate(cohort = names(dataset_list))
      } else {
        if (input$analysis_type == "test") {
          df <- panel_tester(dataset_list,
            group_by = input$mutation_group,
            info_metric = input$metric_type, min_pts = input$min_pts,
            custom_panel = custom_panel, keep_order = F
          )
        } else {
          df <- panel_finder(dataset_list,
            group_by = input$mutation_group,
            info_metric = input$metric_type, min_pts = input$min_pts,
            custom_panel = custom_panel
          )
        }
      }




      return(df)
    }) # isolate
  })


  ########## ========== Mutations frequencies & heat data

  ### Merged dataset
  mutation_des_merged_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)

    isolate({
      if (input$mutation_group == "exon_id") {
        feature_length <- "exon_length"
      } else {
        feature_length <- "tot_exons_length"
      }

      if (is.null(input$panel_list) | input$analysis_type == "optimal") {
        custom_panel <- mutations_raw_data() %>%
          dplyr::select(mutation_id = input$mutation_group) %>%
          unlist() %>%
          as.character()
      } else {
        custom_panel <- unlist(panel_df() %>% dplyr::select(mutation_id)) %>% as.character()
      }

      df <- mutations_raw_data() %>%
        mutate(n_pts = length(unique(patient_id))) %>%
        dplyr::select(mutation_id = input$mutation_group, length = feature_length, everything()) %>%
        filter(mutation_id %in% custom_panel) %>%
        mutate(n_pts = length(unique(patient_id))) %>%
        group_by(mutation_id, length, Chromosome, Start_Position, End_Position, Strand, Variant_Classification, HGVSp) %>%
        summarise(
          count = length(unique(patient_id)),
          n_pts = unique(n_pts),
          n_cohort = length(unique(cohort))
        ) %>%
        mutate(percent_tot = count / n_pts) %>%
        ungroup() %>%
        arrange(desc(percent_tot))

      return(df)
    })
  })

  ### Heatmap dataframe
  mutation_heat_df <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)

    isolate({
      optimal_genes <- informative_merged_df() %>%
        distinct(mutation_id) %>%
        unlist() %>%
        as.character()

      df <- mutations_raw_data() %>%
        dplyr::select(
          mutation_id = input$mutation_group,
          patient_id, Variant_Classification
        ) %>%
        mutate(optimal_mutation = if_else(mutation_id %in% optimal_genes, TRUE, FALSE))

      gene_order <- df %>%
        filter(mutation_id %in% optimal_genes) %>%
        group_by(mutation_id) %>%
        summarise(count = n()) %>%
        arrange(desc(count)) %>%
        distinct(mutation_id) %>%
        unlist() %>%
        as.character()

      no_mut_patients <- df %>%
        filter(optimal_mutation == F) %>%
        distinct(patient_id) %>%
        unlist() %>%
        as.character()
      mut_patients <- df %>%
        filter(optimal_mutation == T) %>%
        distinct(patient_id) %>%
        unlist() %>%
        as.character()

      pts_to_add <- data.frame(expand.grid(optimal_genes, setdiff(no_mut_patients, mut_patients))) %>%
        mutate(mutation = "0")
      colnames(pts_to_add) <- c("mutation_id", "patient_id", "mutation")

      heat_df <- df %>%
        filter(optimal_mutation == T) %>%
        dplyr::select(mutation_id, patient_id) %>%
        mutate(mutation = "1") %>%
        mutate(mutation_id = factor(mutation_id, levels = gene_order)) %>%
        arrange(mutation_id) %>%
        rbind(pts_to_add)

      patient_order <- heat_df$patient_id

      heat_df$patient_id <- factor(heat_df$patient_id, levels = unique(patient_order))

      return(heat_df)
    })
  })

  ########## ========== Length analysis
  graph_length_df <- reactive({
    input$run_length_analysis
    req(input$run_length_analysis >= 1)
    isolate({
      data_list <- list(merged = mutations_raw_data())

      ### Base approach: n pts
      top_genes <- table_mut_stat_df() %>%
        filter(unique_pts >= input$min_pts) %>%
        mutate(cum_length = cumsum(length)) %>%
        filter(cum_length <= 1000 * input$max_length) %>%
        distinct(mutation_id) %>%
        unlist() %>%
        as.character()

      base_pts <- panel_tester(data_list,
        group_by = input$mutation_group,
        info_metric = "unique_pts", min_pts = input$min_pts,
        custom_panel = top_genes
      ) %>%
        mutate(cum_length = cumsum(length), group = "unique_pts") %>%
        filter(cum_length <= 1000 * input$max_length)


      ### Base approach: pts per kb
      top_genes <- table_mut_stat_df() %>%
        filter(unique_pts >= input$min_pts) %>%
        arrange(desc(pts_per_kb)) %>%
        mutate(cum_length = cumsum(length)) %>%
        filter(cum_length <= 1000 * input$max_length) %>%
        distinct(mutation_id) %>%
        unlist() %>%
        as.character()

      base_pts_kb <- panel_tester(data_list,
        group_by = input$mutation_group,
        info_metric = "pts_per_kb", min_pts = input$min_pts,
        custom_panel = top_genes
      ) %>%
        mutate(cum_length = cumsum(length), group = "pts_per_kb") %>%
        filter(cum_length <= 1000 * input$max_length)

      ### PIO: n pts
      PIO_pts <- panel_finder(data_list,
        group_by = input$mutation_group,
        info_metric = "unique_pts", min_pts = input$min_pts,
        custom_panel = NULL
      ) %>%
        mutate(cum_length = cumsum(length), group = "PIO_unique_pts") %>%
        filter(cum_length <= 1000 * input$max_length)


      ### PIO: pts per kb
      PIO_pts_per_kb <- panel_finder(data_list,
        group_by = input$mutation_group,
        info_metric = "pts_per_kb", min_pts = input$min_pts,
        custom_panel = NULL
      ) %>%
        mutate(cum_length = cumsum(length), group = "PIO_pts_per_kb") %>%
        filter(cum_length <= 1000 * input$max_length)


      df <- rbind(base_pts, base_pts_kb, PIO_pts, PIO_pts_per_kb) %>%
        mutate(group = factor(group, levels = c(
          "unique_pts", "pts_per_kb",
          "PIO_unique_pts", "PIO_pts_per_kb"
        )))

      return(df)
    }) # Isolate
  }) # Renderplot

  ########## ========== GRAPHS

  ### === Custom graph parameters

  graph_height <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)
    isolate({
      n_cohorts <- length(unique(informative_df() %>% dplyr::select(cohort) %>% unlist()))

      if (n_cohorts <= 2) {
        height <- 400
      } else if (n_cohorts <= 4) {
        height <- 800
      } else {
        height <- 1200
      }
      return(height)
    })
  })

  graph_height_2 <- reactive({
    input$run_analysis
    req(input$run_analysis >= 1)
    isolate({
      n_genes <- length(unique(informative_df() %>% distinct(mutation_id) %>% unlist() %>% as.character()))

      if (n_genes <= 25) {
        height <- 400
      } else if (n_genes <= 50) {
        height <- 800
      } else {
        height <- 1200
      }
      return(height)
    })
  })

  ### ===Panel graphs

  output$graph_informative_merged <- renderPlot(
    {
      input$run_analysis
      input$apply_param
      req(input$run_analysis >= 1)

      isolate({
        plot_genes <- unique(informative_merged_df() %>% dplyr::select(mutation_id) %>% unlist())

        if (length(plot_genes) > input$max_genes) {
          plot_genes <- plot_genes[1:input$max_genes]
        }


        data_plot <- informative_merged_df() %>%
          filter(
            percent_tot <= input$max_freq,
            mutation_id %in% plot_genes
          ) %>%
          droplevels()



        data_anno <- data.frame(
          x = nth(data_plot$mutation_id, 0.85 * length(data_plot$mutation_id)),
          y = 0.3,
          tot_pts = paste0(
            round(max(data_plot %>% dplyr::select(percent_tot)) * 100, 2), "% pts diagnosed \n with ",
            length(unique(data_plot %>% dplyr::select(mutation_id) %>% unlist())), " ",
            if_else(input$mutation_group == "gene_id", "genes", "exons/introns")
          )
        )


        n_pts <- mutations_raw_data() %>%
          distinct(patient_id) %>%
          unlist() %>%
          length()


        #  test <- panel_data_merged %>% filter(Disease == "AML")

        data_length <- data_plot %>%
          mutate(panel_length = cumsum(length)) %>%
          mutate(percent_length = panel_length / max(panel_length))

        ggplot(data = data_plot, aes(
          x = forcats::fct_reorder(mutation_id, percent_tot),
          y = percent_tot, group = 1
        )) +
          geom_point(aes(color = as.factor(n_cohort))) +
          geom_line(aes(color = as.factor(n_cohort))) +
          scale_color_viridis_d() +
          geom_line(
            data = filter(data_length, mutation_id %in% plot_genes),
            aes(
              x = forcats::fct_reorder(mutation_id, percent_tot),
              y = percent_length, group = 1
            ), color = "red"
          ) +
          geom_text_repel(data = data_anno, aes(x = x, y = y, label = tot_pts), size = 6) +
          custom_theme +
          labs(
            title = paste0(
              "Most informative genes in ", if_else(is.null(input$custom_data), input$disease, "Custom dataset"),
              " (n= ", n_pts, " pts)"
            ),
            x = "", group = "test",
            color = "Cohorts (n)", caption = "Cumulated % of total panel exonic length"
          ) +
          ylim(0, 1) +
          theme(
            axis.text.x = element_text(
              angle = 90, vjust = 0.5, hjust = 1, face = "bold",
              size = input$x_size
            ),
            plot.caption = element_text(color = "red", size = 10),
            legend.position = "right"
          )
      })
    },
    height = 400
  )

  output$graph_informative <- renderPlot(
    {
      input$run_analysis
      input$apply_param
      req(input$run_analysis >= 1)
      isolate({
        plot_list <- list()

        for (study in unique(informative_df() %>% dplyr::select(cohort) %>% unlist())) {
          data_plot <- informative_df() %>%
            filter(cohort == study) %>%
            droplevels()

          plot_genes <- unique(data_plot %>% dplyr::select(mutation_id) %>% unlist())

          if (length(plot_genes) > input$max_genes) {
            plot_genes <- plot_genes[1:input$max_genes]
          }

          data_plot <- data_plot %>%
            filter(
              percent_tot <= input$max_freq,
              mutation_id %in% plot_genes
            ) %>%
            droplevels()


          data_anno <- data.frame(
            x = nth(data_plot$mutation_id, 0.85 * length(data_plot$mutation_id)),
            y = 0.3,
            tot_pts = paste0(
              round(max(data_plot %>% dplyr::select(percent_tot)) * 100, 2), "% pts diagnosed \n with ",
              length(unique(data_plot %>% dplyr::select(mutation_id) %>% unlist())),
              if_else(input$mutation_group == "gene_id", " genes", " exons/introns")
            )
          )


          n_pts <- mutations_raw_data() %>%
            filter(cohort == study) %>%
            distinct(patient_id) %>%
            unlist() %>%
            length()


          #  test <- panel_data_merged %>% filter(Disease == "AML")

          data_length <- data_plot %>%
            mutate(panel_length = cumsum(length)) %>%
            mutate(percent_length = panel_length / max(panel_length))



          plot_list[[study]] <- ggplot(data = data_plot, aes(
            x = forcats::fct_reorder(mutation_id, percent_tot),
            y = percent_tot, group = 1
          )) +
            geom_point() +
            geom_line() +
            scale_color_viridis_d() +
            geom_line(
              data = filter(data_length, mutation_id %in% plot_genes),
              aes(
                x = forcats::fct_reorder(mutation_id, percent_tot),
                y = percent_length, group = 1
              ), color = "red"
            ) +
            geom_text_repel(data = data_anno, aes(x = x, y = y, label = tot_pts), size = 6) +
            custom_theme +
            labs(
              title = paste0(study, " (n= ", n_pts, " pts)"),
              color = "Cohorts (n)", x = ""
            ) +
            ylim(0, 1) +
            theme(
              axis.text.x = element_text(
                angle = 90, vjust = 0.5, hjust = 1, face = "bold",
                size = input$x_size
              ),
              plot.caption = element_text(color = "red", size = 10),
              legend.position = "none"
            )
        }

        gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
      }) # Isolate
    },
    height = graph_height
  ) # Renderplot


  ### === Mutations frequencies

  output$graph_mut_freq <- renderPlot(
    {
      input$apply_param
      input$run_analysis
      req(input$run_analysis >= 1)
      isolate({
        top_n_genes <- mutation_des_merged_df() %>%
          group_by(mutation_id) %>%
          summarise(
            tot_pts = sum(percent_tot),
            tot_count = sum(count)
          ) %>%
          arrange(desc(tot_pts)) %>%
          dplyr::slice(1:input$max_genes) %>%
          distinct(mutation_id) %>%
          unlist() %>%
          as.character()


        data_plot <- mutation_des_merged_df() %>%
          filter(mutation_id %in% top_n_genes) %>%
          mutate(mutation_id = factor(mutation_id, levels = top_n_genes)) %>%
          group_by(mutation_id, Variant_Classification) %>%
          summarise(
            count = n(),
            percent_tot = sum(percent_tot)
          )

        title <- if_else(is.null(input$custom_data) == F, "custom dataset", input$disease)

        data_plot %>%
          ggplot(aes_string(x = "mutation_id", y = "percent_tot", group = 1, fill = "Variant_Classification")) +
          geom_bar(stat = "identity", color = "black") +
          custom_theme +
          # scale_color_viridis_c()+
          labs(
            title = paste0(
              "Most frequent mutations in ", title,
              " (n= ", sum(informative_merged_df() %>% dplyr::select(n_pts)), " pts)"
            ),
            y = "Percent patients", x = "",
            fill = "Mutation type"
          ) +
          theme(
            axis.text.x = element_text(
              angle = 90, vjust = 0.5, hjust = 1,
              size = input$x_size, face = "bold"
            ),
            legend.position = "bottom"
          ) +
          guides(fill = guide_legend(nrow = 6))
      }) # Isolate
    },
    height = 800
  ) # Renderplot

  ### === Heatmap

  output$graph_heat <- renderPlot(
    {
      input$run_analysis
      input$apply_param
      req(input$run_analysis >= 1)
      isolate({
        plot_genes <- unique(mutation_heat_df() %>% dplyr::select(mutation_id) %>% unlist())
        if (length(plot_genes) > input$max_genes) {
          plot_genes <- plot_genes[1:input$max_genes]
        }

        title <- if_else(is.null(input$custom_data) == F, "custom dataset", input$disease)

        fill_color <- c("grey", "#39568CFF")
        names(fill_color) <- c("0", "1")

        alpha <- c(0, 1)
        names(alpha) <- c("0", "1")

        mutation_heat_df() %>%
          filter(mutation_id %in% plot_genes) %>%
          ggplot(aes(x = patient_id, y = fct_rev(mutation_id), fill = mutation, alpha = as.factor(mutation))) +
          geom_tile() +
          custom_theme +
          scale_fill_manual(values = fill_color) +
          scale_alpha_manual(values = alpha) +
          labs(
            title = paste0("Mutations distribution across ", title, " patients"),
            x = "Patients", y = ""
          ) +
          theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = input$x_size),
            legend.position = "none"
          )
      }) # Isolate
    },
    height = graph_height_2
  ) # Renderplot

  ### === Length analysis



  output$graph_length <- renderPlot({
    input$apply_param
    input$run_length_analysis
    req(input$run_length_analysis >= 1)
    isolate({
      n_pts <- mutations_raw_data() %>%
        distinct(patient_id) %>%
        unlist() %>%
        length()

      graph_length_df() %>%
        ggplot(aes(x = cum_length / 1000, y = percent_tot, color = group)) +
        # geom_point() +
        geom_line(size = 1.5) +
        scale_color_viridis_d() +
        custom_theme +
        labs(
          title = paste0(
            "Panel length comparison in ", input$disease,
            " (n= ", n_pts, " pts)"
          ),
          x = "Total panel length (kb)", color = "Approach"
        ) +
        ylim(0, 1) +
        theme(legend.position = "right")
    }) # Isolate
  }) # Renderplot

  ### === Download

  # Download panel
  output$download_panel <- downloadHandler(
    filename = function() {
      paste("Informative_merged.tsv")
    },
    content = function(file) {
      write.table(informative_merged_df(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )

  output$download_sug_genes <- downloadHandler(
    filename = function() {
      paste("Suggested_genes.tsv")
    },
    content = function(file) {
      write.table(suggested_genes(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )


  output$download_panel_indiv <- downloadHandler(
    filename = function() {
      paste("Informative_individual.tsv")
    },
    content = function(file) {
      write.table(informative_df(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )

  # Download table mutation stats
  output$download_mut_stat <- downloadHandler(
    filename = function() {
      paste("Mutation_stats.tsv")
    },
    content = function(file) {
      write.table(table_mut_stat_df(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )

  # Download table mutation stats
  output$download_table_des <- downloadHandler(
    filename = function() {
      paste("Table_des.tsv")
    },
    content = function(file) {
      write.table(table_des_df(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )

  # Download table mutation freq
  output$download_mutfreq <- downloadHandler(
    filename = function() {
      paste("Mutation_types_and_frequencies.tsv")
    },
    content = function(file) {
      write.table(mutation_des_merged_df(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )

  # Download table length
  output$download_length <- downloadHandler(
    filename = function() {
      paste("Length_table.tsv")
    },
    content = function(file) {
      write.table(length_df(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )
}
