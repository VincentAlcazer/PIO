
#' PIO main algorithm: panel finder
#' Select an optimal set of gene to reach maximum informativity
#'
#' @param group_by TO which feature mutations should be grouped by: either gene_id or exon_id
#' @param info_metric Used metric to evaluate informativity: either unique_pts or pts_per_kb
#' @param custom_panel Restrict features selection to a custom set of genes or exons. Can be NULL
#' @param min_pts Minimal number of patients per mutation (on the overall cohort)
#' @export

panel_finder <-
  function(data_list,
           group_by = "gene_id",
           info_metric = "unique_pts",
           custom_panel = NULL,
           min_pts = 1) {
    results_df_list <- list()

    if (group_by == "exon_id") {
      length = "exon_length"
    } else {
      length = "tot_exons_length"
    }


    for (datasets in names(data_list)) {
      message("Processing ", datasets)

      ##### ===== Dataset preparation
      data_mutation <- data_list[[datasets]] %>%
        dplyr::select(mutation_id = group_by, length = length,  everything())

      tot_pts = length(unique(data_mutation$patient_id))

      if (is.null(custom_panel)) {
        custom_panel = unique(data_mutation$mutation_id)
      }

      data_mutation <-
        data_mutation %>% filter(mutation_id %in% custom_panel)

      ordered_mutations <- data_mutation %>%
        group_by(mutation_id) %>%
        summarise(
          count = n(),
          unique_pts = length(unique(patient_id)),
          length = mean(length)
        ) %>%
        filter(unique_pts >= min_pts) %>%
        mutate(pts_per_kb = round(unique_pts * 1000 / length, 4)) %>%
        dplyr::select(metric = all_of(info_metric), everything()) %>%
        arrange(desc(metric))

      data_mutation <- data_mutation %>%
        filter(mutation_id %in% ordered_mutations$mutation_id)

      results_list <- list()
      genes_to_add <- ordered_mutations$mutation_id[1]
      pts_tested <- NULL

      n_added = 10

      ##### ===== Run algorithm

      while (n_added > 0) {
        message("testing ", last(genes_to_add))

        genes_tested <- genes_to_add

        new_pts <- data_mutation %>%
          filter(mutation_id %in% last(genes_to_add) &
                   !patient_id %in% pts_tested) %>%
          distinct(patient_id) %>% unlist() %>% as.character()

        n_added <- length(new_pts)

        results_list[[last(genes_tested)]] <- data.frame(n_pts = n_added,
                                                         mutation_id = last(genes_tested))

        pts_tested <- append(pts_tested, new_pts)

        next_gene <- data_mutation %>%
          filter(!patient_id %in% pts_tested) %>%
          group_by(mutation_id) %>%
          summarise(unique_pts = length(unique(patient_id)),
                    length = mean(length)) %>%
          mutate(pts_per_kb = round(unique_pts * 1000 / length, 4)) %>%
          dplyr::select(metric = all_of(info_metric), everything()) %>%
          arrange(desc(metric))

        genes_to_add <-
          append(genes_to_add, as.character(next_gene[1, 2]))

        if (length(pts_tested) == tot_pts) {
          break
        }

      }

      ##### ===== Bind results
      results_df_list[[datasets]] <- bind_rows(results_list)  %>%
        filter(is.na(mutation_id) == F) %>%
        mutate(
          cum_sum = cumsum(n_pts),
          percent_tot = round(cum_sum / tot_pts, 3),
          cohort = datasets
        )

      ##### ===== Add genes or exon informations
      if (group_by == "gene_id") {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(distinct(
            dplyr::select(ordered_mutations, mutation_id, length),
            mutation_id,
            .keep_all = T
          ),
          by = c("mutation_id")) %>%
          mutate(pts_per_kb = round(n_pts * 1000 / length, 4))

      } else {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(dplyr::select(ordered_mutations, mutation_id, length),
                    by = c("mutation_id")) %>%
          mutate(pts_per_kb = round(n_pts * 1000 / length, 4))

      }


    }

    results_df <- bind_rows(results_df_list)

    return(results_df)


  }



#' Simple panel tester
#' Assess the informativity of a set of genes
#' @param group_by TO which feature mutations should be grouped by: either gene_id or exon_id
#' @param info_metric Used metric to evaluate informativity: either unique_pts or pts_per_kb
#' @param custom_panel Gene or exon list to test (as a character verctor)
#' @param keep_order Should the order of the provided panel be kept? (logical)
#' @param min_pts Minimal number of patients per mutation (on the overall cohort)
#' @export


panel_tester <-
  function(data_list,
           group_by = "gene_id",
           info_metric = "unique_pts",
           custom_panel = NULL,
           keep_order = T,
           min_pts = 1) {
    results_df_list <- list()

    if (group_by == "exon_id") {
      length = "exon_length"
    } else {
      length = "tot_exons_length"
    }



    for (datasets in names(data_list)) {
      ##### ===== Dataset preparation
      data_mutation <- data_list[[datasets]] %>%
        dplyr::select(mutation_id = group_by, length = length, everything()) %>%
        filter(mutation_id %in% custom_panel)

      tot_pts = length(unique(data_list[[datasets]]$patient_id))

      if (keep_order == T) {
        ordered_mutations <- custom_panel

      } else {
        ordered_mutations <- data_mutation %>%
          group_by(mutation_id) %>%
          summarise(
            count = n(),
            unique_pts = length(unique(patient_id)),
            length = mean(length)
          ) %>%
          filter(unique_pts >= min_pts) %>%
          mutate(pts_per_kb = round(unique_pts * 1000 / length, 4)) %>%
          dplyr::select(metric = all_of(info_metric), everything()) %>%
          arrange(desc(metric)) %>%
          distinct(mutation_id) %>% unlist() %>% as.character()

      }


      pts_tested <- NULL

      results_list <- list()

      for (mutations in ordered_mutations) {
        new_pts <- data_mutation %>%
          filter(mutation_id == mutations &
                   !patient_id %in% pts_tested) %>%
          distinct(patient_id) %>% unlist() %>% as.character()

        n_added <- length(new_pts)

        results_list[[mutations]] <- data.frame(n_pts = n_added,
                                                mutation_id = mutations)


        pts_tested <- append(pts_tested, new_pts)

        if (length(pts_tested) == tot_pts) {
          break
        }

      }

      ##### ===== Bind results
      results_df_list[[datasets]] <- bind_rows(results_list)  %>%
        filter(is.na(mutation_id) == F) %>%
        mutate(
          cum_sum = cumsum(n_pts),
          percent_tot = round(cum_sum / tot_pts, 3),
          cohort = datasets
        )

      ##### ===== Add genes or exon informations
      if (group_by == "gene_id") {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(distinct(
            dplyr::select(data_mutation, mutation_id, length),
            mutation_id,
            .keep_all = T
          ),
          by = c("mutation_id")) %>%
          mutate(pts_per_kb = round(n_pts * 1000 / length, 4))

      } else {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(distinct(
            dplyr::select(data_mutation, mutation_id, length),
            mutation_id,
            .keep_all = T
          ),
          by = c("mutation_id")) %>%
          mutate(pts_per_kb = round(n_pts * 1000 / length, 4))

      }


    }


    results_df <- bind_rows(results_df_list)


    return(results_df)


  }
