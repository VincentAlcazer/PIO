
#' PIO main algorithm: panel finder
#' Select an optimal set of gene to reach maximum informativity
#'
#' @param group_by TO which feature mutations should be grouped by: either gene_id or exon_id
#' @param info_metric Used metric to evaluate informativity: either UP (unique pts) or UPKB (unique pts per kb)
#' @param custom_panel Restrict features selection to a custom set of genes or exons. Can be NULL
#' @param min_pts Minimal number of patients per mutation (on the overall cohort)
#' @param min_mut Minimal number of mutations per patient
#' @export

panel_finder <-
  function(data_list,
           group_by = "gene_id",
           info_metric = "UPKB",
           custom_panel = NULL,
           min_pts = 2,
           min_mut = 2) {
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
        dplyr::select(mutation_id = group_by, length = length,  everything()) %>%
        as.data.table()

      tot_pts = length(unique(data_mutation$patient_id))

      if (is.null(custom_panel)) {
        custom_panel = unique(data_mutation$mutation_id)
      }

      data_mutation <-
        data_mutation %>% filter(mutation_id %in% custom_panel)

      ordered_mutations <- data_mutation[, .(count=.N,
                                             UP = length(unique(patient_id)),
                                             length=mean(length)),
                                         by=mutation_id] %>%
        filter(UP >= min_pts) %>%
        mutate(UPKB = round(UP * 1000 / length, 4)) %>%
        dplyr::select(metric = all_of(info_metric), UP = UP, UPKB = UPKB, everything()) %>%
        arrange(desc(metric))

      data_mutation <- data_mutation %>%
        filter(mutation_id %in% ordered_mutations$mutation_id)

      results_list <- list()
      genes_tested <- ordered_mutations$mutation_id[1]
      pts_tested <- NULL
      n_added = 10

      ##### ===== Run algorithm

      while (n_added > 0) {
        message("testing ", last(genes_tested))


        ## Calculate statistics for the current gene & panel & update vector
        new_pts <- data_mutation %>%
          filter(mutation_id %in% last(genes_tested) &
                   !patient_id %in% pts_tested) %>%
          distinct(patient_id) %>% unlist() %>% as.character()

        UP_added <- length(new_pts)

        pts_tested <- append(pts_tested, new_pts)

        n_comut <- data_mutation %>%
          filter(mutation_id %in% genes_tested &
                   patient_id %in% pts_tested) %>%
          group_by(patient_id) %>%
          summarise(n_mut = length(unique(mutation_id)))

        if(exists("comut_eval")){
          UP_comut <- as.numeric(comut_eval$UP[comut_eval$mutation_id == last(genes_tested)])
        } else {
          UP_comut = 0
        }
        if(isEmpty(UP_comut)){
          UP_comut = 0
        }

        ## Save results
        results_list[[last(genes_tested)]] <- data.frame(mutation_id = last(genes_tested),
                                                         step_UP = UP_added,
                                                         step_UP_comut = UP_comut,
                                                         n_comut_1 = sum(n_comut$n_mut >= 1),
                                                         n_comut_2 = sum(n_comut$n_mut >= 2),
                                                         n_comut_3 = sum(n_comut$n_mut >= 3),
                                                         n_comut_4 = sum(n_comut$n_mut >= 4),
                                                         n_comut_5 = sum(n_comut$n_mut >= 5))

        n_added =  results_list[[last(genes_tested)]]$step_UP + results_list[[last(genes_tested)]]$step_UP_comut

        ## Find the next gene to add:

        ### === new patients evaluation

        new_pt_eval <- data_mutation[!patient_id %in% pts_tested, .(count=.N,
                                                                    UP = length(unique(patient_id)),
                                                                    length=mean(length)),
                                     by=mutation_id] %>%
          mutate(UPKB = round(UP * 1000 / length, 4)) %>%
          dplyr::select(metric = all_of(info_metric), UP=UP, UPKB=UPKB,  everything()) %>%
          arrange(desc(metric))

        ### === comutations evaluation

        comut_pt_eval = data_mutation[mutation_id %in% genes_tested &patient_id %in% pts_tested,
                                      .(n_mut = length(unique(mutation_id))),
                                      by=patient_id] %>%
          filter(n_mut < min_mut)

        comut_eval <- data_mutation[patient_id %in% comut_pt_eval$patient_id &
                                      !mutation_id %in% genes_tested,
                                    .(UP = length(unique(patient_id)),
                                      length = mean(length)),
                                    by=mutation_id] %>%
          mutate(UPKB = round(UP * 1000 / length, 4)) %>%
          dplyr::select(metric = all_of(info_metric),UP = UP, UPKB = UPKB, everything()) %>%
          arrange(desc(metric))

        if(nrow(new_pt_eval) == 0 & nrow(comut_eval) == 0 ){
          break
        }
        if(nrow(new_pt_eval) == 0){
          new_pt_eval = data.frame(metric = 0,
                                   mutation_id = "NULL")
        }
        if(nrow(comut_eval) == 0){
          comut_eval = data.frame(metric = 0,
                                   mutation_id = "NULL")
        }


        ## Select next gene
        genes_tested <- append(genes_tested,
                               if_else(new_pt_eval$metric[1] >= comut_eval$metric[1],
                                       new_pt_eval$mutation_id[1],
                                       comut_eval$mutation_id[1]))


        # if (length(pts_tested) == tot_pts) {
        #   break
        # }

      }


      ##### ===== Bind results
      results_df_list[[datasets]] <- bind_rows(results_list)  %>%
        filter(is.na(mutation_id) == F) %>%
        mutate(
          percent_comut_1 = round(n_comut_1 / tot_pts, 3),
          percent_comut_2 = round(n_comut_2 / tot_pts, 3),
          percent_comut_3 = round(n_comut_3 / tot_pts, 3),
          percent_comut_4 = round(n_comut_4 / tot_pts, 3),
          percent_comut_5 = round(n_comut_5 / tot_pts, 3),
          cohort = datasets
        )

      ##### ===== Add genes or exon informations
      if (group_by == "gene_id") {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(distinct(
            dplyr::select(ordered_mutations, mutation_id, length, UP, UPKB),
            mutation_id,
            .keep_all = T
          ),
          by = c("mutation_id")) %>%
          mutate(cum_length = cumsum(length))

      } else {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(dplyr::select(ordered_mutations, mutation_id, length, UP, UPKB),
                    by = c("mutation_id")) %>%
          mutate(cum_length = cumsum(length))

      }


    }

    results_df <- bind_rows(results_df_list) %>%
      mutate(step_UPKB = round((step_UP / length),4),
             step_comut_UPKB = round((step_UP_comut / length),4)) %>%
      dplyr::select(mutation_id, UP, length, UPKB, step_UP, step_UP_comut, step_UPKB, step_comut_UPKB,  everything())

    return(results_df)


  }



#' Simple panel tester
#' Assess the informativity of a set of genes
#' @param group_by TO which feature mutations should be grouped by: either gene_id or exon_id
#' @param info_metric Used metric to evaluate informativity: either UP or UPKB
#' @param custom_panel Gene or exon list to test (as a character verctor)
#' @param keep_order Should the order of the provided panel be kept? (logical)
#' @param min_pts Minimal number of patients per mutation (on the overall cohort)
#' @export


panel_tester <-
  function(data_list,
           group_by = "gene_id",
           info_metric = "UP",
           custom_panel = NULL,
           keep_order = F,
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
            UP = length(unique(patient_id)),
            length = mean(length)
          ) %>%
          filter(UP >= min_pts) %>%
          mutate(UPKB = round(UP * 1000 / length, 4)) %>%
          dplyr::select(metric = all_of(info_metric),UP = UP, UPKB = UPKB, everything()) %>%
          arrange(desc(metric))

      }


      pts_tested <- NULL
      genes_tested <- NULL
      results_list <- list()

      for (mutations in ordered_mutations$mutation_id) {

        genes_tested <- append(genes_tested, mutations)

         new_pts <- data_mutation %>%
          filter(mutation_id == mutations &
                   !patient_id %in% pts_tested) %>%
          distinct(patient_id) %>% unlist() %>% as.character()

         UP_added <- length(new_pts)

        pts_tested <- append(pts_tested, new_pts)

        n_comut <- data_mutation %>%
          filter(mutation_id %in% genes_tested &
                   patient_id %in% pts_tested) %>%
          group_by(patient_id) %>%
          summarise(n_mut = length(unique(mutation_id)))



        results_list[[mutations]] <- data.frame(mutation_id = mutations,
                                                step_UP = UP_added,
                                                n_comut_1 = sum(n_comut$n_mut >= 1),
                                                n_comut_2 = sum(n_comut$n_mut >= 2),
                                                n_comut_3 = sum(n_comut$n_mut >= 3),
                                                n_comut_4 = sum(n_comut$n_mut >= 4),
                                                n_comut_5 = sum(n_comut$n_mut >= 5))



        if (length(pts_tested) == tot_pts) {
          break
        }

      }

      ##### ===== Bind results
      results_df_list[[datasets]] <- bind_rows(results_list)  %>%
        filter(is.na(mutation_id) == F) %>%
        mutate(
          percent_comut_1 = round(n_comut_1 / tot_pts, 3),
          percent_comut_2 = round(n_comut_2 / tot_pts, 3),
          percent_comut_3 = round(n_comut_3 / tot_pts, 3),
          percent_comut_4 = round(n_comut_4 / tot_pts, 3),
          percent_comut_5 = round(n_comut_5 / tot_pts, 3),
          cohort = datasets
        )

      ##### ===== Add genes or exon informations
      if (group_by == "gene_id") {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(distinct(
            dplyr::select(ordered_mutations, mutation_id, length, UP, UPKB),
            mutation_id,
            .keep_all = T
          ),
          by = c("mutation_id")) %>%
          mutate(cum_length = cumsum(length))

      } else {
        results_df_list[[datasets]] <- results_df_list[[datasets]] %>%
          left_join(dplyr::select(ordered_mutations, mutation_id, length, UP, UPKB),
                    by = c("mutation_id")) %>%
          mutate(cum_length = cumsum(length))

      }


    }


    results_df <- bind_rows(results_df_list) %>%
      mutate(step_UPKB = round((step_UP / length),4),
             step_UP_comut = NA,
             step_comut_UPKB = NA) %>%
      dplyr::select(mutation_id, UP, length, UPKB, step_UP, step_UP_comut, step_UPKB, step_comut_UPKB,  everything())


    return(results_df)


  }
