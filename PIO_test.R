

## parameters
input <- data.frame(
  disease = "BLCA",
  mutation_group = "gene_id",
  metric_type = "UPKB",
  min_pts = 2,
  min_mut=2,
  max_freq = 1,
  max_genes=200,
  x_size=6,
  analysis_type = "optimal"
)


### Custom panel
panel_df <- data.table::fread("Custom_panel_example.tsv") %>% as.data.frame()
colnames(panel_df)[1] <- "mutation_id"

custom_panel <- unlist(panel_df %>% dplyr::select(mutation_id)) %>% as.character()

files <- list.files(paste0("raw_mutation_data/",input$disease,"/"))

file_list <- list()

for (f in files) {
  name <- gsub(".rds", "", f)

  file_list[[name]] <- read_rds(paste0("raw_mutation_data/",input$disease,"/", f)) %>%
    mutate(cohort = name)
}

mutations_raw_data <- bind_rows(file_list)

merged_dataset <- list(merged = mutations_raw_data)
split_dataset <- split(mutations_raw_data, mutations_raw_data[,"cohort"])

informative_merged_df <- panel_finder(split_dataset,
             group_by = input$mutation_group,
             info_metric = input$metric_type, min_pts = input$min_pts, min_mut = input$min_mut,
             custom_panel = NULL
)

custom_panel <- read_tsv("Custom_dataset_example_optimized.tsv") %>% unlist() %>% as.character()

## Add the number of independent datasets
independent_sets <- mutations_raw_data %>%
  dplyr::select(mutation_id = input$mutation_group, everything()) %>%
  group_by(mutation_id) %>%
  summarise(n_cohort = length(unique(cohort)))

informative_merged_df <- informative_merged_df %>%
  left_join(independent_sets, by = "mutation_id") %>%
  droplevels()









##### Main plot

plot_genes <- unique(informative_merged_df %>% dplyr::select(mutation_id) %>% unlist())

if (length(plot_genes) > input$max_genes) {
  plot_genes <- plot_genes[1:input$max_genes]
}

data_plot <- informative_merged_df %>%
  filter(
    percent_comut_1 <= input$max_freq,
    mutation_id %in% plot_genes
  ) %>%
  droplevels() %>%
  pivot_longer(starts_with("percent_comut"), names_to = "p_comut_cat", values_to ="p_comut_val")

data_plot$mutation_id <- factor(data_plot$mutation_id, levels = plot_genes)

data_plot$p_comut_cat <- factor(data_plot$p_comut_cat,
                                labels = c(">=1",">=2",">=3",">=4",">=5"))

n_pts <- mutations_raw_data %>%
  distinct(patient_id) %>%
  unlist() %>%
  length()

data_anno <- data_plot %>%
  group_by(p_comut_cat) %>%
  summarise(y = round(max(p_comut_val),2),
            n_pts = round(max(p_comut_val) * n_pts),
            x = max(cum_length)/1000) %>%
  mutate(diff = (y - lead(y))/2) %>%
  mutate(y_pos = if_else(p_comut_cat == ">=5", y/2, y-diff))

ggplot(data = data_plot) +
  geom_line(aes(x = cum_length / 1000, y = p_comut_val,
                color = p_comut_cat, group = p_comut_cat),
            size = 1.2) +
  geom_text(data = data_anno,
                  aes(x = x, y = y_pos, label = paste0(y*100,"%"), color = p_comut_cat),
                  size = 6) +
  custom_theme +
  labs(
    title = paste0(
      "Most informative genes in ", if_else(is.null(input$custom_data), input$disease, "Custom dataset"),
      " (n= ", n_pts, " pts)"
    ),
    x = "Panel size (kb)", group = "test", y = "Total patients (%)",
    color = "Mutations", fill = "Mutations"
  ) +
  ylim(0, 1) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(
    legend.position = "right"
  )







#### Mutation type & freq

if (input$mutation_group == "exon_id") {
  feature_length <- "exon_length"
} else {
  feature_length <- "tot_exons_length"
}

if (is.null(input$panel_list) | input$analysis_type == "optimal") {
  custom_panel <- mutations_raw_data %>%
    dplyr::select(mutation_id = input$mutation_group) %>%
    unlist() %>%
    as.character()
} else {
  custom_panel <- unlist(panel_df %>% dplyr::select(mutation_id)) %>% as.character()
}

mutation_des_merged_df <- mutations_raw_data %>%
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



top_n_genes <- mutation_des_merged_df %>%
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


data_plot <- mutation_des_merged_df %>%
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
      " (n= ", unique(mutation_des_merged_df %>% dplyr::select(n_pts)), " pts)"
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
