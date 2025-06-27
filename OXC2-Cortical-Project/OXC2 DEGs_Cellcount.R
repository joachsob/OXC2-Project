nCells_celltype_treatment$Treatment
deg_data

nCellLabel_map$Condition <- gsub("nCell_(.*?)_.*", "\\1", nCellLabel_map$Label)
nCellLabel_map

oRG_condition <- nCells_celltype_treatment %>%
  filter(Celltype %in% "oRG") %>%
  group_by(Treatment) %>%
  summarise(Total_oRG = sum(Cellcount)) %>%
  arrange(Treatment)
oRG_condition

oRG_condition <- nCells_celltype_treatment %>%
  filter(Celltype == "oRG") %>%
  left_join(nCellLabel_map, by = "Treatment") %>%
  group_by(Condition) %>%
  summarise(Total_oRG = sum(Cellcount)) %>%
  arrange(Condition)
oRG_condition


deg_long <- deg_long %>%
  mutate(Condition = gsub("DEGs", "", Treatment))


# this code is from ChatGPT, but the column names do not correspond to the dataframes
plot_df_oRG <- deg_long %>%
  filter(Celltype == "oRG") %>%
  left_join(total_oRG, by = "Condition") %>%
  mutate(DEG_per_cell = DEG_Count / Total_oRG)
plot_df_oRG


cell_counts_condition <- nCells_celltype_treatment %>%
  left_join(nCellLabel_map, by = "Treatment") %>%
  select(Condition, Celltype, Cellcount)

total_cells <- cell_counts_condition %>%
  group_by(Celltype, Condition) %>%
  summarise(Total_Cellcount = sum(Cellcount), .groups = "drop")
total_cells

deg_normalized <- deg_long %>%
  left_join(total_cells, by = c("Celltype", "Condition")) %>%
  mutate(DEG_per_cell = DEG_Count / Total_Cellcount)
deg_normalized

deg_normalized <- deg_normalized %>%
  select(-Treatment)


library(ggplot2)

celltype_order <- c(
  "NPCs",
  "Proliferating Progenitors",
  "Neuroblasts",
  "vRG",
  "oRG",
  "Granule Cells",
  "GAD1GAD2+ Granule Cells",
  "Migratory Granule Cells",
  "SORCS1+ Immature Ex. Neurons",
  "ARPP21+ Immature Ex. Neurons",
  "Preplate Neurons",
  "RELNGAD2+ Inh. Neurons",
  "GPC5GAD2+ Inh. Neurons",
  "ATP1A2+ Fibroblast-Like"
)

condition_order <- c(
  "Control",
  "Mig",
  "AceLeu",
  "Combined"
)

deg_normalized$Celltype <- factor(deg_normalized$Celltype, levels = celltype_order)
deg_normalized$Condition <- factor(deg_normalized$Condition, levels = condition_order)

deg_long$Celltype <- factor(deg_long$Celltype, levels = celltype_order)
deg_long$Condition <- factor(deg_long$Condition, levels = condition_order)

total_cells$Celltype <- factor(total_cells$Celltype, levels = celltype_order)
total_cells$Condition <- factor(total_cells$Condition, levels = condition_order)

# plotting ----

  ## plotting av cellcount ----
ggplot(total_cells, aes(x = Celltype, y = Total_Cellcount, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  theme_minimal() +
  labs(title = "Number of Cells per Cell Type and Treatment",
       x = "Cell Type", y = "Number of Cells") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    values = c(
      'Control' = "grey",
      'Mig' = "steelblue",
      "AceLeu" = "darkorange",
      "Combined" = "forestgreen"
    ),
    labels = c(
      "Control" = "DMSO (Control)",
      "Mig" = "Miglustat (100uM)",
      "AceLeu" = "Acetyl-Leucine (1mg/ml)",
      "Combined" = "Combined Miglustat (100uM)\n + AceLeu (1mg/ml)"
    )
  )

  ## plotting av DEGcount ----

ggplot(deg_long, aes(x = Celltype, y = DEG_Count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
  theme_minimal() +
  labs(title = "Number of Significant DEGs per Celltype and Treatment",
       x = "Cell Type",
       y = "Number of Significant DEGs (p.adj < 0.05 & abs(avg.log2FC) > 1))") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(hjust = 1, angle = 60),
    axis.title.y = element_text(size = 10, margin = margin(l = 10)),
    axis.title.x = element_text(size = 10, margin = margin(b = 10))
  ) +
  
  scale_fill_manual(values = c(
    'Control' = "grey",
    'Mig' = "steelblue",
    "AceLeu" = "darkorange",
    "Combined" = "forestgreen"
  ),
  labels = c(
    "Control" = "DMSO (Control)",
    "Mig" = "Miglustat (100uM)",
    "AceLeu" = "Acetyl-Leucine (1mg/ml)",
    "Combined" = "Combined Miglustat (100uM)\n + AceLeu (1mg/ml)"
  )
                  )

  ## plotting av DEGs/Cellcount ----
ggplot(deg_normalized, aes(x = Celltype, y = DEG_per_cell, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  theme_minimal() +
  labs(title = "DEGs per Cell Count by Cell Type and Treatment",
       x = "Cell Type", y = "DEGs per Cellcount") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    values = c(
    'Control' = "grey",
    'Mig' = "steelblue",
    "AceLeu" = "darkorange",
    "Combined" = "forestgreen"
    ),
    labels = c(
      "Control" = "DMSO (Control)",
      "Mig" = "Miglustat (100uM)",
      "AceLeu" = "Acetyl-Leucine (1mg/ml)",
      "Combined" = "Combined Miglustat (100uM)\n + AceLeu (1mg/ml)"
    )
  )

  



