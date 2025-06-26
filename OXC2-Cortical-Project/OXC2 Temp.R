nCells_celltype_treatment$Treatment
deg_data

total_oRG <- nCells_celltype_treatment %>%
  filter(Celltype %in% "oRG") %>%
  group_by(Treatment) %>%
  summarise(Total_oRG = sum(Cellcount)) %>%
  arrange(Treatment)
total_oRG

# this code is from ChatGPT, but the column names do not correspond to the dataframes
plot_df <- deg_long %>%
  filter(Celltype == "oRG") %>%
  left_join(total_oRG, by = "Treatment") %>%
  mutate(DEG_per_cell = DEG_Count / total_oRG)



nCellLabel_map$Condition <- gsub("nCell_(.*?)_.*", "\\1", nCellLabel_map$Label)
nCellLabel_map

combined_oRG <- nCells_celltype_treatment %>%
  filter(Celltype %in% "oRG") %>%
  left_join(nCellLabel_map, by ="Treatment") %>%
  group_by(Condition) %>%
  summarise(Total_oRG = sum(Cellcount)) %>%
  arrange(Condition)
combined_oRG

deg_data

deg_data %>%
  filter(Celltype %in% "oRG")
