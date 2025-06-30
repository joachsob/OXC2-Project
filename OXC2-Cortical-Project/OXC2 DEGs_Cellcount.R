library(Seurat)
library(openxlsx)
library(dplyr)
library(readxl)
library(writexl)
library(presto)
library(tidyr)
library(ggplot2)
library(knitr) # for formatted tables (function kable)

seuratObj <- readRDS('oxcwtc_seurat.RDS')

deg_data <- read_excel("./DEGs/completeDEGdata.xlsx")
# create long-format to collect a frame with celltypes, treatments, and DEG count
deg_long <- pivot_longer(deg_data, 
                         cols = -Celltype, 
                         names_to = "Treatment", 
                         values_to = "DEG_Count")
nCells_celltype_treatment <- read_xlsx("./nCells per celltype per treatment.xlsx")

head(nCells_celltype_treatment)
head(deg_data)
head(deg_long)

nCellLabel_map$Condition <- gsub("nCell_(.*?)_.*", "\\1", nCellLabel_map$Label)
nCellLabel_map

# number of a celltype per condition
celltype_condition <- nCells_celltype_treatment %>%
  filter(Celltype == "oRG") %>%
  left_join(nCellLabel_map, by = "Treatment") %>%
  group_by(Treatment, Condition) %>%
  summarise(Total_nCell = sum(Cellcount)) %>%
  arrange(Condition)
celltype_condition


deg_long <- deg_long %>%
  mutate(Condition = gsub("DEGs", "", Treatment))
deg_long

#### Redone to match OXC2 treatment to WTC DMSO ####
nCellLabel_map$Group <- ifelse(grepl("OXC", nCellLabel_map$Label), "OXC2", "WTC")

cell_counts_condition <- nCells_celltype_treatment %>%
  left_join(nCellLabel_map, by = "Treatment") %>%
  select(Treatment, Group, Condition, Celltype, Cellcount)
cell_counts_condition

# Step 1: Get the D group cell counts (they have the true treatment label)
cell_counts_D <- cell_counts_condition %>%
  filter(Group == "OXC2") %>%
  select(Condition, Celltype, Cellcount)

# Step 2: Get the H group control cell counts (shared across treatments)
cell_counts_H_control <- cell_counts_condition %>%
  filter(Group == "WTC", Condition == "Control") %>%
  select(Celltype, Cellcount)

# Step 3: Expand H-Control cell counts to each treatment condition
conditions <- unique(cell_counts_D$Condition)
expanded_H <- expand.grid(
  Condition = conditions,
  Celltype = unique(cell_counts_H_control$Celltype)
) %>%
  left_join(cell_counts_H_control, by = "Celltype")

# Step 4: Combine D and replicated H rows
deg_denominator <- bind_rows(cell_counts_D, expanded_H) %>%
  group_by(Condition, Celltype) %>%
  summarise(Total_Cellcount = sum(Cellcount), .groups = "drop")


deg_normalized_corrected <- deg_long %>%
  left_join(deg_denominator, by = c("Condition", "Celltype")) %>%
  mutate(DEG_per_cell = DEG_Count / Total_Cellcount)


#### ####

total_cells <- cell_counts_condition %>%
  group_by(Celltype, Condition) %>%
  summarise(Total_Cellcount = sum(Cellcount), .groups = "drop")
total_cells

# should I count the number of cells in OXC2 and WTC separately?
total_cells_group <- cell_counts_condition %>%
  group_by(Group, Celltype, Condition) %>%
  summarise(Total_Cellcount = sum(Cellcount), .groups = "drop")
total_cells_group

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

deg_normalized_corrected$Celltype <- factor(deg_normalized$Celltype, levels = celltype_order)
deg_normalized_corrected$Condition <- factor(deg_normalized$Condition, levels = condition_order)

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
ggplot(deg_normalized_corrected, aes(x = Celltype, y = DEG_per_cell, fill = Condition)) +
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

## Plotting of DEG/cell for specific treatments ----
deg_normalized_corrected %>%
  filter(Condition %in% c("Control", "Combined")) %>%
  ggplot(aes(x = Celltype, y = DEG_per_cell, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  theme_minimal() +
  labs(title = "DEGs per Cell by Cell Type and Treatment",
       x = "Cell Type", y = "DEGs per Cell") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    values = c(
      'Control' = "grey",
      "Combined" = "forestgreen"
    ),
    labels = c(
      "Control" = "DMSO (Control)",
      "Combined" = "Combined Miglustat (100uM)\n + AceLeu (1mg/ml)"
    )
  )

  



