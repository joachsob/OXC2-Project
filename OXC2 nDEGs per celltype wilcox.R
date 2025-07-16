library(Seurat)
library(openxlsx)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(Matrix)
library(limma)
library(msigdbr)
library(dplyr)
library(BiocParallel)
register(SerialParam())
library(readxl)
library(writexl)
library(presto)
library(tidyr)
library(ggplot2)
library(knitr) # for formatted talbes (function kable)

seuratObj <- readRDS('oxcwtc_seurat.RDS')
# newAnnot <- read.xlsx("Annotation.xlsx") # Fixed the celltype names containing a slash
# seuratObj$celltype <- newAnnot$celltype[match(seuratObj$seurat_clusters, newAnnot$cluster)] # match the new celltype names to the seuratObj
# saveRDS(seuratObj, "oxcwtc_seurat.RDS") # save the updated seuratObj in the RDS-file

seuratObj <- PrepSCTFindMarkers(seuratObj)

# Creating an arbitrary markers list to check column names for adjusting the DEG function
testMarkers <- FindMarkers(seuratObj, group.by = 'celltype', ident.1 = 'vRG')
colnames(testMarkers)

#### Actually correct code ####

count_DEGs_wilcox <- function(seuratObj, celltype_list, condition_pairs, condition_shortnames, 
                              group_by1, group_by2, assay = "RNA") {
  
  gene_counts <- list()  # Store DEG counts per cell type
  
  for (celltype in celltype_list) {
    full_name <- celltype$full_name
    short_name <- celltype$short_name
    
    celltype_subset <- subset(seuratObj, cells = which(seuratObj[[group_by1]] == full_name))
    
    for (condition_pair in condition_pairs) {
      condition1 <- condition_pair[[1]]
      condition2 <- condition_pair[[2]]
      
      seuratObj_sub <- celltype_subset
      degs_df <- FindMarkers(
        celltype_subset, group.by = group_by2, assay = assay, verbose = TRUE,
        ident.1 = condition1, ident.2 = condition2, min.pct = 0.2
      )
      
      degs_df <- degs_df[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      
      # Filter DEGs by log2 fold change > 1
      degs_filtered <- subset(degs_df, abs(avg_log2FC) > 1 & p_val_adj < 0.05)
      
      # Store the count of filtered DEGs
      gene_counts[[short_name]] <- nrow(degs_filtered)
    }
  }
  
  # Convert the list to a data frame for plotting
  gene_counts_df <- data.frame(CellType = names(gene_counts), DEG_Count = unlist(gene_counts))
  
  # Create the bar chart
  number_degs <- ggplot(gene_counts_df, aes(x = reorder(CellType, DEG_Count), y = DEG_Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Number of DEGs per Cell Type", x = "Cell Type", y = "DEG Count") +
    theme_minimal()
  
  return(number_degs)
}


# full celltype list for DEGs of control, mig, and ace-leu
celltype_list <- list(
  list(full_name = "vRG", short_name = "vRG"),
  list(full_name = "Granule Cells", short_name = "Granule Cells"),
  list(full_name = "NPCs", short_name = "NPCs"),
  list(full_name = "Proliferating Progenitors", short_name = "Proliferating Progenitors"),
  list(full_name = "SORCS1+ Immature Ex. Neurons", short_name = "SORCS1+ Immature Ex. Neurons"),
  list(full_name = "GAD1GAD2+ Granule Cells", short_name = "GAD1GAD2+ Granule Cells"),
  list(full_name = "oRG", short_name = "oRG"),
  list(full_name = "Neuroblasts", short_name = "Neuroblasts"),
  list(full_name = "ARPP21+ Immature Ex. Neurons", short_name = "ARPP21+ Immature Ex. Neurons"),
  list(full_name = "RELNGAD2+ Inh. Neurons", short_name = "RELNGAD2+ Inh. Neurons"),
  list(full_name = "Migratory Granule Cells", short_name = "Migratory Granule Cells"),
  list(full_name = "GPC5GAD2+ Inh. Neurons", short_name = "GPC5GAD2+ Inh. Neurons"),
  list(full_name = "Preplate Neurons", short_name = "Preplate Neurons"),
  list(full_name = "ATP1A2+ Fibroblast-Like", short_name = "ATP1A2+ Fibroblast-Like") # få celler for combined treatment, fjernes for DEG analyse
)

# unique(seuratObj$celltype)

# Define short names for conditions
condition_shortnames <- list(
  "OXC" = "OXC",
  "WTC" = "WTC"
)

#### DEGs ####

# Define condition pairs for differential expression analysis
condition_pairs <- list(
  c("4-OXC2-DMSO-control", "11-WTC-control")
)

# DEGs for control treatment
number_degs_control <- count_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'sample',
                                 assay = 'RNA')

number_degs_control +
  geom_bar(stat = "identity", fill = "#6BAED6") +  # Nature-style muted blue
  coord_flip() +
  labs(title = "Number of DEGs per Cell Type (Control)", x = "Cell Type", y = "DEG Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

# DEGs for Miglustat treatment
condition_pairs <- list(
  c("1-OXC2-Miglustat-100uM", "11-WTC-control")
)

number_degs_Mig <- count_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'sample',
                                 assay = 'RNA')

number_degs_Mig +
  geom_bar(stat = "identity", fill = "#6BAED6") +  # Nature-style muted blue
  coord_flip() +
  labs(title = "Number of DEGs per Cell Type (Miglustat)", x = "Cell Type", y = "DEG Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

# DEGs for ace-leu treatment
condition_pairs <- list(
  c("2-OXC2-Acetyl-leucine", "11-WTC-control")
)

number_degs_ace_leu <- count_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'sample',
                                 assay = 'RNA')

number_degs_ace_leu +
  geom_bar(stat = "identity", fill = "#6BAED6") +  # Nature-style muted blue
  coord_flip() +
  labs(title = "Number of DEGs per Cell Type (ace-leu)", x = "Cell Type", y = "DEG Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

# DEGs for combined treatment.
# This analysis excludes ATP1A2+ Fibroblast-Like from the cellist because these have too few cells
# The code is therefore modified to circumvent this:


# List excluding problematic celltype:
celltype_list <- list(
  list(full_name = "vRG", short_name = "vRG"),
  list(full_name = "Granule Cells", short_name = "Granule Cells"),
  list(full_name = "NPCs", short_name = "NPCs"),
  list(full_name = "Proliferating Progenitors", short_name = "Proliferating Progenitors"),
  list(full_name = "SORCS1+ Immature Ex. Neurons", short_name = "SORCS1+ Immature Ex. Neurons"),
  list(full_name = "GAD1GAD2+ Granule Cells", short_name = "GAD1GAD2+ Granule Cells"),
  list(full_name = "oRG", short_name = "oRG"),
  list(full_name = "Neuroblasts", short_name = "Neuroblasts"),
  list(full_name = "ARPP21+ Immature Ex. Neurons", short_name = "ARPP21+ Immature Ex. Neurons"),
  list(full_name = "RELNGAD2+ Inh. Neurons", short_name = "RELNGAD2+ Inh. Neurons"),
  list(full_name = "Migratory Granule Cells", short_name = "Migratory Granule Cells"),
  list(full_name = "GPC5GAD2+ Inh. Neurons", short_name = "GPC5GAD2+ Inh. Neurons"),
  list(full_name = "Preplate Neurons", short_name = "Preplate Neurons")
  # list(full_name = "ATP1A2+ Fibroblast-Like", short_name = "ATP1A2+ Fibroblast-Like") # få celler for combined treatment, fjernes for DEG analyse
)

condition_pairs <- list(
  c("3-OXC2-Mig-ace-leu", "11-WTC-control")
)


number_degs_combined <- count_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'sample',
                                 assay = 'RNA')

number_degs_combined +
  geom_bar(stat = "identity", fill = "#6BAED6") +  # Nature-style muted blue
  coord_flip() +
  labs(title = "Number of DEGs per Cell Type (Combined)", x = "Cell Type", y = "DEG Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

# Now include the missing celltype and result to the number_DEGs_combined$data dataframe to continue with analysis
missing_cell <- data.frame(CellType = "ATP1A2+ Fibroblast-Like", DEG_Count = NA)
number_degs_combined$data <- rbind(number_degs_combined$data, missing_cell)

#### Manage DEG-data ####
# the next few lines were a one-time formatting of the data frame and excel file
# wb <- loadWorkbook("./DEGs/completeDEGdata.xlsx")
# deg_df <- data.frame(
#   Celltype = c(
#     "vRG",
#     "Granule Cells",
#     "NPCs",
#     "Proliferating Progenitors",
#     "SORCS1+ Immature Ex. Neurons",
#     "GAD1GAD2+ Granule Cells",
#     "oRG",
#     "Neuroblasts",
#     "ARPP21+ Immature Ex. Neurons",
#     "RELNGAD2+ Inh. Neurons",
#     "Migratory Granule Cells",
#     "GPC5GAD2+ Inh. Neurons",
#     "Preplate Neurons",
#     "ATP1A2+ Fibroblast-Like"
#     ),
# 
#   ControlDEGs = number_degs_control$data$DEG_Count,
# 
#   MigDEGs = number_degs_Mig$data$DEG_Count,
# 
#   AceLeuDEGs = number_degs_ace_leu$data$DEG_Count,
# 
#   CombinedDEGs = number_degs_combined$data$DEG_Count
#   )
# 
# write_xlsx(deg_df, path = "./DEGs/completeDEGdata.xlsx")
# 
# deg_data <- read_excel("./DEGs/completeDEGdata.xlsx")
# deg_data$CombinedDEGs <- as.numeric(deg_data$CombinedDEGs) # after adjusting this dataframe the datatype is set to chr, revert to numeric
# 
# # update the excel file with corrected datatype for combined treatment values
# write_xlsx(deg_data, path = "./DEGs/completeDEGdata.xlsx")

deg_data <- read_excel("./DEGs/completeDEGdata.xlsx")

# create long-format to collect a frame with celltypes, treatments, and DEG count
deg_long <- pivot_longer(deg_data, 
                         cols = -Celltype, 
                         names_to = "Treatment", 
                         values_to = "DEG_Count")

#### Plotting the DEGs ####

# Line for ordering the treatments in metadata for the treatment to have correct order in the plot
deg_long$Treatment <- factor(deg_long$Treatment, levels = c("ControlDEGs", "MigDEGs", "AceLeuDEGs", "CombinedDEGs"))

# Setting the order in which the cells appear on the x-axis so similar celltypes are close to eachother
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

deg_long$Celltype <- factor(deg_long$Celltype, levels = celltype_order)

# plotting av DEGs
ggplot(deg_long, aes(x = Celltype, y = DEG_Count, fill = Treatment)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
        theme_minimal() +
        labs(title = "Number of significant DEGs per Celltype and Treatment",
             x = "Cell Type",
             y = "Number of Significant DEGs (p.adj < 0.05 & abs(avg.log2FC) > 1))") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(hjust = 1, angle = 60),
    axis.title.y = element_text(size = 10, margin = margin(l = 10)),
    axis.title.x = element_text(size = 10, margin = margin(b = 10))
    ) +
  
  scale_fill_manual(values = c("ControlDEGs" = "grey",
                               "MigDEGs" = "steelblue",
                               "AceLeuDEGs" = "darkorange",
                               "CombinedDEGs" = "forestgreen"),
                    labels = c(
                      "ControlDEGs" = "DMSO (Control)",
                      "MigDEGs" = "Miglustat (100uM)",
                      "AceLeuDEGs" = "Acetyl-Leucine (1mg/ml)",
                      "CombinedDEGs" = "Combined Miglustat (100uM)\n + AceLeu (1mg/ml)"
                    ))

#### Plotting for specific treatments ####
deg_long %>% 
  filter(Treatment %in% c("ControlDEGs", "CombinedDEGs")) %>% # this line filters which treatments are included in the plot
  ggplot(aes(x = Celltype, y = DEG_Count, fill = Treatment)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
    theme_minimal() +
    labs(title = "Number of significant DEGs per Celltype and Treatment",
        x = "Cell Type",
        y = "Number of Significant DEGs (p.adj < 0.05 & abs(avg.log2FC) > 1))") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(hjust = 1, angle = 60),
      axis.title.y = element_text(size = 10, margin = margin(l = 10)),
      axis.title.x = element_text(size = 10, margin = margin(b = 10))
    ) +
  
    scale_fill_manual(values = c("ControlDEGs" = "grey",
                               "MigDEGs" = "steelblue",
                               "AceLeuDEGs" = "darkorange",
                               "CombinedDEGs" = "forestgreen"),
                    labels = c(
                      "ControlDEGs" = "DMSO (Control)",
                      "MigDEGs" = "Miglustat (100uM)",
                      "AceLeuDEGs" = "Acetyl-Leucine (1mg/ml)",
                      "CombinedDEGs" = "Combined Miglustat (100uM)\n + AceLeu (1mg/ml)"
                    ))

#### Trying to sort for number of cell of each celltypes in each treatment for debug purposes ####

    # nCells_celltype_treatment <- table(seuratObj$sample, seuratObj$celltype)
    # nCells_celltype_treatment <- as.data.frame(nCells_celltype_treatment)
    # nCells_celltype_treatment <- nCells_celltype_treatment %>% rename(Cellcount = Freq)
    # write.xlsx(nCells_celltype_treatment, "./nCells per celltype per treatment.xlsx")

# get the dataframe containing information on cellcount for each celltype and treatment
nCells_celltype_treatment <- read_xlsx("./nCells per celltype per treatment.xlsx")

# set new name for the treatments in the cellcount file/dataframe ----
nCellLabel_map <- data.frame(
  Treatment = c(
    "4-OXC2-DMSO-control",
    "1-OXC2-Miglustat-100uM",
    "2-OXC2-Acetyl-leucine",
    "3-OXC2-Mig-ace-leu",

    "11-WTC-control",
    "6-WTC-Miglustat-100uM",
    "8-WTC-Acetyl-leucine",
    "9-WTC-Mig-ace-leu"
  ),
  Label = c(
    "nCell_Control_OXC",
    "nCell_Mig_OXC",
    "nCell_AceLeu_OXC",
    "nCell_Combined_OXC",

    "nCell_Control_WTC",
    "nCell_Mig_WTC",
    "nCell_AceLeu_WTC",
    "nCell_Combined_WTC"
  ),

  DEG_label = c(
    "ControlDEGs",
    "MigDEGs",
    "AceLeuDEGs",
    "CombinedDEGs"
  )
)
# 
# nCells_celltype_treatment <- nCells_celltype_treatment %>%
#   left_join(nCellLabel_map, by = "Treatment") %>%
#   relocate(Label, .before = Treatment)
# 
# nCells_celltype_treatment
# write_xlsx(nCells_celltype_treatment, "./nCells per celltype per treatment.xlsx")

#### ####

# sort the number of cells in only the OXC2 samples by treatment, disregarding celltype
sort_nCells_OXC2 <- nCells_celltype_treatment %>% 
  filter(grepl("OXC2", Treatment)) %>% #'grepl' for finding samples containing the 'OXC2' phrase
  group_by(Treatment) %>% # grouping the data by treatment
  summarise(Total_Cells = sum(Cellcount)) # %>% # finding the total number of cells in each treatment
  # arrange(Total_Cells) # sort them in ascending order
sort_nCells_OXC2$Treatment <- factor(sort_nCells_OXC2$Treatment, levels = c(
  "4-OXC2-DMSO-control", "1-OXC2-Miglustat-100uM", "2-OXC2-Acetyl-leucine", "3-OXC2-Mig-ace-leu"
))
# sort_nCells_OXC2 <- arrange(sort_nCells_OXC2, Treatment)
sort_nCells_OXC2

# same process as for sort_nCells_OXC2
sort_nCells_WTC <- nCells_celltype_treatment %>% 
  filter(grepl("WTC", Treatment)) %>%
  group_by(Treatment) %>%
  summarise(Total_Cells = sum(Cellcount)) #%>%
  # arrange(Total_Cells)
sort_nCells_WTC$Treatment <- factor(sort_nCells_WTC$Treatment, levels = c(
  "11-WTC-control", "6-WTC-Miglustat-100uM","8-WTC-Acetyl-leucine", "9-WTC-Mig-ace-leu"))
sort_nCells_WTC <- arrange(sort_nCells_WTC, Treatment)
sort_nCells_WTC

# extract the information on combined and control treatment for OXC2 specifically for the celltypes excluded from the GSEA
# to look at the cellcounts.
# filter cellcount for specific celltypes in specific treatments
filtered_nCells <- nCells_celltype_treatment %>% filter(
  # Treatment %in% c("3-OXC2-Mig-ace-leu", "4-OXC2-DMSO-control"),
  Celltype %in% c("Neuroblasts", "ATP1A2+ Fibroblast Like")
  )
filtered_nCells

# for getting a formatted table of the filtered cells (above)
kable(sort_nCells_OXC2, format = "markdown")

#### Cellcount per celltype per group (whole OXC2 and WTC) ####

# getting number of cells for a specific celltype in each treatment
# change the name of the celltypeOI to the celltype of interest
celltypeOI <- "ARPP21+ Immature Ex. Neurons"
nCelltype <- nCells_celltype_treatment %>% 
  filter(grepl("OXC2|WTC", Treatment), # extract line if 'OXC2' or 'WTC' is in the Treatment-name
         Celltype %in% celltypeOI) %>% # extract line if the Celltype is 'Neuroblast'
         group_by(Treatment) %>% # group the result by Treatment (the same treatments will be next to eachother)
         summarise(Total_Cells = sum(Cellcount))
nCelltype

# get number of celltype per group (total over all oxc2 treatments and wtc treatments)
nCelltype <- nCelltype %>%
  mutate(Source = ifelse(grepl("OXC2", Treatment), "OXC2", "WTC")) %>% #add column 'Source', add OXC2 if 'OXC2' is in the Treatment-name, else 'WTC'
  group_by(Source) %>%
  summarise(Total = sum(Total_Cells))
colnames(nCelltype)[2] <- paste0("Total_", celltypeOI) # change column name of cellcount dynamically
nCelltype

#### ####
# Formatted table 
kable(nCelltype, format = "markdown")


# 