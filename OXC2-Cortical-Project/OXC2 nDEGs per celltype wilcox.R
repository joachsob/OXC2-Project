# seuratObj <- PrepSCTFindMarkers(seuratObj)


fcount_DEGs_wilcox <- function(seuratObj, celltype_list, condition_pairs, condition_shortnames, 
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
  list(full_name = "ATP1A2+ Fibroblast-Like", short_name = "ATP1A2+ Fibroblast-Like")
)

unique(seuratObj$celltype)

# Define condition pairs for differential expression analysis
condition_pairs <- list(
  c("4-OXC2-DMSO-control", "11-WTC-control")
)

# Define short names for conditions
condition_shortnames <- list(
  "OXC" = "OXC",
  "WTC" = "WTC"
)

number_degs <- fcount_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'sample',
                                 assay = 'RNA')

number_degs +
  geom_bar(stat = "identity", fill = "#6BAED6") +  # Nature-style muted blue
  coord_flip() +
  labs(title = "Number of DEGs per Cell Type", x = "Cell Type", y = "DEG Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

# Creating an arbitrary markers list to check column names for adjusting the DEG function
# testMarkers <- FindMarkers(seuratObj, group.by = 'celltype', ident.1 = 'vRG')

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
  list(full_name = "ATP1A2+ Fibroblast-Like", short_name = "ATP1A2+ Fibroblast-Like")
)

unique(seuratObj$celltype)

# Define condition pairs for differential expression analysis
condition_pairs <- list(
  c("4-OXC2-DMSO-control", "11-WTC-control")
)

# Define short names for conditions
condition_shortnames <- list(
  "OXC" = "OXC",
  "WTC" = "WTC"
)

number_degs <- count_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'sample',
                                 assay = 'RNA')

number_degs +
  geom_bar(stat = "identity", fill = "#6BAED6") +  # Nature-style muted blue
  coord_flip() +
  labs(title = "Number of DEGs per Cell Type", x = "Cell Type", y = "DEG Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )
