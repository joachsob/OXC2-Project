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
      
      degs_df <- degs_df[, c("gene", "avg_log2FC", "logCPM", "P.Value", "p_val_adj")]
      
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
  list(full_name = "Arterial ECs", short_name = "Arterial ECs"),
  list(full_name = "Atrial ECs", short_name = "Atrial ECs"),
  list(full_name = "Ventricular Cardiomyocytes", short_name = "Ventricular Cardiomyocytes"),
  list(full_name = "Macrophages", short_name = "Macrophages"),
  list(full_name = "Cardiac Fibroblasts", short_name = "Cardiac Fibroblasts"),
  list(full_name = "Schwann Cells", short_name = "Schwann Cells"),
  list(full_name = "Pericytes", short_name = "Pericytes"),
  list(full_name = "T-Lymphocytes", short_name = "T-Lymphocytes"),
  list(full_name = "ECM-High Cardiac Fibroblasts", short_name = "ECM-High Cardiac Fibroblasts"),
  list(full_name = "Myofibroblasts", short_name = "Myofibroblasts"),
  list(full_name = "mt-High Ventricular Cardiomyocytes", short_name = "mt-High Ventricular Cardiomyocytes"),
  list(full_name = "Cardiac Neurons", short_name = "Cardiac Neurons"),
  list(full_name = "Capillary ECs", short_name = "Capillary ECs"),
  list(full_name = "Lymphatic ECs", short_name = "Lymphatic ECs"),
  list(full_name = "Vascular Smooth Muscle Cells", short_name = "VSMC")
)

unique(seuratObj$celltype)

# Define condition pairs for differential expression analysis
condition_pairs <- list(
  c("Patient", "Control")
)

# Define short names for conditions
condition_shortnames <- list(
  "Control" = "CTRL",
  "Patient" = "PTNT"
)


number_degs <- count_DEGs_wilcox(seuratObj, 
                                 celltype_list = celltype_list, 
                                 condition_pairs = condition_pairs, 
                                 condition_shortnames = condition_shortnames, 
                                 group_by1 = 'celltype', group_by2 = 'condition',
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