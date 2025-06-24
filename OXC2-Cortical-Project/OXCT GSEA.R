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

seuratObj <- readRDS('oxcwtc_seurat.RDS')
newAnnot <- read.xlsx("Annotation.xlsx") # Fixed the celltype names containing a slash
seuratObj$celltype <- newAnnot$celltype[match(seuratObj$seurat_clusters, newAnnot$cluster)]

FeaturePlot(seuratObj, "ELOVL7", order = T)

perform_DEG_and_GSEA_analysis <- function(seuratObj, celltype_list, condition_pairs, condition_shortnames, 
                                          group_by1, group_by2, assay = "RNA", output_dir, 
                                          dea_method = c("MAST", "pseudobulk")) {
  
  dea_method <- match.arg(dea_method)
  
  # Iterate over cell types and condition pairs
  for (celltype in celltype_list) {
    skip_celltype <- TRUE
    
    full_name <- celltype$full_name
    short_name <- celltype$short_name
    
    # Subset the Seurat object by cell type
    celltype_subset <- subset(seuratObj, cells = which(seuratObj[[group_by1]] == full_name))
    
    # Count cells per sample
    sample_counts <- table(celltype_subset$sample)
    
    # Keep only samples with at least 20 cells
    valid_samples <- names(sample_counts[sample_counts >= 20])
    valid_cells <- rownames(celltype_subset@meta.data[celltype_subset$sample %in% valid_samples, ])
    celltype_subset <- subset(celltype_subset, cells = valid_cells)
    
    # Create directories
    celltype_folder <- file.path(output_dir, short_name)
    dir.create(celltype_folder, recursive = TRUE, showWarnings = FALSE)
    figures_folder <- file.path(celltype_folder, "Figures")
    dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Create workbooks for DEG and GSEA results
    wb_deg <- createWorkbook()
    wb_gsea <- createWorkbook()
    
    for (condition_pair in condition_pairs) {
      condition1 <- condition_pair[[1]]
      condition2 <- condition_pair[[2]]
      
      # Count samples per condition
      valid_meta <- celltype_subset@meta.data
      samples_per_condition <- table(valid_meta[[group_by2]], valid_meta$sample)
      
      # Ensure both conditions are present
      if (!(condition1 %in% rownames(samples_per_condition)) ||
          !(condition2 %in% rownames(samples_per_condition))) {
        message(paste("Skipping", short_name, "-", condition1, "vs", condition2, 
                      "- one or both conditions missing"))
        next
      }
      
      # Number of valid samples per condition (with ≥ 20 cells)
      condition1_sample_count <- sum(samples_per_condition[condition1, ] >= 20)
      condition2_sample_count <- sum(samples_per_condition[condition2, ] >= 20)
      
      if (condition1_sample_count < 1 || condition2_sample_count < 1) {
        message(paste("Skipping", short_name, "-", condition1, "vs", condition2, 
                      "- not enough valid samples (≥1 required per condition)."))
        next
      }
      
      
      short_condition1 <- condition_shortnames[[condition1]]
      short_condition2 <- condition_shortnames[[condition2]]
      
      if (dea_method == "MAST") {
        # Perform single-cell DEA using MAST
        degs_df <- FindMarkers(
          celltype_subset, group.by = group_by2, assay = assay, verbose = TRUE,
          ident.1 = condition1, ident.2 = condition2, test.use = 'MAST', min.pct = 0.1
        )
      } else {
        # Pseudobulk DEA
        seuratObj_sub <- celltype_subset
        threshold <- 15
        metadata <- seuratObj_sub@meta.data
        conditions <- unique(metadata[[group_by2]])
        original_genes <- rownames(seuratObj_sub)
        genes_to_keep <- c()
        
        for (cond in conditions) {
          cells_in_condition <- rownames(metadata[metadata[[group_by2]] == cond, ])
          gene_pct <- Matrix::rowMeans(seuratObj_sub@assays$RNA$counts[, cells_in_condition] > 0) * 100
          genes_to_keep <- unique(c(genes_to_keep, names(gene_pct[gene_pct >= threshold])))
        }
        
        genes_to_keep <- intersect(genes_to_keep, rownames(seuratObj_sub))
        y <- Seurat2PB(seuratObj_sub, sample = 'sample', cluster = group_by2)
        keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=5, min.total.count=15)
        y <- y[keep.genes, , keep = F]
        keep.genes.manual <- rownames(y) %in% genes_to_keep
        y <- y[keep.genes.manual, , keep = FALSE]
        y <- calcNormFactors(y, method = "TMM")
        sample <- as.factor(y$samples$sample)
        group <- as.factor(y$samples$cluster)
        design <- model.matrix(~ group)
        colnames(design)[1] <- "Int"
        y <- estimateDisp(y, design, robust=TRUE)
        fit <- glmQLFit(y, design, robust=TRUE)
        lrt <- glmLRT(fit, coef = 2)
        results <- topTags(lrt, n = Inf)
        degs_df <- as.data.frame(results)
        degs_df$gene <- rownames(degs_df)
        colnames(degs_df)[colnames(degs_df) == "logFC"] <- "avg_log2FC"
        colnames(degs_df)[colnames(degs_df) == "FDR"] <- "p_val_adj"
        colnames(degs_df)[colnames(degs_df) == "PValue"] <- "P.Value"
      }
      
      skip_celltype <- FALSE
      
      # Add gene column to the results
      degs_df$gene <- rownames(degs_df)
      degs_df <- degs_df[, c("gene", setdiff(names(degs_df), "gene"))]  
      
      if (dea_method == "pseudobulk") {
        degs_df <- degs_df[, c("gene", "avg_log2FC", "logCPM", "P.Value", "p_val_adj")]
      }
      
      # Filter DEGs by log2 fold change > 1
      degs_filtered <- subset(degs_df, abs(avg_log2FC) > 1)
      
      # Write DEGs and filtered DEGs to the workbook
      sheet_name_deg <- paste(short_condition1, "v", short_condition2, sep = "_")
      addWorksheet(wb_deg, sheet_name_deg)
      writeDataTable(wb_deg, sheet = sheet_name_deg, x = degs_df, tableStyle = "TableStyleMedium9")
      
      filtered_sheet_name <- paste(sheet_name_deg, "logFC_filtered", sep = "_")
      addWorksheet(wb_deg, filtered_sheet_name)
      writeDataTable(wb_deg, sheet = filtered_sheet_name, x = degs_filtered, tableStyle = "TableStyleMedium9")
      
      # Generate a volcano plot of non-filtered DEGs
      volcano_plot <- EnhancedVolcano(
        degs_df,
        lab = degs_df$gene,
        x = 'avg_log2FC',
        y = 'p_val_adj',
        title = paste(short_name, short_condition1, "vs", short_condition2),
        pCutoff = 0.05, # Adjust the p-value threshold as needed
        FCcutoff = 1,   # Fold-change cutoff
        xlim = c(-3, 3), # Customize x-axis limits
        colAlpha = 0.8,
        legendLabels = c("NS", "Log2FC", "Adj.P", "Adj.P & Log2FC"),
        legendPosition = 'bottom'
      )
      
      # Save the volcano plot
      volcano_plot_file <- file.path(figures_folder, paste0(sheet_name_deg, "_volcano_plot.png"))
      ggsave(volcano_plot_file, plot = volcano_plot, width = 10, height = 14)
      
      
      # Prepare the gene list for GSEA (filtered genes only)
      gene_symbols <- degs_filtered$gene
      avg_log2FC <- degs_filtered$avg_log2FC
      gene_list <- setNames(avg_log2FC, gene_symbols)
      
      # Convert gene symbols to Entrez IDs
      entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = "ENTREZID", keytype = "SYMBOL")
      
      # Merge Entrez IDs into the gene list
      gene_list_entrez <- sapply(names(gene_list), function(symbol) {
        entrez_id <- entrez_ids$ENTREZID[entrez_ids$SYMBOL == symbol]
        if (length(entrez_id) > 0) {
          return(entrez_id)
        } else {
          return(NA)
        }
      })
      
      # Final gene list with Entrez IDs
      final_gene_list <- setNames(gene_list, gene_list_entrez)
      sorted_gene_list <- sort(final_gene_list, decreasing = TRUE)
      
      # Perform GSEA with GO (BP, MF, CC), WikiPathways, and Reactome
      # Possible to add limits for minimum and maximum gene set size
      gsea_results <- list()
      gsea_results[[paste(short_condition1, short_condition2, "GO_BP", sep = "_")]] <- gseGO(geneList = sorted_gene_list,
                                                                                             OrgDb = org.Hs.eg.db,
                                                                                             ont = "BP",
                                                                                             pvalueCutoff = 1,
                                                                                             verbose = FALSE,
                                                                                             scoreType = 'pos')
      gsea_results[[paste(short_condition1, short_condition2, "GO_MF", sep = "_")]] <- gseGO(geneList = sorted_gene_list,
                                                                                             OrgDb = org.Hs.eg.db,
                                                                                             ont = "MF",
                                                                                             pvalueCutoff = 1,
                                                                                             verbose = FALSE,
                                                                                             scoreType = 'pos')
      gsea_results[[paste(short_condition1, short_condition2, "GO_CC", sep = "_")]] <- gseGO(geneList = sorted_gene_list,
                                                                                             OrgDb = org.Hs.eg.db,
                                                                                             ont = "CC",
                                                                                             pvalueCutoff = 1,
                                                                                             verbose = FALSE,
                                                                                             scoreType = 'pos')
      gsea_results[[paste(short_condition1, short_condition2, "WP", sep = "_")]] <- gseWP(sorted_gene_list,
                                                                                          organism = "Homo sapiens",
                                                                                          pvalueCutoff = 1,
                                                                                          verbose = FALSE,
                                                                                          scoreType = 'pos')
      gsea_results[[paste(short_condition1, short_condition2, "RCTM", sep = "_")]] <- gsePathway(sorted_gene_list, 
                                                                                                 pvalueCutoff = 1,
                                                                                                 pAdjustMethod = "BH", 
                                                                                                 verbose = FALSE,
                                                                                                 scoreType = 'pos')
      # Perform MSigDB enrichment analyses for categories H, C1, C2, and C3
      msig_categories <- c("H", "C1", "C2", "C3")
      for (category in msig_categories) {
        msig_t2g <- msigdbr(species = "Homo sapiens", category = category) %>% 
          dplyr::select(gs_name, entrez_gene)
        gsea_results[[paste(short_condition1, short_condition2, paste0("MSigDB_", category), sep = "_")]] <- 
          GSEA(sorted_gene_list, TERM2GENE = msig_t2g)
      }
      
      # Save GSEA results and create enrichment plots
      for (gsea_name in names(gsea_results)) {
        gsea_result <- gsea_results[[gsea_name]]
        
        if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result@result)) > 0) {
          # Convert the `core_enrichment` column into a list of Entrez IDs
          core_entrez_ids <- strsplit(gsea_result@result$core_enrichment, "/")
          
          # Function to convert Entrez IDs to gene symbols
          convert_to_symbols <- function(entrez_ids) {
            symbols <- unlist(mget(entrez_ids, org.Hs.egSYMBOL, ifnotfound = NA))
            return(symbols)
          }
          
          # Apply the conversion to each list of Entrez IDs
          gsea_result@result$core_enrichment <- lapply(core_entrez_ids, convert_to_symbols)
          
          # Replace the `core_enrichment` column in the dataframe with gene symbols
          gsea_result@result$core_enrichment <- sapply(gsea_result@result$core_enrichment, paste, collapse = "/")
          
          # Save GSEA result to a workbook
          gsea_result_df <- as.data.frame(gsea_result@result)
          addWorksheet(wb_gsea, gsea_name)
          writeDataTable(wb_gsea, sheet = gsea_name, x = gsea_result_df, tableStyle = "TableStyleMedium9")
          
          # Create an enrichment plot
          plot_title <- paste(short_name, gsea_name, sep = ": ")
          enrichment_plot <- dotplot(gsea_result, showCategory = 20) + ggtitle(plot_title)
          plot_file <- file.path(figures_folder, paste0(gsea_name, "_enrichment_plot.png"))
          ggsave(plot_file, plot = enrichment_plot, width = 8, height = 12, units = "in", dpi = 100)
        }
      }
    }
    
    # Save the DEG and GSEA workbooks
    deg_output_file <- file.path(celltype_folder, paste(short_name, "DEG_results.xlsx", sep = "_"))
    gsea_output_file <- file.path(celltype_folder, paste(short_name, "GSEA_results.xlsx", sep = "_"))
    
    if (skip_celltype) {
      message(paste("Deleting output folder for", short_name, "since all comparisons were skipped."))
      unlink(celltype_folder, recursive = TRUE, force = TRUE)
    } else {
      # Save DEG workbook
      saveWorkbook(wb_deg, file.path(celltype_folder, paste0(short_name, "_DEG_results.xlsx")), overwrite = TRUE)
      saveWorkbook(wb_gsea, file.path(celltype_folder, paste0(short_name, "_GSEA_results.xlsx")), overwrite = TRUE)
      cat("DEG and GSEA analysis completed for", full_name, ". Results saved to:", celltype_folder, "\n")
    }
  }
}

# Define cell types to analyze
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



# Define condition pairs for differential expression analysis
condition_pairs <- list(
  c("3-OXC2-Mig-ace-leu", "9-WTC-Mig-ace-leu")
)

# Define short names for conditions
condition_shortnames <- list(
  "OXC" = "OXC",
  "WTC" = "WTC"
)

# Define output directory
output_dir <- "./GSEA results/OXC2 combined vs WTC combined"

# Run the analysis using pseudobulk
perform_DEG_and_GSEA_analysis(
  seuratObj = seuratObj,
  celltype_list = celltype_list,
  condition_pairs = condition_pairs,
  condition_shortnames = condition_shortnames,
  group_by1 = 'celltype',
  group_by2 = "sample",
  assay = "RNA",
  output_dir = output_dir,
  dea_method = "MAST"
)

OXC_combined_sample_subset <- subset(seuratObj, sample == "3-OXC2-Mig-ace-leu")
table(OXC_combined_sample_subset$celltype)

WTC_combined_sample_subset <- subset(seuratObj, sample == "9-WTC-Mig-ace-leu")
table(WTC_combined_sample_subset$celltype)
