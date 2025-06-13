library(Seurat)
library(readxl)

seuratObj <- readRDS('oxcwtc_seurat.RDS')
DimPlot(seuratObj, group.by = 'celltype', label = TRUE)
FeaturePlot(seuratObj, c("SP100", "STAT1", "IFI16"), ncol = 1, split.by = "group", order = TRUE)

#### Metadata ####
# Create metadata dataframe

sample_metadata <- data.frame(
  original_sample = c(
             "HCMV-1", 
             "HCMV-2", 
             "HCMV-3",
             "____-1",
             "MOCK-2",
             "MOCK-3"),
  
  sample = c("HCMV1", 
             "HCMV2", 
             "HCMV3",
             "MOCK1",
             "MOCK2",
             "MOCK3"),
  
  group = c("HCMV", "HCMV", "HCMV", 
            "Control", "Control", "Control")
  
)

#rename sample
seuratObj$sample <- sample_metadata$sample[match(seuratObj$sample, sample_metadata$original_sample)]

# Add metadata to Seurat object
seuratObj$group <- sample_metadata$group[match(seuratObj$sample, sample_metadata$sample)]

#### ####
DimPlot(seuratObj, group.by = 'seurat_clusters', label = TRUE, repel = T)

#### Create and list markers ####
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seuratObj <- CellCycleScoring(seuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(seuratObj, group.by = 'Phase', label = F)

#### Store and view markers (can group by clusters or celltype ####
seuratObj <- PrepSCTFindMarkers(seuratObj)

# create a group in the seuratObj containing markers for each specific cluster. 
# Can compare markers present in one cluster (ident.1) compared to other clusters (ident.2)
markers <- FindMarkers(seuratObj, group.by = 'celltype',ident.1 = "Viral-infected neural progenitor cells", only.pos = T, recorrect_umi = F)

# create list of marker names from the list created based on clusters
writeLines(rownames(head(markers, n = 25)))

#### Annotation ####
annotationdf <- read_excel("cluster annotation.xlsx")
seuratObj$celltype <- annotationdf$celltype[match(seuratObj$seurat_clusters, annotationdf$cluster)]
DimPlot(seuratObj, group.by = "celltype", label = T, repel = T)
DotPlot(seuratObj, group.by = "celltype", features = "SLC17A7")

DimPlot(seuratObj, group.by = "celltype", split.by = "sample", label = F, ncol = 3)

saveRDS(seuratObj, 'OXCT_annotated.RDS')

#### Subset for Mock samples ####
mock_samples <- subset(seuratObj, group == "Control")
DimPlot(mock_samples, group.by = "celltype", label = F, repel = T)
DimPlot(mock_samples, split.by = "sample", label = T)
FeaturePlot(mock_samples, c("SP100", "STAT1", "IFI16"), ncol = 1, split.by = "sample", order = TRUE)

mock_markers <- FindMarkers(mock_samples, group.by = 'celltype',
                            ident.1 = "Viral-infected neural progenitor cells", 
                            only.pos = T, recorrect_umi = F,
                            min.diff.pct = 0.2)
writeLines(rownames(head(mock_markers, n = 25)))


#### Subset for HCMV samples ####
HCMV_samples <- subset(seuratObj, group == "HCMV")
FeaturePlot(HCMV_samples, c("SP100", "STAT1", "IFI16"), ncol = 1, split.by = "sample", order = TRUE)


