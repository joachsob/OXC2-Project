library(Seurat)
library(readxl)
 # THIS SCRIPT IS CURRENTLY USING HCMV SCRIPT AS TEMPLATE, AND IS NOT UPDATED FOR OXC2
seuratObj <- readRDS('oxcwtc_seurat.RDS')

treatment_order <- c(
  "11-WTC-control",
  "6-WTC-Miglustat-100uM",
  "8-WTC-Acetyl-leucine",
  "9-WTC-Mig-ace-leu",
  "4-OXC2-DMSO-control",
  "1-OXC2-Miglustat-100uM",
  "2-OXC2-Acetyl-leucine",
  "3-OXC2-Mig-ace-leu"
)
seuratObj$sample <- factor(seuratObj$sample, levels = treatment_order)

#### Adding metadata for groups for permutation plot ####
sample_metadata <- data.frame(
  group = c(
    "WTC", "WTC", "WTC","WTC",
    "OXC2", "OXC2", "OXC2", "OXC2"
  ),
  
  samples = c(
    "11-WTC-control",
    "6-WTC-Miglustat-100uM",
    "8-WTC-Acetyl-leucine",
    "9-WTC-Mig-ace-leu",
    "4-OXC2-DMSO-control",
    "1-OXC2-Miglustat-100uM",
    "2-OXC2-Acetyl-leucine",
    "3-OXC2-Mig-ace-leu"
  )
)

seuratObj$group <- sample_metadata$group[match(seuratObj$sample, sample_metadata$samples)]
table(seuratObj$group)

#### ####
DimPlot(seuratObj, label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, c("SP100", "STAT1", "IFI16"), ncol = 1, split.by = "group", order = TRUE)

DimPlot(seuratObj, group.by = 'celltype', label = TRUE,repel = T)

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



