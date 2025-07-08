library("scProportionTest")
library(writexl)

prop_test <- sc_utils(seuratObj)

prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WTC", sample_2 = "OXC2",
  sample_identity = "group"
)
permutation_plot(prop_test)

# Getting celltype distribution in each sample ----
  # this is not the same as celltype abundancy in the permutation, but an absolute measurment

celltype_per_sample <- seuratObj@meta.data %>%
  group_by(sample, celltype) %>%
  summarise(Cellcount = n(), .groups = "drop")


celltype_per_sample$celltype <- factor(celltype_per_sample$celltype, levels = celltype_order)
sample_order <- c(
  "11-WTC-control",
  "6-WTC-Miglustat-100uM",
  "8-WTC-Acetyl-leucine",
  "9-WTC-Mig-ace-leu",
  "4-OXC2-DMSO-control",
  "1-OXC2-Miglustat-100uM",
  "2-OXC2-Acetyl-leucine",
  "3-OXC2-Mig-ace-leu"
)

celltype_per_sample$sample <- factor(celltype_per_sample$sample, levels = sample_order)

ggplot(celltype_per_sample, aes(x = celltype, y = Cellcount, fill = celltype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sample, scales = "free_y", ncol = 4, nrow = 2) +
  theme_minimal() +
  labs(title = "Celltype Composition per Sample",
       x = "Celltype", y = "Cell Count") +
  theme(axis.text.x = element_blank())
