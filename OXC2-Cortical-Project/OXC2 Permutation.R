library("scProportionTest")

prop_test <- sc_utils(seuratObj)

prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WTC", sample_2 = "OXC2",
  sample_identity = "group"
)
permutation_plot(prop_test)
