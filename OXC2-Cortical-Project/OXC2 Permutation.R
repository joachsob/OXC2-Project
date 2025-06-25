library("scProportionTest")

prop_test <- sc_utils(seuratObj)

prop_test <- permutation_test(
  prop_test, cluster_identity = "sctype_classification",
  sample_1 = "Control", sample_2 = "Patient",
  sample_identity = "condition"
)
permutation_plot(prop_test)