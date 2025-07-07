

significantGSEAs <- list.files(path = "./GSEA results/WTC Combined vs WTC DMSO/AllSignificantGSEAs", 
                               pattern = "\\.xlsx$", full.names = TRUE)


# function to read and filter sheets in excel files
read_and_filter_sheet <- function(file, sheet) {
  df <- read_excel(file, sheet = sheet)
  
  # filter the results by adjusted p-value and selected columns
  df_filtered <- df %>%
    filter(p.adjust <= 0.05) %>%
    select(c(ID, Description, setSize, NES, p.adjust, core_enrichment))
  
  return(df_filtered)
}

# use the filter function to extract the significant results 
filtered_results <- map(significantGSEAs, function(file) {
  sheets <- excel_sheets(file) # get names of all sheets in the Excel file
  sheet_data <- map(sheets, ~ read_and_filter_sheet(file, .x)) # apply filter function to each sheet
  names(sheet_data) <- sheets # name each element in the list by its sheet name
  return(sheet_data)
})
names(filtered_results) <- basename(significantGSEAs) # name each top-level list item by the file name

head(filtered_results)

# combine every data frame (excel sheet) into one
# containing every significant result from GSEA
all_results_named <- map2_dfr(filtered_results, names(filtered_results), function(sheet_list, file_name) {
  map2_dfr(sheet_list, names(sheet_list), function(df, sheet_name) {
    if (nrow(df) > 0) {
      df %>%
        mutate(Celltype = file_name, Source = sheet_name) %>%
        relocate(Celltype, Source)  
    }
  })
})

write_xlsx(all_results_named, "filtered_GSEA_results.xlsx")

head(all_results_named)
