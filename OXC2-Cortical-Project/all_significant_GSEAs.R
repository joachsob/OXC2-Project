library(readxl)
library(dplyr)
library(purrr)
library(writexl)
library(tidyr)
library(stringr)


path = "./GSEA results/WTC Miglustat vs WTC DMSO/"          # Path for GSEA results
gsea_output_file = "WTCmiglustat_WTCDMSO_filteredGSEA.xlsx"  # Desired file for all significant GSEA results
recurrent_genes_output_file = "WTCmiglustat_WTCDMSO_recurrent genes.xlsx" #Desired file for recurring genes

# Finding all GSEA files in the celltype sub-folders
all_excel_files <- list.files(path = path, pattern = "\\.xlsx$", full.names =  TRUE, recursive = TRUE)
gsea_files <- all_excel_files[grepl("GSEA", basename(all_excel_files), ignore.case = TRUE)]
cell_types <- basename(dirname(gsea_files))
names(cell_types) <- basename(gsea_files)


# function to read and filter sheets in excel files
read_and_filter_sheet <- function(file, sheet) {
  df <- read_excel(file, sheet = sheet)
  
  print(paste("Reading:", basename(file), "| Sheet:", sheet))
  print(colnames(df))
  print(str(df$p.adjust))
  
  # If p.adjust is not numeric, convert or skip
  if (!"p.adjust" %in% colnames(df) || !is.numeric(df$p.adjust)) {
    message("Skipping: p.adjust missing or not numeric")
    return(tibble())  # Return empty tibble
  }
  
  # filter the results by adjusted p-value and selected columns
  df_filtered <- df %>%
    filter(p.adjust <= 0.05) %>%
    select(c(ID, Description, NES, p.adjust, core_enrichment))
  
  return(df_filtered)
}

# use the filter function to extract the significant results 
filtered_results <- map(gsea_files, function(file) {
  sheets <- excel_sheets(file) # get names of all sheets in the Excel file
  sheet_data <- map(sheets, ~ read_and_filter_sheet(file, .x)) # apply filter function to each sheet
  names(sheet_data) <- sheets # name each element in the list by its sheet name
  return(sheet_data)
})
names(filtered_results) <- basename(gsea_files) # name each top-level list item by the file name


# combine every data frame (excel sheet) into one
  # containing every significant result from GSEA
all_results_named <- map2_dfr(filtered_results, names(filtered_results), function(sheet_list, file_name) {
  map2_dfr(sheet_list, names(sheet_list), function(df, sheet_name) {
    if (nrow(df) > 0) {
      df %>%
        mutate(
          Celltype = cell_types[[file_name]], 
          Source = sheet_name
          ) %>%
        relocate(Celltype, Source)  
    }
  })
})

# write all significant GSEA results to an Excel file in desired path/file
write_xlsx(all_results_named, file.path(path, gsea_output_file))


#### Create a list of the most recurring genes throughout the GSEA results ####

# select and split gene list
all_genes <- all_results_named %>% 
  select(core_enrichment) %>%                                     # select gene column
  filter(!is.na(core_enrichment)) %>%                             # Remove NAs
  mutate(core_enrichment = str_split(core_enrichment, "/")) %>%   # Split gene list into vectors
  unnest(core_enrichment) %>%                                     # Flatten the list into rows
  rename(Gene = core_enrichment)

gene_counts <- all_genes %>%
  group_by(Gene) %>%
  summarise(Appearances = n()) %>%
  arrange(desc(Appearances))

recurrent_genes <- gene_counts %>%
  filter(Appearances > 1)

write_xlsx(recurrent_genes, file.path(path, recurrent_genes_output_file))
