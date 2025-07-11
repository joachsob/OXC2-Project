# COPY OF SCRIPT FOR SYNCING TO GITHUB

library(readxl)
library(dplyr)
library(purrr)
library(writexl)
library(tidyr)
library(stringr)

#

path = "./GSEA results/WTC Miglustat vs WTC DMSO/"           # Path where GSEA results are stored
gsea_output_file = "WTCmiglustat_WTCDMSO_filteredGSEA.xlsx"  # Desired file for all significant GSEA results
deg_output_file = file.path(path, "WTCmiglustat_WTCDMSO_filteredDEG.xlsx") # set output file for the filtered result


# Finding all GSEA files in the celltype sub-folders
all_excel_files <- list.files(path = path, pattern = "\\.xlsx$", full.names =  TRUE, recursive = TRUE) # list all excel files in GSEA results folder
gsea_files <- all_excel_files[grepl("GSEA", basename(all_excel_files), ignore.case = TRUE)] # list all GSEA excel files 
cell_types <- basename(dirname(gsea_files)) # list all celltypes in the GSEA results
names(cell_types) <- basename(gsea_files) # name each element in cell_types by celltype name for reference


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
filtered_results <- map(gsea_files, function(file) { #map() -> apply function to each element in gsea_files
  sheets <- excel_sheets(file) # get names of all sheets in the Excel file
  sheet_data <- map(sheets, ~ read_and_filter_sheet(file, .x)) # apply filter function to each sheet
  names(sheet_data) <- sheets # name each element in the list by its sheet name
  return(sheet_data)
})
names(filtered_results) <- basename(gsea_files) # name each top-level list item by the file name


# combine every data frame (excel sheet) into one
# containing every significant result from GSEA
gsea_results <- map2_dfr(filtered_results, names(filtered_results), function(sheet_list, file_name) {
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


# create one sheet in excel file for each celltype for managebility
gsea_split_by_celltype <- split(gsea_results, gsea_results$Celltype)
# clean sheet names cuz Excel is weird abt special characters
names(gsea_split_by_celltype) <- make.names(names(gsea_split_by_celltype))

# Write to file
write_xlsx(gsea_split_by_celltype, file.path(path, gsea_output_file))





#### Create file with all DEG results ####

# uses the code from first part - collecting all files containing DEG information
deg_files <- all_excel_files[grepl("DEG", basename(all_excel_files), ignore.case = TRUE)]
deg_cell_types <- basename(dirname(deg_files))
names(deg_cell_types) <- basename(deg_files)


# adjusted filtering function for DEG files
read_and_filter_deg_sheet <- function(file) {
  df <- read_excel(file)
  
  print(paste("Reading:", basename(file)))
  print(colnames(df))
  print(str(df$p.adjust))
  
  # If p.adjust is not numeric, convert or skip
  if (!"p_val_adj" %in% colnames(df) || !is.numeric(df$p_val_adj)) {
    message("Skipping: p.adjust missing or not numeric")
    return(tibble())  # Return empty tibble
  }
  
  # filter the results by adjusted p-value and selected columns
  df_filtered <- df %>%
    filter(p_val_adj <= 0.05) %>%
    select(c(gene, avg_log2FC, pct.1, pct.2, p_val_adj))
  
  return(df_filtered)
}

deg_results <- map2_dfr(deg_files, basename(deg_files), function(file, file_name) {
  df <- read_and_filter_deg_sheet(file)
  
  if (nrow(df) > 0) {
    df %>%
      mutate(
        Celltype = deg_cell_types[[file_name]],
      ) %>%
      relocate(Celltype)
  }
})

# create one sheet in excel file for each celltype for managebility
deg_split_by_celltype <- split(deg_results, deg_results$Celltype)
# clean sheet names cuz Excel is weird abt special characters
names(deg_split_by_celltype) <- make.names(names(deg_split_by_celltype))

# Write to file
write_xlsx(deg_split_by_celltype, deg_output_file)

