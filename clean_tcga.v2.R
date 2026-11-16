library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
file_meta <- "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_calls_matched.cptac.target.tsv"
metafusion <- fread(file_meta, sep = "\t")

# TCGA
## info
TCGA_info <- fread("/data/data2/Projects/Fusion/fusiondb/Clininfo/TCGA_info.tsv", sep = "\t")
### TCGA_sheet pretreatment
TCGA_sheet <- fread("/data/data2/Projects/Fusion/unclean/TCGA_sample_sheet.2024-05-14.tsv", sep = "\t") |>
              select(c(`Sample ID`, `Sample Type`,`Project ID`))

colnames(TCGA_sheet) <- c("Sample_ID", "Sample_type","Project_id")
TCGA_sheet_pretreatment <- TCGA_sheet  |>
                           mutate(Sample_ID = substr(Sample_ID, 1, 15)) |>
                           unique()
### We use our naming standard mutate a new column cohort.
TCGA_after_processing <- TCGA_sheet_pretreatment |>
  left_join(TCGA_info, by = c("Sample_ID" = "Sample_ID")) |>
  mutate(
    cohort = Project_id,
    Run = Sample_ID,
    Patient_ID = substr(Sample_ID, 1, 12),
  )
### TCGA Sample type
### TCGA_after_processing$Sample_type |> unique()
###   [1] "Primary Tumor"                                   "Solid Tissue Normal"                             "Metastatic"
###   [4] "Additional - New Primary"                        "Primary Blood Derived Cancer - Peripheral Blood" "Additional Metastatic"
###   [7] "Recurrent Tumor"
###
### "Solid Tissue Normal" is normal sample, others is tumor sample

TCGA_NORMAL <- TCGA_after_processing |>
               filter(Sample_type == "Solid Tissue Normal")

NORMAL_id <- TCGA_NORMAL$Run |> unique()

TCGA_BLOOD <- TCGA_after_processing |>
               filter(Sample_type == "Primary Blood Derived Cancer - Peripheral Blood")

BLOOD_id <- TCGA_BLOOD$Run |> unique()

TCGA_after_processing <- TCGA_after_processing |>
  mutate(
    Project_id = cohort,
    cohort = case_when(
      Run %in% NORMAL_id ~ "TCGA-NORMAL",
      TRUE ~ cohort))

### make sure nrow is same between Before and After Processing.
TCGA_after_processing |> nrow()
### [1] 10545
TCGA_sheet_pretreatment |> nrow()
### [1] 10545

### make sure IOBR/TME data just contain sample which tcga_sheet gave.
detected_id_TCGA <- TCGA_after_processing$Run |> unique()

### save info
cohorts <- unique(TCGA_after_processing$cohort)

file_name1 <- paste0("/data/data2/Projects/Fusion/fusiondb_v4/Clininfo/")

for (cohort1 in cohorts) {
  cohort_data <- subset(TCGA_after_processing, cohort == cohort1)
  name <- cohort1
  file_name <- paste0(file_name1, name ,"_info.tsv")
  fwrite(cohort_data , file = file_name, sep = "\t")
}
## Gene

metafusion1 <- metafusion |>
  mutate(
    cohort = case_when(
      Sample_ID %in% TCGA_after_processing$Run & grepl("TCGA", cohort) ~ TCGA_after_processing$cohort[match(Sample_ID, TCGA_after_processing$Run)],
      TRUE ~ cohort
    ),
    sample_type = case_when(
      Sample_ID %in% NORMAL_id ~ "Normal_Adjacent",
      Sample_ID %in% BLOOD_id  ~ "Tumor_Blood" ,
      TRUE ~ sample_type
    ))
  
cp_file_name1 <- "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_gdc_divide.tsv"
fwrite(metafusion1, file = cp_file_name1, sep = "\t")

## IOBR
folder_path <- "/data/data2/Projects/Fusion/fusiondb/IOBR"

file_pattern <- "^TCGA.*_IOBR_signature.tsv.gz$"
file_list <- list.files(path = folder_path, pattern = file_pattern, full.names = TRUE)

combined_data <- data.table()

for (file in file_list) {

  file_data <- fread(file, sep = "\t")
  file_name <- basename(file)
  project_id <- sub("^(TCGA-[^_]+).*", "\\1", file_name)
  file_data <- file_data %>%
               mutate(Project_id = project_id)
  combined_data <- rbind(combined_data, file_data)
}

### Filter TCGA Replicate Samples
### The filter rule following broad institute says:
###
### In many instances there is more than one aliquot for a given combination of individual, platform, and data type. However, only one aliquot may be ingested into Firehose. Therefore, a set of precedence rules are applied to select the most scientifically advantageous one among them. Two filters are applied to achieve this aim: an Analyte Replicate Filter and a Sort Replicate Filter.
###
### Analyte Replicate Filter
### The following precedence rules are applied when the aliquots have differing analytes. For RNA aliquots, T analytes are dropped in preference to H and R analytes, since T is the inferior extraction protocol. If H and R are encountered, H is the chosen analyte. This is somewhat arbitrary and subject to change, since it is not clear at present whether H or R is the better protocol. If there are multiple aliquots associated with the chosen RNA analyte, the aliquot with the later plate number is chosen. For DNA aliquots, D analytes (native DNA) are preferred over G, W, or X (whole-genome amplified) analytes, unless the G, W, or X analyte sample has a higher plate number.
###
### Sort Replicate Filter
### The following precedence rules are applied when the analyte filter still produces more than one sample. The sort filter chooses the aliquot with the highest lexicographical sort value, to ensure that the barcode with the highest portion and/or plate number is selected when all other barcode fields are identical.
###
### Ref Link: <https://confluence.broadinstitute.org/display/GDAC/FAQ#FAQ-sampleTypesQWhatTCGAsampletypesareFirehosepipelinesexecutedupon>
###           <https://gist.github.com/ShixiangWang/33b2f9b49b77eaa8f773b428480f9101>


result_data <- combined_data |>
  mutate(
    prefixes = str_split(ID, "-", simplify = TRUE)[, 1:4] %>% apply(1, paste, collapse = "-"),
    suffixes = str_split(ID, "-", simplify = TRUE)[, 5],
    second_chars = str_extract(suffixes, "[0-9]+"),
    second_chars = as.integer(second_chars))
result_IOBR <- result_data  |>
  group_by(prefixes) |>
  slice(which.max(second_chars)) |>
  ungroup() |>
  mutate(cohort = Project_id,
         ID = substr(prefixes, 1, 15),
         cohort = case_when(ID %in% NORMAL_id ~ "TCGA-NORMAL",
                            TRUE ~ cohort)
  ) |>
  filter(ID %in% detected_id_TCGA)  |>
  select(-c("suffixes", "second_chars","prefixes")) |>
  unique()

cohorts_IOBR <- unique(result_IOBR$cohort)

for (cohort1 in cohorts_IOBR) {
  cohort_data <- subset(result_IOBR, cohort == cohort1) |> select(-c(cohort,Project_id) )
  file_name <- paste0("/data/data2/Projects/Fusion/fusiondb_v4/IOBR/", cohort1 ,"_IOBR_signature.tsv.gz")
  write.table(cohort_data,
              file = gzfile(file_name),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}

## TME

folder_path <- "/data/data2/Projects/Fusion/fusiondb/IOBR"
file_pattern <- "^TCGA.*_TME.tsv.gz$"
file_list <- list.files(path = folder_path, pattern = file_pattern, full.names = TRUE)

combined_data1 <- data.table()

for (file in file_list) {

  file_data <- fread(file, sep = "\t")

  file_name <- basename(file)
  project_id <- sub("^(TCGA-[^_]+).*", "\\1", file_name)

  file_data <- file_data %>%
    mutate(Project_id = project_id)

  combined_data1 <- rbind(combined_data1, file_data)
}

result_data1 <- combined_data1 |>
  mutate(
    prefixes = str_split(ID, "-", simplify = TRUE)[, 1:4] %>% apply(1, paste, collapse = "-"),
    suffixes = str_split(ID, "-", simplify = TRUE)[, 5],
    second_chars = str_extract(suffixes, "[0-9]+"),
    second_chars = as.integer(second_chars))
result_TME <- result_data1  |>
  group_by(prefixes) |>
  slice(which.max(second_chars)) |>
    ungroup() |>
      mutate(cohort = Project_id,
             ID = substr(prefixes, 1, 15),
             cohort = case_when(ID %in% NORMAL_id ~ "TCGA-NORMAL",
                                TRUE ~ cohort)
      ) |>
      filter(ID %in% detected_id_TCGA)  |>
      select(-c("suffixes", "second_chars","prefixes")) |>
      unique()

cohorts_TME <- unique(result_TME$cohort)

for (cohort1 in cohorts_TME) {
  cohort_data <- subset(result_TME, cohort == cohort1) |> select(-c(cohort,Project_id) )
  file_name <- paste0("/data/data2/Projects/Fusion/fusiondb_v4/IOBR/", cohort1 ,"_TME.tsv.gz")
  write.table(cohort_data,
              file = gzfile(file_name),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}

