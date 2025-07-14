library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

file_meta <- "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_calls.tsv"
metafusion <- fread(file_meta, sep = "\t")

# CPTAC
## info
CPTAC_sheet <- fread("/data/data2/Projects/Fusion/unclean/CPTAC_gdc_sample_sheet.2023-11-15.tsv", sep = "\t")

### Composite Sample 
###
### Composite Sample use comma to separate like C3L-03400-05, C3L-03400-02
### 
### It leads one question is whether these composite samples using comma to separate from identical patient?
rows_with_comma <- grepl(", ", CPTAC_sheet$`Sample ID`)
rows_with_comma_data <- CPTAC_sheet[rows_with_comma, ]

check_same_chars <- function(x) {
  parts <- unlist(strsplit(x, ", "))
  return(parts[1] == parts[2])
}
same_chars <- sapply(rows_with_comma_data$`Case ID`, check_same_chars)
if (all(same_chars)) {
  cat("TRUE\n")
} else {
  cat("FALSE\n")
}
### They are from identical patient.

### For the TARGET databases, check whether there is any overlap between composite and normal Sample_ID in the info data. If there is an overlap, retain only the normal Sample_ID. 
check_cptac_overlap <- CPTAC_sheet |>  
                       select(c(`Case ID`, `Sample ID`, `Sample Type`, `Project ID`)) |>
                       unique() |> 
                       separate_rows(`Sample ID`, sep = ",\\s*") 
 
check_cptac_overlap <- check_cptac_overlap |>
                        mutate(Sample_Type = ifelse(grepl(", ", `Sample ID`), "Composite", "Normal"))
                      
overlap_samples <- check_cptac_overlap |>
                        group_by(`Sample ID`) |>
                        summarise(Count = n(), .groups = 'drop') |>
                        filter(Count > 1)
### > overlap_samples 
### A tibble: 0 × 2
### ℹ 2 variables: Sample ID <chr>, Count <int> 
### no this situation

### Due to the less rigorous patient naming conventions(like 604) in the CPTAC2 project, check for any overlapping patient names between CPTAC2 and CPTAC3.
CPTAC_sheet_rigorous  <- CPTAC_sheet |> 
                         select(c(`Case ID`, `Sample ID`, `Sample Type`, `Project ID`)) |>  
                         unique() 
duplicates <- CPTAC_sheet_rigorous |>
                          group_by(`Case ID`, `Sample ID`) |>
                          summarise(n_projects = n(), projects = paste(`Project ID`, collapse = ", ")) |>
                          filter(n_projects > 1)
### nrow(duplicates)
### [1] 0
### no same patient names between CPTAC2 and CPTAC3

### CPTAC_sheet pretreatment
matched_CPTAC_pretreatment <- CPTAC_sheet |> 
  select(c(`Case ID`, `Sample ID`, `Sample Type`, `Project ID`)) |>
  rename("Patient_ID" = `Case ID`,
         "Sample_ID" = `Sample ID`,
         "Sample_type"  =  `Sample Type`,
         "Project_id" = `Project ID`) |> 
  separate_rows(Patient_ID, sep = ",\\s*") |>  #because they are from identical patient.
  unique() 

### CPTAC_info pretreatment
CPTAC_info <- fread("/data/data2/Projects/Fusion/fusiondb/Clininfo/CPTAC_info.tsv", sep = "\t")
# > CPTAC_info$Cancer |> unique()
# [1] "BRCA"       "CCRCC"      "COAD"       "GBM"        "HGSC"       "HNSCC"      "LSCC"       "LUAD"       "PDAC"       "UCEC"       "GTEX-Brain"
### We use our naming standard mutate a new column cohort.
cptac_replacement_table2 <- c(
  "BRCA" = "CPTAC-BRCA",
  "CCRCC" = "CPTAC-KIRC",
  "COAD" =  "CPTAC-COAD",
  "GBM" = "CPTAC-GBM",
  "HGSC" = "CPTAC-OV",
  "HNSCC" = "CPTAC-HNSC",
  "LSCC" = "CPTAC-LUSC" ,
  "LUAD" = "CPTAC-LUAD",
  "PDAC" = "CPTAC-PAAD",
  "UCEC" = "CPTAC-UCEC",
  "GTEX-Brain" = "CPTAC-GTEX-Brain"
)     

matched_CPTAC_info <- CPTAC_info$Cancer %in% names(cptac_replacement_table2)
CPTAC_info_pretreatment <- CPTAC_info %>%
  mutate(cohort = case_when(
    matched_CPTAC_info ~ cptac_replacement_table2[Cancer],
    TRUE ~ as.character(NA)  
  )) |> 
  select(-Sample_ID)
### CPTAC_info_pretreatment  |> filter(is.na(cohort)) |> nrow()
### [1] 0

matched_CPTAC <- matched_CPTAC_pretreatment |>
                 left_join(CPTAC_info_pretreatment, by = c("Patient_ID" = "Patient_ID")) |>
                 mutate(Run = str_replace_all(Sample_ID, ", ", "-"))

### CPTAC Sample type
### matched_CPTAC_pretreatment$Sample_type |> unique()
### [1] "Primary Tumor"                                                             "Solid Tissue Normal"                                                      
### [3] "Primary Tumor, Primary Tumor"                                              "Primary Tumor, Primary Tumor, Primary Tumor"                              
### [5] "Solid Tissue Normal, Solid Tissue Normal, Solid Tissue Normal"             "Solid Tissue Normal, Solid Tissue Normal"                                 
### [7] "Primary Tumor, Primary Tumor, Primary Tumor, Primary Tumor, Primary Tumor" "Primary Tumor, Primary Tumor, Primary Tumor, Primary Tumor"   
### 
### Just two sample type so if grepling "Solid Tissue Normal" is normal sample, grepling "Primary Tumor" is tumor sample
matched_CPTAC_after_processing <- matched_CPTAC |>
  mutate(cohort = case_when(
    grepl("Solid Tissue Normal", Sample_type) ~ "CPTAC-NORMAL",
    TRUE ~ cohort  
  )) 

### make sure nrow is same between Before and After Processing. 
matched_CPTAC_after_processing |> nrow()
### [1] 2558
matched_CPTAC_pretreatment |> nrow() 
### [1] 2558

### make sure IOBR/TME data just contain sample which cptac_sheet gave. 
detected_id_CPTAC <- matched_CPTAC_after_processing$Run |> unique()

### save info
cohorts <- unique(matched_CPTAC_after_processing$cohort)

for (cohort1 in cohorts) {
  cohort_data <- subset(matched_CPTAC_after_processing, cohort == cohort1)
  name <- gsub("-", "_", cohort1)
  file_name <- paste0("/data/data2/Projects/Fusion/fusiondb_v4/Clininfo/", name ,"_info.tsv")
  fwrite(cohort_data , file = file_name, sep = "\t")  
}

## TME
CPTAC2_TME <- fread("/data/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-2_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC2")
CPTAC3_TME <- fread("/data/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-3_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC3")
CPTAC_TME <- rbind(CPTAC2_TME, CPTAC3_TME)

### Match matched_CPTAC_after_processing$Run to CPTAC_TME$ID
### 
### CPTAC_TME can match
matched_can <- CPTAC_TME %>%
               filter(ID %in% matched_CPTAC_after_processing$Run) 

### 
### CPTAC_TME can't match
matched_not <- CPTAC_TME %>%
               filter(!(ID %in% matched_CPTAC_after_processing$Run))

### 
### situtation 1: Due to the barcode naming convention and etc, some portion character are added to the endning of the IOBR/TME sample name.
### str_count(matched_last$ID, "_")
### [1] 1 2 0
### 1.1 "CPT000817-0008" and "CPT000817". If the value of str_count(matched_last$ID, "_") for a row is 0, then the value of the prefixes column for that row should be the characters before the last -, and the value of the suffixes column should be the characters after the last -.
### 1.2 "f6a7f764-2495-4a1d-89a2-e12f51_D7_1" and "f6a7f764-2495-4a1d-89a2-e12f51". If the value of str_count(matched_last$ID, "_") for a row is 1, then the value of the prefixes column for that row should be the characters before the first _, and the value of the suffixes column should be the characters after the first _.
### 1.3 "CPT008756_0001" and "CPT008756". f6a7f764-2495-4a1d-89a2-e12f51_D7_1, .If the value of str_count(matched_last$ID, "_") for a row is 2, then the value of the prefixes column for that row should be the characters before the first _, and the value of the suffixes column should be the characters after the first _.

matched_not <- matched_not %>%
  rowwise() %>%  
  mutate(
    prefixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][1],  
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][1],  
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][1],  
      TRUE ~ NA_character_
    ),
    suffixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][length(str_split(ID, "-")[[1]])],  
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][2],  
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][2],  
      TRUE ~ NA_character_
    ),
    name1 = prefixes
  ) %>%
  ungroup() 

### situtation 1 after processing can match
matched_in_situtation1 <-  matched_not %>%
                           filter((name1 %in% matched_CPTAC_after_processing$Run))%>%
                           mutate(ID = name1)  %>%
                           select(-c(suffixes,prefixes,name1) )     
matched_can <- rbind(matched_can, matched_in_situtation1)

### not in situtation 1 
matched_not_situtation1 <- matched_not %>%
                           filter(!(name1 %in% matched_CPTAC_after_processing$Run))

### situtation 2: Due to IOBR package analyze code, an additional "X" is added to the beginning of the IOBR/TME sample name.
### "X1104806" and "1104806", "X9f905736-f662-41d6-b3ac-16758d" and "9f905736-f662-41d6-b3ac-16758d". 
matched_not_situtation1  <- matched_not_situtation1 %>%
                            mutate(name2 = str_remove(ID, "^X"))  

### situtation 2 after processing can match
matched_in_situtation2 <- matched_not_situtation1 %>%
                          filter((name2 %in% matched_CPTAC_after_processing$Run))    %>%
                          mutate(ID = name2 )  %>%
                          select(-c(name1,prefixes,suffixes,name2)) 
matched_can <- rbind(matched_can, matched_in_situtation2)

### situtation 2 after processing can't match
matched_not_situtation2 <- matched_not_situtation1 %>%
                           filter(!(name2 %in% matched_CPTAC_after_processing$Run))

### situtation 3: Mixed situation1 and situation2.
### X161317f3-1a39-438e-af45-cd2502_D7 and 161317f3-1a39-438e-af45-cd2502
matched_not_situtation2  <- matched_not_situtation2 %>%
                            mutate(name3 = str_remove(prefixes, "^X")) 

### situtation 3 after processing can match    
matched_in_situtation3 <- matched_not_situtation2  %>%
                          filter((name3 %in% matched_CPTAC_after_processing$Run))  %>%
                          mutate(ID = name3) %>%
                          select(-c(name1,name2,prefixes,suffixes, name3))           
matched_can <- rbind(matched_can, matched_in_situtation3)

### situtation 3 after processing can't match
matched_not_situtation3 <- matched_not_situtation2  %>%
                           filter(!(name3 %in% matched_CPTAC_after_processing$Run))
nrow(matched_not_situtation3 )
### [1] 81

### Check the reason why it still can't match after Type 3 processing
### After manual review, these Sample_ID are not included in either the info dataset or the gdc_sheet.

### Check if CPTAC TME data ID are any samples that have been divided into multiple portions.
matched_can |> select(ID) |>  nrow()
### [1] 2558
matched_can |> select(ID) |>  unique() |> nrow()
### [1] 2558
### No this situation.

### Check if CPTAC TME data has same norw with CPTAC info data.
matched_can |> select(ID) |>   unique() |> nrow()
### [1] 2558
matched_CPTAC_after_processing |> select(Run) |>  unique() |> nrow()
### [1] 2558
### Sure.

cohort_data <- matched_CPTAC_after_processing |> 
               select(Run, cohort)   |>
               unique()   

matched_can <- matched_can  |>
               mutate(Run = ID)  |>
               left_join(cohort_data, by = "Run")  |>
               select(-Run)

cohorts_TME <- unique(matched_can$cohort)

for (cohort1 in cohorts_TME) {
  cohort_data <- subset(matched_can, cohort == cohort1) |> select(-cohort)
  file_name <- paste0("/data/data2/Projects/Fusion/fusiondb_v4/IOBR/", cohort1 ,"_TME.tsv.gz")
  write.table(cohort_data, 
              file = gzfile(file_name), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE) 
}  

#IOBR
CPTAC2_Signature <- fread("/data/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-2_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC2")
CPTAC3_Signature <- fread("/data/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-3_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC3")
CPTAC_Signature <- rbind(CPTAC2_Signature, CPTAC3_Signature)

### Match matched_CPTAC_after_processing$Run to CPTAC_TME$ID
### 
### CPTAC_TME can match
matched_can <- CPTAC_Signature %>%
               filter(ID %in% matched_CPTAC_after_processing$Run) 


### 
### CPTAC_TME can't match
matched_not <-  CPTAC_Signature %>%
                filter(!(ID %in% matched_CPTAC_after_processing$Run))

### 
### situtation 1:like TME 
matched_not <- matched_not %>%
  rowwise() %>%  
  mutate(
    prefixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][1],  
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][1],  
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][1], 
      TRUE ~ NA_character_
    ),
    suffixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][length(str_split(ID, "-")[[1]])],  
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][2],  
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][2],  
      TRUE ~ NA_character_
    ),
    name1 = prefixes
  ) %>%
  ungroup()  

### situtation 1 after processing can match
matched_in_situtation1 <- matched_not %>%
                          filter((name1 %in% matched_CPTAC_after_processing$Run))%>%
                          mutate(ID = name1)  %>%
                          select(-c(suffixes,prefixes,name1) )    

matched_can <- rbind(matched_can, matched_in_situtation1)

### situtation 1 after processing can't match
matched_not_situtation1 <- matched_not %>%
                           filter(!(name1 %in% matched_CPTAC_after_processing$Run))

### situtation 2:like TME 
matched_not_situtation1  <- matched_not_situtation1 %>%
                            mutate(name2 = str_remove(ID, "^X"))  

### situtation 2 after processing can match
matched_in_situtation2 <- matched_not_situtation1 %>%
                          filter((name2 %in% matched_CPTAC_after_processing$Run))    %>%
                          mutate(ID = name2 )  %>%
                          select(-c(name1,prefixes,suffixes,name2)) 

matched_can <- rbind(matched_can, matched_in_situtation2)

### situtation 2 after processing can't match
matched_not_situtation2 <- matched_not_situtation1 %>%
                           filter(!(name2 %in% matched_CPTAC_after_processing$Run))

### situtation 3:like TME 
matched_not_situtation2   <- matched_not_situtation2  %>%
                             mutate(name3 = str_remove(prefixes, "^X")) 

### situtation 3 after processing can match   
matched_in_situtation3 <- matched_not_situtation2  %>%
                          filter((name3 %in% matched_CPTAC_after_processing$Run))  %>%
                          mutate(ID = name3) %>%
                          select(-c(name1,name2,prefixes,suffixes, name3))     

matched_can <- rbind(matched_can, matched_in_situtation3)

### situtation 3 after processing can't match  
matched_not_situtation3 <- matched_not_situtation2  %>%
                           filter(!(name3 %in% matched_CPTAC_after_processing$Run))

nrow(matched_not_situtation3 )
### [1] 81
                          
### Check the reason why it still can't match after Type 3 processing
### After manual review, these Sample_ID are not included in either the info dataset or the gdc_sheet.

### Check if CPTAC IOBR data ID are any samples that have been divided into multiple portions.
matched_can |> select(ID) |>  nrow()
### [1] 2558
matched_can |> select(ID) |>  unique() |> nrow()
### [1] 2558
### No this situation.

### Check if CPTAC IOBR data has same norw with CPTAC info data.
matched_can |> select(ID) |>   unique() |> nrow()
### [1] 2558
matched_CPTAC_after_processing |> select(Run) |>  unique() |> nrow()
### [1] 2558
### Sure.

cohort_data   <- matched_CPTAC_after_processing |> 
                 select(Run, cohort)   |>
                 unique()   

matched_can <- matched_can |>
               mutate(Run = ID)  |>
               left_join(cohort_data, by = "Run")  |>
               select(-Run)

cohorts_IOBR <- unique(matched_can$cohort)

for (cohort1 in cohorts_IOBR) {
  cohort_data <- subset(matched_can, cohort == cohort1) |> select(-cohort)
  file_name <- paste0("/data/data2/Projects/Fusion/fusiondb_v4/IOBR/", cohort1 ,"_IOBR_signature.tsv.gz")
  write.table(cohort_data, 
              file = gzfile(file_name), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE) 
}  
 
## Gene

metafusion1 <- metafusion |> 
  mutate(
    Sample_ID = case_when(
      cohort == "CPTAC-3" & grepl("_", Sample_ID) ~  str_replace_all(Sample_ID, "_", "-"),
      TRUE ~ Sample_ID 
    ),
    cohort = case_when(
      Sample_ID %in% matched_CPTAC_after_processing$Run  & grepl("CPTAC", cohort) ~ matched_CPTAC_after_processing$cohort[match(Sample_ID, matched_CPTAC_after_processing$Run)],
      TRUE ~ cohort
    ) 
  )  

### In the past, the metafusion file had an issue with the sample_type being incorrectly labeled.
### After further investigation, it was found that only the composite Sample_ID in CPTAC had this issue, where samples labeled as 'Normal' were mistakenly identified as 'Tumor' in the sample_type column.

normal_id <-   matched_CPTAC_after_processing |> 
               filter(cohort == "CPTAC-NORMAL") |> 
               select(Run) |>
               unique() |>
               pull(Run)

meta_normal <- metafusion1 |>
               filter(Sample_ID %in% normal_id)

unique(meta_normal$sample_type)
#[1] "Normal" "Tumor" 

meta_processing  <- metafusion1 |>
                    mutate(
                      sample_type = ifelse(
                        sample_type != "Normal" & Sample_ID %in% normal_id,
                        "Normal",
                        sample_type
                      )
                    )

cp_file_name1 <- "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_calls_matched.cptac.tsv"
fwrite(meta_processing, file = cp_file_name1, sep = "\t") 




