library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

file_meta <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_calls_matched.cptac.tsv"
metafusion <- fread(file_meta, sep = "\t")

# TARGET
## info
### Because some projects are not included in the target_sheet, we should abandon them.
TARGET_sheet <- fread("/home/data2/Projects/Fusion/unclean/TARGET_gdc_sample_sheet.2023-11-15.tsv", sep = "\t")
project_choose <- substr(TARGET_sheet$`Sample ID`, 1, 9) |> unique()
project_choose
### "TARGET-20" "TARGET-00" "TARGET-21" "TARGET-30"

### TARGET_sheet pretreatment

### Composite Sample 
###
### Composite Sample use comma to separate like TARGET-30-PAPEAV, TARGET-30-PAPTFZ
TARGET_Composite <- TARGET_sheet |> 
  filter(grepl(", ", `Sample ID`)) |>  
  select(`Sample ID`) |>
  unique()
TARGET_Composite
### 1: TARGET-30-PAPEAV-01A, TARGET-30-PAPTFZ-01A
### 2: TARGET-30-PANUKV-01A, TARGET-30-PASUML-01A
### 3: TARGET-30-PANKFE-01A, TARGET-30-PAPUAR-01A
### 4: TARGET-30-PASYPX-01A, TARGET-30-PAIXIF-01A

### For the CPTAC databases, check whether there is any overlap between composite and normal Sample_ID in the info data. 
check_TARGET_overlap <- TARGET_sheet |>  
  select(c(`Case ID`, `Sample ID`, `Sample Type`, `Project ID`)) |>
  unique() |> 
  separate_rows(`Sample ID`, sep = ",\\s*") 

check_TARGET_overlap <- check_TARGET_overlap |>
  mutate(Sample_Type = ifelse(grepl(", ", `Sample ID`), "Composite", "Normal"))

overlap_samples <- check_TARGET_overlap |>
  group_by(`Sample ID`) |>
  summarise(Count = n(), .groups = 'drop') |>
  filter(Count > 1)

overlap_samples
###  A tibble: 4 × 2
###  `Sample ID`          Count
###  <chr>                <int>
###  1 TARGET-30-PANKFE-01A     2
###  2 TARGET-30-PAPTFZ-01A     2
###  3 TARGET-30-PASUML-01A     2
###  4 TARGET-30-PASYPX-01A     2

### So we choose abandon them.

TARGET_pretreatment <- TARGET_sheet |> 
  select(c(`Case ID`, `Sample ID`, `Sample Type`)) |>
  rename("Patient_ID" = `Case ID`,
         "Sample_ID" = `Sample ID`,
         "Sample_type"  =  `Sample Type`) |> 
  filter(!grepl(",", Sample_ID)) |> 
  filter(grepl("^TARGET-(20|00|21|30)", Sample_ID)) |> 
  unique()

### TARGET_xena_info pretreatment
### This data use
TARGET_xena_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_xena_info.tsv", sep = "\t") 
TARGET_xena_pretreatment   <- TARGET_xena_info |>
                              mutate(Patient_ID = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+")) |>
                              filter(grepl("^TARGET-(20|00|21|30)", Sample_ID)) |>
                              mutate(Age = round(Age, digits = 0),
                                     across(c(Age, OS_Time, OS_Status), ~ as.integer(.))) |>
                              rename(Cancer = disease_code)

### TARGET_info pretreatment
### In the TARGET database, rows with the info column marked as '--' are annotated as NA.
filter_lose <- function(x) {
  is.na(x) | x == "" | x == "--" | x == "Unknown" | x == "'--"
 }
TARGET_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_info.tsv", sep = "\t") 
replacement_table2 <- c(
  "Acute lymphocytic leukemia" =  "TARGET-ALL",             
  "'--"   = "TARGET-Unkown",                                  
  "Precursor B-cell lymphoblastic leukemia" = "TARGET-B-ALL", 
  "Not Reported"  = "TARGET-Unkown",                         
  "Acute myeloid leukemia, NOS" = "TARGET-LAML",             
  "Clear cell sarcoma of kidney" = "TARGET-CCSK",           
  "Neuroblastoma, NOS"  = "TARGET-NBL",                     
  "Ganglioneuroblastoma" = "TARGET-GNB",                   
  "Osteosarcoma, NOS"   = "TARGET-OS",                    
  "Malignant rhabdoid tumor" = "TARGET-MRT",               
  "Wilms tumor" = "TARGET-WT"  
)
matched2 <- TARGET_info$primary_diagnosis %in% names(replacement_table2)
TARGET_info_pretreatment <- TARGET_info |>
                            mutate(cohort = case_when(
                                   matched2 ~ replacement_table2[primary_diagnosis],
                                   TRUE ~ as.character(NA)  
                                   )) |>
                            rename(Cancer = primary_diagnosis,
                                   Project_id = project_id) |> 
                            filter(grepl("^TARGET-(20|00|21|30)", Sample_ID)) |> 
                            mutate(across(c("Age", "Sex", "OS_Time", "OS_Status", "race"), ~ ifelse(filter_lose(.x), NA, .x))) |> 
                            mutate(Cancer = ifelse(Cancer == "'--", NA, Cancer))

### Check if the valid rows in TARGET_info match the valid rows in TARGET_xena_info.

### 1.Filter the valid rows.
filter_valid <- function(x) {
  !is.na(x) & x != "" & x != "--" & x != "Unknown" & x == "'--"
}

### Filter the valid rows in TARGET_info.
TARGET_info_valid <- TARGET_info_pretreatment %>%
                     filter(filter_valid(Age) & filter_valid(OS_Time) & filter_valid(OS_Status) & filter_valid(Sex)) |>
                     select(c(Patient_ID, Age, OS_Time, OS_Status, Sex))

### Filter the valid rows in TARGET_xena_info.
TARGET_xena_info_valid <- TARGET_xena_pretreatment %>%
                          filter(filter_valid(Age) & filter_valid(OS_Time) & filter_valid(OS_Status) & filter_valid(Sex)) |>
                          select(c(Patient_ID, Age, OS_Time, OS_Status, Sex))

### 2. Merge the two data frames based on Patient_ID.
merged_data <- TARGET_info_valid %>%
               inner_join(TARGET_xena_info_valid, by = "Patient_ID", suffix = c("_info", "_xena"))

### 3. Check whether the valid values in the merged data frame are consistent.
is_consistent <- merged_data %>%
  rowwise() %>%
  mutate(
    Age_consistent = Age_info == Age_xena,
    Sex_consistent = Sex_info == Sex_xena,
    OS_Time_consistent = OS_Time_info == OS_Time_xena,
    OS_Status_consistent = OS_Status_info == OS_Status_xena
  ) %>%
  ungroup()

### 4. Age not onsistent
Age_not_consistent <- is_consistent |> 
                      filter(Age_consistent == "FALSE") 

### 5. Sex not onsistent
Sex_not_consistent <- is_consistent |> 
                      filter(Sex_consistent == "FALSE") 
### Sex_not_consistent

###  A tibble: 0 × 13

### 6. OS_Time not onsistent
OS_Time_not_consistent <- is_consistent |> 
                          filter(OS_Time_consistent == "FALSE") |>
                          mutate(
                          OS_Time_info_greater = OS_Time_info > OS_Time_xena,  
                          OS_Time_xena_greater = OS_Time_xena > OS_Time_info,  
                          OS_Time_equal = OS_Time_info == OS_Time_xena         
                          )
nrow(OS_Time_not_consistent |> filter(OS_Time_info_greater == "FALSE"))
### [1]0
nrow(OS_Time_not_consistent |> filter(OS_Time_xena_greater == "TRUE"))
### [1]0
nrow(OS_Time_not_consistent |> filter(OS_Time_equal == "TRUE"))
### [1]0

### 7. OS_Status not onsistent
OS_Status_not_consistent <- is_consistent |> 
                            filter(OS_Status_consistent == "FALSE") 

unique(OS_Status_not_consistent$OS_Status_info)
### [1] 1
unique(OS_Status_not_consistent$OS_Status_xena)
### [1] 0

### All evidence shows that the GDC data has been updated more recently than the Xena data.
### Therefore, if there is a conflict between GDC data and Xena data, we should use the GDC data.

### Use the valid values from TARGET_xena_info to fill the missing values in TARGET_info.
### 1: Determine the matching Patient_ID.
matched_ids <- intersect(TARGET_xena_pretreatment$Patient_ID, TARGET_info_pretreatment$Patient_ID)
setDT(TARGET_xena_pretreatment)
setDT(TARGET_info_pretreatment)

### 2: Find the columns that require updates.
cols_to_update <- setdiff(names(TARGET_info_pretreatment), c("Sample_ID", "Patient_ID", "Cancer","project_id","cohort","race"))

### 3.Iterate through each matching Patient_ID and update TARGET_info_pretreatment using TARGET_xena_pretreatment.
### update Age
TARGET_xena_can  <-  TARGET_xena_pretreatment |> 
                     filter(if_any("Age", ~ !filter_lose(.x))) 

TARGET_info_match <- TARGET_info_pretreatment |>
                            mutate(Age = case_when(
                                                   filter_lose(Age) & Patient_ID %in% TARGET_xena_can$Patient_ID ~  TARGET_xena_can$Age[match(Patient_ID, TARGET_xena_can$Patient_ID)],
                                                   TRUE ~ Age))

### update Sex
TARGET_xena_can  <-  TARGET_xena_pretreatment |> 
                     filter(if_any("Sex", ~ !filter_lose(.x))) 

TARGET_info_match1 <- TARGET_info_match |>
                            mutate(Sex = case_when(
                                                   filter_lose(Sex) & Patient_ID %in% TARGET_xena_can$Patient_ID ~  TARGET_xena_can$Sex[match(Patient_ID, TARGET_xena_can$Patient_ID)],
                                                   TRUE ~ Sex))

### update OS_Time
TARGET_xena_can  <-  TARGET_xena_pretreatment |> 
                     filter(if_any("OS_Time", ~ !filter_lose(.x))) 

TARGET_info_match2 <- TARGET_info_match1 |>
                      mutate(OS_Time = case_when(
                                                   filter_lose(OS_Time) & Patient_ID %in% TARGET_xena_can$Patient_ID  & OS_Status %in% TARGET_xena_can$OS_Status ~  TARGET_xena_can$OS_Time[match(Patient_ID, TARGET_xena_can$Patient_ID)],
                                                   TRUE ~ OS_Time))

### update OS_Status
TARGET_xena_can  <-  TARGET_xena_pretreatment |> 
                     filter(if_any("OS_Status", ~ !filter_lose(.x))) 

TARGET_info_after_match <- TARGET_info_match2 |>
                            mutate(OS_Status = case_when(
                                                   filter_lose(OS_Status) & Patient_ID %in% TARGET_xena_can$Patient_ID  & OS_Time %in% TARGET_xena_can$OS_Time ~  TARGET_xena_can$OS_Status[match(Patient_ID, TARGET_xena_can$Patient_ID)],
                                                   TRUE ~ OS_Status)) 
  
### Use TARGET_sheet to perform a left join, ensuring that the Sample_ID is correct.
TARGET_info_after_match <- TARGET_info_after_match |> select(- Sample_ID)
TARGET_correct      <- TARGET_pretreatment |>
                       left_join(TARGET_info_after_match, by = c("Patient_ID" = "Patient_ID")) |>
                       mutate(Run = Sample_ID)

### Replace the cohort column values with Project_ID where the row is TARGET-Unkown.
### TARGET_correct$project_id |> unique()
### [1] "TARGET-AML" "TARGET-NBL"
TARGET_correct_after_pro <- TARGET_correct  %>%
  mutate(
    cohort = case_when(
             cohort == "TARGET-Unkown" ~ Project_id,
             TRUE ~ cohort),
    cohort = case_when(
             cohort == "TARGET-AML" ~ "TARGET-LAML",
             TRUE ~ cohort)
  ) 

### Divide a part of data into a group based on whether the Sample_type is cell

target_cell <- TARGET_correct_after_pro |> filter(`Sample_type` == "Cell Lines")
target_cell_sample <- target_cell$Sample_ID |> unique()
target_cell_patient <- target_cell$Patient_ID |> unique()
### target_cell_sample
### [1] "TARGET-20-MUTZ3-50A"         "TARGET-20-ML1-50A"           "TARGET-20-MV411D1-50A"       "TARGET-20-KasumiD1-50A"    "TARGET-20-HEK293-50A"      "TARGET-20-MOLM14CBFGLIS-50A"     
### [7] "TARGET-20-TF1-50A"           "TARGET-20-IGROV1-50A"        "TARGET-20-CMS-50A"           "TARGET-20-WSUAML-50A"      "TARGET-20-KasumiAZAD5-50A" "TARGET-20-KasumiAZAD11-50A"
### [13]"TARGET-20-REH-50A"           "TARGET-20-MV411AZAD11-50A"   "TARGET-20-MV411AZAD5-50A"    "TARGET-20-THP1-50A"        "TARGET-20-OCIAML2-50A"     "TARGET-20-HL60-50A"  
### target_cell_patient
### [1] "TARGET-20-MUTZ3"         "TARGET-20-ML1"           "TARGET-20-MV411D1"       "TARGET-20-KasumiD1"      "TARGET-20-HEK293"        "TARGET-20-MOLM14CBFGLIS"
### [7] "TARGET-20-TF1"           "TARGET-20-IGROV1"        "TARGET-20-CMS"           "TARGET-20-WSUAML"        "TARGET-20-KasumiAZAD5"   "TARGET-20-KasumiAZAD11" 
### [13]"TARGET-20-REH"           "TARGET-20-MV411AZAD11"   "TARGET-20-MV411AZAD5"    "TARGET-20-THP1"          "TARGET-20-OCIAML2"       "TARGET-20-HL60"  
### reference : https://www.cellosaurus.org/
### MUTZ3：Acute myeloid leukemia
### ML1：	Adult acute myeloid leukemia
### MV411D1：biphenotypic B-myelomonocytic leukemia
### KasumiD1：Acute myeloblastic leukemia
### TF1：Acute erythroid leukemia 
### IGROV1：	Ovarian endometrioid adenocarcinoma
### WSUAML：Acute myeloid leukemia。
### KasumiAZAD5：can't find, so KasumiD1：Acute myeloblastic leukemia
### KasumiAZAD11：can't find, so KasumiD1：Acute myeloblastic leukemia
### REH：Childhood B acute lymphoblastic leukemia
### MV411AZAD11：can't find, so MV411D1：biphenotypic B-myelomonocytic leukemia
### MV411AZAD5：can't find, so MV411D1：biphenotypic B-myelomonocytic leukemia
### THP1：Childhood acute monocytic leukemia
### OCIAML2：Adult acute myeloid leukemia
### HL60：Adult acute myeloid leukemia
### HEK293：Vaccine production cell line.
### CML：Childhood acute megakaryoblastic leukemia
### MOLM14CBFGLIS：can't find, so MOLM14：Adult acute myeloid leukemia
cell_data <- metafusion[metafusion$Sample_ID %in% target_cell_sample]
### > nrow(cell_data)
###  [1] 1853
cell_info <- TARGET_info_pretreatment[TARGET_info_pretreatment$Patient_ID %in% target_cell_patient]
###  nrow(cell_info)
###  [1] 18
### no os.etc
cell_info2 <- TARGET_xena_pretreatment[TARGET_xena_pretreatment$Patient_ID %in% target_cell_patient]
### nrow(cell_info2)
### [1] 0

TARGET_info_after_cell <- TARGET_correct_after_pro |>
  mutate(cohort = case_when(
    `Sample_type` %in% "Cell Lines" ~ "TARGET-CELL",
    TRUE ~ cohort  
  ))

### The criteria for distinguishing normal samples based on Sample_type column in the TARGET databases are as follows:
TARGET_info_after_cell$Sample_type |> unique()
### [1] "Primary Blood Derived Cancer - Bone Marrow"              "Recurrent Blood Derived Cancer - Bone Marrow"           
### [3] "Blood Derived Normal"                                    "Primary Blood Derived Cancer - Peripheral Blood"        
### [5] "Bone Marrow Normal"                                      "Cell Lines"                                             
### [7] "Next Generation Cancer Model"                            "Blood Derived Cancer - Bone Marrow, Post-treatment"     
### [9] "Primary Tumor"                                           "Recurrent Tumor"                                        
### [11] "Recurrent Blood Derived Cancer - Peripheral Blood"       "Blood Derived Cancer - Peripheral Blood, Post-treatment"
### 
### If the sample type is "Blood Derived Normal" or "Bone Marrow Normal", then it is a normal sample. Otherwise, it is a tumor sample.
TARGET_info_after_normal <- TARGET_info_after_cell |>
  mutate(cohort = case_when(
    `Sample_type` %in% c("Blood Derived Normal", "Bone Marrow Normal") ~ "TARGET-NORMAL",
    TRUE ~ cohort  
  ))

### Adjust the column order.
new_order <- c("Patient_ID", "Sample_ID", "Run", "Sample_type", "Age", "Sex", "race", "Cancer", 
               "OS_Time", "OS_Status", "Project_id", "cohort")

TARGET_info_after_processing <- TARGET_info_after_normal |>
                                select(all_of(new_order)) 

### make sure nrow is same between Before and After Processing. 
TARGET_pretreatment |> nrow()
### [1]  3222
TARGET_info_after_processing |> nrow() 
### [1] 3222

### make sure IOBR/TME data just contain sample which target_sheet gave. 
detected_id_TARGET <- TARGET_info_after_processing$Run

### save info
cohorts <- unique(TARGET_info_after_processing$cohort)

for (cohort1 in cohorts) {
  cohort_data <- subset(TARGET_info_after_processing, cohort == cohort1)
  name <- gsub("-", "_", cohort1)
  file_name <- paste0("/home/data2/Projects/Fusion/fusiondb_v3/Clininfo/", name ,"_info.tsv")
  fwrite(cohort_data , file = file_name, sep = "\t")  
}

##Gene
metafusion1 <- metafusion |> 
  mutate(
    cohort = case_when(
      Sample_ID %in% TARGET_info_after_processing$Run & grepl("TARGET", cohort) ~ TARGET_info_after_processing$cohort[match(Sample_ID, TARGET_info_after_processing$Run)],
      TRUE ~ cohort
    )
  )  
  cp_file_name1 <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_calls_matched.cptac.target.tsv"
  fwrite(metafusion1, file = cp_file_name1, sep = "\t") 

##TME
NBL_TME <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-NBL_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-NBL")
AML_TME <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-AML_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-AML")
TARGET_TME <- rbind(NBL_TME,AML_TME) |>
              mutate(dash_count = str_count(ID, "-"))
### Manually check and delete composite samples from the info dataset.(TARGET-30-PAPTFZ-01A-01R-TARGET-30-PAPEAV-01A-01R)
TARGET_TME |> 
filter(dash_count == 9) |> 
select(ID) |>
unique()
### 1: TARGET-30-PASYPX-01A-01R-TARGET-30-PAIXIF-01A-01R
### 2: TARGET-30-PANKFE-01A-01R-TARGET-30-PAPUAR-01A-01R
### 3: TARGET-30-PANUKV-01A-01R-TARGET-30-PASUML-01A-01R
### 4: TARGET-30-PAPTFZ-01A-01R-TARGET-30-PAPEAV-01A-01R

TARGET_TME_pretreatment <- TARGET_TME |>
                           filter(dash_count != 9) |>
                           select(-dash_count)

### Filter TARGET Replicate Samples
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
TARGET_TME_after_replicate <- TARGET_TME_pretreatment |>
                              mutate(
                              prefixes = str_remove(ID, "-[^-]+$"),  
                              suffixes = str_extract(ID, "[^-]+$"), 
                              second_chars = str_extract(suffixes, "[0-9]+"),
                              second_chars = as.integer(second_chars)) |>
                              group_by(prefixes) |>
                              slice(which.max(second_chars)) |>
                              ungroup() |>
                              mutate(ID = prefixes) |> 
                              filter(ID %in% detected_id_TARGET)  |> ### Ensure that all IDs in the TME data can be found in the TARGET info data. 
                              select(-c(prefixes,suffixes,second_chars))

### Using TARGET info data cohort column to annotate TARGET TME data 
cohort_data   <- TARGET_info_after_processing |> 
                 select(Run, cohort)   |>
                 unique()   

TARGET_TME_after_cohort <- TARGET_TME_after_replicate |>
                           mutate(Run = ID)  |>
                           left_join(cohort_data, by = "Run")  |>
                           select(-Run)

### save TME
cohorts_TME <- unique(TARGET_TME_after_cohort$cohort)

for (cohort1 in cohorts_TME) {
  cohort_data <- subset(TARGET_TME_after_cohort, cohort == cohort1) |> select(-cohort)
  file_name <- paste0("/home/data2/Projects/Fusion/fusiondb_v3/IOBR/", cohort1 ,"_TME.tsv.gz")
  write.table(cohort_data, 
              file = gzfile(file_name), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE) 
}  

## IOBR
NBL_IOBR_signature <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-NBL_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-NBL")
AML_IOBR_signature <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-AML_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-AML")
TARGET_IOBR <- rbind(NBL_IOBR_signature,AML_IOBR_signature) |>
               mutate(dash_count = str_count(ID, "-"))

### Manually check and delete composite samples from the info dataset.(TARGET-30-PAPTFZ-01A-01R-TARGET-30-PAPEAV-01A-01R)
TARGET_IOBR |> 
filter(dash_count == 9) |> 
select(ID) |>
unique()
### 1: TARGET-30-PASYPX-01A-01R-TARGET-30-PAIXIF-01A-01R
### 2: TARGET-30-PANKFE-01A-01R-TARGET-30-PAPUAR-01A-01R
### 3: TARGET-30-PANUKV-01A-01R-TARGET-30-PASUML-01A-01R
### 4: TARGET-30-PAPTFZ-01A-01R-TARGET-30-PAPEAV-01A-01R

TARGET_IOBR_pretreatment <- TARGET_IOBR |>
                            filter(dash_count != 9) |>
                            select(-dash_count)

TARGET_IOBR_after_replicate <- TARGET_IOBR_pretreatment |>
                            mutate(
                            prefixes = str_remove(ID, "-[^-]+$"),  
                            suffixes = str_extract(ID, "[^-]+$"), 
                            second_chars = str_extract(suffixes, "[0-9]+"),
                            second_chars = as.integer(second_chars)) |>
                            group_by(prefixes) |>
                            slice(which.max(second_chars)) |>
                            ungroup() |>
                            mutate(ID = prefixes) |> 
                            filter(ID %in% detected_id_TARGET)  |> ### Ensure that all IDs in the TME data can be found in the TARGET info data. 
                            select(-c(prefixes,suffixes,second_chars))

### Using TARGET info data cohort column to annotate TARGET TME data 
cohort_data   <- TARGET_info_after_processing |> 
                 select(Run, cohort)   |>
                 unique()   

TARGET_IOBR_after_cohort <- TARGET_IOBR_after_replicate |>
                            mutate(Run = ID)  |>
                            left_join(cohort_data, by = "Run")  |>
                            select(-Run)
### save IOBR
cohorts_IOBR <- unique(TARGET_IOBR_after_cohort$cohort)

for (cohort1 in cohorts_IOBR) {
  cohort_data <- subset(TARGET_IOBR_after_cohort, cohort == cohort1) |> select(-cohort)
  file_name <- paste0("/home/data2/Projects/Fusion/fusiondb_v3/IOBR/", cohort1 ,"_IOBR_signature.tsv.gz")
  write.table(cohort_data, 
              file = gzfile(file_name), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE) 
}  

