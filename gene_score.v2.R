library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

file_meta <- "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_gdc_divide.tsv"
metafusion <- fread(file_meta , sep = "\t")
cohort1 <- fread("/data/data2/Projects/Fusion/fusiondb_v4/cohorts_meta.csv")
# Stratification Strategy
## Tumor data
### 1. Assign a column num_cohorts based on the number of distinct cohorts corresponding to each disease, +1 score; 
### When num_cohorts > 1, there should be at least one disease that corresponds to at least two different cohorts, and each of these cohorts must have at least one sample_id that has detected the same fusion value, +1 score; 
### When num_cohorts == 1, there must be at least two different sample_ids within the same cohort that have detected the same fusion value, +1 score;
### 2. If one fusion event's max_split_cnt and max_span_cnt are both greater than or equal to 1, +1 score;
### 3. If Arriba and STAR-Fusion detected identical fusion event, +1 score; 
### 4. score >= 3 , high confident fusion, else low confident fusion
Tumor_data <-  metafusion  |>
               mutate(
                Tumor_score = 0,
                fusion = paste0(gene1,"::",gene2),
                disease = case_when(
                  grepl( "TCGA", cohort) ~ gsub("^TCGA-", "", cohort),
                  grepl( "TARGET", cohort) ~ gsub("^TARGET-", "", cohort),
                  grepl( "CPTAC", cohort) ~ gsub("^CPTAC-", "", cohort),
                  !cohort %in% c("TCGA", "TARGET", "CPTAC") ~ cohort1$cancer_type[match(cohort, cohort1$id)],
                  TRUE ~ NA
                )) 
cohort_count <- Tumor_data %>%
               group_by(disease) %>%
               summarise(num_cohorts = n_distinct(cohort)) 
             
### 1. frequence
Tumor_data1 <- Tumor_data %>%
  left_join(cohort_count, by = "disease") %>%
  group_by(disease, fusion) %>%
  mutate(n_cohorts_with_fusion = n_distinct(cohort[Sample_ID != "" & !is.na(Sample_ID)])) %>%
  group_by(disease, fusion, cohort) %>%
  mutate(Tumor_score = case_when(
    (num_cohorts > 1 & n_cohorts_with_fusion >= 2) |  
    (num_cohorts == 1 & n_distinct(Sample_ID) >= 2) ~ Tumor_score + 1,
   TRUE ~ Tumor_score)
  ) %>%
  ungroup() |>
  select(-c(n_cohorts_with_fusion, num_cohorts))

###  2. Reads counts

Tumor_data2 <- Tumor_data1 %>%
  mutate(
    Tumor_score = case_when(
    max_split_cnt >= 1 & max_span_cnt >= 1 ~ Tumor_score + 1,
      TRUE ~ Tumor_score
    )
  )

###  3. Tools detected

Tumor_data3 <- Tumor_data2 %>%
    mutate(
      Tumor_score = case_when(
           tools == "arriba,starfusion" ~ Tumor_score + 1,
           TRUE ~ Tumor_score)
          )



### New metafusion
metafusion_new <- Tumor_data3  |> 
  mutate(Fusion_confidence = case_when(
    Tumor_score == 3 ~ "High Confidence",
    TRUE ~ "Low Confidence")) |> 
               select(-c(disease,fusion, Tumor_score)) 

file_name <- "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_calls.star.tsv"
write.table(metafusion_new, 
                      file = file_name, 
                      sep = "\t", 
                      row.names = FALSE, 
                      quote = FALSE) 

  
  
