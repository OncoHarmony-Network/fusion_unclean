library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

file_meta <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_gdc_divide.tsv"
metafusion <- fread(file_meta , sep = "\t")
# Stratification Strategy
## Normal data
### 1. If one fusion event was found in the past Arriba blocklist, it will be given a score of 0 and will not be changed.
### 2. If one fusion event's max_split_cnt >= 1 and max_span_cnt >= 1, +1 score.
### 3. If Arriba and STAR-Fusion detected identical fusion event, +1 score.
### 4. If past FusionAnnotator annotated one normal fusion event, +2 score(Because it is a strong evidence)
Normal_data <- metafusion |>
               filter(sample_type == "Normal") |>
               mutate(Normal_score = 0) 

### 1. Arriba blocklist
#### reference:https://arriba.readthedocs.io/en/latest/input-files/#blacklist
Normal_data <- Normal_data %>%
  mutate(
    Normal_score = case_when(
      blocklist == TRUE ~ 0,  
      TRUE ~ Normal_score     
    ),
    score0 = case_when(
      blocklist == TRUE ~ TRUE,  
      TRUE ~ FALSE    
    )
  )

###  2. Reads counts
Normal_data <- Normal_data %>%
  mutate(
    Normal_score = case_when(
         max_split_cnt >= 1 & max_span_cnt >= 1  & score0 == FALSE ~ Normal_score + 1 ,
         TRUE ~ Normal_score 
                            )
        )

###  3. Tools detected
Normal_data <- Normal_data %>%
    mutate(
    Normal_score = case_when(
           tools == "arriba,starfusion" & score0 == FALSE ~ Normal_score + 1,
           TRUE ~ Normal_score
                            )
          )
            
### 4. Fusionanntator annotate 
normal_markers <- c("GTEx_recurrent_StarF2019", "DGD_PARALOGS", "HGNC_GENEFAM", "BodyMap", "Greger_Normal", "ConjoinG", "Babiceanu_Normal")

#### Past Red Herrings summary didn't contain TCGA, TARGET, CPTAC(except for GTEx-Brain) normal sample data.
#### reference:https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki
#### BodyMap	40
#### Babiceanu_Normal	239
#### ConjoinG	698
#### DGD_PARALOGS	14463
#### Greger_Normal	113
#### GTEx_recurrent_StarF2019	4047
#### HGNC_GENEFAM	Unknown

Normal_data <- Normal_data %>%
  mutate(
    Normal_score = case_when(
      rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) >= 1 & score0 == FALSE ~ Normal_score + 2,
      TRUE ~ Normal_score   )
        )

### Nomral true 
Normal_data_true <- Normal_data |>
                          group_by(gene1, gene2) %>%
                          slice_max(Normal_score) %>%
                          ungroup() %>%
                          filter(Normal_score >= 2) 

## Tumor data
### 1. If one fusion event was found in the past Arriba blocklist, it will be given a score of 0 and will not be changed.
### 2. If fusion event was in our detcted normal fusion（score>=2), it will be given a score of 0 and will not be changed.
### 3. If one fusion event's max_split_cnt and max_span_cnt are both greater than or equal to 1, and at least one of them is less than 5, +1 score.
###    If one fusion event's max_split_cnt and max_span_cnt are both greater than or equal to 5, +2 score.
### 4. If Arriba detected exclusive fusion event, +0 score; 
###    If STAR-Fusion detected exclusive fusion event, +0 score; 
###    If Arriba and STAR-Fusion detected identical fusion event, +1 score; 
### 5. If past FusionAnnotator annotated one tumor fusion event, +2 score(Because it is a strong evidence)
Tumor_data <-  metafusion |>
               filter(sample_type == "Tumor" | sample_type == "Blood") |>
               mutate(Tumor_score = 0) 

### 1. Arriba blocklist
#### reference:https://arriba.readthedocs.io/en/latest/input-files/#blacklist
Tumor_data <- Tumor_data %>%
  mutate(
    Tumor_score = case_when(
      blocklist == TRUE ~ 0,  
      TRUE ~ Tumor_score   
    ),
    score0 = case_when(
      blocklist == TRUE ~ TRUE,  
      TRUE ~ FALSE    
    )
  )

### 2. Filter not in normal data.
Tumor_data <- Tumor_data %>%
              mutate(gene_pair = paste0(gene1, "::", gene2))

Normal_data_true <- Normal_data_true %>%
                    mutate(gene_pair = paste0(gene1, "::", gene2))

Tumor_data <- Tumor_data %>%
  mutate(
    Tumor_score = case_when(
      gene_pair %in% Normal_data_true$gene_pair & score0 == FALSE ~ 0,
      TRUE ~ Tumor_score   
    ),
    score0 = case_when(
      gene_pair %in% Normal_data_true$gene_pair & score0 == FALSE ~ TRUE,
      TRUE ~ FALSE    
    )
  )

###  3. Reads counts

Tumor_data <- Tumor_data %>%
  mutate(
    Tumor_score = case_when(
      max_split_cnt >= 5 & max_span_cnt >= 5  & score0 == FALSE ~ Tumor_score + 2,
      (max_split_cnt >= 1 & max_span_cnt >= 1) & ((max_split_cnt < 5 & max_span_cnt >= 5)| (max_split_cnt >= 5 & max_span_cnt < 5)| (max_split_cnt < 5 & max_span_cnt < 5))~ Tumor_score + 1,
      TRUE ~ Tumor_score
    )
  )

###  4. Tools detected

Tumor_data <- Tumor_data %>%
    mutate(
      Tumor_score = case_when(
           tools == "arriba,starfusion" & score0 == FALSE ~ Tumor_score + 1,
           tools == "arriba" & score0 == FALSE ~ Tumor_score,
           tools == "starfusion" & score0 == FALSE ~ Tumor_score,
           TRUE ~ Tumor_score
                            )
          )

### 5. Fusionanntator annotate 
tumor_markers <- c("Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic", "YOSHIHARA_TCGA", "ChimerKB", "ChimerPub", "ChimerSeq", "Klijn_CellLines", "Larsson_TCGA", "CCLE_StarF2019", "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018", "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023")

Tumor_data <- Tumor_data %>%
    mutate(
      Tumor_score = case_when(
        rowSums(sapply(tumor_markers, grepl, x = cancer_db_hits)) >= 1 ~ Tumor_score + 2,
        TRUE ~ Tumor_score)
          ) 

### New metafusion
Tumor_data1 <- Tumor_data  |> 
               select(-c(gene_pair, score0)) |>
               rename(gene_score = Tumor_score)
Normal_data1 <- Normal_data |> 
               select(-c( score0)) |>
               rename(gene_score = Normal_score)

metafusion_new <- rbind(Normal_data1,Tumor_data1)
file_name <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_calls.star.tsv"
write.table(metafusion_new, 
                      file = file_name, 
                      sep = "\t", 
                      row.names = FALSE, 
                      quote = FALSE) 

  
  
