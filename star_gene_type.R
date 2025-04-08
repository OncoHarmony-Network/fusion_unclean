library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggVennDiagram)
library(sf)
library(stringi)
library(stringr)
library(plotly)
library(reticulate)
library(RColorBrewer)
library(plotly)
library(paletteer)
library(ComplexHeatmap)


file_meta <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_gdc_divide.tsv"
metafusion <- fread(file_meta , sep = "\t")

# Stratification Strategy
##  
##  Level num
##  
##  GDC samples num 
### TCGA_gdc have 10545 samples, normal cohort have 729 samples.
### TARGET_gdc have 3222 samples, normal cohort have 425 samples.
### CPTAC_gdc have 2558 samples, normal cohort have 653 samples.
### GDC have 16325 sample, normal cohort have 1807 samples.

## NORMAL
### 1. If one fusion event was found in the past Arriba blocklist, it will be given a score of 0 and will not be changed.
### 2. If one fusion event's max_split_cnt >= 1 and max_span_cnt >= 1, +1 score.
### 3. If Arriba and STAR-Fusion detected identical fusion event, +1 score.
### 4. If past FusionAnnotator annotated one normal fusion event, +2 score(Because it is a strong evidence)
Normal_data <- metafusion |>
               filter(sample_type == "Normal") |>
               mutate(Normal_score = 0) 

### 1. metafusion$blocklist
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

#### Normal_score Tools detected code
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
#### DGD_PARALOGS	3543
#### Greger_Normal	113
#### GTEx_recurrent_StarF2019	4047
#### HGNC_GENEFAM	3543

Normal_data <- Normal_data %>%
  mutate(
    Normal_score = case_when(
      rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) >= 1 & score0 == FALSE ~ Normal_score + 2,
      TRUE ~ Normal_score   )
        )

### Because Metafusion merge rows.
Normal_data1 <- Normal_data %>%
               separate_rows(gene1, sep = ",") %>%
               separate_rows(gene2, sep = ",")

Normal_high_level <- Normal_data1 |>
                      filter(Normal_score >= 2) 

file_name <- "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Normal_high_level.tsv"
write.table(Normal_high_level , 
            file = file_name, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 


### Just choose the Highest gene pair for fusion event in different situation
Normal_score_max <- Normal_data1 %>%
                    group_by(gene1, gene2) %>%
                    summarise(
                            Max_Normal_score = max(Normal_score),
                            .groups = 'drop'
                          )
### Nomral proportion                                         
Normal_score_proportion <- Normal_score_max %>%
                          group_by(Max_Normal_score) %>%
                          summarise(
                            Gene_Pairs = n(),
                            .groups = 'drop'
                          ) %>%
                          mutate(Proportion = Gene_Pairs / sum(Gene_Pairs))


### Tumor data
#### 1. If one fusion event was found in the past Arriba blocklist, it will be given a score of 0 and will not be changed.
#### 2. If fusion event was in our detcted normal fusion（score>=2), it will be given a score of 0 and will not be changed.
#### 3. If one fusion event's max_split_cnt and max_span_cnt are both greater than or equal to 1, and at least one of them is less than 5, +1 score.
####    If one fusion event's max_split_cnt and max_span_cnt are both greater than or equal to 5, +2 score.
#### 4. If Arriba detected exclusive fusion event, +0 score; 
####    If STAR-Fusion detected exclusive fusion event, +0 score; 
####    If Arriba and STAR-Fusion detected identical fusion event, +1 score; 
#### 5. If past FusionAnnotator annotated one tumor fusion event, +2 score(Because it is a strong evidence)
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

Normal_data_true <- Normal_high_level |>
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

### Just choose the Highest gene pair for fusion event in different situation
Tumor_data_max <- Tumor_data %>%
                  group_by(gene1, gene2) %>%
                  summarise(
                            Max_Tumor_score = max(Tumor_score),
                            .groups = 'drop'
                          )
### Tumor proportion                                         
Tumor_score_proportion <- Tumor_data_max %>%
                           group_by(Max_Tumor_score) %>%
                           summarise(
                            Gene_Pairs = n(),
                            .groups = 'drop'
                          ) %>%
                          mutate(Proportion = Gene_Pairs / sum(Gene_Pairs))
### Because Metafusion merge rows.
Tumor_data1 <-  Tumor_data %>%
                separate_rows(gene1, sep = ",") %>%
                separate_rows(gene2, sep = ",")
### Tumor high level 
Tumor_high_level      <- Tumor_data1 |>
                         filter(Tumor_score>= 2) 

file_name <- "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Tumor_high_level.tsv"
write.table(Tumor_high_level , 
            file = file_name, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 

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



### Overlap part1 
### Figure1 原有的改为https://mp.weixin.qq.com/s/gaP9EtBlmAbq-Hj-8AFJNw的Abstract

### Figure2a
reticulate::use_python("/home/data2/home/zyk/.virtualenvs/myenv/bin/python")
reticulate::use_virtualenv("myenv")
tumor_high_level <- fread("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Tumor_high_level.tsv")
normal_high_level <- fread("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Normal_high_level.tsv")
cohort <-  fread("/home/data2/Projects/Fusion/fusiondb_v3/cohorts_meta.csv") 
directory <- "/home/data2/Projects/Fusion/fusiondb_v3/Clininfo"
files <- list.files(directory, pattern = "_info.tsv$", full.names = TRUE)
exclude_files <- c("TCGA_info.tsv", "TARGET_xena_info.tsv", "TARGET_info.tsv", "CPTAC_info.tsv")
files <- files[!basename(files) %in% exclude_files]
all_data_Run <- NULL
for (file in files) {
  data <- fread(file)
  if ("Run" %in% colnames(data)) {
    data <- data %>% select(Run)
    cohort_id <- tools::file_path_sans_ext(basename(file))
    cohort_id <- sub("_info$", "", cohort_id)
    data$cohort_id <- cohort_id
    all_data_Run <- bind_rows(all_data_Run, data)
  } else {
    warning(paste("File", basename(file), "does not contain a 'Run' column."))
  }
}
cohort1 <- cohort |> 
           select(id,cancer_type) |>
           rename(cancer_type1 = cancer_type)
all_data_Run_pro <- all_data_Run |>    
  mutate(
    cancer_type = ifelse(
      grepl("TCGA|CPTAC|TARGET", cohort_id),
      sub("^(TCGA|CPTAC|TARGET)_", "", cohort_id),
      NA 
    )
  )|>
    left_join(cohort1, by = c("cohort_id" = "id")) |>
    mutate(
      cancer_type = ifelse(is.na(cancer_type), cancer_type1, cancer_type)
    ) |>
      select(-cancer_type1) |>
        mutate(
          cohort_id = str_replace_all(cohort_id, "TCGA_", "TCGA-"),
          cohort_id = str_replace_all(cohort_id, "CPTAC_", "CPTAC-"),
          cohort_id = str_replace_all(cohort_id, "TARGET_", "TARGET-")
        )

### normal_high_level 
all_data_Run_pro_normal <- all_data_Run_pro |>
                           filter(cancer_type == "NORMAL")
normal_high_level_choose <-normal_high_level |> 
                           select(c(gene1,gene2,cohort,Sample_ID)) |>
                           pivot_longer(
                           cols = c(gene1, gene2),  
                           names_to = "gene_type",  
                           values_to = "gene_name" ) |> 
                           select(-gene_type) |>
                           distinct()

Run_normal <- all_data_Run_pro_normal |>
                            left_join(normal_high_level_choose, by = c("cohort_id" = "cohort", "Run" = "Sample_ID")) |>
                            mutate(
                              gene_name = ifelse(is.na(gene_name), "No high-level fusion events detected", gene_name)
                            ) 
### TCGA-NORMAL, TARGET-NORMAL, CPTAC-NORMAL
independent_normal_run_counts <- Run_normal %>%
                                 group_by(cohort_id) %>%
                                 summarise(total_runs = n_distinct(Run)) 

Run_normal <- Run_normal %>%
              left_join(independent_normal_run_counts, by = "cohort_id")                            
result_normal <-          Run_normal %>%
                          group_by(cohort_id, gene_name,total_runs) %>%
                          reframe(
                            detected_runs = n(),  
                            detection_rate = detected_runs / total_runs
                          ) |>
                          distinct() 
no_high_level_fusion_normal <- result_normal %>%
                          filter(gene_name == "No high-level fusion events detected")
top_3_per_cohort_normal <- result_normal %>%
               filter(gene_name != "No high-level fusion events detected") %>%
               group_by(cohort_id) %>%
               arrange(desc(detection_rate)) %>%
               slice_head(n = 3)

### GDC-NORMAL

all_normal_run_counts <- Run_normal %>%
                         summarise(total_runs = n_distinct(Run))   |>
                         pull()
all_result_normal <- Run_normal %>%
  mutate(total_runs = all_normal_run_counts) |>
  group_by(gene_name,total_runs) |>
  reframe(
    detected_runs = n(),  
    detection_rate = detected_runs / total_runs
  ) |>
  distinct() 
all_no_high_level_fusion_normal <- all_result_normal |>
                                   filter(gene_name == "No high-level fusion events detected") |>
                                   mutate(cohort_id = "GDC-NORMMAL")
all_top_3_per_cohort_normal <- all_result_normal |>
                               filter(gene_name != "No high-level fusion events detected") |>
                               arrange(desc(detection_rate)) |>
                               slice_head(n = 3) |>
                               mutate(cohort_id = "GDC-NORMMAL")

all_data_Run_circos_normal <- bind_rows(all_no_high_level_fusion_normal,all_top_3_per_cohort_normal,no_high_level_fusion_normal, top_3_per_cohort_normal)
saveRDS(all_data_Run_circos_normal , "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Run_circos_normal.rds")  
all_data_Run_circos_normal <- readRDS("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Run_circos_normal.rds")  |> 
                              arrange(cohort_id) 
cohort_id_normal <- all_data_Run_circos_normal$cohort_id |> unique() 
all_data_Run_circos_normal1 <- all_data_Run_circos_normal
all_data_Run_circos_normal1$detection_rate <- round(all_data_Run_circos_normal$detection_rate, digits = 5) 
data <- all_data_Run_circos_normal1
for(idn in cohort_id_normal){
  
  ids <- c("Root","nodetected","detected",
           "Top1_ndsdf","Top1_dsdf","Top1",
           "Top2_ndsdf","Top2_dsdf","Top2_top1","Top2",
           "Top3_ndsdf","Top3_dsdf","Top3_top1","Top3_top2","Top3")
  parents <- c("", "Root","Root", 
               "nodetected","detected","detected",
               "Top1_ndsdf","Top1_dsdf","Top1","Top1",
               "Top2_ndsdf","Top2_dsdf","Top2_top1","Top2","Top2")
  colors <- c("#ffffff", "#808080", "#32037d", 
              "#ffffff", "#ffffff", "#0001a1",
              "#ffffff", "#ffffff", "#ffffff", "#037f77",
              "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#e90f44")
  nodetected_rate <- data |>
    filter(cohort_id == idn) |> 
    filter(gene_name == "No high-level fusion events detected") |> 
    select(detection_rate) |> 
    pull()
  detected_rate <- 1 - nodetected_rate
  label_DSPF  <- paste0("Fusion+",
                        "<br>", round(detected_rate * 100, digits = 2),"%")
  label_NDSPF <- ""
  detectingtop <- data  |> 
    filter(cohort_id == idn) |> 
    filter(gene_name != "No high-level fusion events detected")  
  detecting_ids <-  detectingtop  |>
    arrange(desc(detection_rate)) |> 
    select(gene_name) |> 
    pull()
  detecting_rates <-  detectingtop  |>
    select(detection_rate) |> 
    arrange(desc(detection_rate)) |> 
    pull()
  label_top1 <- paste0(detecting_ids[1],
                       "<br>",  round(detecting_rates[1] * 100, digits = 2),"%")
  label_top2 <- paste0(detecting_ids[2],
                       "<br>", round(detecting_rates[2] * 100, digits = 2),"%")
  label_top3 <- paste0(detecting_ids[3],
                       "<br>", round(detecting_rates[3] * 100, digits = 2),"%")
  labels_all <- c("", label_NDSPF, label_DSPF, 
                  "", "", label_top1, 
                  "", "", "", label_top2, 
                  "", "", "", "",label_top3)
  nodetected_counts  <- data |>
    filter(cohort_id == idn) |> 
    filter(gene_name == "No high-level fusion events detected") |> 
    select(detected_runs) |> 
    pull()
  total_counts     <-   data |>
    filter(cohort_id == idn) |> 
    filter(gene_name == "No high-level fusion events detected") |> 
    select(total_runs) |> 
    pull()
  detected_counts  <-   total_counts - nodetected_counts
  detecting_counts <-   detectingtop  |>
    select(detected_runs) |> 
    arrange(desc(detected_runs)) |> 
    pull()
  counts_all <- c( total_counts, nodetected_counts, detected_counts,
                   nodetected_counts, detected_counts - detecting_counts[1], detecting_counts[1],                    
                   nodetected_counts, detected_counts - detecting_counts[1], detecting_counts[1] - detecting_counts[2], detecting_counts[2],
                   nodetected_counts, detected_counts - detecting_counts[1], detecting_counts[1] - detecting_counts[2], detecting_counts[2] - detecting_counts[3], detecting_counts[3]) 

  Figure2a_part <- data.frame(ids= ids, labels = labels_all, parents = parents, colors = colors, counts = counts_all) 


Figure2a_part1 <- Figure2a_part
Figure2a_part1$labels <- ""
Figure2a_part2 <-  Figure2a_part1[1:3,]
Figure2a_part3 <-  Figure2a_part1[1,]
Figure2a_part2[1,3] <- "circos"
Figure2a_part3[1,1] <- "circos"
Figure2a_part4<- rbind(Figure2a_part3,Figure2a_part2 )
p_figure2a_part1 <- plot_ly(Figure2a_part4, 
                           ids = ~ids, 
                           labels = ~labels, 
                           parents = ~parents, 
                           type = 'sunburst',
                           #textfont = list(color = "#ffffff",size = 166),
                           insidetextorientation = "radial",
                           values = ~counts, 
                           sort = FALSE,
                           branchvalues = 'total',
                           rotation = 90,
                           marker = list(colors = ~colors))
save_image(p_figure2a_part1, file = paste0("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/",idn,"_normal_Figure2a.pdf"), height = 1000, width = 1000)
save_image(p_figure2a_part1, file = paste0("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/",idn,"_normal_Figure2a.png"), height = 1000, width = 1000)
}
### Tumor 
all_data_Run_pro_tumor   <- all_data_Run_pro |>
                            distinct() |>
                            filter(cancer_type != "NORMAL")

tumor_high_level_choose  <- tumor_high_level |> 
                            select(c(gene1,gene2,cohort,Sample_ID)) |>
                            pivot_longer(
                              cols = c(gene1, gene2),  
                              names_to = "gene_type",  
                              values_to = "gene_name"  
                            ) |>
                            select(-gene_type) |>
                            distinct()

Run_tumor <- all_data_Run_pro_tumor |>
             left_join(tumor_high_level_choose, by = c("cohort_id" = "cohort", "Run" = "Sample_ID")) |>
             mutate(gene_name = ifelse(is.na(gene_name), "No high-level fusion events detected", gene_name)) 

independent_tumor_run_counts <- Run_tumor %>%
                                group_by(cancer_type) %>%
                                summarise(total_runs = n_distinct(Run))
Run_tumor <- Run_tumor %>%
             left_join(independent_tumor_run_counts, by = "cancer_type")   

result_tumor <- Run_tumor %>%
                group_by(cancer_type, gene_name,total_runs) %>%
                reframe(detected_runs = n(),  
                        detection_rate = detected_runs / total_runs) |>
                distinct() 
no_high_level_fusion_tumor <- result_tumor  %>%
                              filter(gene_name == "No high-level fusion events detected")
top_3_per_cohort_tumor <- result_tumor  %>%
                          filter(gene_name != "No high-level fusion events detected") %>%
                          group_by(cancer_type) %>%
                          arrange(desc(detection_rate)) %>%
                          slice_head(n = 3)
all_data_Run_circos_tumor  <- rbind(no_high_level_fusion_tumor, top_3_per_cohort_tumor)
saveRDS(all_data_Run_circos_tumor , "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Run_circos_tumor.rds")      
all_data_Run_circos_tumor <- readRDS("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Run_circos_tumor.rds") |> 
  arrange(cancer_type)
cancer_type_tumor <- all_data_Run_circos_tumor$cancer_type |> unique() 
all_data_Run_circos_tumor1 <- all_data_Run_circos_tumor
all_data_Run_circos_tumor1$detection_rate <- round(all_data_Run_circos_tumor$detection_rate, digits = 5) 
data <- all_data_Run_circos_tumor1
for(idn in cancer_type_tumor){
  data_nrow <- data |> 
    filter(cancer_type == idn) |>
    nrow()
  if(data_nrow == 3){
    ids <- c("Root","detected",
             "Top1_dsdf","Top1",
             "Top2_dsdf","Top2_top1","Top2",
             "Top3_dsdf","Top3_top1","Top3_top2","Top3")
    parents <- c("","Root", 
                 "detected","detected",
                 "Top1_dsdf","Top1","Top1",
                 "Top2_dsdf","Top2_top1","Top2","Top2")
    colors <- c("#ffffff", "#32037d", 
                "#ffffff", "#0001a1",
                "#ffffff", "#ffffff", "#037f77",
                "#ffffff", "#ffffff", "#ffffff", "#e90f44")
    detected_rate <- 1 
    label_DSPF <- paste0("Fusion+",
                         "<br>",round(detected_rate * 100, digits = 2),"%")
    detectingtop <- data  |> 
      filter(cancer_type == idn) 
    detecting_ids <-  detectingtop  |>
      arrange(desc(detection_rate)) |> 
      select(gene_name) |> 
      pull()
    detecting_rates <-  detectingtop  |>
      select(detection_rate) |> 
      arrange(desc(detection_rate)) |> 
      pull()
    label_top1 <- paste0(detecting_ids[1],
                         "<br>", round(detecting_rates[1] * 100, digits = 2),"%")
    label_top2 <- paste0(detecting_ids[2],
                         "<br>", round(detecting_rates[2] * 100, digits = 2),"%")
    label_top3 <- paste0( detecting_ids[3],
                          "<br>",  round(detecting_rates[3] * 100, digits = 2),"%")
    labels_all <- c("", label_DSPF, 
                    "", label_top1, 
                    "", "", label_top2, 
                    "", "", "",label_top3)
    total_counts     <-   data |>
      filter(cancer_type == idn) |> 
      filter(gene_name != "No high-level fusion events detected") |> 
      select(total_runs) |>  
      distinct() |>
      pull()
    detected_counts  <-   total_counts
    detecting_counts <-   detectingtop  |>
      select(detected_runs) |> 
      arrange(desc(detected_runs)) |> 
      pull()
    counts_all <- c( total_counts,  detected_counts,
                     detected_counts - detecting_counts[1], detecting_counts[1],                    
                     detected_counts - detecting_counts[1], detecting_counts[1] - detecting_counts[2], detecting_counts[2],
                     detected_counts - detecting_counts[1], detecting_counts[1] - detecting_counts[2], detecting_counts[2] - detecting_counts[3], detecting_counts[3]) 
    Figure2a_part <- data.frame(ids= ids, labels = labels_all, parents = parents, colors = colors, counts = counts_all) 
    Figure2a_part1 <- Figure2a_part
    Figure2a_part1$labels <- ""
    Figure2a_part2 <-  Figure2a_part1[1:2,]
    Figure2a_part3 <-  Figure2a_part1[1,]
    Figure2a_part2[1,3] <- "circos"
    Figure2a_part3[1,1] <- "circos"
    Figure2a_part4<- rbind(Figure2a_part3,Figure2a_part2 )
  }else if(data_nrow == 4){
    ids <- c("Root","nodetected","detected",
             "Top1_ndsdf","Top1_dsdf","Top1",
             "Top2_ndsdf","Top2_dsdf","Top2_top1","Top2",
             "Top3_ndsdf","Top3_dsdf","Top3_top1","Top3_top2","Top3")
    parents <- c("", "Root","Root", 
                 "nodetected","detected","detected",
                 "Top1_ndsdf","Top1_dsdf","Top1","Top1",
                 "Top2_ndsdf","Top2_dsdf","Top2_top1","Top2","Top2")
    colors <- c("#ffffff", "#808080", "#32037d", 
                "#ffffff", "#ffffff", "#0001a1",
                "#ffffff", "#ffffff", "#ffffff", "#037f77",
                "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#e90f44")
    nodetected_rate <- data |>
      filter(cancer_type == idn) |> 
      filter(gene_name == "No high-level fusion events detected") |> 
      select(detection_rate) |> 
      pull()
    detected_rate <- 1 - nodetected_rate
    label_DSPF  <- paste0("Fusion+",
                          "<br>", round(detected_rate * 100, digits = 2),"%")
    label_NDSPF <- paste0("")
    detectingtop <- data  |> 
      filter(cancer_type == idn) |> 
      filter(gene_name != "No high-level fusion events detected")  
    detecting_ids <-  detectingtop  |>
      arrange(desc(detection_rate)) |> 
      select(gene_name) |> 
      pull()
    detecting_rates <-  detectingtop  |>
      select(detection_rate) |> 
      arrange(desc(detection_rate)) |> 
      pull()
    label_top1 <- paste0(detecting_ids[1],
                         "<br>",  round(detecting_rates[1] * 100, digits = 2),"%")
    label_top2 <- paste0(detecting_ids[2],
                         "<br>", round(detecting_rates[2] * 100, digits = 2),"%")
    label_top3 <- paste0(detecting_ids[3],
                         "<br>", round(detecting_rates[3] * 100, digits = 2),"%")
    labels_all <- c("", label_NDSPF, label_DSPF, 
                    "", "", label_top1, 
                    "", "", "", label_top2, 
                    "", "", "", "", label_top3)
    nodetected_counts  <- data |>
      filter(cancer_type == idn) |> 
      filter(gene_name == "No high-level fusion events detected") |> 
      select(detected_runs) |> 
      pull()
    total_counts     <-   data |>
      filter(cancer_type == idn) |> 
      filter(gene_name == "No high-level fusion events detected") |> 
      select(total_runs) |> 
      pull()
    detected_counts  <-   total_counts - nodetected_counts
    detecting_counts <-   detectingtop  |>
      select(detected_runs) |> 
      arrange(desc(detected_runs)) |> 
      pull()
    counts_all <- c( total_counts, nodetected_counts, detected_counts,
                     nodetected_counts, detected_counts - detecting_counts[1], detecting_counts[1],                    
                     nodetected_counts, detected_counts - detecting_counts[1], detecting_counts[1] - detecting_counts[2], detecting_counts[2],
                     nodetected_counts, detected_counts - detecting_counts[1], detecting_counts[1] - detecting_counts[2], detecting_counts[2] - detecting_counts[3], detecting_counts[3]) 
    Figure2a_part <- data.frame(ids= ids, labels = labels_all, parents = parents, colors = colors, counts = counts_all) 
  }
  Figure2a_part1 <- Figure2a_part
  Figure2a_part1$labels <- ""
  Figure2a_part2 <-  Figure2a_part1[1:3,]
  Figure2a_part3 <-  Figure2a_part1[1,]
  Figure2a_part2[1,3] <- "circos"
  Figure2a_part3[1,1] <- "circos"
  Figure2a_part4<- rbind(Figure2a_part3,Figure2a_part2 )
  p_figure2a_part1 <- plot_ly(Figure2a_part4, 
                             ids = ~ids, 
                             labels = ~labels, 
                             parents = ~parents, 
                             type = 'sunburst',
                            # textfont = list(color = "#ffffff",size = 166),
                            # insidetextorientation = "radial",
                             values = ~counts, 
                             sort = FALSE,
                             branchvalues = 'total',
                             rotation = 90,
                             marker = list(colors = ~colors))
  save_image(p_figure2a_part1, file = paste0("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/",idn,"_tumor_Figure2a.pdf"), height = 1000, width = 1000)
  save_image(p_figure2a_part1, file = paste0("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/",idn,"_tumor_Figure2a.png"), height = 1000, width = 1000)
}

colnames(all_data_Run_circos_normal1)[5] <- "cancer_type"
tabelS1 <- rbind(all_data_Run_circos_normal1,all_data_Run_circos_tumor1)
tabelS1 <- tabelS1 |>
           mutate(gene_name = case_when(
                  gene_name == "No high-level fusion events detected" ~ "Fusion-",  
                  TRUE ~ gene_name))
colnames(tabelS1)[1] <- "Fusion_events"
colnames(tabelS1)[2] <- "Total_Sample_Num"
colnames(tabelS1)[3] <- "Detected_Sample_Num"
fwrite(tabelS1, "/home/data2/Projects/Fusion/paper_plot/Figure2/tabelS1.tsv")

### Figure2b
### no-GDC
reticulate::use_python("/home/data2/home/zyk/.virtualenvs/myenv/bin/python")
reticulate::use_virtualenv("myenv")
file_name <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_calls.star.tsv"
metafusion <- fread(file_name)
cohort <-  fread("/home/data2/Projects/Fusion/fusiondb_v3/cohorts_meta.csv") |> arrange(cancer_type, year)
num_years  <- length(unique(cohort$year))
cohort_gdc <- metafusion |>
              select(cohort) |>
              distinct() |>
              filter(grepl("TCGA|CPTAC|TARGET", cohort)) |>
              mutate(cancer_type = sub("^(TCGA|CPTAC|TARGET)-", "", cohort)) |> arrange(cancer_type)
cancer_type_gdc_nrow <- cohort_gdc$cancer_type |> unique() |> length()
### 39
color_cancer_type_all <- c(
  "COAD"  = "#BDC314",                       "UCEC" = "#C2900F", "LUAD" = "#77DCDD", "PDAC"  = "#848385", "LUSC"   = "#3981AF", "LAML"       = "#5F97C6", 
  "NBL"   = "#60966D", "GNB"  = "#E7308B", "ACC"  = "#B31761", "CESC" = "#338FEB", "CHOL"  = "#A5CC5B", "DLBC"   = "#088c00", "CELL"       = "#D2DC38", 
  "ESCA"  = "#5B3660", "KICH" = "#d40d8c", "KIRP" = "#00b3b0", "LGG"  = "#c2d0ea", "LIHC"  = "#c85e62", "MESO"   = "#C0845C", "OV"         = "#C3B0D2", 
  "PAAD"  = "#C99BFF", "PCPG" = "#A4E048", "PRAD" = "#037F77", "TGCT" = "#BABABA", "THCA"  = "#F4A99B", "THYM"   = "#BD9AAD", "UCS"        = "#E8D2B3",  
  "UVM"   = "#F6631C", "BLCA" = "#5861AC", "KIRC" = "#EF767B", "SKCM" = "#FFC080", "NSCLC" = "#AA66EB", "SCLC"   = "#39B086", "SARC"       = "#FF6100", 
  "BRCA"  = "#6AAF0A", "STAD" = "#6B7EB9", "HNSC" = "#AED6D5", "SGC"  = "#96A4C7", "GBM"   = "#965478", "NORMAL" = "#F1B656", "Pan-Cancer" = "#B58581", 
  "READ"  = "#592E13")

color_year <-  c("2012" = "#FEB44E", "2014" = "#FDA446", "2016" = "#FD943F", 
                 "2017" = "#FC7F38", "2018" = "#FA4B29", "2019" = "#E51F1D",
                 "2020" = "#D7111F", "2021" = "#C70623", "2022" = "#B30026", 
                 "2023" = "#800026")

color_id <- paletteer_c("scico::berlin", n = 100)
### id
years <- cohort$year
cancer_types_io <- cohort$cancer_type
names <- cohort$id
ids    <- c(as.character(cancer_types_io), as.character(names), as.character(years))
labels <- c(as.character(cancer_types_io), as.character(names),as.character(years))
Donut_Chart_data <- data.frame(ids = ids , labels = labels)

Donut_Chart_data$parents <- ""
Donut_Chart_data$parents[30:58] <- as.character(cancer_types_io)
Donut_Chart_data$parents[59:87] <- as.character(names)

Donut_Chart_data$colors <- ""

selected_ids_cancer <- Donut_Chart_data$ids[1:29]
Donut_Chart_data$colors[1:29]  <- color_cancer_type_all[selected_ids_cancer]

Donut_Chart_data$colors[30:58] <- color_id[1:29]

selected_ids_year <- Donut_Chart_data$ids[59:87]
Donut_Chart_data$colors[59:87] <- color_year[selected_ids_year]
Donut_Chart_data1 <- Donut_Chart_data[1:29,] |> distinct()
Donut_Chart_data <- bind_rows(Donut_Chart_data1, Donut_Chart_data[30:87, ])
Donut_Chart_data <- Donut_Chart_data |>
  mutate(
    id_count = row_number(),  
    ids_new = ifelse(
      duplicated(ids) | duplicated(ids, fromLast = TRUE),  
      paste0(ids, "_", ave(id_count, ids, FUN = seq_along)),  
      ids  
    ),
    ids = ids_new
  ) |>
  select(-c(id_count,ids_new) )
saveRDS(Donut_Chart_data, "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2b_data/Donut_Chart_data.rds")
p_figure2b <- plot_ly(Donut_Chart_data , 
  ids = ~ids, 
  labels = ~labels, 
  parents = ~parents, 
  type = 'sunburst', 
  sort = FALSE, 
  textfont = list(color = 'white',size = 180),
  marker = list(colors = ~colors))
p_figure2b
save_image(p_figure2b, file = "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2b.pdf", height = 1000, width = 1000)
save_image(p_figure2b, file = "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2b.png", height = 1000, width = 1000)
### Figure2c
### gdc
cancer_types_gdc <- cohort_gdc$cancer_type
names_gdc  <- cohort_gdc$cohort
ids_gdc    <- c(as.character(cancer_types_gdc), as.character(names_gdc))
labels_gdc <- c(as.character(cancer_types_gdc), as.character(names_gdc))
Donut_Chart_data_gdc <- data.frame(ids = ids_gdc , labels = labels_gdc )
Donut_Chart_data_gdc$parents <- ""
Donut_Chart_data_gdc$parents[51:100] <- Donut_Chart_data_gdc$ids[1:50] 
Donut_Chart_data_gdc$colors[51:100] <- color_id[42:91]
selected_ids_cancer_gdc <- Donut_Chart_data_gdc$labels[1:50]
Donut_Chart_data_gdc$colors[1:50]  <- color_cancer_type_all[selected_ids_cancer_gdc]
Donut_Chart_data_gdc1 <- Donut_Chart_data_gdc  |>
                         slice(1:50) |>
                         distinct(labels, .keep_all = TRUE)
Donut_Chart_data_gdc <- bind_rows(Donut_Chart_data_gdc1, Donut_Chart_data_gdc[51:100, ]) 
saveRDS(Donut_Chart_data_gdc,"/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2c_data/Donut_Chart_data_gdc.rds")
p_figure2c <- plot_ly(Donut_Chart_data_gdc , 
  ids = ~ids, 
  labels = ~labels, 
  parents = ~parents, 
  type = 'sunburst', 
  sort = FALSE, 
  textfont = list(color = 'white',size = 150),
  insidetextorientation = "radial",
  marker = list(colors = ~colors))
p_figure2c
save_image(p_figure2c, file = "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2c.pdf", height = 6000, width = 6000)
save_image(p_figure2c, file = "/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2c.png", height = 6000, width = 6000)

### Figure2d 
library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggh4x)
library(dplyr)
library(patchwork)
#Normal data
Normal_data <- fread("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Normal_high_level.tsv")
group_GTEx        <- c("GTEx_recurrent_StarF2019")   
group_DGD  <- c("DGD_PARALOGS")
group_HGNC    <- c("HGNC_GENEFAM")
group_Greger      <- c("Greger_Normal")
normal_markers <- c("GTEx_recurrent_StarF2019", "DGD_PARALOGS", "HGNC_GENEFAM",  "Greger_Normal")
gene_up_normal <- Normal_data  |>
  mutate(gene_pair = paste0(gene1,"--",gene2)) |>
  select(c(gene_pair,cancer_db_hits)) |>
  distinct()|>
  mutate(
    Group_GTEx = ifelse(rowSums(sapply(group_GTEx, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_DGD = ifelse(rowSums(sapply(group_DGD, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_HGNC = ifelse(rowSums(sapply(group_HGNC, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Greger = ifelse(rowSums(sapply(group_Greger, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Not_Reported = ifelse(rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) == 0, 1, 0),
  )  
GTEx <- gene_up_normal |> 
select(gene_pair, Group_GTEx) |>
filter(Group_GTEx == 1) |>
rename(Group = Group_GTEx)|>
mutate(Group = "Group_GTEx") 

DGD <- gene_up_normal |> 
select(gene_pair, Group_DGD)|>
filter(Group_DGD == 1) |>
rename(Group = Group_DGD) |>
mutate(Group = "Group_DGD") 

HGNC <- gene_up_normal |> 
select(gene_pair, Group_HGNC)|>
filter(Group_HGNC == 1) |>
rename(Group = Group_HGNC) |>
mutate(Group = "Group_HGNC") 

Greger <- gene_up_normal |> 
select(gene_pair, Group_Greger)|>
filter(Group_Greger == 1) |>
rename(Group = Group_Greger) |>
mutate(Group = "Group_Greger") 

Not_ReportedN <- gene_up_normal |> 
select("gene_pair", "Group_Not_Reported") |>
filter(Group_Not_Reported == 1) |>
rename(Group = Group_Not_Reported) |>
mutate(Group = "Group_Not_Reported") 

Group_All <- rbind(GTEx, DGD, HGNC, Greger, Not_ReportedN)
lists <- lapply(split(Group_All, Group_All$Group), function(x)x$gene_pair)
m <- make_comb_mat(lists)
cs = comb_size(m)
labels <- c("Group_GTEx", "Group_DGD","Group_HGNC", "Group_Greger","Group_Not_Reported")
colors <- c("#fcdaba", "#a7d2ba", "#d0cade", "#abc6e4", "#c39398")

pdf("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2e.pdf", width = 8, height = 4)
p_normal = draw(ComplexHeatmap::UpSet(m,                  
                      top_annotation = upset_top_annotation(m,  
                                                            ylim = c(0, max(cs)*1.5),
                                                            annotation_name_rot = 90,                                                      
                                                            annotation_name_side = "left",height = unit(4, "cm")),               
                                                            bg_col = rev(colors),              
                                                            set_order = labels)) 
od = column_order(p_normal)
decorate_annotation("intersection_size", {  
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)})

dev.off()
png("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2e.png", width = 8, height = 4, units = "in", res = 600)
p_normal = draw(ComplexHeatmap::UpSet(m,                  
  top_annotation = upset_top_annotation(m,  
                                        ylim = c(0, max(cs)*1.5),
                                        annotation_name_rot = 90,                                                      
                                        annotation_name_side = "left",height = unit(4, "cm")),               
                                        bg_col = rev(colors),              
                                        set_order = labels)) 
od = column_order(p_normal)
decorate_annotation("intersection_size", {  
grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)})

dev.off()

### Tumor data
Tumor_data  <- fread("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2a_data/Tumor_high_level.tsv")
group_chimer      <- c("ChimerKB", "ChimerPub", "ChimerSeq", "chimerdb_omim", "chimerdb_pubmed")
group_tcga        <- c("DEEPEST2019", "YOSHIHARA_TCGA", "Larsson_TCGA", "TCGA_StarF2019", "TumorFusionsNAR2018", "GUO2018CR_TCGA")
group_ccle        <- c("Klijn_CellLines", "CCLE_StarF2019", "DepMap2023")
group_karyotype   <- "Mitelman"
group_reviews     <- c("Cosmic", "HaasMedCancer")
tumor_markers <- c("Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic", "YOSHIHARA_TCGA", 
                   "ChimerKB", "ChimerPub", "ChimerSeq", "Klijn_CellLines", "Larsson_TCGA", 
                   "CCLE_StarF2019", "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018", 
                   "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023")
gene_up_tumor <- Tumor_data |>
  mutate(gene_pair = paste0(gene1,"--",gene2)) |>
  select(c(gene_pair,cancer_db_hits)) |>
  distinct()|>
  mutate(
    Group_TCGA = ifelse(rowSums(sapply(group_tcga, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_CCLE = ifelse(rowSums(sapply(group_ccle, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Karyotype = ifelse(rowSums(sapply(group_karyotype, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Chimer = ifelse(rowSums(sapply(group_chimer, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Reviews = ifelse(rowSums(sapply(group_reviews, grepl, x = cancer_db_hits)) > 0, 1, 0),
    Group_Not_Reported = ifelse(rowSums(sapply(tumor_markers, grepl, x = cancer_db_hits)) == 0, 1, 0),
  )  
TCGA <- gene_up_tumor |> 
              select(gene_pair, Group_TCGA) |>
              filter(Group_TCGA == 1) |>
              rename(Group = Group_TCGA)|>
              mutate(Group = "Group_TCGA") 
CCLE <- gene_up_tumor |> 
  select(gene_pair, Group_CCLE)|>
  filter(Group_CCLE == 1) |>
  rename(Group = Group_CCLE) |>
  mutate(Group = "Group_CCLE") 
Karyotype <- gene_up_tumor |> 
  select(gene_pair, Group_Karyotype)|>
  filter(Group_Karyotype == 1) |>
  rename(Group = Group_Karyotype) |>
  mutate(Group = "Group_Karyotype") 
Chimer <- gene_up_tumor |> select(gene_pair, Group_Chimer )|>
  filter(Group_Chimer == 1) |>
  rename(Group = Group_Chimer) |>
  mutate(Group = "Group_Chimer") 
Reviews <- gene_up_tumor |> select(gene_pair, Group_Reviews)|>
  filter(Group_Reviews == 1) |>
  rename(Group = Group_Reviews) |>
  mutate(Group = "Group_Reviews") 
Not_ReportedT <- gene_up_tumor |> 
  select("gene_pair", "Group_Not_Reported") |>
  filter(Group_Not_Reported == 1) |>
  rename(Group = Group_Not_Reported) |>
  mutate(Group = "Group_Not_Reported") 

Group_All <- rbind(TCGA, CCLE, Karyotype, Chimer, Reviews, Not_ReportedT)
lists <- lapply(split(Group_All, Group_All$Group), function(x)x$gene_pair)
m <- make_comb_mat(lists)
cs = comb_size(m)
labels <- c("Group_TCGA", "Group_CCLE", "Group_Karyotype","Group_Chimer", "Group_Reviews","Group_Not_Reported")
colors <- c("#f3a59a",  "#a6ddea",  "#80d0c3", "#9eaac4", "#c2c8da", "#d8cdc1")
pdf("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2d.pdf", width = 8, height = 4)
p_tumor = draw(ComplexHeatmap::UpSet(m,                  
                      top_annotation = upset_top_annotation(m,  
                                                            ylim = c(0, max(cs)*1.5),
                                                            annotation_name_rot = 90,                                                      
                                                            annotation_name_side = "left",height = unit(4, "cm")),               
                                                            bg_col = rev(colors),              
                                                            set_order = labels)) 
od = column_order(p_tumor)
decorate_annotation("intersection_size", {  
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)})

dev.off()
png("/home/data2/Projects/Fusion/paper_plot/Figure2/Figure2d.png", width = 8, height = 4, units = "in", res = 600)
p_tumor = draw(ComplexHeatmap::UpSet(m,                  
  top_annotation = upset_top_annotation(m,  
                                        ylim = c(0, max(cs)*1.5),
                                        annotation_name_rot = 90,                                                      
                                        annotation_name_side = "left",height = unit(4, "cm")),               
                                        bg_col = rev(colors),              
                                        set_order = labels)) 
od = column_order(p_tumor)
decorate_annotation("intersection_size", {  
grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)})

dev.off()
 

### Method 

### Figure 3网站介绍和截图
gene1 <- data.frame()
gene2 <- data.frame()
gene3 <- data.frame()
for (a in unique(combined_data$id)) {
  current_data <- combined_data |> filter(id == a)
  cancer_name <- current_data$cancer_type |> unique() 
  sample_ids   <- current_data$Sample_ID
  patient_ids <- data.frame(Sample_ID = current_data$Sample_ID, Patient_ID = current_data$Patient_ID)
  matched_data <- immunotherapy[immunotherapy$Sample_ID %in% sample_ids, ]
  matched_data[,Patient_ID := ifelse(Sample_ID %in% patient_ids$Sample_ID, patient_ids$Patient_ID)]
  matched_gene1a <- matched_data |> select(gene1,tools)
  while (any(grepl(",", matched_gene1a$gene1))) {
    matched_gene1a <- matched_gene1a %>%
      separate_rows(gene1, sep = ",")
  }
  gene_count_per_tool1 <- matched_gene1a |> 
    group_by(tools) |> 
    summarize(n = n_distinct(gene1))
  matched_gene1 <- data.frame(id = a, n = gene_count_per_tool1$n, tools = gene_count_per_tool1$tools,cancer = cancer_name)
  gene1 <- rbind(gene1,matched_gene1)
  
  matched_gene2a <- matched_data |> select(gene2,tools) 
  while (any(grepl(",", matched_gene2a$gene2))) {
    matched_gene2a <- matched_gene2a %>%
      separate_rows(gene2, sep = ",")
  }
  gene_count_per_tool2 <- matched_gene2a |> 
    group_by(tools) |> 
    summarize(n = n_distinct(gene2))
  matched_gene2 <- data.frame(id = a, n = gene_count_per_tool2$n, tools = gene_count_per_tool2$tools,cancer = cancer_name)
  gene2 <- rbind(gene2,matched_gene2)
  
  matched_gene3a <- matched_data |> select(gene2,gene1,tools) 
  while (any(grepl(",", matched_gene3a$gene2))) {
    matched_gene3a <- matched_gene3a %>%
      separate_rows(gene2, sep = ",")
  }
  while (any(grepl(",", matched_gene3a$gene1))) {
    matched_gene3a <- matched_gene3a %>%
      separate_rows(gene1, sep = ",")
  }
  matched_gene3a <- matched_gene3a %>%
    mutate(Fusion_name = paste(gene1, gene2, sep = "::"))
  
  gene_count_per_tool3 <- matched_gene3a |> 
    group_by(tools) |> 
    summarize(n = n_distinct(Fusion_name))
  matched_gene3 <- data.frame(id = a, n = gene_count_per_tool3$n, tools = gene_count_per_tool3$tools,cancer = cancer_name)
  gene3 <- rbind(gene3,matched_gene3)
  
}

gene1$id <- factor(gene1$id, levels = unique(gene1$id))
gene2$id <- factor(gene2$id, levels = unique(gene2$id))
gene3$id <- factor(gene3$id, levels = unique(gene3$id))

p1 <- ggplot(data=gene1,aes(x=id,y=n,fill=tools))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), 
        legend.title = element_blank()
  )+
  scale_fill_manual(values = c("arriba" = "#1c3c63","starfusion" = "#d3d5d4", "arriba,starfusion" = "#93c8c0")) +
  guides(fill = "none") +
  labs(x=NULL,y="Detected 3' gene counts") + 
  annotate(geom = "text",x=1.5,y=2200,label="BLCA")+
  geom_vline(xintercept = 2.5,lty="dashed")+   
  annotate(geom = "text",x=4,y=2200,label="KIRC")+
  geom_vline(xintercept = 5.5,lty="dashed")+
  annotate(geom = "text",x=11,y=2200,label="SKCM")+
  geom_vline(xintercept = 16.5,lty="dashed")+
  annotate(geom = "text",x=18,y=2200,label="NSCLC")+
  geom_vline(xintercept = 19.5,lty="dashed")+
  annotate(geom = "text",x=20,y=2200,label="SCLC")+
  geom_vline(xintercept = 20.5,lty="dashed")+
  annotate(geom = "text",x=21,y=2200,label="SARC")+
  geom_vline(xintercept = 21.5,lty="dashed")+
  annotate(geom = "text",x=22,y=2200,label="BRCA")+
  geom_vline(xintercept = 22.5,lty="dashed")+
  annotate(geom = "text",x=23,y=2200,label="STAD")+
  geom_vline(xintercept = 23.5,lty="dashed")+
  annotate(geom = "text",x=24,y=2200,label="HNSC")+
  geom_vline(xintercept = 24.5,lty="dashed")+
   annotate(geom = "text",x=25,y=2200,label="SGC")+
  geom_vline(xintercept = 25.5,lty="dashed")+
  annotate(geom = "text",x=26.5,y=2200,label="GBM")
  geom_vline(xintercept = 27.5,lty="dashed")   
#p2  
  p2 <-  ggplot(data=gene2,aes(x=id,y=n,fill=tools))+
    geom_bar(stat="identity",position = "dodge")+
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(), 
          legend.title = element_blank())+
    scale_fill_manual(values = c("arriba" = "#1c3c63","starfusion" = "#d3d5d4", "arriba,starfusion" = "#93c8c0")) +
    guides(fill = "none") +
    labs(x=NULL,y="Detected 5' gene counts") + 
    geom_vline(xintercept = 2.5,lty="dashed")+   
    geom_vline(xintercept = 5.5,lty="dashed")+
    geom_vline(xintercept = 16.5,lty="dashed")+
    geom_vline(xintercept = 19.5,lty="dashed")+
    geom_vline(xintercept = 20.5,lty="dashed")+
    geom_vline(xintercept = 21.5,lty="dashed")+
    geom_vline(xintercept = 22.5,lty="dashed")+
    geom_vline(xintercept = 23.5,lty="dashed")+
    geom_vline(xintercept = 24.5,lty="dashed")+
    geom_vline(xintercept = 25.5,lty="dashed")+
  geom_vline(xintercept = 27.5,lty="dashed")  
 # p3  
  p3 <- ggplot(data=gene3,aes(x=id,y=n,fill=tools))+
    geom_bar(stat="identity",position = "dodge")+
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle=60,hjust = 1,vjust = 1),
          legend.title = element_blank())+
    scale_fill_manual(values = c("arriba" = "#1c3c63","starfusion" = "#d3d5d4", "arriba,starfusion" = "#93c8c0")) +
    labs(x=NULL,y="Detected Fusion counts") + 
    geom_vline(xintercept = 2.5,lty="dashed")+   
    geom_vline(xintercept = 5.5,lty="dashed")+
    geom_vline(xintercept = 16.5,lty="dashed")+
    geom_vline(xintercept = 19.5,lty="dashed")+
    geom_vline(xintercept = 20.5,lty="dashed")+
    geom_vline(xintercept = 21.5,lty="dashed")+
    geom_vline(xintercept = 22.5,lty="dashed")+
    geom_vline(xintercept = 23.5,lty="dashed")+
    geom_vline(xintercept = 24.5,lty="dashed")+
    geom_vline(xintercept = 25.5,lty="dashed")+
  geom_vline(xintercept = 27.5,lty="dashed")
  
  # Arrange the plots using patchwork
  final_plot <- p1 / p2 / p3
  
  # Show the final plot
  final_plot  
  
  ggsave("C:/Users/Administrator/Desktop/write_paper/PAPER/2c.pdf", 
         plot = final_plot,         
         width = 14,                
         height = 7,                
         dpi = 1200,                
         units = "in",              
         device = "pdf",           
         bg = "white")  

### Figure 4摘一些分析



