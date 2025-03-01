library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

metafusion <- fread("/home/data2/Projects/Fusion/fusiondb/metafusion_calls_cohort_update.tsv", sep = "\t")

#TARGET
#cli
##确认GDC里有几个TARGET项目类型
TARGET_sheet <- fread("/home/data2/Projects/Fusion/unclean/TARGET_gdc_sample_sheet.2023-11-15.tsv", sep = "\t")
cohort_choose <- substr(TARGET_sheet$`Sample ID`, 1, 9) |> unique()
#"TARGET-20" "TARGET-00" "TARGET-21" "TARGET-30"
##TARGET_xena_info预处理匹配
#TODO:重合的是否重合
TARGET_xena_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_xena_info.tsv", sep = "\t") 
matched_data1   <- TARGET_xena_info |>
                   mutate(Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+")) |>
                   filter(grepl("^TARGET-(20|00|21|30)", Sample_ID_short))
  
##TARGET_info预处理匹配
TARGET_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_info.tsv", sep = "\t") 
replacement_table2 <- c(
  "Acute lymphocytic leukemia" =  "TARGET-ALL",             
  "'--"   = "TARGET-Unkown",                                  
  "Precursor B-cell lymphoblastic leukemia" = "TARGET-B-ALL", 
  "Not Reported"  = "TARGET-Unkown",                         
  "Acute myeloid leukemia, NOS" = "TARGET-AML",             
  "Clear cell sarcoma of kidney" = "TARGET-CCSK",           
  "Neuroblastoma, NOS"  = "TARGET-NBL",                     
  "Ganglioneuroblastoma" = "TARGET-GNB",                   
  "Osteosarcoma, NOS"   = "TARGET-OS",                    
  "Malignant rhabdoid tumor" = "TARGET-MRT",               
  "Wilms tumor" = "TARGET-WT"  
)
matched2 <- TARGET_info$primary_diagnosis %in% names(replacement_table2)
TARGET_info <- TARGET_info %>%
  mutate(cohort = case_when(
    matched2 ~ replacement_table2[primary_diagnosis],
    TRUE ~ as.character(NA)  # 如果不匹配，则设置为NA
  ))
matched_data2 <- TARGET_info[matched2, ]

matched_data2 <- TARGET_info |> 
                 mutate(Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+")) |>
                 filter(grepl("^TARGET-(20|00|21|30)", Sample_ID_short))

## TARGET_sheet预处理
matched_data3 <- TARGET_sheet |> 
                 select(c(`Case ID`, `Sample ID`, `Sample Type`)) |>
                 rename("Sample_ID_short" = `Case ID`,
                        "Sample_ID" = `Sample ID`,
                        "Sample_type"  =  `Sample Type`) |> 
                 separate_rows(c(Sample_ID, Sample_type), sep = ",\\s*")|> 
                 mutate(Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+"))|> 
                 unique()
            
## 使用TARGET_sheet更新Sample_ID,(原先也只是Patient_ID复制)并添加sample_type

# Step 1: 找到匹配的 Sample_ID_short
matched_ids <- intersect(matched_data2$Sample_ID_short, matched_data3$Sample_ID_short)
#setdiff(matched_data3$Sample_ID_short,matched_ids)
#character(0)
# 将数据框转换为 data.table
setDT(matched_data2)
setDT(matched_data3)
# Step 2: 获取需要更新的列名（除了 Sample_ID 和 Sample_ID_short）
cols_to_update <- setdiff(names(matched_data2), c("Sample_ID", "Sample_ID_short"))

# Step 3: 确保 matched_data3 中存在这些列，如果不存在则创建并填充为 NA
for (col in cols_to_update) {
  if (!col %in% names(matched_data3)) {  
    matched_data3[[col]] <- "Unkown"  
  }
}

# Step 4: 遍历每个匹配的 Sample_ID_short，并更新 matched_data3
matched_data3[matched_data2, on = .(Sample_ID_short), 
  (cols_to_update) := mget(paste0("i.", cols_to_update))]

matched_data3 <- matched_data3 |>
                 mutate(across(c(Age, OS_Time, OS_Status), ~ as.integer(.)))
## 使用TARGET_xena_info对TARGET_info的Sex, Age, OS_Time, OS_Status列进行补充
valid_data1 <- matched_data1 %>%
  filter(!is.na(Sex) & !is.na(Age) & !is.na(OS_Time) & !is.na(OS_Status))
# 将valid_data1转换为一个查找表，方便后续匹配
lookup_table <- valid_data1 %>%
  select(Sample_ID_short, Sex, Age, OS_Time, OS_Status) %>%
  distinct()
# 使用left_join将matched_data3与lookup_table合并，并用coalesce填充缺失值
matched_data4 <- matched_data3 %>%
  left_join(lookup_table, by = "Sample_ID_short", suffix = c("", "_valid")) %>%
  mutate(
    Sex = coalesce(na_if(Sex, ""), Sex_valid),  
    Age = coalesce(na_if(Age, NA_integer_), Age_valid),  
    OS_Time = coalesce(na_if(OS_Time, NA_integer_), OS_Time_valid),  
    OS_Status = coalesce(na_if(OS_Status, NA_integer_), OS_Status_valid)  
  ) %>%
  select(-c(Sex_valid, Age_valid, OS_Time_valid, OS_Status_valid))

##替换TARGET-Unkown为project_id
matched_data4 <- matched_data4 %>%
  mutate(cohort = case_when(
    cohort == "TARGET-Unkown" ~ project_id,
    TRUE ~ cohort
  )) %>%
  mutate(cohort = case_when(
    project_id %in% c("TARGET-ALL-P3", "TARGET-ALL-P1", "TARGET-ALL-P2") ~ "TARGET-ALL",
    TRUE ~ cohort
  )) %>%
    mutate(cohort = case_when(
      project_id %in% c("TARGET-ALL-P3") ~ "TARGET-ALAL",
      TRUE ~ cohort
    ))

##拆出细胞系

target_cell <- TARGET_sheet |> filter(`Sample Type` == "Cell Lines")
target_cell_sample_id_short <- sub("^([^-\n]+-[^-\n]+-[^-\n]+).*", "\\1", target_cell$`Sample ID`) |> unique()
# target_cell_sample_id_short
# [1] "TARGET-20-MUTZ3"         "TARGET-20-ML1"           "TARGET-20-MV411D1"       "TARGET-20-KasumiD1"      "TARGET-20-HEK293"       
# [6] "TARGET-20-MOLM14CBFGLIS" "TARGET-20-TF1"           "TARGET-20-IGROV1"        "TARGET-20-CMS"           "TARGET-20-WSUAML"       
# [11] "TARGET-20-KasumiAZAD5"   "TARGET-20-KasumiAZAD11"  "TARGET-20-REH"           "TARGET-20-MV411AZAD11"   "TARGET-20-MV411AZAD5"   
# [16] "TARGET-20-THP1"          "TARGET-20-OCIAML2"       "TARGET-20-HL60" 
#https://www.cellosaurus.org/
#MUTZ3：Acute myeloid leukemia
#ML1：	Adult acute myeloid leukemia
#MV411D1：biphenotypic B-myelomonocytic leukemia
#KasumiD1：Acute myeloblastic leukemia
#TF1：Acute erythroid leukemia 
#IGROV1：	Ovarian endometrioid adenocarcinoma
#WSUAML：Acute myeloid leukemia。
#KasumiAZAD5：can't find, so KasumiD1：Acute myeloblastic leukemia
#KasumiAZAD11：can't find, so KasumiD1：Acute myeloblastic leukemia
#REH：Childhood B acute lymphoblastic leukemia
#MV411AZAD11：can't find, so MV411D1：biphenotypic B-myelomonocytic leukemia
#MV411AZAD5：can't find, so MV411D1：biphenotypic B-myelomonocytic leukemia
#THP1：Childhood acute monocytic leukemia
#OCIAML2：Adult acute myeloid leukemia
#HL60：Adult acute myeloid leukemia
#HEK293：Vaccine production cell line.
#CML：Childhood acute megakaryoblastic leukemia
#MOLM14CBFGLIS：can't find, so MOLM14：Adult acute myeloid leukemia
cell_data <- metafusion[metafusion$Sample_ID %in% target_cell_sample_id_short]
#> nrow(cell_data)
# [1] 1853
cell_info <- TARGET_info[TARGET_info$Patient_ID %in% target_cell_sample_id_short]
# nrow(cell_info)
# [1] 18
#no os.etc
cell_info2 <- TARGET_xena_info[TARGET_xena_info$Sample_ID %in% unique(target_cell$`Sample ID`)]
#nrow(cell_info2)
#[1] 0

# 拆出"TARGET-CELL"
matched_data4 <- matched_data4 |>
  mutate(cohort = case_when(
    `Sample_type` %in% "Cell Lines" ~ "TARGET-CELL",
    TRUE ~ cohort  
  ))

##拆出"TARGET-NORMAL"
matched_data4 <- matched_data4 |>
  mutate(cohort = case_when(
    `Sample_type` %in% c("Blood Derived Normal", "Bone Marrow Normal") ~ "TARGET-NORMAL",
    TRUE ~ cohort  
  ))
##重新排列列顺序
new_order <- c("Patient_ID", "Sample_ID", "Age", "Sex", "race", "primary_diagnosis", 
               "OS_Time", "OS_Status", "Sample_type", "project_id", "cohort","Sample_ID_short")
matched_data4 <- matched_data4 |>
  select(all_of(new_order))

matched_data5 <- matched_data4 |>
                 mutate(Sample_ID_short = Sample_ID)
##保存info
cohorts <- unique(matched_data5$cohort)
file_name1 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TARGET")

if (!dir.exists(file_name1)) {
  dir.create(file_name1, recursive = TRUE)
} 
for (cohort1 in cohorts) {
cohort_data <- subset(matched_data5, cohort == cohort1)
name <- gsub("-", "_", cohort1)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TARGET/", name ,"_info.tsv")
fwrite(cohort_data , file = file_name, sep = "\t")  
}
#其km时删掉这种数据？
#TARGET-30-PAPEAV-01A, TARGET-30-PAPTFZ-01A
##Gene
metafusion1 <- metafusion |> 
               mutate(
                     cohort = case_when(
                     Sample_ID %in% matched_data5$Sample_ID_short ~ matched_data5$cohort[match(Sample_ID, matched_data5$Sample_ID_short)],
                     TRUE ~ cohort
                                       )
                    )  


##TME
NBL_TME <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-NBL_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-NBL")
AML_TME <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-AML_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-AML")
TARGET_TME <- rbind(NBL_TME,AML_TME)
#TODO:实际纳入数据库的只有几个Project_ID,纳入排除
##由于TARGET-00-BM3897-14A-01R和TARGET-00-BM3897-14A-02R，检查这样的有多少行

#suffix    <-   TARGET_TME |> 
#               select(ID) |>
#               mutate(ID = str_extract(ID , "^[^-]+-[^-]+-[^-]+-[^-]+")，) |> 
#               unique() 
#TARGET_TME|> nrow() 
#[1] 3344  
#suffix |> nrow() 
#[1] 3166         
#TARGET_TME|> nrow() - suffix |> nrow()
#[1] 178
#讨论后决定馆子取样靠后面会稳定点, 比如TARGET-00-BM3897-14A-02R和TARGET-00-BM3897-14A-01R取01R
TARGET_TME <- TARGET_TME |>
              mutate( prefixes = str_extract(ID, ".*(?=-[^-]+$)"),  
                      suffixes = str_extract(ID, "[^-]+$"),
                      second_char = str_sub(suffixes, 2, 2),
                      second_char = as.integer(second_char))

result_TME <- TARGET_TME |>
          group_by(prefixes) |>
          slice(which.max(second_char)) |>
          ungroup() |>
          rename(Sample_ID_short = prefixes) |>
          select(-c("suffixes", "second_char")) |>
          unique() 
          
#deleted_rows <- TARGET_TME |>
#  anti_join(result, by = c("ID", "prefixes", "suffixes", "second_char"))
# 删除数据框里里面没有TARGET-30-PAPTFZ-01A-01R-TARGET-30-PAPEAV-01A-01R这种情况
#找到特殊的Sample_id号匹配
target_NORMAL_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-NORMAL")  
NORMAL_id <- target_NORMAL_sample_id_short$Sample_ID |> unique()

target_GNB_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-GNB") 
GNB_id <- target_GNB_sample_id_short$Sample_ID|> unique()

target_CELL_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-CELL") 
CELL_id <- target_CELL_sample_id_short$Sample_ID|> unique()

result_TME <- result_TME |> 
  mutate(cohort = Project_id) |>
  mutate(cohort = case_when(
    Sample_ID_short %in% NORMAL_id ~ "TARGET-NORMAL",
    TRUE ~ cohort  
  )) |>
  mutate(cohort = case_when(
    Sample_ID_short %in% CELL_id ~ "TARGET-CELL",
    TRUE ~ cohort  
  ))|>
  mutate(cohort = case_when(
    Sample_ID_short %in% GNB_id ~ "TARGET-GNB",
    TRUE ~ cohort  
  )) 

cohorts_TME <- unique(result_TME$cohort)
file_name3 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TARGET/")

if (!dir.exists(file_name3)) {
  dir.create(file_name3, recursive = TRUE)
} 
for (cohort1 in cohorts_TME) {
cohort_data <- subset(result_TME, cohort == cohort1) |> select(-cohort)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TARGET/", cohort1 ,"_TME.tsv.gz")
write.table(cohort_data, 
            file = gzfile(file_name), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 
}  

##IOBR
NBL_IOBR_signature <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-NBL_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-NBL")
AML_IOBR_signature <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/TARGET-AML_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "TARGET-AML")
#拆出CELL和NORMAL和GNB,ALAL过滤掉了
TARGET_IOBR <- rbind(NBL_IOBR_signature,AML_IOBR_signature)

TARGET_IOBR <- TARGET_IOBR |>
               mutate(prefixes = str_extract(ID, ".*(?=-[^-]+$)"),  
                      suffixes = str_extract(ID, "[^-]+$"),
                      second_char = str_sub(suffixes, 2, 2),
                      second_char = as.integer(second_char))

result_IOBR <- TARGET_IOBR |>
               group_by(prefixes) |>
               slice(which.max(second_char)) |>
               ungroup() |>
               rename(Sample_ID_short = prefixes) |>
               select(-c("suffixes", "second_char")) |>
               unique() 

target_NORMAL_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-NORMAL")  
NORMAL_id <- target_NORMAL_sample_id_short$Sample_ID |> unique()

target_GNB_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-GNB") 
GNB_id <- target_GNB_sample_id_short$Sample_ID |> unique()

target_CELL_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-CELL") 
CELL_id <- target_CELL_sample_id_short$Sample_ID |> unique()

result_IOBR <- result_IOBR |> 
  mutate(cohort = Project_id) |>
  mutate(cohort = case_when(
    Sample_ID_short %in% NORMAL_id ~ "TARGET-NORMAL",
    TRUE ~ cohort  
  )) |>
  mutate(cohort = case_when(
    Sample_ID_short %in% CELL_id ~ "TARGET-CELL",
    TRUE ~ cohort  
  ))|>
  mutate(cohort = case_when(
    Sample_ID_short %in% GNB_id ~ "TARGET-GNB",
    TRUE ~ cohort  
  )) 

cohorts_IOBR <- unique(result_IOBR$cohort)
file_name2 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TARGET/")

if (!dir.exists(file_name2)) {
  dir.create(file_name2, recursive = TRUE)
} 
for (cohort1 in cohorts_IOBR) {
cohort_data <- subset(result_IOBR, cohort == cohort1) |> select(-cohort)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TARGET/", cohort1 ,"_IOBR_signature.tsv.gz")
write.table(cohort_data, 
            file = gzfile(file_name), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 
}  

#CPTAC
CPTAC_sheet <- fread("/home/data2/Projects/Fusion/GeneFusionPipelines/data_cff/CPTAC/CPTAC_gdc_sample_sheet.2023-11-15.tsv", sep = "\t")

##出现C3L-03262-04, C3L-03262-05和C3L-03262-01, C3L-03262-03俩行这种情况
##检查其", "字符俩边的是否一样，一样的话就是一个患者多个样本直接混一起了
rows_with_comma <- grepl(", ", CPTAC_sheet$`Case ID`)
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
#TRUE就是一个患者多个样本
##CPTAC_sheet预处理
matched_CPTAC <- CPTAC_sheet |> 
                 select(c(`Case ID`, `Sample ID`, `Sample Type`)) |>
                 rename("Sample_ID_short" = `Case ID`,
                        "Sample_ID" = `Sample ID`,
                        "Sample_type"  =  `Sample Type`) |> 
                 separate_rows(Sample_ID_short, sep = ",\\s*")|> 
                 unique() 
 
#matched_CPTAC$Sample_type |> unique()
#[1] "Primary Tumor"                                                            
#[2] "Solid Tissue Normal"                                                      
#[3] "Primary Tumor, Primary Tumor"                                             
#[4] "Primary Tumor, Primary Tumor, Primary Tumor"                              
#[5] "Solid Tissue Normal, Solid Tissue Normal, Solid Tissue Normal"            
#[6] "Solid Tissue Normal, Solid Tissue Normal"                                 
#[7] "Primary Tumor, Primary Tumor, Primary Tumor, Primary Tumor, Primary Tumor"
#[8] "Primary Tumor, Primary Tumor, Primary Tumor, Primary Tumor"
#"Solid Tissue Normal", "Primary Tumor"俩种，所以直接"Primary Tumor" , "Solid Tissue Normal"
##CPTAC_info预处理
CPTAC_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/CPTAC_info.tsv", sep = "\t")
 cptac_replacement_table2 <- c(
  "BRCA" = "CPTAC-BRCA",
  "CCRCC" = "CPTAC-KIRC",
  "COAD" =  "CPTAC-COAD",
  "GBM" = "CPTAC-GBM",
  "HGSC" = "CPTAC-HGSC",
  "HNSCC" = "CPTAC-HNSC",
  "LSCC" = "CPTAC-LUSC" ,
  "LUAD" = "CPTAC-LUAD",
  "PDAC" = "CPTAC-PDAC",
  "UCEC" = "CPTAC-UCEC",
  "GTEX-Brain" = "CPTAC-GTEX-Brain"
)     
 cptac_replacement_table3 <- c(
  "BRCA" = "CPTAC-BRCA",
  "CCRCC" = "CPTAC-CCRCC",
  "COAD" =  "CPTAC-COAD",
  "GBM" = "CPTAC-GBM",
  "HGSC" = "CPTAC-HGSC",
  "HNSCC" = "CPTAC-HNSCC",
  "LSCC" = "CPTAC-LUSC" ,
  "LUAD" = "CPTAC-LUAD",
  "PDAC" = "CPTAC-PDAC",
  "UCEC" = "CPTAC-UCEC",
  "GTEX-Brain" = "CPTAC-GTEX-Brain"
) 
matched_CPTAC_info <- CPTAC_info$Cancer %in% names(cptac_replacement_table2)
CPTAC_info <- CPTAC_info %>%
  mutate(cohort = case_when(
    matched_CPTAC_info ~ cptac_replacement_table2[Cancer],
    TRUE ~ as.character(NA)  # 如果不匹配，则设置为NA
  ))
CPTACinfo <- CPTAC_info[matched_CPTAC_info, ] |> select(-Sample_ID)

## 使用CPTAC_info更新TARGET_sheet,(CPTAC_info$Sample_ID原先也只是Patient_ID复制)并添加sample_type,
#由于如下原因，直接left_join
#matched_CPTAC |> nrow()
#[1] 2558
#matched_CPTAC |> select(Sample_ID_short) |> nrow()
#[1] 2558    

matched_ids <- intersect(CPTACinfo$Patient_ID, matched_CPTAC$Sample_ID_short)
#setdiff(matched_CPTAC$Sample_ID_short,matched_ids)
#character(0)

matched_CPTAC <- left_join(matched_CPTAC, CPTACinfo, by = c("Sample_ID_short" = "Patient_ID")) |>
                 rename("Patient_ID" = "Sample_ID_short") 
##拆出CPTAC-NORMAL                 
              
matched_CPTAC <- matched_CPTAC |>
                 mutate(cohort = case_when(
                 grepl("Solid Tissue Normal", Sample_type) ~ "CPTAC-NORMAL",
                 TRUE ~ cohort  
                 )) |>
                 mutate(Sample_ID_short = str_replace_all(Sample_ID, ", ", "-"))

  
##保存info
cohorts <- unique(matched_CPTAC$cohort)
file_name1 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/CPTAC")

if (!dir.exists(file_name1)) {
  dir.create(file_name1, recursive = TRUE)
} 
for (cohort1 in cohorts) {
cohort_data <- subset(matched_CPTAC, cohort == cohort1)
name <- gsub("-", "_", cohort1)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/CPTAC/", name ,"_info.tsv")
fwrite(cohort_data , file = file_name, sep = "\t")  
}
##Gene
##Gene
metafusion1 <- metafusion1 |> 
  mutate(
    cohort = case_when(
      Sample_ID %in% matched_CPTAC$Sample_ID_short ~ matched_CPTAC$cohort[match(Sample_ID, matched_CPTAC$Sample_ID_short)],
      TRUE ~ cohort
    )
  )  
##TME
CPTAC2_TME <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-2_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC2")
CPTAC3_TME <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-3_TME.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC3")
CPTAC_TME <- rbind(CPTAC2_TME, CPTAC3_TME)
#nrow(CPTAC_TME)
#[1] 2639
##matched_CPTAC$Sample_ID_short和CPTAC_TME$ID号匹配
# 匹配到的CPTAC_TME行
matched_first <- CPTAC_TME %>%
  filter(ID %in% matched_CPTAC$Sample_ID_short) |>
  mutate(Sample_ID_short = ID)
# 未匹配到的CPTAC_TME行
matched_last <- CPTAC_TME %>%
  filter(!(ID %in% matched_CPTAC$Sample_ID_short))
##第一种情况CPT008756_0001和CPT008756
#如果该行的str_count(matched_last$ID, "_") 为0则要求最后一个-之后的字符，
#为1则该行prefixes列的值为_之前的字符，该行suffixes列的值为第一个_之后的字符
#为2则该行prefixes列的值为第一个_之前的字符，该行suffixes列的值为第一个_之后的字符
#str_count(matched_last$ID, "_")
#[1] 1 2 0
matched_last <- matched_last %>%
  rowwise() %>%  # 确保逐行处理
  mutate(
    prefixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][1],  # 没有下划线，取最后一个"-"之前的部分
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][1],  # 有一个下划线，取第一个"_"之前的部分
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][1],  # 有两个下划线，取第一个"_"之前的部分
      TRUE ~ NA_character_
    ),
    suffixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][length(str_split(ID, "-")[[1]])],  # 没有下划线，取最后一个"-"之后的部分
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][2],  # 有一个下划线，取第一个"_"之后的部分
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][2],  # 有两个下划线，取第一个"_"之后的部分
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()  
#情况1
matched_lastf1 <- matched_last %>%
                  rename(Sample_ID_short = prefixes)  %>%
                  select(-suffixes)                   %>%
                  filter((Sample_ID_short %in% matched_CPTAC$Sample_ID_short))
matched_first <- rbind(matched_first, matched_lastf1)
#非情况1
matched_last1 <- matched_last %>%
  filter(!(prefixes %in% matched_CPTAC$Sample_ID_short))

##第二种情况能匹配上只是前面多了个X，（X1104806和1104806，X9f905736-f662-41d6-b3ac-16758d和9f905736-f662-41d6-b3ac-16758d）
matched_last1  <- matched_last1 %>%
                  mutate(prefixe_name = str_remove(ID, "^X"))  
#情况2
matched_lastf2 <- matched_last1 %>%
                  rename(Sample_ID_short = prefixe_name) %>%
                  select(-c(prefixes,suffixes)) %>%
                  filter((Sample_ID_short %in% matched_CPTAC$Sample_ID_short)) 
matched_first <- rbind(matched_first, matched_lastf2)
#非情况2
matched_last2 <- matched_last1 %>%
                 filter(!(prefixe_name %in% matched_CPTAC$Sample_ID_short))

##第三种情况其实就是一二混合（X161317f3-1a39-438e-af45-cd2502_D7和161317f3-1a39-438e-af45-cd2502）
matched_last2  <- matched_last2 %>%
                  mutate(prefixe_name = str_remove(prefixes, "^X")) 
#情况3    
matched_lastf3 <- matched_last2 %>%
                  rename(Sample_ID_short = prefixe_name) %>%
                  select(-c(prefixes,suffixes)) %>%
                  filter((Sample_ID_short %in% matched_CPTAC$Sample_ID_short))              
matched_first <- rbind(matched_first, matched_lastf3)
#非情况3
matched_last3 <- matched_last2 %>%
                 filter(!(prefixe_name %in% matched_CPTAC$Sample_ID_short))
#手动排除后发现没有出现C3L-02662-05-C3L-02662-03-C3L-02662-02颠倒顺序情况，TME的info对照列筛选完成
#nrow(matched_first)
#[1] 2558
#nrow(matched_CPTAC)
#[1] 2558
##匹配
#根据matched_first$Sample_ID_short 和 matched_CPTAC$Sample_ID_short的匹配，给matched_first添加新列cohort,值为matched_CPTAC$cohort,
cohort_data   <- matched_CPTAC %>% select(Sample_ID_short, cohort)
matched_first <- matched_first %>%
                 left_join(cohort_data, by = "Sample_ID_short")

cohorts_TME <- unique(matched_first$cohort)
file_name3 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/CPTAC/")

if (!dir.exists(file_name3)) {
  dir.create(file_name3, recursive = TRUE)
} 
for (cohort1 in cohorts_TME) {
cohort_data <- subset(matched_first, cohort == cohort1) |> select(-cohort)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/CPTAC/", cohort1 ,"_TME.tsv.gz")
write.table(cohort_data, 
            file = gzfile(file_name), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 
}  
               
#匹配链metafusion$Sample_ID == info$Sample_ID_short == iobr$Sample_ID_short
#IOBR
CPTAC2_Signature <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-2_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC2")
CPTAC3_Signature <- fread("/home/data2/Projects/Fusion/fusiondb/IOBR/CPTAC-3_IOBR_signature.tsv.gz", sep = "\t") |> mutate(Project_id = "CPTAC3")
CPTAC_Signature <- rbind(CPTAC2_Signature, CPTAC3_Signature)

##matched_CPTAC$Sample_ID_short和CPTAC_Signature$ID号匹配
# 匹配到的CPTAC_Signature行
matched_first <- CPTAC_Signature %>%
  filter(ID %in% matched_CPTAC$Sample_ID_short) |>
  mutate(Sample_ID_short = ID)
# 未匹配到的CPTAC_Signature行
matched_last <- CPTAC_Signature %>%
  filter(!(ID %in% matched_CPTAC$Sample_ID_short))
##第一种情况CPT008756_0001和CPT008756
#拆出prefixes列和suffixes列
#如果该行的str_count(matched_last$ID, "_") 为0则要求最后一个-之后的字符，
#为1则该行prefixes列的值为_之前的字符，该行suffixes列的值为第一个_之后的字符
#为2则该行prefixes列的值为第一个_之前的字符，该行suffixes列的值为第一个_之后的字符
#str_count(matched_last$ID, "_")
#[1] 1 2 0
matched_last <- matched_last %>%
  rowwise() %>%  # 确保逐行处理
  mutate(
    prefixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][1],  # 没有下划线，取最后一个"-"之前的部分
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][1],  # 有一个下划线，取第一个"_"之前的部分
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][1],  # 有两个下划线，取第一个"_"之前的部分
      TRUE ~ NA_character_
    ),
    suffixes = case_when(
      str_count(ID, "_") == 0 ~ str_split(ID, "-")[[1]][length(str_split(ID, "-")[[1]])],  # 没有下划线，取最后一个"-"之后的部分
      str_count(ID, "_") == 1 ~ str_split(ID, "_")[[1]][2],  # 有一个下划线，取第一个"_"之后的部分
      str_count(ID, "_") == 2 ~ str_split(ID, "_")[[1]][2],  # 有两个下划线，取第一个"_"之后的部分
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()  
#情况1
matched_lastf1 <- matched_last %>%
                  rename(Sample_ID_short = prefixes)  %>%
                  select(-suffixes)                   %>%
                  filter((Sample_ID_short %in% matched_CPTAC$Sample_ID_short))
matched_first <- rbind(matched_first, matched_lastf1)
#非情况1
matched_last1 <- matched_last %>%
  filter(!(prefixes %in% matched_CPTAC$Sample_ID_short))

##第二种情况能匹配上只是前面多了个X，（X1104806和1104806，X9f905736-f662-41d6-b3ac-16758d和9f905736-f662-41d6-b3ac-16758d）
matched_last1  <- matched_last1 %>%
                  mutate(prefixe_name = str_remove(ID, "^X"))  
#情况2
matched_lastf2 <- matched_last1 %>%
                  rename(Sample_ID_short = prefixe_name) %>%
                  select(-c(prefixes,suffixes)) %>%
                  filter((Sample_ID_short %in% matched_CPTAC$Sample_ID_short)) 
matched_first <- rbind(matched_first, matched_lastf2)
#非情况2
matched_last2 <- matched_last1 %>%
                 filter(!(prefixe_name %in% matched_CPTAC$Sample_ID_short))

##第三种情况其实就是一二混合（X161317f3-1a39-438e-af45-cd2502_D7和161317f3-1a39-438e-af45-cd2502）
matched_last2  <- matched_last2 %>%
                  mutate(prefixe_name = str_remove(prefixes, "^X")) 
#情况3    
matched_lastf3 <- matched_last2 %>%
                  rename(Sample_ID_short = prefixe_name) %>%
                  select(-c(prefixes,suffixes)) %>%
                  filter((Sample_ID_short %in% matched_CPTAC$Sample_ID_short))              
matched_first <- rbind(matched_first, matched_lastf3)
#非情况3
matched_last3 <- matched_last2 %>%
                 filter(!(prefixe_name %in% matched_CPTAC$Sample_ID_short))
#手动排除后发现没有出现C3L-02662-05-C3L-02662-03-C3L-02662-02颠倒顺序情况，TME的info对照列筛选完成
#nrow(matched_first)
#[1] 2558
#nrow(matched_CPTAC)
#[1] 2558
##匹配
#根据matched_first$Sample_ID_short 和 matched_CPTAC$Sample_ID_short的匹配，给matched_first添加新列cohort,值为matched_CPTAC$cohort,
cohort_data   <- matched_CPTAC %>% select(Sample_ID_short, cohort)
matched_first <- matched_first %>%
                 left_join(cohort_data, by = "Sample_ID_short")

cohorts_IOBR <- unique(matched_first$cohort)
file_name2 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/CPTAC/")

if (!dir.exists(file_name2)) {
  dir.create(file_name2, recursive = TRUE)
} 
for (cohort1 in cohorts_IOBR) {
cohort_data <- subset(matched_first, cohort == cohort1) |> select(-cohort)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/CPTAC/", cohort1 ,"_IOBR_signature.tsv.gz")
write.table(cohort_data, 
            file = gzfile(file_name), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 
}  

#TCGA 只用拆info

TCGA_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TCGA_info.tsv", sep = "\t")
#TCGA_sheet <- fread("/home/data2/Projects/Fusion/unclean/TCGA_sample_sheet.2024-05-14.tsv", sep = "\t")
TCGA_info <- TCGA_info |>
             mutate(cohort = paste0("TCGA-",Cancer),
                    Sample_ID_short = Sample_ID
             )
TCGA_NORMAL <- metafusion |> 
              filter(sample_type == "Normal" & startsWith(cohort, "TCGA")) 
            
NORMAL_id <- TCGA_NORMAL$Sample_ID |> unique()
 

TCGA_info <- TCGA_info |> 
             mutate(cohort = case_when(
             Sample_ID_short %in% NORMAL_id ~ "TCGA-NORMAL",
             TRUE ~ cohort))    

cohorts <- unique(TCGA_info$cohort)
file_name1 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TCGA")

if (!dir.exists(file_name1)) {
  dir.create(file_name1, recursive = TRUE)
} 
for (cohort1 in cohorts) {
cohort_data <- subset(TCGA_info, cohort == cohort1)
name <- gsub("-", "_", cohort1)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TCGA/", name ,"_info.tsv")
fwrite(cohort_data , file = file_name, sep = "\t")  
}  
#Gene
metafusion1 <- metafusion1 |> 
               mutate(
                     cohort = case_when(
                     Sample_ID %in% TCGA_info$Sample_ID_short ~ TCGA_info$cohort[match(Sample_ID, TCGA_info$Sample_ID_short)],
                     TRUE ~ cohort
                                       )
                    ) 
fwrite(metafusion1, file = "/home/data2/Projects/Fusion/unclean/gdc_split/fusion/metafusion_calls_matched.tsv", sep = "\t")
      
#IOBR
folder_path <- "/home/data2/Projects/Fusion/fusiondb/IOBR"

# 查询所有符合条件的文件
file_pattern <- "^TCGA.*_IOBR_signature.tsv.gz$"
file_list <- list.files(path = folder_path, pattern = file_pattern, full.names = TRUE)

# 初始化一个空的数据表
combined_data <- data.table()

# 遍历文件列表，逐个读取并处理
for (file in file_list) {
  # 读取文件
  file_data <- fread(file, sep = "\t")
  
  # 提取 Project_id
  file_name <- basename(file)  # 获取文件名
  project_id <- sub("^(TCGA-[^_]+).*", "\\1", file_name)  # 提取 TCGA-ACC 部分
  
  # 添加 Project_id 列
  file_data <- file_data %>%
    mutate(Project_id = project_id)
  
  # 将当前文件的数据追加到总数据表中
  combined_data <- rbind(combined_data, file_data)
}

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
               rename(Sample_ID_short = prefixes) |>
               mutate(cohort = Project_id,
                      Sample_ID_short = substr(Sample_ID_short, 1, nchar(Sample_ID_short) - 1),
                      cohort = case_when(Sample_ID_short %in% NORMAL_id ~ "TCGA-NORMAL",
                                         TRUE ~ cohort)
                      )  |>
               select(-c("suffixes", "second_chars")) |>
               unique() 
               
cohorts_IOBR <- unique(result_IOBR$cohort)
file_name2 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TCGA/")

if (!dir.exists(file_name2)) {
  dir.create(file_name2, recursive = TRUE)
} 
for (cohort1 in cohorts_IOBR) {
cohort_data <- subset(result_IOBR, cohort == cohort1) |> select(-cohort)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TCGA/", cohort1 ,"_IOBR_signature.tsv.gz")
write.table(cohort_data, 
            file = gzfile(file_name), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 
}  

#TME

folder_path <- "/home/data2/Projects/Fusion/fusiondb/IOBR"

# 查询所有符合条件的文件
file_pattern <- "^TCGA.*_TME.tsv.gz$"
file_list <- list.files(path = folder_path, pattern = file_pattern, full.names = TRUE)

# 初始化一个空的数据表
combined_data1 <- data.table()

# 遍历文件列表，逐个读取并处理
for (file in file_list) {
  # 读取文件
  file_data <- fread(file, sep = "\t")
  
  # 提取 Project_id
  file_name <- basename(file)  # 获取文件名
  project_id <- sub("^(TCGA-[^_]+).*", "\\1", file_name)  # 提取 TCGA-ACC 部分
  
  # 添加 Project_id 列
  file_data <- file_data %>%
    mutate(Project_id = project_id)
  
  # 将当前文件的数据追加到总数据表中
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
               rename(Sample_ID_short = prefixes) |>
               mutate(cohort = Project_id,
                      Sample_ID_short = substr(Sample_ID_short, 1, nchar(Sample_ID_short) - 1),
                      cohort = case_when(Sample_ID_short %in% NORMAL_id ~ "TCGA-NORMAL",
                                         TRUE ~ cohort)
                      )  |>
               select(-c("suffixes", "second_chars")) |>
               unique() 
               
cohorts_TME <- unique(result_TME$cohort)
file_name2 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TCGA/")

if (!dir.exists(file_name2)) {
  dir.create(file_name2, recursive = TRUE)
} 
for (cohort1 in cohorts_TME) {
cohort_data <- subset(result_TME, cohort == cohort1) |> select(-cohort)
file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TCGA/", cohort1 ,"_TME.tsv.gz")
write.table(cohort_data, 
            file = gzfile(file_name), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE) 
}  