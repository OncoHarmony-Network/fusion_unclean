# NOTE: 该脚本部分内容已在我与彦昆讨论后弃用
# 相关内容已被拆分到 clean_tcga/cptac/target.R & gene_exon.R 等脚本中进行分解优化
# 王诗翔，2025年3月


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

filter_valid <- function(x) {
  !is.na(x) & x != "" & x != "--" & x != "Unknown"
}
# 筛选 TARGET_info 中的有效行
TARGET_info_valid <- TARGET_info %>%
  filter(filter_valid(Age) & filter_valid(OS_Time) & filter_valid(OS_Status) & filter_valid(Sex))
# 筛选 TARGET_xena_info 中的有效行(TODO:check 都是unique(),Sample重叠)
TARGET_xena_info_valid <- TARGET_xena_info %>%
  filter(filter_valid(Age) & filter_valid(OS_Time) & filter_valid(OS_Status) & filter_valid(Sex))
# 2. 根据 Sample_ID 合并两个数据框(NOTE:合并空，Sample_ID未相同)
merged_data <- TARGET_info_valid %>%
  inner_join(TARGET_xena_info_valid, by = "Sample_ID", suffix = c("_info", "_xena"))
# 3. 检查合并后的数据框中有效值是否一致
# 创建一个逻辑向量，检查每行的对应列是否一致
is_consistent <- merged_data %>%
  rowwise() %>%
  mutate(
    Age_consistent = Age_info == Age_xena,
    Sex_consistent = Sex_info == Sex_xena,
    OS_Time_consistent = OS_Time_info == OS_Time_xena,
    OS_Status_consistent = OS_Status_info == OS_Status_xena
  ) %>%
  ungroup()
# 4. 汇总结果
summary <- is_consistent %>%
  summarise(
    Age_consistent = all(Age_consistent, na.rm = TRUE),
    Sex_consistent = all(Sex_consistent, na.rm = TRUE),
    OS_Time_consistent = all(OS_Time_consistent, na.rm = TRUE),
    OS_Status_consistent = all(OS_Status_consistent, na.rm = TRUE),
  )
# > summary
# A tibble: 1 × 4
#  Age_consistent Sex_consistent OS_Time_consistent OS_Status_consistent
#  <lgl>          <lgl>          <lgl>              <lgl>               
#1 TRUE           TRUE           TRUE               TRUE 
## TARGET_sheet预处理(TODO:Sample_ID_short变量命名,TARGET复合型删除)
matched_data3 <- TARGET_sheet |> 
  select(c(`Case ID`, `Sample ID`, `Sample Type`)) |>
  rename("Sample_ID_short" = `Case ID`,
         "Sample_ID" = `Sample ID`,
         "Sample_type"  =  `Sample Type`) |> 
  separate_rows(c(Sample_ID, Sample_type), sep = ",\\s*")|> 
  mutate(Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+"))|> 
  unique()
#3226nrow
#TODO数据对准
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
    matched_data3[[col]] <- "Unknown"  
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
    project_id %in% c( "TARGET-ALL-P1", "TARGET-ALL-P2") ~ "TARGET-ALL",
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

detected_id_TARGET <- matched_data5$Sample_ID_short
##保存info
cohorts <- unique(matched_data5$cohort)
file_name1 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TARGET")
##TODO:race '--' 改为Unkown
if (!dir.exists(file_name1)) {
  dir.create(file_name1, recursive = TRUE)
} 
for (cohort1 in cohorts) {
  cohort_data <- subset(matched_data5, cohort == cohort1)
  name <- gsub("-", "_", cohort1)
  file_name <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TARGET/", name ,"_info.tsv")
  fwrite(cohort_data , file = file_name, sep = "\t")  
}

##Gene(TODO:不同cohort的id重叠)
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
#讨论后决定管子取样靠后面会稳定点, 比如TARGET-00-BM3897-14A-02R和TARGET-00-BM3897-14A-01R取01R
TARGET_TME <- TARGET_TME |>
  mutate( prefixes = str_extract(ID, ".*(?=-[^-]+$)"),  
          suffixes = str_extract(ID, "[^-]+$"),
          second_char = str_sub(suffixes, 2, 2),
          second_char = as.integer(second_char))
#TODO:重合样本未拆开
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
target_NORMAL_id <- target_NORMAL_sample_id_short$Sample_ID |> unique()

target_GNB_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-GNB") 
GNB_id <- target_GNB_sample_id_short$Sample_ID|> unique()

target_CELL_sample_id_short <- matched_data4 |> filter(cohort == "TARGET-CELL") 
CELL_id <- target_CELL_sample_id_short$Sample_ID|> unique()

result_TME <- result_TME |> 
  mutate(cohort = Project_id) |>
  mutate(cohort = case_when(
    Sample_ID_short %in% target_NORMAL_id ~ "TARGET-NORMAL",
    TRUE ~ cohort  
  )) |>
  mutate(cohort = case_when(
    Sample_ID_short %in% CELL_id ~ "TARGET-CELL",
    TRUE ~ cohort  
  ))|>
  mutate(cohort = case_when(
    Sample_ID_short %in% GNB_id ~ "TARGET-GNB",
    TRUE ~ cohort  
  )) |>
  filter(Sample_ID_short %in% detected_id_TARGET)   #确保都在GDC-SHEET里

cohorts_TME <- unique(result_TME$cohort)
file_name3 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/IOBR/TARGET/")
#ID换位Sample_ID,PROJECT去掉
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
  )) |>
  filter(Sample_ID_short %in% detected_id_TARGET)  #有些样本没有检出fusion

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
#TODO:GDC所有数据均应Sample_ID匹配cli,fusion,tme/IOBR。
#CLI:RUN,iobr:ID,fusion:Sample_ID
#TODO:有复合和普通样本重叠，选普通
##出现C3L-03262-04, C3L-03262-05和C3L-03262-01, C3L-03262-03俩行这种情况
##检查其", "字符俩边的是否一样，一样的话就是一个患者多个样本直接混一起了
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
#TRUE就是一个患者多个样本
##CPTAC_sheet预处理
matched_CPTAC1 <- CPTAC_sheet |> 
  select(c(`Case ID`, `Sample ID`, `Sample Type`, `Project ID`)) |>
  rename("Sample_ID_short" = `Case ID`,
         "Sample_ID" = `Sample ID`,
         "Sample_type"  =  `Sample Type`,
         "project_id" = `Project ID`) |> 
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

matched_CPTAC <- matched_CPTAC1 |>
                left_join(CPTACinfo, by = c("Sample_ID_short" = "Patient_ID")) |>
                rename("Patient_ID" = "Sample_ID_short") 

#matched_ids <- intersect(CPTACinfo$Patient_ID, matched_CPTAC$Sample_ID_short)
#setdiff(matched_CPTAC$Sample_ID_short,matched_ids)
#character(0)

##拆出CPTAC-NORMAL                 
matched_CPTAC <- matched_CPTAC |>
  mutate(cohort = case_when(
    grepl("Solid Tissue Normal", Sample_type) ~ "CPTAC-NORMAL",
    TRUE ~ cohort  
  )) |>
  mutate(Sample_ID_short = str_replace_all(Sample_ID, ", ", "-"))
CPTAC_NORMAL1 <- matched_CPTAC |>
                   filter(grepl("Solid Tissue Normal", Sample_type))
CPTAC_NORMAL_id <- CPTAC_NORMAL1$Sample_ID_short |> unique()
detected_id_CPTAC <- matched_CPTAC$Sample_ID_short |> unique()
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
  left_join(cohort_data, by = "Sample_ID_short")  |>
  filter(Sample_ID_short %in% detected_id_CPTAC) |> #有些样本没有检出fusion
  unique()
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
  left_join(cohort_data, by = "Sample_ID_short") |>
  filter(Sample_ID_short %in% detected_id_CPTAC) |> #有些样本没有检出fusion
  unique()

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
TCGA_sheet <- fread("/home/data2/Projects/Fusion/unclean/TCGA_sample_sheet.2024-05-14.tsv", sep = "\t") |> 
              select(c(`Sample ID`, `Sample Type`,`Project ID`)) 
colnames(TCGA_sheet) <- c("Sample_ID", "Sample_type","Project_id")

TCGA_sheet <- TCGA_sheet  |>
              mutate(Sample_ID = substr(Sample_ID, 1, 15)) |> 
              unique()

TCGA_info <- TCGA_info |>
  mutate(cohort = paste0("TCGA-",Cancer) )
 

TCGA_info1 <- TCGA_sheet |>
  left_join(TCGA_info, by = c("Sample_ID" = "Sample_ID")) |>
  mutate(
    cohort = Project_id,
    Sample_ID_short = Sample_ID,
    Patient_ID = substr(Sample_ID, 1, nchar(Sample_ID) - 3),  
    Cancer = gsub("TCGA-", "", Project_id)  
  )
  

TCGA_NORMAL <- TCGA_info1 |> 
               filter(Sample_type == "Solid Tissue Normal") 

NORMAL_id <- TCGA_NORMAL$Sample_ID |> unique()

detected_id_TCGA <- TCGA_info1$Sample_ID |> unique()

TCGA_info1 <- TCGA_info1 |> 
  mutate(
    Project_id = cohort,
    cohort = case_when(
      Sample_ID_short %in% NORMAL_id ~ "TCGA-NORMAL",
      TRUE ~ cohort)) 
cohorts <- unique(TCGA_info1$cohort)

file_name1 <- paste0("/home/data2/Projects/Fusion/unclean/gdc_split/cli/TCGA")

if (!dir.exists(file_name1)) {
  dir.create(file_name1, recursive = TRUE)
} 
for (cohort1 in cohorts) {
  cohort_data <- subset(TCGA_info1, cohort == cohort1)
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
         Sample_ID_short = substr(Sample_ID_short, 1, 15),
         cohort = case_when(Sample_ID_short %in% NORMAL_id ~ "TCGA-NORMAL",
                            TRUE ~ cohort)
  ) |>
  filter(Sample_ID_short %in% detected_id_TCGA)  |> #有些样本没有检出fusion
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
  )|>
  filter(Sample_ID_short %in% detected_id_TCGA)  |> #有些样本没有检出fusion
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
#metafusion_all
combined_id <- c(NORMAL_id, target_NORMAL_id, CPTAC_NORMAL_id)

meta_normal <- metafusion1 |>
               filter(Sample_ID %in% combined_id)
#> unique(meta_normal$sample_type)
#[1] "Normal"
#nrow(meta_normal)
#[1] 101141
#meta_new |>
#filter(sample_type == "Normal" ) |>
#nrow()
#[1] 114007

#meta_normal1 <-    meta_new |>
#                   filter(sample_type != "Normal" &  Sample_ID %in% combined_id) 
#> > meta_normal1$sample_type |> unique()
#[1] "Tumor" past not point right ,they should be normal sample
meta_new  <- metafusion1 |>
                    mutate(
                      sample_type = ifelse(
                        sample_type != "Normal" & Sample_ID %in% combined_id,
                        "Normal",
                        sample_type
                      )
                    )
cp_file_name <- "/home/data2/Projects/Fusion/unclean/gdc_split/fusion/metafusion_calls_matched.tsv"
write.table(meta_new, 
                      file = cp_file_name, 
                      sep = "\t", 
                      row.names = FALSE, 
                      quote = FALSE) 
#TCGA Fix :because of TCGA_sheet <- TCGA_sheet  |>
#              mutate(Sample_ID = str_sub(Sample_ID, 1, -2)) wrong
#should use mutate(Sample_ID = substr(Sample_ID, 1, nchar(Sample_ID) - 1))
#CPTAC fread bug? fwrite bug? fread bug.. use cp_file_name <- "/home/data2/Projects/Fusion/unclean/gdc_split/fusion/metafusion_calls_matched.tsv"
#meta_new <- fread(cp_file_name , sep = "\t")



#Not in Red Herrings
#metafusion
#TODO:检查适合注释的行就是少，数据库来源那个？肿瘤注释和正常注释？我先找一下俩reads标准。肿瘤分析多少没被注释
#Overlap假如分析之后，如果reads少的甩了，注释就多了，就star标明一下fusion检出可信标准，卡阈值。
#注释源来自？
cp_file_name <- "/home/data2/Projects/Fusion/unclean/gdc_split/fusion/metafusion_calls_matched.tsv"
meta_new <- fread(cp_file_name , sep = "\t")
Herring_Supplement <- meta_new |>
                      filter(sample_type == "Normal" & blocklist == "FALSE") #blocklist mean false positives
#> nrow(Herring_Supplement)
#[1] 109968
Herring_Supplement1 <- Herring_Supplement |>
                      filter(cancer_db_hits == "") 
#nrow(Herring_Supplement1)
#[1] 102449
Herring_Supplement2 <- Herring_Supplement |>
  filter(grepl("GTEx_recurrent_StarF2019", cancer_db_hits) |
  grepl("DGD_PARALOGS", cancer_db_hits) |
  grepl("HGNC_GENEFAM", cancer_db_hits) |
  grepl("BodyMap", cancer_db_hits) |
  grepl("Greger_Normal", cancer_db_hits) |
  grepl("ConjoinG", cancer_db_hits) |
  grepl("Babiceanu_Normal", cancer_db_hits))
#nrow(Herring_Supplement2)
#[1] 81
Herring_Supplement3 <-  Herring_Supplement |>
  filter(!grepl("GTEx_recurrent_StarF2019", cancer_db_hits) &
  !grepl("DGD_PARALOGS", cancer_db_hits) &
  !grepl("HGNC_GENEFAM", cancer_db_hits) &
  !grepl("BodyMap", cancer_db_hits) &
  !grepl("Greger_Normal", cancer_db_hits) &
  !grepl("ConjoinG", cancer_db_hits) &
  !grepl("Babiceanu_Normal", cancer_db_hits) & 
  cancer_db_hits != "")
#nrow(Herring_Supplement3)
#[1] 7438
#  meta_new |>
#  filter(sample_type == "Normal" & blocklist == "FALSE" & cancer_db_hits == "") |>
#  nrow()
#[1] 102449
#  meta_new |>
#  filter(sample_type == "Normal" & blocklist == "FALSE") |>
#  nrow()
#[1] 109968
#TODO 癌旁和normal
Herring_file_name <- "/home/data2/Projects/Fusion/fusiondb/Herring_Supplement.tsv"
write.table(Herring_Supplement, 
            file = Herring_file_name, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


#gene annotation
library(rtracklayer)
#metafusion
#meta_new <- fread("/home/data2/Projects/Fusion/fusiondb/metafusion_calls_matched.tsv", sep = "\t") will lead lose C3N-03849-03 after's -03
cp_file_name <- "/home/data2/Projects/Fusion/unclean/gdc_split/fusion/metafusion_calls_matched.tsv"
meta_new <- fread(cp_file_name , sep = "\t")  
meta_new_with_commas <- meta_new |>
                        filter(grepl(",", exon1) | grepl(",", exon2)) |>
                        separate(col = exon1, into = c("exon1_part1", "exon1_part2"), sep = ",", remove = FALSE) |>
                        separate(col = chr1, into = c("chr11", "chr12"), sep = ",", remove = FALSE,fill = "right") |>
                        mutate(
                            chr12 = ifelse(is.na(chr12), chr11, chr12),
                            exon1_part1_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon1_part1)),  
                            exon1_part1_end = as.integer(sub(".*-(\\d+)", "\\1", exon1_part1)),      
                            exon1_part2_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon1_part2)),  
                            exon1_part2_end = as.integer(sub(".*-(\\d+)", "\\1", exon1_part2))      
                          ) |>
                        separate(col = exon2, into = c("exon2_part1", "exon2_part2"), sep = ",", remove = FALSE) |>
                        separate(col = chr2, into = c("chr21", "chr22"), sep = ",", remove = FALSE,fill = "right") |>
                        mutate(
                            chr22 = ifelse(is.na(chr22), chr21, chr22),
                            exon2_part1_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon2_part1)),  
                            exon2_part1_end = as.integer(sub(".*-(\\d+)", "\\1", exon2_part1)),      
                            exon2_part2_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon2_part2)),  
                            exon2_part2_end = as.integer(sub(".*-(\\d+)", "\\1", exon2_part2))      
                          ) 
meta_new_without_commas <- meta_new |>
                        filter(!(grepl(",", exon1) | grepl(",", exon2))) |>
                        mutate(
                               exon1_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon1)),  
                               exon1_end = as.integer(sub(".*-(\\d+)", "\\1", exon1)),      
                               exon2_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon2)),  
                               exon2_end = as.integer(sub(".*-(\\d+)", "\\1", exon2))      
                              )
     
#Transcript
Cof <- fread("/home/data2/Projects/Fusion/fusiondb/TranscriptionFactor/Homo_sapiens_Cof.txt", sep = "\t")
TF <- fread("/home/data2/Projects/Fusion/fusiondb/TranscriptionFactor/Homo_sapiens_TF.txt", sep = "\t")

refer_Homo <- fread("/home/data2/Projects/Fusion/MetaFusion/reference_files/Homo_sapiens.gene_info")
refer1 <- refer_Homo |> 
          select(type_of_gene,Symbol) |> 
          dplyr::distinct()
meta_ref_gene_name_to_type <- setNames(refer1$type_of_gene, refer1$Symbol) 
refer2 <- refer_Homo |> 
          select(Synonyms,type_of_gene) |> 
          separate_rows(Synonyms, sep = "\\|") |>
          filter(Synonyms != "-")|> 
          dplyr::distinct()
meta_ref_gene_name_to_type2 <- setNames(refer2$type_of_gene, refer2$Synonyms)

#use refer_bed to annote gene
refer_bed <- fread("/home/data2/Projects/Fusion/MetaFusion/reference_files/hg38_gene.bed") |>
             mutate(V2 = as.integer(V2 + 1) ) |>
             select(c(V1,V2,V3,V8)) 

colnames(refer_bed) <- c("chr1", "exon1_start", "exon1_end", "gene_name")
result1 <- refer_bed |>
                dplyr::group_by(chr1, exon1_start, exon1_end) |>
                dplyr::distinct() |>
                dplyr::mutate(
                  detected_name1 = ifelse(n() > 1, paste(gene_name, collapse = ","), gene_name),
                  chr1 = as.character(chr1)
                )  |>
                dplyr::ungroup() |>
                dplyr::select(-gene_name) |>
                dplyr::distinct()
# found gene detected region gene
#without_commas
colnames(result1) <- c("chr1", "exon1_start", "exon1_end", "gene1_name")
meta_new_without_commas1 <- meta_new_without_commas |>
                           left_join(result1 , by = c("chr1", "exon1_start", "exon1_end")) 


colnames(result1) <- c("chr2","exon2_start", "exon2_end","gene2_name")
meta_new_without_commas1  <- meta_new_without_commas1 |>
                            left_join(result1 , by = c("chr2", "exon2_start", "exon2_end")) 

meta_new_without_commas1  <- meta_new_without_commas1 |>
                             select(-c("exon1_start", "exon1_end",
                                       "exon2_start", "exon2_end"))

#with_commas
colnames(result1) <- c("chr11","exon1_part1_start", "exon1_part1_end","detected_name1_part1")
meta_new_with_commas1  <- meta_new_with_commas |>
                         left_join(result1 , by = c("chr11","exon1_part1_start", "exon1_part1_end"))

colnames(result1) <- c("chr12","exon1_part2_start", "exon1_part2_end","detected_name1_part2")
meta_new_with_commas1  <- meta_new_with_commas1 |>
                         left_join(result1 , by = c("chr12","exon1_part2_start", "exon1_part2_end"))  

colnames(result1) <- c("chr21", "exon2_part1_start", "exon2_part1_end","detected_name2_part1")
meta_new_with_commas1  <- meta_new_with_commas1 |>
                         left_join(result1 , by = c("chr21", "exon2_part1_start", "exon2_part1_end"))

colnames(result1) <- c("chr22","exon2_part2_start", "exon2_part2_end","detected_name2_part2")
meta_new_with_commas1  <- meta_new_with_commas1 |>
                          left_join(result1 , by = c("chr22", "exon2_part2_start", "exon2_part2_end")) 

meta_new_with_commas1  <- meta_new_with_commas1 |>
                          select(-c("chr21", "exon2_part1_start", "exon2_part1_end",
                                    "chr22", "exon2_part2_start", "exon2_part2_end",
                                    "chr12","exon1_part2_start", "exon1_part2_end",
                                    "exon1_part1","exon1_part2","exon2_part1","exon2_part2",
                                    "chr11","exon1_part1_start", "exon1_part1_end")) 

meta_new_with_commas1 <- meta_new_with_commas1 |>
                                      mutate(
                                        gene1_name = ifelse(detected_name1_part1 == detected_name1_part2, detected_name1_part1, paste(detected_name1_part1, detected_name1_part2, sep = ",")),
                                        gene2_name = ifelse(detected_name2_part1 == detected_name2_part2, detected_name2_part1, paste(detected_name2_part1, detected_name2_part2, sep = ","))
                                      )   
meta_new_with_commas1 <- meta_new_with_commas1 |>
                         select(-c(detected_name1_part1, detected_name1_part2,detected_name2_part1, detected_name2_part2 )) 
 
#gene1

meta_new1 <- meta_new1 |>
  mutate(
    gene1_type = case_when(
      detected_name1 %in% names(gene_name_to_type) ~ gene_name_to_type[detected_name1],
      TRUE ~ NA_character_
    )
  ) |>
    mutate(
      gene1_type = case_when(
        is.na(gene1_type) & detected_name1 %in% names(meta_ref_gene_name_to_type) ~ meta_ref_gene_name_to_type[detected_name1],
        TRUE ~ gene1_type  
      )
    ) 
#gene2

meta_new1 <- meta_new1 |>
  mutate(
    gene2_type = case_when(
      detected_name2 %in% names(gene_name_to_type) ~ gene_name_to_type[detected_name2],
      TRUE ~ NA_character_
    )
  ) |>
    mutate(
      gene2_type = case_when(
        is.na(gene2_type) & detected_name2 %in% names(meta_ref_gene_name_to_type) ~ meta_ref_gene_name_to_type[detected_name2],
        TRUE ~ gene2_type  
      )
    )


    


#use refer_bed to annote gene



