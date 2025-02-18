library(data.table)
library(dplyr)
library(stringr)
#TARGET临床信息处理
metafusion <- fread("/home/data2/Projects/Fusion/fusiondb/metafusion_calls_cohort_update.tsv", sep = "\t")
TARGET_xena_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_xena_info.tsv", sep = "\t") 
TARGET_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_info.tsv", sep = "\t")  
#TARGET_xena_info 替换为缩写
replacement_table1 <- c(
"Non cancerous tissue"                        =   "TARGET-Normal",
"ALL"                                         =   "TARGET-ALL",                                
"AML"                                         =   "TARGET-AML",                            
"AML-IF"                                      =   "TARGET-AML-Ind-Fail",                             
"NBL"                                         =   "TARGET-NBL",                           
"OS"                                          =   "TARGET-OS",                             
"Ewing sarcoma"                               =   "TARGET-ES",                     
"WT"                                          =   "TARGET-WT",                             
"CCSK"                                        =   "TARGET-CCSK",                          
"RT"                                          =   "TARGET-RT", 
"CNS, ependymoma"                             =   "TARGET-CNS-EPN",                  
"CNS, glioblastoma (GBM)"                     =   "TARGET-CNS-GBM",              
"CNS, rhabdoid tumor"                         =   "TARGET-CNS-RT",                 
"CNS, low grade glioma (LGG)"                 =   "TARGET-CNS-LGG",                 
"CNS, medulloblastoma"                        =   "TARGET-CNS-MB",                    
"CNS, other"                                  =   "TARGET-CNS-other",          
"NHL, anaplastic large cell lymphoma"         =   "TARGET-ALCL",      
"NHL, Burkitt lymphoma (BL)"                  =   "TARGET-BL", 
"Rhabdomyosarcoma"                            =   "TARGET-RMS",   
"Soft tissue sarcoma, non-rhabdomyosarcoma"   =   "TARGET-NRMSSTS"  
)

matched1 <- TARGET_xena_info$disease_code %in% names(replacement_table1)
TARGET_xena_info <- TARGET_xena_info |>
  mutate(disease_code = case_when(
    matched1 ~ replacement_table1[disease_code],
    TRUE ~ disease_code
  ))
matched_data <- TARGET_xena_info[matched1, ]
matched_data1   <- matched_data |>
                      mutate(
                      Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+")
                      )
#TARGET_info 替换为缩写
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
TARGET_info <- TARGET_info |>
  mutate(primary_diagnosis = case_when(
    matched2 ~ replacement_table2[primary_diagnosis],
    TRUE ~ primary_diagnosis
  ))
matched_data2 <- TARGET_info[matched2, ] |> 
                   mutate(
                   Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+"))
  

#使用TARGET_xena_info第三个-之前  匹配   

#用TARGET_info第三个-之前  匹配 
meta_target <- metafusion |>
  filter(cohort %in% c("TARGET-AML", "TARGET-NBL"))

meta_target_normal <- meta_target |> 
  mutate(
    Sample_ID_short = str_extract(Sample_ID, "^[^-]+-[^-]+-[^-]+"))

matched_indices2 <- match(meta_target_normal$Sample_ID_short, matched_data2$Sample_ID_short) 
meta_target_normal2 <-  meta_target_normal |>
                        mutate(cohort = matched_data2$primary_diagnosis[matched_indices2[!is.na(matched_indices2)]]) |>
                        select(-Sample_ID_short)
meta_target_normal_final <- meta_target_normal2
#meta_target_normal_final <-  rbind(meta_target_normal21,meta_target_normal1,meta_target_normal22)              
#CPTAC
#"CPTAC-2","CPTAC-3"
CPTAC_info <- fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/CPTAC_info.tsv", sep = "\t")
CPTAC_combined <- fread("CPTAC_combined_with_cancer_annotated.tsv", sep = "\t")
CPTAC3_non_included <-  fread("CPTAC3_non_included.tsv", sep = "\t")
CPTAC_cli_paper_combined  <- fread("CPTAC_cli_paper_combined.tsv", sep = "\t")

#CPTAC3_non_included primary_diagnosis
cptac_replacement_table4 <- c(
  "Adenocarcinoma, NOS" = "CPTAC-ADC",
  "Renal cell carcinoma, NOS" = "CPTAC-RCC",
  "Glioblastoma" = "CPTAC-GBM",
  "Endometrioid adenocarcinoma, NOS" = "CPTAC-EAC",
  "Infiltrating duct carcinoma, NOS" = "CPTAC-IDC",
  "Glioblastoma multiforme" = "CPTAC-GBM-MF",
  "Oligodendroglioma, anaplastic" = "CPTAC-ODG-A",
  "'--" = "CPTAC-Unknown",
  "Squamous cell carcinoma, NOS" = "CPTAC-SCC",
  "Oligodendroglioma, NOS" = "CPTAC-ODG",
  "Unknown" = "CPTAC-Unknown"
)
meta_cptac  <- metafusion |>
  filter(cohort %in% c("CPTAC-2", "CPTAC-3"))
#sample_sheet切换为case_id才匹配的上info
sample_sheet <- fread("/home/data2/Projects/Fusion/GeneFusionPipelines/data_cff/CPTAC/CPTAC_gdc_sample_sheet.2023-11-15.tsv", sep = "\t") |>
                select(c(`Case ID`,`Sample ID`))
sample_sheet1 <- sample_sheet |>
                 mutate(`Sample ID` = substr(`Sample ID`, 1, 9)) |>
                 distinct()
replacement_sample <- setNames(sample_sheet$`Case ID`, sample_sheet$`Sample ID`)
meta_cptac <- meta_cptac |>
              mutate(Sample_ID_A = replacement_sample[Sample_ID])
meta_cptac1 <- meta_cptac |>
               filter(!is.na(Sample_ID_A)) |>
               mutate(Sample_ID = Sample_ID_A)  |>
               select(-Sample_ID_A)            
meta_cptac2 <- meta_cptac |>
               filter(is.na(Sample_ID_A)) |> 
               mutate(Sample_ID_A = substr(Sample_ID, 1, 9))
replacement_sample2 <- setNames(sample_sheet1$`Case ID`, sample_sheet1$`Sample ID`)
meta_cptac2 <- meta_cptac2 |>
               mutate(Sample_ID_B = replacement_sample2[Sample_ID_A]) |>
               mutate(Sample_ID = Sample_ID_B)  |>
               select(-c(Sample_ID_B,Sample_ID_A)) 
meta_cptac_sheet <- rbind(meta_cptac1,meta_cptac2)

#CPTAC_info Cancer
cptac_replacement_table2 <- c(
  "BRCA" = "CPTAC-BRCA",
  "CCRCC" = "CPTAC-CCRCC",
  "COAD" =  "CPTAC-COAD",
  "GBM" = "CPTAC-GBM",
  "HGSC" = "CPTAC-HGSC",
  "HNSCC" = "CPTAC-HNSCC",
  "LSCC" = "CPTAC-LSCC" ,
  "LUAD" = "CPTAC-LUAD",
  "PDAC" = "CPTAC-PDAC",
  "UCEC" = "CPTAC-UCEC",
  "GTEX-Brain" = "CPTAC-GTEX-Brain"
)

#替换CPTAC里的cohort行值
matched_CPTAC_info <- CPTAC_info$Cancer %in% names(cptac_replacement_table2)
CPTAC_info <- CPTAC_info |>
  mutate(Cancer = case_when(
    matched_CPTAC_info ~ cptac_replacement_table2[Cancer],
    TRUE ~ Cancer
  ))
matched_data_info <- CPTAC_info[matched_CPTAC_info, ]
#CPTAC_combined  primary_diagnosis
cptac_replacement_table1 <- c(
  "Infiltrating duct and lobular carcinoma" = "CPTAC-IDLC",
  "Serous adenocarcinoma, NOS" = "CPTAC-SAC",
  "Infiltrating duct carcinoma, NOS" = "CPTAC-IDC",
  "Adenocarcinoma, NOS" = "CPTAC-ADC",
  "Infiltrating lobular carcinoma, NOS" = "CPTAC-IFLC",
  "Infiltrating duct and mucinous carcinoma" = "CPTAC-IDMC",
  "Intraductal papillary carcinoma" = "CPTAC-IDPC",
  "Invasive lobular carcinoma" = "CPTAC-IVLC",
  "Mucinous carcinoma" = "CPTAC-MC",
  "Squamous cell carcinoma, NOS" = "CPTAC-SCC",
  "Renal cell carcinoma, NOS" = "CPTAC-RCC",
  "Glioblastoma" = "CPTAC-GBM",
  "Endometrioid adenocarcinoma, NOS" = "CPTAC-EAC",
  "Glioblastoma multiforme" = "CPTAC-GBM-MF",
  "Oligodendroglioma, anaplastic" = "CPTAC-ODG-A",
  "Oligodendroglioma, NOS" = "CPTAC-ODG",
  "'--" = "CPTAC-Unknown",
  "Not Reported" = "CPTAC-Unknown",
  "Unknown" = "CPTAC-Unknown"
)
matched_CPTAC_info <- CPTAC_combined$primary_diagnosis %in% names(cptac_replacement_table1)
CPTAC_combined  <- CPTAC_combined  |>
  mutate(primary_diagnosis = case_when(
    matched_CPTAC_info ~ cptac_replacement_table1[primary_diagnosis],
    TRUE ~ primary_diagnosis
  ))
matched_data_combined <-CPTAC_combined[matched_CPTAC_info, ]
#匹配
matched_indices_cptac <- match(meta_cptac_sheet$Sample_ID, matched_data_info$Sample_ID)
meta_cptac_sheet_normal <- meta_cptac_sheet[!is.na(matched_indices_cptac), ]
meta_cptac_sheet_normal <- meta_cptac_sheet_normal|>
                           mutate(cohort = matched_data_info$Cancer[matched_indices_cptac[!is.na(matched_indices_cptac)]])
meta_cptac_sheet_unnormal <- meta_cptac_sheet[is.na(matched_indices_cptac), ] |>
                             mutate(Sample_ID_A = sub(",.*", "", Sample_ID))
#contains_comma <- grepl(",", meta_cptac_sheet_unnormal$Sample_ID)
#meta_cptac_sheet_no_comma <- meta_cptac_sheet_unnormal |>
#                             filter(!contains_comma) 
#nrow(meta_cptac_sheet_no_comma)
#[1] 0
matched_indices_cptac2 <- match(meta_cptac_sheet_unnormal$Sample_ID_A, matched_data_info$Sample_ID)
meta_cptac_sheet_unnormal <- meta_cptac_sheet_unnormal |>
                             mutate(cohort = matched_data_info$Cancer[matched_indices_cptac2[!is.na(matched_indices_cptac2)]]) |>
                             select(-Sample_ID_A)
meta_cptac_last <- rbind(meta_cptac_sheet_normal,meta_cptac_sheet_unnormal)
#由于combined更多的是补充，所以要分的话直接手动匹配了
IDLC <- matched_data_combined |> filter(primary_diagnosis == "CPTAC-IDLC")

#metafusion_last
meta_remain <- metafusion |>
  filter(!(cohort %in% c("CPTAC-2", "CPTAC-3", "TARGET-AML", "TARGET-NBL")))
meta_new <- rbind(meta_remain,meta_cptac_last,meta_target_normal_final)
fwrite(meta_new, file = "/home/data2/Projects/Fusion/fusiondb/metafusion_calls_matched.tsv", sep = "\t")