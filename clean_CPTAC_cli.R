library(data.table)


# Clinical data from CPTAC paper ------------------------------------------


cli_list = fs::dir_ls("CPTAC_Clinical_meta_data_v1/")

cli = purrr::map(cli_list, fread) |> rbindlist(fill = TRUE, use.names = TRUE)
colnames(cli)
# 问题，这个第 2 行要跳过才行
# 另外这个同时包含 CPTAC-2/3 的样本，但总数少一些

cli2 = purrr::map(cli_list, function(x) {
  # 由于列名也不同，这里读两次，一次读 colname，一次读数据
  data1 = fread(x)
  data2 = fread(x, skip = 2, header = FALSE)
  # 更新列名
  colnames(data2) = colnames(data1)
  data2$Cancer = stringr::str_remove(basename(x), "_meta.txt")
  data2 |> as.data.table()
}) |> rbindlist(fill = TRUE, use.names = TRUE)
#colnames(cli2) = colnames(cli)
colnames(cli2)

cli2[Tumor == "Yes" & Normal == "Yes"] # 这个是标记患者是否同时 tumor / normal sample


cli_list2 = fs::dir_ls("meta_data_BCM/")

cli3 = purrr::map(cli_list2, function(x) {
  # 由于列名也不同，这里读两次，一次读 colname，一次读数据
  data1 = fread(x)
  data2 = fread(x, skip = 2, header = FALSE)
  # 更新列名
  colnames(data2) = colnames(data1)
  data2$Cancer = stringr::str_remove(basename(x), "_meta.txt")
  data2 |> as.data.table()
}) |> rbindlist(fill = TRUE, use.names = TRUE)
colnames(cli3)

all.equal(cli2, cli3)
all.equal(sort(colnames(cli2)), sort(colnames(cli3)))

cli4 = cli3[, colnames(cli2), with = FALSE]
all.equal(cli2, cli4)
identical(cli2, cli4)

dplyr::setdiff(cli2, cli4) |> View()
dplyr::setdiff(cli4, cli2) |> nrow()

# rbind(cli2[idx == "01OV007"],
#       cli4[idx == "01OV007"]) |> View()

# 2 个数据源目录的文件基本没有区别

cli.final = cli2 |>
  dplyr::rename(
    Patient_ID = idx,
    OS_Status = OS_event,
    OS_Time = OS_days,
    PFS_Status = PFS_event,
    PFS_Time = PFS_days
  ) |>
  dplyr::mutate(
    Sample_ID = Patient_ID,
    Sex = ifelse(Sex == "Female", "F", "M")
  ) |>
  dplyr::select(
    Patient_ID, Sample_ID, Cancer, Age, Sex, OS_Time, OS_Status, PFS_Time, PFS_Status, dplyr::everything()
  )
  # dplyr::select(
  #   - dplyr::starts_with("CIBERSORT"), - dplyr::starts_with("xCell")
  # )


fwrite(cli.final, file = "CPTAC_cli_paper_combined.tsv", sep = "\t")

# Clinical data from GDC --------------------------------------------------


system("
cd /home/data2/Projects/Fusion/unclean/gdc_cli/CPTAC
mkdir CPTAC-3 CPTAC-2
tar zxvf clinical.project-cptac-2.2024-08-07.tar.gz -C CPTAC-2
tar zxvf clinical.project-cptac-3.2024-08-07.tar.gz -C CPTAC-3
       ")

cptac2 = fread("gdc_cli/CPTAC/CPTAC-2/clinical.tsv")
cptac3 = fread("gdc_cli/CPTAC/CPTAC-3/clinical.tsv")

which(colnames(cptac2) %in% "residual_disease")
which(colnames(cptac3) %in% "residual_disease")
all.equal(cptac2[[119]], cptac2[[192]])
all.equal(cptac3[[119]], cptac3[[192]])
cptac2[[192]] = NULL
cptac3[[192]] = NULL

setdiff(cptac2$case_submitter_id, cli.final$Sample_ID)
setdiff(cptac3$case_submitter_id, cli.final$Sample_ID)
# > setdiff(cptac2$case_submitter_id, cli.final$Sample_ID)
# [1] "01OV049"   "26OV010"   "100004028" "05BR058"   "05BR055"   "17OV019"   "11BR069"   "02OV042"   "100004012" "05BR051"   "03BR012"   "02OV045"   "04OV041"
# [14] "05BR052"   "05BR031"   "1488"      "11OV009"   "02OV040"   "01OV002"   "100002921" "22OV001"   "13OV004"   "100003304" "02OV035"   "01OV046"   "01OV045"
# [27] "17OV034"
# > setdiff(cptac3$case_submitter_id, cli.final$Sample_ID)
# [1] "C3N-02672"                    "C3L-01929"                    "C3N-02682"                    "C3N-02012"                    "C3N-02070"
# [6] "C3N-02696"                    "C3L-02354"                    "C3N-03026"                    "C3L-01739"                    "C3N-02088"
# [11] "C3L-04213"                    "C3N-03205"                    "C3L-04081"                    "C3N-03791"                    "C3L-03679"
# [16] "C3L-02643"                    "C3N-04686"                    "C3L-01558"                    "C3L-03268"                    "C3L-01633"
# [21] "C3L-02220"                    "C3N-02439"                    "C3N-02028"                    "C3N-03800"                    "C3L-04392"
# [26] "C3L-02556"                    "C3N-02721"                    "C3L-01672"                    "C3L-00938"                    "C3L-02201"
# [31] "C3N-02146"                    "C3N-02639"                    "C3N-03788"                    "C3L-03733"                    "C3L-04037"
# [36] "C3N-01091"                    "C3L-02747"                    "C3N-01879"                    "C3N-02996"                    "C3N-02296"
# [41] "C3L-02553"                    "C3L-03632"                    "C3N-03755"                    "GTEX-NPJ7-0011-R10A-SM-HAKXW" "C3L-02746"
length(setdiff(cptac3$case_submitter_id, cli.final$Sample_ID))  # 586


View(cptac2[case_submitter_id %in% setdiff(cptac2$case_submitter_id, cli.final$Sample_ID)])
View(cptac3[case_submitter_id %in% setdiff(cptac3$case_submitter_id, cli.final$Sample_ID)])

data.table::fwrite(
  cptac3[case_submitter_id %in% setdiff(cptac3$case_submitter_id, cli.final$Sample_ID)],
  file = "CPTAC3_non_included.tsv", sep = "\t"
)

cptac3


# 根据 GDC 已有信息补全 paper 缺少的样本信息 ---------------------------------------------

# 可能标记了 cancer type 的列
# primary_diagnosis site_of_resection_or_biopsy tissue_or_organ_of_origin

cptac = rbind(cptac2, cptac3)[, .(case_submitter_id, primary_diagnosis, site_of_resection_or_biopsy, tissue_or_organ_of_origin)] |> unique()
cptac = dplyr::left_join(cptac, cli.final[, .(Sample_ID, Cancer)], by = c("case_submitter_id"="Sample_ID"))

setdiff(cli.final$Sample_ID, cptac$case_submitter_id)
# 有 29 例不在 GDC 数据中

table(cptac$Cancer, cptac$primary_diagnosis) |> pheatmap::pheatmap()
table(cptac$Cancer, cptac$site_of_resection_or_biopsy) |> pheatmap::pheatmap()
table(cptac$Cancer, cptac$tissue_or_organ_of_origin) |> pheatmap::pheatmap()

table(cptac$Cancer,
      paste(
        cptac$primary_diagnosis,
        cptac$site_of_resection_or_biopsy,
        cptac$tissue_or_organ_of_origin, sep = "_"
      )) |> pheatmap::pheatmap(show_colnames = FALSE)

# 一个个处理
cat(sub("_meta.txt", "", basename(names(cli_list))))

typelist = sub("_meta.txt", "", basename(names(cli_list)))
# BRCA CCRCC COAD GBM HGSC HNSCC LSCC LUAD PDAC UCEC
for (i in typelist) {
  message("check type ", i)
  cptac_part = cptac |> dplyr::filter(Cancer == i)
  message("  primary_diagnosis")
  table(cptac_part$Cancer, cptac_part$primary_diagnosis) |> print()
  message("  site_of_resection_or_biopsy")
  table(cptac_part$Cancer, cptac_part$site_of_resection_or_biopsy) |> print()
  message("  tissue_or_organ_of_origin")
  table(cptac_part$Cancer, cptac_part$tissue_or_organ_of_origin) |> print()
}
# BRCA (122): tissue_or_organ_of_origin "Breast, NOS"
# CCRCC (103): primary_diagnosis "Renal cell carcinoma, NOS"
#        tissue_or_organ_of_origin "Kidney, NOS"
# COAD (106): tissue_or_organ_of_origin "Colon, NOS", "Rectum, NOS"
# GBM (99): primary_diagnosis "Glioblastoma"
# HGSC (87): primary_diagnosis "Serous adenocarcinoma, NOS" AND
#  tissue_or_organ_of_origin Fallopian tube Ovary Peritoneum, NOS
# HNSCC (108): 其他情况 primary_diagnosis "Squamous cell carcinoma, NOS"
# Base of tongue, NOS Cheek mucosa Floor of mouth, NOS Gum, NOS
# HNSCC                   1            2                  18        3
#
# Head, face or neck, NOS Larynx, NOS Lip, NOS Oropharynx, NOS
# HNSCC                       2          47        4               4
#
# Overlapping lesion of lip, oral cavity and pharynx Tongue, NOS Tonsil, NOS
# HNSCC
# LSCC (108): primary_diagnosis "Squamous cell carcinoma, NOS" AND
#   tissue_or_organ_of_origin "lung" pattern
# LUAD (110): primary_diagnosis "Adenocarcinoma, NOS" AND
#   tissue_or_organ_of_origin "lung" pattern
# PDAC (105): tissue_or_organ_of_origin "pancreas" pattern
# UCEC (95): primary_diagnosis "Endometrioid adenocarcinoma, NOS"
sum(!is.na(cptac$Cancer))

# rbind(cptac2, cptac3)[,192] # 有重复列，需要回到前面处理 (fixed)

# 实际测试和更新
cptac_updated = rbind(cptac2, cptac3) |>
  unique() |>
  dplyr::mutate(
    Cancer = dplyr::case_when(
      tissue_or_organ_of_origin %in% "Breast, NOS" ~ "BRCA",
      primary_diagnosis %in% "Renal cell carcinoma, NOS" ~ "CCRCC",
      tissue_or_organ_of_origin %in% c("Colon, NOS", "Rectum, NOS") ~ "COAD",
      grepl("Glioblastoma", primary_diagnosis, ignore.case = TRUE) | grepl("Oligodendroglioma", primary_diagnosis, ignore.case = TRUE) ~ "GBM",
      grepl("GTEX", case_submitter_id, ignore.case = TRUE) ~ "GTEX-Brain",
      case_submitter_id %in% c("C3L-06912", "C3L-07212") ~ "GBM",
      #primary_diagnosis %in% "Glioblastoma" ~ "GBM",
      primary_diagnosis %in% "Serous adenocarcinoma, NOS" & tissue_or_organ_of_origin %in% c("Fallopian tube", "Ovary", "Peritoneum, NOS") ~ "HGSC",
      primary_diagnosis %in% "Squamous cell carcinoma, NOS" & grepl("lung", tissue_or_organ_of_origin, ignore.case = TRUE) ~ "LSCC",
      primary_diagnosis %in% "Adenocarcinoma, NOS" & grepl("lung", tissue_or_organ_of_origin, ignore.case = TRUE) ~ "LUAD",
      grepl("pancreas", tissue_or_organ_of_origin, ignore.case = TRUE) ~ "PDAC",
      primary_diagnosis %in% "Endometrioid adenocarcinoma, NOS" ~ "UCEC",
      TRUE ~ "HNSCC"
    )
  ) |> as.data.table()

for (i in "HNSCC") { # typelist
  message("check type ", i)
  cptac_part = cptac_updated |> dplyr::filter(Cancer == i)
  message("  primary_diagnosis")
  table(cptac_part$Cancer, cptac_part$primary_diagnosis) |> print()
  message("  site_of_resection_or_biopsy")
  table(cptac_part$Cancer, cptac_part$site_of_resection_or_biopsy) |> print()
  message("  tissue_or_organ_of_origin")
  table(cptac_part$Cancer, cptac_part$tissue_or_organ_of_origin) |> print()
}

table(cptac$Cancer)
table(cptac_updated$Cancer)

cptac_updated |>
  dplyr::filter(
    primary_diagnosis == "'--" | site_of_resection_or_biopsy == "'--" |
    tissue_or_organ_of_origin %in% c("'--", "Unknown")
  )
# 1-9 from GTEX-Brain
# C3L-06912 C3L-07212 Gliomas # https://portal.gdc.cancer.gov/cases/e6027623-a5b2-4d82-a9ca-f37c948dd3ed
cptac_updated$case_submitter_id[grepl("GTEX", cptac_updated$case_submitter_id)]

nrow(cptac_updated)

fwrite(cptac_updated, file = "CPTAC_combined_with_cancer_annotated.tsv", sep = "\t")

cptac_unmatched = cptac_updated[! case_submitter_id %in% cli.final$Sample_ID]
colnames(cptac_unmatched)
colnames(cli.final)[1:20]

cptac_unmatched2 = cptac_unmatched |>
  dplyr::select(
    case_submitter_id, age_at_diagnosis, gender,
    days_to_death, vital_status,
    dplyr::starts_with("ajcc"),
    dplyr::starts_with("days"),
    last_known_disease_status, progression_or_recurrence, residual_disease,
    tumor_grade, Cancer
  )

sum(cptac_unmatched2$days_to_last_follow_up != "'--")
sum(cptac_unmatched2$days_to_death != "'--")
sum(cptac_unmatched2$days_to_recurrence != "'--")
sum(cptac_unmatched2$days_to_last_known_disease_status != "'--")

sum(cptac_unmatched2$progression_or_recurrence != "'--")


sum(!is.na(cli.final$OS_Time))
sum(!is.na(cli.final$PFS_Time))
sum(cli.final$OS_Time >= cli.final$PFS_Time, na.rm = TRUE)
cli.final[cli.final$OS_Time < cli.final$PFS_Time] |> View()

sum(as.integer(cptac_unmatched2$days_to_last_follow_up) >= as.integer(cptac_unmatched2$days_to_last_known_disease_status), na.rm = TRUE)

cptac_unmatched3 = cptac_unmatched2 |>
  dplyr::rename(
    Sample_ID = case_submitter_id,
    Age = age_at_diagnosis,
    Sex = gender,
    OS_Time = days_to_last_follow_up,
    OS_Status = vital_status,
    PFS_Time = days_to_last_known_disease_status,
    PFS_Status = progression_or_recurrence,
    Histologic_Grade = tumor_grade,
    Path_Stage_pT = ajcc_pathologic_t,
    Path_Stage_pN = ajcc_pathologic_n,
    Stage = ajcc_pathologic_stage
  ) |>
  dplyr::select(
    Sample_ID, Age, Sex, OS_Time, OS_Status, PFS_Time, PFS_Status,
    Path_Stage_pN:Path_Stage_pT,
    Histologic_Grade, Cancer
  ) |>
  dplyr::mutate(
    Patient_ID = Sample_ID,
    Age = as.integer(round(as.integer(Age) / 365)),
    Sex = fcase(Sex %in% "female", "F",
                Sex %in% "male", "M", default = NA_character_),
    OS_Time = as.integer(OS_Time),
    OS_Status = fcase(OS_Status %in% "Alive", 0L,
                      OS_Status %in% "Dead", 1L, default = NA_integer_),
    PFS_Time = as.integer(PFS_Time),
    PFS_Status = fcase(PFS_Status %in% "no", 0L,
                       PFS_Status %in% "yes", 1L, default = NA_integer_),
    OS_Time = ifelse(OS_Time < 0, NA_integer_, OS_Time),
    Path_Stage_pN = fcase(
      Path_Stage_pN %in% "N0", "pN0",
      Path_Stage_pN %in% c("N1", "N1a"), "pN1",
      Path_Stage_pN %in% c("N2", "N2a"), "pN2",
      default = NA_character_
    ),
    Path_Stage_pT = fcase(
      Path_Stage_pT %in% "T0", "pT0",
      Path_Stage_pT %in% c("T1", "T1a", "T1b", "T1c"), "pT1",
      Path_Stage_pT %in% c("T2", "T2a", "T2b"), "pT2",
      Path_Stage_pT %in% c("T3", "T3a", "T3b"), "pT3",
      Path_Stage_pT %in% "T4", "pT4",
      default = NA_character_
    ),
    Stage = fcase(
      Stage %in% c("Stage I", "Stage IA", "Stage IA1", "Stage IA2", "Stage IA3", "Stage IB"), "Stage I",
      Stage %in% c("Stage II", "Stage IIA", "Stage IIB"), "Stage II",
      Stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IIIC1", "Stage IIIC2"), "Stage III",
      Stage %in% c("Stage IV", "Stage IVB"), "Stage IV",
      default = NA_character_
    ),
    Histologic_Grade = fcase(
      Histologic_Grade %in% "G1", "G1 Well differentiated",
      Histologic_Grade %in% "G2", "G2 Moderately differentiated",
      Histologic_Grade %in% "G3", "G3 Poorly differentiated",
      Histologic_Grade %in% "G4", "G4 Undifferentiated",
      default = NA_character_
    )
  )

cptac_unmatched3

# X: cannot be assessed

table(cptac_unmatched3$Path_Stage_pN)
table(cli.final$Path_Stage_pN)

table(cptac_unmatched3$Path_Stage_pT)
table(cli.final$Path_Stage_pT)

table(cptac_unmatched3$Stage)
table(cli.final$Stage)

table(cptac_unmatched3$Histologic_Grade)
table(cli.final$Histologic_Grade)



# Final output ------------------------------------------------------------

cli.final.2 = dplyr::bind_rows(
  cli.final, cptac_unmatched3
) |> as.data.table()

table(cli.final.2$Cancer)
sum(table(cli.final.2$Cancer))

# Race 全部加进去？
cli.final.2 = dplyr::left_join(
  cli.final.2, unique(rbind(cptac2, cptac3)[, .(case_submitter_id, race)]),
  by = c("Sample_ID"="case_submitter_id")
) |>
  dplyr::mutate(
    race = ifelse(race %in% c("Unknown", "not reported"), NA_character_,
                  race)
  ) |>
  as.data.table()

sum(is.na(cli.final.2$race))

fwrite(cli.final.2, file = "/home/data2/Projects/Fusion/fusiondb/Clininfo/CPTAC_info.tsv", sep = "\t")

