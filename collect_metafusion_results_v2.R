library(fs)
library(data.table)

# 注意: CPTAC 和 TARGET 的 id 与 cli 后面可能需要
# 在包内进行匹配处理，此处没有处理。

idmap = fread("EGAD00001003977_idmap.tsv")
idmap_vec = idmap$id
names(idmap_vec) = sub("GO29293_RNAseq_", "", idmap$fileid)

# GDC fusion sample sheets ------------------------------------------------
tcga = fread("TCGA_sample_sheet.2024-05-14.tsv")
tcga = tcga |>
  dplyr::mutate(`Sample ID` = substr(`Sample ID`, 1, 15))

cptac = fread("CPTAC_gdc_sample_sheet.2023-11-15.tsv")
target = fread("TARGET_gdc_sample_sheet.2023-11-15.tsv")

gdc = rbind(tcga, cptac, target)
gdc = gdc |> dplyr::rename(Patient_ID = `Case ID`, Sample_ID = `Sample ID`, Sample_Type = `Sample Type`, Project_ID = `Project ID`)

fwrite(gdc, file = "/data/data2/Projects/Fusion/fusiondb_v4/gdc_fusion_sample_sheet.tsv", sep = "\t")

# metafusion calls. final ----------------------------------------------

final_cluster_dt = fread("/data/data2/Projects/Fusion/MetaFusion/RUNS_merged_Results_202506/RUNS_merged_v2_unconvert/final.n1.cluster.CANCER_FUSIONS")
colnames(final_cluster_dt)[1] = 'gene1'
final_cluster_dt$fusion_IDs = NULL

final_cluster_dt = final_cluster_dt |>
  tidyr::separate_longer_delim("samples", delim = ",") |>
  dplyr::rename(Sample_ID = samples) |>
  dplyr::mutate(disease = NULL) |>
  tidyr::separate("Sample_ID", into = c("cohort", "Sample_ID"), sep = "__") |>
  dplyr::mutate(Sample_ID = ifelse(startsWith(Sample_ID, "TCGA"), substr(Sample_ID, 1, 15), Sample_ID)) |>
  dplyr::mutate(cohort = ifelse(startsWith(Sample_ID, "TCGA"), paste0("TCGA-", cohort), cohort)) |>
  dplyr::mutate(tools = ifelse(tools %in% "starfusion,arriba", "arriba,starfusion", tools)) |>
  as.data.table()


table(final_cluster_dt$cohort)
table(final_cluster_dt$tools)

final_cluster_dt[, cohort := fcase(
  startsWith(cohort, "CPTAC"), cohort,
  cohort %in% c("AML", "NBL"), paste0("TARGET-", cohort),
  startsWith(cohort, "RNA_"), toupper(substring(cohort, 5)),
  !(startsWith(cohort, "CPTAC") | startsWith(cohort, "TARGET") | startsWith(cohort, "RNA_")), cohort
)]

table(final_cluster_dt$cohort)

# 更新 sample_type
final_cluster_dt$sample_type = "Tumor"
table(gdc$Sample_Type)
gdc_normal = gdc[grepl("Normal", Sample_Type)]

table(final_cluster_dt$cohort)
final_cluster_dt[cohort == "TARGET-AML"]
final_cluster_dt[cohort == "CPTAC3"]

table(final_cluster_dt[Sample_ID %in% gdc_normal$Sample_ID]$cohort)

any(grepl(",", gdc_normal$Sample_ID))
gdc_normal[grepl(",", gdc_normal$Sample_ID)]
gdc_normal[, Sample_ID := gsub(", ", "_", Sample_ID)]
final_cluster_dt[grepl("_", final_cluster_dt$Sample_ID) & cohort %in% c("CPTAC2", "CPTAC3")]

final_cluster_dt[, sample_type := ifelse(Sample_ID %in% gdc_normal$Sample_ID,
                                         "Normal", sample_type)]
# CPTAC有一些混合样本，使用的_分隔，在gdc sheet中使用的, 分割
final_cluster_dt[grepl("_", final_cluster_dt$Sample_ID) & cohort %in% c("CPTAC2", "CPTAC3") & sample_type == "Normal"]

table(final_cluster_dt$sample_type)

final_cluster_IMvigor210 = final_cluster_dt[cohort == "IMvigor210"]
final_cluster_IMvigor210[, Sample_ID := as.character(
  idmap_vec[Sample_ID]
)]
final_cluster_IMvigor210

final_cluster_dt = rbind(
  final_cluster_dt[cohort != "IMvigor210"],
  final_cluster_IMvigor210
)[order(cohort, gene1, gene2)]

final_cluster_dt

table(final_cluster_dt$inferred_fusion_type)
final_cluster_dt$Sample_ID[grepl("_TARGET", final_cluster_dt$Sample_ID)] |> unique()

# target 数据存在几例数据合并了来自不同病人的数据，不合理，和彦昆讨论后过滤去除。
final_cluster_dt = final_cluster_dt[!grepl("_TARGET", final_cluster_dt$Sample_ID), ]

# Update cohort names for better mapping
PHS000452_PD1 <- data.table::fread("/data/data2/Projects/Fusion/fusiondb/Clininfo/PHS000452_PD1_info.tsv")$Run

final_cluster_dt[Sample_ID %in% PHS000452_PD1]$cohort = "PHS000452_PD1"
final_cluster_dt[cohort == "PHS000452"]$cohort = "PHS000452_CTLA4"

COHORT_MAP <- c(
  "CPTAC3" = "CPTAC-3",
  "CPTAC2" = "CPTAC-2",
  "EGAD00001006282" = "EGAD00001006282",
  "HRA000524" = "HRA000524",
  "PRJEB23709" = "PRJEB23709",
  "PRJEB25780" = "PRJEB25780",
  "PRJNA312948" = "PRJNA312948",
  "PRJNA356761" = "PRJNA356761",
  "PRJNA476140" = "PRJNA476140",
  "PRJNA744780" = "PRJNA744780",
  "PRJNA940989" = "PRJNA940989",
  "IMmotion150" = "EGAD00001004183",
  "OAK" = "EGAD00001008549",
  "POPLAR" = "EGAD00001008548",
  "IMmotion151" = "EGAD00001006619",
  "IMvigor010" = "EGAD00001007575",
  "PHS000452" = "PHS000452",
  "PHS001038" = "PHS001038",
  "PHS001427" = "PHS001427",
  "PHS001464" = "PHS001464",
  "PHS001493" = "PHS001493",
  "PHS002176" = "PHS002176",
  "PHS003284" = "PHS003284",
  "PHS003316" = "PHS003316",
  "PHS001919" = "PHS001919",
  "IMvigor210" = "EGAD00001003977",
  "PHS001994" = "PHS001994",
  "PHS002270" = "PHS002270",
  "PRJNA482620" = "PRJNA482620",
  "PRJNA498500" = "PRJNA498500",
  "PRJNA923698" = "PRJNA923698",
  "PRJNA795330" = "PRJNA795330"
)

final_cluster_dt[, cohort := ifelse(cohort %in% names(COHORT_MAP), COHORT_MAP[cohort], cohort)]
table(final_cluster_dt$cohort)
sum(is.na(final_cluster_dt$cohort))
sum(final_cluster_dt$cohort == "")

# update PRJNA795330 sample types
final_cluster_dt[cohort == "PRJNA795330"]$sample_type = "Blood"
table(final_cluster_dt$sample_type)

final_cluster_dt$cohort[startsWith(final_cluster_dt$cohort, "TARGET")] |> unique()
target = data.table::fread("/data/data2/Projects/Fusion/GeneFusionPipelines/data_cff/TARGET/TARGET_gdc_sample_sheet.2023-11-15.tsv")
target
table(target$`Project ID`)
#table(target$`Sample Type`)

final_cluster_dt2 = copy(final_cluster_dt)
# Create Fusion ID based on positions
final_cluster_dt2[, Pos_Index := paste(gene1, chr1, breakpoint_1, gene2, chr2, breakpoint_2, sep = "_")]
fusion_idx = final_cluster_dt2[ , .N, by = Pos_Index][order(Pos_Index)]
fusion_idx[, Fusion_ID := .I]

final_cluster_dt3 = final_cluster_dt2 |>
  dplyr::full_join(fusion_idx, by = "Pos_Index") |>
  dplyr::mutate(N = NULL, Pos_Index = NULL) |>
  dplyr::arrange(Fusion_ID)

data.table::fwrite(final_cluster_dt3, file = "/data/data2/Projects/Fusion/fusiondb_v4/metafusion_calls.tsv", sep = "\t")

# Generating data for database app ---------------------------------------------
# TODO: 等彦昆更新一个 score 的 calls 版本再切分存储

# system("rm -rf /home/data/Fusion/fusiondb_v4/fusion_calls_arrow")
# library(arrow)
# calls = data.table::fread("/data/data2/Projects/Fusion/fusiondb_v4/metafusion_calls.tsv")
#
# table(calls$cohort)
# table(calls$tools)
# table(calls$sample_type)
#
# calls |>
#   dplyr::group_by(cohort, inferred_fusion_type, sample_type) |>
#   write_dataset(file.path("/home/data/Fusion/fusiondb_v4", "fusion_calls_arrow"), format = "arrow")



# Fullfill unmatched samples -----------------------------------------------
# https://gitea.zhoulab.ac.cn/OncoHarmony-Network/GeneFusionIO/issues/93
pro_cohorts = c("PRJEB23709", "PRJEB25780", "PRJNA356761") # "PRJNA498500" (SRR8112658)
for (i in pro_cohorts) {
  message(i)
  info = fread(
    glue::glue("/data/data2/Projects/Fusion/fusiondb_v4/Clininfo/{i}_info.tsv")
  )
  iobr = fread(
    glue::glue("/data/data2/Projects/Fusion/fusiondb_v4/IOBR/{i}_TME.tsv.gz")
  )
  # iobr2 = fread(
  #   glue::glue("/data/data2/Projects/Fusion/fusiondb/IOBR/{i}_TME.tsv.gz")
  # )
  deal_samples = iobr$ID[!iobr$ID %in% info$Run]
  if (length(deal_samples) > 1) {
    # PRJEB23709 IOBR 出现了重复 "ERR2208912.1" "ERR2208913.1"
    message("unmatched samples detected")
    print(deal_samples)
    deal_samples = setdiff(deal_samples, c("ERR2208912.1", "ERR2208913.1"))
    # back up
    file.copy(glue::glue("/data/data2/Projects/Fusion/fusiondb_v4/Clininfo/{i}_info.tsv"),
              glue::glue("/data/data2/Projects/Fusion/fusiondb_v4/Clininfo/{i}_info.tsv.bk"))
    # create new one
    info2 = info |> dplyr::add_row(
      Patient_ID = deal_samples,
      Sample_ID = deal_samples,
      Run = deal_samples
    )
    message("updating file")
    data.table::fwrite(
      info2,
      glue::glue("/data/data2/Projects/Fusion/fusiondb_v4/Clininfo/{i}_info.tsv"),
      sep = "\t")

  }
}

