library(fs)
library(data.table)

# 注意: CPTAC 和 TARGET 的 id 与 cli 后面可能需要
# 在包内进行匹配处理，此处没有处理。

# DONE:
# 0. IMvigor210 的 id 进行匹配修改
# ~~1. cohort id 与其他保持一致~~ 在R包中进行处理
# 2. (done) 标记分开 tumor 和 normal sample
# 3. (done) fusion id 需要根据 position 设置，因为数据是按不同队列 run 的
#   根据讨论，这一点按 metafusion 结果5'和3'的上下游正确数据，相同breakpoint不算同一个fusion
# 4. (不做) 处理过滤掉 fusion records 不在 metafusion 注释文件中的情况，由下游针对性处理

idmap = fread("EGAD00001003977_idmap.tsv")
idmap_vec = idmap$id
names(idmap_vec) = sub("GO29293_RNAseq_", "", idmap$fileid)

# # Ref genes in metafusion
# ref = data.table::fread("/home/data2/Projects/Fusion/MetaFusion/reference_files/Homo_sapiens.gene_info", header = TRUE)[`#tax_id` == 9606]
# ref_genes1 = ref$Symbol
# ref_genes2 = sapply(ref$Synonyms, stringr::str_split_1, pattern = "\\|") |> unlist() |> unique()
# ref_genes = c(ref_genes1, ref_genes2) |> na.omit() |> sort() |> unique()
# ref_genes = setdiff(ref_genes, "-")

# GDC fusion sample sheets ------------------------------------------------
tcga = fread("TCGA_sample_sheet.2024-05-14.tsv")
tcga = tcga |>
  dplyr::mutate(`Sample ID` = substr(`Sample ID`, 1, 15))

cptac = fread("CPTAC_gdc_sample_sheet.2023-11-15.tsv")
target = fread("TARGET_gdc_sample_sheet.2023-11-15.tsv")

gdc = rbind(tcga, cptac, target)
gdc = gdc |> dplyr::rename(Patient_ID = `Case ID`, Sample_ID = `Sample ID`, Sample_Type = `Sample Type`, Project_ID = `Project ID`)

fwrite(gdc, file = "/home/data2/Projects/Fusion/fusiondb/gdc_fusion_sample_sheet.tsv", sep = "\t")

# metafusion calls. final ----------------------------------------------
# 去掉 caller 本身的 3 个 filter
# https://gitea.zhoulab.ac.cn/OncoHarmony-Network/GeneFusionIO/issues/89

final_cluster.list = fs::dir_ls("/home/data2/Projects/Fusion/MetaFusion/RUNS_Results/", recurse = TRUE, glob = "*cluster.CANCER_FUSIONS")
final_cluster.list = final_cluster.list[!stringr::str_detect(final_cluster.list, "phs001994|phs002270|EGAD00001005500")]

names(final_cluster.list) = basename(dirname(final_cluster.list))

use_unfilter_data = function(f) {
  message("handling ", f)
  # f: path to metafusion final cancer fusion file
  f1 = list.files(dirname(f), pattern = "BLOCKLIST$",
                  full.names = TRUE) # BLOCKLIST file
  f2 = list.files(dirname(f), pattern = "cff.reformat.renamed.reann.WITH_SEQ.exons.cluster$",
                  full.names = TRUE)
  data = fread(f)
  data1 = fread(f1, header = FALSE)
  data2 = fread(f2)
  colnames(data1) = colnames(data2)
  data1$blocklist = TRUE

  data$breakpoint_1 = as.character(data$breakpoint_1)
  data1$breakpoint_1 = as.character(data1$breakpoint_1)
  data2$breakpoint_1 = as.character(data2$breakpoint_1)
  data$breakpoint_2 = as.character(data$breakpoint_2)
  data1$breakpoint_2 = as.character(data1$breakpoint_2)
  data2$breakpoint_2 = as.character(data2$breakpoint_2)

  data2 =
    dplyr::left_join(
      data2,
      dplyr::select(data1, gene1, gene2,
                    chr1, breakpoint_1, chr2, breakpoint_2, samples, blocklist),
      by = c("gene1", "gene2", "chr1", "breakpoint_1", "chr2", "breakpoint_2", "samples")
    ) |>
    dplyr::mutate(
      blocklist = ifelse(is.na(blocklist), FALSE, blocklist)
    ) |>
    data.table::as.data.table()

  data2$`#cluster_type` = NULL
  data2$cancer_db_hits = NULL
  data2$gene1_close_to_bnd = NULL
  data2$gene1_on_bnd = NULL
  data2$gene2_close_to_bnd = NULL
  data2$gene2_on_bnd = NULL
  colnames(data)[1] = "gene1"

  data_final =
  dplyr::left_join(
    data2,
    dplyr::select(data, gene1, gene2,
                  chr1, breakpoint_1, chr2, breakpoint_2, samples, cancer_db_hits),
    by = c("gene1", "gene2", "chr1", "breakpoint_1", "chr2", "breakpoint_2", "samples")
  ) |>
    dplyr::select(c(colnames(data), c("exon1", "exon2"), "blocklist")) |>
    data.table::as.data.table()
}

#final_cluster_dt = purrr::map(final_cluster.list, .f = fread) |> rbindlist(use.names = TRUE, fill = TRUE, idcol = "cohort")
#f = final_cluster.list[5]
final_cluster_dt = purrr::map(final_cluster.list, .f = use_unfilter_data) |> rbindlist(use.names = TRUE, fill = TRUE, idcol = "cohort")
final_cluster_dt

final_cluster_dt$disease = NULL
final_cluster_dt[, cohort := fcase(
  startsWith(cohort, "CPTAC"), substr(cohort, 7, 12),
  startsWith(cohort, "TARGET"), substr(cohort, 1, 10),
  startsWith(cohort, "RNA_"), toupper(substring(cohort, 5)),
  !(startsWith(cohort, "CPTAC") | startsWith(cohort, "TARGET") | startsWith(cohort, "RNA_")), cohort
)]

table(final_cluster_dt$cohort)
table(final_cluster_dt$tools, final_cluster_dt$cohort)

table(final_cluster_dt$sample_type)
final_cluster_dt[sample_type == "Normal", ]$cohort |> table()

final_cluster_dt$fusion_IDs = NULL

#colnames(final_cluster_dt)[2] = "gene1"

final_cluster_dt
table(final_cluster_dt$cohort)

final_cluster_dt2 = copy(final_cluster_dt)


# Create Fusion ID based on positions
final_cluster_dt2[, Pos_Index := paste(gene1, chr1, breakpoint_1, gene2, chr2, breakpoint_2, sep = "_")]
fusion_idx = final_cluster_dt2[ , .N, by = Pos_Index][order(Pos_Index)]
fusion_idx[, Fusion_ID := paste0("FID", .I)]

final_cluster_dt3 = final_cluster_dt2 |>
  dplyr::full_join(fusion_idx, by = "Pos_Index") |>
  tidyr::separate_longer_delim("samples", delim = ",") |>
  dplyr::rename(Sample_ID = samples) |>
  dplyr::mutate(Sample_ID = ifelse(startsWith(Sample_ID, "TCGA"), substr(Sample_ID, 1, 15), Sample_ID)) |>
  as.data.table()

nrow(final_cluster_dt)
nrow(final_cluster_dt2)
nrow(final_cluster_dt3)

final_cluster_dt3$N = NULL
final_cluster_dt3$Pos_Index = NULL

# 更新 sample_type
final_cluster_dt3$sample_type = "Tumor"
table(gdc$Sample_Type)
gdc_normal = gdc[grepl("Normal", Sample_Type)]

table(final_cluster_dt3$cohort)
final_cluster_dt3[cohort == "TARGET-AML"]
final_cluster_dt3[cohort == "CPTAC3"]

table(final_cluster_dt3[Sample_ID %in% gdc_normal$Sample_ID]$cohort)

final_cluster_dt3[, sample_type := ifelse(Sample_ID %in% gdc_normal$Sample_ID, "Normal", sample_type)]


table(final_cluster_dt3$sample_type)

final_cluster_IMvigor210 = final_cluster_dt3[cohort == "IMvigor210"]
final_cluster_IMvigor210[, Sample_ID := as.character(
  idmap_vec[Sample_ID]
)]
final_cluster_IMvigor210

final_cluster_dt4 = rbind(
  final_cluster_dt3[cohort != "IMvigor210"],
  final_cluster_IMvigor210
)[order(cohort, gene1, gene2)]

final_cluster_dt4

fwrite(final_cluster_dt4, file = "/home/data2/Projects/Fusion/fusiondb/metafusion_calls.tsv", sep = "\t")

table(final_cluster_dt4$inferred_fusion_type)

final_cluster_dt4$Sample_ID[grepl("_", final_cluster_dt4$Sample_ID)]
final_cluster_dt4$Sample_ID[grepl("_TARGET", final_cluster_dt4$Sample_ID)] |> unique()

any(grepl("_", final_cluster_dt4[cohort == "CPTAC2"]$Sample_ID))
sum(grepl("_", final_cluster_dt4[cohort == "CPTAC3"]$Sample_ID))
grep("_", final_cluster_dt4[cohort == "CPTAC3"]$Sample_ID, value = TRUE)

# target 数据存在几例数据合并了来自不同病人的数据，不合理，和彦昆讨论后过滤去除。
final_cluster_dt4 = final_cluster_dt4[!grepl("_TARGET", final_cluster_dt4$Sample_ID), ]
final_cluster_dt4[, Sample_ID := ifelse(
  cohort == "CPTAC3" & grepl("_", Sample_ID),
  gsub("_", "-", Sample_ID),
  Sample_ID
)]

fwrite(final_cluster_dt4, file = "/home/data2/Projects/Fusion/fusiondb/metafusion_calls.tsv", sep = "\t")

# Update cohort names for better mapping -------------------------------------
final_dt = data.table::fread("/home/data2/Projects/Fusion/fusiondb/metafusion_calls.tsv")
PHS000452_PD1 <- data.table::fread("/home/data2/Projects/Fusion/fusiondb/Clininfo/PHS000452_PD1_info.tsv")$Run

final_dt[Sample_ID %in% PHS000452_PD1]$cohort = "PHS000452_PD1"
final_dt[cohort == "PHS000452"]$cohort = "PHS000452_CTLA4"

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

final_dt[, cohort := ifelse(cohort %in% names(COHORT_MAP), COHORT_MAP[cohort], cohort)]
table(final_dt$cohort)
sum(is.na(final_dt$cohort))
sum(final_dt$cohort == "")


# update PRJNA795330 sample types
final_dt[cohort == "PRJNA795330"]$sample_type = "Blood"
table(final_dt$sample_type)

final_dt$cohort[startsWith(final_dt$cohort, "TARGET")] |> unique()
target = data.table::fread("/home/data2/Projects/Fusion/GeneFusionPipelines/data_cff/TARGET/TARGET_gdc_sample_sheet.2023-11-15.tsv")
target
table(target$`Project ID`)
table(target$`Sample Type`)

data.table::fwrite(final_dt, file = "/home/data2/Projects/Fusion/fusiondb/metafusion_calls_cohort_update.tsv", sep = "\t")

# Generating data for app --------------------------------------------------
system("rm -rf /home/data/Fusion/fusiondb_v2/fusion_calls_arrow")
library(arrow)
calls = data.table::fread("/home/data/Fusion/fusiondb_v2/metafusion_calls.tsv")

table(calls$cohort)
table(calls$tools)
table(calls$blocklist)
table(calls$sample_type)

calls |>
  dplyr::group_by(cohort, inferred_fusion_type, sample_type) |>
  write_dataset(file.path("/home/data/Fusion/fusiondb_v2", "fusion_calls_arrow"), format = "arrow")

