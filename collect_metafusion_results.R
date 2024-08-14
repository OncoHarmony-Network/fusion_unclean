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

fwrite(gdc, file = "/home/data2/Projects/Fusion/fusiondb/gdc_fusion_sample_sheet.tsv", sep = "\t")

# metafusion calls. final ----------------------------------------------
# call filters 要不要调整？filters 源代码写的不好调整，除非自己修正

final_cluster.list = fs::dir_ls("/home/data2/Projects/Fusion/MetaFusion/RUNS_Results/", recurse = TRUE, glob = "*cluster.CANCER_FUSIONS")
final_cluster.list = final_cluster.list[!stringr::str_detect(final_cluster.list, "phs001994|phs002270|EGAD00001005500")]

names(final_cluster.list) = basename(dirname(final_cluster.list))
final_cluster_dt = purrr::map(final_cluster.list, .f = fread) |> rbindlist(use.names = TRUE, fill = TRUE, idcol = "cohort")
final_cluster_dt

final_cluster_dt$disease = NULL
final_cluster_dt[, cohort := fcase(
  startsWith(cohort, "CPTAC"), substr(cohort, 7, 12),
  startsWith(cohort, "TARGET"), substr(cohort, 1, 10),
  startsWith(cohort, "RNA_"), toupper(substring(cohort, 5)),
  !(startsWith(cohort, "CPTAC") | startsWith(cohort, "TARGET") | startsWith(cohort, "RNA_")), cohort
)]

table(final_cluster_dt$cohort)

table(final_cluster_dt$sample_type)
final_cluster_dt[sample_type == "Normal", ]$cohort |> table()

final_cluster_dt$fusion_IDs = NULL

colnames(final_cluster_dt)[2] = "gene1"

final_cluster_dt
table(final_cluster_dt$cohort)

final_cluster_dt2 = copy(final_cluster_dt)

# # 检查下是否存在完全 gene1/2 position 对调的情况
#
# # Same positioned fusion detected in gene KRT8
# for (i in unique(final_cluster_dt2$gene1)) {
#   cat("gene ", i, "\t")
#   #i = "NBAS"
#   #dt_gene[c(25, 34)]
#   #unique(final_cluster_dt2[gene1 %in% i | gene2 %in% i]) |> View()
#
#   #i = "KRT8"
#   #dt_gene[c(27, 186)]
#
#   dt_gene = final_cluster_dt2[gene1 %in% i | gene2 %in% i,
#                               list(brk1 = paste0(gene1, "_", chr1, "_", breakpoint_1), # gene1, "_",
#                                    brk2 =  paste0(gene2, "_", chr2, "_", breakpoint_2))] |> unique() # gene2, "_",
#   if (nrow(dt_gene) < 2) next()
#   for (j in 1:(nrow(dt_gene)-1)) {
#     for (k in (j+1):nrow(dt_gene)) {
#       message("j:", j, ", k:", k)
#       if (dt_gene$brk1[j] == dt_gene$brk2[k] & dt_gene$brk2[j] == dt_gene$brk1[k]) {
#         message("Same positioned fusion detected in gene ", i)
#         stop()
#       }
#     }
#   }
# }
#
# # >   dt_gene[c(25, 34)]
# # brk1               brk2
# # <char>             <char>
# #   1: NBAS_chr2_15558634 DDX1_chr2_15591553
# # 2: DDX1_chr2_15591553 NBAS_chr2_15558634
# #
# # > dt_gene[c(27, 186)]
# # brk1                            brk2
# # <char>                          <char>
# #   1:  KRT8_chr12_52898771 AC121757.2,KRT86_chr12_52306206
# # 2: KRT86_chr12_52306206             KRT8_chr12_52898771
# #
# # 结论：存在！得先处理统一 fusion 的情况
# # 根据 https://academic.oup.com/view-large/figure/405597626/btab249f3.tif 是通过聚类基于 breakpoint 和 name 做了合并
# # 但队列间（队列内应该处理了？）还是不同

# 用基因+位点标记ID

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

# > table(final_cluster_dt2$inferred_fusion_type)
#
# CodingFusion         NoHeadGene           SameGene    TruncatedCoding TruncatedNoncoding
# 574101               1965              15561             123047              87671
#
# ReadThrough 已经被过滤掉了。

# metafusion calls. cis-sage ----------------------------------------------
# 注意它和metafusion的最终结果应该不是互斥的！！！
# 这个数据后续考虑单独做分析
# diff=100,000
# fusion.chr1 == fusion.chr2  and abs(bp1 - bp2) < diff and fusion.inferred_fusion_type != "CodingFusion" or fusion.inferred_fusion_type == "ReadThrough" or fusion.inferred_fusion_type == "SameGene":
# ai 解释：
#  这个条件表达式的逻辑可以这样理解：
# 如果chr1和chr2相同，并且两个断点的位置差小于diff，并且inferred_fusion_type不是"CodingFusion"，则条件为真。
# 或者，如果inferred_fusion_type是"ReadThrough"或"SameGene"，条件也为真。

cis_cluster.list = fs::dir_ls("/home/data2/Projects/Fusion/MetaFusion/RUNS_Results/", recurse = TRUE, glob = "*cis-sage*")
cis_cluster.list = cis_cluster.list[!stringr::str_detect(cis_cluster.list, "phs001994|phs002270|phs001464|EGAD00001005500")]


names(cis_cluster.list) = basename(dirname(cis_cluster.list))
cis_cluster_dt = purrr::map(cis_cluster.list, .f = fread) |> rbindlist(use.names = TRUE, fill = TRUE, idcol = "cohort")
cis_cluster_dt

cis_cluster_dt$disease = NULL
cis_cluster_dt[, cohort := fcase(
  startsWith(cohort, "CPTAC"), substr(cohort, 7, 12),
  startsWith(cohort, "TARGET"), substr(cohort, 1, 10),
  startsWith(cohort, "RNA_"), toupper(substring(cohort, 5)),
  !(startsWith(cohort, "CPTAC") | startsWith(cohort, "TARGET") | startsWith(cohort, "RNA_")), cohort
)]

cis_cluster_dt$fusion_IDs = NULL
cis_cluster_dt$`#cluster_type` = NULL

cis_cluster_dt
table(cis_cluster_dt$cohort)

any(grepl(",", cis_cluster_dt$samples))
all(is.na(cis_cluster_dt$cancer_db_hits))

cis_cluster_dt$cancer_db_hits = NULL

# 可能也有用，加个辅助 id 标记 Cis-SAGe
cis_cluster_dt[, Pos_Index := paste(gene1, chr1, breakpoint_1, gene2, chr2, breakpoint_2, sep = "_")]
fusion_idx = cis_cluster_dt[ , .N, by = Pos_Index][order(Pos_Index)]
fusion_idx[, Fusion_ID := paste0("CS", .I)]

cis_cluster_dt2 = cis_cluster_dt |>
  dplyr::full_join(fusion_idx, by = "Pos_Index") |>
  tidyr::separate_longer_delim("samples", delim = ",") |>
  dplyr::rename(Sample_ID = samples) |>
  dplyr::mutate(Sample_ID = ifelse(startsWith(Sample_ID, "TCGA"), substr(Sample_ID, 1, 15), Sample_ID)) |>
  as.data.table()


cis_cluster_dt2$N = NULL
cis_cluster_dt2$Pos_Index = NULL

# 更新 sample_type
cis_cluster_dt2$sample_type = "Tumor"
cis_cluster_dt2[, sample_type := ifelse(Sample_ID %in% gdc_normal$Sample_ID, "Normal", sample_type)]

table(cis_cluster_dt2$sample_type)

cis_cluster_IMvigor210 = cis_cluster_dt2[cohort == "IMvigor210"]
cis_cluster_IMvigor210[, Sample_ID := as.character(
  idmap_vec[Sample_ID]
)]
cis_cluster_IMvigor210

cis_cluster_dt3 = rbind(
  cis_cluster_dt2[cohort != "IMvigor210"],
  cis_cluster_IMvigor210
)[order(cohort, gene1, gene2)]

cis_cluster_dt3

fwrite(cis_cluster_dt3, file = "/home/data2/Projects/Fusion/fusiondb/metafusion_calls.cis-sage.tsv", sep = "\t")

# # CFF merged calls --------------------------------------------------------
# 暂时不管这个数据了，避免数据多套导致的复杂性，先关注核心数据和分析！
#
# cff.list = fs::dir_ls("/home/data2/Projects/Fusion/GeneFusionPipelines/data_cff_gene_name_all/", recurse = TRUE, glob = "*cff")
# cff.list = cff.list[!stringr::str_detect(cff.list, "phs001994|EGAD00001005500")]
#
# cff_dt = purrr::map(cff.list, .f = fread, header = FALSE) |> rbindlist()
# cff_dt
#
# table(cff_dt$V10)
#
# cff_dt = cff_dt |>
#   dplyr::mutate(V8 = ifelse(startsWith(V8, "TCGA"), substr(V8, 1, 15), V8)) |>
#   as.data.table()
#
# fwrite(cff_dt, file = "/home/data2/Projects/Fusion/fusiondb/all_fusion_calls.tsv", sep = "\t")
# # V15 是被 metafusion 统一过么？还是来自各自的 caller？
#
# table(cff_dt[V11 == "arriba"]$V15)
# table(cff_dt[V11 != "arriba"]$V15)
# # 看来是拼起来的
#
# # 所以对这种数据的标注信息都需要注意
# table(cff_dt$V9)
