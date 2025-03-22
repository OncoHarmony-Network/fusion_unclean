# This script is used to recreate arrow file list 
# based on recent updates from zyk

system("rm -rf /home/data2/Projects/Fusion/fusiondb_v2/fusion_calls_arrow")
library(arrow)
calls = data.table::fread("/home/data2/Projects/Fusion/fusiondb_v2/metafusion_calls_v2.tsv")

table(calls$cohort)
table(calls$tools)
table(calls$blocklist)
table(calls$sample_type)

calls |>
  dplyr::group_by(cohort, inferred_fusion_type, sample_type) |>
  write_dataset(file.path("/home/data2/Projects/Fusion/fusiondb_v2", "fusion_calls_arrow"), format = "arrow")
