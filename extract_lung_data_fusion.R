calls = data.table::fread("/home/data2/Projects/Fusion/fusiondb/metafusion_calls_cohort_update.tsv")

table(calls$cohort)

calls_lung = calls |>
    dplyr::filter(cohort %in% c("PHS002176", "PHS001464", "EGAD00001008548", "EGAD00001008549", "TCGA-LUAD", "TCGA-LUSC"))

nrow(calls_lung)

table(calls_lung$cohort)

data.table::fwrite(calls_lung, file = "/home/data2/Projects/Fusion/unclean/metafusion_calls_lung_cohorts.tsv.gz", sep = "\t")
