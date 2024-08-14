library(data.table)
library(dplyr)

# load("/home/data2/Projects/Fusion/unclean/TCGA_ALL_clinical.rdata")
#
# TCGA_ALL_clinical |>
#   rename(Patient_ID = bcr_patient_barcode,
#          Age = age_at_initial_pathologic_diagnosis,
#          Sex = gender,
#          Race = race,
#          )

# [Discrepancy] [Not Applicable]  [Not Available]
# [Discrepancy] [Not Available]       [Unknown]
# # [Discrepancy]                                        [Not Applicable]
# 29                                                     105
# [Not Available]                                         [Not Evaluated]
# 5140                                                       4
# [Unknown]


library(UCSCXenaShiny)

table(tcga_clinical_fine$Stage_ajcc)
tcga_clinical

tcga_clinical_fine
tcga_genome_instability
tcga_purity
tcga_subtypes
tcga_surv

nrow(tcga_clinical_fine)
nrow(unique(tcga_clinical_fine))

tcga_clinical = unique(tcga_clinical_fine) |>
  left_join(tcga_surv, by = c("Sample"="sample")) |>
  left_join(tcga_purity, by = c("Sample"="sample"))

tcga_clinical[551, ]$Sample

sum(startsWith(tcga_genome_instability$sample, "TCGA"))

tcga_clinical = tcga_clinical |>
  left_join(tcga_genome_instability[startsWith(tcga_genome_instability$sample, "TCGA") & !duplicated(tcga_genome_instability$sample), ],  by = c("Sample"="sample")) |>
  #left_join(tcga_genome_instability[duplicated(tcga_genome_instability$sample), ],  by = c("Sample"="sample")) |>
  left_join(tcga_subtypes,  by = c("Sample"="sampleID"))


sum(is.na(tcga_clinical$Sample))
sum(duplicated(tcga_clinical$Sample))

tcga_clinical[duplicated(tcga_clinical$Sample), ]$Sample

tcga_clinical.final = tcga_clinical[!duplicated(tcga_clinical$Sample), ] |>
  rename(Sample_ID = Sample, Sex = Gender, OS_Time = OS.time, OS_Status = OS, DSS_Time = DSS.time, DSS_Status = DSS, DFS_Time = DFI.time, DFS_Status = DFI, PFS_Time = PFI.time, PFS_Status = PFI) |>
  select(-cancer_type) |>
  mutate(Sex = ifelse(Sex == "FEMALE", "F", "M"),
         Patient_ID = substr(Sample_ID, 1, 12))

sum(is.na(tcga_clinical.final$Sample_ID))
sum(is.na(tcga_clinical.final$Cancer))
sum(duplicated(tcga_clinical.final$Sample_ID))

tcga_clinical.final = tcga_clinical.final |>
  select(Patient_ID, Sample_ID, everything()) |>
  as.data.table()

tcga_clinical.final

fwrite(tcga_clinical.final, file = "/home/data2/Projects/Fusion/fusiondb/Clininfo/TCGA_info.tsv", sep = "\t")
