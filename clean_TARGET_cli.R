library(data.table)

cli1 = fread("TARGET_donor_allprojects_transfer_to_sample.gz")
cli2 = fread("TARGET_phenotype.gz")

target_cli = cli2 |>
  dplyr::select(-`_cohort`) |>
  dplyr::rename(Patient_ID = `_PATIENT`,
                Sample_ID = sample_id,
                disease_code = primary_disease_code,
                disease = `_primary_disease`,
                sample_type = `_sample_type`) |>
  dplyr::full_join(
    cli1 |> dplyr::select(
      -`_OS_UNIT`, -`_EVENT`, -`_TIME_TO_EVENT`
    ) |>
      dplyr::rename(
        Sample_ID = xena_sample,
        Sex = `_gender`,
        Age = `_age_at_diagnosis`,
        OS_Status = `_OS_IND`,
        OS_Time = `_OS`
      ), by = "Sample_ID"
  ) |>
  dplyr::mutate(
    Sex = ifelse(Sex == "Male", "M", "F"),
    OS_Status = ifelse(OS_Status == "Alive", 0L, 1L)
  ) |>
  dplyr::select(Patient_ID, Sample_ID, Age, Sex, OS_Time, OS_Status, dplyr::everything()) |>
  as.data.table()

cli1

target_cli

fwrite(target_cli, file = "/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_xena_info.tsv", sep = "\t")


# Using cli from GDC ------------------------------------------------------

system("
cd /home/data2/Projects/Fusion/unclean/gdc_cli/TARGET
mkdir ALL_P1 ALL_P2 ALL_P3 AML CCSK NBL OS RT WT
tar zxvf clinical.project-target-all-p1.2024-08-07.tar.gz -C ALL_P1
tar zxvf clinical.project-target-all-p2.2024-08-07.tar.gz -C ALL_P2
tar zxvf clinical.project-target-all-p3.2024-08-07.tar.gz -C ALL_P3
tar zxvf clinical.project-target-aml.2024-08-07.tar.gz -C AML
tar zxvf clinical.project-target-ccsk.2024-08-07.tar.gz -C CCSK
tar zxvf clinical.project-target-nbl.2024-08-07.tar.gz -C NBL
tar zxvf clinical.project-target-os.2024-08-07.tar.gz -C OS
tar zxvf clinical.project-target-rt.2024-08-07.tar.gz -C RT
tar zxvf clinical.project-target-wt.2024-08-07.tar.gz -C WT
       ")

gdc_target = rbindlist(
  list(
    fread("gdc_cli/TARGET/ALL_P1/clinical.tsv"),
    fread("gdc_cli/TARGET/ALL_P2/clinical.tsv"),
    fread("gdc_cli/TARGET/ALL_P3/clinical.tsv"),
    fread("gdc_cli/TARGET/AML/clinical.tsv"),
    fread("gdc_cli/TARGET/CCSK/clinical.tsv"),
    fread("gdc_cli/TARGET/NBL/clinical.tsv"),
    fread("gdc_cli/TARGET/OS/clinical.tsv"),
    fread("gdc_cli/TARGET/RT/clinical.tsv"),
    fread("gdc_cli/TARGET/WT/clinical.tsv")
  ),
  use.names = TRUE, fill = TRUE
)

table(gdc_target$project_id)
any(duplicated(gdc_target))
gdc_target = gdc_target[!duplicated(gdc_target)]

all.equal(gdc_target[[119]], gdc_target[[192]])
gdc_target[[192]] = NULL

sum(gdc_target$days_to_last_follow_up != "'--")
sum(gdc_target$days_to_death != "'--")

gdc_target2 = gdc_target |>
  dplyr::rename(
    Sample_ID = case_submitter_id,
    Age = age_at_diagnosis,
    Sex = gender,
    OS_Time = days_to_death,
    OS_Status = vital_status
  ) |>
  dplyr::select(
    Sample_ID, Age, Sex, OS_Time, OS_Status,
    project_id, race, primary_diagnosis
  ) |>
  dplyr::mutate(
    Patient_ID = Sample_ID,
    Age = as.integer(round(as.integer(Age) / 365)),
    Sex = fcase(Sex %in% "female", "F",
                Sex %in% "male", "M", default = NA_character_),
    OS_Time = as.integer(OS_Time),
    OS_Status = fcase(OS_Status %in% "Alive", 0L,
                      OS_Status %in% "Dead", 1L, default = NA_integer_),
    OS_Time = ifelse(OS_Time < 0, NA_integer_, OS_Time)
  ) |>
  as.data.table()

#unique()

# 存在大量 case level 重复
gdc_target2[duplicated(gdc_target2), ]
duplicated(gdc_target[case_submitter_id == "TARGET-20-PANTIV"])
gdc_target[case_submitter_id == "TARGET-20-PANTIV"] |> View() # 存在少量记录的更新？

gdc_target2 = gdc_target2[!duplicated(gdc_target2)]

#table(gdc_target2$PFS_Time)
any(duplicated(gdc_target2))
any(duplicated(gdc_target2$Sample_ID))

gdc_target2[duplicated(gdc_target2$Sample_ID)]
# 会存在不同 ALL 的 Project 中
gdc_dup = gdc_target2[Sample_ID %in% Sample_ID[duplicated(Sample_ID)]]
table(gdc_dup$project_id)
# TARGET-ALL-P1 TARGET-ALL-P2 TARGET-ALL-P3
# 10            57            47
#
# 这种先不处理，如果 ALL 合并分析考虑去重

colnames(gdc_target2)
gdc_target2 = gdc_target2[, c("Patient_ID", "Sample_ID", "Age", "Sex", "race", "primary_diagnosis", "OS_Time", "OS_Status", "project_id")]

fwrite(gdc_target2, file = "/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_info.tsv", sep = "\t")

sum(!is.na(gdc_target2$OS_Time))
sum(!is.na(cli1$`_OS`))

nrow(na.omit(unique(target_cli[, .(Patient_ID, OS_Time, OS_Status)])))
nrow(na.omit(unique(gdc_target2[, .(Patient_ID, OS_Time, OS_Status)])))
