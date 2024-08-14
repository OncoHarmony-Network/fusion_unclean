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

fwrite(target_cli, file = "/home/data2/Projects/Fusion/fusiondb/Clininfo/TARGET_info.tsv", sep = "\t")
