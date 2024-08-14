# TODO: to check this
# Check ID between IOBR/phenotype info ---------------

map = c("" = "")

for (j in seq_along(map)) {
  i = map[j]
  path = file.path("../CircleRNADb", i)
  # Check if has sample ID issue
  f1 = paste0(names(i), "_TME.tsv.gz")
  f2 = file.path(path, paste0(i, "_TME.tsv.gz"))
  f11 = paste0(names(i), "_IOBR_signature.tsv.gz")
  f22 = file.path(path, paste0(i, "_IOBR_signature.tsv.gz"))
  d1 = data.table::fread(file.path(path, paste0(i, "_info.tsv")))
  d2 = data.table::fread(f1)
  if (any(d1$Run %in% d2$ID)) {
    message("ID can be found and matched, copying files...")
    file.copy(f1, f2)
    file.copy(f11, f22)
  } else {
    message("ID cannot be matched, cohort is ", i)
  }
}
