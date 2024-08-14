outpath = "/home/data2/Projects/Fusion/fusiondb/IOBR"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(data.table)
library(fs)
library(dplyr)

fl = fs::dir_ls("./IOBR.all.signature/", glob = "*Rdata")
cohorts = sub(".IOBR.*", "", basename(fl)) |>
  stringr::str_remove("IOBR.all.signature_integration_tpm.") |>
  stringr::str_remove(".Rdata")

for (t in cohorts) {
  message("Processing ", t)
  f = fl[grepl(t, fl)]
  load(f)

  expr.all.signature_integration$Index = NULL

  rd = expr.all.signature_integration |>
    mutate(ID = gsub("\\.", "-", ID)) |>
    as.data.table()

  t2 = sub("_", "-", t)
  f2 = file.path(outpath, paste0(t2, "_IOBR_signature.tsv.gz"))

  message("Writing to ", f2)
  fwrite(rd, f2, sep = "\t")
}
