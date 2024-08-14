outpath = "/home/data2/Projects/Fusion/fusiondb/IOBR"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(data.table)
library(fs)
library(dplyr)

fl = fs::dir_ls("./IOBR.expr.tme_combine", glob = "*Rdata")
cohorts = sub(".IOBR.*", "", basename(fl)) |>
  stringr::str_remove("IOBR.expr.tme_combine.") |>
  stringr::str_remove(".Rdata")

for (t in cohorts) {
  message("Processing ", t)
  f = fl[grepl(t, fl)]
  load(f)
  colnames(expr.cibersort_abs)[-1] = paste0(colnames(expr.cibersort_abs)[-1], "-ABS")

  exprs = ls(pattern = "expr\\.")
  exprs = setdiff(exprs, "expr.tme_combine")
  exprs_list = purrr::map(exprs, ~eval(parse(text = .)))
  rd = purrr::reduce(exprs_list, full_join, by = "ID")
  rd = rd |>
    mutate(ID = gsub("\\.", "-", ID)) |>
    as.data.table()

  t2 = sub("_", "-", t)
  f2 = file.path(outpath, paste0(t2, "_TME.tsv.gz"))

  message("Writing to ", f2)
  fwrite(rd, f2, sep = "\t")
}
