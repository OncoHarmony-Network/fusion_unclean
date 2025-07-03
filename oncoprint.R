
library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
all_data_Run_pro <- fread("/data/data2/Projects/Fusion/paper_plot.old/sample_id.tsv")
tumor1 <- fread("/data/data2/Projects/Fusion/paper_plot.old/oncoprint.tsv") |> distinct()
normal1 <- fread("/data/data2/Projects/Fusion/paper_plot.old/oncoprint_normal.tsv")
cancer_type <- all_data_Run_pro$cancer_type |> unique()
for(cancer_type1 in cancer_type){
  Run <- all_data_Run_pro |>
         filter(cancer_type == cancer_type1)
  cohort_id <- Run$cohort_id |> unique()
  if(cancer_type1 == "NORMAL"){
    tumor_onco <- normal1 |>
      filter(cohort %in% cohort_id)
  }else{
    tumor_onco <- tumor1 |>
      filter(cohort %in% cohort_id)
  }
  gene1 <- data.frame(gene = tumor_onco$gene1,
                   Sample_ID = tumor_onco$Sample_ID)
  gene2 <- data.frame(gene = tumor_onco$gene2,
                   Sample_ID = tumor_onco$Sample_ID)
  gene <- rbind(gene1, gene2) |>
         distinct()
  gene_count <- gene  |>
    group_by(gene) |>
    summarise(count = n()) |>
    arrange(desc(count))
  top_20_genes <- gene_count |>
                  head(20)
  matrix <- data.frame(
    Sample_ID = rep(Run$Run, each = nrow(top_20_genes)),
    gene = rep(top_20_genes$gene, times = length(Run$Run)),
    Fusion = ""
  )
  for (gene_top in top_20_genes$gene) {
    sample_genes <-  gene  |>
                     filter(gene == gene_top)
    sample <-  Run |>
               filter(Run %in% sample_genes$Sample_ID)
    matrix <- matrix %>%
      mutate(Fusion = case_when(
        Sample_ID %in% sample_genes$Sample_ID & gene == gene_top ~ "Fusion+",
        TRUE ~ Fusion
      ))
  }

  matrix_result <- matrix |>
                   pivot_wider(names_from = Sample_ID, values_from = Fusion, values_fill = "")|>
                   column_to_rownames(var = "gene")
  col = c("Fusion+" = "#01665E")
  alter_fun = list(
    background = alter_graphic("rect", fill = "#CCCCCC"),
    `Fusion+` = alter_graphic("rect", fill = col["Fusion+"])
  )

  column_title = paste0("Top20 Fusions For ", cancer_type1)
  output_file <- paste0("/data/data2/Projects/Fusion/paper_plot.old/oncoprint/", cancer_type1, ".png")
  # 打开 PNG 图形设备，设置分辨率为 600
  plot.new()
  png(output_file, width = 10, height = 8, units = "in", res = 600)
  oncoPrint(matrix_result,
            alter_fun = alter_fun, col = col,
            remove_empty_columns = TRUE,
            column_title = column_title)
  dev.off()
}
