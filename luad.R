
library(dplyr)
library(data.table)
library(stringr)
cff_format <-
  c(
    "chr1",
    "pos1",
    "strand1",
    "chr2",
    "pos2",
    "strand2",
    "library",
    "sample_name",
    "sample_type",
    "disease",
    "tool",
    'split_cnt',
    "span_cnt",
    "t_gene1",
    "t_area1",
    "t_gene2",
    "t_area2"
  )

opt <- commandArgs(TRUE)
ensure_dir_exists <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

make_arriba <- function(sample_file) {
  df <- as.data.frame(matrix(ncol = 0, nrow = nrow(sample_file)))
  
  df$t_gene1 <- sample_file$`#gene1`
  df$t_gene2 <- sample_file$gene2
  df$chr1 <- str_split_fixed(sample_file$breakpoint1, ":", 2)[, 1]
  df$pos1 <- str_split_fixed(sample_file$breakpoint1, ":", 2)[, 2]
  df$chr2 <-  str_split_fixed(sample_file$breakpoint2, ":", 2)[, 1]
  df$pos2 <- str_split_fixed(sample_file$breakpoint2, ":", 2)[, 2]
  df$strand1 <-
    str_split_fixed(sample_file$`strand1(gene/fusion)`, "/", 2)[, 1]
  df$strand2 <-
    str_split_fixed(sample_file$`strand2(gene/fusion)`, "/", 2)[, 1]
  df$tool <- "arriba"
  df$split_cnt <-
    ifelse(!is.na(sample_file$split_reads1),
           sample_file$split_reads1,
           -1)
  df$span_cnt <-
    ifelse(!is.na(sample_file$discordant_mates),
           sample_file$discordant_mates,
           -1)
  df$t_area1 <- sample_file$site1
  df$t_area2 <- sample_file$site2
  df$sample_name <- sample_file$sample
  
  return(df)
}


make_starfusion <- function(sample_file) {
  df <- as.data.frame(matrix(ncol = 0, nrow = nrow(sample_file)))
  df$t_gene1 <- str_split_fixed(sample_file$LeftGene, "\\^", 2)[, 1]
  df$t_gene2 <- str_split_fixed(sample_file$RightGene, "\\^", 2)[, 1]
  df$chr1 <-
    str_replace(str_split_fixed(sample_file$LeftBreakpoint, ":", 3)[, 1],
                "chr",
                "")
  df$pos1 <- str_split_fixed(sample_file$LeftBreakpoint, ":", 3)[, 2]
  df$chr2 <-
    str_replace(str_split_fixed(sample_file$RightBreakpoint, ":", 3)[, 1],
                "chr",
                "")
  df$pos2 <-  str_split_fixed(sample_file$RightBreakpoint, ":", 3)[, 2]
  df$strand1 <-
    str_split_fixed(sample_file$LeftBreakpoint, ":", 3)[, 3]
  df$strand2 <-
    str_split_fixed(sample_file$RightBreakpoint, ":", 3)[, 3]
  df$tool <- "starfusion"
  df$split_cnt <-
    ifelse(!is.na(sample_file$JunctionReadCount),
           sample_file$JunctionReadCount,
           -1)
  df$span_cnt <-
    ifelse(!is.na(sample_file$SpanningFragCount),
           sample_file$SpanningFragCount,
           -1)
  df$t_area1 <- sample_file$SpliceType
  df$t_area2 <- sample_file$SpliceType
  df$sample_name <- sample_file$sample
  return(df)
  
}
file_path1 <- "/data/data2/Projects/Fusion/GeneFusionPipelines/data_cff"
dir <- "OAK"
dir_path_arriba <- file.path(file_path1, dir, "arriba1")
arriba_sample <- list.dirs(path = dir_path_arriba, full.names = FALSE, recursive = FALSE)
cleaned_arriba_sample <- gsub("\\.output", "", arriba_sample)

arriba_sample1 <- list.dirs(path = dir_path_arriba, full.names = TRUE, recursive = FALSE)

results_list1 <- lapply(seq_along(arriba_sample1), function(i) {
  dir_path2 <- file.path(arriba_sample1[[i]], "fusions.tsv")
  
  if (file.exists(dir_path2)) {
    
    if (file.info(dir_path2)$size > 0) {
      data <- fread(dir_path2) |>
              filter(confidence == "high") 
      data[, sample := cleaned_arriba_sample[[i]]]
      return(data)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
  return(data)
})
dir1 <- rbindlist(results_list1, use.names = TRUE, fill = TRUE)
tool_cff1 <- make_arriba(dir1)

if (nrow(tool_cff1) > 0) {
  tool_cff1$library <- "RNA"
  tool_cff1$sample_type <- "Tumor"
  tool_cff1$disease <- NA
}
tool_cff1$strand1[tool_cff1$strand1 == "."] <- NA
tool_cff1$strand2[tool_cff1$strand2 == "."] <- NA
tool_cff1 <- tool_cff1[, cff_format]

dir_path_starfusion <- file.path(file_path1, dir, "star_fusion_outdir1")

starfusion_sample <- list.dirs(path = dir_path_starfusion, full.names = FALSE, recursive = FALSE)

starfusion_sample1 <- list.dirs(path = dir_path_starfusion, full.names = TRUE, recursive = FALSE)

results_list2 <- lapply(seq_along(starfusion_sample1), function(i) {
  dir_path2 <- file.path(starfusion_sample1[[i]], "star-fusion.fusion_predictions.tsv")
  if (file.exists(dir_path2)) {
    
    if (file.info(dir_path2)$size > 0) {
      data <- fread(dir_path2) |>
              filter(FFPM >= 0.5 & LargeAnchorSupport == "YES_LDAS")
      data[, sample := starfusion_sample[[i]]]
      return(data)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
  return(data)
})
dir2 <- rbindlist(results_list2, use.names = TRUE, fill = TRUE)
tool_cff2 <- make_starfusion(dir2)

if (nrow(tool_cff2) > 0) {
  tool_cff2$library <- "RNA"
  tool_cff2$sample_type <- "Tumor"
  tool_cff2$disease <- NA
}
tool_cff2$strand1[tool_cff2$strand1 == "."] <- NA
tool_cff2$strand2[tool_cff2$strand2 == "."] <- NA
tool_cff2 <- tool_cff2[, cff_format]
tool_cff <- rbind(tool_cff1,tool_cff2)
tool_cff$chr1 <- gsub("chr", "", tool_cff$chr1)
tool_cff$chr2 <- gsub("chr", "", tool_cff$chr2)
tool_cff$t_gene1 = convert_unmatched_symbol(tool_cff$t_gene1, data_sub)
tool_cff$t_gene2 = convert_unmatched_symbol(tool_cff$t_gene2, data_sub)
