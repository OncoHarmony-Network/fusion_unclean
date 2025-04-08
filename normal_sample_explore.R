library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

file_meta <- "/home/data2/Projects/Fusion/fusiondb_v3/metafusion_gdc_divide.tsv"
metafusion <- fread(file_meta , sep = "\t")

#Not in Red Herrings
#metafusion
#TODO:检查适合注释的行就是少，数据库来源那个？肿瘤注释和正常注释？我先找一下俩reads标准。肿瘤分析多少没被注释
#Overlap假如分析之后，如果reads少的甩了，注释就多了，就star标明一下fusion检出可信标准，卡阈值。
#注释源来自？

Herring_Supplement <- metafusion |>
                      filter(sample_type == "Normal" & blocklist == "FALSE") #blocklist mean false positives
#> nrow(Herring_Supplement)
#[1] 109968
Herring_Supplement2 <- Herring_Supplement |>
  filter(grepl("GTEx_recurrent_StarF2019", cancer_db_hits) |
  grepl("DGD_PARALOGS", cancer_db_hits) |
  grepl("HGNC_GENEFAM", cancer_db_hits) |
  grepl("BodyMap", cancer_db_hits) |
  grepl("Greger_Normal", cancer_db_hits) |
  grepl("ConjoinG", cancer_db_hits) |
  grepl("Babiceanu_Normal", cancer_db_hits))
#nrow(Herring_Supplement2)
#[1] 81

Herring_Supplement1 <- Herring_Supplement |>
                      filter(cancer_db_hits == "") 
#nrow(Herring_Supplement1)
#[1] 102449
Herring_Supplement2 <- Herring_Supplement |>
  filter(grepl("GTEx_recurrent_StarF2019", cancer_db_hits) |
  grepl("DGD_PARALOGS", cancer_db_hits) |
  grepl("HGNC_GENEFAM", cancer_db_hits) |
  grepl("BodyMap", cancer_db_hits) |
  grepl("Greger_Normal", cancer_db_hits) |
  grepl("ConjoinG", cancer_db_hits) |
  grepl("Babiceanu_Normal", cancer_db_hits))
#nrow(Herring_Supplement2)
#[1] 81
Herring_Supplement3 <-  Herring_Supplement |>
  filter(!grepl("GTEx_recurrent_StarF2019", cancer_db_hits) &
  !grepl("DGD_PARALOGS", cancer_db_hits) &
  !grepl("HGNC_GENEFAM", cancer_db_hits) &
  !grepl("BodyMap", cancer_db_hits) &
  !grepl("Greger_Normal", cancer_db_hits) &
  !grepl("ConjoinG", cancer_db_hits) &
  !grepl("Babiceanu_Normal", cancer_db_hits) & 
  cancer_db_hits != "")
#nrow(Herring_Supplement3)
#[1] 7438
#  meta_new |>
#  filter(sample_type == "Normal" & blocklist == "FALSE" & cancer_db_hits == "") |>
#  nrow()
#[1] 102449
#  meta_new |>
#  filter(sample_type == "Normal" & blocklist == "FALSE") |>
#  nrow()
#[1] 109968
#TODO 癌旁和normal
Herring_file_name <- "/home/data2/Projects/Fusion/fusiondb/Herring_Supplement.tsv"
write.table(Herring_Supplement, 
            file = Herring_file_name, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


            Herring_Supplement_type_high <- metafusion  |>
              +                                 filter(sample_type == "Normal" & blocklist == "FALSE" & tools == "arriba,starfusion")


Tumor_type_high <- metafusion  |>
                   filter(sample_type == "Tumor" & blocklist == "FALSE" & tools == "arriba,starfusion")
nrow(Tumor_type_high)
### [1] 95097
Tumor_type_high2 <- metafusion  |>
                   filter(sample_type == "Tumor" & blocklist == "FALSE" & tools == "arriba,starfusion"& cancer_db_hits == "")
nrow(Tumor_type_high2)    

## Normal Data
### Firstly, Basic data
Normal_data <- metafusion  |>
               filter(sample_type == "Normal")

### Sample num
Normal_data |> select(Sample_ID) |> unique() |> nrow()
### [1] 1781

### Fusion pair num
nrow(Normal_data)
### [1] 114007

### Fusion pair unique num
Normal_data |> select(c(gene1,gene2)) |> unique() |> nrow()
### [1] 25339

### Reads Tier
### split_cnt≥1 & span_cnt≥1
Normal_data_1 <- Normal_data  |>
                 filter(max_split_cnt >= 1 & max_span_cnt>= 1)
### split_cnt≥3 & span_cnt≥3
Normal_data_3 <- Normal_data  |>
                 filter(max_split_cnt >= 3 & max_span_cnt >= 3)
### split_cnt≥5 & span_cnt≥5
Normal_data_5 <- Normal_data  |>
                 filter(max_split_cnt >= 5 & max_span_cnt >= 5)
### split_cnt≥10 & span_cnt≥10
Normal_data_10 <- Normal_data  |>
                 filter(max_split_cnt >= 10 & max_span_cnt >= 10)

### Tools Tier
Normal_data_arriba <- Normal_data  |>
                      filter(tools == "arriba")
nrow(Normal_data_arriba)
### [1] 110685

Normal_data_starfusion <- Normal_data  |>
                      filter(tools == "starfusion")
nrow(Normal_data_starfusion)
### [1] 259

Normal_data_two <- Normal_data  |>
                      filter(tools == "arriba,starfusion")
nrow(Normal_data_two)
### [1] 3063

### FusionAnnotator Tier
normal_markers <- c("GTEx_recurrent_StarF2019", "DGD_PARALOGS", "HGNC_GENEFAM", "BodyMap", "Greger_Normal", "ConjoinG", "Babiceanu_Normal")
tumor_markers <- c( "Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic", "YOSHIHARA_TCGA", "ChimerKB", "ChimerPub", "ChimerSeq", "Klijn_CellLines", "Larsson_TCGA", "CCLE_StarF2019", "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018", "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023")

### Past data not contain TCGA, TARGET, CPTAC(Except GTEx-brain) normal sample data  
###            BodyMap	40
###            Babiceanu_Normal	239
###            ConjoinG	698
###            DGD_PARALOGS	3543
###            Greger_Normal	113
###            GTEx_recurrent_StarF2019	4047
###            HGNC_GENEFAM	3543

Normal_data_annormal0 <- Normal_data |>
                         filter(rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) == 0)
nrow(Normal_data_annormal0)
### [1] 113926

Normal_data_annormal1 <- Normal_data |>
                         filter(rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) >= 1)
nrow(Normal_data_annormal1 )
### [1] 81

Normal_data_annormal2 <- Normal_data |>
                          filter(rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) >= 2)
nrow(Normal_data_annormal2)
### [1] 20

Normal_data_annormal3 <- Normal_data |>
                         filter(rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) >= 3)
nrow(Normal_data_annormal3)
### [1] 0

## Tumor Data
Tumor_data <- metafusion  |>
              filter(sample_type == "Tumor")