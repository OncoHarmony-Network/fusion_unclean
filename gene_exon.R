library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

#gene annotation
library(rtracklayer)
#metafusion
#meta_new <- fread("/home/data2/Projects/Fusion/fusiondb/metafusion_calls_cohort_update.tsv", sep = "\t") will lead lose C3N-03849-03 after's -03
cp_file_name <- "/home/data2/Projects/Fusion/unclean/gdc_split/fusion/metafusion_calls_matched.tsv"
meta_new <- fread(cp_file_name , sep = "\t")  
meta_new_with_commas <- meta_new |>
                        filter(grepl(",", exon1) | grepl(",", exon2)) |>
                        separate(col = exon1, into = c("exon1_part1", "exon1_part2"), sep = ",", remove = FALSE) |>
                        separate(col = chr1, into = c("chr11", "chr12"), sep = ",", remove = FALSE,fill = "right") |>
                        mutate(
                            chr12 = ifelse(is.na(chr12), chr11, chr12),
                            exon1_part1_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon1_part1)),  
                            exon1_part1_end = as.integer(sub(".*-(\\d+)", "\\1", exon1_part1)),      
                            exon1_part2_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon1_part2)),  
                            exon1_part2_end = as.integer(sub(".*-(\\d+)", "\\1", exon1_part2))      
                          ) |>
                        separate(col = exon2, into = c("exon2_part1", "exon2_part2"), sep = ",", remove = FALSE) |>
                        separate(col = chr2, into = c("chr21", "chr22"), sep = ",", remove = FALSE,fill = "right") |>
                        mutate(
                            chr22 = ifelse(is.na(chr22), chr21, chr22),
                            exon2_part1_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon2_part1)),  
                            exon2_part1_end = as.integer(sub(".*-(\\d+)", "\\1", exon2_part1)),      
                            exon2_part2_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon2_part2)),  
                            exon2_part2_end = as.integer(sub(".*-(\\d+)", "\\1", exon2_part2))      
                          ) 
meta_new_without_commas <- meta_new |>
                        filter(!(grepl(",", exon1) | grepl(",", exon2))) |>
                        mutate(
                               exon1_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon1)),  
                               exon1_end = as.integer(sub(".*-(\\d+)", "\\1", exon1)),      
                               exon2_start = as.integer(sub(".*:(\\d+)-.*", "\\1", exon2)),  
                               exon2_end = as.integer(sub(".*-(\\d+)", "\\1", exon2))      
                              )
     
#Transcript
Cof <- fread("/home/data2/Projects/Fusion/fusiondb/TranscriptionFactor/Homo_sapiens_Cof.txt", sep = "\t")
TF <- fread("/home/data2/Projects/Fusion/fusiondb/TranscriptionFactor/Homo_sapiens_TF.txt", sep = "\t")

refer_Homo <- fread("/home/data2/Projects/Fusion/MetaFusion/reference_files/Homo_sapiens.gene_info")
refer1 <- refer_Homo |> 
          select(type_of_gene,Symbol) |> 
          dplyr::distinct()
meta_ref_gene_name_to_type <- setNames(refer1$type_of_gene, refer1$Symbol) 
refer2 <- refer_Homo |> 
          select(Synonyms,type_of_gene) |> 
          separate_rows(Synonyms, sep = "\\|") |>
          filter(Synonyms != "-")|> 
          dplyr::distinct()
meta_ref_gene_name_to_type2 <- setNames(refer2$type_of_gene, refer2$Synonyms)

#use refer_bed to annote gene
refer_bed <- fread("/home/data2/Projects/Fusion/MetaFusion/reference_files/hg38_gene.bed") |>
             mutate(V2 = as.integer(V2 + 1) ) |>
             select(c(V1,V2,V3,V8)) 

colnames(refer_bed) <- c("chr1", "exon1_start", "exon1_end", "gene_name")
result1 <- refer_bed |>
                dplyr::group_by(chr1, exon1_start, exon1_end) |>
                dplyr::distinct() |>
                dplyr::mutate(
                  detected_name1 = ifelse(n() > 1, paste(gene_name, collapse = ","), gene_name),
                  chr1 = as.character(chr1)
                )  |>
                dplyr::ungroup() |>
                dplyr::select(-gene_name) |>
                dplyr::distinct()
# found gene detected region gene
#without_commas
colnames(result1) <- c("chr1", "exon1_start", "exon1_end", "gene1_name")
meta_new_without_commas1 <- meta_new_without_commas |>
                           left_join(result1 , by = c("chr1", "exon1_start", "exon1_end")) 


colnames(result1) <- c("chr2","exon2_start", "exon2_end","gene2_name")
meta_new_without_commas1  <- meta_new_without_commas1 |>
                            left_join(result1 , by = c("chr2", "exon2_start", "exon2_end")) 

meta_new_without_commas1  <- meta_new_without_commas1 |>
                             select(-c("exon1_start", "exon1_end",
                                       "exon2_start", "exon2_end"))

#with_commas
colnames(result1) <- c("chr11","exon1_part1_start", "exon1_part1_end","detected_name1_part1")
meta_new_with_commas1  <- meta_new_with_commas |>
                         left_join(result1 , by = c("chr11","exon1_part1_start", "exon1_part1_end"))

colnames(result1) <- c("chr12","exon1_part2_start", "exon1_part2_end","detected_name1_part2")
meta_new_with_commas1  <- meta_new_with_commas1 |>
                         left_join(result1 , by = c("chr12","exon1_part2_start", "exon1_part2_end"))  

colnames(result1) <- c("chr21", "exon2_part1_start", "exon2_part1_end","detected_name2_part1")
meta_new_with_commas1  <- meta_new_with_commas1 |>
                         left_join(result1 , by = c("chr21", "exon2_part1_start", "exon2_part1_end"))

colnames(result1) <- c("chr22","exon2_part2_start", "exon2_part2_end","detected_name2_part2")
meta_new_with_commas1  <- meta_new_with_commas1 |>
                          left_join(result1 , by = c("chr22", "exon2_part2_start", "exon2_part2_end")) 

meta_new_with_commas1  <- meta_new_with_commas1 |>
                          select(-c("chr21", "exon2_part1_start", "exon2_part1_end",
                                    "chr22", "exon2_part2_start", "exon2_part2_end",
                                    "chr12","exon1_part2_start", "exon1_part2_end",
                                    "exon1_part1","exon1_part2","exon2_part1","exon2_part2",
                                    "chr11","exon1_part1_start", "exon1_part1_end")) 

meta_new_with_commas1 <- meta_new_with_commas1 |>
                                      mutate(
                                        gene1_name = ifelse(detected_name1_part1 == detected_name1_part2, detected_name1_part1, paste(detected_name1_part1, detected_name1_part2, sep = ",")),
                                        gene2_name = ifelse(detected_name2_part1 == detected_name2_part2, detected_name2_part1, paste(detected_name2_part1, detected_name2_part2, sep = ","))
                                      )   
meta_new_with_commas1 <- meta_new_with_commas1 |>
                         select(-c(detected_name1_part1, detected_name1_part2,detected_name2_part1, detected_name2_part2 )) 
 
#gene1

meta_new1 <- meta_new1 |>
  mutate(
    gene1_type = case_when(
      detected_name1 %in% names(gene_name_to_type) ~ gene_name_to_type[detected_name1],
      TRUE ~ NA_character_
    )
  ) |>
    mutate(
      gene1_type = case_when(
        is.na(gene1_type) & detected_name1 %in% names(meta_ref_gene_name_to_type) ~ meta_ref_gene_name_to_type[detected_name1],
        TRUE ~ gene1_type  
      )
    ) 
#gene2

meta_new1 <- meta_new1 |>
  mutate(
    gene2_type = case_when(
      detected_name2 %in% names(gene_name_to_type) ~ gene_name_to_type[detected_name2],
      TRUE ~ NA_character_
    )
  ) |>
    mutate(
      gene2_type = case_when(
        is.na(gene2_type) & detected_name2 %in% names(meta_ref_gene_name_to_type) ~ meta_ref_gene_name_to_type[detected_name2],
        TRUE ~ gene2_type  
      )
    )


    


#use refer_bed to annote gene



