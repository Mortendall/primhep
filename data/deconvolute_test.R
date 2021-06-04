library(tidyverse)
library(here)
library(openxlsx)
library(immunedeconv)

count_matrix <- openxlsx::read.xlsx(here::here("data_raw/count_matrix.xlsx"))

setup <- openxlsx::read.xlsx(here::here("data_raw/metadata.xlsx")) %>%
  dplyr::filter(!Sample == "558L") %>%
  dplyr::filter(Material == "Liver") %>%
  dplyr::arrange(desc(Genotype))
count_matrix <- count_matrix %>%
  dplyr::select(setup$Sample, Geneid)
count_matrix$Geneid<-stringr::str_remove_all(count_matrix$Geneid, "\\..*")

key <- clusterProfiler::bitr(count_matrix$Geneid, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
count_matrix <- count_matrix %>%
  dplyr::filter(Geneid %in% key$ENSEMBL)

count_matrix <- left_join(count_matrix, key, by = c("Geneid" = "ENSEMBL"))

matrix_key <- count_matrix %>%
  dplyr::select(Geneid, SYMBOL)

count_matrix <- count_matrix %>%
  dplyr::select(-Geneid, -SYMBOL)
count_matrix <- as.matrix(count_matrix)
rownames(count_matrix)<-toupper(matrix_key$SYMBOL)


set_cibersort_binary(here::here("CIBERSORT.R"))
set_cibersort_mat(here::here("LM22.txt"))

res_cibersort <- deconvolute(count_matrix, "cibersort")
res_xcell <- deconvolute(count_matrix, "xcell")

meta_heat_map <- setup %>%
  dplyr::select(Genotype)
rownames(meta_heat_map)<- setup$Sample

res_hm <- res_cibersort %>%
  dplyr::select(-cell_type)
rownames(res_hm)<-res_cibersort$cell_type
pheatmap::pheatmap(res_hm,
                   # treeheight_col = 0,
                   # treeheight_row = 0,
                   # legend = T,
                   # na_col = "white",
                   # Colv = NA,
                   # na.rm = T,
                    cluster_cols = F,
                   # fontsize_row = 8,
                   # fontsize_col = 11,
                   # cellwidth = 12,
                   # cellheight = 7,
                   annotation_col = meta_heat_map
)

#####statistical analysis of data#####

res_cibersort_long <- tidyr::pivot_longer(res_cibersort, cols = -"cell_type", names_to = "Sample", values_to = "percentage")
res_cibersort_long <- dplyr::left_join(res_cibersort_long, setup, by = "Sample") %>% dplyr::select(cell_type, Sample, percentage, Genotype)

sum_stats <- res_cibersort_long %>% group_by(Genotype, cell_type) %>% rstatix::get_summary_stats()

res_cibersort_long %>% group_by(Genotype, cell_type) %>% rstatix::identify_outliers(percentage) %>%  View()
