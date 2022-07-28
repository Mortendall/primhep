edgeR_data <- list(Liver = NA,
                   Cell_susp = NA,
                   Prim_hep = NA,
                   Liv_vs_CS_WT = NA,
                   Liv_vs_PH_WT = NA,
                   CS_vs_PH_WT = NA,
                   Liv_vs_CS_KO = NA,
                   Liv_vs_PH_KO = NA,
                   CS_vs_PH_KO = NA)
  for (i in 1:9){
    edgeR_data[[i]]<- openxlsx::read.xlsx(here("data/edgeR.xlsx"),sheet = i)
  }

metadata <- load_metadata("metadata.xlsx")


#Quality_control_plots(counts, metadata)
#QC reveals that 558L looks very odd. Remove and rerun. QCplots saved as QCplots_before_filtering
metadata <- metadata %>%
  dplyr::filter(!Sample == "558L")

design <- Generate_design_matrix(metadata)

ctrsts <- limma::makeContrasts(
  Liver = HNKO_L - WT_L,
  Cell_susp = HNKO_CS - WT_CS,
  Prim_hep = HNKO_PH - WT_PH,
  Liv_vs_CS_WT = WT_L - WT_CS,
  Liv_vs_PH_WT = WT_L - WT_PH,
  CS_vs_PH_WT = WT_CS - WT_PH,
  Liv_vs_CS_KO = HNKO_L - HNKO_CS,
  Liv_vs_PH_KO = HNKO_L - HNKO_PH,
  CS_vs_PH_KO = HNKO_CS - HNKO_PH,
  levels = design)



cpm_matrix <- readRDS(here("data/cpm_matrix.rds"))
####camera with reactome####

#rLst <- fread("https://reactome.org/download/current/Ensembl2Reactome.txt", header = FALSE)
rLst <- openxlsx::read.xlsx("C:/Users/tvb217/Documents/R/tmp/Reactomedatabase.xlsx")
rLst <- data.table::as.data.table(rLst, keep.rownames = T)
rLst <- rLst %>%
  dplyr::filter(V6 == "Mus musculus")
reactomeName <- data.table::data.table(ID = unique(rLst$V2), TERM = unique(rLst$V4))
reactomeList <- tapply(rLst$V1, rLst$V2, list)
reactomeList <- Filter(. %>% length %>% is_greater_than(4), reactomeList) # Remove small categories
reactomeList <- Filter(. %>% length %>% is_less_than(501), reactomeList) # Remove small categories

reactome_data <- camera_reactome(reactomeList, reactomeName, count_matrix = cpm_matrix, design_matrix = design, contrast_matrix = ctrsts)

####camera with GO####
org.db <- "org.Mm.eg.db"

keysGO <- AnnotationDbi::keys(GO.db)
termGO <- AnnotationDbi::select(GO.db, keys=keysGO, columns=c("TERM", "ONTOLOGY")) %>% data.table
termGO <- termGO %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  dplyr::select(-ONTOLOGY)
setnames(termGO, "GOID", "ID")

cyt.go.genes <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
cyt.go.genes<- cyt.go.genes[names(cyt.go.genes) %in% termGO$ID]
cyt.go.genes <- Filter(. %>% length %>% is_greater_than(4), cyt.go.genes) # Remove small categories
cyt.go.genes <- Filter(. %>% length %>% is_less_than(501), cyt.go.genes) # Remove large categories

go_data <- camera_go(org.db, cpm_matrix = cpm_matrix, design_matrix = design, contrast_matrix = ctrsts)

#write.xlsx(camera_test, here::here("data/Reactome_data.xlsx"), asTable = TRUE)
#write.xlsx(GO_test, here::here("data/GO_data.xlsx"), asTable = TRUE)


####GO on signficant genes####
go_sig_genes <- goAnalysis(edgeR_data, "MF", "")
go_sig_genes_cc <- goAnalysis(edgeR_data, "CC", "")
#printGOterms(go_data)

for (i in 1:length(go_sig_genes)){
  go_sig_genes[[i]] <- clusterProfiler::setReadable(go_sig_genes[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
}
#saveRDS(go_sig_genes, file = here::here("data/go_sig.data.rds"))

#####annotate GO data#####
cameraGoAnnotated <- annotateWithGenes(go_data, cyt.go.genes, fromType = 'ENTREZID')
cameraReactomeAnnotated <- annotateWithGenes(reactome_data, reactomeList, fromType = "ENSEMBL")
#saveRDS(cameraReactomeAnnotated, here("data/annotatedReactome.rds"))
#saveRDS(cameraGoAnnotated, here("data/annotatedGO.rds"))

#####Extraction of significant genes from GO list and figures for paper#####
#Note: test is actually BP data
test <- readRDS(here("data/go_sig.data.rds"))
Liver <- dotplot(test[[1]], font.size = 14)+
  ggtitle("Liver - Genotype effect")

tiff("LiverGO.tif", units = "cm", width = 25, height = 15, res = 300)
Liver
dev.off()

PH <- dotplot(test[[3]], font.size = 14)+
  ggtitle("Primary Hepatocyte - Genotype effect")

tiff("PHGO.tif", units = "cm", width = 25, height = 15, res = 300)
PH
dev.off()

dotplot(test[[5]], font.size = 14)+
  ggtitle("Liver vs Prim hep WT")
dotplot(test[[8]], font.size = 14)+
  ggtitle("Liver vs Prim hep HNKO")


mito_genes <- test[[1]]@result
View(mito_genes)

gene_list <- mito_genes %>%
  filter(Description == "cellular respiration"|Description =="NADH dehydrogenase complex assembly") %>%
  dplyr::select(geneID)

gene_list_unsplit <- c(unlist(str_split(gene_list[1,], "/")),unlist(str_split(gene_list[2,], "/")))

gene_list_unique <- unique(gene_list_unsplit)
cpm_key <-   clusterProfiler::bitr(
  gene_list_unique,
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db"
)


cpm_annotated <- as.data.frame(cpm_matrix)
cpm_annotated <- cpm_annotated %>%
  dplyr::filter(rownames(cpm_matrix) %in% cpm_key$ENSEMBL)

columns <- colnames(cpm_annotated)
column_order <- c("544L", "548L", "554L", "556L", "564L", "540L", "542L", "546L", "552L", "562L", "566L", "544CS", "548CS", "554CS", "556CS", "558CS", "564CS", "540CS", "542CS", "546CS", "552CS", "562CS", "566CS", "544PH", "548PH", "554PH", "556PH", "558PH", "564PH", "540PH", "542PH", "546PH", "552PH", "562PH",  "566PH")

cpm_test <- cpm_annotated %>%
  dplyr::select(column_order)

conv <- clusterProfiler::bitr(rownames(cpm_test),
                                            fromType = "ENSEMBL",
                                            toType = "SYMBOL",
                                            OrgDb = "org.Mm.eg.db")
rownames(cpm_test) <- conv$SYMBOL
#generate and organize metadata for heatmap
meta_heat_map <- metadata %>%
 dplyr::arrange(match(Sample, column_order)) %>%
  dplyr::filter(!Sample == "558L") %>%
  dplyr::select(Sample, Group)
rownames(meta_heat_map)<-meta_heat_map$Sample
meta_heat_map <- meta_heat_map %>%
  dplyr::select(-Sample)
meta_heat_map <- meta_heat_map %>%
  dplyr::mutate(Group = case_when(
    Group == "WT_L"~"WT L",
    Group =="WT_CS"~"WT CS",
    Group == "WT_PH"~"WT PH",
    Group =="HNKO_L"~"HNKO L",
    Group =="HNKO_CS"~"HNKO CS",
    Group =="HNKO_PH"~"HNKO PH"
  ))

meta_heat_map$Group <- factor(meta_heat_map$Group, levels = c("WT L", "HNKO L", "WT CS", "HNKO CS", "WT PH", "HNKO PH"))

OxPhos <- pheatmap(cpm_test,
         treeheight_col = 0,
         treeheight_row = 0,
         scale = "row",
         legend = T,
         na_col = "white",
         Colv = NA,
         na.rm = T,
         cluster_cols = F,
         fontsize_row = 8,
         fontsize_col = 11,
         cellwidth = 12,
         cellheight = 7,
         annotation_col = meta_heat_map
)

# tiff("Oxphos.tif", units = "cm", width = 25, height = 20, res = 300)
# OxPhos
# dev.off()

#####Repeat analysis with MF####
go_sig_genes_MF <- goAnalysis(edgeR_data, "MF","")
#printGOterms(go_sig_genes_MF)

for (i in 1:length(go_sig_genes_MF)){
  go_sig_genes_MF[[i]] <- clusterProfiler::setReadable(go_sig_genes_MF[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
}
#saveRDS(go_sig_genes_MF, file = here::here("data/go_sig_MF.data.rds"))
test <- readRDS(here::here("data/go_sig_MF.data.rds"))
Liver <- clusterProfiler::dotplot(test[[1]], font.size = 14)+
  ggtitle("Liver - Genotype effect")
#
 tiff("LiverGO_MF.tif", units = "cm", width = 30, height = 20, res = 300)
  Liver
  dev.off()

PH <- dotplot(test[[3]], font.size = 14)+
  ggtitle("Primary Hepatocyte - Genotype effect")

  tiff("PHGO_MF.tif", units = "cm", width = 30, height = 20, res = 300)
  PH
  dev.off()

dotplot(test[[5]], font.size = 14)+
  ggtitle("Liver vs Prim hep WT")
dotplot(test[[8]], font.size = 14)+
  ggtitle("Liver vs Prim hep HNKO")

#####Repeat analysis with CC####
go_sig_genes_CC <- goAnalysis(edgeR_data, "CC","")
#printGOterms(go_sig_genes_CC,"CC_go")

for (i in 1:length(go_sig_genes_CC)){
  go_sig_genes_CC[[i]] <- clusterProfiler::setReadable(go_sig_genes_CC[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
}
saveRDS(go_sig_genes_CC, file = here::here("data/go_sig_CC.data.rds"))
test <- readRDS(here::here("data/go_sig_CC.data.rds"))
Liver <- clusterProfiler::dotplot(test[[1]], font.size = 14)+
  ggtitle("Liver - Genotype effect")
#
 tiff("LiverGO_CC.tif", units = "cm", width = 25, height = 15, res = 300)
  Liver
  dev.off()

PH <- dotplot(test[[3]], font.size = 14)+
  ggtitle("Primary Hepatocyte - Genotype effect")

 tiff("PHGO_CC.tif", units = "cm", width = 25, height = 15, res = 300)
 PH
 dev.off()

dotplot(test[[5]], font.size = 14)+
  ggtitle("Liver vs Prim hep WT")
dotplot(test[[8]], font.size = 14)+
  ggtitle("Liver vs Prim hep HNKO")

#####Investigate GO analysis for WT L vs PH for up- and down-regulated genes####
pos_genes<- edgeR_data[[5]] %>%
  dplyr::filter(logFC >0) |>
  dplyr::filter(FDR<0.05)
bg <- edgeR_data[[1]]
bg_list <- clusterProfiler::bitr(
  bg$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Mm.eg.db",
  drop = T
)


  eg <- clusterProfiler::bitr(
    pos_genes$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
  )
  goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                         universe = bg_list$ENTREZID,
                                         OrgDb = org.Mm.eg.db,
                                         ont = "MF")
  neg_genes<- edgeR_data[[5]] %>%
    dplyr::filter(logFC <0)|>
    dplyr::filter(FDR<0.05)
eg_neg <- clusterProfiler::bitr(
  neg_genes$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Mm.eg.db",
  drop = T
)
goResults_neg <- clusterProfiler::enrichGO(gene = eg_neg$ENTREZID,
                                       universe = bg_list$ENTREZID,
                                       OrgDb = org.Mm.eg.db,
                                       ont = "MF")

plot_neg <- enrichplot::dotplot(goResults_neg, font.size = 14)+
  ggtitle("Liver vs PH- Downregulated Genes")
plot_pos <- enrichplot::dotplot(goResults, font.size = 14)+
  ggtitle("Liver vs PH- Upregulated Genes")

tiff("L_vs_PH_MF.tif", units = "cm", width = 60, height = 20, res = 300)
plot_pos+plot_neg
dev.off()
#repeat for CC
goResults_CC <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                       universe = bg_list$ENTREZID,
                                       OrgDb = org.Mm.eg.db,
                                       ont = "CC")
goResults_neg_CC <- clusterProfiler::enrichGO(gene = eg_neg$ENTREZID,
                                           universe = bg_list$ENTREZID,
                                           OrgDb = org.Mm.eg.db,
                                           ont = "CC")
plot_neg_CC <- enrichplot::dotplot(goResults_neg_CC, font.size = 20)+
  ggtitle("Genes Upregulated in PH  \n compared to L")+
  ggplot2::theme(plot.title = ggplot2::element_text(size = 24))
plot_pos_CC <- enrichplot::dotplot(goResults_CC, font.size = 20)+
  ggtitle("Genes Upregulated in L \n compared to PH")+
  ggplot2::theme(plot.title = ggplot2::element_text(size = 24))

# tiff(here::here("L_vs_PH_CC.tif"), units = "cm", width = 40, height = 20, res = 300)
# plot_pos_CC+plot_neg_CC
# dev.off()

goResults_CC <- clusterProfiler::setReadable(goResults_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
goResults_neg_CC <- clusterProfiler::setReadable(goResults_neg_CC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
collected_data <- list("negative_cc"=goResults_neg_CC@result,
                       "positive_cc"=goResults_CC@result)
openxlsx::write.xlsx(collected_data, here::here("data/Liver_vs_Hep_CC_seperated.xlsx"))
collected_data_MF <- list("negative_mf"=goResults_neg@result,
                       "positive_mf"=goResults@result)
openxlsx::write.xlsx(collected_data_MF, here::here("data/Liver_vs_Hep_MF_seperated.xlsx"))


#####Generate heatmap for calcium signaling proteins#####
calcium_genes <- test[[3]]@result

gene_list <- calcium_genes %>%
  filter(Description == "calcium ion binding") %>%
  dplyr::select(geneID)

gene_list_unsplit <- c(unlist(str_split(gene_list[1,], "/")),unlist(str_split(gene_list[2,], "/")))

gene_list_unique <- unique(gene_list_unsplit)
cpm_key <-   clusterProfiler::bitr(
  gene_list_unique,
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db"
)
#find signficant genes from L group
calcium_genes_L <- calcium_genes <- test[[1]]@result

gene_list_L <- calcium_genes_L %>%
  filter(Description == "calcium ion binding") %>%
  dplyr::select(geneID)

gene_list_unsplit_L <- c(unlist(str_split(gene_list_L[1,], "/")),unlist(str_split(gene_list_L[2,], "/")))

gene_list_unique_L <- unique(gene_list_unsplit_L)
cpm_key_L <-   clusterProfiler::bitr(
  gene_list_unique_L,
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db"
)


cpm_annotated <- as.data.frame(cpm_matrix)
cpm_annotated <- cpm_annotated %>%
  dplyr::filter(rownames(cpm_matrix) %in% cpm_key$ENSEMBL | rownames(cpm_matrix) %in% cpm_key_L$ENSEMBL)

columns <- colnames(cpm_annotated)
column_order <- c("544L", "548L", "554L", "556L", "564L", "540L", "542L", "546L", "552L", "562L", "566L", "544CS", "548CS", "554CS", "556CS", "558CS", "564CS", "540CS", "542CS", "546CS", "552CS", "562CS", "566CS", "544PH", "548PH", "554PH", "556PH", "558PH", "564PH", "540PH", "542PH", "546PH", "552PH", "562PH",  "566PH")

cpm_test <- cpm_annotated %>%
  dplyr::select(column_order)

conv <- clusterProfiler::bitr(rownames(cpm_test),
                              fromType = "ENSEMBL",
                              toType = "SYMBOL",
                              OrgDb = "org.Mm.eg.db")
rownames(cpm_test) <- conv$SYMBOL
cpm_test <- cpm_test %>%
  dplyr::arrange(rownames(cpm_test))
#generate and organize metadata for heatmap
meta_heat_map <- metadata %>%
  dplyr::arrange(match(Sample, column_order)) %>%
  dplyr::filter(!Sample == "558L") %>%
  dplyr::select(Sample, Group)
rownames(meta_heat_map)<-meta_heat_map$Sample
meta_heat_map <- meta_heat_map %>%
  dplyr::select(-Sample)
meta_heat_map <- meta_heat_map %>%
  dplyr::mutate(Group = case_when(
    Group == "WT_L"~"WT L",
    Group =="WT_CS"~"WT CS",
    Group == "WT_PH"~"WT PH",
    Group =="HNKO_L"~"HNKO L",
    Group =="HNKO_CS"~"HNKO CS",
    Group =="HNKO_PH"~"HNKO PH"
  ))

meta_heat_map$Group <- factor(meta_heat_map$Group, levels = c("WT L", "HNKO L", "WT CS", "HNKO CS", "WT PH", "HNKO PH"))
#mark by star in name if things were signficant in L

rownames(cpm_test) <-case_when(
   (rownames(cpm_test) %in% cpm_key_L$SYMBOL) & (rownames(cpm_test) %in% cpm_key$SYMBOL) ~ paste(rownames(cpm_test), "#*", sep = ""),
    (rownames(cpm_test) %in% cpm_key_L$SYMBOL) & (!rownames(cpm_test) %in% cpm_key$SYMBOL) ~ paste(rownames(cpm_test), "*", sep = ""),
  TRUE ~ rownames(cpm_test))

CalciumHm <- pheatmap(cpm_test,
                   treeheight_col = 0,
                   treeheight_row = 0,
                   scale = "row",
                   legend = T,
                   na_col = "white",
                   Colv = NA,
                   na.rm = T,
                   cluster_cols = F,
                   cluster_rows = F,
                   fontsize_row = 8,
                   fontsize_col = 11,
                   cellwidth = 12,
                   cellheight = 7,
                   annotation_col = meta_heat_map,
                   main = "Calcium Ion Binding"
)

# tiff("Calcium.tif", units = "cm", width = 25, height = 30, res = 300)
# CalciumHm
# dev.off()

#####Create OxRed heatmap for PH####
OxRed_genes <- test[[3]]@result

gene_list <- OxRed_genes %>%
  filter(Description == "steroid hydroxylase activity") %>%
  dplyr::select(geneID)

gene_list_unsplit <- c(unlist(str_split(gene_list[1,], "/")),unlist(str_split(gene_list[2,], "/")))

gene_list_unique <- unique(gene_list_unsplit)
cpm_key <-   clusterProfiler::bitr(
  gene_list_unique,
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db"
)


cpm_annotated <- as.data.frame(cpm_matrix)
cpm_annotated <- cpm_annotated %>%
  dplyr::filter(rownames(cpm_matrix) %in% cpm_key$ENSEMBL)

columns <- colnames(cpm_annotated)
column_order <- c("544L", "548L", "554L", "556L", "564L", "540L", "542L", "546L", "552L", "562L", "566L", "544CS", "548CS", "554CS", "556CS", "558CS", "564CS", "540CS", "542CS", "546CS", "552CS", "562CS", "566CS", "544PH", "548PH", "554PH", "556PH", "558PH", "564PH", "540PH", "542PH", "546PH", "552PH", "562PH",  "566PH")

cpm_test <- cpm_annotated %>%
  dplyr::select(column_order)

conv <- clusterProfiler::bitr(rownames(cpm_test),
                              fromType = "ENSEMBL",
                              toType = "SYMBOL",
                              OrgDb = "org.Mm.eg.db")
rownames(cpm_test) <- conv$SYMBOL
cpm_test <- cpm_test %>%
  dplyr::arrange(rownames(cpm_test))
#generate and organize metadata for heatmap
meta_heat_map <- metadata %>%
  dplyr::arrange(match(Sample, column_order)) %>%
  dplyr::filter(!Sample == "558L") %>%
  dplyr::select(Sample, Group)
rownames(meta_heat_map)<-meta_heat_map$Sample
meta_heat_map <- meta_heat_map %>%
  dplyr::select(-Sample)
meta_heat_map <- meta_heat_map %>%
  dplyr::mutate(Group = case_when(
    Group == "WT_L"~"WT L",
    Group =="WT_CS"~"WT CS",
    Group == "WT_PH"~"WT PH",
    Group =="HNKO_L"~"HNKO L",
    Group =="HNKO_CS"~"HNKO CS",
    Group =="HNKO_PH"~"HNKO PH"
  ))

meta_heat_map$Group <- factor(meta_heat_map$Group, levels = c("WT L", "HNKO L", "WT CS", "HNKO CS", "WT PH", "HNKO PH"))


OxRedHm <- pheatmap(cpm_test,
                      treeheight_col = 0,
                      treeheight_row = 0,
                      scale = "row",
                      legend = T,
                      na_col = "white",
                      Colv = NA,
                      na.rm = T,
                      cluster_cols = F,
                      cluster_rows = F,
                      fontsize_row = 8,
                      fontsize_col = 11,
                      cellwidth = 12,
                      cellheight = 7,
                      annotation_col = meta_heat_map,
                      main = "steroid hydroxylase activity"
)

 tiff("SteroidHydroxy.tif", units = "cm", width = 25, height = 30, res = 300)
 OxRedHm
 dev.off()

 #####Repeat figure generation for CC#####

 #####Test ReactomePA####

#####Upsetplot####
 sig_genes_names <- names(edgeR_data[4:6])
 sig_genes <-vector(mode = "list", length = 3)
 names(sig_genes)<-names(edgeR_data[4:6])
 for (i in 1:3){
   sig_genes[[i]]<- edgeR_data[[i+3]]
   sig_genes[[i]]<- sig_genes[[i]] %>%
     dplyr::filter(FDR < 0.05)
   sig_genes[[i]]<-sig_genes[[i]]$SYMBOL
 }


 names(sig_genes) <- c("Liver vs CS","Liver vs PH", "CS vs PH")

#  upsetPlot <- UpSetR::upset(UpSetR::fromList(sig_genes),
#                             # sets = order_upset,
#                             order.by = "freq",
#                             keep.order = T,
#                             text.scale = c(3, 2, 2, 1.2, 3, 3),
#                             point.size = 4
#  )
# upsetPlot

#make same upsetplot with ComplexUpset
order_upset <- c("Liver vs CS","Liver vs PH", "CS vs PH")
upset_data <- data.frame("Genes"= sig_genes[[1]],
                         "Group"= names(sig_genes[1]))

for (i in 2:length(sig_genes)){
  upset_data <- dplyr::add_row(upset_data,"Genes"=sig_genes[[i]],
                               "Group"=names(sig_genes)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
  dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
  dplyr::filter(!is.na(Genes))

gene_names<-upset_wide$Genes

upset_wide <- upset_wide |> dplyr::select(-Genes)
rownames(upset_wide)<-gene_names

upsetRNA <- ComplexUpset::upset(upset_wide, order_upset, name = "", sort_sets = F, themes = ComplexUpset::upset_modify_themes(list(
  'intersections_matrix'=ggplot2::theme(text = ggplot2::element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetRNA <- ggplotify::as.ggplot(upsetRNA)
upsetRNA <- upsetRNA +
  # ggplot2::ggtitle("Upsetplot") +
  ggplot2::theme(plot.title = ggplot2::element_text(
    size = 18,
    hjust = 0.5,
    vjust = 0.95
  ))
upsetRNA
  # tiff(here::here("Upset_WT_subset.tif"), units = "cm", width = 25, height = 20, res = 300)
  # upsetPlot
  # grid::grid.text("Differentially Expressed Genes", x=0.65, y =0.98, gp=grid::gpar(fontsize = 24))
  # dev.off()

 upsetPlot

 #####comparecluster####
 PHSubset <- edgeR_data[4:6]
 for (i in 1:3){
   PHSubset[[i]]<-PHSubset[[i]] |>
     dplyr::filter(FDR<0.05)
   PHSubset[i] <- clusterProfiler::bitr(PHSubset[[i]]$SYMBOL,
                                     fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = "org.Mm.eg.db",
                                     drop = T) |>
   dplyr::select(ENTREZID)
 }
 bg <-clusterProfiler::bitr(edgeR_data[[i]]$SYMBOL,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = "org.Mm.eg.db")

 PHCompare <- clusterProfiler::compareCluster(geneClusters = PHSubset,
                                              fun = "enrichGO",
                                              universe = bg$ENTREZID,
                                              OrgDb = org.Mm.eg.db,
                                              ont = "MF",
                                              readable = T)
 clusterProfiler::dotplot(PHCompare, showCategory = 10)
 PHCompare_CC <- clusterProfiler::compareCluster(geneClusters = PHSubset,
                                              fun = "enrichGO",
                                              universe = bg$ENTREZID,
                                              OrgDb = org.Mm.eg.db,
                                              ont = "CC",
                                              readable = T)
 cc_plot <- clusterProfiler::dotplot(PHCompare_CC, showCategory = 10)+ggplot2::ggtitle("Cellular Component - 10 categories")

 # tiff("CompareClusterCC.tif", units = "cm", width = 30, height = 30, res = 300)
 #  cc_plot
 #  dev.off()

 PHCompare_R <- clusterProfiler::compareCluster(geneClusters = PHSubset,
                                                 fun = "ReactomePA::enrichPathway",
                                                organism = "mouse",
                                                universe = bg$ENTREZID,
                                                readable = T)
 clusterProfiler::dotplot(PHCompare_R, showCategory = 10)


 #####Treeplot####
 go_MF_up <- goAnalysis(edgeR_data[4:6], "MF", "Upregulated")
 go_MF_down <- goAnalysis(edgeR_data[4:6], "MF", "Downregulated")
 go_CC_up <- goAnalysis(edgeR_data[4:6], "CC", "Upregulated")
 go_CC_down <- goAnalysis(edgeR_data[4:6], "CC", "Downregulated")

enrichplot::dotplot(go_MF_down[[3]])
#Treeplots from CC. This code could be made more elegant or turned into a function
treeplot <- enrichplot::pairwise_termsim(go_CC_down$Liv_vs_CS_WT)
enrichplot::treeplot(treeplot)
# tiff("Down_L_vs_CS.tif", units = "cm", width = 25, height = 20, res = 200)
# enrichplot::treeplot(treeplot)+ggtitle("Genes downregulated in liver vs CS")
# dev.off()


treeplot <- enrichplot::pairwise_termsim(go_CC_up$Liv_vs_CS_WT)
enrichplot::treeplot(treeplot)
#no terms

treeplot <- enrichplot::pairwise_termsim(go_CC_up$Liv_vs_PH_WT)
# tiff("Up_L_vs_PH.tif", units = "cm", width = 25, height = 20, res = 200)
# enrichplot::treeplot(treeplot)+ggtitle("Genes upregulated in liver vs PH")
# dev.off()



treeplot <- enrichplot::pairwise_termsim(go_CC_down$Liv_vs_PH_WT)
# tiff("Down_L_vs_PH.tif", units = "cm", width = 25, height = 20, res = 200)
# enrichplot::treeplot(treeplot)+ggtitle("Genes upregulated in PH vs Liver")
# dev.off()

 # tiff(here::here("data/treeplotsBP.tiff"), height = 25, width = 70, res = 200, units = "cm")
 # treeplot_up + treeplot_down
 #  dev.off()



treeplot <- enrichplot::pairwise_termsim(go_CC_down$CS_vs_PH_WT)
enrichplot::treeplot(treeplot)

treeplot <- enrichplot::pairwise_termsim(go_CC_up$CS_vs_PH_WT)
enrichplot::treeplot(treeplot)

#For MF
treeplot <- enrichplot::pairwise_termsim(go_MF_down$Liv_vs_CS_WT)
enrichplot::treeplot(treeplot)

treeplot <- enrichplot::pairwise_termsim(go_MF_up$Liv_vs_CS_WT)
enrichplot::treeplot(treeplot)

treeplot <- enrichplot::pairwise_termsim(go_MF_down$Liv_vs_PH_WT)
enrichplot::treeplot(treeplot)

treeplot <- enrichplot::pairwise_termsim(go_MF_up$Liv_vs_PH_WT)
enrichplot::treeplot(treeplot)

treeplot <- enrichplot::pairwise_termsim(go_MF_down$CS_vs_PH_WT)
enrichplot::treeplot(treeplot)

treeplot <- enrichplot::pairwise_termsim(go_MF_up$CS_vs_PH_WT)
enrichplot::treeplot(treeplot)

#Add reactome analysis?

#####ComplexUpset with colors####
sig_genes_names <- names(edgeR_data[4:6])
sig_genes <-vector(mode = "list", length = 3)
names(sig_genes)<-names(edgeR_data[4:6])
for (i in 1:3){
  sig_genes[[i]]<- edgeR_data[[i+3]]
  sig_genes[[i]]<- sig_genes[[i]] %>%
    dplyr::filter(FDR < 0.05) |>
    dplyr::mutate(Direction =
                    dplyr::case_when(
      logFC > 0 ~ "Up",
      logFC < 0 ~ "Down"
    ))
  sig_genes[[i]]<-sig_genes[[i]] |>
    dplyr::select(SYMBOL, Direction)
}


names(sig_genes) <- c("Liver vs CS","Liver vs PH", "CS vs PH")


#make same upsetplot with ComplexUpset
order_upset <- c("Liver vs CS","Liver vs PH", "CS vs PH")
upset_data <- data.frame("Genes"= sig_genes[[1]]$SYMBOL,
                         "Group"= names(sig_genes[1]),
                         "Direction"= sig_genes[[1]]$Direction)

for (i in 2:length(sig_genes)){
  upset_data <- dplyr::add_row(upset_data,"Genes"=sig_genes[[i]]$SYMBOL,
                               "Group"=names(sig_genes)[i],
                               "Direction" = sig_genes[[i]]$Direction)
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
  dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
  dplyr::filter(!is.na(Genes))

gene_names<-upset_wide$Genes

upset_wide <- upset_wide |> dplyr::select(-Genes)
rownames(upset_wide)<-gene_names

upsetRNA <- ComplexUpset::upset(upset_wide, order_upset, name = "", sort_sets = F, themes = ComplexUpset::upset_modify_themes(list(
  'intersections_matrix'=ggplot2::theme(text = ggplot2::element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetRNA <- ggplotify::as.ggplot(upsetRNA)
upsetRNA <- upsetRNA +
  # ggplot2::ggtitle("Upsetplot") +
  ggplot2::theme(plot.title = ggplot2::element_text(
    size = 18,
    hjust = 0.5,
    vjust = 0.95
  ))
upsetRNA
# tiff(here::here("Upset_WT_subset.tif"), units = "cm", width = 25, height = 20, res = 300)
# upsetPlot
# grid::grid.text("Differentially Expressed Genes", x=0.65, y =0.98, gp=grid::gpar(fontsize = 24))
# dev.off()

upsetPlot

#####Make enrichGO plot####
L_vs_PH <- vector(mode = "list", length = 2)
L_vs_PH[[1]]<- edgeR_data[[5]] |>
  dplyr::filter(logFC >0) |>
  dplyr::filter(FDR<0.05) |>
  dplyr::select(SYMBOL)

L_vs_PH[[2]]<- edgeR_data[[5]] |>
  dplyr::filter(logFC <0) |>
  dplyr::filter(FDR<0.05)|>
  dplyr::select(SYMBOL)
names(L_vs_PH)<- (c("Upregulated in L vs PH","Upregulated in PH vs L"))

bg <- edgeR_data[[1]]
bg <- clusterProfiler::bitr(
  bg$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Mm.eg.db",
  drop = T
)

for (i in 1:2){
  L_vs_PH[[i]] <- clusterProfiler::bitr(
    L_vs_PH[[i]]$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
  )[2]
  L_vs_PH[[i]]<-L_vs_PH[[i]]$ENTREZID
}
  clusterCompare <- clusterProfiler::compareCluster(gene = L_vs_PH,
                                                    fun = "enrichGO",
                                                    universe = bg$ENTREZID,
                                                    key = "ENTREZID",
                                                    OrgDb = "org.Mm.eg.db",
                                                    ont = "CC")
  CompareClusterFigure <- clusterProfiler::dotplot(clusterCompare)

 tiff(here::here("data/figures/compareclusterLvsPH.tiff"), width = 25, height = 20, res = 200, units = "cm")
 CompareClusterFigure+
     ggplot2::ggtitle(label = "Gene Ontology Enrichment",
                      subtitle = "Cellular Component")+
     ggplot2::theme(title = ggplot2::element_text(size = 24,
                                                  hjust = 0.5))
 dev.off()






plot_neg <- enrichplot::dotplot(goResults_neg, font.size = 14)+
  ggtitle("Liver vs PH- Downregulated Genes")
plot_pos <- enrichplot::dotplot(goResults, font.size = 14)+
  ggtitle("Liver vs PH- Upregulated Genes")

# tiff("L_vs_PH_MF.tif", units = "cm", width = 60, height = 20, res = 300)
# plot_pos+plot_neg
# dev.off()
