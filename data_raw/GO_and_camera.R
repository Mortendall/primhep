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
go_sig_genes <- goAnalysis(edgeR_data, "MF")
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
go_sig_genes_MF <- goAnalysis(edgeR_data, "MF")
#printGOterms(go_data)

for (i in 1:length(go_sig_genes_MF)){
  go_sig_genes_MF[[i]] <- clusterProfiler::setReadable(go_sig_genes_MF[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
}
#saveRDS(go_sig_genes_MF, file = here::here("data/go_sig_MF.data.rds"))
test <- readRDS(here::here("data/go_sig_MF.data.rds"))
Liver <- clusterProfiler::dotplot(test[[1]], font.size = 14)+
  ggtitle("Liver - Genotype effect")
#
# tiff("LiverGO_MF.tif", units = "cm", width = 25, height = 15, res = 300)
#  Liver
#  dev.off()

PH <- dotplot(test[[3]], font.size = 14)+
  ggtitle("Primary Hepatocyte - Genotype effect")

 # tiff("PHGO_MF.tif", units = "cm", width = 25, height = 15, res = 300)
 # PH
 # dev.off()

dotplot(test[[5]], font.size = 14)+
  ggtitle("Liver vs Prim hep WT")
dotplot(test[[8]], font.size = 14)+
  ggtitle("Liver vs Prim hep HNKO")

#####Repeat analysis with CC####
go_sig_genes_CC <- goAnalysis(edgeR_data, "CC")
#printGOterms(go_data)

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
