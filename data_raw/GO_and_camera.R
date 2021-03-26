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
    edgeR_data[[i]]<- open.xlsx::read.xlsx(here("data/edgeR.xlsx"),sheet = i)
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
go_sig_genes <- goAnalysis(edgeR_data)
#printGOterms(go_data)

cameraGoAnnotated <- annotateWithGenes(go_data, cyt.go.genes, fromType = 'ENTREZID')
cameraReactomeAnnotated <- annotateWithGenes(reactome_data, reactomeList, fromType = "ENSEMBL")
#saveRDS(cameraReactomeAnnotated, here("data/annotatedReactome.rds"))
#saveRDS(cameraGoAnnotated, here("data/annotatedGO.rds"))
