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


cpm_matrix <- readRDS(here("data/cpm_matrix.rds"))
#rLst <- fread("https://reactome.org/download/current/Ensembl2Reactome.txt", header = FALSE)
rLst <- openxlsx::read.xlsx("C:/Users/tvb217/Documents/R/tmp/Reactomedatabase.xlsx")
organism <- "Mus musculus"
reactome_data <- camera_reactome(rLst = rLst, organism = organism, count_matrix = cpm_matrix, design_matrix = design, contrast_matrix = ctrsts)
go_data <- camera_go(org.database = "org.Mm.eg.db", cpm_matrix = cpm_matrix, design_matrix = design, contrast_matrix = ctrsts)

#write.xlsx(camera_test, here::here("data/Reactome_data.xlsx"), asTable = TRUE)
#write.xlsx(GO_test, here::here("data/GO_data.xlsx"), asTable = TRUE)



go_sig_genes <- goAnalysis(edgeR_data)
printGOterms(go_data)
