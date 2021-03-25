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

cpm_data <- read.xlsx(here("data/CPM_matrix.xlsx"),rowNames = T)

go_data <- goAnalysis(edgeR_data)
printGOterms(go_data)
