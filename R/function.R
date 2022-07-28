library(fs)
library(vroom)
library(Biobase)
library(here)
library(magrittr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)
library(openxlsx)
library(pheatmap)
library(gridExtra)
library(PoiClaClu)
library(RColorBrewer)
library(limma)
library(GO.db)
library(stringr)
library(dplyr)
library(ggplot2)
library(ComplexUpset)

#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly <- function(file_type){
  count_file <- fs::dir_ls(here::here("data_raw/"),
                           regexp = file_type,
                           recurse = TRUE)
  count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file,sheet = 1)
  count_matrix$Geneid<-stringr::str_remove_all(count_matrix$Geneid, "\\..*")
  rownames(count_matrix) <- count_matrix$Geneid
  count_matrix <- count_matrix %>%
    dplyr::select(-Geneid)

  return(count_matrix)
}

#' Load metadata and sort them according to count matrix input
#'
#' @param file_name the name of teh metadata file (default "metadata.csv")
#'
#' @return metadata file sorted according to count matrix order

 load_metadata <- function(file_name) {
   data_file <- fs::dir_ls(here::here("data_raw/"),
                           regexp = file_name,
                           recurse = T)
   metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
   return(metadata)
 }

#' Quality control generator
#'
#' @param count_matrix a count matrix generated through the count_matrix function
#' @param setup setup data.frame
#'
#' @return

Quality_control_plots <- function(count_matrix, setup) {
  group <- as.matrix(setup$Group)
  RNAseq <- edgeR::DGEList(counts = count_matrix, group = group)

  pD <-
    reshape2::melt(cpm(RNAseq, normalized.lib.sizes = TRUE, log = TRUE))
  p <- ggplot(pD, aes(value)) +
    geom_density() +
    facet_wrap( ~ Var2)
  dir.create(here("data/figures"), showWarnings = F)
  dir.create(here("data/figures/QCplots"), showWarnings = F)
  ggplot2::ggsave(
    p,
    filename = here("data/figures/QCplots/Density_plot.png"),
    width = 12,
    height = 12,
    units = "cm",
    scale = 2.5
  )

  #Create mdPlots

  oldpar <- par()$mfrow
  pdf(file.path(here("data/figures/QCplots"), "beforeFiltering_MD.pdf"), width = 4, height = 4)
  par(mfrow = c(2, 2))
  for (i in seq_len(ncol(RNAseq))) {
    plotMD(RNAseq, column = i)
    abline(h = 0)
  }
  par(mfrow = oldpar)
  dev.off()

  #create Possion heatmap
  poisd <- PoissonDistance(t(RNAseq$counts))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- colnames(cpm(RNAseq))
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  heatmap <- pheatmap(samplePoisDistMatrix,  clustering_distance_rows=poisd$dd,clustering_distance_cols=poisd$dd, col=colors)
  ggsave(heatmap, filename = here("data/figures/QCplots/Poisson_heatmap.png"),
         width = 12,
         height = 12,
         units = "cm",
         scale = 2.5)
  #crete mdsPlots
  mdsData <- plotMDS(RNAseq, ndim = 3, plot = FALSE)
  mdsData <-
    mdsData$cmdscale.out %>% data.table(keep.rownames = TRUE) %>%
    mutate(ID = rownames(RNAseq$samples)) %>%
    dplyr::select(-rn) %>%
    mutate(Group = setup$Group)

  setnames(mdsData,
           c("V1", "V2", "V3", "ID", "Group"),
           c("dim1", "dim2", "dim3", "ID", "Group"))
  plotMDS(RNAseq, ndim = 3)
  pBase <-
    ggplot(mdsData, aes(x = dim1, y = dim2, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()
  pBase2 <-
    ggplot(mdsData, aes(x = dim1, y = dim3, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()

  pBase3 <-
    ggplot(mdsData, aes(x = dim2, y = dim3, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()

  pdf(file.path(here("data/figures/QCplots"), "MDSplots.pdf"), width = 4, height = 4)
  par(mfrow = c(1, 1))
  plot(pBase)
  plot(pBase2)
  plot(pBase3)
  plotMDS(RNAseq, ndim = 3)
  par(mfrow = oldpar)
  dev.off()


  #check mds

  #calc norm factors
  keep <- edgeR::filterByExpr(RNAseq)
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- edgeR::calcNormFactors(RNAseq)



  #check MDplots and density after filtering

  pdf(file.path(here("data/figures/QCplots"), "afterFiltering_MD.pdf"), width = 4, height = 4)
  par(mfrow = c(2, 2))
  for (i in seq_len(ncol(RNAseq))) {
    plotMD(RNAseq, column = i)
    abline(h = 0)
  }
  par(mfrow = oldpar)
  dev.off()

  pD <-
    reshape2::melt(cpm(RNAseq, normalized.lib.sizes = TRUE, log = TRUE))
  p <- ggplot(pD, aes(value)) +
    geom_density() +
    facet_wrap( ~ Var2)
  dir.create(here("data/figures/QCplots"), showWarnings = F)
  ggplot2::ggsave(
    p,
    filename = here("data/figures/QCplots/Density_plot_post_filtering.png"),
    width = 12,
    height = 12,
    units = "cm",
    scale = 2.5
  )

  #crete mdsPlots
  mdsData <- plotMDS(RNAseq, ndim = 3, plot = FALSE)
  mdsData <-
    mdsData$cmdscale.out %>% data.table(keep.rownames = TRUE) %>%
    mutate(ID = rownames(RNAseq$samples)) %>%
    dplyr::select(-rn) %>%
    mutate(Group = setup$Group)

  setnames(mdsData,
           c("V1", "V2", "V3", "ID", "Group"),
           c("dim1", "dim2", "dim3", "ID", "Group"))
  plotMDS(RNAseq, ndim = 3)
  pBase <-
    ggplot(mdsData, aes(x = dim1, y = dim2, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()
  pBase2 <-
    ggplot(mdsData, aes(x = dim1, y = dim3, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()

  pBase3 <-
    ggplot(mdsData, aes(x = dim2, y = dim3, colour = Group)) +
    geom_point(size = 5) +
    #geom_label(show.legend = FALSE, size = 5) +
    theme_bw()

  pdf(file.path(here("data/figures/QCplots"), "MDSplots_after_filter.pdf"), width = 4, height = 4)
  par(mfrow = c(1, 1))
  plot(pBase)
  plot(pBase2)
  plot(pBase3)
  plotMDS(RNAseq, ndim = 3)
  par(mfrow = oldpar)
  dev.off()

  print("All your plots can be found in the Figures/QCplots folder")
}

#' RNAseq_processing
#'
#' @param count_matrix generated with the count matrix assembly function
#' @param metadata loaded in with the load_metadata function
#' @param design Generated with Generate_design_matrix function
#' @param ctrsts defined in the WATanalysis script
#'
#' @return a dgeResults list object

RNAseq_processing <- function(count_matrix, metadata, design, ctrsts) {
  group <- as.matrix(metadata$Group)
  RNAseq <- edgeR::DGEList(counts = count_matrix, group = group)
  keep <- edgeR::filterByExpr(RNAseq, design = design, min.count = 20)
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- edgeR::calcNormFactors(RNAseq)
  RNAseq <- edgeR::estimateDisp(RNAseq,design)
  efit <- edgeR::glmQLFit(RNAseq, design)
  dgeResults <- apply(ctrsts, 2, . %>%
                        edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                        edgeR::topTags(n = Inf, p.value = 1) %>%
                        magrittr::extract2("table") %>%
                        data.table::as.data.table(keep.rownames = TRUE))
  return(dgeResults)
}

#' Generate design matrix
#'
#' @param metadata a metadata object generated through the load_metadata function
#'
#' @return a design matrix file


Generate_design_matrix <- function(metadata){
  design <- stats::model.matrix( ~0+Group, metadata)
  colnames(design) <-
    stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")
  return(design)
}


#' camera_reactome
#'
#' @param rLst the reactome database downloaded from https://reactome.org/download/current/Ensembl2Reactome.txt
#' @param organism eg. "Mus musculus"
#' @param count_matrix a cpm matrix generated through edgeR's "cpm" function
#' @param design_matrix design matrix generated thriough stats::niodel.matrix
#' @param contrast_matrix contrast matrix generated with limma::makeContrasts
#'
#' @return a sheet with reactome analysis

camera_reactome <- function(reactomeList,reactomeName, count_matrix, design_matrix, contrast_matrix){
  camera_test <- apply(contrast_matrix, 2, limma::camera, index = reactomeList, y = count_matrix, design = design_matrix)
  camera_test <- lapply(camera_test, data.table, keep.rownames = T)
  camera_test <- lapply(camera_test, setnames, old = "rn", new = "ID")
  camera_test <- lapply(camera_test, magrittr::extract, !is.na(PValue))
  camera_test <- lapply(camera_test, magrittr::extract, reactomeName, on = "ID", nomatch = FALSE)
  camera_test <- lapply(camera_test, magrittr::extract, order(PValue, decreasing = FALSE))
  return(camera_test)
}

#### Gene Ontology (only BP)
#' Title
#'
#' @param org.database ex Org.Mm.eg.db
#' @param cpm_matrix a cpm matrix generated through edgeR's "cpm" function
#' @param design_matrix design matrix generated thriough stats::model.matrix
#' @param contrast_matrix contrast matrix generated with limma::makeContrasts
#'
#' @return a GO analysis

camera_go <- function(org.database, cpm_matrix, design_matrix, contrast_matrix) {

  entrez_matrix <- cpm_matrix
  conv <- clusterProfiler::bitr(rownames(entrez_matrix), fromType='ENSEMBL', toType='ENTREZID', OrgDb = org.database) %>%
    data.table::data.table(key = "ENSEMBL")
  entrez_matrix <- as.data.frame(entrez_matrix) %>%
    dplyr::mutate(Ensembl = rownames(entrez_matrix))
  entrez_matrix <- dplyr::left_join(entrez_matrix, conv, by = c("Ensembl" = "ENSEMBL"))
  entrez_matrix <- dplyr::filter(entrez_matrix,!is.na(entrez_matrix$ENTREZID))
  entrez_matrix <- dplyr::distinct(entrez_matrix, ENTREZID, .keep_all = T)
  rownames(entrez_matrix)<-entrez_matrix$ENTREZID

  entrez_matrix <- dplyr::select(entrez_matrix, -c(ENTREZID,Ensembl))

  GO_test <- apply(ctrsts, 2, camera, index = cyt.go.genes, y = entrez_matrix, design = design)
  GO_test <- lapply(GO_test, data.table, keep.rownames = T)
  GO_test <- lapply(GO_test, setnames, old = "rn", new = "ID")
  GO_test <- lapply(GO_test, magrittr::extract, !is.nan(PValue))
  GO_test <- lapply(GO_test, magrittr::extract, termGO, on = "ID", nomatch = FALSE)
  GO_test <- lapply(GO_test, magrittr::extract, order(PValue, decreasing = FALSE))
  return(GO_test)
}

#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list list of data.tables generated from edgeR. Must be data.table and contain a SYMBOL annotation column
#'
#' @return a list containing enrichresults for each element in the results file list

goAnalysis <- function(result_list, ontology, direction){
  bg <- result_list[[1]]
  bg_list <- clusterProfiler::bitr(
    bg$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
  )

  goResult_list <- vector(mode = "list", length = length(result_list))
  if (direction == "Upregulated"){
    print("Running analysis on upregulated genes")
    for(i in 1:length(result_list)){
      sig_list<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05) |>
        dplyr::filter(logFC>0)
      eg <- clusterProfiler::bitr(
        sig_list$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
      )
      goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                             universe = bg_list$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = ontology)
      goResult_list[[i]]<- goResults
    }


  }

  if(direction == "Downregulated"){
    print("Running analysis on downregulated genes")
    for(i in 1:length(result_list)){
      sig_list<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05) |>
        dplyr::filter(logFC<0)
      eg <- clusterProfiler::bitr(
        sig_list$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
      )
      goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                             universe = bg_list$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = ontology)
      goResult_list[[i]]<- goResults
    }

  }
  if (direction == ""){
    for(i in 1:length(result_list)){
      sig_list<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05)
      eg <- clusterProfiler::bitr(
        sig_list$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
      )
      goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                             universe = bg_list$ENTREZID,
                                             OrgDb = org.Mm.eg.db,
                                             ont = ontology)
      goResult_list[[i]]<- goResults
    }
  }


  for (i in 1:length(goResult_list)){
    names(goResult_list)[i]<-names(result_list)[i]
  }
  return(goResult_list)

}


#' File exporter - exports GOresults as an excel sheet, and prints dotplot and cnet plots
#'
#' @param goList a list object containing one or more enrichResults
#'
#' @return

printGOterms <- function(goList, fileName){
  goSheets<- vector(mode = "list", length = length(goList))
  for (i in 1:length(goSheets)){
    goSheets[[i]] <- goList[[i]]@result
    names(goSheets)[i]<-names(goList)[i]
  }
  openxlsx::write.xlsx(goSheets, file = paste(here::here("data"),"/",fileName,".xlsx",sep = ""), asTable = TRUE)

  for (i in 1:length(goList)){
    dotplot <- enrichplot::dotplot(goList[[i]], title = names(goList)[i])
    ggplot2::ggsave(dotplot, filename = paste(here::here("data/figures"),"/dotplot_",names(goList[i]),".png", sep = ""),width = 12, height = 12, units = "cm", scale = 2.5)
    goList_anno <- clusterProfiler::setReadable(goList[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    cnetplot <- enrichplot::cnetplot(goList_anno, title = names(goList)[i], size = 1)
    ggplot2::ggsave(cnetplot, filename = paste(here("data/figures"),"/cnetplot_",names(goList[i]),".png", sep = ""),scale = 2.5)
  }
}

#
#' code from Lars to extract genes
#'
#' @param tests
#' @param termList
#' @param fromType
#'
#' @return

annotateWithGenes <- function(tests, termList, fromType){
  conv <- bitr(unique(unlist(termList)), fromType = fromType, toType = 'SYMBOL', OrgDb = "org.Mm.eg.db") %>%
    data.table(key = fromType)
  pasteSymbols <- function(x) conv[termList[[x]], paste0(SYMBOL, collapse = ", ")]
  termSymbols <- lapply(names(termList), pasteSymbols) %>% data.table
  termSymbols[, id:=names(termList)]
  setnames(termSymbols, c("symbolList", "id"))
  setkey(termSymbols, id)

  tests <- lapply(tests, copy)
  for (i in names(tests)){
    tests[[i]][, genesInTerm:=termSymbols[ID, symbolList]]
  }
  tests
}

#' reactomeAnalysis
#'
#' @param result_list
#' @param direction
#'
#' @return
#' @export
#'
#' @examples
reactomeAnalysis <- function(result_list, direction){
  bg <- result_list[[1]]
  bg_list <- clusterProfiler::bitr(
    bg$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
  )
  reactomeResult_list <- vector(mode = "list", length = length(result_list))
  sig_list <- vector(mode = "list", length = length(result_list))
  if (direction == "Upregulated"){
    for(i in 1:length(result_list)){
      sig_list[[i]]<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05) |>
        dplyr::filter(logFC>0)
    }

  }

  if(direction == "Downregulated"){
    for(i in 1:length(result_list)){
      sig_list[[i]]<- result_list[[i]] %>%
        dplyr::filter(FDR<0.05) |>
        dplyr::filter(logFC<0)
    }

  }
  else{for(i in 1:length(result_list)){
    sig_list[[i]]<- result_list[[i]] %>%
      dplyr::filter(FDR<0.05)
  }

  }
  for(i in 1:length(result_list)){
    eg <- clusterProfiler::bitr(
      sig_list[[i]]$SYMBOL,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db",
      drop = T
    )
    reactomeResults <- ReactomePA::enrichPathway(gene = eg$ENTREZID,
                                                 universe = bg_list$ENTREZID,
                                                 organism ="mouse",
                                                 readable = T)
    reactomeResult_list[[i]]<- reactomeResults
  }

  for (i in 1:length(reactomeResult_list)){
    names(reactomeResult_list)[i]<-names(result_list)[i]
  }
  return(reactomeResult_list)

}

#' Heatmap Generator clustered
#'
#' @param input_genes list of genes selected from enrichresult
#' @param cpm_matrix count matrix
#' @param setup metadta
#' @param heatmap_title title for heatmap
#'
#'

heatmap_generator_clustered <- function(input_genes, cpm_matrix, setup, heatmap_title){
  #Generate Gene list
  if(length(input_genes==1)){
    gene_list <- input_genes
  gene_list <- unlist(str_split(gene_list, "/"))
  }
  else{
    gene_list <- input_genes
  }

  #Select Candidates in Gene List
  trimmed_cpm <- cpm_matrix |>
    dplyr::filter(SYMBOL %in% gene_list) |>
    dplyr::distinct(SYMBOL, .keep_all = T)
  trimmed_cpm <- trimmed_cpm |>
    dplyr::filter(!is.na(SYMBOL))
  rownames(trimmed_cpm)<-trimmed_cpm$SYMBOL
  trimmed_cpm <- trimmed_cpm |>
    dplyr::select(-SYMBOL)

  #Order metadata

  order <- c("Con_0", "Con_4", "Con_8", "Con_12", "Con_16", "Con_20",
             "Low_0", "Low_4", "Low_8", "Low_12", "Low_16", "Low_20",
             "High_0", "High_4", "High_8", "High_12", "High_16", "High_20")
  setup_ordered <- setup%>%
    dplyr::arrange(desc(Condition1),desc(Condition2))

  #arrange data based on setup sheet
  trimmed_cpm <- trimmed_cpm |>
    dplyr::select(setup_ordered$Sample)

  #create annotation key for heatmap
  key <- as.data.frame(setup_ordered)
  key <- key |>
    dplyr::select(Group)
  rownames(key) <- setup_ordered$ID
  key$Group <-factor(key$Group, c("Con_0", "Con_4", "Con_8", "Con_12", "Con_16", "Con_20",
                                  "Low_0", "Low_4", "Low_8", "Low_12", "Low_16", "Low_20",
                                  "High_0", "High_4", "High_8", "High_12", "High_16", "High_20"))

  #create heatmap
  Heatmap <- pheatmap::pheatmap(trimmed_cpm,
                                treeheight_col = 0,
                                treeheight_row = 0,
                                scale = "row",
                                legend = T,
                                na_col = "white",
                                Colv = NA,
                                na.rm = T,
                                cluster_cols = F,
                                fontsize_row = 5,
                                fontsize_col = 8,
                                cellwidth = 7,
                                cellheight = 5,
                                annotation_col = key,
                                show_colnames = F,
                                show_rownames = T,
                                main = heatmap_title
  )

}

