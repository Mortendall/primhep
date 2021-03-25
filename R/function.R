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
  count_matrix$Geneid<-str_remove_all(count_matrix$Geneid, "\\..*")
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
