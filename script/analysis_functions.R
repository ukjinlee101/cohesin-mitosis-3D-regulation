################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
################################################################################
# THIS IS A COMPREHENSIVE TOOLSET FOT GENOMICS ANALYSIS
# AUTHOR: UK JIN LEE
# APOSTOLOU LAB
# WEILL CORNELL MEDICINE
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
################################################################################
# LIST OF REQUIRED PACKAGES
pkgs <- c("tidyverse", "here", "svglite", "ggplot2", "ggrepel", "stringr",
         "BiocManager", "RColorBrewer", "viridis", "magick", "nlme", "ashr",
         "factoextra", "plotly", "remotes", "cowplot", "gridExtra", "eulerr", 
         "data.table", "clusterProfiler", "ggpubr", "ggsci", "Rtsne", "circlize")
bio_pkgs <- c("biomaRt", "DESeq2", "ComplexHeatmap", "apeglm", "vsn", "sva",
             "limma", "mixOmics", "plotgardener", "org.Mm.eg.db")
# INSTALL PACKAGES IF NECESSARY
#install.packages(pkgs)
#BiocManager::install(bio_pkgs)
#remotes::install_github("EvaYiwenWang/PLSDAbatch")

# LOAD PACKAGES
lapply(pkgs, require, character.only = TRUE)
lapply(bio_pkgs, require, character.only = TRUE)
require(PLSDAbatch)
require("DOSE")
library(rtracklayer)

# CLEANING
rm(pkgs, bio_pkgs)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
################################################################################
# ENSEMBL DATABASE
# v102 is for the latest mm10 version
# ensembl.v102 <- useMart(host = "https://nov2020.archive.ensembl.org",
#                        biomart = "ENSEMBL_MART_ENSEMBL",
#                        dataset = "mmusculus_gene_ensembl")
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
################################################################################
# CREATE DIRECTORY STRUCTURE FOR ANALYSIS
dir.create(here("data"), showWarnings = FALSE)
dir.create(here("r_data"), showWarnings = FALSE)
dir.create(here("figure"), showWarnings = FALSE)
################################################################################
# RNA-seq analysis
################################################################################
rnaseq_makeCountMatrix <- function(fileList, expName){
  # PURPOSE: make a merged count matrix with annotation & TPM normalization
  # INPUT:
  #   fileList: list of "*_featureCounts.txt" files generated from RNA-seq pipeline
  #   expName: experiment name which will be used for naming outputs
  
  sampleList = fileList %>%
    dplyr::mutate(value = str_replace(value, "_featureCounts.txt", ""))
  
  # GENERATE SEED DATABASE
  i = 1
  sample = sampleList$value[[i]]
  sample.tb = as_tibble(read.table(
    file = here("raw_data", fileList$value[[i]]),
    header = TRUE,
    sep = "\t"
  )) %>%
    dplyr::select(c(1, 6, 7))
  colnames(sample.tb) = c("ensembl_gene_id", "length", sample)
  sample.tb = sample.tb %>%
    dplyr::mutate(ensembl_gene_id = str_replace(ensembl_gene_id, "\\.[0-9]*$", ""))
  
  # ADD FOLLOWING DATABASE
  for (i in 2:length(fileList$value)) {
    # MAKE TEMP DATABASE
    sample = sampleList$value[[i]]
    temp.tb = as_tibble(read.table(
      file = here("raw_data", fileList$value[[i]]),
      header = TRUE,
      sep = "\t"
    )) %>%
      dplyr::select(c(1, 7))
    colnames(temp.tb) = c("ensembl_gene_id", sample)
    temp.tb = temp.tb %>%
      dplyr::mutate(ensembl_gene_id = str_replace(ensembl_gene_id, "\\.[0-9]*$", ""))
    
    # MERGE DATABASE
    sample.tb = sample.tb %>% full_join(temp.tb, by = "ensembl_gene_id")
  }
  
  # CONVERTING ensembl_gene_id TO external_gene_name
  mm10.ensembl2gene = as_tibble(
    getBM(
      mart = ensembl.v102,
      values = sample.tb$ensembl_gene_id,
      filters = "ensembl_gene_id",
      attributes = c("ensembl_gene_id", "external_gene_name")
    )
  )
  
  sample.tb = sample.tb %>%
    dplyr::full_join(mm10.ensembl2gene, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
    dplyr::rename(ensembl = ensembl_gene_id, gene = external_gene_name) %>%
    dplyr::relocate(gene, .after = ensembl)
  
  fwrite(sample.tb, here("data", paste0("readCount.", expName, ".tsv")), sep = '\t')

  # RPM NORMALIZATION
  temp.tb = sample.tb %>% dplyr::select(-c(1, 2, 3))
  normFactor = 1000000 / colSums(temp.tb)
  sample.RPM.tb = as_tibble(sweep(temp.tb, 2, normFactor, "*"))
  temp = sample.tb %>% dplyr::select(c(1, 2))
  sample.RPM.tb = bind_cols(temp, sample.RPM.tb)
  
  fwrite(sample.RPM.tb, here("data", paste0("readCount.RPM.", expName, ".tsv")), sep = '\t')
  
  # TPM NORMALIZATION
  temp.tb = sample.tb %>% dplyr::select(-c(1, 2, 3))
  
  ## 1_LENGTH NORMALIZATION
  length.tb = sample.tb %>% dplyr::select(c(3)) %>%
    dplyr::mutate(length = length / 1000)
  temp.tb = as_tibble(mapply("/", temp.tb, length.tb))
  
  ## 2_DEPTH NORMALIZATION
  normFactor = 1000000 / colSums(temp.tb)
  sample.TPM.tb = as_tibble(sweep(temp.tb, 2, normFactor, "*"))
  temp = sample.tb %>% dplyr::select(c(1, 2))
  sample.TPM.tb = bind_cols(temp, sample.TPM.tb)
  
  fwrite(sample.TPM.tb, here("data", paste0("readCount.TPM.", expName, ".tsv")), sep = '\t')
}

rnaseq_filterLowRPM <- function(expName, RPMCutoff, sampleNumCutoff, ncol = 3,
                                parentDir = "figure",
                                dataDir = "data"){
  # PURPOSE: Use genes that has RPM > RPMCutoff in at least sampleNumCutoff samples 
  # to remove improve quality and reliability of downstream analysis
  # INPUT:
  #   expName: experiment name which were used for naming outputs
  #   RPMCutoff: minimum RPM number required
  #   sampleNumCutoff: minimum number of samples over RPMCutoff
  #   width: width of histogram showing readcount distribution
  #   height: height of histogram showing readcount distribution
  
  dir.create(here(parentDir, "histogram"), showWarnings = FALSE, recursive = TRUE)
  
  sample.tb <- fread(here(dataDir, paste0("readCount.", expName, ".tsv")))
  sample.RPM.tb <- fread(here(dataDir, paste0("readCount.RPM.", expName, ".tsv")))
  sample.TPM.tb <- fread(here(dataDir, paste0("readCount.TPM.", expName, ".tsv")))
  ##############################################################################
  ### Plotting the distribution of raw readCounts
  fileName <- here(parentDir, "histogram", paste0("Histogram_readCounts.", expName))
  title <- "Raw readcount distribution"
  
  options(scipen=999)
  sampleNum <- ncol(sample.tb)-3
  
  plotList <- list()
  for(i in seq(1, sampleNum)){
    data <- data.frame(sample.tb[[c(i+3)]])
    colnames(data) <- "value"
    sampleName <- colnames(sample.tb)[i+3]
    plot <- ggplot(data, aes(x = value)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white") +
      scale_x_continuous(trans = "log10", limits = c(1, NA)) + theme_classic() +
      geom_density(alpha=.2, fill="#FF6666") +
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 6))+
      labs(x = "readCount", title = sampleName)
    plotList[[i]] <- plot
  }
  
  width <- 2*ncol
  height <- 1.5*ceiling(sampleNum/ncol)
  
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  
  svglite(paste0(fileName, ".svg"), width = width, height = height)
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  dev.off()
  
  png(paste0(fileName, ".png"), width = width, height = height, units = "in", res = 600)
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  dev.off()
  
  ##############################################################################
  # Plotting the distribution of RPM before filtering 
  fileName <- here(parentDir, "histogram", paste0("Histogram_RPM.", expName))
  title <- "RPM distribution"
  
  options(scipen=999)
  sampleNum <- ncol(sample.RPM.tb)-2
  
  plotList <- list()
  for(i in seq(1, sampleNum)){
    data <- data.frame(sample.RPM.tb[[c(i+2)]])
    colnames(data) <- "value"
    sampleName <- colnames(sample.RPM.tb)[i+2]
    plot <- ggplot(data, aes(x = value)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white") +
      scale_x_continuous(trans = "log10", limits = c(1, NA)) + theme_classic() +
      geom_density(alpha=.2, fill="#FF6666") +
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 6))+
      labs(x = "RPM", title = sampleName)
    plotList[[i]] <- plot
  }
  
  width <- 2*ncol
  height <- 1.5*ceiling(sampleNum/ncol)
  
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  
  svglite(paste0(fileName, ".svg"), width = width, height = height)
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  dev.off()
  
  png(paste0(fileName, ".png"), width = width, height = height, units = "in", res = 600)
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  dev.off()
  
  ##############################################################################
  #RPMCutoff <- 1
  #sampleNumCutoff <- 2
  # Plotting the distribution of RPM after filtering 
  fileName <- here(parentDir, "histogram", paste0("Histogram_RPM_filtered.", expName))
  title <- paste0("RPM distribution after filtering...RPMCutoff: ", RPMCutoff, " / sampleNumCutoff: ", sampleNumCutoff)
  
  data.tb <- sample.RPM.tb %>% dplyr::select(seq(3, 2+sampleNum))
  
  
  
  # FILTERING
  filter1 <- apply(data.tb, 1, function(x)
    length(x[x > RPMCutoff]) >= sampleNumCutoff)
  filter2 <- apply(data.tb, 1, function(x)
    length(x[x > RPMCutoff]) < sampleNumCutoff)
  data.tb.filtered1 <- data.tb[filter1, ]
  data.tb.filtered2 <- data.tb[filter2, ]
  sample.filtered.tb <- sample.tb[filter1, ]
  sample.RPM.filtered.tb <- sample.RPM.tb[filter1, ]
  sample.RPM.filteredOut.tb <- sample.RPM.tb[filter2, ]
  sample.TPM.filtered.tb <- sample.TPM.tb[filter1, ]
  
  fwrite(sample.filtered.tb, here(dataDir, paste0("readCount.filtered.", expName, ".tsv")), sep = '\t')
  fwrite(sample.RPM.filtered.tb, here(dataDir, paste0("readCount.filtered.RPM.", expName, ".tsv")), sep = '\t')
  fwrite(sample.TPM.filtered.tb, here(dataDir, paste0("readCount.filtered.TPM.", expName, ".tsv")), sep = '\t')
  
  # PLOTTING
  options(scipen=999)
  sampleNum <- ncol(sample.RPM.tb)-2
  
  plotList <- list()
  for(i in seq(1, sampleNum)){
    data <- data.frame(sample.RPM.filtered.tb[[c(i+2)]])
    colnames(data) <- "value"
    sampleName <- colnames(sample.RPM.filtered.tb)[i+2]
    plot <- ggplot(data, aes(x = value)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white") +
      scale_x_continuous(trans = "log10", limits = c(1, NA)) + theme_classic() +
      geom_density(alpha=.2, fill="#FF6666") +
      theme(axis.title = element_text(size = 8), plot.title = element_text(size = 6))+
      labs(x = "RPM", title = sampleName)
    plotList[[i]] <- plot
  }
  
  width <- 2*ncol
  height <- 1.5*ceiling(sampleNum/ncol)
  
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  
  svglite(paste0(fileName, ".svg"), width = width, height = height)
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  dev.off()
  
  png(paste0(fileName, ".png"), width = width, height = height, units = "in", res = 600)
  grid.arrange(grobs = plotList, ncol = 3, top = title)
  dev.off()
  
  print(paste0(nrow(sample.RPM.filtered.tb), " genes left after filtering"))
}

rnaseq_batchCorrect_SVA <- function(sample.tb, batch, covar_mat){
  # PURPOSE: To do batch correction
  # to remove improve quality and reliability of downstream analysis
  # INPUT:
  #   expName: experiment name which were used for naming outputs
  #   RPMCutoff: minimum RPM number required
  #   sampleNumCutoff: minimum number of samples over RPMCutoff
  #   width: width of histogram showing readcount distribution
  #   height: height of histogram showing readcount distribution
  
  
  count_matrix = as.matrix(sample.tb %>% dplyr::select(-c(1, 2, 3)))
  
  adjusted_counts.cov = ComBat_seq(counts = count_matrix, batch = batch, group = NULL, covar_mod = covar_mat)
  adjusted_counts.cov = cbind(sample.tb %>% dplyr::select(c(1, 2, 3)),
                              as_tibble(adjusted_counts.cov))
  
  fwrite(adjusted_counts.cov, file = here("data", paste0("readCount.filtered.ComBatSeq.", expName, ".tsv")), sep = "\t")
  
  sample.tb = as_tibble(adjusted_counts.cov)
  
  # RPM NORMALIZATION
  temp.tb = sample.tb %>% dplyr::select(-c(1, 2, 3))
  normFactor = 1000000 / colSums(temp.tb)
  sample.RPM.tb = as_tibble(sweep(temp.tb, 2, normFactor, "*"))
  temp = sample.tb %>% dplyr::select(c(1, 2))
  sample.RPM.tb = bind_cols(temp, sample.RPM.tb)
  
  fwrite(sample.RPM.tb, here("data", paste0("readCount.filtered.ComBatSeq.RPM.", expName, ".tsv")), sep = '\t')
  
  # TPM NORMALIZATION
  temp.tb = sample.tb %>% dplyr::select(-c(1, 2, 3))
  
  ## 1_LENGTH NORMALIZATION
  length.tb = sample.tb %>% dplyr::select(c(3)) %>%
    dplyr::mutate(length = length / 1000)
  temp.tb = as_tibble(mapply("/", temp.tb, length.tb))
  
  ## 2_DEPTH NORMALIZATION
  normFactor = 1000000 / colSums(temp.tb)
  sample.TPM.tb = as_tibble(sweep(temp.tb, 2, normFactor, "*"))
  temp = sample.tb %>% dplyr::select(c(1, 2))
  sample.TPM.tb = bind_cols(temp, sample.TPM.tb)
  
  fwrite(sample.TPM.tb, here("data", paste0("readCount.filtered.ComBatSeq.TPM.", expName, ".tsv")), sep = '\t')
}


rnaseq_runDESeq <- function(data.tb, expName, note, colData, design, dataDir = "data"){
  # PURPOSE: Use filtered readCount table to run differential analysis 
  # INPUT:
  #   data.tb: filtered raw readCount table from `rnaseq_filterLowRPM`
  #   expName: experiment name which will be used for naming outputs
  #   colData: colData for DESeq2
  #   design: design for DESeq2
  
  #########################################################
  # STEP 4: Construct DESeq2 object
  #########################################################
  dds <- DESeqDataSetFromMatrix(data.tb,
                                colData,
                                design = design)
  #########################################################
  # STEP 5: Run DESeq2
  #########################################################
  dds <- DESeq(dds)
  readCount.filtered.DESeqNorm <- counts(dds, normalized = TRUE)
  temp <- as_tibble(readCount.filtered.DESeqNorm, rownames = "ensembl")
  fwrite(temp, here(dataDir, paste0("readCount.filtered.DESeqNorm.", expName, ".tsv")), sep = '\t')
  dir.create(here(dataDir, "r_data"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(dds, file = here(dataDir, "r_data", paste0("dds_", expName, "_", note, ".rds")))
  return(dds)
}

rnaseq_runDiff <- function(dds, c1, c2, alpha, foldchange, log2FC_cutoff, width, height, dataDir = "data", figDir = "figure"){
  contrast <- c("group", c1, c2)
  output_prefix <- paste0(expName, "_", c1, "_vs_", c2)
  rnaseq_pairwiseDiff(dds, contrast, output_prefix, alpha = alpha, foldchange = foldchange, dataDir = dataDir)
  create_volcanoMA(dds, output_prefix, alpha = alpha, foldchange = foldchange,
                   log2FC_cutoff = log2FC_cutoff, width = width, height = height, dataDir = dataDir, figDir = figDir)
  
}

rnaseq_pairwiseDiff <- function(dds, contrast, output_prefix, alpha = 0.05, foldchange = 0.5, dataDir = "data") {
  readCount.filtered.DESeqNorm <- counts(dds, normalized = TRUE)
  data.tb <- as_tibble(readCount.filtered.DESeqNorm, rownames = "ensembl")
  
  res <- results(dds, alpha = alpha, contrast = contrast)
  res.tb <- as_tibble(res)
  res.tb$ensembl_gene_id <- rownames(res)
  
  resLFC <- lfcShrink(dds, contrast = contrast, type = "ashr") # Log fold shrinkage, better for gene ranking
  resLFC.tb <- as_tibble(resLFC)
  resLFC.tb$ensembl_gene_id <- rownames(resLFC)
  temp <- resLFC.tb %>% dplyr::mutate(shrinked_log2FC = log2FoldChange) %>%
    dplyr::select(ensembl_gene_id, shrinked_log2FC)
  
  mergedData <- full_join(res.tb, temp, by = join_by("ensembl_gene_id"))
  
  resLFCSort.tb <- mergedData %>% dplyr::arrange(padj)
  
  # mm10.ensembl2gene <- as_tibble(
  #   getBM(
  #     mart = ensembl.v102,
  #     values = resLFCSort.tb$ensembl_gene_id,
  #     filters = "ensembl_gene_id",
  #     attributes = c("ensembl_gene_id", "external_gene_name")
  #   )
  # )
  
  gtf <- import("/Volumes/UKJIN_SSD/p300/reference/mm10/Mus_musculus.GRCm38.102.gtf")

  mm10.ensembl2gene <- gtf %>%
    as_tibble() %>%
    filter(type == "gene") %>%
    dplyr::select(
      ensembl_gene_id = gene_id,
      external_gene_name = gene_name
    ) %>%
    distinct() %>%
    filter(ensembl_gene_id %in% resLFCSort.tb$ensembl_gene_id)

  resLFCSort.tb <- resLFCSort.tb %>%
    dplyr::full_join(mm10.ensembl2gene, by = c("ensembl_gene_id")) %>%
    dplyr::relocate(ensembl_gene_id, .before = baseMean) %>%
    dplyr::relocate(external_gene_name, .before = baseMean) %>%
    dplyr::relocate(shrinked_log2FC, .after = log2FoldChange)
  
  resLFCSort.tb$diffExpressed <- "NO"
  resLFCSort.tb$diffExpressed[resLFCSort.tb$shrinked_log2FC > foldchange &
                                resLFCSort.tb$padj < alpha] <- "UP"
  resLFCSort.tb$diffExpressed[resLFCSort.tb$shrinked_log2FC < -foldchange &
                                resLFCSort.tb$padj < alpha] <- "DOWN"
  
  out.tb <- resLFCSort.tb
  out.tb <- left_join(out.tb, data.tb, 
                      by = c("ensembl_gene_id" = "ensembl"))
  
  # Save data
  saveRDS(resLFCSort.tb, file = here(dataDir, "r_data", paste0("resLFCSort.tb_", output_prefix, ".rds")))
  fwrite(out.tb, file = here(dataDir, "diff", paste0("diff_", output_prefix, ".tsv")), sep = "\t")
  
  # Save summary
  sink(file = here(dataDir, paste0("diff_", output_prefix, "_summary.txt")))
  summary(res)
  sink(file = NULL)
}

create_PCA_DESeq2 <- function(dds, colData, varThreshold = 0.25, title, label = "group", deseq2BatchCorrect = FALSE, width = 6, height = 4,
                              figDir = "figure"){
  # PURPOSE: data used here is `rld` transformed DESeq2 normalized reads
  # INPUT:
  #   dds: rld transformed DESeq2 normalized reads
  #   colData: colData from DESeq2
  #   varThreshold: percentage of top variable genes to use for PCA
  
  rld <- rlog(dds, blind = FALSE) # log transform option 2 with variance shrinkage
  data <- assay(rld)
  
  if(deseq2BatchCorrect){
    # APPLY BATCH CORRECTION FROM DESIGN
    data <- assay(rld)
    mm <- model.matrix(~group, colData(rld))
    data <- limma::removeBatchEffect(data, batch = rld$batch, design = mm)
    assay(rld) <- data
  }
  variances <- rowVars(data)
  threshold <- quantile(variances, 1-varThreshold) 
  selected_genes <- variances > threshold
  data_selected <- data[selected_genes, ]
  
  PCA_abundance_table <- prcomp(t(data_selected), scale. = FALSE)
  std_devs <- PCA_abundance_table$sdev
  variances <- std_devs^2
  proportion_variance <- variances / sum(variances)
  cumulative_variance <- cumsum(proportion_variance)
  pca_contributions <- data.frame(
    PC = paste0("PC", 1:length(proportion_variance)),
    Variance = variances,
    Proportion_Variance = proportion_variance,
    Cumulative_Variance = cumulative_variance
  )
  
  PCA_res <- as.data.frame(PCA_abundance_table$x) %>% rownames_to_column(var = "Sample")
  colData <- as.data.frame(colData) %>% rownames_to_column(var = "Sample")
  merged_data <- left_join(PCA_res, colData, by = "Sample")
  if(label == "group"){
    p <- ggplot(data = merged_data, aes(x = PC1, y = PC2, color = group)) +
      geom_point(size = 3) +
      scale_color_npg() +
      geom_text_repel(aes(label = group)) +
      labs(title = paste0(title, " (using top ", varThreshold*100, "% variable genes)"), 
           x = paste0("PC1: ", round(pca_contributions$Proportion_Variance[1]*100, 1), "% variance"),
           y = paste0("PC2: ", round(pca_contributions$Proportion_Variance[2]*100, 1), "% variace")) +
      theme_pubr() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
      )
    fileName <- paste0("PCA_", title)
    
  }else{
    p <- ggplot(data = merged_data, aes(x = PC1, y = PC2, color = group)) +
      geom_point(size = 3) +
      scale_color_npg() +
      geom_text_repel(aes(label = label)) +
      labs(title = paste0(title, " (using top ", varThreshold*100, "% variable genes)"), 
           x = paste0("PC1: ", round(pca_contributions$Proportion_Variance[1]*100, 1), "% variance"),
           y = paste0("PC2: ", round(pca_contributions$Proportion_Variance[2]*100, 1), "% variace")) +
      theme_pubr() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
      )
    fileName <- paste0("PCA_", title, "_label")
    
  }
  print(p)
  
  dir.create(here(figDir, "PCA"), showWarnings = FALSE, recursive = TRUE)
  svglite(here(figDir, "PCA", paste0(fileName, ".svg")), width = width, height = height)
  print(p)
  dev.off()
  
  png(here(figDir, "PCA", paste0(fileName, ".png")), width = width, height = height, units = "in", res = 600)
  print(p)
  dev.off()
}

create_tSNE_DESeq2 <- function(dds, colData, varThreshold = 0.25, perplexity = 1, title, label = "group", deseq2BatchCorrect = FALSE, width = 4, height = 4){
  # PURPOSE: data used here is `rld` transformed DESeq2 normalized reads
  # INPUT:
  #   data: rld transformed DESeq2 normalized reads
  #   colData: colData from DESeq2
  #   varThreshold: percentage of top variable genes to use for PCA
  rld <- rlog(dds, blind = FALSE) # log transform option 2 with variance shrinkage
  
  if(deseq2BatchCorrect){
    # APPLY BATCH CORRECTION FROM DESIGN
    data <- assay(rld)
    mm <- model.matrix(~group, colData(rld))
    data <- limma::removeBatchEffect(data, batch = rld$batch, design = mm)
    assay(rld) <- data
  }
  
  data <- assay(rld)
  variances <- rowVars(data)
  threshold <- quantile(variances, 1-varThreshold) 
  selected_genes <- variances > threshold
  data_selected <- data[selected_genes, ]
  
  colData <- as.data.frame(colData) %>% rownames_to_column(var = "Sample")
  set.seed(123)
  tsne_model_1 <- Rtsne(t(as.matrix(data_selected)), check_duplicates = FALSE, pca = TRUE, perplexity = perplexity, theta = 0.1, dims = 2)
  d_tsne_1 <- as.data.frame(tsne_model_1$Y) %>% dplyr::mutate(Sample = colData$Sample)
  
  merged_data <- left_join(d_tsne_1, colData, by = "Sample")
  if(label == "group"){
    
    p <- ggplot(data = merged_data, aes(x = V1, y = V2, color = group)) +
      geom_point(size = 3) +
      scale_color_npg() +
      geom_text_repel(aes(label = group)) +
      labs(title = paste0(title, " (using top ", varThreshold*100, "% variable genes, perplexity: ", perplexity, ")"), 
           x = "V1",
           y = "V2") +
      theme_pubr() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
      )
    fileName <- paste0("tSNE_", title)
    
  }else{
    p <- ggplot(data = merged_data, aes(x = V1, y = V2, color = group)) +
      geom_point(size = 3) +
      scale_color_npg() +
      geom_text_repel(aes(label = label)) +
      labs(title = paste0(title, " (using top ", varThreshold*100, "% variable genes, perplexity: ", perplexity, ")"), 
           x = "V1",
           y = "V2") +
      theme_pubr() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
      )
    fileName <- paste0("tSNE_", title, "_label")
    
  }
  print(p)
  
  dir.create(here("figure", "tSNE"), showWarnings = FALSE)
  svglite(here("figure", "tSNE", paste0(fileName, ".svg")), width = width, height = height)
  print(p)
  dev.off()
  
  png(here("figure", "tSNE", paste0(fileName, ".png")), width = width, height = height, units = "in", res = 600)
  print(p)
  dev.off()
  
}

create_distHeatmap_DESeq2 <- function(dds, title, deseq2BatchCorrect = FALSE, width = 7, height = 4,
                                      figDir = "figure"){
  # PURPOSE: Calculate sample-to-sample distance and draw heatmap
  # INPUT:
  #   dds: DESeq2 output
  
  rld <- rlog(dds)
  fileName <- paste0("heatmap_", title)
  
  if(deseq2BatchCorrect){
    # APPLY BATCH CORRECTION FROM DESIGN
    data <- assay(rld)
    mm <- model.matrix(~group, colData(rld))
    data <- limma::removeBatchEffect(data, batch = rld$batch, design = mm)
    assay(rld) <- data
  }
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colData <- colData(rld)
  rownames(sampleDistMatrix) <-(as_tibble(colData) %>% dplyr::select(label) %>% pull(1))
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  p <- pheatmap(sampleDistMatrix,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                col=colors,
                main = title)
  
  dir.create(here(figDir, "heatmap"), showWarnings = FALSE, recursive = TRUE)
  svglite(here(figDir, "heatmap", paste0(fileName, ".svg")), width = width, height = height)
  print(p)
  dev.off()
  
  png(here(figDir, "heatmap", paste0(fileName, ".png")), width = width, height = height, units = "in", res = 600)
  print(p)
  dev.off()
  
  print(p)
}

create_volcanoMA <- function(dds, output_prefix, 
                             alpha = 0.05, foldchange = 0.5, log2FC_cutoff = 5,
                             width = 4, height = 3.5, dataDir = "data", figDir = "figure"){
  dir.create(here(figDir, "MA"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here(figDir, "volcano"), showWarnings = FALSE, recursive = TRUE)
  
  # MA plot
  resLFCSort.tb <- readRDS(here(dataDir,
    "r_data",
    paste0("resLFCSort.tb_", output_prefix, ".rds")
  ))
  
  resLFCSort.tb$label <- NA
  resLFCSort.tb$label[resLFCSort.tb$diffExpressed != "NO"] <- resLFCSort.tb$external_gene_name[resLFCSort.tb$diffExpressed != "NO"]
  resLFCSort.tb <- resLFCSort.tb %>% dplyr::filter(!is.na(padj))
  
  mycolors <- c("blue", "red", "grey")
  names(mycolors) <- c("DOWN", "UP", "NO")
  
  down.num <- table(resLFCSort.tb$diffExpressed)["DOWN"]
  up.num <- table(resLFCSort.tb$diffExpressed)["UP"]
  
  # Apply limit to log2FC
  resLFCSort.tb <- resLFCSort.tb %>% dplyr::mutate(log2FoldChange_MAX = ifelse(abs(log2FoldChange) >= log2FC_cutoff,
                                                                               sign(log2FoldChange)*log2FC_cutoff,
                                                                               log2FoldChange),
                                                   shrinked_log2FC_MAX = ifelse(abs(shrinked_log2FC) >= log2FC_cutoff,
                                                                                sign(shrinked_log2FC)*log2FC_cutoff,
                                                                                shrinked_log2FC))
  
  plot <-  ggplot(data = resLFCSort.tb, aes(x = baseMean, y = log2FoldChange_MAX,
                                      col = diffExpressed, label = label)) +
    geom_point(
      size = 0.7,
      fill = 'black',
      alpha = 1,
      shape = ifelse(abs(resLFCSort.tb$log2FoldChange) >= log2FC_cutoff, 2, 16)
    ) + 
    geom_text_repel(
      nudge_x = 0.1,
      nudge_y = 0.1,
      size = 1.5,
      hjust = 0, force_pull = 20
    )+
    scale_x_log10() +
    scale_y_continuous(breaks = seq(-log2FC_cutoff, log2FC_cutoff)) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "black") +
    ggtitle(paste0(
      "[MA Plot]", "\n",
      output_prefix, "\n",
      "padj < ", alpha, ", abs(log2FC) > ", foldchange, "\n",
      "down (", down.num, "), up (", up.num, ")"
    )) +
    xlab(paste0("baseMean")) +
    ylab(paste0("log2FC")) +
    theme_classic() +
    theme(
      aspect.ratio = 0.8,
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title = element_text(size = 8),
      text = element_text(size = 8)
    ) + scale_colour_manual(values = mycolors)
  
  
  fileName = here(figDir, "MA", paste0("MA_", output_prefix, "_", alpha, "_", foldchange))
  
  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot)
  dev.off()
  
  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot)
  dev.off()
  print(plot)
  
  ##############################################################################
  # Volcano
  plot <- ggplot(data = resLFCSort.tb,
                 aes(
                   x = shrinked_log2FC_MAX,
                   y = -log10(padj),
                   col = diffExpressed,
                   label = label
                 )) + xlim(-log2FC_cutoff, log2FC_cutoff) +
    geom_point(size = 0.7,
               alpha = 1,
               fill = 'black',
               shape = ifelse(abs(resLFCSort.tb$shrinked_log2FC) >= log2FC_cutoff, 2, 16)) +
    geom_text_repel(
      nudge_x = 0.1,
      nudge_y = 0.1,
      size = 1.5,
      hjust = 0, force_pull = 20
    ) +
    geom_vline(
      xintercept = c(-foldchange, foldchange),
      linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "black"
    ) +
    scale_colour_manual(values = mycolors) +
    ggtitle(
      paste0(
        "[Volcano Plot]", "\n",
        output_prefix, "\n",
        "padj < ", alpha, ", abs(log2FC) > ", foldchange, "\n",
        "down (", down.num, "), up (", up.num, ")"
      )) +
    theme_classic() +
    theme(
      aspect.ratio = 0.8,
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title = element_text(size = 8),
      text = element_text(size = 8)
    ) + coord_cartesian(clip = "off") +
    xlab("shrinked log2fc")
  
  plot2 <- ggplot(data = resLFCSort.tb,
                  aes(
                    x = log2FoldChange_MAX,
                    y = -log10(padj),
                    col = diffExpressed,
                    label = label
                  )) + xlim(-log2FC_cutoff, log2FC_cutoff) +
    geom_point(size = 0.7,
               alpha = 1,
               fill = 'black',
               shape = ifelse(abs(resLFCSort.tb$log2FoldChange) >= log2FC_cutoff, 2, 16)) +
    geom_text_repel(
      nudge_x = 0.1,
      nudge_y = 0.1,
      size = 1.5,
      hjust = 0, force_pull = 20
    ) +
    geom_vline(
      xintercept = c(-foldchange, foldchange),
      linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "black"
    ) +
    scale_colour_manual(values = mycolors) +
    ggtitle(
      paste0(
        "[Volcano Plot]", "\n",
        output_prefix, "\n",
        "padj < ", alpha, ", abs(log2FC) > ", foldchange, "\n",
        "down (", down.num, "), up (", up.num, ")"
      )) +
    theme_classic() +
    theme(
      aspect.ratio = 0.8,
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title = element_text(size = 8),
      text = element_text(size = 8)
    ) + coord_cartesian(clip = "off") +
    xlab("log2fc")
  fileName = here(figDir, "MA", paste0("MA_", output_prefix, "_", alpha, "_", foldchange))
  
  fileName = here(figDir, "volcano", paste0("volcano_", output_prefix, "_", alpha, "_", foldchange, "_shrlog2FC"))
  
  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot)
  dev.off()
  
  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot)
  dev.off()
  
  fileName = here(figDir, "volcano", paste0("volcano_", output_prefix, "_", alpha, "_", foldchange, "_log2FC"))
  
  
  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot2)
  dev.off()
  
  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot2)
  dev.off()
  print(plot)
  print(plot2)
}

create_volcanoHighlighted <- function(dds, output_prefix, 
                                      alpha = 0.05, foldchange = 0.5, log2FC_cutoff = 5,
                                      width = 4, height = 3.5, geneList, note, color = "purple"){
  
  dir.create(here("figure", "volcano"), showWarnings = FALSE, recursive = TRUE)
  
  resLFCSort.tb <- readRDS(here(
    "r_data",
    paste0("resLFCSort.tb_", output_prefix, ".rds")
  ))
  
  resLFCSort.tb <- resLFCSort.tb %>%
    dplyr::mutate(flag =if_else(external_gene_name %in% geneList, "YES", "NO"))
  resLFCSort.tb$flag <- factor(resLFCSort.tb$flag, levels = c("YES", "NO"))
  resLFCSort.tb$label <- NA
  
  resLFCSort.tb$label[resLFCSort.tb$flag != "NO"] <- resLFCSort.tb$external_gene_name[resLFCSort.tb$flag != "NO"]
  resLFCSort.tb <- resLFCSort.tb %>% dplyr::filter(!is.na(padj)) %>% dplyr::arrange(desc(flag))
  
  
  resLFCSort.tb = resLFCSort.tb %>% dplyr::filter(!is.na(padj)) %>% dplyr::arrange(desc(flag))
  mycolors = c(color, "grey")
  names(mycolors) = c("YES", "NO")     
  
  # Apply limit to log2FC
  resLFCSort.tb <- resLFCSort.tb %>% dplyr::mutate(log2FoldChange_MAX = ifelse(abs(log2FoldChange) >= log2FC_cutoff,
                                                                               sign(log2FoldChange)*log2FC_cutoff,
                                                                               log2FoldChange),
                                                   shrinked_log2FC_MAX = ifelse(abs(shrinked_log2FC) >= log2FC_cutoff,
                                                                                sign(shrinked_log2FC)*log2FC_cutoff,
                                                                                shrinked_log2FC))
  
  plot <- ggplot(data = resLFCSort.tb,
                 aes(
                   x = shrinked_log2FC_MAX,
                   y = -log10(padj),
                   col = flag,
                   label = label
                 )) + xlim(-log2FC_cutoff, log2FC_cutoff) +
    geom_point(size = 0.7,
               alpha = 1,
               fill = 'black',
               shape = ifelse(abs(resLFCSort.tb$shrinked_log2FC) >= log2FC_cutoff, 2, 16)) +
    geom_text_repel(
      nudge_x = 0.1,
      nudge_y = 0.1,
      size = 1.5,
      hjust = 0, force_pull = 20
    ) +
    geom_vline(
      xintercept = c(-foldchange, foldchange),
      linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "black"
    ) +
    scale_colour_manual(values = mycolors) +
    ggtitle(
      paste0(
        "[Volcano Plot]", "\n",
        output_prefix, "\n",
        "padj < ", alpha, ", abs(log2FC) > ", foldchange, "\n"
      )) +
    theme_classic() +
    theme(
      aspect.ratio = 0.8,
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title = element_text(size = 8),
      text = element_text(size = 8)
    ) + coord_cartesian(clip = "off") +
    xlab("shrinked log2fc")
  
  plot2 <- ggplot(data = resLFCSort.tb,
                  aes(
                    x = log2FoldChange_MAX,
                    y = -log10(padj),
                    col = flag,
                    label = label
                  )) + xlim(-log2FC_cutoff, log2FC_cutoff) +
    geom_point(size = 0.7,
               alpha = 1,
               fill = 'black',
               shape = ifelse(abs(resLFCSort.tb$log2FoldChange) >= log2FC_cutoff, 2, 16)) +
    geom_text_repel(
      nudge_x = 0.1,
      nudge_y = 0.1,
      size = 1.5,
      hjust = 0, force_pull = 20
    ) +
    geom_vline(
      xintercept = c(-foldchange, foldchange),
      linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "black"
    ) +
    scale_colour_manual(values = mycolors) +
    ggtitle(
      paste0(
        "[Volcano Plot]", "\n",
        output_prefix, "\n",
        "padj < ", alpha, ", abs(log2FC) > ", foldchange, "\n"
      )) +
    theme_classic() +
    theme(
      aspect.ratio = 0.8,
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title = element_text(size = 8),
      text = element_text(size = 8)
    ) + coord_cartesian(clip = "off") +
    xlab("log2fc")

  fileName = here("figure", "volcano", paste0("volcano_", output_prefix, "_", alpha, "_", foldchange, "_shrlog2FC_", note))
  
  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot)
  dev.off()
  
  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot)
  dev.off()
  
  fileName = here("figure", "volcano", paste0("volcano_", output_prefix, "_", alpha, "_", foldchange, "_log2FC_", note))
  
  svglite(paste0(fileName, ".svg"),
          width = width,
          height = height)
  print(plot2)
  dev.off()
  
  png(
    paste0(fileName, ".png"),
    width = width,
    height = height,
    res = 600,
    units = "in"
  )
  print(plot2)
  dev.off()
  print(plot)
  print(plot2)  
  
  
}



getGO = function(sigListName, sigList, universe, categoryNum = 15, width = 7, height = 10){
  GO = enrichGO(gene = sigList, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP")
  fwrite(as.data.frame(GO), here("data", paste0("GO_", sigListName, ".tsv")), sep = "\t")
  dir.create(here("figure", "GO"), showWarnings = FALSE, recursive = TRUE)
  fileName <- paste0("GO_", sigListName)
  svglite(here("figure", "GO", paste0(fileName, ".svg")), height = height, width = width)
  print(clusterProfiler::dotplot(GO, showCategory = categoryNum) + scale_color_continuous(limits = c(0, 0.05),
                                                                low = "red", high = "black"))
  dev.off()
  png(here("figure", "GO", paste0(fileName, ".png")), width = width, height = height, res = 600, units = "in")
  print(clusterProfiler::dotplot(GO, showCategory = categoryNum) + scale_color_continuous(limits = c(0, 0.05),
                                                                         low = "red", high = "black"))
  dev.off()
  
  print(clusterProfiler::dotplot(GO, showCategory = categoryNum) + scale_color_continuous(limits = c(0, 0.05),
                                                                         low = "red", high = "black"))
}

