setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas")
{
  library(ArchR)
  library(ggalt)
  library(Seurat)
  library(BSgenome.OSativa.NCBI.IRGSPv1.0)
  library(parallel)
  library(clustree)
  library(dplyr)
  library(patchwork)
  library(universalmotif)
  library(TFBSTools)
  library(genomation)
  library(GenomicRanges)
  library(scales)
  library(ComplexHeatmap)
  library(pheatmap)
  library(Hmisc)
  library(Rmagic)
  library(AUCell)
  library(tidyverse)
  library(stringr)
  library(VennDiagram)
  library(RColorBrewer)
  library(grDevices)
  library(ggrepel)
  library(ggplot2)
  library(ggExtra)
  library(VennDiagram)
  library(ggpubr)
  
  genome_annotation <- createGenomeAnnotation(genome = BSgenome.OSativa.NCBI.IRGSPv1.0)
  blacklist_file <- "./Ref/irgsp1_repeat_unit.bed"
  bl <- read.table(blacklist_file)
  colnames(bl) <- c("chr", "start", "end")
  bl$strand <- "+"
  blacklist <- makeGRangesFromDataFrame(bl, keep.extra.columns = T)
  genome_annotation$blacklist <- blacklist
  gene_file <- "./Ref/rice_IRGSPv1.0.gene.txt"
  trans_file <- "./Ref/rice_IRGSPv1.0.transcript.txt"
  exon_file <- "./Ref/rice_IRGSPv1.0.exon.txt"
  archr_utils <- "./Ref/archr_utils.R"
  source(archr_utils)
  gene_annotation <- loadGeneinfo(gene = gene_file, exon = exon_file, trans = trans_file)
  gene_annotation$exons@elementMetadata$gene_id <- unlist(lapply(gene_annotation$exons@elementMetadata$gene_id,
                                                                 function(x){
                                                                   unlist(strsplit(x, "-"))[1]
                                                                 }))
  gene_annotation$exons@elementMetadata$symbol <- unlist(lapply(gene_annotation$exons@elementMetadata$symbol,
                                                                function(x){
                                                                  unlist(strsplit(x, "-"))[1]
                                                                }))
  gene_annotation$TSS@elementMetadata$tx_id <- unlist(lapply(gene_annotation$TSS@elementMetadata$tx_id,
                                                             function(x){
                                                               unlist(strsplit(x, "-"))[1]
                                                             }))
  gene_annotation$TSS@elementMetadata$tx_name <- unlist(lapply(gene_annotation$TSS@elementMetadata$tx_name,
                                                               function(x){
                                                                 unlist(strsplit(x, "-"))[1]
                                                               }))
  addArchRThreads(threads = 30)
  rscripts <- list.files(path = "./Ref/ArchR_R", pattern = ".R$", recursive = F, full.names = F)
  for (i in rscripts) {
    source(paste0("./Ref/ArchR_R/",i))
  }
  
  proj_pass_filter <- readRDS("3.proj_pass_filter_withPeaks.rds")
  
  celltype_color <- c("#faa818", # BM
                      "#d5a0e3", # Cortex
                      "#fbdf72", # cb/b
                      "#367d7d", # Endodermis
                      "#FF88C2", "#FF8888", # Epidermal
                      "#8ff0a4", # Fiber
                      "#5555FF", "#7744FF", "#B94FFF", # IM[1:3]
                      "#76dcb0", "#6ec9bb", "#5eb5c7", "#4c9fd4", "#3d8ede", # Mesophyll
                      "#df9cbf", # Pericycle
                      "#5c7ada", # Phloem
                      "#62589d", # Procambium
                      "#a4a5ee", # Rachis
                      "#c85c6c", # Root hair
                      "#b4446c", # Spikelet meristem(SM)
                      "#b7d62d" # Stele
  )
  names(celltype_color) <- sort(unique(proj_pass_filter$Celltype))
  
  cellnames <- getCellNames(proj_pass_filter)
  cellnames <- cellnames[which(proj_pass_filter$TissuesSub %in% c("10dRoot",
                                                                  "10dLeaf",
                                                                  "IM0.5cm",
                                                                  "IM1cm",
                                                                  "60dRootTip"))]
  ATAC_common_tissues <- proj_pass_filter[cellnames,]
  
  RNA <- readRDS("./Rice-snRNA/all_data_annotated.rds")
  RNA_common_tissues <- subset(RNA, subset = tissues %in% c("10dRoot",
                                                            "10dLeaf",
                                                            "IM0.5cm",
                                                            "IM1cm",
                                                            "60dRootTip"))
  unique(RNA_common_tissues$Celltype)
  RNA_common_tissues$Celltype[which(RNA_common_tissues$Celltype == "Branch meristems(BM)")] <- "Branch meristems (BM)"
  RNA_common_tissues$Celltype[which(RNA_common_tissues$Celltype == "Inflorescence meristem(IM)")] <- "Inflorescence meristem (IM)"
  RNA_common_tissues$Celltype[which(RNA_common_tissues$Celltype == "Spikelet meristem(SM)")] <- "Spikelet meristem (SM)"
  RNA_common_tissues$Celltype[which(RNA_common_tissues$Celltype == "Lemma(le)")] <- "Lemma (le)"
  RNA_common_tissues$Celltype[which(RNA_common_tissues$Celltype == "Cryptic bract/bract(cb/b)")] <- "Cryptic bract/bract (cb/b)"
  
  
  common_celltypes <- intersect(unique(RNA_common_tissues$Celltype),
                                unique(ATAC_common_tissues$Celltype))
  setdiff(unique(RNA_common_tissues$Celltype),
          unique(ATAC_common_tissues$Celltype))
  
  RNA_common_tissues <- subset(RNA_common_tissues, subset = Celltype %in% common_celltypes)
  cellnames <- getCellNames(ATAC_common_tissues)
  cellnames <- cellnames[which(ATAC_common_tissues$Celltype %in% common_celltypes)]
  ATAC_common_tissues <- ATAC_common_tissues[cellnames,]
  
  gene_id_map <- readRDS("gene_id_map.rds")
  rownames(gene_id_map) <- gene_id_map$gene_id
  RNA_id_gene <- data.frame(id = rownames(RNA_common_tissues),
                            symbol = gene_id_map[rownames(RNA_common_tissues),2],
                            row.names = rownames(RNA_common_tissues))
  RNA_id_gene <- as.data.frame(na.omit(RNA_id_gene))
  RNA_common_tissues <- RNA_common_tissues[RNA_id_gene$id,]
  RNA_common_tissues@assays$RNA@counts@Dimnames[[1]] <- RNA_id_gene$symbol
  RNA_common_tissues@assays$RNA@data@Dimnames[[1]] <- RNA_id_gene$symbol
  rownames(RNA_common_tissues@assays$RNA@scale.data) <- RNA_id_gene[rownames(RNA_common_tissues@assays$RNA@scale.data),"symbol"]
  RNA_common_tissues <- CreateSeuratObject(counts = RNA_common_tissues@assays$RNA@counts,
                                           meta.data = RNA_common_tissues@meta.data)
  RNA_common_tissues <- NormalizeData(object = RNA_common_tissues, verbose = FALSE)
  
  
  main_celltype <- data.frame(celltype = sort(unique(ATAC_common_tissues$Celltype)),
                              main = c("Meristem", "RootGroundTissue", "Meristem",
                                       "RootGroundTissue", "Epidermis", "Fiber",
                                       "Meristem", "Mesophyll", "Mesophyll",
                                       "Vasculature", "Vasculature", "Rachis",
                                       "Meristem"),
                              row.names = sort(unique(ATAC_common_tissues$Celltype)))
  ATAC_common_tissues$Mainclusters <- main_celltype[ATAC_common_tissues$Celltype, "main"]
  
  peaksets <- getPeakSet(ATAC_common_tissues)
  pr1 <- paste(as.character(seqnames(peaksets)),as.character(start(peaksets)),as.character(end(peaksets)),sep="_")
  peaksets <- as.data.frame(peaksets@elementMetadata)
  rownames(peaksets) <- pr1
  ###### celltype marker peaks
  {
    markersPeaks_Celltype <- getMarkerFeatures(
      ArchRProj = ATAC_common_tissues, 
      useMatrix = "PeakMatrix", 
      groupBy = "Celltype",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
    )
    
    markerPeaksList_Celltype <- getMarkers(markersPeaks_Celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
    
    markerPeaksDF <- c()
    for (i in names(markerPeaksList_Celltype@listData)) {
      temp <- as.data.frame(markerPeaksList_Celltype@listData[[i]])
      temp$Celltype <- i
      temp <- temp[order(temp$Log2FC, decreasing = T),]
      temp$Chr_start_end <- paste0(temp$seqnames,"_",temp$start,"_",temp$end)
      temp$nearestGene <- peaksets[temp$Chr_start_end, "nearestGene"]
      markerPeaksDF <- as.data.frame(rbind(markerPeaksDF,
                                           temp))
    }
    
  }
  
  RNA_common_tissues$Mainclusters <- main_celltype[RNA_common_tissues$Celltype, "main"]
  
}
########## celltype cluster
### ATAC
getAvailableMatrices(ATAC_common_tissues)
{
  GeneScoreMatrix <- getMatrixFromProject(ATAC_common_tissues,
                                          useMatrix = "GeneScoreMatrix")
  GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]@Dimnames[[1]] <- GeneScoreMatrix@elementMetadata@listData[["name"]]
  
  GeneScoreMatrix_mean <- aggregate.data.frame(t(GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]),
                                               by = list(GeneScoreMatrix@colData@listData[["Celltype"]]),
                                               FUN = mean)
  rownames(GeneScoreMatrix_mean) <- GeneScoreMatrix_mean$Group.1
  GeneScoreMatrix_mean <- GeneScoreMatrix_mean[,-1]
  GeneScoreMatrix_mean <- as.data.frame(t(GeneScoreMatrix_mean))
  GeneScoreMatrix_mean_cor <- rcorr(as.matrix(GeneScoreMatrix_mean))
  
  Maincluster <- c("Meristem", "RootGroundTissue", "Epidermis", "Fiber",
                   "Mesophyll", "Vasculature", "Rachis")
  Maincluster_color <- RColorBrewer::brewer.pal(n = 7, "Dark2")
  names(Maincluster_color) <- Maincluster
  cellanno_col <- columnAnnotation(df = data.frame(Maincluster = c("Meristem", "RootGroundTissue", "Meristem",
                                                                   "RootGroundTissue", "Epidermis", "Fiber",
                                                                   "Meristem", "Mesophyll", "Mesophyll",
                                                                   "Vasculature", "Vasculature", "Rachis",
                                                                   "Meristem")),
                                   col = list(Maincluster = Maincluster_color))
  cellanno_row <- rowAnnotation(df = data.frame(Maincluster = c("Meristem", "RootGroundTissue", "Meristem",
                                                                "RootGroundTissue", "Epidermis", "Fiber",
                                                                "Meristem", "Mesophyll", "Mesophyll",
                                                                "Vasculature", "Vasculature", "Rachis",
                                                                "Meristem")),
                                col = list(Maincluster = Maincluster_color))
  pdf("Figure6_4_ATAC_GeneScoreMatrix_celltype_cor.pdf", height = 6, width = 9)
  Heatmap(GeneScoreMatrix_mean_cor[["r"]],
          col = viridis::viridis(n = 100, option = "D"),
          name = "Pearson Coefficient\nGeneScoreMatrix", width = unit(80, "mm"),
          height = unit(80, "mm"), bottom_annotation = cellanno_col, right_annotation = cellanno_row)
  dev.off()
  rm(GeneScoreMatrix)
  
  
  ATAC_common_tissues <- addGeneIntegrationMatrix(
    ArchRProj = ATAC_common_tissues, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = RNA_common_tissues,
    addToArrow = TRUE,
    groupRNA = "Celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE,
    sampleCellsATAC = 15000,
    sampleCellsRNA = 15000
  )
  Gene_integrated <- getMatrixFromProject(ATAC_common_tissues,
                                          useMatrix = "GeneIntegrationMatrix")
  Gene_integrated@assays@data@listData[["GeneIntegrationMatrix"]]@Dimnames[[1]] <- Gene_integrated@elementMetadata@listData[["name"]]
  
  Gene_integrated_mean <- aggregate.data.frame(t(Gene_integrated@assays@data@listData[["GeneIntegrationMatrix"]]),
                                               by = list(Gene_integrated@colData@listData[["Celltype"]]),
                                               FUN = mean)
  rownames(Gene_integrated_mean) <- Gene_integrated_mean$Group.1
  Gene_integrated_mean <- Gene_integrated_mean[,-1]
  Gene_integrated_mean <- as.data.frame(t(Gene_integrated_mean))
  Gene_integrated_mean_cor <- rcorr(as.matrix(Gene_integrated_mean))
  
  pdf("Figure6_4_ATAC_GeneIntegrationMatrix_celltype_cor.pdf", height = 6, width = 9)
  Heatmap(Gene_integrated_mean_cor[["r"]],
          col = viridis::viridis(n = 100, option = "D"),
          name = "Pearson Coefficient\nIntegrationMatrix", width = unit(80, "mm"),
          height = unit(80, "mm"), bottom_annotation = cellanno_col, right_annotation = cellanno_row)
  dev.off()
  rm(Gene_integrated)
  
}
### RNA
{
  RNA_mean <- c()
  for (i in unique(RNA_common_tissues$Celltype)) {
    print(i)
    temp <- subset(RNA_common_tissues, subset = Celltype == i)
    temp_sum <- rowSums(temp@assays$RNA@data)
    temp_mean <- temp_sum / ncol(temp)
    temp_mean <- as.data.frame(temp_mean)
    colnames(temp_mean) <- i
    RNA_mean <- c(RNA_mean, list(temp_mean))
  }
  RNA_mean <- reduce(RNA_mean, cbind)
  RNA_cor <- rcorr(as.matrix(RNA_mean))
  colnames(RNA_cor$r)
  cellanno_col <- columnAnnotation(df = data.frame(Maincluster = c("Mesophyll", "Fiber", "Vasculature", "Mesophyll",
                                                                   "Vasculature", "Epidermis", "RootGroundTissue", "RootGroundTissue",
                                                                   "Meristem", "Meristem", "Meristem", "Rachis",
                                                                   "Meristem")),
                                   col = list(Maincluster = Maincluster_color))
  cellanno_row <- rowAnnotation(df = data.frame(Maincluster = c("Mesophyll", "Fiber", "Vasculature", "Mesophyll",
                                                                "Vasculature", "Epidermis", "RootGroundTissue", "RootGroundTissue",
                                                                "Meristem", "Meristem", "Meristem", "Rachis",
                                                                "Meristem")),
                                col = list(Maincluster = Maincluster_color))
  pdf("Figure6_4_RNA_celltype_cor.pdf", height = 6, width = 9)
  Heatmap(RNA_cor[["r"]],
          col = viridis::viridis(n = 100, option = "D"),
          name = "Pearson Coefficient\nGeneScoreMatrix", width = unit(80, "mm"),
          height = unit(80, "mm"), bottom_annotation = cellanno_col, right_annotation = cellanno_row)
  dev.off()
  rm(temp)
}
gc()

######### peak gene correlation
{
  peaksets <- getPeakSet(ATAC_common_tissues)
  ranges <- as.data.frame(peaksets@ranges)
  seqname <- as.data.frame(peaksets@seqnames)
  peaksets <- as.data.frame(peaksets@elementMetadata@listData)
  peaksets <- as.data.frame(cbind(peaksets, ranges))
  peaksets$Chr <- seqname$value
  peaksets$ID <- paste0(peaksets$Chr, "_", peaksets$start, "_", peaksets$end)
  rownames(peaksets) <- peaksets$ID
  
  RNA_mean
  peak_matrix <- getMatrixFromProject(ATAC_common_tissues,
                                      useMatrix = "PeakMatrix")
  peak_matrix_meta <- as.data.frame(peak_matrix@rowRanges)
  peak_matrix_meta$ID <- paste0(peak_matrix_meta$seqnames, "_", peak_matrix_meta$start, "_", peak_matrix_meta$end)
  rownames(peak_matrix_meta) <- peak_matrix_meta$ID
  identical(peaksets$ID, peak_matrix_meta$ID) # TRUE
  peak_matrix@assays@data@listData[["PeakMatrix"]]@Dimnames[[1]] <- peak_matrix_meta$ID
  peak_matrix_log <- peak_matrix
  peak_matrix_log <- t(peak_matrix@assays@data@listData[["PeakMatrix"]])/colSums(peak_matrix@assays@data@listData[["PeakMatrix"]])
  rowSums(peak_matrix_log)
  peak_matrix_log <- log2(t(peak_matrix_log) * 10^4 + 1)
  
  identical(colnames(peak_matrix), colnames(peak_matrix_log))
  Peak_mean <- c()
  for (i in unique(peak_matrix$Celltype)) {
    print(i)
    temp <- peak_matrix_log[,which(peak_matrix@colData@listData[["Celltype"]] == i)]
    temp_sum <- rowSums(temp)
    temp_mean <- temp_sum / ncol(temp)
    temp_mean <- as.data.frame(temp_mean)
    colnames(temp_mean) <- i
    Peak_mean <- c(Peak_mean, list(temp_mean))
  }
  Peak_mean <- reduce(Peak_mean, cbind)
  Peak_mean
  
  RNA_mean <- RNA_mean[, sort(colnames(RNA_mean))]
  RNA_mean <- t(RNA_mean)
  Peak_mean <- Peak_mean[, sort(colnames(Peak_mean))]
  Peak_mean <- t(Peak_mean)
  
  
  Peak_gene_pearson_cor <- cor(RNA_mean,Peak_mean, method = "pearson")
  Peak_gene_pearson_cor[1:5,1:5]
  Peak_gene_spearman_cor <- cor(RNA_mean,Peak_mean, method = "spearman")
  Peak_gene_spearman_cor[1:5,1:5]
  gc()
  
  peaksets_2 <- peaksets[which(peaksets$nearestGene %in% colnames(RNA_mean)),]
  pearson <- c()
  spearman <- c()
  for (i in 1:nrow(peaksets_2)) {
    print(i)
    pearson <- c(pearson, Peak_gene_pearson_cor[peaksets_2[i,"nearestGene"],
                                                peaksets_2[i,"ID"]])
    spearman <- c(spearman, Peak_gene_spearman_cor[peaksets_2[i,"nearestGene"],
                                                   peaksets_2[i,"ID"]])
  }
  peaksets_2$pearson <- pearson
  peaksets_2$spearman <- spearman
  peaksets_2 <- peaksets_2[!is.na(peaksets_2$pearson),]
  pdf("Figure6_5_PG_cor_pearson_distribution.pdf", width = 5.5, height = 5.5)
  p1 <- ggplot(data = peaksets_2, aes(x = pearson)) +
    geom_density(linewidth = 0.8, color = "black", fill = "#6495ED") +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(size =10, color = "black"),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position = "none") +
    # geom_vline(xintercept = quantile(peaksets_2$pearson, 0.95),
    #            color = "gray20") +
    geom_vline(xintercept = 0.6,
               color = "gray20") +
    annotate(geom = "text", x = 0.6, y = 1,
             label = paste0("Threshold = ", 0.6,
                            "\n", sum(peaksets_2$pearson > 0.6), " Peak-Gene pairs",
                            "\n", length(unique(peaksets_2[which(peaksets_2$pearson > 0.6), "nearestGene"])), " genes"),
             vjust = 1, hjust = 0.5, size = 3.3) +
    # annotate(geom = "text", x = quantile(peaksets_2$pearson, 0.95), y = 1,
    #          label = paste0("Threshold = ", round(quantile(peaksets_2$pearson, 0.95), 3),
    #                         "\n", "95th percentile", "\n", sum(peaksets_2$pearson > quantile(peaksets_2$pearson, 0.95)), " Peak-Gene pairs",
    #                         "\n", length(unique(peaksets_2[which(peaksets_2$pearson > quantile(peaksets_2$pearson, 0.95)), "nearestGene"])), " genes"),
    #          vjust = 1, hjust = 0.5, size = 3.3) +
    labs(x = "Pearson Coefficient", y = "Density")
  p1
  dev.off()
  
  pdf("Figure6_5_PG_cor_spearman_distribution.pdf", width = 5.5, height = 5.5)
  p1 <- ggplot(data = peaksets_2, aes(x = spearman)) +
    geom_density(linewidth = 0.8, color = "black", fill = "#6495ED") +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(size =10, color = "black"),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position = "none") +
    # geom_vline(xintercept = quantile(peaksets_2$spearman, 0.95),
    #            color = "gray20") +
    geom_vline(xintercept = 0.6,
               color = "gray20") +
    annotate(geom = "text", x = 0.6, y = 1,
             label = paste0("Threshold = ", 0.6,
                            "\n", sum(peaksets_2$spearman > 0.6), " Peak-Gene pairs",
                            "\n", length(unique(peaksets_2[which(peaksets_2$spearman > 0.6), "nearestGene"])), " genes"),
             vjust = 1, hjust = 0.5, size = 3.3) +
    # annotate(geom = "text", x = quantile(peaksets_2$spearman, 0.95), y = 1,
    #          label = paste0("Threshold = ", round(quantile(peaksets_2$spearman, 0.95), 3),
    #                         "\n", "95th percentile", "\n", sum(peaksets_2$spearman > quantile(peaksets_2$spearman, 0.95)), " Peak-Gene pairs",
    #                         "\n", length(unique(peaksets_2[which(peaksets_2$spearman > quantile(peaksets_2$spearman, 0.95)), "nearestGene"])), " genes"),
    #          vjust = 1, hjust = 0.5, size = 3.3) +
    labs(x = "Spearman Coefficient", y = "Density")
  p1
  dev.off()
  
  cor_threshold <- 0.6
  pearson_spearman_common_peak <- peaksets_2[which(peaksets_2$pearson > cor_threshold &
                                                   peaksets_2$spearman > cor_threshold), "ID"]
  pearson_spearman_common_gene <- unique(peaksets_2[which(peaksets_2$pearson > cor_threshold &
                                                     peaksets_2$spearman > cor_threshold), "nearestGene"])
  pearson_gene <- unique(peaksets_2[which(peaksets_2$pearson > cor_threshold), "nearestGene"])
  spearman_gene <- unique(peaksets_2[which(peaksets_2$spearman > cor_threshold), "nearestGene"])

  
  gene_linked_peaks <- peaksets_2[which(peaksets_2$pearson > cor_threshold), "ID"]
  Not_gene_linked_peaks <- peaksets_2[which(peaksets_2$pearson <= cor_threshold), "ID"]
  
}

#### RNA DEGs
{
  Idents(RNA_common_tissues) <- RNA_common_tissues$Celltype
  celltype_DEGs <- FindAllMarkers(RNA_common_tissues, only.pos = T, logfc.threshold = 0.25)
  celltype_DEGs <- celltype_DEGs[which(celltype_DEGs$p_val_adj < 0.05),]
  length(unique(celltype_DEGs$gene)) # 6140
  peaksets_2$DEGs <- ifelse(peaksets_2$nearestGene %in% celltype_DEGs$gene,
                            "YES", "NO")
  table(peaksets_2$DEGs)
  
  temp <- peaksets_2[order(peaksets_2$pearson,decreasing = F),]
  temp <- data.frame(row.names = 1:nrow(temp),
                     Type = temp$DEGs)
  pdf("Figure6_4_peak_DEG.pdf", width = 4.5, height = 2)
  Heatmap(t(temp), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
          col = c("white", "red"))
  dev.off()
  
  # sum(peaksets_2[which(peaksets_2$DEGs == "YES"), "ID"] %in% Global_Peaks)
  # length(Global_Peaks)
  # sum(peaksets_2[which(peaksets_2$DEGs == "YES"), "ID"] %in% Lineage_Peaks)
  # length(Lineage_Peaks)
  # sum(peaksets_2[which(peaksets_2$DEGs == "YES"), "ID"] %in% Other_Peaks)
  # length(Other_Peaks)
  # sum(peaksets_2[which(peaksets_2$DEGs == "YES"), "ID"] %in% Gene_Linked_Peaks)
  # length(Other_Peaks)
  
  Idents(RNA_common_tissues) <- RNA_common_tissues$Mainclusters
  lineage_DEGs <- FindAllMarkers(RNA_common_tissues, only.pos = T, logfc.threshold = 0.25)
}

######### global constitutive peaks
#### average insertion count rank
{
  # {
  #   Peak_mean <- t(Peak_mean)
  #   Peak_mean_2 <- Peak_mean[Not_gene_linked_peaks,]
  #   Peak_mean_melt <- c()
  #   for (i in colnames(Peak_mean_2)) {
  #     print(i)
  #     temp <- as.data.frame(Peak_mean_2[,i])
  #     colnames(temp) <- "Value"
  #     temp$Celltype <- i
  #     temp <- temp[order(temp$Value, decreasing = T),]
  #     temp$Rank <- 1:nrow(temp)
  #     temp$PeakType <- "Other Peaks"
  #     temp[1:ceiling(nrow(temp)*0.05), "PeakType"] <- "Celltype Constitutive Peaks"
  #     temp$PeakID <- rownames(temp)
  #     Peak_mean_melt <- as.data.frame(rbind(Peak_mean_melt,
  #                                           temp))
  #   }
  #   pdf("Figure6_3_Plan_A_celltype_constitutive_peaks.pdf", width = 12, height = 10)
  #   ggplot(data = Peak_mean_melt, aes(x = Rank, y = Value, color = PeakType)) +
  #     geom_point(size = 0.3) +
  #     theme_bw() +
  #     theme(panel.grid.major.y = element_blank(),
  #           panel.grid.minor.y = element_blank(),
  #           axis.text = element_text(size =10, color = "black"),
  #           axis.title = element_text(size = 10),
  #           legend.text = element_text(size = 10),
  #           legend.title = element_text(size = 10)) +
  #     labs(y = "Average insertion count", x = "Rank") +
  #     facet_wrap(.~Celltype, )
  #   dev.off()
  #   
  #   table(Peak_mean_melt$PeakType)
  #   Peak_mean_melt$Binary <- ifelse(Peak_mean_melt$PeakType == "Celltype Constitutive Peaks",
  #                                   1, 0)
  #   Celltype_index <- data.frame(celltype = sort(unique(Peak_mean_melt$Celltype)),
  #                                index = 1:13,
  #                                row.names = sort(unique(Peak_mean_melt$Celltype)))
  #   Peaks_index <- data.frame(peak = rownames(Peak_mean_2),
  #                             index = 1:nrow(Peak_mean_2),
  #                             row.names = rownames(Peak_mean_2))
  #   PeakIndex <- Peaks_index[Peak_mean_melt$PeakID, "index"]
  #   CelltypeIndex <- Celltype_index[Peak_mean_melt$Celltype, "index"]
  #   
  #   Peak_mean_melt$PeakIndex <- PeakIndex
  #   Peak_mean_melt$CelltypeIndex <- CelltypeIndex
  #   Peak_mean_melt <- Peak_mean_melt[order(Peak_mean_melt$Celltype,
  #                                          Peak_mean_melt$PeakIndex),]
  #   nrow(Peak_mean_melt) == nrow(Peaks_index) * nrow(Celltype_index)
  #   Peak_type_binary <- Peak_mean_melt$Binary
  #   dim(Peak_type_binary) <- c(nrow(Peaks_index), nrow(Celltype_index))
  #   rownames(Peak_type_binary) <- Peak_mean_melt$PeakID[1:nrow(Peaks_index)]
  #   colnames(Peak_type_binary) <- unique(Peak_mean_melt$Celltype)
  #   
  #   table(rowSums(Peak_type_binary))
  #   
  #   Other_Peak_Matrix <- Peak_type_binary[which(rowSums(Peak_type_binary) == 0),]
  #   nrow(Other_Peak_Matrix)
  #   Global_Peak_Matrix <- Peak_type_binary[which(rowSums(Peak_type_binary) == 13),]
  #   nrow(Global_Peak_Matrix)
  #   Cell_Peak_Matrix <- Peak_type_binary[which(rowSums(Peak_type_binary) %in% c(1:12)),]
  #   nrow(Cell_Peak_Matrix)
  #   
  #   nrow(Peak_type_binary) == nrow(Other_Peak_Matrix)+nrow(Global_Peak_Matrix)+nrow(Cell_Peak_Matrix)
  #   
  #   
  #   ht <- pheatmap(Cell_Peak_Matrix, show_rownames = F, show_colnames = T)
  #   ht
  #   # rownames(Cell_Peak_Matrix[ht$tree_row[["order"]],])
  #   peak_order <- rownames(Cell_Peak_Matrix[ht$tree_row[["order"]],])
  #   Cell_Peak_Matrix <- Cell_Peak_Matrix[peak_order,c("Spikelet meristem (SM)", "Cryptic bract/bract (cb/b)",
  #                                           "Branch meristems (BM)", "Inflorescence meristem (IM)",
  #                                           "Rachis", "Mesophyll precursor", "Mesophyll (MO)",
  #                                           "Endodermis", "Cortex", "Phloem", "Procambium", "Fiber",
  #                                           "Epidermis")]
  #   Heatmap(Cell_Peak_Matrix, show_row_names = F, show_column_names = T,
  #           cluster_rows = F, cluster_columns = F, raster_by_magick = F,
  #           show_row_dend = F)
  #   
  #   Peak_type_binary <- rbind(Global_Peak_Matrix[,colnames(Cell_Peak_Matrix)],
  #                             Cell_Peak_Matrix)
  #   ht1 <- Heatmap(Global_Peak_Matrix[,colnames(Cell_Peak_Matrix)], show_row_names = F, show_column_names = T,
  #                  cluster_rows = F, cluster_columns = F, raster_by_magick = F,
  #                  show_row_dend = F)
  #   cellanno_col <- columnAnnotation(df = data.frame(Maincluster = c("Meristem", "Meristem", "Meristem", "Meristem",
  #                                                                    "Rachis", "Mesophyll", "Mesophyll", "RootGroundTissue",
  #                                                                    "RootGroundTissue", "Vasculature", "Vasculature",
  #                                                                    "Fiber", "Epidermis")),
  #                                    col = list(Maincluster = Maincluster_color))
  #   ht2 <- Heatmap(Cell_Peak_Matrix, show_row_names = F, show_column_names = T,
  #                  cluster_rows = F, cluster_columns = F, raster_by_magick = F,
  #                  show_row_dend = F, bottom_annotation = cellanno_col,
  #                  col = c("#00acc1", "#e53935"))
  #   ht3 <- Heatmap(Other_Peak_Matrix[,colnames(Cell_Peak_Matrix)], show_row_names = F, show_column_names = T,
  #                  cluster_rows = F, cluster_columns = F, raster_by_magick = F,
  #                  show_row_dend = F)
  #   
  #   pdf("Figure6_3_Plan_A_celltype_constitutive_peaks_heatmap.pdf", width = 4.5, height = 6)
  #   ht2
  #   dev.off()
  #   
  #   PlanA_Other_Peaks <- rownames(Other_Peak_Matrix)
  #   PlanA_Celltype_Peaks <- rownames(Cell_Peak_Matrix)
  #   PlanA_Global_Peaks <- rownames(Global_Peak_Matrix)
  # }
  }
#### accessible in > 5% cells & Not celltype marker peaks
{
  
  length(unique(markerPeaksDF$Chr_start_end)) # 43689 celltype marker peaks
  Peak_percent_celltype <- c()
  for (i in unique(peak_matrix@colData@listData[["Celltype"]])) {
    print(i)
    temp <- peak_matrix[,which(peak_matrix@colData@listData[["Celltype"]] == i)]
    temp <- temp@assays@data@listData[["PeakMatrix"]]
    temp <- temp[Not_gene_linked_peaks,]
    temp <- temp > 0
    temp <- temp + 0
    ncell <- rowSums(temp)
    percent <- ncell / ncol(temp)
    temp <- data.frame(PeakID = rownames(temp),
                       Celltype = i,
                       Percent = percent
    )
    temp <- temp[order(temp$Percent, decreasing = T),]
    temp$Rank <- 1:nrow(temp)
    temp$PeakType <- "Other Peaks"
    temp[which(temp$Percent > 0.05), "PeakType"] <- "Celltype Constitutive Peaks"
    
    temp2 <- markerPeaksDF[which(markerPeaksDF$Celltype == i),]
    
    temp$MarkerPeaks <- "Not Cell-type Marker Peak"
    temp$MarkerPeaks[which(temp$PeakID %in% temp2$Chr_start_end)] <- "Cell-type Marker Peak"
    
    Peak_percent_celltype <- as.data.frame(rbind(Peak_percent_celltype,
                                                 temp))
  }
  pdf("Figure6_5_celltype_constitutive_peaks.pdf", width = 12, height = 10)
  ggplot(data = Peak_percent_celltype, aes(x = Rank, y = Percent, color = PeakType)) +
    geom_point(size = 0.3) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(size =10, color = "black"),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Percentage of cells with accessible peak in cell-type", x = "Peak rank") +
    facet_wrap(.~Celltype)
  dev.off()
  
  temp <- Peak_percent_celltype[which(Peak_percent_celltype$PeakType == "Celltype Constitutive Peaks"),]
  temp_2 <- as.data.frame(table(temp$Celltype, temp$MarkerPeaks))
  colnames(temp_2) <- c("Celltype", "MarkerPeaks", "Count")
  pdf("Figure6_5_celltype_accessible_marker_peak_count.pdf", width = 6, height = 8)
  ggplot() +
    geom_bar(data = temp_2,
             aes(x = MarkerPeaks, y = Count, fill = MarkerPeaks),
             stat = "identity") +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(size =10, color = "black"),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    labs(y = "Count", x = "") +
    geom_text(data = temp_2, aes(x = MarkerPeaks, y = Count, label = Count),
              vjust = -0.5, size = 3.3) +
    scale_fill_manual(values = c("#FF6347", "#4682B4")) +
    scale_y_continuous(limits = c(0, 32000)) +
    facet_wrap(.~Celltype)
  dev.off()
  
  length(unique(Peak_percent_celltype$PeakID)) == length(Not_gene_linked_peaks)
  table(Peak_percent_celltype$MarkerPeaks)
  Peak_percent_celltype <- Peak_percent_celltype[which(Peak_percent_celltype$MarkerPeaks == "Not Cell-type Marker Peak"),]
  
  unique(Peak_percent_celltype$PeakType)
  Peak_percent_celltype$Binary <- ifelse(Peak_percent_celltype$PeakType == "Celltype Constitutive Peaks",
                                         1, 0)
  Celltype_index <- data.frame(celltype = sort(unique(Peak_percent_celltype$Celltype)),
                               index = 1:13,
                               row.names = sort(unique(Peak_percent_celltype$Celltype)))
  Peaks_index <- data.frame(peak = unique(Peak_percent_celltype$PeakID),
                            index = 1:length(unique(Peak_percent_celltype$PeakID)),
                            row.names = unique(Peak_percent_celltype$PeakID))
  PeakIndex <- Peaks_index[Peak_percent_celltype$PeakID, "index"]
  CelltypeIndex <- Celltype_index[Peak_percent_celltype$Celltype, "index"]
  
  Peak_percent_celltype$PeakIndex <- PeakIndex
  max(Peak_percent_celltype$PeakIndex)
  Peak_percent_celltype$CelltypeIndex <- CelltypeIndex
  nrow(Peak_percent_celltype) == nrow(Peaks_index) * nrow(Celltype_index)
  Peak_percent_celltype_2 <- Peak_percent_celltype[which(Peak_percent_celltype$Binary == 1),]
  maxindex <- max(Peak_percent_celltype_2$PeakIndex)
  nrow(Peaks_index)
  removed_peaks <- Peaks_index[which(Peaks_index$index %in% c((maxindex+1):nrow(Peaks_index))),"peak"]
  Peak_percent_celltype[which(Peak_percent_celltype$PeakID %in% removed_peaks),]
  Peak_type_binary <- matrix(0, nrow = maxindex, ncol = length(unique(Celltype_index$celltype)))
  colnames(Peak_type_binary) <- unique(Celltype_index$celltype)
  rownames(Peak_type_binary) <- Peaks_index$peak[1:maxindex]
  nrow(Peak_percent_celltype_2)
  for (i in 1:nrow(Peak_percent_celltype_2)) {
    print(i)
    Peak_type_binary[Peak_percent_celltype_2[i, "PeakID"],
                     Peak_percent_celltype_2[i, "Celltype"]] <- 1
  }
  
  # Peak_type_binary <- Peak_percent_celltype_2$Binary
  # Peak_type_binary <- sparseMatrix(i = Peak_percent_celltype_2$PeakIndex, j = Peak_percent_celltype_2$CelltypeIndex, x = Peak_type_binary)
  # nrow(Peak_type_binary)
  # Peak_type_binary <- as.matrix(Peak_type_binary)
  # rownames(Peak_type_binary) <- Peaks_index$peak[1:nrow(Peak_type_binary)]
  # colnames(Peak_type_binary) <- Celltype_index$celltype
  
  
  table(rowSums(Peak_type_binary))
  
  Other_Peak_Matrix <- Peak_type_binary[which(rowSums(Peak_type_binary) == 0),]
  nrow(Other_Peak_Matrix)
  Global_Peak_Matrix <- Peak_type_binary[which(rowSums(Peak_type_binary) == 13),]
  nrow(Global_Peak_Matrix)
  Cell_Peak_Matrix <- Peak_type_binary[which(rowSums(Peak_type_binary) %in% c(1:12)),]
  nrow(Cell_Peak_Matrix)
  
  nrow(Peak_type_binary) == nrow(Other_Peak_Matrix) + nrow(Global_Peak_Matrix) + nrow(Cell_Peak_Matrix)
  
  
  ht <- pheatmap(Cell_Peak_Matrix, show_rownames = F, show_colnames = T)
  # ht
  peak_order <- rownames(Cell_Peak_Matrix[ht$tree_row[["order"]],])
  Cell_Peak_Matrix <- Cell_Peak_Matrix[peak_order,c("Spikelet meristem (SM)", "Cryptic bract/bract (cb/b)",
                                                    "Branch meristems (BM)", "Inflorescence meristem (IM)",
                                                    "Rachis", "Mesophyll precursor", "Mesophyll (MO)",
                                                    "Endodermis", "Cortex", "Phloem", "Procambium", "Fiber",
                                                    "Epidermis")]
  # 
  # ht1 <- Heatmap(Global_Peak_Matrix[,colnames(Cell_Peak_Matrix)], show_row_names = F, show_column_names = T,
  #                cluster_rows = F, cluster_columns = F, raster_by_magick = F,
  #                show_row_dend = F)
  cellanno_col <- columnAnnotation(df = data.frame(Maincluster = c("Meristem", "Meristem", "Meristem", "Meristem",
                                                                   "Rachis", "Mesophyll", "Mesophyll", "RootGroundTissue",
                                                                   "RootGroundTissue", "Vasculature", "Vasculature",
                                                                   "Fiber", "Epidermis")),
                                   col = list(Maincluster = Maincluster_color))
  ht2 <- Heatmap(Cell_Peak_Matrix, show_row_names = F, show_column_names = T,
                 cluster_rows = F, cluster_columns = F, raster_by_magick = F,
                 show_row_dend = F, bottom_annotation = cellanno_col,
                 col = c("#00acc1", "#e53935"))
  # ht3 <- Heatmap(Other_Peak_Matrix[,colnames(Cell_Peak_Matrix)], show_row_names = F, show_column_names = T,
  #                cluster_rows = F, cluster_columns = F, raster_by_magick = F,
  #                show_row_dend = F)
  
  pdf("Figure6_5_celltype_constitutive_peaks_heatmap.pdf", width = 4.5, height = 6)
  ht2
  dev.off()
  
  PlanB_Other_Peaks <- rownames(Other_Peak_Matrix)
  PlanB_Other_Peaks <- c(PlanB_Other_Peaks, removed_peaks)
  PlanB_Celltype_Peaks <- rownames(Cell_Peak_Matrix)
  PlanB_Global_Peaks <- rownames(Global_Peak_Matrix)
  
  length(c(PlanB_Other_Peaks, PlanB_Celltype_Peaks, PlanB_Global_Peaks)) == length(Not_gene_linked_peaks)
  
  temp_peaks <- intersect(PlanB_Global_Peaks, unique(markerPeaksDF$Chr_start_end)) # 0
  
  PlanB_Global_Peaks_genes <- peaksets_2[PlanB_Global_Peaks, "nearestGene"]
  PlanB_Global_Peaks_genes <- unique(PlanB_Global_Peaks_genes)
}

######### lineage constitutive peaks
{
  Lineage_Epidermis <- c("Epidermis")
  Lineage_Fiber <- c("Fiber")
  Lineage_Meristem <- c("Inflorescence meristem (IM)",
                        "Cryptic bract/bract (cb/b)",
                        "Branch meristems (BM)",
                        "Spikelet meristem (SM)")
  Lineage_Rachis <- c("Rachis")
  Lineage_RootGroundTissue <- c("Cortex", "Endodermis")
  Lineage_Vasculature <- c("Phloem", "Procambium")
  Lineage_Mesophyll <- c("Mesophyll (MO)", "Mesophyll precursor")
  # Lineage_Epidermis
  Lineage_Epidermis_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_Epidermis])
  colnames(Lineage_Epidermis_peaks_matrix) <- "Epidermis"
  Lineage_Epidermis_peaks <- rownames(Lineage_Epidermis_peaks_matrix)
  Lineage_Epidermis_peaks <- Lineage_Epidermis_peaks[which(rowSums(Lineage_Epidermis_peaks_matrix) == length(Lineage_Epidermis))]
  # Lineage_Fiber
  Lineage_Fiber_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_Fiber])
  colnames(Lineage_Fiber_peaks_matrix) <- "Fiber"
  Lineage_Fiber_peaks <- rownames(Lineage_Fiber_peaks_matrix)
  Lineage_Fiber_peaks <- Lineage_Fiber_peaks[which(rowSums(Lineage_Fiber_peaks_matrix) == length(Lineage_Fiber))]
  # Lineage_Meristem
  Lineage_Meristem_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_Meristem])
  Lineage_Meristem_peaks <- rownames(Lineage_Meristem_peaks_matrix)
  Lineage_Meristem_peaks <- Lineage_Meristem_peaks[which(rowSums(Lineage_Meristem_peaks_matrix) == length(Lineage_Meristem))]
  # Lineage_Rachis
  Lineage_Rachis_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_Rachis])
  colnames(Lineage_Rachis_peaks_matrix) <- "Rachis"
  Lineage_Rachis_peaks <- rownames(Lineage_Rachis_peaks_matrix)
  Lineage_Rachis_peaks <- Lineage_Rachis_peaks[which(rowSums(Lineage_Rachis_peaks_matrix) == length(Lineage_Rachis))]
  # Lineage_RootGroundTissue
  Lineage_RootGroundTissue_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_RootGroundTissue])
  Lineage_RootGroundTissue_peaks <- rownames(Lineage_RootGroundTissue_peaks_matrix)
  Lineage_RootGroundTissue_peaks <- Lineage_RootGroundTissue_peaks[which(rowSums(Lineage_RootGroundTissue_peaks_matrix) == length(Lineage_RootGroundTissue))]
  # Lineage_Vasculature
  Lineage_Vasculature_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_Vasculature])
  Lineage_Vasculature_peaks <- rownames(Lineage_Vasculature_peaks_matrix)
  Lineage_Vasculature_peaks <- Lineage_Vasculature_peaks[which(rowSums(Lineage_Vasculature_peaks_matrix) == length(Lineage_Vasculature))]
  # Lineage_Mesophyll
  Lineage_Mesophyll_peaks_matrix <- as.data.frame(Cell_Peak_Matrix[,Lineage_Mesophyll])
  Lineage_Mesophyll_peaks <- rownames(Lineage_Mesophyll_peaks_matrix)
  Lineage_Mesophyll_peaks <- Lineage_Mesophyll_peaks[which(rowSums(Lineage_Mesophyll_peaks_matrix) == length(Lineage_Mesophyll))]
  
  library(VennDiagram)
  library(venn)
  
  Lineage_color <- c("#7570B3", "#E7298A", "#1B9E77",
                     "#A6761D", "#D95F02", "#E6AB02", "#66A61E")
  names(Lineage_color) <- c("Epidermis", "Fiber", "Meristem", "Rachis",
                            "RootGroundTissue", "Vasculature", "Mesophyll")
  venn_list <- list(Epidermis = Lineage_Epidermis_peaks,
                    Fiber = Lineage_Fiber_peaks,
                    Meristem = Lineage_Meristem_peaks,
                    Rachis = Lineage_Rachis_peaks,
                    RootGroundTissue = Lineage_RootGroundTissue_peaks,
                    Vasculature = Lineage_Vasculature_peaks,
                    Mesophyll = Lineage_Mesophyll_peaks)
  pdf("Figure6_5_Lineage_constitutive_peaks_venn", width = 7, height = 7)
  venn(venn_list,
       zcolor = c("#7570B3", "#E7298A", "#1B9E77",
                  "#A6761D", "#D95F02", "#E6AB02", "#66A61E"), # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
       opacity = 0.6,  # 调整颜色透明度
       box = F,        # 是否添加边框
       ilcs = 1,     # 数字大小
       ellipse = T,
       sncs = 1.2        # 组名字体大小
  )
  dev.off()
  
  inter <- get.venn.partitions(venn_list)
  for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
  inter <- subset(inter, select = -..values.. )
  inter <- subset(inter, select = -..set.. )
  write.table(inter, "Figure6_5_Lineage_constitutive_peaks_venn_result.csv", row.names = FALSE, sep = ',', quote = FALSE)
  
  PlanB_Lineage_peaks <- unique(c(Lineage_Epidermis_peaks, Lineage_Fiber_peaks,
                                  Lineage_Meristem_peaks, Lineage_Rachis_peaks,
                                  Lineage_RootGroundTissue_peaks,
                                  Lineage_Vasculature_peaks, Lineage_Mesophyll_peaks))
  sum(PlanB_Lineage_peaks %in% PlanB_Celltype_Peaks) == length(PlanB_Lineage_peaks)
  PlanB_Other_Peaks
  PlanB_Celltype_Peaks
  PlanB_Global_Peaks
  length(setdiff(PlanB_Celltype_Peaks, PlanB_Lineage_peaks)) == length(PlanB_Celltype_Peaks) - length(PlanB_Lineage_peaks)

}

######### 4 peak type count
{
  Other_Peaks <- c(PlanB_Other_Peaks, 
                   setdiff(PlanB_Celltype_Peaks, PlanB_Lineage_peaks))
  length(Other_Peaks)
  
  Lineage_Peaks <- PlanB_Lineage_peaks
  length(Lineage_Peaks)
  
  Global_Peaks <- PlanB_Global_Peaks
  length(Global_Peaks)
  
  Gene_Linked_Peaks <- gene_linked_peaks
  length(Gene_Linked_Peaks)
  
  length(Other_Peaks) +
    length(Global_Peaks) +
    length(Lineage_Peaks) +
    length(Gene_Linked_Peaks) == nrow(peaksets_2)
  
  peak_group <- data.frame(peakType = c(rep("Gene-linked peaks", length(Gene_Linked_Peaks)),
                                        rep("Lineage constitutive peaks", length(Lineage_Peaks)),
                                        rep("Global constitutive peaks", length(Global_Peaks)),
                                        rep("Other peaks", length(Other_Peaks))),
                           peakID = c(Gene_Linked_Peaks,
                                      Lineage_Peaks,
                                      Global_Peaks,
                                      Other_Peaks))
  rownames(peak_group) <- peak_group$peakID
  saveRDS(peak_group, "Figure6_peak_group.rds")
  
  Gene_Linked_Peaks_DEGs <- Gene_Linked_Peaks[which(Gene_Linked_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "YES")])]
  Gene_Linked_Peaks_DEGs_not <- Gene_Linked_Peaks[which(Gene_Linked_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "NO")])]
  Lineage_Peaks_DEGs <- Lineage_Peaks[which(Lineage_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "YES")])]
  Lineage_Peaks_DEGs_not <- Lineage_Peaks[which(Lineage_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "NO")])]
  Global_Peaks_DEGs <- Global_Peaks[which(Global_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "YES")])]
  Global_Peaks_DEGs_not <- Global_Peaks[which(Global_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "NO")])]
  Other_Peaks_DEGs <- Other_Peaks[which(Other_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "YES")])]
  Other_Peaks_DEGs_not <- Other_Peaks[which(Other_Peaks %in% peaksets_2$ID[which(peaksets_2$DEGs == "NO")])]
  
  
  temp <- data.frame(Peak_Type = c("Gene-linked peaks",
                                   "Gene-linked peaks",
                                   "Lineage constitutive peaks",
                                   "Lineage constitutive peaks",
                                   "Global constitutive peaks",
                                   "Global constitutive peaks",
                                   "Other peaks",
                                   "Other peaks"),
                     Count = c(length(Gene_Linked_Peaks_DEGs),
                               length(Gene_Linked_Peaks_DEGs_not),
                               length(Lineage_Peaks_DEGs),
                               length(Lineage_Peaks_DEGs_not),
                               length(Global_Peaks_DEGs),
                               length(Global_Peaks_DEGs_not),
                               length(Other_Peaks_DEGs),
                               length(Other_Peaks_DEGs_not)),
                     Group = c("Nearest with DEGs", "Nearest with other genes",
                               "Nearest with DEGs", "Nearest with other genes",
                               "Nearest with DEGs", "Nearest with other genes",
                               "Nearest with DEGs", "Nearest with other genes"
                               ))
  
  temp$Peak_Type <- factor(temp$Peak_Type,
                           levels = c("Gene-linked peaks",
                                          "Lineage constitutive peaks",
                                          "Global constitutive peaks",
                                          "Other peaks"))
  temp$Group <- factor(temp$Group,
                       levels = c("Nearest with DEGs",
                                  "Nearest with other genes"))
  
  peak_type_color <- c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD")
  names(peak_type_color) <- c("Gene-linked peaks",
                              "Lineage constitutive peaks",
                              "Global constitutive peaks",
                              "Other peaks")
  
  temp <- aggregate.data.frame(temp$Count, by = list(temp$Peak_Type), FUN = sum)
  colnames(temp) <- c("Peak_Type", "Count")
  temp$Peak_Type <- factor(temp$Peak_Type,
                           levels = c("Gene-linked peaks",
                                      "Lineage constitutive peaks",
                                      "Global constitutive peaks",
                                      "Other peaks"))
  pdf("Figure6_5_Peak_type_count.pdf", width = 5.5, height = 5.5)
  ggplot() +
    geom_bar(data = temp, aes(x = Peak_Type,
                              y = log2(Count),
                              # y = Count,
                              fill = Peak_Type),
             stat = "identity", position = "stack") +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    scale_fill_manual(values = peak_type_color) +
    # scale_fill_manual(values = c("red", "gray")) +
    # geom_text(data = temp, aes(x = Peak_Type, y = log2(Count), label = Count),
    #           hjust = 1.3) +
    labs(x = "", y = "Log2(Count)", fill = "PeakType")
    # coord_flip()
  dev.off()
  
}

length(unique(peaksets_2[which(peaksets_2$DEGs == "YES"), "nearestGene"]))
temp <- peaksets_2[which(peaksets_2$pearson > 0.6),]
length(unique(temp[which(temp$DEGs == "YES"), "nearestGene"]))

# 绘图
pdf("Figure6_5_Peak_type_pie.pdf")
pie(c(length(Gene_Linked_Peaks),
      length(Lineage_Peaks),
      length(Global_Peaks),
      length(Other_Peaks)),
      labels = c("Gene-linked peaks",
                 "Lineage constitutive peaks",
                 "Global constitutive peaks",
                 "Other peaks"),
      col = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD"))
dev.off()

saveRDS(GCPs, "Figure6_GCPs.rds")
saveRDS(LCPs, "Figure6_LCPs.rds")

###### expression level
{
  peaksets_2$peakGroup <- "Other peaks"
  peaksets_2$peakGroup[which(peaksets_2$ID %in% Gene_Linked_Peaks)] <- "Gene-linked peaks"
  peaksets_2$peakGroup[which(peaksets_2$ID %in% Lineage_Peaks)] <- "Lineage constitutive peaks"
  peaksets_2$peakGroup[which(peaksets_2$ID %in% Global_Peaks)] <- "Global constitutive peaks"
  
  gene_peak <- as.data.frame.array(table(peaksets_2$nearestGene, peaksets_2$peakGroup))
  gene_peak <- as.matrix(gene_peak)
  gene_peak[which(gene_peak > 0)] <- 1
  
  # No Gene_linked peaks 不考虑
  type1_genes <- rownames(gene_peak)[which(gene_peak[,1] == 0)]
  with_gene_linked_peaks <- rownames(gene_peak)[which(gene_peak[,1] == 1)]
  
  # Gene_linked peaks & Lineage constitutive peaks & Global constitutive peaks
  type2_genes <- rownames(gene_peak)[which(rowSums(gene_peak[,1:3]) == 3)]
  sum(type2_genes %in% RNA_Celltye_DEG$gene)
  gene_peak2 <- gene_peak[-which(rownames(gene_peak) %in% c(type1_genes, type2_genes)),]
  
  # Gene_linked peaks & Lineage constitutive peaks
  type3_genes <- intersect(rownames(gene_peak2)[which(rowSums(gene_peak2[,c(1,3)]) == 2)],
                           rownames(gene_peak2)[which(gene_peak2[,c(2)] == 0)])
  intersect(type2_genes, type3_genes)
  gene_peak2 <- gene_peak[-which(rownames(gene_peak) %in% c(type1_genes, type2_genes, type3_genes)),]
  
  # Gene_linked peaks & Global constitutive peaks
  type4_genes <- intersect(rownames(gene_peak2)[which(rowSums(gene_peak2[,c(1,2)]) == 2)],
                           rownames(gene_peak2)[which(gene_peak2[,c(3)] == 0)])
  gene_peak2 <- gene_peak[-which(rownames(gene_peak) %in% c(type1_genes, type2_genes, type3_genes, type4_genes)),]
  
  # Gene_linked peaks
  type5_genes <- rownames(gene_peak2)
  
  length(type1_genes) + length(type2_genes) + length(type3_genes) + length(type4_genes) + length(type5_genes)
  nrow(gene_peak)
  
  
  gplot2 <- function(temp, temp2, ncol){
    ggplot(data = temp, aes(x = type, y =  temp[,ncol], fill = type)) +
      geom_boxplot(outlier.shape = NA) +
      # geom_point(data = temp2, aes(x = Group.1, y = x), color = "red", alpha = 0.6) +
      # geom_line(data = temp2, aes(x = Group.1, y = x, group = 1), color = "red", alpha = 0.6) +
      theme_bw() +
      ggtitle(colnames(temp)[ncol]) +
      coord_cartesian(ylim = c(0,0.55)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
            legend.position = "none",
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black", size = 10),
            plot.title = element_text(hjust = 0.5, size = 10, color = "black", face= "bold")) +
      labs(y = "Average Log-normalized RNA expression", x = "") +
      stat_compare_means(bracket.size = 0.3, label = "p.signif",
                         tip.length = 0.002,
                         comparisons = list(c("Only GLPs", "No GLPs"),
                                            # c("GLPs & LCPs", "No GLPs"),
                                            # c("GLPs & GCPs", "No GLPs"),
                                            # c("GLPs & LCPs & GCPs", "No GLPs"),
                                            c("GLPs & LCPs", "Only GLPs"),
                                            c("GLPs & GCPs", "Only GLPs"),
                                            c("GLPs & LCPs & GCPs", "Only GLPs")
                                            # c("GLPs & GCPs", "GLPs & LCPs"),
                                            # c("GLPs & LCPs & GCPs", "GLPs & LCPs"),
                                            # c("GLPs & LCPs & GCPs", "GLPs & GCPs")
                         ),
                         label.y = seq(from = 0.25, to = 0.45, length.out = 9),
                         vjust = 0.6)
  }
  # RNA_mean
  temp <- RNA_mean[,c(type1_genes, type5_genes, type3_genes, type4_genes, type2_genes)]
  temp <- as.data.frame(na.omit(temp))
  temp <- data.frame(t(temp), check.names = F)
  temp <- data.frame(type = c(rep("No GLPs", length(type1_genes)),
                              rep("Only GLPs", length(type5_genes)),
                              rep("GLPs & LCPs", length(type3_genes)),
                              rep("GLPs & GCPs", length(type4_genes)),
                              rep("GLPs & LCPs & GCPs", length(type2_genes))
                              ),
                     temp, check.names = F)
  temp$type <- factor(temp$type, levels = c("No GLPs",
                                            "Only GLPs",
                                            "GLPs & LCPs",
                                            "GLPs & GCPs",
                                            "GLPs & LCPs & GCPs"
                                            ))
  temp2 <- aggregate.data.frame(temp[,-1], by = list(temp$type), FUN = mean)
  rownames(temp2) <- temp2$Group.1
  temp2 <- temp2[,-1]
  temp2 <- t(temp2)
  # temp2 <- t(scale(t(temp2)))
  # temp2[which(temp2 > 1.5)] <- 1.5
  # temp2[which(temp2 < -0.6)] <- -0.6
  pdf("Figure6_peak_combination_expression_level_heatmap_2.pdf", height = 7, width = 4.5)
  pheatmap(temp2, cluster_cols = F, scale = "none",
           cluster_rows = F, border_color = NA, angle_col = "90")
  dev.off()
  
  
  plot_list <- lapply(2:ncol(temp), function(x){
    temp2 <- aggregate.data.frame(temp[,x], by = list(temp$type), FUN = median)
    gplot2(temp, temp2, x)
  })
  pdf("Figure6_5_peak_combination_expression_level.pdf", width = 12, height = 12)
  cowplot::plot_grid(plotlist = plot_list, nrow = 3)
  dev.off()
  
  
  RNA_mean_all <- rowSums(RNA_common_tissues@assays[["RNA"]]@data) / ncol(RNA_common_tissues@assays[["RNA"]]@data)
  temp <- RNA_mean_all[c(type1_genes, type5_genes, type3_genes, type4_genes, type2_genes)]
  temp <- data.frame(type = c(rep("No GLPs", length(type1_genes)),
                              rep("Only GLPs", length(type5_genes)),
                              rep("GLPs & LCPs", length(type3_genes)),
                              rep("GLPs & GCPs", length(type4_genes)),
                              rep("GLPs & LCPs & GCPs", length(type2_genes))
                              ),
                     value = temp, check.names = F)
  temp$type <- factor(temp$type, levels = c("No GLPs",
                                            "Only GLPs",
                                            "GLPs & LCPs",
                                            "GLPs & GCPs",
                                            "GLPs & LCPs & GCPs"))
  pdf("Figure6_5_peak_combination_expression_level_all.pdf", width = 3, height = 5)
  ggplot(data = temp, aes(x = type, y = value, fill = type)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    coord_cartesian(ylim = c(0,0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", size = 10),
          plot.title = element_text(hjust = 0.5, size = 10, color = "black", face= "bold")) +
    labs(y = "Average Log-normalized RNA expression", x = "") +
    stat_compare_means(bracket.size = 0.3, label = "p.signif",
                       tip.length = 0.002,
                       comparisons = list(c("Only GLPs", "No GLPs"),
                                          # c("GLPs & LCPs", "No GLPs"),
                                          # c("GLPs & GCPs", "No GLPs"),
                                          # c("GLPs & LCPs & GCPs", "No GLPs"),
                                          c("GLPs & LCPs", "Only GLPs"),
                                          c("GLPs & GCPs", "Only GLPs"),
                                          c("GLPs & LCPs & GCPs", "Only GLPs")
                                          ),
                       label.y = seq(from = 0.25, to = 0.4, length.out = 9),
                       vjust = 0.6)
  dev.off()
}

{
  # TF Motif
  MSU_RAP <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RAP-MSU_2023-03-15.txt",
                      sep = "\t", header = F)
  MSU_RAP <- MSU_RAP[-which(MSU_RAP$V2 == "None"),]
  MSU_RAP <- MSU_RAP[-which(MSU_RAP$V1 == "None"),]
  MSU <- c()
  row <- c()
  for (i in 1:nrow(MSU_RAP)) {
    temp <- unlist(strsplit(MSU_RAP[i,2], ",", fixed = T))
    row <- c(row, rep(i, length(temp)))
    MSU <- c(MSU, temp)
  }
  MSU_RAP <- data.frame(RAP = MSU_RAP[row,1],
                        MSU = MSU)
  MSU_RAP$MSU <- unlist(lapply(MSU_RAP$MSU, function(x){
    unlist(strsplit(x, ".", fixed = T))[1]
  }))
  MSU_RAP <- MSU_RAP[!duplicated(MSU_RAP),]
  
  motif <- read_meme("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs.meme")
  motif_all <- convert_motifs(motif, class = "TFBSTools-PWMatrix")
  motif_all <- do.call(PWMatrixList, motif_all)
  MSU_ID <- unlist(lapply(motif_all@listData, function(x){
    x@name
  }))
  names(motif_all) <- MSU_ID
  RAP_ID <- MSU_RAP[match(MSU_ID, MSU_RAP$MSU),]
  RAP_ID <- na.omit(RAP_ID)
  motif_all <- motif_all[RAP_ID$MSU]
  names(motif_all) <- RAP_ID$RAP
  motif_all <- motif_all[!duplicated(names(motif_all))] # ID duplicated !!!
  rownames(gene_id_map) <- gene_id_map$gene_id
  motif_all <- motif_all[names(motif_all) %in% gene_id_map$gene_id]
  Symbol <- match(names(motif_all), gene_id_map$gene_id)
  Symbol <- gene_id_map[Symbol,]
  Annotation_genes <- read.csv("./Ref/IRGSP-1.0_representative_annotation_2023-03-15.tsv",
                               sep = "\t", header = T, check.names = F)
  Annotation_genes <- Annotation_genes[,c("Locus_ID", "CGSNL Gene Symbol", "CGSNL Gene Name")]
  Annotation_genes <- Annotation_genes[Annotation_genes$Locus_ID %in% names(motif_all),]
  Annotation_genes <- Annotation_genes[!duplicated(Annotation_genes$Locus_ID),]
  Annotation_genes$`CGSNL Gene Symbol`[168] <- "EREBP96"
  Annotation_genes$`CGSNL Gene Symbol`[23] <- "EREBP127"
  Annotation_genes$`CGSNL Gene Symbol`[244] <- "RID1"
  Annotation_genes$`CGSNL Gene Symbol` <- paste0("Os",Annotation_genes$`CGSNL Gene Symbol`)
  empty_row <- which(Annotation_genes$`CGSNL Gene Symbol` == "Os")
  Annotation_genes$`CGSNL Gene Symbol`[empty_row] <- Annotation_genes$Locus_ID[empty_row]
  Annotation_genes$DataSets_Symbol <- gene_id_map[Annotation_genes$Locus_ID, "symbol"]
  rownames(Annotation_genes) <- Annotation_genes$Locus_ID
  names(motif_all) <- Annotation_genes[names(motif_all), "CGSNL Gene Symbol"]
  
  motifPositions <- motifmatchr::matchMotifs(
    pwms = motif_all,
    subject = peakSet,
    genome = BSgenome.OSativa.NCBI.IRGSPv1.0, 
    out = "positions", 
    p.cutoff = 5e-05, 
    w = 7
  )
  nO <- lapply(motifPositions, function(x)length(x)) %>% unlist
  mF <- names(which(nO == 0))
  if(length(mF) > 0){
    #Filter
    motifPositions <- motifPositions[nO > 0]
    motifSummary <- motifSummary[names(motifPositions),,drop=FALSE]
    motifs <- motifs[names(motifPositions)]
  } 
  allPositions <- unlist(motifPositions)
  overlapMotifs <- findOverlaps(peakSet, allPositions, ignore.strand = TRUE)
  overlapMotifs <- as.data.frame(overlapMotifs)
  peak_id <- peakSet$Peak_id
  allPositions$TF <- allPositions@ranges@NAMES
  allPositions@ranges@NAMES <- NULL
  allPositions_df <- as.data.frame(allPositions)
  peak_TF <- data.frame(Peak_nearestGene = peakSet$nearestGene[overlapMotifs$queryHits],
                            Peak_id = peakSet$Peak_id[overlapMotifs$queryHits],
                            Motif_chr = allPositions_df$seqnames[overlapMotifs$subjectHits],
                            Motif_start = allPositions_df$start[overlapMotifs$subjectHits],
                            Motif_end = allPositions_df$end[overlapMotifs$subjectHits],
                            Motif_TF = allPositions_df$TF[overlapMotifs$subjectHits],
                            Motif_strand = allPositions_df$strand[overlapMotifs$subjectHits],
                            Motif_score = allPositions_df$score[overlapMotifs$subjectHits],
                            Motif_width = allPositions_df$width[overlapMotifs$subjectHits]
  )
  gene_df <- as.data.frame(gene_annotation@listData[["genes"]])
  rownames(gene_df) <- gene_df$symbol
  peak_TF$Peak_nearestGene_strand <- gene_df[peak_TF$Peak_nearestGene, "strand"]
  peak_TF <- peak_TF[which(peak_TF$Motif_strand == peak_TF$Peak_nearestGene_strand),]
  peak_TF$Peak_TF <- paste0(peak_TF$Peak_id, "_", peak_TF$Motif_TF)
  sort(table(peak_TF$Peak_TF), decreasing = T)
  
  peak_tf_matches <- Matrix::sparseMatrix(
    i = match(peak_TF$Peak_id, peakSet$Peak_id), # peaks
    j = match(peak_TF$Motif_TF, names(motif_all)), # motifs
    x = rep(1, nrow(peak_TF)),
    dims = c(length(peakSet), length(motif_all)),
    dimnames = list(peakSet$Peak_id, names(motif_all))
  )
  peak_tf_matches["7_19533282_19533782", "OsERF19"]
  peak_tf_matches["3_11031118_11031618", "OsERF19"]
  matches <- t(peak_tf_matches)
  matches@x <- rep(1, length(matches@x))
  matches <- as.matrix(matches)
  matches <- t(matches)
  matches <- as.data.frame(matches)
  
  # GLPs & LCPs
  GLPs_LCPs <- c()
  for (i in type3_genes) {
    print(i)
    temp_gene_peaks <- peaksets_2[which(peaksets_2$nearestGene == i),]
    temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakGroup %in% c("Gene-linked peaks",
                                                                              "Lineage constitutive peaks",
                                                                              "Global constitutive peaks")),]
    temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakType %in% c("Promoter", "Distal")),]
    if (length(table(temp_gene_peaks$peakGroup)) < 2) {
      GLPs_LCPs <- c(GLPs_LCPs, -1)
    } else {
      if (nrow(temp_gene_peaks) > 1) {
        temp_gene_peaks_match <- matches[rownames(temp_gene_peaks),]
        temp_gene_peaks_match <- aggregate.data.frame(temp_gene_peaks_match,
                                                      by = list(temp_gene_peaks$peakGroup),
                                                      FUN = max)
        temp_gene_peaks_match <- temp_gene_peaks_match[,-1]
        if (length(which(colSums(temp_gene_peaks_match) == nrow(temp_gene_peaks_match))) > 0) {
          GLPs_LCPs <- c(GLPs_LCPs, 1)
        } else {
          GLPs_LCPs <- c(GLPs_LCPs, 0)
        }
      } else {
        GLPs_LCPs <- c(GLPs_LCPs, -1)
      }
    }
    
  }
  table(GLPs_LCPs)
  
  # GLPs & GCPs
  GLPs_GCPs <- c()
  for (i in type4_genes) {
    print(i)
    temp_gene_peaks <- peaksets_2[which(peaksets_2$nearestGene == i),]
    temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakGroup %in% c("Gene-linked peaks",
                                                                              "Lineage constitutive peaks",
                                                                              "Global constitutive peaks")),]
    temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakType %in% c("Promoter", "Distal")),]
    if (length(table(temp_gene_peaks$peakGroup)) < 2) {
      GLPs_GCPs <- c(GLPs_GCPs, -1)
    } else {
      if (nrow(temp_gene_peaks) > 1) {
        temp_gene_peaks_match <- matches[rownames(temp_gene_peaks),]
        temp_gene_peaks_match <- aggregate.data.frame(temp_gene_peaks_match,
                                                      by = list(temp_gene_peaks$peakGroup),
                                                      FUN = max)
        temp_gene_peaks_match <- temp_gene_peaks_match[,-1]
        if (length(which(colSums(temp_gene_peaks_match) == nrow(temp_gene_peaks_match))) > 0) {
          GLPs_GCPs <- c(GLPs_GCPs, 1)
        } else {
          GLPs_GCPs <- c(GLPs_GCPs, 0)
        }
      } else {
        GLPs_GCPs <- c(GLPs_GCPs, -1)
      }
    }
  }
  table(GLPs_GCPs)
  
  # GLPs & LCPs & GCPs
  GLPs_LCPs_GCPs <- c()
  for (i in type2_genes) {
    print(i)
    temp_gene_peaks <- peaksets_2[which(peaksets_2$nearestGene == i),]
    temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakGroup %in% c("Gene-linked peaks",
                                                                              "Lineage constitutive peaks",
                                                                              "Global constitutive peaks")),]
    temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakType %in% c("Promoter", "Distal")),]
    if (length(table(temp_gene_peaks$peakGroup)) < 3) {
      GLPs_LCPs_GCPs <- c(GLPs_LCPs_GCPs, -1)
    } else {
      if (nrow(temp_gene_peaks) > 2) {
        temp_gene_peaks_match <- matches[rownames(temp_gene_peaks),]
        temp_gene_peaks_match <- aggregate.data.frame(temp_gene_peaks_match,
                                                      by = list(temp_gene_peaks$peakGroup),
                                                      FUN = max)
        temp_gene_peaks_match <- temp_gene_peaks_match[,-1]
        if (length(which(colSums(temp_gene_peaks_match) == nrow(temp_gene_peaks_match))) > 0) {
          GLPs_LCPs_GCPs <- c(GLPs_LCPs_GCPs, 1)
        } else {
          GLPs_LCPs_GCPs <- c(GLPs_LCPs_GCPs, 0)
        }
      } else {
        GLPs_LCPs_GCPs <- c(GLPs_LCPs_GCPs, -1)
      }
    }
  }
  table(GLPs_LCPs_GCPs)
  
  #### random sequence, 保持GLP不变
  {
    # 产生BSgenome
    {
      set.seed(123)
      random_sequence <- lapply(1:10, function(x){
        paste0(sample(c("A","T","G","C"), size = 501*nrow(peaksets_2), replace = TRUE, prob = c(0.25,0.25,0.25,0.25)),
               collapse = "")
      })
      random_sequence <- unlist(random_sequence)
      names(random_sequence) <- paste0("random_", 1:10)
      random_sequence <- DNAStringSet(random_sequence)
      writeXStringSet(random_sequence,
                      filepath = paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Figure6_Permutation_random_genome/fasta_sequences_random.fasta"),
                      format = "fasta")
      library(BSgenome)
      setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Figure6_Permutation_random_genome")
      forgeBSgenomeDataPkg("./BSgenome.Oryza.sativa.IRGSP.1.0.random-seed",
                           seqs_srcdir = getwd(), destdir = getwd(),
                           verbose = T)
      system(command = "R CMD build BSgenome.Oryza.sativa.IRGSP.1.0.random")
      system(command = "R CMD check --no-manual BSgenome.Oryza.sativa.IRGSP.1.0.random_1.0.0.tar.gz")
      system(command = "R CMD INSTALL BSgenome.Oryza.sativa.IRGSP.1.0.random_1.0.0.tar.gz")
    }
    library(BSgenome.Oryza.sativa.IRGSP.1.0.random)
    random_GLPs_LCPs <- c()
    random_GLPs_GCPs <- c()
    random_GLPs_LCPs_GCPs <- c()
    
    for (r in 1:10) {
      print(r)
      start <- unlist(lapply(1:nrow(peaksets_2), function(x){
        (x-1)*501 + 1
      }))
      end <- start + 500
      random_seq_region <- data.frame(chr = paste0("random_", r),
                                      start = start,
                                      end = end,
                                      strand = "*")
      random_seq_region <- makeGRangesFromDataFrame(random_seq_region,
                                                    seqinfo = SeqinfoForBSGenome(BSgenome.Oryza.sativa.IRGSP.1.0.random))
      random_seq_region$Peak_id <- peaksets_2$ID
      temp_match <- motifmatchr::matchMotifs(
        pwms = motif_all,
        subject = random_seq_region,
        genome = BSgenome.Oryza.sativa.IRGSP.1.0.random,
        out = "positions", 
        p.cutoff = 5e-5, 
        w = 7
      )
      nO <- lapply(temp_match, function(x)length(x)) %>% unlist
      mF <- names(which(nO == 0))
      motifs <- motif_all
      if(length(mF) > 0){
        #Filter
        temp_match <- temp_match[nO > 0]
        motifSummary <- motifSummary[names(temp_match),,drop=FALSE]
        motifs <- motif_all[names(temp_match)]
      }
      
      temp_allPositions <- unlist(temp_match)
      temp_overlapMotifs <- findOverlaps(random_seq_region, temp_allPositions, ignore.strand = TRUE)
      temp_overlapMotifs <- as.data.frame(temp_overlapMotifs)
      
      temp_allPositions$TF <- temp_allPositions@ranges@NAMES
      temp_allPositions@ranges@NAMES <- NULL
      temp_allPositions_df <- as.data.frame(temp_allPositions)
      temp_random_peak_TF <- data.frame(Peak_id = peaksets_2$ID[temp_overlapMotifs$queryHits],
                                        Motif_chr = temp_allPositions_df$seqnames[temp_overlapMotifs$subjectHits],
                                        Motif_TF = temp_allPositions_df$TF[temp_overlapMotifs$subjectHits],
                                        Motif_strand = temp_allPositions_df$strand[temp_overlapMotifs$subjectHits]
      )
      temp_random_peak_TF$Peak_nearestGene <- peaksets_2[temp_random_peak_TF$Peak_id, "nearestGene"]
      gene_df <- as.data.frame(gene_annotation@listData[["genes"]])
      rownames(gene_df) <- gene_df$symbol
      temp_random_peak_TF$Peak_nearestGene_strand <- gene_df[temp_random_peak_TF$Peak_nearestGene, "strand"]
      temp_random_peak_TF <- temp_random_peak_TF[which(temp_random_peak_TF$Motif_strand == temp_random_peak_TF$Peak_nearestGene_strand),]
      
      temp_peak_tf_matches <- Matrix::sparseMatrix(
        i = match(temp_random_peak_TF$Peak_id, random_seq_region$Peak_id), # peaks
        j = match(temp_random_peak_TF$Motif_TF, names(motifs)), # motifs
        x = rep(1, nrow(temp_random_peak_TF)),
        dims = c(length(random_seq_region), length(motifs)),
        dimnames = list(random_seq_region$Peak_id, names(motifs))
      )
      temp_peak_tf_matches@x <- rep(1, length(temp_peak_tf_matches@x))
      temp_peak_tf_matches <- as.data.frame(temp_peak_tf_matches)
      
      # GLPs & LCPs
      temp_GLPs_LCPs <- c()
      for (i in type3_genes) {
        print(i)
        temp_gene_peaks <- peaksets_2[which(peaksets_2$nearestGene == i),]
        temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakGroup %in% c("Gene-linked peaks",
                                                                                  "Lineage constitutive peaks",
                                                                                  "Global constitutive peaks")),]
        temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakType %in% c("Promoter", "Distal")),]
        GLP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Gene-linked peaks"), "ID"]
        LCP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Lineage constitutive peaks"), "ID"]
        if (length(table(temp_gene_peaks$peakGroup)) < 2) {
          temp_GLPs_LCPs <- c(temp_GLPs_LCPs, -1)
        } else {
          if (nrow(temp_gene_peaks) > 1) {
            common_tf <- intersect(colnames(temp_peak_tf_matches), colnames(matches))
            temp_gene_peaks_match <- as.data.frame(rbind(matches[GLP,common_tf],
                                                         temp_peak_tf_matches[LCP,common_tf]))
            temp_gene_peaks_match <- aggregate.data.frame(temp_gene_peaks_match,
                                                          by = list(c(rep("Gene-linked peaks", length(GLP)),
                                                                      rep("Lineage constitutive peaks", length(LCP)))),
                                                          FUN = max)
            temp_gene_peaks_match <- temp_gene_peaks_match[,-1]
            if (length(which(colSums(temp_gene_peaks_match) == nrow(temp_gene_peaks_match))) > 0) {
              temp_GLPs_LCPs <- c(temp_GLPs_LCPs, 1)
            } else {
              temp_GLPs_LCPs <- c(temp_GLPs_LCPs, 0)
            }
          } else {
            temp_GLPs_LCPs <- c(temp_GLPs_LCPs, -1)
          }
        }
        
      }
      table(temp_GLPs_LCPs)
      
      # GLPs & GCPs
      temp_GLPs_GCPs <- c()
      for (i in type4_genes) {
        print(i)
        temp_gene_peaks <- peaksets_2[which(peaksets_2$nearestGene == i),]
        temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakGroup %in% c("Gene-linked peaks",
                                                                                  "Lineage constitutive peaks",
                                                                                  "Global constitutive peaks")),]
        temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakType %in% c("Promoter", "Distal")),]
        GLP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Gene-linked peaks"), "ID"]
        GCP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Global constitutive peaks"), "ID"]
        
        if (length(table(temp_gene_peaks$peakGroup)) < 2) {
          temp_GLPs_GCPs <- c(temp_GLPs_GCPs, -1)
        } else {
          if (nrow(temp_gene_peaks) > 1) {
            common_tf <- intersect(colnames(temp_peak_tf_matches), colnames(matches))
            temp_gene_peaks_match <- as.data.frame(rbind(matches[GLP,common_tf],
                                                         temp_peak_tf_matches[GCP,common_tf]))
            temp_gene_peaks_match <- aggregate.data.frame(temp_gene_peaks_match,
                                                          by = list(c(rep("Gene-linked peaks", length(GLP)),
                                                                      rep("Global constitutive peaks", length(GCP)))),
                                                          FUN = max)
            temp_gene_peaks_match <- temp_gene_peaks_match[,-1]
            if (length(which(colSums(temp_gene_peaks_match) == nrow(temp_gene_peaks_match))) > 0) {
              temp_GLPs_GCPs <- c(temp_GLPs_GCPs, 1)
            } else {
              temp_GLPs_GCPs <- c(temp_GLPs_GCPs, 0)
            }
          } else {
            temp_GLPs_GCPs <- c(temp_GLPs_GCPs, -1)
          }
        }
      }
      table(temp_GLPs_GCPs)
      
      # GLPs & LCPs & GCPs
      temp_GLPs_LCPs_GCPs <- c()
      for (i in type2_genes) {
        print(i)
        temp_gene_peaks <- peaksets_2[which(peaksets_2$nearestGene == i),]
        temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakGroup %in% c("Gene-linked peaks",
                                                                                  "Lineage constitutive peaks",
                                                                                  "Global constitutive peaks")),]
        GLP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Gene-linked peaks"), "ID"]
        GCP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Global constitutive peaks"), "ID"]
        LCP <- temp_gene_peaks[which(temp_gene_peaks$peakGroup == "Lineage constitutive peaks"), "ID"]
        
        temp_gene_peaks <- temp_gene_peaks[which(temp_gene_peaks$peakType %in% c("Promoter", "Distal")),]
        if (length(table(temp_gene_peaks$peakGroup)) < 3) {
          temp_GLPs_LCPs_GCPs <- c(temp_GLPs_LCPs_GCPs, -1)
        } else {
          if (nrow(temp_gene_peaks) > 2) {
            common_tf <- intersect(colnames(temp_peak_tf_matches), colnames(matches))
            temp_gene_peaks_match <- as.data.frame(rbind(matches[GLP,common_tf],
                                                         temp_peak_tf_matches[GCP,common_tf]))
            temp_gene_peaks_match <- as.data.frame(rbind(temp_gene_peaks_match,
                                                         temp_peak_tf_matches[LCP,common_tf]))
            
            temp_gene_peaks_match <- aggregate.data.frame(temp_gene_peaks_match,
                                                          by = list(c(rep("Gene-linked peaks", length(GLP)),
                                                                      rep("Global constitutive peaks", length(GCP)),
                                                                      rep("Lineage constitutive peaks", length(LCP)))),
                                                          FUN = max)
            temp_gene_peaks_match <- temp_gene_peaks_match[,-1]
            if (length(which(colSums(temp_gene_peaks_match) == nrow(temp_gene_peaks_match))) > 0) {
              temp_GLPs_LCPs_GCPs <- c(temp_GLPs_LCPs_GCPs, 1)
            } else {
              temp_GLPs_LCPs_GCPs <- c(temp_GLPs_LCPs_GCPs, 0)
            }
          } else {
            temp_GLPs_LCPs_GCPs <- c(temp_GLPs_LCPs_GCPs, -1)
          }
        }
      }
      
      temp_GLPs_LCPs <- as.data.frame(table(temp_GLPs_LCPs))
      temp_GLPs_LCPs <- temp_GLPs_LCPs[3,2] / sum(c(temp_GLPs_LCPs[2,2], temp_GLPs_LCPs[3,2]))
      temp_GLPs_GCPs <- as.data.frame(table(temp_GLPs_GCPs))
      temp_GLPs_GCPs <- temp_GLPs_GCPs[3,2] / sum(c(temp_GLPs_GCPs[2,2], temp_GLPs_GCPs[3,2]))
      temp_GLPs_LCPs_GCPs <- as.data.frame(table(temp_GLPs_LCPs_GCPs))
      temp_GLPs_LCPs_GCPs <- temp_GLPs_LCPs_GCPs[3,2] / sum(c(temp_GLPs_LCPs_GCPs[2,2], temp_GLPs_LCPs_GCPs[3,2]))
      random_GLPs_LCPs <- c(random_GLPs_LCPs,
                            temp_GLPs_LCPs)
      random_GLPs_GCPs <- c(random_GLPs_GCPs,
                            temp_GLPs_GCPs)
      random_GLPs_LCPs_GCPs <- c(random_GLPs_LCPs_GCPs,
                                 temp_GLPs_LCPs_GCPs)
    }
    
    
    
  }
  
  
  
  #### random region
  {
    peaksets_2_temp <- peaksets_2
    peaksets_2_temp <- peaksets_2_temp[,c("nearestGene", "peakType", "peakGroup", "ID")]
    GLPs <- peaksets_2_temp[which(peaksets_2_temp$peakGroup == "Gene-linked peaks"),]
    GCPs <- peaksets_2_temp[which(peaksets_2_temp$peakGroup == "Global constitutive peaks"),]
    LCPs <- peaksets_2_temp[which(peaksets_2_temp$peakGroup == "Lineage constitutive peaks"),]
    Other <- peaksets_2_temp[which(peaksets_2_temp$peakGroup == "Other peaks"),]
    
    set.seed(123)
    library(GenomicRanges)
    chr_size <- genome_annotation@listData[["chromSizes"]]@seqinfo@seqlengths[1:12]
    names(chr_size) <- genome_annotation@listData[["chromSizes"]]@seqinfo@seqnames[1:12]
    
    chr_random <- sample(c(1:12), (nrow(GCPs)+nrow(LCPs)) * 1000, replace = T,
                         prob = chr_size / sum(chr_size))
    chr_random <- sort(chr_random)
    chr_num <- table(chr_random)
    
    start_random <- c()
    for (i in names(chr_num)) {
      start_random <- c(start_random, sample(1:(chr_size[i]-500), chr_num[i]))
    }
    
    length(start_random) == length(chr_random)
    end_random <- start_random + 500
    
    
    random_range <- data.frame(chr = chr_random, start = start_random, end = end_random,
                               strand = "*")
    for (i in 1:20) {
      print(i)
      random_range <- random_range[sample(1:nrow(random_range),
                                          nrow(random_range), replace = F), ]
    }
    
    random_range <- makeGRangesFromDataFrame(random_range)
    
    random_range_temp <- as.data.frame(random_range)
    random_peaks <- paste0(random_range_temp$seqnames, "_", random_range_temp$start, "_", random_range_temp$end)
    rm(random_range_temp)
    random_range@elementMetadata@listData$PeakID <- random_peaks
    saveRDS(random_range, "Figure6_5_random_range.rds")
    
    random_range <- readRDS("Figure6_5_random_range.rds")
    
    temp_list <- list(GLPs = GLPs, GCPs = GCPs, LCPs = LCPs, Other = Other)
    temp_df <- temp_list$GLPs
    temp_df <- temp_df[,rep("ID", 1000)]
    colnames(temp_df) <- paste0("random_", 1:1000)
    temp_list$GLPs <- as.data.frame(cbind(temp_list$GLPs,
                                          temp_df))
    temp_df <- temp_list$GCPs
    temp_df <- temp_df[,rep("ID", 1000)]
    colnames(temp_df) <- paste0("random_", 1:1000)
    temp_list$GCPs <- as.data.frame(cbind(temp_list$GCPs,
                                          temp_df))
    temp_df <- temp_list$LCPs
    temp_df <- temp_df[,rep("ID", 1000)]
    colnames(temp_df) <- paste0("random_", 1:1000)
    temp_list$LCPs <- as.data.frame(cbind(temp_list$LCPs,
                                          temp_df))
    
    
    
    
    # TRUE
    (nrow(temp_list$GCPs) + nrow(temp_list$LCPs)) * 1000 == length(random_peaks)
    
    temp_df <- random_peaks[1:(1000*nrow(temp_list$GCPs))]
    sum(duplicated(temp_df))
    temp_df <- matrix(temp_df, nrow = nrow(temp_list$GCPs), ncol = 1000)
    colnames(temp_df) <- paste0("random_", 1:1000)
    temp_list$GCPs <- as.data.frame(cbind(temp_list$GCPs,
                                          temp_df))
    
    temp_df <- random_peaks[((1000*nrow(temp_list$GCPs))+1):length(random_peaks)]
    sum(duplicated(temp_df))
    temp_df <- matrix(temp_df, nrow = nrow(temp_list$LCPs), ncol = 1000)
    colnames(temp_df) <- paste0("random_", 1:1000)
    temp_list$LCPs <- as.data.frame(cbind(temp_list$LCPs,
                                          temp_df))
    
    # random match
    motif_match_list <- c()
    for (i in paste0("random_", 1:1000)) {
      print(i)
      temp_peak <- c(temp_list$GCPs[,i], temp_list$LCPs[,i])
      temp_random_range <- random_range[match(temp_peak,
                                              random_range@elementMetadata@listData[["PeakID"]]),]
      motifPositions <- motifmatchr::matchMotifs(
        pwms = motif_all,
        subject = temp_random_range,
        # subject = peaksets,
        genome = genome_annotation$genome, 
        out = "positions", 
        p.cutoff = 5e-5, 
        w = 7
      )
      nO <- lapply(motifPositions, length) %>% unlist
      mF <- names(which(nO == 0))
      
      if(length(mF) > 0){
        #Filter
        motifPositions <- motifPositions[nO > 0]
        motifs_temp <- motif_all[names(motifPositions)]
      } else {
        motifs_temp <- motif_all
      }
      
      allPositions <- unlist(motifPositions)
      # overlapMotifs <- findOverlaps(temp_random_range, allPositions, ignore.strand = TRUE)
      overlapMotifs <- findOverlaps(peaksets, allPositions, ignore.strand = TRUE)
      
      motifMat <- Matrix::sparseMatrix(
        i = queryHits(overlapMotifs),
        j = match(names(allPositions),names(motifPositions))[subjectHits(overlapMotifs)],
        x = rep(TRUE, length(overlapMotifs)),
        # dims = c(length(temp_random_range), length(motifPositions))
        dims = c(length(peaksets), length(motifPositions))
      )
      colnames(motifMat) <- names(motifPositions)
      rownames(motifMat) <- temp_peak
      motif_match_list <- c(motif_match_list, list(motifMat))
    }
    
    quantile(rowSums(motifMat))
    
    # for (i in 1:1000) {
    #   rownames(motif_match_list[[i]]) <- c(temp_list$GCPs[,paste0("random_", i)], temp_list$LCPs[,paste0("random_", i)])
    # }
    
    sum(temp_list$GCPs[,paste0("random_", i)] %in% rownames(motif_match_list[[i]]))
    sum(temp_list$GCPs[,paste0("random_", 3)] %in% random_peaks)
    motif_match_list[[i]][temp_list$GCPs[,paste0("random_", i)],]
    
    saveRDS(motif_match_list, "Figure6_5_motif_match_list.rds")
    
    # merged <- rlist::list.rbind(motif_match_list)
    # saveRDS(merged, "merged_motif_match.rds")
    
    random_result <- data.frame(type = c("GLPs & LCPs", "GLPs & GCPs", "GLPs & LCPs & GCPs"))
    
    save(random_peaks,
         temp_list,
         motif_match_list,
         peaksets_2,
         random_result,
         type1_genes,
         type2_genes,
         type3_genes,
         type4_genes,
         type5_genes,
         file = "./Figure6_Imputation_P/GLP_LCP_GCP.RData")
    
    
    # 
    # matches_random <- getMatches(ATAC_common_tissues, name = "TF-random_peak")
    # 
    # r1 <- SummarizedExperiment::rowRanges(matches_random)
    # pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
    # 
    # matches <- matches@assays@data@listData[["matches"]]
    # rownames(matches) <- pr1
    # matches <- matches + 0
    # matches <- as.matrix(matches)
    # matches <- matches[peaksets_2$ID,]
    
    
    # temp_list <- list(GLPs = GLPs, GCPs = GCPs, LCPs = LCPs, Other = Other)
    
    # for (i in names(temp_list)) {
    #   print(i)
    #   if (i %in% c("GLPs","Other")) {
    #     temp_df <- temp_list[[i]]
    #     temp_df_2 <- temp_df[,rep("ID", 1000)]
    #     colnames(temp_df_2) <- paste0("random_", 1:1000)
    #     temp_list[[i]] <- as.data.frame(cbind(temp_df, temp_df_2))
    #   } else {
    #     temp_df <- temp_list[[i]]
    #     unique(temp_df$peakType)
    #     promoter <- temp_df[which(temp_df$peakType == "Promoter"),]
    #     exonic <- temp_df[which(temp_df$peakType == "Exonic"),]
    #     distal <- temp_df[which(temp_df$peakType == "Distal"),]
    #     intronic <- temp_df[which(temp_df$peakType == "Intronic"),]
    #     temp_df_2 <- list(promoter = promoter, exonic = exonic, distal = distal, intronic = intronic)
    #     for (j in names(temp_df_2)) {
    #       print(j)
    #       temp_df_3 <- temp_df_2[[j]]
    #       random_time <- 0
    #       random_peak_list <- c()
    #       while (random_time < 1000) {
    #         random_peak <- sample(temp_df_3$ID, nrow(temp_df_3), replace = F)
    #         random_gene <- peaksets_2[random_peak, "nearestGene"]
    #         if (sum(random_gene == temp_df_3$nearestGene) > 0) {
    #           random_time <- random_time
    #         } else {
    #           # print(random_time)
    #           random_time <- random_time + 1
    #           random_peak_list <- c(random_peak_list, list(as.data.frame(random_peak)))
    #         }
    #       }
    #       random_peak <- purrr::list_cbind(random_peak_list)
    #       colnames(random_peak) <- paste0("random_", 1:1000)
    #       temp_df_3 <- as.data.frame(cbind(temp_df_3,
    #                                        random_peak))
    #       temp_df_2[[j]] <- temp_df_3
    #     }
    #     temp_df <- purrr::list_rbind(temp_df_2)
    #     temp_list[[i]] <- temp_df
    #   }
      
      
      # random_peak <- lapply(1:1000, function(x){
      #   data.frame(sample(temp_df$ID, nrow(temp_df), replace = F))
      # })
      # random_peak <- list_cbind(random_peak)
      # colnames(random_peak) <- paste0("random_", 1:1000)
      # temp_df <- as.data.frame(cbind(temp_df,
      #                                random_peak))
      # temp_list[[i]] <- temp_df
    }
    
    peaksets_2_random <- purrr::list_rbind(temp_list)
    identical(peaksets_2[peaksets_2_random$random_1, "peakType"],
              peaksets_2_random$peakType
    )
    
    random_result <- data.frame(type = c("GLPs & LCPs", "GLPs & GCPs", "GLPs & LCPs & GCPs"))
    save(matches,
         peaksets_2,
         peaksets_2_random,
         random_result,
         type1_genes,
         type2_genes,
         type3_genes,
         type4_genes,
         type5_genes,
         file = "./Figure6_Imputation_P/GLP_LCP_GCP.RData")
    
  }
  
#### random result
{
  random_result <- c()
  for (i in 1:100) {
    random_result_temp <- readRDS(paste0("./Figure6_Imputation_P/random_result_", i, ".rds"))
    colnames(random_result_temp) <- random_result_temp[1,]
    random_result_temp <- random_result_temp[-1,]
    random_result <- as.data.frame(rbind(random_result,
                                         random_result_temp))
  }
  random_result$`GLPs & LCPs` <- as.numeric(random_result$`GLPs & LCPs`)
  random_result$`GLPs & GCPs` <- as.numeric(random_result$`GLPs & GCPs`)
  random_result$`GLPs & LCPs & GCPs` <- as.numeric(random_result$`GLPs & LCPs & GCPs`)
  sum(random_result$`GLPs & LCPs` > 0.9172) / nrow(random_result) # 0 / 0
  sum(random_result$`GLPs & GCPs` > 0.8995) / nrow(random_result) # 0.149 / 0.226 / 0
  sum(random_result$`GLPs & LCPs & GCPs` > 0.5879) / nrow(random_result) # 0.479 / 0.405 / 0
  colSums(random_result) / nrow(random_result)
  
  random_result_melt <- data.frame(GeneType = c(rep("GLPs & LCPs", 1000),
                                                rep("GLPs & GCPs", 1000),
                                                rep("GLPs & LCPs & GCPs", 1000)),
                                   Percent = c(random_result$`GLPs & LCPs`,
                                               random_result$`GLPs & GCPs`,
                                               random_result$`GLPs & LCPs & GCPs`))
  
  random_result_melt <- data.frame(GeneType = c(rep("GLPs & LCPs", 10),
                                                rep("GLPs & GCPs", 10),
                                                rep("GLPs & LCPs & GCPs", 10)),
                                   Percent = c(random_GLPs_LCPs,
                                               random_GLPs_GCPs,
                                               random_GLPs_LCPs_GCPs))
  table(GLPs_GCPs)
  pdf("Figure6_5_genes_peaktype_same_TF_percent_sig.pdf", width = 5, height = 5)
  ggplot() +
    geom_boxplot(data = random_result_melt[which(random_result_melt$GeneType == "GLPs & LCPs"),],
                 aes(x = 0.9, y = Percent),
                 width = 0.3) +
    geom_bar(aes(x = 0.5, y = 0.6986688), # GLPs & LCPs, 1732 / (1732 + 747)
             stat = "identity", width = 0.3, color = "#00BF7D", fill = "#00BF7D") +
    geom_boxplot(data = random_result_melt[which(random_result_melt$GeneType == "GLPs & GCPs"),],
                 aes(x = 1.9, y = Percent),
                 width = 0.3) +
    geom_bar(aes(x = 1.5, y = 0.6890459), # GLPs & GCPs, 195 / (195 + 88)
             stat = "identity", width = 0.3, color = "#00B0F6", fill = "#00B0F6") +
    geom_boxplot(data = random_result_melt[which(random_result_melt$GeneType == "GLPs & LCPs & GCPs"),],
                 aes(x = 2.9, y = Percent),
                 width = 0.3) +
    geom_bar(aes(x = 2.5, y = 0.2533333), # GLPs & LCPs & GCPs, 114 / (114 + 336)
             stat = "identity", width = 0.3, color = "#E76BF3", fill = "#E76BF3") +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black")) +
    labs(y = "Percentage", x = "")
  dev.off()
  hue_pal()(5)
}
  
#### gene browse
{
  peak_group <- data.frame(peakType = c(rep("Gene-linked peaks", length(Gene_Linked_Peaks_DEGs)),
                                        rep("Lineage constitutive peaks", length(Lineage_Peaks_DEGs)),
                                        rep("Global constitutive peaks", length(Global_Peaks_DEGs)),
                                        rep("Other peaks", length(Other_Peaks_DEGs))),
                           peakID = c(Gene_Linked_Peaks_DEGs,
                                      Lineage_Peaks_DEGs,
                                      Global_Peaks_DEGs,
                                      Other_Peaks_DEGs))
  
  peak_group$nearestGene <- peaksets[peak_group$peakID, "nearestGene"]
  peak_group$peakRegion <- peaksets[peak_group$peakID, "peakType"]
  
  Other_Peaks_PD <- peaksets[Other_Peaks_DEGs,]
  Other_Peaks_PD <- Other_Peaks_PD[which(Other_Peaks_PD$peakType %in% c("Promoter",
                                                                        "Distal")),]
  Other_Peaks_PD <- Other_Peaks_PD$ID
  
  Lineage_Peaks_PD <- peaksets[Lineage_Peaks_DEGs,]
  Lineage_Peaks_PD <- Lineage_Peaks_PD[which(Lineage_Peaks_PD$peakType %in% c("Promoter",
                                                                              "Distal")),]
  Lineage_Peaks_PD <- Lineage_Peaks_PD$ID
  
  Global_Peaks_PD <- peaksets[Global_Peaks_DEGs,]
  Global_Peaks_PD <- Global_Peaks_PD[which(Global_Peaks_PD$peakType %in% c("Promoter",
                                                                           "Distal")),]
  Global_Peaks_PD <- Global_Peaks_PD$ID
  
  Gene_Linked_Peaks_PD <- peaksets[Gene_Linked_Peaks_DEGs,]
  Gene_Linked_Peaks_PD <- Gene_Linked_Peaks_PD[which(Gene_Linked_Peaks_PD$peakType %in% c("Promoter",
                                                                                          "Distal")),]
  Gene_Linked_Peaks_PD <- Gene_Linked_Peaks_PD$ID
  peak_group_PD <- data.frame(peakType = c(rep("Gene-linked peaks", length(Gene_Linked_Peaks_PD)),
                                           rep("Lineage constitutive peaks", length(Lineage_Peaks_PD)),
                                           rep("Global constitutive peaks", length(Global_Peaks_PD)),
                                           rep("Other peaks", length(Other_Peaks_PD))),
                              peakID = c(Gene_Linked_Peaks_PD,
                                         Lineage_Peaks_PD,
                                         Global_Peaks_PD,
                                         Other_Peaks_PD))
  
  peak_group_PD$nearestGene <- peaksets[peak_group_PD$peakID, "nearestGene"]
  peak_group_temp <- as.data.frame(table(peak_group_PD$peakType, peak_group_PD$nearestGene))
  peak_group_temp <- peak_group_temp[which(peak_group_temp$Freq > 0),]
  temp <- as.data.frame.array(table(peak_group_PD$peakType, peak_group_PD$nearestGene))
  temp <- apply(temp, 2, function(x){identical(as.numeric(x[1:3]), c(1,1,1))})
  temp <- names(temp)[which(temp)]
  genes <- names(table(peak_group_temp$Var2))[which(table(peak_group_temp$Var2) %in% c(4))]
  
  temp_genes <- temp[temp %in% genes]
  
  # RNA Lineage DEGS
  Idents(RNA_common_tissues) <- RNA_common_tissues$Mainclusters
  RNA_Lineage_DEGs <- FindAllMarkers(RNA_common_tissues, only.pos = T, logfc.threshold = 0.5)
  
  genes_DF <- RNA_Lineage_DEGs[which(RNA_Lineage_DEGs$gene %in% genes),]
  
  peakSet <- getPeakSet(ATAC_common_tissues)
  Chr <- as.data.frame(seqnames(peakSet))
  Range <- as.data.frame(peakSet@ranges)
  peakSet$Peak_id <- paste0(Chr$value, "_", Range$start, "_", Range$end)
  
  # peak_group_selected <- peak_group[which(peak_group$nearestGene %in% genes_DF$gene),]
  # peak_group_selected <- peak_group[which(peak_group$nearestGene == "OsRLCK85"),]
  peak_group_selected <- peak_group[which(peak_group$nearestGene %in% "SAP11"),]
  
  peaksets[peak_group_selected$peakID,]
  
  Gene_Linked_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Gene-linked peaks"),"peakID"])]
  Global_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Global constitutive peaks"),"peakID"])]
  Lineage_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Lineage constitutive peaks"),"peakID"])]
  Other_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Other peaks"),"peakID"])]
  
  featureList <- SimpleList("Gene-linked peaks" = Gene_Linked_Peaks_PD_RG,
                            "Global constitutive peaks" = Global_Peaks_PD_RG,
                            "Lineage constitutive peaks" = Lineage_Peaks_PD_RG,
                            "Other peaks" = Other_Peaks_PD_RG)
  peak_type_color
  
  
  
  p_Maincluster <- plotBrowserTrack(ArchRProj = ATAC_common_tissues,
                        groupBy = "Mainclusters",
                        geneSymbol = unique(peak_group_selected$nearestGene),
                        upstream = 15000, downstream = 15000,
                        features = featureList,
                        peakPal = peak_type_color,
                        groupPal = Maincluster_color,
                        facetbaseSize = 10)
  dev.off()
  
  # OsMT3a AHA7 Os01g0711400 Os01g0930200 Lip19
  # SAP11
  
  pdf("Figure6_all_genes_genebrowse.pdf", width = 8.5, height = 6.5)
  for (i in unique(genes_DF$gene)) {
    grid::grid.newpage()
    grid::grid.draw(p_Maincluster[[i]])
  }
  dev.off()
  temp_genes
  grid::grid.newpage()
  grid::grid.draw(p_Maincluster$OsMT3a)
  temp <- data.frame(Maincluster = RNA_common_tissues$Mainclusters,
                     expression = RNA_common_tissues@assays$RNA@data["OsBGLU1",])
  ggplot() +
    geom_boxplot(data = temp, aes(x = Maincluster, y = expression, fill = Maincluster)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = Maincluster_color) +
    labs(x = "", y = "RNA Expression")
  dev.off()
  
  pdf("Figure6_5_SAP11_genomeBrowse.pdf", width = 6, height = 5)
  grid::grid.newpage()
  grid::grid.draw(p$SAP11)
  dev.off()
  
  
  temp <- data.frame(Maincluster = RNA_common_tissues$Mainclusters,
                     expression = RNA_common_tissues@assays$RNA@data["SAP11",])
  
  pdf("Figure6_5_SAP11_RNA.pdf", width = 5, height = 3.5)
  ggplot() +
    geom_boxplot(data = temp, aes(x = Maincluster, y = expression, fill = Maincluster)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = Maincluster_color) +
    labs(x = "", y = "RNA Expression")
  dev.off()
  
}


######### peak type characterize
### gene expression level
{
  RNA_mean_all <- rowSums(RNA_common_tissues@assays[["RNA"]]@data) / ncol(RNA_common_tissues@assays[["RNA"]]@data)
  temp <- RNA_mean_all[c(type1_genes, type5_genes, type3_genes, type4_genes, type2_genes)]
  temp <- data.frame(type = c(rep("No GLPs", length(type1_genes)),
                              rep("Only GLPs", length(type5_genes)),
                              rep("GLPs & LCPs", length(type3_genes)),
                              rep("GLPs & GCPs", length(type4_genes)),
                              rep("GLPs & LCPs & GCPs", length(type2_genes))
  ),
  value = temp, check.names = F)
  temp$type <- factor(temp$type, levels = c("No GLPs",
                                            "Only GLPs",
                                            "GLPs & LCPs",
                                            "GLPs & GCPs",
                                            "GLPs & LCPs & GCPs"))
  temp$expression_level <- ""
  temp <- temp[order(temp$value, decreasing = F),]
  temp$rank <- 1:nrow(temp)
  
  pdf("Figure6_5_gene_type_expression_rank_density.pdf", height = 5, width = 7)
  ggplot(data = temp,
         aes(x = rank, group = type, fill = type, color = type)) +
    geom_density(alpha = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "Genes ranked by expression level", y = "Density", fill = "Gene type")
  
  ggplot(data = temp,
         aes(x = value, group = type, color = type)) +
    stat_ecdf(geom = "step")
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "Genes ranked by expression level", y = "Density", fill = "Gene type")
  dev.off()
  
}

{
  RNA_Celltye_DEG
  Meristem_DEG <- RNA_Celltye_DEG[which(RNA_Celltye_DEG$cluster %in% 
                                          c("Inflorescence meristem (IM)",
                                            "Spikelet meristem (SM)",
                                            "Branch meristems (BM)",
                                            "Cryptic bract/bract (cb/b)")),]
  
  Meristem_DEG <- as.data.frame(table(Meristem_DEG$gene))
  Meristem_DEG <- Meristem_DEG[which(Meristem_DEG[,2]>2),]
  
  peaksets_2$peakGroup <- "Other peaks"
  peaksets_2$peakGroup[which(peaksets_2$ID %in% Gene_Linked_Peaks)] <- "Gene-linked peaks"
  peaksets_2$peakGroup[which(peaksets_2$ID %in% Lineage_Peaks)] <- "Lineage constitutive peaks"
  peaksets_2$peakGroup[which(peaksets_2$ID %in% Global_Peaks)] <- "Global constitutive peaks"
  table(peaksets_2$peakGroup)
  
  temp_gene_peak_type <- as.data.frame(table(peaksets_2$nearestGene,
                                             peaksets_2$peakGroup))
  temp_gene_peak_type <- temp_gene_peak_type[which(temp_gene_peak_type$Freq > 0),]
  temp_gene_peak_type <- temp_gene_peak_type[which(temp_gene_peak_type$Var1 %in%
                                                     Meristem_DEG$Var1),]
  main_celltype <- main_celltype[order(main_celltype$main, decreasing = T),]
}

### 不止局限于DEGs
### Peak region percent
{
  peaksets_2_DEG <- peaksets_2
  temp <- as.data.frame(table(peaksets_2_DEG$peakGroup, peaksets_2_DEG$peakType))
  colnames(temp) <- c("group", "peak_region", "Count")
  temp$group <- factor(temp$group,
                       levels = c("Gene-linked peaks",
                                  "Lineage constitutive peaks",
                                  "Global constitutive peaks",
                                  "Other peaks"))
  temp$peak_region <- factor(temp$peak_region,
                             levels = c("Exonic", "Intronic", "Distal", "Promoter"))
  peak_region_color <- c("#ADC698", "#A0C1D1",
                         "#B18ED7", "#D676A1")
  names(peak_region_color) <- c("Exonic", "Intronic", "Distal", "Promoter")
  pdf("Figure6_5_peak_type_region.pdf", width = 4, height = 5.5)
  ggplot() +
    geom_bar(data = temp, aes(x = group, y = Count, fill = peak_region),
             stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(size = 10)) +
    labs(x = "", y = "Percentage", fill = "Peak Region") +
    scale_fill_manual(values = peak_region_color)
  dev.off()

}

### Peak conserve score
{
  peaks_group <- peaksets_2_DEG
  chr <- unlist(lapply(peaksets_2$ID, function(x){
    unlist(strsplit(x, "_"))[1]
  }))
  chr <- paste0("Chr", chr)
  start <- unlist(lapply(peaksets_2$ID, function(x){
    unlist(strsplit(x, "_"))[2]
  }))
  end <- unlist(lapply(peaksets_2$ID, function(x){
    unlist(strsplit(x, "_"))[3]
  }))
  peaks_bed <- data.frame(chr = chr,
                          start = start,
                          end = end,
                          meta = peaksets_2$ID)
  write.table(peaks_bed, "peaks_bed.bed", sep = "\t", row.names = F, col.names = F, quote = F)
  
  peak_mean_PhastCons <- read.table("peak_meanscore.tab",
                                    sep = "\t")
  rownames(peak_mean_PhastCons) <- peak_mean_PhastCons$V1
  peaks_group$Mean_PhastCons <- peak_mean_PhastCons[peaks_group$ID,5]
  sum(is.na(peaks_group$Mean_PhastCons))
  peaks_group_PD <- peaks_group[which(peaks_group$peakType %in% c("Promoter", "Distal")),]
  
  pdf("Figure6_5_Mean_PhastCons_diff.pdf", width = 6, height = 6)
  ggplot() +
    geom_violin(data = peaks_group,
                aes(x = peakGroup, y = Mean_PhastCons, fill = peakGroup, color = peakGroup)) +
    geom_boxplot(data = peaks_group,
                 aes(x = peakGroup, y = Mean_PhastCons),
                 width = 0.3, alpha = 0.6, linewidth = 0.8,
                 outlier.size = 0.7) +
    stat_summary(data = peaks_group,
                 aes(x = peakGroup, y = Mean_PhastCons, fill = peakGroup),
                 fun = "mean", geom = "point", shape = 18, size = 2.5,
                 color = "black", fill = "black", alpha = 0.7) +
    stat_compare_means(data = peaks_group,
                       aes(x = peakGroup, y = Mean_PhastCons)) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = peak_type_color) +
    scale_color_manual(values = peak_type_color) +
    labs(x = "", y = "Mean PhastCons Score")
  dev.off()
  
  pdf("Figure6_5_Mean_PhastCons_diff_PD.pdf", width = 6, height = 6)
  ggplot() +
    geom_violin(data = peaks_group_PD,
                aes(x = peakGroup, y = Mean_PhastCons, fill = peakGroup, color = peakGroup)) +
    geom_boxplot(data = peaks_group_PD,
                 aes(x = peakGroup, y = Mean_PhastCons),
                 width = 0.3, alpha = 0.6, linewidth = 0.8,
                 outlier.size = 0.7) +
    stat_summary(data = peaks_group_PD,
                 aes(x = peakGroup, y = Mean_PhastCons, fill = peakGroup),
                 fun = "mean", geom = "point", shape = 18, size = 2.5,
                 color = "black", fill = "black", alpha = 0.7) +
    stat_compare_means(data = peaks_group_PD,
                       aes(x = peakGroup, y = Mean_PhastCons)) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = peak_type_color) +
    scale_color_manual(values = peak_type_color) +
    labs(x = "", y = "Mean PhastCons Score")
  dev.off()
  
  CE <- read.csv("./Ref/Conserved elements (CEs, comparative genomics)_Osj.gtf",
                 header = F, sep = "\t")
  table(CE$V1)
  CE <- CE[-which(CE$V1 %in% c("ChrSy", "ChrUn")),]
  chr <- data.frame(chr = 1:12,
                    row.names = paste0("Chr",1:12))
  CE$V1 <- chr[CE$V1,1]
  CE <- makeGRangesFromDataFrame(CE, keep.extra.columns = TRUE,
                                 seqnames.field = "V1",
                                 start.field = "V4",
                                 end.field = "V5",
                                 strand.field = "V7")
  ATAC_common_tissues <- addPeakAnnotations(ATAC_common_tissues,
                                            list(CE),
                                            name = "CE",
                                            force = TRUE)
  CE_peak_anno <- getPeakAnnotation(ATAC_common_tissues, name = "CE")
  CE_peak_anno_match <- readRDS(CE_peak_anno[["Matches"]])
  CE_peak_anno_match <- CE_peak_anno_match@assays@data@listData[["matches"]]
  CE_peak_anno_match <- as.matrix(CE_peak_anno_match)
  CE_peak_anno_match <- CE_peak_anno_match + 0
  CE_peak_anno_match <- as.data.frame(CE_peak_anno_match)
  rownames(CE_peak_anno_match) <- rownames(peaksets)
  
  peaks_group$CE <- CE_peak_anno_match[peaks_group$ID, 1]
  
  peaks_group_PD <- peaks_group[which(peaks_group$peakType %in% c("Promoter", "Distal")),]
  
  CE_DF <- as.data.frame(table(peaks_group_PD$peakGroup, peaks_group_PD$CE))
  group_num <- as.data.frame(table(peaks_group_PD$peakGroup))
  rownames(group_num) <- group_num$Var1
  CE_DF$Freq <- CE_DF$Freq / group_num[CE_DF$Var1,2]
  
  CE_DF$Var2 <- ifelse(CE_DF$Var2 == 0, "Not-matched CE", "Matched CE")
  
  pdf("Figure6_5_conserved_element_percent.pdf", height = 7, width = 6.5)
  ggplot() +
    geom_bar(data = CE_DF, aes(x = Var2, y = Freq, fill = Var1),
             stat = "identity", position = position_dodge(), width = 0.8) +
    theme_bw(base_size = 10) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = peak_type_color) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 10),
          legend.title = element_text(color = "black", size = 10)) +
    labs(x = "", y = "Percentage", fill = "Peak Type")
  dev.off()
  
  CE_DF_temp <- as.data.frame.array(table(peaks_group_PD$peakGroup, peaks_group_PD$CE))
  
  chisq.test(CE_DF_temp[c(1,4),], correct = F) # GLPs & Other
  chisq.test(CE_DF_temp[c(2,4),], correct = F) # GCPs & Other
  chisq.test(CE_DF_temp[c(3,4),], correct = F) # LCPs & Other
  chisq.test(CE_DF_temp[c(4,4),], correct = F) # GLPs & Other
}

#### JSD
{
  
  sum(colnames(peak_matrix_log) %in% ATAC_common_tissues$cellNames)
  
  peak_celltype_mean_scaled <- c()
  for (i in unique(ATAC_common_tissues$Celltype)) {
    print(i)
    cells <- getCellNames(ATAC_common_tissues)
    cells <- cells[which(ATAC_common_tissues$Celltype == i)]
    temp <- peak_matrix_log[,cells]
    temp <- rowSums(temp) / ncol(temp)
    peak_celltype_mean_scaled <- as.data.frame(cbind(peak_celltype_mean_scaled,
                                                     temp))
  }
  colnames(peak_celltype_mean_scaled) <- unique(ATAC_common_tissues$Celltype)
  rm(temp)
  sums <- rowSums(peak_celltype_mean_scaled)
  peak_celltype_mean_scaled <- peak_celltype_mean_scaled / sums
  
  JSD <- c()
  for (i in 1:nrow(peak_celltype_mean_scaled)) {
    print(i)
    temp <- data.frame(x = t(peak_celltype_mean_scaled[i,])[,1], y = rep(1/13, 13))
    temp <- as.matrix(temp)
    temp <- philentropy::JSD(t(temp))
    JSD <- c(JSD, temp)
  }
  names(JSD) <- rownames(peak_celltype_mean_scaled)
  
  peaks_group$JSD <- JSD[peaks_group$ID]
  pdf("Figure6_5_JSD_peakGroup_2.pdf", width = 5, height = 5.5)
  ggplot(data = peaks_group, aes(x = peakGroup, y = JSD, fill = peakGroup)) +
    geom_boxplot(data = peaks_group, aes(x = peakGroup, y = JSD, fill = peakGroup),
                 outlier.shape = NA) +
    coord_cartesian(ylim = c(0, 0.3)) +
    scale_fill_manual(values = peak_type_color) +
    labs(fill = "PeakType", x = "", y = "Cell-type specific score (JSD)") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggpubr::stat_compare_means(comparisons = list(c("Gene-linked peaks", "Global constitutive peaks"),
                                                c("Gene-linked peaks", "Lineage constitutive peaks"),
                                                c("Other peaks", "Lineage constitutive peaks"),
                                                c("Other peaks", "Global constitutive peaks")),
                             label = "p.signif",
                             label.y = c(0.1, 0.15, 0.2, 0.25))
  dev.off()
  
}

### tissue insertion count
{
  unique(ATAC_common_tissues$TissuesSub)
  temp_tissue <- c()
  peak_type_color
  for (i in unique(ATAC_common_tissues$TissuesSub)) {
    print(i)
    cells <- getCellNames(ATAC_common_tissues)
    cells <- cells[which(ATAC_common_tissues$TissuesSub == i)]
      temp <- peak_matrix_log[Gene_Linked_Peaks_DEGs,cells]
      temp <- rowSums(temp) / ncol(temp)
      temp_tissue <- as.data.frame(rbind(temp_tissue,
                                         data.frame(tissue = i,
                                                    peakGroup = "Gene-linked peaks",
                                                    value = temp)))
      temp <- peak_matrix_log[Lineage_Peaks_DEGs,cells]
      temp <- rowSums(temp) / ncol(temp)
      temp_tissue <- as.data.frame(rbind(temp_tissue,
                                         data.frame(tissue = i,
                                                    peakGroup = "Lineage constitutive peaks",
                                                    value = temp)))
      temp <- peak_matrix_log[Global_Peaks_DEGs,cells]
      temp <- rowSums(temp) / ncol(temp)
      temp_tissue <- as.data.frame(rbind(temp_tissue,
                                         data.frame(tissue = i,
                                                    peakGroup = "Global constitutive peaks",
                                                    value = temp)))
      temp <- peak_matrix_log[Other_Peaks_DEGs,cells]
      temp <- rowSums(temp) / ncol(temp)
      temp_tissue <- as.data.frame(rbind(temp_tissue,
                                         data.frame(tissue = i,
                                                    peakGroup = "Other peaks",
                                                    value = temp)))
    
  }
  
  
  pdf("Figure6_5_peakType_tissues_boxplot.pdf", width = 7, height = 5.5)
  ggplot(data = temp_tissue,
         aes(x = tissue, y = value, fill = peakGroup)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = peak_type_color) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "", fill = "PeakType", y = "Average Log-normalized Insertion Count")
  dev.off()
  
  unique(temp_tissue$tissue)
  ggplot(data = temp_tissue[which(temp_tissue$tissue == "IM1cm"),],
         aes(x = peakGroup, y = value, fill = peakGroup)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = peak_type_color) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "", fill = "PeakType", y = "Average Log-normalized Insertion Count") +
    ggpubr::stat_compare_means(comparisons = list(c("Gene-linked peaks", "Global constitutive peaks"),
                                                  c("Gene-linked peaks", "Lineage constitutive peaks"),
                                                  c("Other peaks", "Lineage constitutive peaks"),
                                                  c("Other peaks", "Global constitutive peaks")),
                               label = "p.signif")
  
}

#### peak TFBS
{
# 
#   peaks_group_PD <- peak_group[which(peak_group$peakRegion %in% c("Promoter", "Distal")),]
# 
#   peak_group$CE <- CE_peak_anno_match[peak_group$peakID, 1]
# 
#   peaks_group_PD <- peak_group[which(peak_group$peakRegion %in% c("Promoter", "Distal")),]
#   peaks_group_PD$CE_type <- ifelse(peaks_group_PD$CE == 0, "Not-matched CE", "Matched CE")
# 
#   CE_DF <- as.data.frame(table(peaks_group_PD$peakType, peaks_group_PD$CE_type))
#   group_num <- as.data.frame(table(peaks_group_PD$peakType))
#   rownames(group_num) <- group_num$Var1
#   CE_DF$Freq <- CE_DF$Freq / group_num[CE_DF$Var1,2]
# 
#   pdf("Figure6_3_conserved_element_percent_PD.pdf", height = 7, width = 5.5)
#   ggplot() +
#     geom_bar(data = CE_DF, aes(x = Var2, y = Freq, fill = Var1),
#              stat = "identity", position = position_dodge(), width = 0.8) +
#     theme_bw(base_size = 10) +
#     scale_y_continuous(labels = scales::percent) +
#     scale_fill_manual(values = peak_type_color) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#           axis.text = element_text(color = "black", size = 10),
#           axis.title = element_text(color = "black", size = 10)) +
#     labs(x = "", y = "Percent")
#   dev.off()
  # TF Motif
  MSU_RAP <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RAP-MSU_2023-03-15.txt",
                      sep = "\t", header = F)
  MSU_RAP <- MSU_RAP[-which(MSU_RAP$V2 == "None"),]
  MSU_RAP <- MSU_RAP[-which(MSU_RAP$V1 == "None"),]
  MSU <- c()
  row <- c()
  for (i in 1:nrow(MSU_RAP)) {
    temp <- unlist(strsplit(MSU_RAP[i,2], ",", fixed = T))
    row <- c(row, rep(i, length(temp)))
    MSU <- c(MSU, temp)
  }
  MSU_RAP <- data.frame(RAP = MSU_RAP[row,1],
                        MSU = MSU)
  MSU_RAP$MSU <- unlist(lapply(MSU_RAP$MSU, function(x){
    unlist(strsplit(x, ".", fixed = T))[1]
  }))
  MSU_RAP <- MSU_RAP[!duplicated(MSU_RAP),]
  
  motif <- read_meme("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs.meme")
  motif_all <- convert_motifs(motif, class = "TFBSTools-PWMatrix")
  motif_all <- do.call(PWMatrixList, motif_all)
  MSU_ID <- unlist(lapply(motif_all@listData, function(x){
    x@name
  }))
  names(motif_all) <- MSU_ID
  RAP_ID <- MSU_RAP[match(MSU_ID, MSU_RAP$MSU),]
  RAP_ID <- na.omit(RAP_ID)
  motif_all <- motif_all[RAP_ID$MSU]
  names(motif_all) <- RAP_ID$RAP
  motif_all <- motif_all[!duplicated(names(motif_all))] # ID duplicated !!!
  rownames(gene_id_map) <- gene_id_map$gene_id
  motif_all <- motif_all[names(motif_all) %in% gene_id_map$gene_id]
  Symbol <- match(names(motif_all), gene_id_map$gene_id)
  Symbol <- gene_id_map[Symbol,]
  Annotation_genes <- read.csv("./Ref/IRGSP-1.0_representative_annotation_2023-03-15.tsv",
                               sep = "\t", header = T, check.names = F)
  Annotation_genes <- Annotation_genes[,c("Locus_ID", "CGSNL Gene Symbol", "CGSNL Gene Name")]
  Annotation_genes <- Annotation_genes[Annotation_genes$Locus_ID %in% names(motif_all),]
  Annotation_genes <- Annotation_genes[!duplicated(Annotation_genes$Locus_ID),]
  Annotation_genes$`CGSNL Gene Symbol`[168] <- "EREBP96"
  Annotation_genes$`CGSNL Gene Symbol`[23] <- "EREBP127"
  Annotation_genes$`CGSNL Gene Symbol`[244] <- "RID1"
  Annotation_genes$`CGSNL Gene Symbol` <- paste0("Os",Annotation_genes$`CGSNL Gene Symbol`)
  empty_row <- which(Annotation_genes$`CGSNL Gene Symbol` == "Os")
  Annotation_genes$`CGSNL Gene Symbol`[empty_row] <- Annotation_genes$Locus_ID[empty_row]
  Annotation_genes$DataSets_Symbol <- gene_id_map[Annotation_genes$Locus_ID, "symbol"]
  rownames(Annotation_genes) <- Annotation_genes$Locus_ID
  names(motif_all) <- Annotation_genes[names(motif_all), "CGSNL Gene Symbol"]
  
  ATAC_common_tissues <- addMotifAnnotations(ArchRProj = ATAC_common_tissues, motifPWMs = motif_all,
                                             annoName = "TF-Motif", force = T)
  matches <- getMatches(ATAC_common_tissues, name = "TF-Motif")
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  
  matches <- matches@assays@data@listData[["matches"]]
  rownames(matches) <- pr1
  matches <- matches + 0
  matches <- as.matrix(matches)
  TFBS <- rowSums(matches)
  quantile(TFBS)
  
  peaks_group_PD$TFBS <- TFBS[peaks_group_PD$ID]
  # peaks_group_PD$CE_type <- factor(peaks_group_PD$CE_type, levels = c("Matched CE", "Not-matched CE"))
  temp <- lapply(c("Other peaks", "Lineage constitutive peaks", "Gene-linked peaks", "Global constitutive peaks"), function(x){
    temp <- peaks_group_PD[which(peaks_group_PD[,"peakGroup"] == as.character(x)),]
    quan <- quantile(temp[, "TFBS"])
    IQR <- quan[4] - quan[2]
    Q1 <- quan[2]
    Q2 <- quan[3]
    Q3 <- quan[4]
    low <- Q1-1.5*IQR
    if (low <= min(temp[, "TFBS"])) {
      low <- min(temp[, "TFBS"])
    }
    high <- Q3+1.5*IQR
    if (high >= max(temp[, "TFBS"])) {
      high <- max(temp[, "TFBS"])
    }
    data.frame(group = as.character(x),
               low = low,
               Q1 = Q1,
               Q2 = Q2,
               Q3 = Q3,
               high = high
    )
  })
  temp <- rlist::list.rbind(temp)
  
  pdf("Figure6_5_PD_TFBS.pdf", height = 6, width = 5.5)
  ggplot() +
    geom_violin(data = peaks_group_PD, aes(x = peakGroup, y = TFBS, fill = peakGroup)) +
    scale_fill_manual(values = peak_type_color) +
    geom_segment(aes(x = group, y = low, xend = group, yend = high), linewidth = 0.5, colour = "black", data = temp) +
    geom_segment(aes(x = group, y = Q1, xend = group, yend = Q3), linewidth = 2, colour = "black", data = temp) +
    geom_point(aes(x = group, y = Q2), data = temp, colour = "white", size = 2) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    stat_compare_means(comparisons = list(c(names(peak_type_color)[3], names(peak_type_color)[1]),
                                          c(names(peak_type_color)[2], names(peak_type_color)[1]),
                                          # c(names(peak_type_color)[3], names(peak_type_color)[2]),
                                          c(names(peak_type_color)[3], names(peak_type_color)[4]),
                                          c(names(peak_type_color)[2], names(peak_type_color)[4])), label = "p.signif") +
    labs(x = "", y = "Number of TF binding sites", fill = "PeakType")
  
  dev.off()
}


#### TF enrichment
{
  Other_Peaks_PD <- peaksets[Other_Peaks,]
  Other_Peaks_PD <- Other_Peaks_PD[which(Other_Peaks_PD$peakType %in% c("Promoter",
                                                                        "Distal")),]
  Other_Peaks_PD <- Other_Peaks_PD$ID
  
  Lineage_Peaks_PD <- peaksets[Lineage_Peaks,]
  Lineage_Peaks_PD <- Lineage_Peaks_PD[which(Lineage_Peaks_PD$peakType %in% c("Promoter",
                                                                              "Distal")),]
  Lineage_Peaks_PD <- Lineage_Peaks_PD$ID
  
  Global_Peaks_PD <- peaksets[Global_Peaks,]
  Global_Peaks_PD <- Global_Peaks_PD[which(Global_Peaks_PD$peakType %in% c("Promoter",
                                                                           "Distal")),]
  Global_Peaks_PD <- Global_Peaks_PD$ID
  
  Gene_Linked_Peaks_PD <- peaksets[Gene_Linked_Peaks,]
  Gene_Linked_Peaks_PD <- Gene_Linked_Peaks_PD[which(Gene_Linked_Peaks_PD$peakType %in% c("Promoter",
                                                                                          "Distal")),]
  Gene_Linked_Peaks_PD <- Gene_Linked_Peaks_PD$ID
  
  
  peak_type_color
  Global_Peaks_PD
  Lineage_Peaks_PD
  Gene_Linked_Peaks_PD
  Other_Peaks_PD
  
  
  ATAC_common_tissues <- addBgdPeaks(ATAC_common_tissues, force = T)
  bgdPeaks <- getBgdPeaks(ATAC_common_tissues)
  bgdPeaks_range <- rowRanges(bgdPeaks)
  bgdPeaks_peak <- paste(seqnames(bgdPeaks_range),
                         start(bgdPeaks_range),
                         end(bgdPeaks_range),
                         sep="_")
  bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ATAC_common_tissues))
  rownames(bgdPeaks) <- bgdPeaks_peak
  sum(bgdPeaks_peak %in% peaksets$ID) == length(bgdPeaks_peak)
  
  
  matches <- getMatches(ATAC_common_tissues, name = "TF-Motif")
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  rownames(matches) <- pr1
  
  
  peaksets_2_pd <- peaksets_2[which(peaksets_2$peakType %in% c("Promoter",
                                                               "Distal")),]
  
  sum(Gene_Linked_Peaks_PD %in% peaksets_2_pd$ID) == length(Gene_Linked_Peaks_PD)
  sum(Lineage_Peaks_PD %in% peaksets_2_pd$ID) == length(Lineage_Peaks_PD)
  sum(Global_Peaks_PD %in% peaksets_2_pd$ID) == length(Global_Peaks_PD)
  sum(Other_Peaks_PD %in% peaksets_2_pd$ID) == length(Other_Peaks_PD)
  (length(Gene_Linked_Peaks_PD) + length(Lineage_Peaks_PD) + length(Global_Peaks_PD) + length(Other_Peaks_PD)) == nrow(peaksets_2_pd)
  
  temp <- matrix(FALSE, nrow = nrow(peaksets_2_pd), ncol = 4)
  rownames(temp) <- peaksets_2_pd$ID
  colnames(temp) <- names(peak_type_color)
  temp[Gene_Linked_Peaks_PD, 1] <- TRUE
  temp[Lineage_Peaks_PD, 2] <- TRUE
  temp[Global_Peaks_PD, 3] <- TRUE
  temp[Other_Peaks_PD, 4] <- TRUE
  
  matches <- matches[which(rownames(matches) %in% rownames(temp)),]
  matches <- matches[rownames(temp),]
  identical(rownames(temp), rownames(matches))

  Other_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,4]), 1:nrow(matches))
  Other_Peaks_Enrich$Enrichment_log <- log2(Other_Peaks_Enrich$Enrichment)
  Other_Peaks_Enrich$Enrichmented <- "NO"
  Other_Peaks_Enrich$Enrichmented[which(Other_Peaks_Enrich$Enrichment_log >= 0.25 & Other_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Other_Peaks_Enrich$Enrichmented)
  
  colnames(temp)
  
  Global_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,3]), 1:nrow(matches))
  Global_Peaks_Enrich$Enrichment_log <- log2(Global_Peaks_Enrich$Enrichment)
  Global_Peaks_Enrich$Enrichmented <- "NO"
  Global_Peaks_Enrich$Enrichmented[which(Global_Peaks_Enrich$Enrichment_log >= 0.25 & Global_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Global_Peaks_Enrich$Enrichmented)
  
  Lineage_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,2]), 1:nrow(matches))
  Lineage_Peaks_Enrich$Enrichment_log <- log2(Lineage_Peaks_Enrich$Enrichment)
  Lineage_Peaks_Enrich$Enrichmented <- "NO"
  Lineage_Peaks_Enrich$Enrichmented[which(Lineage_Peaks_Enrich$Enrichment_log >= 0.25 & Lineage_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Peaks_Enrich$Enrichmented)
  
  Gene_Linked_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,1]), 1:nrow(matches))
  Gene_Linked_Peaks_Enrich$Enrichment_log <- log2(Gene_Linked_Peaks_Enrich$Enrichment)
  Gene_Linked_Peaks_Enrich$Enrichmented <- "NO"
  Gene_Linked_Peaks_Enrich$Enrichmented[which(Gene_Linked_Peaks_Enrich$Enrichment_log >= 0.25 & Gene_Linked_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Gene_Linked_Peaks_Enrich$Enrichmented)
  
  gplot1 <- function(DATA) {
    DATA$Enrichmented <- factor(DATA$Enrichmented,
                                levels = c("YES", 'NO'))
    temp <- NULL
    if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") > 10) {
      temp <- DATA[which(DATA$Enrichmented == "YES"),]
      temp <- temp[order(temp$Enrichment_log, decreasing = T),]
      temp <- temp[1:10,]
    }
    if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") <= 10) {
      temp <- DATA[which(DATA$Enrichmented == "YES"),]
      temp <- temp[order(temp$Enrichment_log, decreasing = T),]
    }
    
    g <- ggplot() +
      geom_point(data = DATA, aes(x = Enrichment_log, y = mlog10Padj, color = Enrichmented)) +
      theme_bw() +
      theme(axis.title = element_text(size = 10, colour = "black"),
            axis.text = element_text(size = 10, color = "black")) +
      labs(y = "-Log10(adjusted P-value)", x = "Log2(Enrichment)") +
      geom_hline(yintercept = 5) +
      geom_vline(xintercept = 0.25)
    if (is.null(temp)) {
      g <- g + 
        scale_color_manual(values = c("gray"))
    }
    if (!is.null(temp)) {
      g <- g +
        scale_color_manual(values = c("red", "gray")) + 
        geom_text_repel(data = temp, aes(x = Enrichment_log, y = mlog10Padj, 
                                            label = temp$feature),
                      size = 3,
                      box.padding = unit(0.8, "lines"),
                      point.padding = unit(0, "lines"), 
                      min.segment.length = 0,
                      segment.color = "black",
                      colour="#000000",
                      show.legend = FALSE,
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
    }
    
    if (max(DATA$Enrichment_log) < 0.25) {
      g <- g +
        scale_x_continuous(limits = c(min(DATA$Enrichment_log), 0.5))
    }
    g <- g + theme(legend.position = "none")
    return(g)
  }
  
  Other_Peaks_Enrich_ggplot <- gplot1(Other_Peaks_Enrich)
  Other_Peaks_Enrich_ggplot
  Gene_Linked_Peaks_Enrich_ggplot <- gplot1(Gene_Linked_Peaks_Enrich)
  Gene_Linked_Peaks_Enrich_ggplot
  Lineage_Peaks_Enrich_ggplot <- gplot1(Lineage_Peaks_Enrich)
  Lineage_Peaks_Enrich_ggplot
  Global_Peaks_Enrich_ggplot <- gplot1(Global_Peaks_Enrich)
  Global_Peaks_Enrich_ggplot
  
  pdf("Figure6_5_PeakGroup_TF_enrich.pdf", width = 13.33, height = 4)
  Global_Peaks_Enrich_ggplot +
    Lineage_Peaks_Enrich_ggplot +
    Other_Peaks_Enrich_ggplot +
    Gene_Linked_Peaks_Enrich_ggplot +
    plot_layout(nrow = 1, byrow = FALSE)
  dev.off()
  
  library(gplots)
  library(VennDiagram)
  venn_list_2 <- list("Gene-linked peaks\nenriched TFs" = Gene_Linked_Peaks_Enrich[which(Gene_Linked_Peaks_Enrich$Enrichmented == "YES"),"feature"],
                    "Global constitutive peaks\nenriched TFs" = Global_Peaks_Enrich[which(Global_Peaks_Enrich$Enrichmented == "YES"),"feature"],
                    "Lineage constitutive peaks\nenriched TFs" = Lineage_Peaks_Enrich[which(Lineage_Peaks_Enrich$Enrichmented == "YES"),"feature"])
  pdf("Figure6_5_TF_enrich_venn.pdf", width = 5, height = 5)
  venn.plot <- venn.diagram(
    x = venn_list_2,
    euler.d = TRUE,
    cex = 2.5,
    cat.cex = 2.5,
    cat.pos = 0,
    filename = NULL
  )
  grid.draw(venn.plot)
  
  
  Tvenn <- venn.diagram(venn_list_2,
                        filename=NULL,
                        lwd=1,lty=2,
                        col = c(peak_type_color[c(1,3,2)]) ,
                        fill = c(peak_type_color[c(1,3,2)]),
                        cat.col = c(peak_type_color[c(1,3,2)]),
                        reverse = TRUE)
  grid.draw(Tvenn)
  
  dev.off()
  
  global_tf <- Global_Peaks_Enrich[which(Global_Peaks_Enrich$Enrichmented == "YES"),"feature"]
  lineage_tf <- Lineage_Peaks_Enrich[which(Lineage_Peaks_Enrich$Enrichmented == "YES"),"feature"]
  setdiff(global_tf, lineage_tf)
  gene_linked_tf <- Gene_Linked_Peaks_Enrich[which(Gene_Linked_Peaks_Enrich$Enrichmented == "YES"),"feature"]
}

#### TF family
{
  TF_Family <- read.csv("./Ref/Rice_TF_Family.txt", sep = "\t", header = F)
  TF_Family <- as.data.frame(na.omit(TF_Family))
  colnames(TF_Family) <- c("MSU", "TF_Family")
  TF_Family <- TF_Family[!duplicated(TF_Family),]
  TF_Family$MSU <- unlist(lapply(TF_Family$MSU, function(x){
    unlist(strsplit(x, ".", fixed = T))[1]
  }))
  TF_Family$TF_Family <- unlist(lapply(TF_Family$TF_Family, function(x){
    unlist(strsplit(x, " ", fixed = T))[1]
  }))
  sort(unique(TF_Family$TF_Family))
  length(unique(TF_Family$TF_Family))
  
  TF_Family$RAP_ID <- MSU_RAP$RAP[match(TF_Family$MSU, MSU_RAP$MSU)]
  TF_Family <- as.data.frame(na.omit(TF_Family))
  rownames(gene_id_map) <- gene_id_map$gene_id
  TF_Family$symbol <- gene_id_map[TF_Family$RAP_ID, "symbol"]
  TF_Family[which(is.na(TF_Family$symbol)), "symbol"] <- TF_Family[which(is.na(TF_Family$symbol)), "RAP_ID"]
  TF_Family <- TF_Family[!duplicated(TF_Family$RAP_ID),]
  
  
  Annotation_genes$MSU <- MSU_RAP$MSU[match(Annotation_genes$Locus_ID, MSU_RAP$RAP)]
  Annotation_genes$Gene_family <- TF_Family$TF_Family[match(Annotation_genes$MSU,
                                                            TF_Family$MSU)]
  rownames(Annotation_genes) <- Annotation_genes$`CGSNL Gene Symbol`
  
  global_only_tf <- setdiff(global_tf, gene_linked_tf)
  
  global_only_tf <- table(Annotation_genes[global_only_tf, "Gene_family"])
  global_only_tf <- as.data.frame(global_only_tf)
  colnames(global_only_tf) <- c("TF Family", "Count")
  global_only_tf$Group <- "A"
  global_only_tf$Name <- paste0(global_only_tf$`TF Family`, " (", global_only_tf$Count, ")")
  global_only_tf$Name2 <- paste0(global_only_tf$`TF Family`, " (", global_only_tf$Count, " TFs)")
  
  TF_total_number <- 1750
  P <- c()
  for (i in 1:nrow(global_only_tf)) {
    common_num <- global_only_tf[i, "Count"]
    p <- phyper(common_num-1, sum(global_only_tf[,"Count"]), TF_total_number-sum(global_only_tf[,"Count"]), as.numeric(table(TF_Family$TF_Family)[global_only_tf[i, "TF Family"]]), lower.tail=F)
    P <- c(P, p)
  }
  
  global_only_tf$P <- P
  global_only_tf$Sig_if <- ""
  for (i in 1:nrow(global_only_tf)) {
    if (global_only_tf$P[i] < 0.05) {
      global_only_tf$Sig_if[i] <- "*"
    }
    if (global_only_tf$P[i] < 0.01) {
      global_only_tf$Sig_if[i] <- "**"
    }
    if (global_only_tf$P[i] < 0.001) {
      global_only_tf$Sig_if[i] <- "***"
    }
    if (global_only_tf$P[i] < 0.0001) {
      global_only_tf$Sig_if[i] <- "****"
    }
  }
  for (i in 1:nrow(global_only_tf)) {
    if (global_only_tf$P[i] < 0.05) {
      global_only_tf$Name[i] <- paste0(global_only_tf$Name[i], "\n", global_only_tf$Sig_if[i])
    }
  }
  
  global_only_tf$P2 <- -log10(global_only_tf$P)
  global_only_tf$P2[which(global_only_tf$P2 == Inf)] <- 15
  library(treemapify)
  pdf("Figure6_5_global_tf_family.pdf", width = 8, height = 5)
  ggplot(global_only_tf, aes(area = Count, fill = P2, label = Name)) +
    geom_treemap() +
    geom_treemap_text(colour = "black", place = "centre",
                      grow = FALSE, min.size = 5, padding.x = grid::unit(4, "mm"),
                      padding.y = grid::unit(4, "mm")) +
    scale_fill_gradientn(colours = colorRampPalette(c("#EEE9E9", "#EEB4B4", "#FF6A6A", "#EE3B3B"))(20))
  dev.off()
  
  global_only_tf <- global_only_tf[order(global_only_tf$P2, decreasing = T),]
  global_only_tf$P2[2] <- global_only_tf$P2[2] + 1
  global_only_tf <- global_only_tf[order(global_only_tf$P2, decreasing = T),]
  global_only_tf$`TF Family` <- factor(global_only_tf$`TF Family`, levels = rev(global_only_tf$`TF Family`))
  global_only_tf$Name2 <- factor(global_only_tf$Name2, levels = rev(global_only_tf$Name2))
  pdf("Figure6_5_global_tf_family_barplot.pdf", width = 4, height = 4)
  ggplot(global_only_tf, aes(x = P2, y = Name2)) +
    geom_bar(stat = "identity", fill = "#7AC5CD") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black")) +
    geom_vline(xintercept = -log10(0.05), color = "red") +
    labs(x = "-Log10(P-value)", y = "TF family")
  dev.off()
  
  
  # pdf("Figure6_3_global_only_tf_family.pdf", width = 6, height = 4)
  # treemap::treemap(global_only_tf, #数据
  #         index = "Name",#分类变量
  #         vSize = "Count",#分类变量对应数据值
  #         vColor="P",#颜色深浅的对应列
  #         type = "value",#"颜色映射方式,"index"、"value"、"comp"、"dens"、"depth"、"categorical"、"color"、"manual"
  #         title = '79 global constitutive peaks enriched TFs',#标题
  #         border.col = "black",#边框颜色
  #         border.lwds = 2,#边框线宽度
  #         fontsize.labels = c(10),#标签大小
  #         fontcolor.labels = 'black',#标签颜色
  #         align.labels = list(c("center", "center")),#标签位置
  #         fontface.labels = 1, #标签字体：1,2,3,4 表示正常、粗体、斜体、粗斜体
  #         lowerbound.cex.labels = .3,
  #         palette = RColorBrewer::brewer.pal(n = 5, "RdYlBu"))
  # dev.off()
  #
  
  gene_linked_tf <- table(Annotation_genes[gene_linked_tf, "Gene_family"])
  gene_linked_tf <- as.data.frame(gene_linked_tf)
  colnames(gene_linked_tf) <- c("TF Family", "Count")
  gene_linked_tf$Group <- "A"
  gene_linked_tf$Name <- paste0(gene_linked_tf$`TF Family`, " (", gene_linked_tf$Count, ")")
  
  TF_total_number <- 1750
  P <- c()
  for (i in 1:nrow(gene_linked_tf)) {
    common_num <- gene_linked_tf[i, "Count"]
    p <- phyper(common_num-1, sum(gene_linked_tf[,"Count"]), TF_total_number-sum(gene_linked_tf[,"Count"]), as.numeric(table(TF_Family$TF_Family)[global_only_tf[i, "TF Family"]]), lower.tail=F)
    P <- c(P, p)
  }
  
  gene_linked_tf$P <- P
  gene_linked_tf$Sig_if <- ""
  for (i in 1:nrow(gene_linked_tf)) {
    if (gene_linked_tf$P[i] < 0.05) {
      gene_linked_tf$Sig_if[i] <- "*"
    }
    if (gene_linked_tf$P[i] < 0.01) {
      gene_linked_tf$Sig_if[i] <- "**"
    }
    if (gene_linked_tf$P[i] < 0.001) {
      gene_linked_tf$Sig_if[i] <- "***"
    }
    if (gene_linked_tf$P[i] < 0.0001) {
      gene_linked_tf$Sig_if[i] <- "****"
    }
  }
  for (i in 1:nrow(gene_linked_tf)) {
    if (gene_linked_tf$P[i] < 0.05) {
      gene_linked_tf$Name[i] <- paste0(gene_linked_tf$Name[i], "\n", gene_linked_tf$Sig_if[i])
    }
  }
  
  gene_linked_tf$P2 <- -log10(gene_linked_tf$P)
  library(treemapify)
  pdf("Figure6_5_gene_linked_tf_family.pdf", width = 8, height = 5)
  ggplot(gene_linked_tf, aes(area = Count, fill = P2, label = Name)) +
    geom_treemap() +
    geom_treemap_text(colour = "black", place = "centre",
                      grow = TRUE, min.size = 1, padding.x = grid::unit(4, "mm"),
                      padding.y = grid::unit(4, "mm")) +
    scale_fill_gradient(low = "#DEF5E5FF", high = "#440154FF")
  dev.off()
  
  
}

##### Global TF GO
{
  DAVID_MF <- read.csv("./Ref/DAVIDKnowledgebase/ENSEMBL_GENE_ID2GOTERM_MF_ALL.txt", sep = "\t", header = F)
  unique(str_count(DAVID_MF$V2,"~")) 
  DAVID_MF$Ontology <- "MF"
  DAVID_MF$GO_term <- unlist(lapply(DAVID_MF$V2, function(x){
    unlist(strsplit(x, "~"))[1]
  }))
  DAVID_MF$Description <- unlist(lapply(DAVID_MF$V2, function(x){
    unlist(strsplit(x, "~"))[2]
  }))
  DAVID_BP <- read.csv("./Ref/DAVIDKnowledgebase/ENSEMBL_GENE_ID2GOTERM_BP_ALL.txt", sep = "\t", header = F)
  unique(str_count(DAVID_BP$V2,"~")) 
  DAVID_BP$Ontology <- "BP"
  DAVID_BP$GO_term <- unlist(lapply(DAVID_BP$V2, function(x){
    unlist(strsplit(x, "~"))[1]
  }))
  DAVID_BP$Description <- unlist(lapply(DAVID_BP$V2, function(x){
    unlist(strsplit(x, "~"))[2]
  }))
  DAVID_CC <- read.csv("./Ref/DAVIDKnowledgebase/ENSEMBL_GENE_ID2GOTERM_CC_ALL.txt", sep = "\t", header = F)
  unique(str_count(DAVID_CC$V2,"~")) 
  DAVID_CC$Ontology <- "CC"
  DAVID_CC$GO_term <- unlist(lapply(DAVID_CC$V2, function(x){
    unlist(strsplit(x, "~"))[1]
  }))
  DAVID_CC$Description <- unlist(lapply(DAVID_CC$V2, function(x){
    unlist(strsplit(x, "~"))[2]
  }))
  DAVID <- as.data.frame(rbind(DAVID_BP[,c(1,3,4,5)],
                               DAVID_CC[,c(1,3,4,5)],
                               DAVID_MF[,c(1,3,4,5)]))
  colnames(DAVID)[1] <- "Gene_ID"
  DAVID <- DAVID[!duplicated(DAVID),]
  DAVID$gene_count <- as.numeric(table(DAVID$GO_term)[DAVID$GO_term])
  DAVID_passfilter <- DAVID[which(DAVID$gene_count >= 5),]
  
  ATAC_common_tissues
  markersGenes_celltype <- getMarkerFeatures(
    ArchRProj = ATAC_common_tissues, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerGenesList_celltype <- getMarkers(markersGenes_celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  markerGenesList_celltype <- markerGenesList_celltype@listData
  markerGenesList_celltype_TF <- markerGenesList_celltype
  rownames(gene_id_map) <- gene_id_map$symbol
  for (i in names(markerGenesList_celltype_TF)) {
    temp <- markerGenesList_celltype_TF[[i]]
    temp <- as.data.frame(temp)
    temp <- temp$name
    temp <- temp[which(temp %in% TF_Family$symbol)]
    markerGenesList_celltype_TF[[i]] <- gene_id_map[temp, "gene_id"]
  }
  
  global_only_tf <- setdiff(global_tf, gene_linked_tf)
  rownames(Annotation_genes) <- Annotation_genes$`CGSNL Gene Symbol`
  global_only_tf <- Annotation_genes[global_only_tf, "Locus_ID"]
  genes_list <- list(
    "Only global constitutive peaks enriched TFs" = global_only_tf
  )
  genes_list <- c(genes_list, markerGenesList_celltype_TF)
  
  for (i in names(genes_list)) {
    print(i)
    genes <- genes_list[[i]]
    GO_result <- c()
    for (z in unique(DAVID_passfilter$GO_term)) {
      temp <- DAVID_passfilter[which(DAVID_passfilter$GO_term==z),]
      gene_count <- intersect(genes,
                              temp$Gene_ID)
      if (length(gene_count) > 0) {
        gene_ratio <- length(gene_count) / nrow(temp)
        ontology <- unique(temp$Ontology)
        description <- unique(temp$Description)
        x <- length(gene_count) - 1
        m <- length(unique(temp$Gene_ID))
        n <- 37000 - m
        k <- length(genes)
        p <- phyper(x, m, n, k, lower.tail = F)
        GO_result <- rbind.data.frame(GO_result,
                                      data.frame(GO_term = z,
                                                 Ontology = ontology,
                                                 Description = description,
                                                 Gene_count = length(gene_count),
                                                 Gene_ratio = gene_ratio,
                                                 P_value = p))
      } else {
        next
      }
    } # GO term enrichment
    GO_result <- GO_result[which(GO_result$Ontology == "BP"),]
    if (is.null(GO_result)) {
      GO_result <- paste0("There is only ", length(genes), " TFs ",", as a result, no GO pathways were enriched.")   
      write.table(GO_result, quote = F, row.names = F, col.names = F, sep = "\t",
                  paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result.tsv"))
      GO_result_Sig <- paste0("There is no GO pathways were significantly enriched.") 
      write.table(GO_result_Sig, quote = F, row.names = F, col.names = F, sep = "\t",
                  paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result_Sig.tsv"))
    } else {
      GO_result$FDR <- p.adjust(GO_result$P_value, method = "fdr")
      GO_result <- GO_result[order(GO_result$Ontology, GO_result$P_value),]
      GO_result_Sig <- GO_result[which(GO_result$P_value<0.05),]
      if (nrow(GO_result_Sig) < 1) {
        write.table(GO_result, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result.tsv"))
        GO_result_Sig <- paste0("There is no GO pathways were significantly enriched.")
        write.table(GO_result_Sig, quote = F, row.names = F, col.names = F, sep = "\t",
                    paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result_Sig.tsv"))
      } else {
        write.table(GO_result, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result.tsv"))
        write.table(GO_result_Sig, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result_Sig.tsv"))
      }
    }
  }
  
  # barplot for Sig-enrichment
  Enrichment_data <- c()
  for (i in names(genes_list)) {
    print(i)
    temp_enrich <- read.table(paste0("./Figure6_5_",gsub(pattern = "/", replacement = "_", fixed = T, x = i),"_enrich_GO_result_Sig.tsv"), sep = "\t", header = T, row.names = 1)
    temp_enrich$Group <- unlist(strsplit(i, "_genes"))[1]
    Enrichment_data <- as.data.frame(rbind(Enrichment_data,
                                           temp_enrich[1:20,]))
  }
  length(unique(Enrichment_data$Description))
  
  
  library(Hmisc)
  Enrichment_data$Description2 <- paste(Enrichment_data$Description, Enrichment_data$Group, sep = "_")
  Enrichment_data$Description2 <- capitalize(as.character(Enrichment_data$Description2))
  Enrichment_data$Description <- capitalize(as.character(Enrichment_data$Description))
  Enrichment_data$Description2 <- factor(Enrichment_data$Description2,
                                         levels = Enrichment_data$Description2,
                                         labels = Enrichment_data$Description2)
  
  Enrichment_data$FDR_2 <- -log10(Enrichment_data$FDR)
  Enrichment_data$P_value_2 <- -log10(Enrichment_data$P_value)
  
  pdf("Figure6_5_Global constitutive peaks enriched TFs_Top20.pdf", width = 10, height = 4)
  ggplot(Enrichment_data, aes(x = fct_reorder(Description2, P_value_2, .desc = F), y = P_value_2)) +
    geom_bar(aes(fill = P_value_2,
                 color = P_value_2),
             stat = "identity",
             # fill = "#87CEEB",
             # color = "black",
             width = 0.7) +
    labs(x = "", y = "-Log10(P value)") +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,3.3)) +
    coord_flip() +
    theme_classic() +
    facet_grid(Group ~ ., scale = "free") +
    scale_color_viridis(option = "C") +
    scale_fill_viridis(option = "C") +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))
  dev.off()
  
  
}

##### Lineage TF
{
  nrow(peaks_group_PD)
  table(peaks_group_PD$peakType)
  ### lineages
  Lineage_Epidermis_peaks_PD <- Lineage_Epidermis_peaks[which(Lineage_Epidermis_peaks %in% peaks_group_PD$ID)]
  Lineage_Fiber_peaks_PD <- Lineage_Fiber_peaks[which(Lineage_Fiber_peaks %in% peaks_group_PD$ID)]
  Lineage_Meristem_peaks_PD <- Lineage_Meristem_peaks[which(Lineage_Meristem_peaks %in% peaks_group_PD$ID)]
  Lineage_Rachis_peaks_PD <- Lineage_Rachis_peaks[which(Lineage_Rachis_peaks %in% peaks_group_PD$ID)]
  Lineage_RootGroundTissue_peaks_PD <- Lineage_RootGroundTissue_peaks[which(Lineage_RootGroundTissue_peaks %in% peaks_group_PD$ID)]
  Lineage_Vasculature_peaks_PD <- Lineage_Vasculature_peaks[which(Lineage_Vasculature_peaks %in% peaks_group_PD$ID)]
  Lineage_Mesophyll_peaks_PD <- Lineage_Mesophyll_peaks[which(Lineage_Mesophyll_peaks %in% peaks_group_PD$ID)]
  
  matches <- getMatches(ATAC_common_tissues, name = "TF-Motif")
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  rownames(matches) <- pr1
  matches <- matches[which(rownames(matches) %in% c(Lineage_Epidermis_peaks_PD,
                                                    Lineage_Fiber_peaks_PD,
                                                    Lineage_Meristem_peaks_PD,
                                                    Lineage_Rachis_peaks_PD,
                                                    Lineage_RootGroundTissue_peaks_PD,
                                                    Lineage_Vasculature_peaks_PD,
                                                    Lineage_Mesophyll_peaks_PD)),]
  
  temp <- matrix(FALSE, nrow = nrow(matches), ncol = 7)
  rownames(temp) <- rownames(matches)
  colnames(temp) <- names(Lineage_color)
  temp[Lineage_Epidermis_peaks_PD, 1] <- TRUE
  temp[Lineage_Fiber_peaks_PD, 2] <- TRUE
  temp[Lineage_Meristem_peaks_PD, 3] <- TRUE
  temp[Lineage_Rachis_peaks_PD, 4] <- TRUE
  temp[Lineage_RootGroundTissue_peaks_PD, 5] <- TRUE
  temp[Lineage_Vasculature_peaks_PD, 6] <- TRUE
  temp[Lineage_Mesophyll_peaks_PD, 7] <- TRUE
  identical(rownames(temp), rownames(matches))
  
  colnames(temp)
  Lineage_Epidermis_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,1]), 1:nrow(matches))
  Lineage_Epidermis_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_Epidermis_peaks_PD_Enrich$Enrichment)
  Lineage_Epidermis_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_Epidermis_peaks_PD_Enrich$Enrichmented[which(Lineage_Epidermis_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_Epidermis_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Epidermis_peaks_PD_Enrich$Enrichmented)
  
  colnames(temp)
  Lineage_Fiber_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,2]), 1:nrow(matches))
  Lineage_Fiber_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_Fiber_peaks_PD_Enrich$Enrichment)
  Lineage_Fiber_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_Fiber_peaks_PD_Enrich$Enrichmented[which(Lineage_Fiber_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_Fiber_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Fiber_peaks_PD_Enrich$Enrichmented)
  
  colnames(temp)
  Lineage_Meristem_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,3]), 1:nrow(matches))
  Lineage_Meristem_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_Meristem_peaks_PD_Enrich$Enrichment)
  Lineage_Meristem_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_Meristem_peaks_PD_Enrich$Enrichmented[which(Lineage_Meristem_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_Meristem_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Meristem_peaks_PD_Enrich$Enrichmented)
  
  colnames(temp)
  Lineage_Rachis_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,4]), 1:nrow(matches))
  Lineage_Rachis_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_Rachis_peaks_PD_Enrich$Enrichment)
  Lineage_Rachis_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_Rachis_peaks_PD_Enrich$Enrichmented[which(Lineage_Rachis_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_Rachis_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Rachis_peaks_PD_Enrich$Enrichmented)
  
  colnames(temp)
  Lineage_RootGroundTissue_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,5]), 1:nrow(matches))
  Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichment)
  Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichmented[which(Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_RootGroundTissue_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichmented)
  
  colnames(temp)
  Lineage_Vasculature_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,6]), 1:nrow(matches))
  Lineage_Vasculature_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_Vasculature_peaks_PD_Enrich$Enrichment)
  Lineage_Vasculature_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_Vasculature_peaks_PD_Enrich$Enrichmented[which(Lineage_Vasculature_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_Vasculature_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Vasculature_peaks_PD_Enrich$Enrichmented)
  
  colnames(temp)
  Lineage_Mesophyll_peaks_PD_Enrich <- .computeEnrichment(matches, which(temp[,7]), 1:nrow(matches))
  Lineage_Mesophyll_peaks_PD_Enrich$Enrichment_log <- log2(Lineage_Mesophyll_peaks_PD_Enrich$Enrichment)
  Lineage_Mesophyll_peaks_PD_Enrich$Enrichmented <- "NO"
  Lineage_Mesophyll_peaks_PD_Enrich$Enrichmented[which(Lineage_Mesophyll_peaks_PD_Enrich$Enrichment_log >= 0.25 & Lineage_Mesophyll_peaks_PD_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Lineage_Mesophyll_peaks_PD_Enrich$Enrichmented)
  
  venn_list_3 <- list(Epidermis = Lineage_Epidermis_peaks_PD_Enrich[which(Lineage_Epidermis_peaks_PD_Enrich$Enrichmented == "YES"),"feature"],
                    # Fiber = Lineage_Fiber_peaks_PD_Enrich[which(Lineage_Fiber_peaks_PD_Enrich$Enrichmented == "YES"),"feature"],
                    Meristem = Lineage_Meristem_peaks_PD_Enrich[which(Lineage_Meristem_peaks_PD_Enrich$Enrichmented == "YES"),"feature"],
                    Rachis = Lineage_Rachis_peaks_PD_Enrich[which(Lineage_Rachis_peaks_PD_Enrich$Enrichmented == "YES"),"feature"],
                    RootGroundTissue = Lineage_RootGroundTissue_peaks_PD_Enrich[which(Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichmented == "YES"),"feature"],
                    Vasculature = Lineage_Vasculature_peaks_PD_Enrich[which(Lineage_Vasculature_peaks_PD_Enrich$Enrichmented == "YES"),"feature"],
                    Mesophyll = Lineage_Mesophyll_peaks_PD_Enrich[which(Lineage_Mesophyll_peaks_PD_Enrich$Enrichmented == "YES"),"feature"])
  pdf("Figure6_5_Lineage_TF_venn.pdf", width = 7, height = 7)
  venn(venn_list_3,
       zcolor = c("#7570B3",
                  # "#E7298A",
                  "#1B9E77",
                  "#A6761D", "#D95F02",
                  "#E6AB02", "#66A61E"), # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
       opacity = 0.6,  # 调整颜色透明度
       box = F,        # 是否添加边框
       ilcs = 1,     # 数字大小
       ellipse = T,
       sncs = 1.2        # 组名字体大小
  )
  dev.off()
  
  inter <- get.venn.partitions(venn_list_3)
  for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
  inter <- subset(inter, select = -..values.. )
  inter <- subset(inter, select = -..set.. )
  write.table(inter, "Figure6_5_Lineage_TF_venn_result.csv", row.names = FALSE, sep = ',', quote = FALSE)
  
  Lineage_TFs <- unique(c(venn_list_3$Epidermis,
                          venn_list_3$Fiber,
                          venn_list_3$Meristem,
                          venn_list_3$Rachis,
                          venn_list_3$RootGroundTissue,
                          venn_list_3$Vasculature,
                          venn_list_3$Mesophyll))
  
  Lineage_Epidermis_peaks_PD_Enrich$Lineage <- "Epidermis"
  Lineage_Fiber_peaks_PD_Enrich$Lineage <- "Fiber"
  Lineage_Meristem_peaks_PD_Enrich$Lineage <- "Meristem"
  Lineage_Rachis_peaks_PD_Enrich$Lineage <- "Rachis"
  Lineage_RootGroundTissue_peaks_PD_Enrich$Lineage <- "RootGroundTissue"
  Lineage_Vasculature_peaks_PD_Enrich$Lineage <- "Vasculature"
  Lineage_Mesophyll_peaks_PD_Enrich$Lineage <- "Mesophyll"
  
  all_enrich <- list_rbind(list(Lineage_Epidermis_peaks_PD_Enrich[which(Lineage_Epidermis_peaks_PD_Enrich$Enrichmented == "YES"),],
                                Lineage_Fiber_peaks_PD_Enrich[which(Lineage_Fiber_peaks_PD_Enrich$Enrichmented == "YES"),],
                                Lineage_Meristem_peaks_PD_Enrich[which(Lineage_Meristem_peaks_PD_Enrich$Enrichmented == "YES"),],
                                Lineage_Rachis_peaks_PD_Enrich[which(Lineage_Rachis_peaks_PD_Enrich$Enrichmented == "YES"),],
                                Lineage_RootGroundTissue_peaks_PD_Enrich[which(Lineage_RootGroundTissue_peaks_PD_Enrich$Enrichmented == "YES"),],
                                Lineage_Vasculature_peaks_PD_Enrich[which(Lineage_Vasculature_peaks_PD_Enrich$Enrichmented == "YES"),],
                                Lineage_Mesophyll_peaks_PD_Enrich[which(Lineage_Mesophyll_peaks_PD_Enrich$Enrichmented == "YES"),]))
  all_enrich_filtered <- matrix(0, nrow = length(Lineage_TFs), ncol = 7)
  rownames(all_enrich_filtered) <- Lineage_TFs
  colnames(all_enrich_filtered) <- names(Lineage_color)
  nrow(all_enrich) == sum(all_enrich$Lineage %in% names(Lineage_color))
  for (i in 1:nrow(all_enrich)) {
    all_enrich_filtered[all_enrich$feature[i], all_enrich$Lineage[i]] <- all_enrich$Enrichment_log[i]
  }
  max(all_enrich_filtered)
  
  
  min(all_enrich$Enrichment_log)
  col_fun <- circlize::colorRamp2(
    breaks = c(seq(from = 0.25, to = 0.9, length.out = 25)), 
    colors = c(viridis::viridis(n = 25, option = "H"))
  )
  row_anno <- rowAnnotation(df = data.frame(Lineage = names(Lineage_color)),
                            col = list(Lineage = Lineage_color))
  pdf("Figure6_5_Lineage_enriched_TF_heatmap.pdf", width = 13, height = 3)
  Heatmap(t(all_enrich_filtered), cluster_rows = T, cluster_columns = T,
          col = col_fun, na_col = "#E0E0E0", row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param = list(title = "Log2(Enrichment)",
                                      at = c(0.25, 0.9)),
          column_names_gp = gpar(fontsize = 10), right_annotation = row_anno
  )
  dev.off()
  
  
  Epidermis_peaks_Enrich_ggplot <- gplot1(Lineage_Epidermis_peaks_PD_Enrich)
  Epidermis_peaks_Enrich_ggplot
  Fiber_peaks_Enrich_ggplot <- gplot1(Lineage_Fiber_peaks_PD_Enrich)
  Fiber_peaks_Enrich_ggplot
  Meristem_peaks_Enrich_ggplot <- gplot1(Lineage_Meristem_peaks_PD_Enrich)
  Meristem_peaks_Enrich_ggplot
  Rachis_peaks_Enrich_ggplot <- gplot1(Lineage_Rachis_peaks_PD_Enrich)
  Rachis_peaks_Enrich_ggplot
  RootGroundTissue_peaks_Enrich_ggplot <- gplot1(Lineage_RootGroundTissue_peaks_PD_Enrich)
  RootGroundTissue_peaks_Enrich_ggplot
  Vasculature_peaks_Enrich_ggplot <- gplot1(Lineage_Vasculature_peaks_PD_Enrich)
  Vasculature_peaks_Enrich_ggplot
  Mesophyll_peaks_Enrich_ggplot <- gplot1(Lineage_Mesophyll_peaks_PD_Enrich)
  Mesophyll_peaks_Enrich_ggplot
  
  pdf("Figure6_4_Lineage_peak_TF_enrich.pdf", width = 14, height = 7)
  Epidermis_peaks_Enrich_ggplot +
    Fiber_peaks_Enrich_ggplot +
    Meristem_peaks_Enrich_ggplot +
    Rachis_peaks_Enrich_ggplot +
    RootGroundTissue_peaks_Enrich_ggplot +
    Vasculature_peaks_Enrich_ggplot +
    Mesophyll_peaks_Enrich_ggplot +
    plot_layout(nrow = 2, byrow = FALSE)
  dev.off()
  
  
  
}

rownames(RNA)
RNA_mean

RNA_mean_all
median(RNA_mean_all)

temp <- setdiff(global_tf, gene_linked_tf)
temp <- Annotation_genes[temp, "DataSets_Symbol"]
temp <- temp[temp %in% names(RNA_mean_all)]
length(temp)
sum(RNA_mean_all[temp] > quantile(RNA_mean_all)[2]) # 65个表达水平大于25%
sum(RNA_mean_all[temp] > quantile(RNA_mean_all)[3]) # 45个表达水平大于50%
sum(RNA_mean_all[temp] > quantile(RNA_mean_all)[4]) # 27个表达水平大于27%



######### 补充

# OsRLCK85_peaks <- peaksets_2[which(peaksets_2$nearestGene == "OsRLCK85"),]
# OsRLCK85_peaks <- OsRLCK85_peaks[order(OsRLCK85_peaks$start, decreasing = F),]
# OsRLCK85_peaks$order <- 1:nrow(OsRLCK85_peaks)
# pdf("Figure6_补充_OsRLCK85_peaks_pearson_barplot.pdf", width = 4, height = 2)
# ggplot(data = OsRLCK85_peaks, aes(x = order, y = pearson, fill = peakGroup)) +
#   geom_bar(stat = "identity", width = 0.5) +
#   scale_fill_manual(values = peak_type_color) +
#   theme_bw() +
#   geom_hline(yintercept = 0.777, colour = "red", size = 0.15) +
#   theme(axis.text = element_text(size = 10, color = "black"),
#         legend.position = "none")
# dev.off()
# 
# pdf("Figure6_补充_OsRLCK85_peaks_pearson_point_plot.pdf", width = 4, height = 2)
# ggplot(data = OsRLCK85_peaks, aes(x = order, y = pearson, color = peakGroup)) +
#   geom_point(shape = 15, size = 3) +
#   scale_color_manual(values = peak_type_color) +
#   theme_bw() +
#   geom_hline(yintercept = 0.777, colour = "red", size = 0.15) +
#   theme(axis.text = element_text(size = 10, color = "black"),
#         legend.position = "none")
# dev.off()


{
  table(RNA_Celltye_DEG$cluster)
  sum(Gene_Linked_Peaks_DEGs %in% peaksets_2$ID)
  sum(peaksets_2$nearestGene %in% rownames(RNA_common_tissues)) == nrow(peaksets_2)
  Gene_Linked_Peaks_DF <- peaksets_2[Gene_Linked_Peaks_DEGs,]
  Gene_Linked_Peaks_DF <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$peakType %in% c("Promoter", "Distal")),]
  
  GL_MO <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                        RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Mesophyll (MO)")]),]
  GL_Fiber <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                           RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Fiber")]),]
  GL_Phloem <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                            RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Phloem")]),]
  GL_MP <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                        RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Mesophyll precursor")]),]
  GL_Procambium <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                                RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Procambium")]),]
  GL_Epidermis <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                               RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Epidermis")]),]
  GL_Cortex <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                            RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Cortex")]),]
  GL_Endodermis <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                                RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Endodermis")]),]
  GL_IM <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                        RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Inflorescence meristem (IM)")]),]
  GL_SM <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                        RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Spikelet meristem (SM)")]),]
  GL_BM <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                        RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Branch meristems (BM)")]),]
  GL_Rachis <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                            RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Rachis")]),]
  GL_CBB <- Gene_Linked_Peaks_DF[which(Gene_Linked_Peaks_DF$nearestGene %in%
                                         RNA_Celltye_DEG$gene[which(RNA_Celltye_DEG$cluster == "Cryptic bract/bract (cb/b)")]),]
  
  matches <- getMatches(ATAC_common_tissues, name = "TF-Motif")
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  rownames(matches) <- pr1
  matches <- matches[which(rownames(matches) %in% c(GL_MO$ID,
                                                    GL_Fiber$ID,
                                                    GL_Phloem$ID,
                                                    GL_MP$ID,
                                                    GL_Procambium$ID,
                                                    GL_Epidermis$ID,
                                                    GL_Cortex$ID,
                                                    GL_Endodermis$ID,
                                                    GL_IM$ID,
                                                    GL_SM$ID,
                                                    GL_BM$ID,
                                                    GL_Rachis$ID,
                                                    GL_CBB$ID
  )),]
  
  temp <- matrix(FALSE, nrow = nrow(matches), ncol = 13)
  rownames(temp) <- rownames(matches)
  colnames(temp) <- as.character(unique(RNA_Celltye_DEG$cluster))
  temp[GL_MO$ID, 1] <- TRUE
  temp[GL_Fiber$ID, 2] <- TRUE
  temp[GL_Phloem$ID, 3] <- TRUE
  temp[GL_MP$ID, 4] <- TRUE
  temp[GL_Procambium$ID, 5] <- TRUE
  temp[GL_Epidermis$ID, 6] <- TRUE
  temp[GL_Cortex$ID, 7] <- TRUE
  temp[GL_Endodermis$ID, 8] <- TRUE
  temp[GL_IM$ID, 9] <- TRUE
  temp[GL_SM$ID, 10] <- TRUE
  temp[GL_BM$ID, 11] <- TRUE
  temp[GL_Rachis$ID, 12] <- TRUE
  temp[GL_CBB$ID, 13] <- TRUE
  
  identical(rownames(temp), rownames(matches))
  
  colnames(temp)
  MO_Enrich <- .computeEnrichment(matches, which(temp[,1]), 1:nrow(matches))
  MO_Enrich$Enrichment_log <- log2(MO_Enrich$Enrichment)
  MO_Enrich$Enrichmented <- "NO"
  MO_Enrich$Enrichmented[which(MO_Enrich$Enrichment_log >= 0.25 & MO_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(MO_Enrich$Enrichmented)
  
  colnames(temp)
  Fiber_Enrich <- .computeEnrichment(matches, which(temp[,2]), 1:nrow(matches))
  Fiber_Enrich$Enrichment_log <- log2(Fiber_Enrich$Enrichment)
  Fiber_Enrich$Enrichmented <- "NO"
  Fiber_Enrich$Enrichmented[which(Fiber_Enrich$Enrichment_log >= 0.25 & Fiber_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Fiber_Enrich$Enrichmented)
  
  colnames(temp)
  Phloem_Enrich <- .computeEnrichment(matches, which(temp[,3]), 1:nrow(matches))
  Phloem_Enrich$Enrichment_log <- log2(Phloem_Enrich$Enrichment)
  Phloem_Enrich$Enrichmented <- "NO"
  Phloem_Enrich$Enrichmented[which(Phloem_Enrich$Enrichment_log >= 0.25 & Phloem_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Phloem_Enrich$Enrichmented)
  
  colnames(temp)
  MP_Enrich <- .computeEnrichment(matches, which(temp[,4]), 1:nrow(matches))
  MP_Enrich$Enrichment_log <- log2(MP_Enrich$Enrichment)
  MP_Enrich$Enrichmented <- "NO"
  MP_Enrich$Enrichmented[which(MP_Enrich$Enrichment_log >= 0.25 & MP_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(MP_Enrich$Enrichmented)
  
  colnames(temp)
  Procambium_Enrich <- .computeEnrichment(matches, which(temp[,5]), 1:nrow(matches))
  Procambium_Enrich$Enrichment_log <- log2(Procambium_Enrich$Enrichment)
  Procambium_Enrich$Enrichmented <- "NO"
  Procambium_Enrich$Enrichmented[which(Procambium_Enrich$Enrichment_log >= 0.25 & Procambium_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Procambium_Enrich$Enrichmented)
  
  colnames(temp)
  Epidermis_Enrich <- .computeEnrichment(matches, which(temp[,6]), 1:nrow(matches))
  Epidermis_Enrich$Enrichment_log <- log2(Epidermis_Enrich$Enrichment)
  Epidermis_Enrich$Enrichmented <- "NO"
  Epidermis_Enrich$Enrichmented[which(Epidermis_Enrich$Enrichment_log >= 0.25 & Epidermis_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Epidermis_Enrich$Enrichmented)
  
  colnames(temp)
  Cortex_Enrich <- .computeEnrichment(matches, which(temp[,7]), 1:nrow(matches))
  Cortex_Enrich$Enrichment_log <- log2(Cortex_Enrich$Enrichment)
  Cortex_Enrich$Enrichmented <- "NO"
  Cortex_Enrich$Enrichmented[which(Cortex_Enrich$Enrichment_log >= 0.25 & Cortex_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Cortex_Enrich$Enrichmented)
  
  colnames(temp)
  Endodermis_Enrich <- .computeEnrichment(matches, which(temp[,8]), 1:nrow(matches))
  Endodermis_Enrich$Enrichment_log <- log2(Endodermis_Enrich$Enrichment)
  Endodermis_Enrich$Enrichmented <- "NO"
  Endodermis_Enrich$Enrichmented[which(Endodermis_Enrich$Enrichment_log >= 0.25 & Endodermis_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Endodermis_Enrich$Enrichmented)
  
  colnames(temp)
  IM_Enrich <- .computeEnrichment(matches, which(temp[,9]), 1:nrow(matches))
  IM_Enrich$Enrichment_log <- log2(IM_Enrich$Enrichment)
  IM_Enrich$Enrichmented <- "NO"
  IM_Enrich$Enrichmented[which(IM_Enrich$Enrichment_log >= 0.25 & IM_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(IM_Enrich$Enrichmented)
  
  
  colnames(temp)
  SM_Enrich <- .computeEnrichment(matches, which(temp[,10]), 1:nrow(matches))
  SM_Enrich$Enrichment_log <- log2(SM_Enrich$Enrichment)
  SM_Enrich$Enrichmented <- "NO"
  SM_Enrich$Enrichmented[which(SM_Enrich$Enrichment_log >= 0.25 & SM_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(SM_Enrich$Enrichmented)
  
  colnames(temp)
  BM_Enrich <- .computeEnrichment(matches, which(temp[,11]), 1:nrow(matches))
  BM_Enrich$Enrichment_log <- log2(BM_Enrich$Enrichment)
  BM_Enrich$Enrichmented <- "NO"
  BM_Enrich$Enrichmented[which(BM_Enrich$Enrichment_log >= 0.25 & BM_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(BM_Enrich$Enrichmented)
  
  colnames(temp)
  Rachis_Enrich <- .computeEnrichment(matches, which(temp[,12]), 1:nrow(matches))
  Rachis_Enrich$Enrichment_log <- log2(Rachis_Enrich$Enrichment)
  Rachis_Enrich$Enrichmented <- "NO"
  Rachis_Enrich$Enrichmented[which(Rachis_Enrich$Enrichment_log >= 0.25 & Rachis_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(Rachis_Enrich$Enrichmented)
  
  colnames(temp)
  CBB_Enrich <- .computeEnrichment(matches, which(temp[,13]), 1:nrow(matches))
  CBB_Enrich$Enrichment_log <- log2(CBB_Enrich$Enrichment)
  CBB_Enrich$Enrichmented <- "NO"
  CBB_Enrich$Enrichmented[which(CBB_Enrich$Enrichment_log >= 0.25 & CBB_Enrich$mlog10p >= -log10(0.05))] <- "YES"
  table(CBB_Enrich$Enrichmented)
  
  
  venn_list_3 <- list(MO = MO_Enrich[which(MO_Enrich$Enrichmented == "YES"),"feature"],
                      Fiber = Fiber_Enrich[which(Fiber_Enrich$Enrichmented == "YES"),"feature"],
                      Phloem = Phloem_Enrich[which(Phloem_Enrich$Enrichmented == "YES"),"feature"],
                      MP = MP_Enrich[which(MP_Enrich$Enrichmented == "YES"),"feature"],
                      Procambium = Procambium_Enrich[which(Procambium_Enrich$Enrichmented == "YES"),"feature"],
                      Epidermis = Epidermis_Enrich[which(Epidermis_Enrich$Enrichmented == "YES"),"feature"],
                      Cortex = Cortex_Enrich[which(Cortex_Enrich$Enrichmented == "YES"),"feature"],
                      Endodermis = Endodermis_Enrich[which(Endodermis_Enrich$Enrichmented == "YES"),"feature"],
                      IM = IM_Enrich[which(IM_Enrich$Enrichmented == "YES"),"feature"],
                      SM = SM_Enrich[which(SM_Enrich$Enrichmented == "YES"),"feature"],
                      BM = BM_Enrich[which(BM_Enrich$Enrichmented == "YES"),"feature"],
                      Rachis = Rachis_Enrich[which(Rachis_Enrich$Enrichmented == "YES"),"feature"],
                      CBB = CBB_Enrich[which(CBB_Enrich$Enrichmented == "YES"),"feature"]
  )
  
  Celltype_TFs <- unique(c(venn_list_3$MO,
                           venn_list_3$Fiber,
                           venn_list_3$Phloem,
                           venn_list_3$MP,
                           venn_list_3$Procambium,
                           venn_list_3$Epidermis,
                           venn_list_3$Cortex,
                           venn_list_3$Endodermis,
                           venn_list_3$IM,
                           venn_list_3$SM,
                           venn_list_3$BM,
                           venn_list_3$Rachis,
                           venn_list_3$CBB))
  
  colnames(temp)
  
  MO_Enrich$Celltype <- "Mesophyll (MO)"
  Fiber_Enrich$Celltype <- "Fiber"
  Phloem_Enrich$Celltype <- "Phloem"
  MP_Enrich$Celltype <- "Mesophyll precursor"
  Procambium_Enrich$Celltype <- "Procambium"
  Epidermis_Enrich$Celltype <- "Epidermis"
  Cortex_Enrich$Celltype <- "Cortex"
  Endodermis_Enrich$Celltype <- "Endodermis"
  IM_Enrich$Celltype <- "Inflorescence meristem (IM)"
  SM_Enrich$Celltype <- "Spikelet meristem (SM)"
  BM_Enrich$Celltype <- "Branch meristems (BM)"
  Rachis_Enrich$Celltype <- "Rachis"
  CBB_Enrich$Celltype <- "Cryptic bract/bract (cb/b)"
  
  
  all_enrich <- list_rbind(list(MO_Enrich[which(MO_Enrich$Enrichmented == "YES"),],
                                Fiber_Enrich[which(Fiber_Enrich$Enrichmented == "YES"),],
                                Phloem_Enrich[which(Phloem_Enrich$Enrichmented == "YES"),],
                                MP_Enrich[which(MP_Enrich$Enrichmented == "YES"),],
                                Procambium_Enrich[which(Procambium_Enrich$Enrichmented == "YES"),],
                                Epidermis_Enrich[which(Epidermis_Enrich$Enrichmented == "YES"),],
                                Cortex_Enrich[which(Cortex_Enrich$Enrichmented == "YES"),],
                                Endodermis_Enrich[which(Endodermis_Enrich$Enrichmented == "YES"),],
                                IM_Enrich[which(IM_Enrich$Enrichmented == "YES"),],
                                SM_Enrich[which(SM_Enrich$Enrichmented == "YES"),],
                                BM_Enrich[which(BM_Enrich$Enrichmented == "YES"),],
                                Rachis_Enrich[which(Rachis_Enrich$Enrichmented == "YES"),],
                                CBB_Enrich[which(CBB_Enrich$Enrichmented == "YES"),]
  ))
  all_enrich_filtered <- matrix(0, nrow = length(Celltype_TFs), ncol = 13)
  rownames(all_enrich_filtered) <- Celltype_TFs
  colnames(all_enrich_filtered) <- colnames(temp)
  nrow(all_enrich) == sum(all_enrich$Celltype %in% colnames(temp))
  for (i in 1:nrow(all_enrich)) {
    all_enrich_filtered[all_enrich$feature[i], all_enrich$Celltype[i]] <- all_enrich$Enrichment_log[i]
  }
  max(all_enrich_filtered)
  
  
  min(all_enrich$Enrichment_log)
  col_fun <- circlize::colorRamp2(
    breaks = c(seq(from = 0.25, to = 2, length.out = 25)), 
    colors = c(viridis::viridis(n = 25, option = "H"))
  )
  row_anno <- rowAnnotation(df = data.frame(Celltype = colnames(temp)),
                            col = list(Celltype = celltype_color[colnames(temp)]))
  pdf("Figure6_补充_celltype_enriched_TF_heatmap.pdf", width = 12, height = 3.5)
  Heatmap(t(all_enrich_filtered), cluster_rows = T, cluster_columns = T,
          col = col_fun, na_col = "#E0E0E0", row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param = list(title = "Log2(Enrichment)",
                                      at = c(0.25, 2)),
          column_names_gp = gpar(fontsize = 10), right_annotation = row_anno
  )
  dev.off()
  
  RNA_celltype_regulon <- readRDS("Figure3_RNA_celltype_regulon.rds")
  rownames(Annotation_genes) <- Annotation_genes$Locus_ID
  RNA_celltype_regulon$TF_name <- Annotation_genes[RNA_celltype_regulon$TF, "CGSNL Gene Symbol"]
  RNA_celltype_regulon <- as.data.frame(na.omit(RNA_celltype_regulon))
  
  rownames(Annotation_genes) <- Annotation_genes$`CGSNL Gene Symbol`
  
  nrow(all_enrich)
  all_enrich$ID <- Annotation_genes[all_enrich$feature, "Locus_ID"]
  sum(all_enrich$feature %in% RNA_celltype_regulon$TF_name)
  sum(all_enrich$ID %in% RNA_celltype_regulon$TF)
  
  Celltype_TFs_intersect_RNA <- data.frame(celltype = colnames(temp))
  RNA_TF_num <- c()
  ATAC_TF_num <- c()
  common_num <- c()
  for (i in Celltype_TFs_intersect_RNA$celltype) {
    ATAC_TF_num_temp <- all_enrich[which(all_enrich$Celltype == i),]
    ATAC_TF_num_temp <- ATAC_TF_num_temp$ID
    ATAC_TF_num <- c(ATAC_TF_num, length(ATAC_TF_num_temp))
    RNA_TF_num_temp <- RNA_celltype_regulon[which(RNA_celltype_regulon$Celltype == i),]
    RNA_TF_num <- c(RNA_TF_num, nrow(RNA_TF_num_temp))
    common_num_temp <- intersect(ATAC_TF_num_temp,
                                 RNA_TF_num_temp$TF)
    common_num <- c(common_num, length(common_num_temp))
  }
  
  Celltype_TFs_intersect_RNA$RNA_TF_num <- RNA_TF_num
  Celltype_TFs_intersect_RNA$ATAC_TF_num <- ATAC_TF_num
  Celltype_TFs_intersect_RNA$common_num <- common_num
}

{
  RNA_Celltye_DEG_DF <- matrix(0, nrow = length(unique(RNA_Celltye_DEG$gene)), ncol = 13)
  rownames(RNA_Celltye_DEG_DF) <- unique(RNA_Celltye_DEG$gene)
  colnames(RNA_Celltye_DEG_DF) <- unique(RNA_Celltye_DEG$cluster)
  for (i in 1:nrow(RNA_Celltye_DEG)) {
    
  }
  
  markerPeaksDF
  
  peak_group_selected <- peak_group[which(peak_group$nearestGene == "OsRLCK85"),]
  Gene_Linked_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Gene-linked peaks"),"peakID"])]
  Global_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Global constitutive peaks"),"peakID"])]
  Lineage_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Lineage constitutive peaks"),"peakID"])]
  Other_Peaks_PD_RG <- peakSet[which(peakSet$Peak_id %in% peak_group_selected[which(peak_group_selected$peakType == "Other peaks"),"peakID"])]
  
  featureList <- SimpleList("Gene-linked peaks" = Gene_Linked_Peaks_PD_RG,
                            "Global constitutive peaks" = Global_Peaks_PD_RG,
                            "Lineage constitutive peaks" = Lineage_Peaks_PD_RG,
                            "Other peaks" = Other_Peaks_PD_RG)
  peak_type_color
  
  "OsCCD1.1"
  "Oshox9"
  "OsARF1"
  "OsC3H53"
  
  ATAC_common_tissues$Celltype <- factor(ATAC_common_tissues$Celltype,
                                         levels = main_celltype$celltype)
  p <- plotBrowserTrack(ArchRProj = ATAC_common_tissues,
                        groupBy = "Celltype",
                        geneSymbol = c("OsRLCK85"),
                        upstream = 10500, downstream = 12000,
                        features = featureList,
                        peakPal = peak_type_color,
                        groupPal = celltype_color,
                        facetbaseSize = 10)
  dev.off()
  # pdf("Figure6_4_OsRLCK85_genomeBrowse.pdf", width = 6, height = 5)
  grid::grid.newpage()
  grid::grid.draw(p$OsRLCK85)
  dev.off()
  
  
  temp <- data.frame(Maincluster = RNA_common_tissues$Mainclusters,
                     expression = RNA_common_tissues@assays$RNA@data["OsRLCK85",])
  
  pdf("Figure6_4_OsRLCK85_RNA.pdf", width = 5, height = 3.5)
  ggplot() +
    geom_boxplot(data = temp, aes(x = Maincluster, y = expression, fill = Maincluster)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = Maincluster_color) +
    labs(x = "", y = "RNA Expression")
  dev.off()
  
  
}

#### DEG peak percent 
{
  peak_group_degs <- unique(peak_group$nearestGene)
  lineage_percent <- c()
  global_percent <- c()
  pg_percent <- c()
  other_percent <- c()
  for (i in peak_group_degs) {
    print(i)
    peak_group_temp <- peak_group[which(peak_group$nearestGene == i),]
    lineage_percent <- c(lineage_percent, sum(peak_group_temp$peakType == "Lineage constitutive peaks") / nrow(peak_group_temp))
    global_percent <- c(global_percent, sum(peak_group_temp$peakType == "Global constitutive peaks") / nrow(peak_group_temp))
    pg_percent <- c(pg_percent, sum(peak_group_temp$peakType == "Gene-linked peaks") / nrow(peak_group_temp))
    other_percent <- c(other_percent, sum(peak_group_temp$peakType == "Other peaks") / nrow(peak_group_temp))
  }
  DEG_peak_percent <- data.frame(DEGs = peak_group_degs,
                                 lineage_percent = lineage_percent,
                                 global_percent = global_percent,
                                 pg_percent = pg_percent,
                                 other_percent = other_percent)
  
  
  
}


####################
#################### 弃用代码
####################
### OsRLCK85 peaks
{
  # 
  # "OsRLCK85" %in% peaksets_2$nearestGene
  # OsRLCK85_peaks <- peaksets_2[which(peaksets_2$nearestGene == "OsRLCK85"),]
  # OsRLCK85_peaks <- OsRLCK85_peaks[order(OsRLCK85_peaks$peakGroup, decreasing = F),]
  # celltype_order <- c("Branch meristems (BM)", "Cryptic bract/bract (cb/b)", "Inflorescence meristem (IM)",
  #                     "Spikelet meristem (SM)", "Rachis", "Mesophyll precursor", "Mesophyll (MO)", "Endodermis",
  #                     "Cortex", "Procambium", "Phloem", "Fiber", "Epidermis")
  # temp_peak_matrix <- Peak_mean[celltype_order, OsRLCK85_peaks$ID]
  # # pheatmap(temp_peak_matrix, scale = "column")
  # temp_RNA_matrix <- data.frame(OsRLCK85 = RNA_mean[celltype_order, "OsRLCK85"])
  # temp <- as.data.frame(cbind(temp_RNA_matrix,
  #                             temp_peak_matrix[rownames(temp_RNA_matrix),]))
  # cellanno_row <- rowAnnotation(df = data.frame(Maincluster = main_celltype[celltype_order, "main"]),
  #                                  col = list(Maincluster = Maincluster_color))
  # temp <- temp[,c(1:6,8,7,13:15, 9:12,16:18)]
  # peak_region_color <- c("#ADC698", "#A0C1D1",
  #                        "#B18ED7", "#D676A1")
  # names(peak_region_color) <- c("Exonic", "Intronic", "Distal", "Promoter")
  # col_color <- c("#008B8B", peak_type_color)
  # names(col_color)[1] <- "RNA"
  # cellanno_col <- columnAnnotation(PG_cor = peaksets_2[colnames(temp)[2:18], "pearson"],
  #                                  df = data.frame(
  #                                       # Type = c("RNA", peaksets_2[colnames(temp)[2:18], "peakGroup"])
  #                                       Type = c(peaksets_2[colnames(temp)[2:18], "peakGroup"]),
  #                                       Region = c(peaksets_2[colnames(temp)[2:18], "peakType"])
  #                                       ),
  #                                  col = list(Type = col_color,
  #                                             PG_cor = circlize::colorRamp2(c(-1, 0, 1),
  #                                                                           c("#1f78b4", "white", "#CD5C5C")),
  #                                             Region = peak_region_color
  #                                             ))
  # 
  # pdf("Figure6_4_OsRLCK85_peak_notScaled.pdf", height = 6, width = 10)
  # Heatmap(temp[,-1], cluster_rows = F, cluster_columns = F,
  #         col = viridis::viridis(n = 100, option = "D"),
  #         name = "Accessibility", width = unit(80, "mm"),
  #         height = unit(45, "mm"), bottom_annotation = cellanno_col, right_annotation = cellanno_row)
  # dev.off()
  # 
  # pdf("Figure6_4_OsRLCK85_RNA_notScaled.pdf", height = 6, width = 10)
  # temp2 <- as.data.frame(temp[,1])
  # colnames(temp2) <- "OsRLCK85"
  # rownames(temp2) <- rownames(temp)
  # Heatmap(temp2, cluster_rows = F, cluster_columns = F,
  #         col = viridis::viridis(n = 100, option = "E"),
  #         name = "Accessibility", width = unit(80, "mm"),
  #         height = unit(45, "mm"), right_annotation = cellanno_row)
  # dev.off()
}

{
  # "OsRLCK85" %in% peaksets_2$nearestGene
  # OsRLCK85_peaks <- peaksets_2[which(peaksets_2$nearestGene == "OsRLCK85"),]
  # OsRLCK85_peaks <- OsRLCK85_peaks[order(OsRLCK85_peaks$peakGroup, decreasing = F),]
  # celltype_order <- c("Branch meristems (BM)", "Cryptic bract/bract (cb/b)", "Inflorescence meristem (IM)",
  #                     "Spikelet meristem (SM)", "Rachis", "Mesophyll precursor", "Mesophyll (MO)", "Endodermis",
  #                     "Cortex", "Procambium", "Phloem", "Fiber", "Epidermis")
  # temp_peak_matrix <- Peak_mean[celltype_order, OsRLCK85_peaks$ID]
  # temp_RNA_matrix <- data.frame(OsRLCK85 = RNA_mean[celltype_order, "OsRLCK85"])
  # temp <- as.data.frame(cbind(temp_RNA_matrix,
  #                             temp_peak_matrix[rownames(temp_RNA_matrix),]))
  # temp <- scale(temp)
  # max(temp)
  # min(temp)
  # temp[which(temp > 2)] <- 2
  # temp[which(temp < -2)] <- -2
  # cellanno_row <- rowAnnotation(df = data.frame(Maincluster = main_celltype[celltype_order, "main"]),
  #                               col = list(Maincluster = Maincluster_color))
  # temp <- temp[,c(1:6,8,7,13:15, 9:12,16:18)]
  # peak_region_color <- c("#ADC698", "#A0C1D1",
  #                        "#B18ED7", "#D676A1")
  # names(peak_region_color) <- c("Exonic", "Intronic", "Distal", "Promoter")
  # col_color <- c("#008B8B", peak_type_color)
  # names(col_color)[1] <- "RNA"
  # cellanno_col <- columnAnnotation(PG_cor = peaksets_2[colnames(temp)[2:18], "pearson"],
  #                                  df = data.frame(
  #                                    # Type = c("RNA", peaksets_2[colnames(temp)[2:18], "peakGroup"]),
  #                                    Type = c(peaksets_2[colnames(temp)[2:18], "peakGroup"]),
  #                                    Region = c(peaksets_2[colnames(temp)[2:18], "peakType"])
  #                                  ),
  #                                  col = list(Type = col_color,
  #                                             PG_cor = circlize::colorRamp2(c(-1, 0, 1),
  #                                                                           c("#1f78b4", "white", "#CD5C5C")),
  #                                             Region = peak_region_color
  #                                  ))
  # 
  # pdf("Figure6_4_OsRLCK85_peak_Scaled.pdf", height = 6, width = 10)
  # Heatmap(temp[,-1], cluster_rows = F, cluster_columns = F,
  #         col = viridis::viridis(n = 100, option = "D"),
  #         name = "Accessibility", width = unit(80, "mm"),
  #         height = unit(45, "mm"), bottom_annotation = cellanno_col, right_annotation = cellanno_row)
  # dev.off()
  # 
  # pdf("Figure6_4_OsRLCK85_RNA_Scaled.pdf", height = 6, width = 10)
  # temp2 <- as.data.frame(temp[,1])
  # colnames(temp2) <- "OsRLCK85"
  # rownames(temp2) <- rownames(temp)
  # Heatmap(temp2, cluster_rows = F, cluster_columns = F,
  #         col = viridis::viridis(n = 100, option = "D"),
  #         name = "Accessibility", width = unit(80, "mm"),
  #         height = unit(45, "mm"), right_annotation = cellanno_row)
  # dev.off()
}

{
  #   OsC3H53_peaks <- peaksets_2[which(peaksets_2$nearestGene == "OsC3H53"),]
  #   temp_peak_matrix <- Peak_mean[, OsC3H53_peaks$ID]
  #   pheatmap(temp_peak_matrix, scale = "none")
  #   temp_RNA_matrix <- data.frame(OsC3H53 = RNA_mean[, "OsC3H53"])
  #   temp <- as.data.frame(cbind(temp_RNA_matrix,
  #                               temp_peak_matrix[rownames(temp_RNA_matrix),]))
  #   temp <- scale(temp)
  #   cellanno_row <- rowAnnotation(df = data.frame(Maincluster = c("Meristem", "RootGroundTissue", "Meristem",
  #                                                                 "RootGroundTissue", "Epidermis", "Fiber",
  #                                                                 "Meristem", "Mesophyll", "Mesophyll",
  #                                                                 "Vasculature", "Vasculature", "Rachis",
  #                                                                 "Meristem")),
  #                                 col = list(Maincluster = Maincluster_color))
  #   col_color <- c("#008B8B", peak_type_color)
  #   names(col_color)[1] <- "RNA"
  #   cellanno_col <- columnAnnotation(df = data.frame(
  #     # Type = c("RNA", peaksets_2[colnames(temp)[2:18], "peakGroup"])
  #     Type = c(peaksets_2[colnames(temp)[2:ncol(temp)], "peakGroup"])
  #   ),
  #   col = list(Type = col_color))
  #   pdf("Figure6_4_ATAC_GeneScoreMatrix_celltype_cor.pdf", height = 6, width = 9)
  #   Heatmap(temp[,-1],
  #           col = viridis::viridis(n = 100, option = "D"),
  #           name = "Pearson Coefficient\nGeneScoreMatrix", width = unit(80, "mm"),
  #           height = unit(80, "mm"), bottom_annotation = cellanno_col, right_annotation = cellanno_row)
  #   dev.off()
  #   pheatmap(temp, scale = "none")
  # }
  # 
  # 
  # {
  #   for (i in 1:nrow(OsRLCK85_peaks)) {
  #     group <- i
  #     Peak_id <- OsRLCK85_peaks$ID[i]
  #     peakmean <- Peak_mean[,Peak_id]
  #     celltype <- rownames(Peak_mean)
  #     type <- OsRLCK85_peaks$peakGroup[i]
  #     temp <- as.data.frame(rbind(temp,
  #                                 data.frame(celltype = celltype,
  #                                            value = peakmean,
  #                                            type = type,
  #                                            group = group,
  #                                            peak_id = Peak_id)))
  #   }
  #   
  #   temp <- as.data.frame(rbind(temp,
  #                               data.frame(celltype = rownames(RNA_mean),
  #                                          value = RNA_mean[,which(colnames(RNA_mean) == "OsRLCK85")],
  #                                          type = "OsRLCK85",
  #                                          group = 11,
  #                                          peak_id = "OsRLCK85")))
  #   
  #   
  #   
  #   celltype_order <- sort(RNA_mean[,which(colnames(RNA_mean) == "OsRLCK85")], decreasing = T)
  #   temp$celltype <- factor(temp$celltype,
  #                           levels = names(celltype_order)[c(1:3,5,6,4,7:13)])
  #   temp$type <- factor(temp$type,
  #                       levels = c("Gene-linked peaks", "Lineage constitutive peaks",
  #                                  "Global constitutive peaks", "Other peaks", "OsRLCK85"))
  #   temp <- temp[which(temp$peak_id %in%
  #                        c("2_32890766_32891266",
  #                          "2_32878843_32879343",
  #                          "2_32892619_32893119",
  #                          "2_32884668_32885168",
  #                          "OsRLCK85")),]
  #   scale <- TRUE
  #   if (scale == TRUE) {
  #     temp[which(temp$peak_id == "OsRLCK85"), "value"] <- temp[which(temp$peak_id == "OsRLCK85"), "value"] /
  #       max(temp[which(temp$peak_id == "OsRLCK85"), "value"])
  #     temp[which(temp$peak_id != "OsRLCK85"), "value"] <- temp[which(temp$peak_id != "OsRLCK85"), "value"] /
  #       max(temp[which(temp$peak_id != "OsRLCK85"), "value"])
  #   }
  #   
  #   pdf("Figure6_3_OsRLCK85_RNA_Peak_mean.pdf", width = 3.5, height = 3.5)
  #   ggplot(data = temp[which(temp$peak_id == "OsRLCK85"),], aes(x = celltype, y = value, group = group, color = type)) +
  #     geom_line() +
  #     geom_point() +
  #     scale_y_continuous(limits = c(0,1)) +
  #     theme_classic() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #           axis.text = element_text(size = 10, color = "black"),
  #           axis.title = element_text(size = 10, color = "black"),
  #           panel.grid = element_blank()) +
  #     scale_color_manual(values = c("black")) +
  #     labs(x = "", y = "snRNA\nScaled Expression Level")
  #   ggplot(data = temp[-which(temp$peak_id == "OsRLCK85"),], aes(x = celltype, y = value, group = group, color = type)) +
  #     geom_line() +
  #     geom_point() +
  #     scale_y_continuous(limits = c(0,1)) +
  #     theme_classic() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #           axis.text = element_text(size = 10, color = "black"),
  #           axis.title = element_text(size = 10, color = "black"),
  #           panel.grid = element_blank()) +
  #     scale_color_manual(values = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD", "black")) +
  #     labs(x = "", y = "snATAC\nScaled Insertion Count")
  #   dev.off()
  #   
  #   pdf("Figure6_3_OsRLCK85_RNA_Peak_mean_NoLegend.pdf", width = 3.5, height = 3.5)
  #   ggplot(data = temp[which(temp$peak_id == "OsRLCK85"),], aes(x = celltype, y = value, group = group, color = type)) +
  #     geom_line() +
  #     geom_point() +
  #     scale_y_continuous(limits = c(0,1)) +
  #     theme_classic() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #           axis.text = element_text(size = 10, color = "black"),
  #           axis.title = element_text(size = 10, color = "black"),
  #           panel.grid = element_blank(),
  #           legend.position = "none") +
  #     scale_color_manual(values = c("black")) +
  #     labs(x = "", y = "snRNA\nScaled Expression Level")
  #   ggplot(data = temp[-which(temp$peak_id == "OsRLCK85"),], aes(x = celltype, y = value, group = group, color = type)) +
  #     geom_line() +
  #     geom_point() +
  #     scale_y_continuous(limits = c(0,1)) +
  #     theme_classic() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #           axis.text = element_text(size = 10, color = "black"),
  #           axis.title = element_text(size = 10, color = "black"),
  #           panel.grid = element_blank(),
  #           legend.position = "none") +
  #     scale_color_manual(values = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD", "black")) +
  #     labs(x = "", y = "snATAC\nScaled Insertion Count")
  #   dev.off()
}

### Histone modification
{
  #### Add peak annotation
  H3K23ac_Osj_leaf_normal <- readBed(file = paste0("./Ref/",
                                                   "Histone modifications (HMs, ChIP-seq)_H3K23ac_Osj_leaf_normal.bed"))
  H3K23ac_Osj_leaf_normal@seqinfo@seqnames <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
  H3K23ac_Osj_leaf_normal@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                                    levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  names(H3K23ac_Osj_leaf_normal@elementMetadata@listData) <- c("int(m10Log10_P)",
                                                               "name",
                                                               "peak_foldchange",
                                                               "mLog10_P",
                                                               "mLog10_Adj")
  
  H4K16ac_Osj_leaf_normal <- readBed(file = paste0("./Ref/",
                                                   "Histone modifications (HMs, ChIP-seq)_H4K16ac_Osj_leaf_normal.bed"))
  H4K16ac_Osj_leaf_normal@seqinfo@seqnames <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
  H4K16ac_Osj_leaf_normal@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                                    levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  names(H4K16ac_Osj_leaf_normal@elementMetadata@listData) <- c("int(m10Log10_P)",
                                                               "name",
                                                               "peak_foldchange",
                                                               "mLog10_P",
                                                               "mLog10_Adj")
  H4K16ac_Osj_leaf_normal@elementMetadata@listData[["Type"]] <- rep("H4K16ac_leaf",length(H4K16ac_Osj_leaf_normal@elementMetadata@listData[["mLog10_Adj"]]))
  
  
  H3K36me3_Osj_seedling_normal <- readBed(file = paste0("./Ref/",
                                                        "Histone modifications (HMs, ChIP-seq)_H3K36me3_Osj_seedling_normal.bed"))
  H3K36me3_Osj_seedling_normal@seqinfo@seqnames <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
  H3K36me3_Osj_seedling_normal@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                                         levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  names(H3K36me3_Osj_seedling_normal@elementMetadata@listData) <- c("int(m10Log10_P)",
                                                                    "name",
                                                                    "peak_foldchange",
                                                                    "mLog10_P",
                                                                    "mLog10_Adj")
  H3K36me3_Osj_seedling_normal@elementMetadata@listData[["Type"]] <- rep("H3K36me3_seedling",length(H3K36me3_Osj_seedling_normal@elementMetadata@listData[["mLog10_Adj"]]))
  
  
  H3K4me2_Osj_seedling_normal <- readBed(file = paste0("./Ref/",
                                                       "Histone modifications (HMs, ChIP-seq)_H3K4me2_Osj_seedling_normal.bed"))
  H3K4me2_Osj_seedling_normal@seqinfo@seqnames <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
  H3K4me2_Osj_seedling_normal@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                                        levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  names(H3K4me2_Osj_seedling_normal@elementMetadata@listData) <- c("int(m10Log10_P)",
                                                                   "name",
                                                                   "peak_foldchange",
                                                                   "mLog10_P",
                                                                   "mLog10_Adj")
  H3K4me2_Osj_seedling_normal@elementMetadata@listData[["Type"]] <- rep("H3K4me2_seedling",length(H3K4me2_Osj_seedling_normal@elementMetadata@listData[["mLog10_Adj"]]))
  
  
  H4K12ac_Osj_seedling_normal <- readBed(file = paste0("./Ref/",
                                                       "Histone modifications (HMs, ChIP-seq)_H4K12ac_Osj_seedling_normal.bed"))
  H4K12ac_Osj_seedling_normal@seqinfo@seqnames <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
  H4K12ac_Osj_seedling_normal@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                                        levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  names(H4K12ac_Osj_seedling_normal@elementMetadata@listData) <- c("int(m10Log10_P)",
                                                                   "name",
                                                                   "peak_foldchange",
                                                                   "mLog10_P",
                                                                   "mLog10_Adj")
  H4K12ac_Osj_seedling_normal@elementMetadata@listData[["Type"]] <- rep("H4K12ac_seedling",length(H4K12ac_Osj_seedling_normal@elementMetadata@listData[["mLog10_Adj"]]))
  
  
  ATAC_common_tissues <- addPeakAnnotations(
    ArchRProj = ATAC_common_tissues,
    regions = list("H3K23ac_leaf" = H3K23ac_Osj_leaf_normal,
                   "H4K16ac_leaf" = H4K16ac_Osj_leaf_normal,
                   "H3K36me3_seedling" = H3K36me3_Osj_seedling_normal,
                   "H3K4me2_seedling" = H3K4me2_Osj_seedling_normal,
                   "H4K12ac_seedling" = H4K12ac_Osj_seedling_normal),
    name = "Histone_modifications",
    force = TRUE
  )
  
  His_peak_anno <- getPeakAnnotation(ATAC_common_tissues, name = "Histone_modifications")
  His_peak_anno_match <- readRDS(His_peak_anno[["Matches"]])
  His_peak_anno_match <- His_peak_anno_match@assays@data@listData[["matches"]]
  His_peak_anno_match <- as.matrix(His_peak_anno_match)
  His_peak_anno_match <- His_peak_anno_match + 0
  His_peak_anno_match <- as.data.frame(His_peak_anno_match)
  rownames(His_peak_anno_match) <- rownames(peaksets)
  
  peaks_group <- as.data.frame(cbind(peaks_group,
                                     His_peak_anno_match[peaks_group$ID,]))
  
  colnames(peaks_group)
  
  peaks_group_his <- peaks_group[,c("peakGroup", "peakType", "H3K23ac_leaf", "H4K16ac_leaf",
                                    "H3K36me3_seedling", "H3K4me2_seedling", "H4K12ac_seedling")]
  peaks_group_his <- peaks_group_his[order(rowSums(peaks_group_his[,3:ncol(peaks_group_his)])),]
  
  table(peaks_group_meth$peakGroup)
  
  pdf("Figure6_3_histo_anno_heatmap_Gene_linked.pdf", width = 2.5, height = 5.5)
  length(intersect(which(peaks_group_his$peakGroup == "Gene-linked peaks"),
                   which(peaks_group_his$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_his[intersect(which(peaks_group_his$peakGroup == "Gene-linked peaks"),
                                    which(peaks_group_his$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_his)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#A9C456", "#C25568"))(5))
  dev.off()
  
  pdf("Figure6_3_histo_anno_heatmap_Other_peaks.pdf", width = 2.5, height = 5.5)
  length(intersect(which(peaks_group_his$peakGroup == "Other peaks"),
                   which(peaks_group_his$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_his[intersect(which(peaks_group_his$peakGroup == "Other peaks"),
                                    which(peaks_group_his$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_his)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#A9C456", "#C25568"))(5))
  dev.off()
  
  pdf("Figure6_3_histo_anno_heatmap_Global_constitutive_peaks.pdf", width = 2.5, height = 5.5)
  length(intersect(which(peaks_group_his$peakGroup == "Global constitutive peaks"),
                   which(peaks_group_his$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_his[intersect(which(peaks_group_his$peakGroup == "Global constitutive peaks"),
                                    which(peaks_group_his$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_his)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#A9C456", "#C25568"))(5))
  dev.off()
  
  pdf("Figure6_3_histo_anno_heatmap_Celltype_constitutive_peaks.pdf", width = 2.5, height = 5.5)
  length(intersect(which(peaks_group_his$peakGroup == "Lineage constitutive peaks"),
                   which(peaks_group_his$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_his[intersect(which(peaks_group_his$peakGroup == "Lineage constitutive peaks"),
                                    which(peaks_group_his$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_his)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#A9C456", "#C25568"))(5))
  dev.off()
  
  
  peaks_group_his_2 <- peaks_group_his[which(peaks_group_his$peakType %in% c("Promoter", "Distal")),]
  group_num <- as.data.frame(table(peaks_group_his_2$peakGroup))
  rownames(group_num) <- group_num$Var1
  peaks_group_his_2_sum <- rowSums(peaks_group_his_2[,3:ncol(peaks_group_his_2)])
  peaks_group_his_DF <- as.data.frame(table(peaks_group_his_2$peakGroup, peaks_group_his_2_sum))
  peaks_group_his_DF <- peaks_group_his_DF[order(peaks_group_his_DF$Var1),]
  peaks_group_his_DF$Freq <- peaks_group_his_DF$Freq / group_num[peaks_group_his_DF$Var1,2]
  
  pdf("Figure6_3_hit_barplot_diff.pdf", width = 7, height = 5.5)
  ggplot() +
    geom_bar(data = peaks_group_his_DF, aes(x = peaks_group_his_2_sum, y = Freq, fill = Var1),
             stat = "identity", position = "dodge", color = "black", width = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = peak_type_color) +
    labs(x = "Annotated histone modification type count", y = "Percent", fill = "Peak Type")+ 
    theme_bw(base_size = 10)+  
    theme(axis.text = element_text(colour = 'black'))
  dev.off()
  
  
}

### Methylation modification
{
  readGeneric2 <- function(file, chr=1,start=2,end=3,strand=NULL,meta.cols=NULL, 
                           keep.all.metadata=FALSE, zero.based=FALSE, 
                           remove.unusual=FALSE, header=FALSE, 
                           skip=0, sep="\t"){
    
    # reads the bed files
    df=readTableFast(file, header=header, skip=skip, sep=sep, chr=chr)                    
    
    # make a list of new column names, and their column numbers
    col.names1=list(chr=chr,start=start,end=end,strand=strand)
    col.names=c(col.names1,meta.cols) # put the meta colums if any
    
    # check if col number exceeds dimensions of the original df.
    if( max(unlist(col.names)) > ncol(df) ) 
      stop("Number of columns is lower than designated number of columns by ",
           "meta.cols,chr,start,end or strand arguments\n")
    
    # change the col names to the ones given by meta.cols and chr,str,end,strand
    colnames(df)[unlist(col.names)] = names(unlist(col.names))
    
    # converts the . strand character to *
    sind = grepl('strand',colnames(df))
    if(any(sind) & !is.null(strand))
      df[, sind] = sub('\\.', '*', df[,sind])
    
    # removes nonstandard chromosome names
    if(remove.unusual)
      df = df[grep("_", as.character(df$chr),invert=TRUE),]
    
    g = makeGRangesFromDataFrame(
      df, 
      keep.extra.columns = FALSE, 
      starts.in.df.are.0based = zero.based,
      ignore.strand = is.null(strand),
      na.rm = TRUE)
    
    # this names can not be column names in meta-data
    black.names=c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                  "isCircular", "start", "end", "width", "element")
    
    
    if(keep.all.metadata){
      my.mcols = df[,-unlist(col.names1),drop=FALSE]
      mcols(g) = my.mcols[, !colnames(my.mcols) %in% black.names]
    }else if(!is.null(meta.cols)){
      my.mcols=df[,unlist(meta.cols),drop=FALSE]
      values(g) = my.mcols[, !colnames(my.mcols) %in% black.names, drop=FALSE]
    }
    
    return(g)
  } # 修改readGeneric函数，在makeGRangesFromDataFrame函数设置'na.rm=TRUE'
  readBed2 <- function(file,track.line=FALSE,
                       remove.unusual=FALSE,zero.based=TRUE){
    
    meta.cols=list(score=5,name=4,thickStart=7,  
                   thickEnd=8, 
                   itemRgb=9,  
                   blockCount=10, 
                   blockSizes=11,  
                   blockStarts=12 )
    
    file <- compressedAndUrl2temp(file)
    if(is.numeric(track.line)){
      skip = track.line
    }else if(track.line=="auto"){
      skip = detectUCSCheader(file)
    }else{
      skip = 0
    }
    
    df=read_delim(file,skip=skip,n_max=2,col_names=FALSE, delim="\t")
    numcol=ncol(df)
    
    if(numcol==3){
      df=readGeneric2(file, chr = 1, start = 2, end = 3, strand = NULL,
                      meta.cols = NULL,   zero.based = zero.based,
                      remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")
    }else if(numcol==4){
      df=readGeneric2(file, chr = 1, start = 2, end = 3, strand = NULL,
                      meta.cols = meta.cols[1],   zero.based = zero.based,
                      remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")    
    }else if(numcol==5){
      df=readGeneric2(file, chr = 1, start = 2, end = 3, strand = NULL,
                      meta.cols = meta.cols[1:2],   zero.based = zero.based,
                      remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")    
    }else if(numcol == 6){
      df=readGeneric2(file, chr = 1, start = 2, end = 3, strand = 6,
                      meta.cols = meta.cols[1:2],   zero.based = zero.based,
                      remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t") 
    }else if(numcol > 6){
      df=readGeneric2(file, chr = 1, start = 2, end = 3, strand = 6,
                      meta.cols = meta.cols[c(1:2,(3:8)[1:(numcol-6)])],   zero.based = zero.based,
                      remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t") 
    }
    df
    
  } 
  Oryza_sativa_Leaf_3w_WT2003 <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                            "Oryza_sativa_Leaf_3w_WT2003.bed.gz"),
                                              skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                     depth = 6,
                                                                                     mCdep = 7,
                                                                                     level = 8))
  Oryza_sativa_Leaf_3w_WT2007 <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                            "Oryza_sativa_Leaf_3w_WT2007.bed.gz"),
                                              skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                     depth = 6,
                                                                                     mCdep = 7,
                                                                                     level = 8))
  Oryza_sativa_Leaf_3w_WT2011 <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                            "Oryza_sativa_Leaf_3w_WT2011.bed.gz"),
                                              skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                     depth = 6,
                                                                                     mCdep = 7,
                                                                                     level = 8))
  Oryza_sativa_Seed_12h <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                      "Oryza_sativa_Seed_12h.bed.gz"),
                                        skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                               depth = 6,
                                                                               mCdep = 7,
                                                                               level = 8))
  Oryza_sativa_Seed_DRY <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                      "Oryza_sativa_Seed_DRY.bed.gz"),
                                        skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                               depth = 6,
                                                                               mCdep = 7,
                                                                               level = 8))
  Oryza_sativa_Seed_SDLG <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                       "Oryza_sativa_Seed_SDLG.bed.gz"),
                                         skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                depth = 6,
                                                                                mCdep = 7,
                                                                                level = 8))
  Oryza_sativa_Shoot_11d_WT <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                          "Oryza_sativa_Shoot_11d_WT.bed.gz"),
                                            skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                   depth = 6,
                                                                                   mCdep = 7,
                                                                                   level = 8))
  Oryza_sativa_Leaf_3w_WT2011_replicate <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                                      "Oryza_sativa_Leaf_3w_WT2011_replicate.bed.gz"),
                                                        skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                               depth = 6,
                                                                                               mCdep = 7,
                                                                                               level = 8))
  Oryza_sativa_Leaf_3w_WT_regenerated <- readGeneric2(file = paste0("./Ref/MethBank/Methylomes/",
                                                                    "Oryza_sativa_Leaf_3w_WT_regenerated.bed.gz"),
                                                      skip = 1, strand = 4, meta.cols = list(mCtype = 5,
                                                                                             depth = 6,
                                                                                             mCdep = 7,
                                                                                             level = 8))
  Oryza_sativa_Leaf_3w_WT_regenerated <- Oryza_sativa_Leaf_3w_WT_regenerated[which(Oryza_sativa_Leaf_3w_WT_regenerated$level == 1),]
  Oryza_sativa_Leaf_3w_WT_regenerated@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Leaf_3w_WT_regenerated@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                                levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Leaf_3w_WT_regenerated_CG <- Oryza_sativa_Leaf_3w_WT_regenerated[which(Oryza_sativa_Leaf_3w_WT_regenerated$mCtype == "CG"),]
  Oryza_sativa_Leaf_3w_WT_regenerated_CHG <- Oryza_sativa_Leaf_3w_WT_regenerated[which(Oryza_sativa_Leaf_3w_WT_regenerated$mCtype == "CHG"),]
  Oryza_sativa_Leaf_3w_WT_regenerated_CHH <- Oryza_sativa_Leaf_3w_WT_regenerated[which(Oryza_sativa_Leaf_3w_WT_regenerated$mCtype == "CHH"),]
  
  Oryza_sativa_Leaf_3w_WT2011_replicate <- Oryza_sativa_Leaf_3w_WT2011_replicate[which(Oryza_sativa_Leaf_3w_WT2011_replicate$level == 1),]
  Oryza_sativa_Leaf_3w_WT2011_replicate@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Leaf_3w_WT2011_replicate@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                                  levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Leaf_3w_WT2011_replicate_CG <- Oryza_sativa_Leaf_3w_WT2011_replicate[which(Oryza_sativa_Leaf_3w_WT2011_replicate$mCtype == "CG"),]
  Oryza_sativa_Leaf_3w_WT2011_replicate_CHG <- Oryza_sativa_Leaf_3w_WT2011_replicate[which(Oryza_sativa_Leaf_3w_WT2011_replicate$mCtype == "CHG"),]
  Oryza_sativa_Leaf_3w_WT2011_replicate_CHH <- Oryza_sativa_Leaf_3w_WT2011_replicate[which(Oryza_sativa_Leaf_3w_WT2011_replicate$mCtype == "CHH"),]
  
  Oryza_sativa_Shoot_11d_WT <- Oryza_sativa_Shoot_11d_WT[which(Oryza_sativa_Shoot_11d_WT$level == 1),]
  Oryza_sativa_Shoot_11d_WT@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Shoot_11d_WT@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                      levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Shoot_11d_WT_CG <- Oryza_sativa_Shoot_11d_WT[which(Oryza_sativa_Shoot_11d_WT$mCtype == "CG"),]
  Oryza_sativa_Shoot_11d_WT_CHG <- Oryza_sativa_Shoot_11d_WT[which(Oryza_sativa_Shoot_11d_WT$mCtype == "CHG"),]
  Oryza_sativa_Shoot_11d_WT_CHH <- Oryza_sativa_Shoot_11d_WT[which(Oryza_sativa_Shoot_11d_WT$mCtype == "CHH"),]
  
  Oryza_sativa_Seed_SDLG <- Oryza_sativa_Seed_SDLG[which(Oryza_sativa_Seed_SDLG$level == 1),]
  Oryza_sativa_Seed_SDLG@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Seed_SDLG@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                   levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Seed_SDLG_CG <- Oryza_sativa_Seed_SDLG[which(Oryza_sativa_Seed_SDLG$mCtype == "CG"),]
  Oryza_sativa_Seed_SDLG_CHG <- Oryza_sativa_Seed_SDLG[which(Oryza_sativa_Seed_SDLG$mCtype == "CHG"),]
  Oryza_sativa_Seed_SDLG_CHH <- Oryza_sativa_Seed_SDLG[which(Oryza_sativa_Seed_SDLG$mCtype == "CHH"),]
  
  Oryza_sativa_Seed_DRY <- Oryza_sativa_Seed_DRY[which(Oryza_sativa_Seed_DRY$level == 1),]
  Oryza_sativa_Seed_DRY@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Seed_DRY@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                  levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Seed_DRY_CG <- Oryza_sativa_Seed_DRY[which(Oryza_sativa_Seed_DRY$mCtype == "CG"),]
  Oryza_sativa_Seed_DRY_CHG <- Oryza_sativa_Seed_DRY[which(Oryza_sativa_Seed_DRY$mCtype == "CHG"),]
  Oryza_sativa_Seed_DRY_CHH <- Oryza_sativa_Seed_DRY[which(Oryza_sativa_Seed_DRY$mCtype == "CHH"),]
  
  Oryza_sativa_Seed_12h <- Oryza_sativa_Seed_12h[which(Oryza_sativa_Seed_12h$level == 1),]
  Oryza_sativa_Seed_12h@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Seed_12h@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                  levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Seed_12h_CG <- Oryza_sativa_Seed_12h[which(Oryza_sativa_Seed_12h$mCtype == "CG"),]
  Oryza_sativa_Seed_12h_CHG <- Oryza_sativa_Seed_12h[which(Oryza_sativa_Seed_12h$mCtype == "CHG"),]
  Oryza_sativa_Seed_12h_CHH <- Oryza_sativa_Seed_12h[which(Oryza_sativa_Seed_12h$mCtype == "CHH"),]
  
  Oryza_sativa_Leaf_3w_WT2011 <- Oryza_sativa_Leaf_3w_WT2011[which(Oryza_sativa_Leaf_3w_WT2011$level == 1),]
  Oryza_sativa_Leaf_3w_WT2011@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Leaf_3w_WT2011@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                        levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Leaf_3w_WT2011_CG <- Oryza_sativa_Leaf_3w_WT2011[which(Oryza_sativa_Leaf_3w_WT2011$mCtype == "CG"),]
  Oryza_sativa_Leaf_3w_WT2011_CHG <- Oryza_sativa_Leaf_3w_WT2011[which(Oryza_sativa_Leaf_3w_WT2011$mCtype == "CHG"),]
  Oryza_sativa_Leaf_3w_WT2011_CHH <- Oryza_sativa_Leaf_3w_WT2011[which(Oryza_sativa_Leaf_3w_WT2011$mCtype == "CHH"),]
  
  Oryza_sativa_Leaf_3w_WT2007 <- Oryza_sativa_Leaf_3w_WT2007[which(Oryza_sativa_Leaf_3w_WT2007$level == 1),]
  Oryza_sativa_Leaf_3w_WT2007@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Leaf_3w_WT2007@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                        levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Leaf_3w_WT2007_CG <- Oryza_sativa_Leaf_3w_WT2007[which(Oryza_sativa_Leaf_3w_WT2007$mCtype == "CG"),]
  Oryza_sativa_Leaf_3w_WT2007_CHG <- Oryza_sativa_Leaf_3w_WT2007[which(Oryza_sativa_Leaf_3w_WT2007$mCtype == "CHG"),]
  Oryza_sativa_Leaf_3w_WT2007_CHH <- Oryza_sativa_Leaf_3w_WT2007[which(Oryza_sativa_Leaf_3w_WT2007$mCtype == "CHH"),]
  
  Oryza_sativa_Leaf_3w_WT2003 <- Oryza_sativa_Leaf_3w_WT2003[which(Oryza_sativa_Leaf_3w_WT2003$level == 1),]
  Oryza_sativa_Leaf_3w_WT2003@seqinfo@seqnames <- c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  Oryza_sativa_Leaf_3w_WT2003@seqnames@values <- factor(c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "C", "M"),
                                                        levels = c("C","10","11", "12", "M", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  Oryza_sativa_Leaf_3w_WT2003_CG <- Oryza_sativa_Leaf_3w_WT2003[which(Oryza_sativa_Leaf_3w_WT2003$mCtype == "CG"),]
  Oryza_sativa_Leaf_3w_WT2003_CHG <- Oryza_sativa_Leaf_3w_WT2003[which(Oryza_sativa_Leaf_3w_WT2003$mCtype == "CHG"),]
  Oryza_sativa_Leaf_3w_WT2003_CHH <- Oryza_sativa_Leaf_3w_WT2003[which(Oryza_sativa_Leaf_3w_WT2003$mCtype == "CHH"),]
  
  
  CHH <- rbindlist(list(as.data.frame(Oryza_sativa_Leaf_3w_WT2003_CHH),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT2007_CHH),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT2011_CHH),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT2011_replicate_CHH),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT_regenerated_CHH),
                        as.data.frame(Oryza_sativa_Seed_12h_CHH),
                        as.data.frame(Oryza_sativa_Seed_DRY_CHH),
                        as.data.frame(Oryza_sativa_Seed_SDLG_CHH),
                        as.data.frame(Oryza_sativa_Shoot_11d_WT_CHH)
  ))
  CHH <- makeGRangesFromDataFrame(CHH,
                                  keep.extra.columns = TRUE, seqnames.field = "seqnames",
                                  start.field = "start", end.field = "end", strand.field = "strand")
  CHG <- rbindlist(list(as.data.frame(Oryza_sativa_Leaf_3w_WT2003_CHG),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT2007_CHG),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT2011_CHG),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT2011_replicate_CHG),
                        as.data.frame(Oryza_sativa_Leaf_3w_WT_regenerated_CHG),
                        as.data.frame(Oryza_sativa_Seed_12h_CHG),
                        as.data.frame(Oryza_sativa_Seed_DRY_CHG),
                        as.data.frame(Oryza_sativa_Seed_SDLG_CHG),
                        as.data.frame(Oryza_sativa_Shoot_11d_WT_CHG)
  ))
  CHG <- makeGRangesFromDataFrame(CHG,
                                  keep.extra.columns = TRUE, seqnames.field = "seqnames",
                                  start.field = "start", end.field = "end", strand.field = "strand")
  CG <- rbindlist(list(as.data.frame(Oryza_sativa_Leaf_3w_WT2003_CG),
                       as.data.frame(Oryza_sativa_Leaf_3w_WT2007_CG),
                       as.data.frame(Oryza_sativa_Leaf_3w_WT2011_CG),
                       as.data.frame(Oryza_sativa_Leaf_3w_WT2011_replicate_CG),
                       as.data.frame(Oryza_sativa_Leaf_3w_WT_regenerated_CG),
                       as.data.frame(Oryza_sativa_Seed_12h_CG),
                       as.data.frame(Oryza_sativa_Seed_DRY_CG),
                       as.data.frame(Oryza_sativa_Seed_SDLG_CG),
                       as.data.frame(Oryza_sativa_Shoot_11d_WT_CG)
  ))
  CG <- makeGRangesFromDataFrame(CG,
                                 keep.extra.columns = TRUE, seqnames.field = "seqnames",
                                 start.field = "start", end.field = "end", strand.field = "strand")
  ATAC_common_tissues <- addPeakAnnotations(
    ArchRProj = ATAC_common_tissues,
    regions = list("CHH" = CHH,
                   "CG" = CG,
                   "CHG" = CHG),
    name = "Meth_modifications",
    force = TRUE
  )
  
  Meth_peak_anno <- getPeakAnnotation(ATAC_common_tissues, name = "Meth_modifications")
  Meth_peak_anno_match <- readRDS(Meth_peak_anno[["Matches"]])
  Meth_peak_anno_match <- Meth_peak_anno_match@assays@data@listData[["matches"]]
  Meth_peak_anno_match <- as.matrix(Meth_peak_anno_match)
  Meth_peak_anno_match <- Meth_peak_anno_match + 0
  Meth_peak_anno_match <- as.data.frame(Meth_peak_anno_match)
  rownames(Meth_peak_anno_match) <- rownames(peaksets)
  
  peaks_group <- as.data.frame(cbind(peaks_group,
                                     Meth_peak_anno_match[peaks_group$ID,]))
  
  colnames(peaks_group)
  peaks_group_meth <- peaks_group[,c("peakGroup", "peakType", "CHH", "CG", "CHG")]
  peaks_group_meth <- peaks_group_meth[order(rowSums(peaks_group_meth[,3:ncol(peaks_group_meth)])),]
  
  pdf("Figure6_3_meth_anno_heatmap_Gene_linked.pdf", width = 1.5, height = 4.5)
  length(intersect(which(peaks_group_meth$peakGroup == "Gene-linked peaks"),
                   which(peaks_group_meth$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_meth[intersect(which(peaks_group_meth$peakGroup == "Gene-linked peaks"),
                                     which(peaks_group_meth$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_meth)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#00ACC1", "#E53935"))(5))
  dev.off()
  
  pdf("Figure6_3_meth_anno_heatmap_Other_peaks.pdf", width = 1.5, height = 4.5)
  length(intersect(which(peaks_group_meth$peakGroup == "Other peaks"),
                   which(peaks_group_meth$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_meth[intersect(which(peaks_group_meth$peakGroup == "Other peaks"),
                                     which(peaks_group_meth$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_meth)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#00ACC1", "#E53935"))(5))
  dev.off()
  
  pdf("Figure6_3_meth_anno_heatmap_Global_constitutive_peaks.pdf", width = 1.5, height = 4.5)
  length(intersect(which(peaks_group_meth$peakGroup == "Global constitutive peaks"),
                   which(peaks_group_meth$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_meth[intersect(which(peaks_group_meth$peakGroup == "Global constitutive peaks"),
                                     which(peaks_group_meth$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_meth)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#00ACC1", "#E53935"))(5))
  dev.off()
  
  pdf("Figure6_3_meth_anno_heatmap_Celltype_constitutive_peaks.pdf", width = 1.5, height = 4.5)
  length(intersect(which(peaks_group_meth$peakGroup == "Lineage constitutive peaks"),
                   which(peaks_group_meth$peakType %in% c("Promoter","Distal"))))
  Heatmap(peaks_group_meth[intersect(which(peaks_group_meth$peakGroup == "Lineage constitutive peaks"),
                                     which(peaks_group_meth$peakType %in% c("Promoter","Distal"))),3:ncol(peaks_group_meth)],
          cluster_rows = F, use_raster = T, raster_quality = 20,
          show_row_names = F, raster_by_magick = F, col = colorRampPalette(c("#00ACC1", "#E53935"))(5))
  dev.off()
  
  peaks_group_meth_2 <- peaks_group_meth[which(peaks_group_meth$peakType %in% c("Promoter", "Distal")),]
  table(peaks_group_meth_2$peakType)
  group_num <- as.data.frame(table(peaks_group_meth_2$peakGroup))
  rownames(group_num) <- group_num$Var1
  peaks_group_meth_2_sum <- rowSums(peaks_group_meth_2[,3:ncol(peaks_group_meth_2)])
  peaks_group_meth_DF <- as.data.frame(table(peaks_group_meth_2$peakGroup, peaks_group_meth_2_sum))
  peaks_group_meth_DF <- peaks_group_meth_DF[order(peaks_group_meth_DF$Var1),]
  peaks_group_meth_DF$Freq <- peaks_group_meth_DF$Freq / group_num[peaks_group_meth_DF$Var1,2]
  
  pdf("Figure6_3_meth_barplot_diff.pdf", width = 6, height = 5.5)
  ggplot() +
    geom_bar(data = peaks_group_meth_DF, aes(x = peaks_group_meth_2_sum, y = Freq, fill = Var1),
             stat = "identity", position = "dodge", color = "black", width = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = peak_type_color) +
    labs(x = "Annotated methylated modification type count", y = "Percent", fill = "Peak Type")+ 
    theme_bw(base_size = 10)+  
    theme(axis.text = element_text(colour = 'black'))
  dev.off()
  
}



### Gene_Linked_Peaks gene & celltype DEG
{
  Idents(RNA_common_tissues) <- RNA_common_tissues$Celltype
  RNA_Celltye_DEG <- FindAllMarkers(RNA_common_tissues, logfc.threshold = 0.5, only.pos = T)
  sum(RNA_Celltye_DEG$gene %in% proj_pass_filter@geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  # RNA_Celltye_DEG <- RNA_Celltye_DEG[which(RNA_Celltye_DEG$gene %in% proj_pass_filter@geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]]),]
  RNA_Celltye_DEG <- RNA_Celltye_DEG[which(RNA_Celltye_DEG$p_val_adj < 0.05),]
  length(unique(RNA_Celltye_DEG$gene))
  
  Gene_Linked_Peaks_genes <- peaksets_2[Gene_Linked_Peaks, "nearestGene"]
  Gene_Linked_Peaks_genes <- unique(Gene_Linked_Peaks_genes)
  
  sum(Gene_Linked_Peaks_genes %in% RNA_Celltye_DEG$gene) # 1775
}




