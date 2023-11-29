setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas")
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
library(monocle)
library(tidyverse)
library(patchwork)

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


RNA <- readRDS("./Rice-snRNA/all_data_annotated.rds")
unique(RNA$Celltype)
RNA$Celltype[which(RNA$Celltype == "Inflorescence meristem(IM)")] <- "Inflorescence meristem (IM)"
RNA$Celltype[which(RNA$Celltype == "Spikelet meristem(SM)")] <- "Spikelet meristem (SM)"
RNA$Celltype[which(RNA$Celltype == "Branch meristems(BM)")] <- "Branch meristems (BM)"
RNA$Celltype[which(RNA$Celltype == "Lemma(le)")] <- "Lemma (le)"
RNA$Celltype[which(RNA$Celltype == "Cryptic bract/bract(cb/b)")] <- "Cryptic bract/bract (cb/b)"
unique(RNA$Celltype)
sum(is.na(RNA$Celltype))

Idents(RNA) <- RNA$Celltype

gene_id_map <- readRDS("gene_id_map.rds")
rownames(gene_id_map) <- gene_id_map$gene_id
RNA_id_gene <- data.frame(id = rownames(RNA),
                          symbol = gene_id_map[rownames(RNA),2],
                          row.names = rownames(RNA))
RNA_id_gene <- as.data.frame(na.omit(RNA_id_gene))
RNA <- RNA[RNA_id_gene$id,]
RNA@assays$RNA@counts@Dimnames[[1]] <- RNA_id_gene$symbol
RNA@assays$RNA@data@Dimnames[[1]] <- RNA_id_gene$symbol
rownames(RNA@assays$RNA@scale.data) <- RNA_id_gene[rownames(RNA@assays$RNA@scale.data),"symbol"]
RNA <- CreateSeuratObject(counts = RNA@assays$RNA@counts,
                          meta.data = RNA@meta.data)
RNA <- NormalizeData(object = RNA)


proj_pass_filter <- readRDS("3.proj_pass_filter_withPeaks.rds")
unique(proj_pass_filter$Celltype)
unique(proj_pass_filter$TissuesSub)

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
momr::plotCol(celltype_color)
names(celltype_color) <- sort(unique(proj_pass_filter$Celltype))

########## IM

markersPeaks_Celltype <- getMarkerFeatures(
  ArchRProj = proj_pass_filter, 
  useMatrix = "PeakMatrix", 
  groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerPeaksList_Celltype <- getMarkers(markersPeaks_Celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
BM_marker_peaks <- markerPeaksList_Celltype@listData[["Branch meristems (BM)"]]
BM_marker_peaks <- as.data.frame(BM_marker_peaks)
rownames(BM_marker_peaks) <- paste(BM_marker_peaks$seqnames,
                                   BM_marker_peaks$start,
                                   BM_marker_peaks$end,
                                   sep="_")
BM_marker_peaks$ID <- rownames(BM_marker_peaks)

IM_marker_peaks <- markerPeaksList_Celltype@listData[["Inflorescence meristem (IM)"]]
IM_marker_peaks <- as.data.frame(IM_marker_peaks)
rownames(IM_marker_peaks) <- paste(IM_marker_peaks$seqnames,
                                   IM_marker_peaks$start,
                                   IM_marker_peaks$end,
                                   sep="_")
IM_marker_peaks$ID <- rownames(IM_marker_peaks)

peaksets <- getPeakSet(proj_pass_filter)
pr1 <- paste(seqnames(peaksets),start(peaksets),end(peaksets),sep="_")
peaksets <- as.data.frame(peaksets@elementMetadata)
rownames(peaksets) <- pr1

BM_marker_peaks$peakType <- peaksets[BM_marker_peaks$ID, "peakType"]
IM_marker_peaks$peakType <- peaksets[IM_marker_peaks$ID, "peakType"]

{
  peak_region_color <- c("#ADC698", "#A0C1D1",
                         "#B18ED7", "#D676A1")
  names(peak_region_color) <- c("Exonic", "Intronic", "Distal", "Promoter")
  temp <- as.data.frame(table(BM_marker_peaks$peakType))
  colnames(temp) <- c("PeakRegion", "Count")
  temp$Color <- peak_region_color[temp$PeakRegion]
  temp$Label <- c("Distal(n=2,983)", "Exonic(n=284)", "Intronic(n=407)", "Promoter(n=1,658)")
  
  pdf("Figure5_2_BM_marker_peaks_percent.pdf", width = 4.5, height = 4.5)
  pie(temp$Count, labels = temp$Label, col = temp$Color)
  dev.off()
  
  table(IM_marker_peaks$peakType)
  temp <- as.data.frame(table(IM_marker_peaks$peakType))
  colnames(temp) <- c("PeakRegion", "Count")
  temp$Color <- peak_region_color[temp$PeakRegion]
  temp$Label <- c("Distal(n=3,301)", "Exonic(n=417)", "Intronic(n=376)", "Promoter(n=1,923)")
  
  pdf("Figure5_2_IM_marker_peaks_percent.pdf", width = 4.5, height = 4.5)
  pie(temp$Count, labels = temp$Label, col = temp$Color)
  dev.off()
}

BM_marker_peaks_PD <- BM_marker_peaks[which(BM_marker_peaks$peakType %in% c("Promoter",
                                                                            "Distal")),]
IM_marker_peaks_PD <- IM_marker_peaks[which(IM_marker_peaks$peakType %in% c("Promoter",
                                                                            "Distal")),]

proj_pass_filter <- addMotifAnnotations(ArchRProj = proj_pass_filter, motifPWMs = motif_all,
                                        annoName = "TF-Motif", force = T)
matches <- getMatches(proj_pass_filter, name = "TF-Motif")
r1 <- SummarizedExperiment::rowRanges(matches)
pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
rownames(matches) <- pr1

peaksets_pd <- peaksets[which(peaksets$peakType %in% c("Promoter",
                                                       "Distal")),]

temp <- matrix(FALSE, nrow = nrow(peaksets_pd), ncol = 2)
rownames(temp) <- rownames(peaksets_pd)
colnames(temp) <- c("Branch meristems (BM)", "Inflorescence meristem (IM)")
temp[rownames(BM_marker_peaks_PD), 1] <- TRUE
temp[rownames(IM_marker_peaks_PD), 2] <- TRUE

matches_pd <- matches[rownames(temp),]
identical(rownames(temp), rownames(matches_pd))
BM_Peaks_Enrich <- .computeEnrichment(matches_pd, which(temp[,1]), 1:nrow(matches_pd))
BM_Peaks_Enrich$Enrichment_log <- log2(BM_Peaks_Enrich$Enrichment)
BM_Peaks_Enrich$Enrichmented <- "NO"
BM_Peaks_Enrich$Enrichmented[which(BM_Peaks_Enrich$Enrichment_log >= 0.25 & BM_Peaks_Enrich$mlog10Padj >= 3)] <- "YES"
table(BM_Peaks_Enrich$Enrichmented)

IM_Peaks_Enrich <- .computeEnrichment(matches_pd, which(temp[,2]), 1:nrow(matches_pd))
IM_Peaks_Enrich$Enrichment_log <- log2(IM_Peaks_Enrich$Enrichment)
IM_Peaks_Enrich$Enrichmented <- "NO"
IM_Peaks_Enrich$Enrichmented[which(IM_Peaks_Enrich$Enrichment_log >= 0.25 & IM_Peaks_Enrich$mlog10Padj >= 3)] <- "YES"
table(IM_Peaks_Enrich$Enrichmented)

gplot1 <- function(DATA) {
  DATA$Enrichmented <- factor(DATA$Enrichmented,
                              levels = c("YES", 'NO'))
  temp <- NULL
  if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") > 10) {
    temp <- DATA[which(DATA$Enrichmented == "YES"),]
    temp <- temp[order(temp$Enrichment_log, decreasing = T),]
    # temp <- temp[1:10,]
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
    geom_hline(yintercept = 3) +
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

pdf("Figure5_2_IM_Peaks_Enrich_TF.pdf", width = 3.5, height = 4)
gplot1(IM_Peaks_Enrich)
dev.off()

pdf("Figure5_2_BM_Peaks_Enrich_TF.pdf", width = 3.5, height = 4)
gplot1(BM_Peaks_Enrich)
dev.off()

#### IM vs BM
{
  cells <- proj_pass_filter$cellNames
  cells <- cells[which(proj_pass_filter$TissuesSub %in% c("IM1cm","IM0.5cm"))]
  subsetArchRProject(ArchRProj = proj_pass_filter, cells = cells,
                     outputDirectory = "./snATAC-IM")
  IM_ATAC <- loadArchRProject(path = "./snATAC-IM")
  
  table(IM_ATAC$Celltype)
  IM_ATAC_IB <- IM_ATAC[which(IM_ATAC$Celltype %in% c("Branch meristems (BM)",
                                                      "Inflorescence meristem (IM)")),]
  table(IM_ATAC_IB$Celltype)
  
  IM_ATAC_IB <- addIterativeLSI(ArchRProj = IM_ATAC_IB, useMatrix = "TileMatrix",
                                name = "IterativeLSI", force = TRUE)
  IM_ATAC_IB <- addHarmony(
    ArchRProj = IM_ATAC_IB,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = T
  )
  
  IM_ATAC_IB <- addUMAP(
    ArchRProj = IM_ATAC_IB, 
    reducedDims = "IterativeLSI", 
    name = "UMAP_withBatchEffect", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
  )
  IM_ATAC_IB <- addUMAP(
    ArchRProj = IM_ATAC_IB, 
    reducedDims = "Harmony", 
    name = "UMAP_Harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
  )
  
  IM_ATAC_IB <- addTSNE(
    ArchRProj = IM_ATAC_IB, 
    reducedDims = "IterativeLSI", 
    name = "TSNE_withBatchEffect", 
    perplexity = 30,
    force = T
  )
  IM_ATAC_IB <- addTSNE(
    ArchRProj = IM_ATAC_IB, 
    reducedDims = "Harmony", 
    name = "TSNE_Harmony", 
    perplexity = 30,
    force = T
  )
  
  IM_ATAC_IB <- addMotifAnnotations(ArchRProj = IM_ATAC_IB, motifPWMs = motif_all,
                                    annoName = "TF-Motif", force = T)
  IM_ATAC_IB <- addBgdPeaks(IM_ATAC_IB)
  IM_ATAC_IB <- addDeviationsMatrix(
    ArchRProj = IM_ATAC_IB,
    peakAnnotation = "TF-Motif",
    force = TRUE,
    matrixName = "TF-deviation"
  )
  # IM_IB_TF_dev <- getMatrixFromProject(IM_ATAC_IB, useMatrix = "TF-deviation")
  # IM_IB_TF_dev_matrix <- IM_IB_TF_dev@assays@data@listData[["deviations"]]
  # 
}

########### RNA trajectory
unique(RNA$tissues)
IM_RNA <- subset(RNA, subset = tissues %in% c("IM1cm", "IM0.5cm"))
Meristem <- subset(IM_RNA, subset = Celltype %in% c("Inflorescence meristem (IM)", "Branch meristems (BM)"))
Meristem <- NormalizeData(Meristem)
Meristem <- FindVariableFeatures(Meristem)
Meristem <- ScaleData(Meristem)
Meristem <- RunPCA(Meristem)
Meristem <- RunHarmony(Meristem, reduction = "pca", dims.use = 1:30,
                       group.by.vars = "Batch")
Meristem <- RunUMAP(Meristem, dims = 1:30, reduction = "harmony")
Idents(Meristem) <- Meristem$Celltype
Meristem_markers <- FindAllMarkers(Meristem, only.pos = T)
Meristem_markers_filtered <- Meristem_markers[which(Meristem_markers$avg_log2FC > 0.5),]

DimPlot(Meristem, group.by = "Celltype", pt.size = 1,
        cols = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)")])
DimPlot(Meristem, group.by = "tissues")

library(monocle)
library(tidyverse)
library(patchwork)

data <- as.matrix(Meristem@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Meristem@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size()) # lowerDetectionLimit = 0.1

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 4)
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.2 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)
mycds <- reduceDimension(mycds, max_components = 2, reduction_method = 'DDRTree')
mycds <- orderCells(mycds)
plot_cell_trajectory(mycds, color_by = "State") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
mycds <- orderCells(mycds, root_state = "5")
plot_cell_trajectory(mycds, color_by = "Celltype") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
plot_cell_trajectory(mycds, color_by = "Pseudotime") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

pdf("Figure5_2_cell_trajectory", height = 5, width = 7.3)
plot_cell_trajectory(mycds, color_by = "Celltype") + 
  theme_bw() +
  scale_colour_manual(values = c("#faa818", "#5555ff")) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figure5_2_cell_trajectory_tissue", height = 5, width = 6.3)
plot_cell_trajectory(mycds, color_by = "tissues") +
  theme_bw() +
  labs(color = "Tissue stage") +
  scale_color_manual(values = c("#eea29a","#f7786b")) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figure5_2_Pseudotime.pdf", width = 7, height = 6)
plot_cell_trajectory(mycds, color_by = "Pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradientn(colours = c("#0000CD","#F08080","#FF4500"))
dev.off()

pdf("Figure5_2_State.pdf", width = 6, height = 5)
plot_cell_trajectory(mycds, color_by = "State") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pseudo_celltype <- data.frame(State = mycds$State,
                              Celltype = mycds$Celltype,
                              Pseudo = mycds$Pseudotime,
                              Cells = colnames(mycds))
rownames(pseudo_celltype) <- pseudo_celltype$Cells

Meristem$State <- pseudo_celltype[colnames(Meristem), "State"]

pdf("Figure5_2_UMAP_state.pdf", width = 5.3, height = 5)
DimPlot(Meristem, group.by = "State", label = T) +
  ggtitle("Trajectory State")
dev.off()

pdf("Figure5_2_UMAP_celltype.pdf", width = 7, height = 5)
DimPlot(Meristem, group.by = "Celltype") +
  ggtitle("Celltype")
dev.off()

table(Meristem$State)
State_Celltype <- as.data.frame(table(Meristem$State, Meristem$Celltype))
colnames(State_Celltype) <- c("State", "Celltype", "Count")
State_Celltype <- State_Celltype[order(State_Celltype$State),]
State_Celltype$Celltype_percent <- State_Celltype$Count / rep(c(868, 383, 398, 15, 2358, 207, 1299), each = 2)
State_Celltype$State <- factor(State_Celltype$State,
                               levels = c("5", "6", "4", "3", "2", "1", "7"))
pdf("Figure5_2_state_celltype_percent.pdf", height = 4.5, width = 6)
ggplot() +
  geom_bar(data = State_Celltype, aes(x = State, y = Celltype_percent, fill = Celltype),
           stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_color) +
  scale_y_continuous(labels = percent) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10), # angle = 45), # hjust = 1, vjust = 1),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  labs(fill = "Cell-type", x = "Trajectory state", y = "Percent")
dev.off()

#### state 1 vs 7
Idents(Meristem) <- Meristem$State
state_markers <- FindMarkers(Meristem, ident.1 = "1", ident.2 = "7",
                             logfc.threshold = -Inf)

state_markers$DataSet_symbol <- rownames(state_markers)
state_markers$Gene_id <- gene_id_map[match(state_markers$DataSet_symbol, gene_id_map$symbol),"gene_id"]
state_markers[which(is.na(state_markers$Gene_id)),]
state_markers <- as.data.frame(na.omit(state_markers))

state_1_row <- intersect(which(state_markers$p_val_adj < 0.05), which(state_markers$avg_log2FC > 0.5))
state_7_row <- intersect(which(state_markers$p_val_adj < 0.05), which(state_markers$avg_log2FC < -0.5))
state_1 <- state_markers[state_1_row,]
state_7 <- state_markers[state_7_row,]
state_1 <- state_1[order(state_1$avg_log2FC,decreasing = T),]
state_7 <- state_7[order(state_7$avg_log2FC,decreasing = F),]

color <- rep("gray", nrow(state_markers))
cex <- rep(1.2, nrow(state_markers))
cex[state_1_row] <- 2
cex[state_7_row] <- 2
cex <- as.character(cex)
color[state_1_row] <- "blue"
color[state_7_row] <- "red"
color <- factor(color, levels = c("blue", "red", "gray"))
# sigdiff <- as.data.frame(rbind(state_1[1:5,], state_7[1:5,]))
sigdiff <- state_markers[c("Os01g0140100", "LP1", "OsSPL18", "RBPA"),]
pdf("Figure5_2_state_1_vs_7_volcano.pdf", height = 3.5, width = 5)
vio_plot <- ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  labs(x = "log2 FoldChange", y = "-log10(adjusted P value)") +
  geom_point(data = state_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), size = cex, colour = color, alpha = "0.3")) +
  geom_point(data = sigdiff[c("RBPA", "Os01g0140100"),], aes(x = avg_log2FC, y = -log10(p_val_adj)),
             color = "black", shape = 21, fill = "#F05662", size = 2) +
  geom_point(data = sigdiff[c("OsSPL18", "LP1"),], aes(x = avg_log2FC, y = -log10(p_val_adj)),
             color = "black", shape = 21, fill = "#7990C8", size = 2) +
  geom_hline(aes(yintercept=-log10(0.05)), colour = "black", linetype="dashed") +
  geom_vline(aes(xintercept= -0.5), colour="black", linetype="dashed") +
  geom_vline(aes(xintercept= 0.5), colour="black", linetype="dashed") +
  scale_color_manual(name = "",
                     values = c("blue" = '#F05662', 'red' = '#7990C8', "gray" = "gray"),
                     labels = c('\nState 1\nN = 506', "\nState 7\nN = 564", "\nNot significant\n")) +
  scale_size_manual(values = c('1.2' = 1.2, "2" = 2)) +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  scale_alpha_manual(values = c("0.3" = 0.7, "1" = 1)) +
  guides(size = "none", alpha = "none")
vio_plot +
  geom_text_repel(data = sigdiff, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                      label = rownames(sigdiff),
                                      alpha = rep("1",nrow(sigdiff))),
                  size = 3.3,
                  box.padding = unit(0.8, "lines"),
                  point.padding = unit(0, "lines"), 
                  min.segment.length = 0,
                  segment.color = "black",
                  colour="#000000",
                  show.legend = FALSE,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30)) 
dev.off()


# GO
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
  
  genes_list <- list(state_1_genes = state_1$Gene_id,
                     state_7_genes = state_7$Gene_id)
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
      GO_result <- paste0("There is only ", length(genes), " DEGs of tissue ",i,", as a result, no GO pathways were enriched.")   
      write.table(GO_result, quote = F, row.names = F, col.names = F, sep = "\t",
                  paste0("./Figure5_2_",i,"_enrich_GO_result.tsv"))
      GO_result_Sig <- paste0("There is no GO pathways were significantly enriched.") 
      write.table(GO_result_Sig, quote = F, row.names = F, col.names = F, sep = "\t",
                  paste0("./Figure5_2_",i,"_enrich_GO_result_Sig.tsv"))
    } else {
      GO_result$FDR <- p.adjust(GO_result$P_value, method = "fdr")
      GO_result <- GO_result[order(GO_result$Ontology, GO_result$P_value),]
      GO_result_Sig <- GO_result[which(GO_result$P_value<0.05),]
      if (nrow(GO_result_Sig) < 1) {
        write.table(GO_result, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure5_2_",i,"_enrich_GO_result.tsv"))
        GO_result_Sig <- paste0("There is no GO pathways were significantly enriched.")
        write.table(GO_result_Sig, quote = F, row.names = F, col.names = F, sep = "\t",
                    paste0("./Figure5_2_",i,"_enrich_GO_result_Sig.tsv"))
      } else {
        write.table(GO_result, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure5_2_",i,"_enrich_GO_result.tsv"))
        write.table(GO_result_Sig, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure5_2_",i,"_enrich_GO_result_Sig.tsv"))
      }
    }
  }
  
}

# barplot for Sig-enrichment
Enrichment_data <- c()
# for (i in names(genes_list)) {
for (i in c("state_1_genes", "state_7_genes")) {
  print(i)
  temp_enrich <- read.table(paste0("./Figure5_2_",i,"_enrich_GO_result_Sig.tsv"), sep = "\t", header = T, row.names = 1)
  temp_enrich$Group <- unlist(strsplit(i, "_genes"))[1]
  Enrichment_data <- as.data.frame(rbind(Enrichment_data,
                                         temp_enrich[1:20,]))
}
library(Hmisc)
# Enrichment_data$Description2 <- paste(Enrichment_data$Description, Enrichment_data$Group, sep = "_")
Enrichment_data$Description2 <- Enrichment_data$Description
c<- capitalize(as.character(Enrichment_data$Description2))
Enrichment_data$Description <- capitalize(as.character(Enrichment_data$Description))
Enrichment_data$Description2 <- factor(Enrichment_data$Description2,
                                       levels = Enrichment_data$Description2,
                                       labels = Enrichment_data$Description2)
Enrichment_data$FDR_2 <- -log10(Enrichment_data$FDR)
pdf("Figure5_2_Branch_GOBP_enrich_Top20.pdf", width = 10, height = 6.5)
ggplot(Enrichment_data, aes(x = forcats::fct_reorder(Description2, FDR_2, .desc = T), y = FDR_2)) +
  geom_bar(aes(fill = Gene_count), stat = "identity") +
  # coord_flip() +
  theme_bw() +
  facet_wrap(Group ~ ., scales = "free", nrow = 1) +
  scale_fill_viridis(option = "C") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black", angle = 65, hjust = 1, vjust = 1),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  labs(color = "-Log10(FDR)")
dev.off()

##### TF list BEAM
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

TF_Family_2 <- TF_Family[which(TF_Family$symbol %in% rownames(mycds)),]
TF_Family_2$mean_exp <- disp_table[match(TF_Family_2$symbol, disp_table$gene_id),"mean_expression"]
TF_Family_2 <- as.data.frame(na.omit(TF_Family_2))
TF_Family_2 <- TF_Family_2[which(TF_Family_2$mean_exp > 0.1),]

BEAM_res <- BEAM(mycds[TF_Family_2$symbol,], branch_point = 3, cores = 1,
                 progenitor_method = "duplicate", verbose = TRUE)
saveRDS(BEAM_res, file = "Figure_5_2_TFs_BEAM_res.rds")
BEAM_res <- BEAM_res[, c("gene_short_name","pval","qval")]
BEAM_res <- as.data.frame(na.omit(BEAM_res))
BEAM_res_filter <- BEAM_res[which(BEAM_res$qval < 0.05),]
branched_heatmap_TFs <- plot_genes_branched_heatmap(mycds[BEAM_res_filter$gene_short_name,],
                                                branch_point = 3,
                                                num_clusters = 3, #这些基因被分成几个group
                                                cores = 1,
                                                branch_labels = c("Branch 1", "Branch 2"),
                                                hmcols = hmcols,
                                                branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                                use_gene_short_name = T,
                                                show_rownames = F,
                                                return_heatmap = T) #是否返回一些重要信息
pdf("Figure5_3_monocle_branch_heatmap.pdf", height = 5, width = 4.5)
branched_heatmap_TFs$ph_res
dev.off()

clusters <- cutree(branched_heatmap_TFs$ph_res$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$Genes <- rownames(clustering)
table(clustering$Gene_Clusters)
clustering$ID <- gene_id_map[match(clustering$Genes, gene_id_map$symbol), "gene_id"]
clustering$TF_family <- TF_Family[match(clustering$ID, TF_Family$RAP_ID), "TF_Family"]
table(clustering$Gene_Clusters, clustering$TF_family)

rownames(clustering) <- clustering$ID

# meristem/panicle/grain related genes
meristem_related_genes <- xlsx::read.xlsx("Figure5_meristem_genes.xlsx",
                                          sheetIndex = 1)
meristem_related_genes <- meristem_related_genes[which(meristem_related_genes$Locus.ID %in% clustering$ID),]
meristem_related_genes$heatmap_clusters <- clustering[match(meristem_related_genes$Locus.ID, clustering$ID), "Gene_Clusters"]
meristem_related_genes$Symbol <- clustering[meristem_related_genes$Locus.ID, "Genes"]
# panicle related genes
panicle_related_genes <- xlsx::read.xlsx("Figure5_panicle_genes.xlsx",
                                         sheetIndex = 1)
panicle_related_genes <- panicle_related_genes[which(panicle_related_genes$Locus.ID %in% clustering$ID),]
panicle_related_genes$heatmap_clusters <- clustering[match(panicle_related_genes$Locus.ID, clustering$ID), "Gene_Clusters"]
panicle_related_genes$Symbol <- clustering[panicle_related_genes$Locus.ID, "Genes"]
# grain related genes
grain_related_genes <- xlsx::read.xlsx("Figure5_grain_genes.xlsx",
                                       sheetIndex = 1)
grain_related_genes <- grain_related_genes[which(grain_related_genes$Locus.ID %in% clustering$ID),]
grain_related_genes$heatmap_clusters <- clustering[match(grain_related_genes$Locus.ID, clustering$ID), "Gene_Clusters"]
grain_related_genes$Symbol <- clustering[grain_related_genes$Locus.ID, "Genes"]

grain_related_genes_selected <- c("OsbZIP60", "OsSPL18")
panicle_related_genes_selected <- c("OsSPL18", "Os02g0284500")
meristem_related_genes_selected <- c("GRF10", "MADS17", "MADS18")
pdf("Figure5_2_psudotime_branch_expression.pdf", width = 4.5, height = 4)
plot_genes_branched_pseudotime(mycds[c("MADS18", "MADS17", "GRF10",
                                       "Os02g0284500", "OsSPL18", "OsbZIP60"),],
                               branch_point = 3,
                               color_by = "State",
                               cell_size = 0.5,
                               ncol = 2, min_expr = NULL) +
  scale_color_manual(values = state_color)
dev.off()

###############  非限制性整合
IM_ATAC_IB <- addGeneIntegrationMatrix(
  ArchRProj = IM_ATAC_IB,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_State_Un",
  reducedDims = "IterativeLSI",
  seRNA = Meristem,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "State",
  nameCell = "predictedCell_State_Un",
  nameGroup = "predictedGroup_State_Un",
  nameScore = "predictedScore_State_Un"
)


{ # Motif_all RAPDB symbol
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
  motif_all_RAPsymbol <- convert_motifs(motif, class = "TFBSTools-PWMatrix")
  motif_all_RAPsymbol <- do.call(PWMatrixList, motif_all_RAPsymbol)
  MSU_ID <- unlist(lapply(motif_all_RAPsymbol@listData, function(x){
    x@name
  }))
  names(motif_all_RAPsymbol) <- MSU_ID
  RAP_ID <- MSU_RAP[match(MSU_ID, MSU_RAP$MSU),]
  RAP_ID <- na.omit(RAP_ID)
  motif_all_RAPsymbol <- motif_all_RAPsymbol[RAP_ID$MSU]
  names(motif_all_RAPsymbol) <- RAP_ID$RAP
  motif_all_RAPsymbol <- motif_all_RAPsymbol[!duplicated(names(motif_all_RAPsymbol))] # ID duplicated !!!
  rownames(gene_id_map) <- gene_id_map$gene_id
  motif_all_RAPsymbol <- motif_all_RAPsymbol[names(motif_all_RAPsymbol) %in% gene_id_map$gene_id]
  Symbol <- match(names(motif_all_RAPsymbol), gene_id_map$gene_id)
  Symbol <- gene_id_map[Symbol,]
  names(motif_all_RAPsymbol) <- Symbol[names(motif_all_RAPsymbol), "symbol"]
  IM_ATAC_IB <- addMotifAnnotations(ArchRProj = IM_ATAC_IB, motifPWMs = motif_all_RAPsymbol,
                                    annoName = "TF-Motif-RAPsymbol", force = T)
  IM_ATAC_IB <- addDeviationsMatrix(
    ArchRProj = IM_ATAC_IB, 
    peakAnnotation = "TF-Motif-RAPsymbol",
    force = TRUE,
    matrixName = "MotifMatrix-RAPsymbol"
  )
  seGroupMotif <- getGroupSE(ArchRProj = IM_ATAC_IB, useMatrix = "MotifMatrix-RAPsymbol", groupBy = "Celltype")
  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
  }) %>% Reduce("cbind", .) %>% rowMaxs
  
  TF_dev_RAPsymbol <- getMatrixFromProject(IM_ATAC_IB, useMatrix = "MotifMatrix-RAPsymbol")
  TF_dev_RAPsymbol <- TF_dev_RAPsymbol@assays@data@listData[["deviations"]]
  TF_dev_RAPsymbol <- as.data.frame(TF_dev_RAPsymbol)
  
  corGIM_MM <- correlateMatrices(IM_ATAC_IB, useMatrix1 = "GeneIntegrationMatrix_State_Un",
                                 useMatrix2 = "MotifMatrix-RAPsymbol",
                                 log2Norm1 = TRUE, log2Norm2 = FALSE,
                                 removeFromName1 = NULL, removeFromName2 = NULL)
  corGIM_MM <- as.data.frame(corGIM_MM)
  corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$GeneIntegrationMatrix_State_Un_name, rowData(seZ)$name), "maxDelta"]
  corGIM_MM$TFRegulator <- "NO"
  corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0 & corGIM_MM$padj < 0.05 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.6))] <- "YES"
  sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
  corGIM_MM_YES <- corGIM_MM[which(corGIM_MM$TFRegulator == "YES"),]
  
  min(corGIM_MM_YES$cor) # 0.15
  
  rownames(Annotation_genes) <- Annotation_genes$RAP_DB_Symbol
  corGIM_MM_YES$Symbol <- Annotation_genes[corGIM_MM_YES$MotifMatrix.RAPsymbol_name, "CGSNL Gene Symbol"]
  p <- ggplot() +
    geom_point(data = data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black")) +
    labs(color = "Positive TF-Regulator") +
    geom_vline(xintercept = 0.15, lty = "dashed") + 
    geom_hline(yintercept = quantile(corGIM_MM$maxDelta, 0.6), lty = "dashed") + 
    scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
    xlab("Correlation to integrated gene expression") +
    ylab("Max delta TF ChromVAR\ndeviation scores between IM and BM") +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(0, max(corGSM_MM$maxDelta)*1.05)
    ) +
    geom_text_repel(data = data.frame(corGIM_MM_YES), aes(cor, maxDelta, label = Symbol),
                    size = 3.3,
                    box.padding = unit(0.4, "lines"),
                    point.padding = unit(0, "lines"), 
                    min.segment.length = 0,
                    segment.color = "black",
                    colour="#000000",
                    show.legend = FALSE,
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 100))
  pdf("Figure5_3_ATAC_positive_TF_regulator.pdf", width = 6, height = 5)
  p
  dev.off()
  
  IM_ATAC_IB <- addGroupCoverages(ArchRProj = IM_ATAC_IB,
                                  groupBy = "Celltype", force = T)
  motifPositions_IB <- getPositions(IM_ATAC_IB, name = "TF-Motif-RAPsymbol")
  seFoot_IB <- getFootprints(
    ArchRProj = IM_ATAC_IB, 
    positions = motifPositions_IB, 
    groupBy = "Celltype"
  )
  plotFootprints(
    seFoot = seFoot_IB,
    ArchRProj = IM_ATAC_IB, 
    normMethod = "Subtract",
    plotName = "IB_TF",
    addDOC = FALSE,
    smoothWindow = 5,
    pal = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)")]
  )
}

{
  markersPeaks_Celltype_IB <- getMarkerFeatures(
    ArchRProj = IM_ATAC_IB,
    useMatrix = "PeakMatrix",
    groupBy = "Celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  markerPeaksList_Celltype_IB <- getMarkers(markersPeaks_Celltype_IB, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  BM_marker_peaks <- markerPeaksList_Celltype_IB@listData[["Branch meristems (BM)"]]
  BM_marker_peaks <- as.data.frame(BM_marker_peaks)
  rownames(BM_marker_peaks) <- paste(BM_marker_peaks$seqnames,
                                     BM_marker_peaks$start,
                                     BM_marker_peaks$end,
                                     sep="_")
  BM_marker_peaks$ID <- rownames(BM_marker_peaks)
  
  IM_marker_peaks <- markerPeaksList_Celltype_IB@listData[["Inflorescence meristem (IM)"]]
  IM_marker_peaks <- as.data.frame(IM_marker_peaks)
  rownames(IM_marker_peaks) <- paste(IM_marker_peaks$seqnames,
                                       IM_marker_peaks$start,
                                       IM_marker_peaks$end,
                                       sep="_")
  IM_marker_peaks$ID <- rownames(IM_marker_peaks)
  
  peaksets <- getPeakSet(IM_ATAC_IB)
  pr1 <- paste(seqnames(peaksets),start(peaksets),end(peaksets),sep="_")
  peaksets <- as.data.frame(peaksets@elementMetadata)
  rownames(peaksets) <- pr1
  
  BM_marker_peaks$peakType <- peaksets[BM_marker_peaks$ID, "peakType"]
  table(BM_marker_peaks$peakType)
  # Distal   Exonic Intronic Promoter 
  # 961       88      160      491 
  IM_marker_peaks$peakType <- peaksets[IM_marker_peaks$ID, "peakType"]
  table(IM_marker_peaks$peakType)
  # Distal   Exonic Intronic Promoter 
  # 952      115      127      666
  
  BM_marker_peaks$Celltype <- "Branch meristems (BM)"
  IM_marker_peaks$Celltype <- "Inflorescence meristem (IM)"
  temp_marker_peaks <- as.data.frame(rbind(BM_marker_peaks,
                                           IM_marker_peaks))
  pdf("Figure5_3_IM_BM_marker_peaktype_percent.pdf", width = 3, height = 5.5)
  ggplot() +
    geom_bar(data = temp_marker_peaks, aes(x = Celltype, fill = peakType),
             stat = "count", position = "stack") +
    scale_fill_manual(values = peak_region_color) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    labs(fill = "Peak Region", x = "", y = "Marker peak count")
  dev.off()
  
  
  BM_marker_peaks_PD <- BM_marker_peaks[which(BM_marker_peaks$peakType %in% c("Promoter",
                                                                              "Distal")),]
  IM_marker_peaks_PD <- IM_marker_peaks[which(IM_marker_peaks$peakType %in% c("Promoter",
                                                                              "Distal")),]
  
  temp <- matrix(FALSE, nrow = length(c(IM_marker_peaks_PD$ID,
                                        BM_marker_peaks_PD$ID)), ncol = 2)
  rownames(temp) <- c(IM_marker_peaks_PD$ID,
                      BM_marker_peaks_PD$ID)
  colnames(temp) <- c("Branch meristems (BM)", "Inflorescence meristem (IM)")
  temp[rownames(BM_marker_peaks_PD), 1] <- TRUE
  temp[rownames(IM_marker_peaks_PD), 2] <- TRUE
  # IM_ATAC_IB <- addMotifAnnotations(ArchRProj = IM_ATAC_IB, motifPWMs = motif_all,
  #                                   annoName = "TF-Motif", force = T)
  matches <- getMatches(IM_ATAC_IB, name = "TF-Motif-RAPsymbol")
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  rownames(matches) <- pr1
  matches <- matches[rownames(temp),]
  identical(rownames(temp), rownames(matches))
  
  BM_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,1]), 1:nrow(matches))
  BM_Peaks_Enrich$Enrichment_log <- log2(BM_Peaks_Enrich$Enrichment)
  BM_Peaks_Enrich$Enrichmented <- "NO"
  BM_Peaks_Enrich$Enrichmented[which(BM_Peaks_Enrich$Enrichment_log >= 0.25 & BM_Peaks_Enrich$mlog10Padj >= -log10(0.05))] <- "YES"
  table(BM_Peaks_Enrich$Enrichmented)
  
  IM_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,2]), 1:nrow(matches))
  IM_Peaks_Enrich$Enrichment_log <- log2(IM_Peaks_Enrich$Enrichment)
  IM_Peaks_Enrich$Enrichmented <- "NO"
  IM_Peaks_Enrich$Enrichmented[which(IM_Peaks_Enrich$Enrichment_log >= 0.25 & IM_Peaks_Enrich$mlog10Padj >= -log10(0.05))] <- "YES"
  table(IM_Peaks_Enrich$Enrichmented)
  
  gplot1_2 <- function(DATA) {
    DATA$Enrichmented <- factor(DATA$Enrichmented,
                                levels = c("YES", 'NO'))
    temp <- DATA[which(DATA$Enrichmented == "YES"),]
    # if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") > 10) {
    #   temp <- DATA[which(DATA$Enrichmented == "YES"),]
    #   temp <- temp[order(temp$Enrichment_log, decreasing = T),]
    #   temp <- temp[1:10,]
    # }
    # if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") <= 10) {
    #   temp <- DATA[which(DATA$Enrichmented == "YES"),]
    #   temp <- temp[order(temp$Enrichment_log, decreasing = T),]
    # }
    
    g <- ggplot() +
      geom_point(data = DATA, aes(x = Enrichment_log, y = mlog10Padj, color = Enrichmented)) +
      theme_bw() +
      theme(axis.title = element_text(size = 10, colour = "black"),
            axis.text = element_text(size = 10, color = "black")) +
      labs(y = "-Log10(adjusted P-value)", x = "Log2(Enrichment)") +
      geom_hline(yintercept = -log10(0.05)) +
      geom_vline(xintercept = 0.25)
    if (is.null(temp)) {
      g <- g + 
        scale_color_manual(values = c("gray"))
    }
    if (!is.null(temp)) {
      g <- g +
        scale_color_manual(values = c("red", "gray")) + 
        geom_text_repel(data = temp, aes(x = Enrichment_log, y = mlog10Padj, 
                                         label = temp$feature2),
                        size = 3.3,
                        box.padding = unit(0.8, "lines"),
                        point.padding = unit(0, "lines"), 
                        min.segment.length = 0,
                        segment.color = "black",
                        colour="#000000",
                        show.legend = FALSE,
                        max.overlaps = getOption("ggrepel.max.overlaps", default = 50))
    }
    
    if (max(DATA$Enrichment_log) < 0.25) {
      g <- g +
        scale_x_continuous(limits = c(min(DATA$Enrichment_log), 0.5))
    }
    g <- g + theme(legend.position = "none")
    return(g)
  }
  
  BM_Peaks_Enrich$feature2 <- Annotation_genes[BM_Peaks_Enrich$feature, "CGSNL Gene Symbol"]
  IM_Peaks_Enrich$feature2 <- Annotation_genes[IM_Peaks_Enrich$feature, "CGSNL Gene Symbol"]
  pdf("Figure5_3_IM_Peaks_Enrich_TF(IM_vs_BM).pdf", width = 3.5, height = 4)
  gplot1_2(IM_Peaks_Enrich)
  dev.off()
  
  pdf("Figure5_3_BM_Peaks_Enrich_TF(IM_vs_BM).pdf", width = 3.5, height = 4)
  gplot1_2(BM_Peaks_Enrich)
  dev.off()
  
  BM_Peaks_Enrich_TF <- BM_Peaks_Enrich[which(BM_Peaks_Enrich$Enrichmented == "YES"), "feature"]
  IM_Peaks_Enrich_TF <- IM_Peaks_Enrich[which(IM_Peaks_Enrich$Enrichmented == "YES"), "feature"]
  
  BM_Peak_TF_positive <- intersect(BM_Peaks_Enrich_TF, corGIM_MM_YES$MotifMatrix.RAPsymbol_name)
  IM_Peak_TF_positive <- intersect(IM_Peaks_Enrich_TF, corGIM_MM_YES$MotifMatrix.RAPsymbol_name)
  
  intersect(c(BM_Peak_TF_positive,
              IM_Peak_TF_positive),
            BEAM_res_filter$gene_short_name)
  
  library(ggVennDiagram)
  pdf("Figure5_3_TF_venn.pdf")
  ggVennDiagram(list("Positive TFs" = corGIM_MM_YES$MotifMatrix.RAPsymbol_name,
                     "BM marker peaks enriched TFs" = BM_Peaks_Enrich_TF,
                     "IM marker peaks enriched TFs" = IM_Peaks_Enrich_TF),
                category.names = c("Positive TFs", "BM marker peaks enriched TFs",
                                   "IM marker peaks enriched TFs"),
                show_intersect = FALSE,
                set_color = "black",
                set_size = NA,
                label = c("both", "count", "percent", "none"),
                label_alpha = 0.5,
                label_geom = c("label", "text"),
                label_color = "black",
                label_size = NA,
                label_percent_digit = 0,
                label_txtWidth = 40,
                edge_lty = "solid",
                edge_size = 1)
  dev.off()
  
  Annotation_genes[BM_Peak_TF_positive,]
  Annotation_genes[IM_Peak_TF_positive,]
  sum(corGIM_MM_YES$MotifMatrix.RAPsymbol_name %in% rownames(mycds))
  mycds_2 <- mycds[corGIM_MM_YES$MotifMatrix.RAPsymbol_name,]
  
  rownames(Annotation_genes) <- Annotation_genes$DataSets_Symbol
  mycds_2@featureData@data[["gene_short_name"]] <- Annotation_genes[mycds_2@featureData@data[["gene_short_name"]], "CGSNL Gene Symbol"]
  rownames(mycds_2) <- Annotation_genes[rownames(mycds_2), "CGSNL Gene Symbol"]
  
  heatmap_matrix <- plot_pseudotime_heatmap(mycds_3,
                                            num_clusters = 3,
                                            cores = 2,
                                            hmcols = hmcols,
                                            show_rownames = T, return_heatmap = T,
                                            use_gene_short_name = T)
  dev.off()
  
  pdf("Figure5_3_monocle_heatmap_Peak_positive.pdf", height = 13, width = 4.5)
  mycds_3 <- mycds_2[corGIM_MM_YES$Symbol,]
  ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
                     cluster_rows = T, show_rownames = T, 
                     show_colnames = F, #clustering_distance_rows = row_dist, 
                     clustering_method = "ward.D",
                     cutree_rows = 3, 
                     # annotation_row = annotation_row,
                     # annotation_col = annotation_col, 
                     treeheight_row = 20, # breaks = bks,
                     fontsize = 6, color = hmcols, 
                     border_color = NA, silent = TRUE, filename = NA)
  ph_res
  dev.off()
  
  
  
  pdf("Figure5_3_pseudotime_celltype_density.pdf")
  ggplot(data = mycds@phenoData@data, aes(x = Pseudotime, group = Celltype, color = Celltype)) +
    geom_density() +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black")) + 
    scale_color_manual(values = c("Branch meristems (BM)" = "#faa818",
                                  "Inflorescence meristem (IM)" = "#5555ff")) +
    labs(x = "Pseudotime", y = "Density")
  dev.off()
  
  branched_heatmap_Peak_positive <- plot_pseudotime_heatmap(mycds[c(BM_Peak_TF_positive,
                                                              IM_Peak_TF_positive),],
                                                      branch_point = 3,
                                                      num_clusters = 3, #这些基因被分成几个group
                                                      cores = 1,
                                                      branch_labels = c("Branch 1", "Branch 2"),
                                                      hmcols = hmcols,
                                                      branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                                      use_gene_short_name = T,
                                                      show_rownames = T,
                                                      return_heatmap = T) #是否返回一些重要信息
  pdf("Figure5_3_monocle_branch_heatmap_Peak_positive.pdf", height = 5, width = 4.5)
  branched_heatmap_Peak_positive$ph_res
  dev.off()
  
  enrich_TF_heatmap <- data.frame(row.names = BM_Peaks_Enrich$feature,
                                  tf = BM_Peaks_Enrich$feature,
                                  BM_Enrichment_log = BM_Peaks_Enrich$Enrichment_log,
                                  IM_Enrichment_log = IM_Peaks_Enrich$Enrichment_log,
                                  BM_adj = BM_Peaks_Enrich$mlog10Padj,
                                  IM_adj = IM_Peaks_Enrich$mlog10Padj)
  ggplot(data = enrich_TF_heatmap, aes(x = BM_Enrichment_log,
                                       y = IM_Enrichment_log)) +
    geom_point()
  
  
}

{
  table(IM_ATAC_IB$predictedGroup_State_Un, IM_ATAC_IB$Celltype)
  table(IM_ATAC_IB$predictedGroup_State_Un)
  
  ATAC_State_Celltype <- as.data.frame(table(IM_ATAC_IB$predictedGroup_State_Un, IM_ATAC_IB$Celltype))
  colnames(ATAC_State_Celltype) <- c("Predicted_State", "Celltype", "Count")
  ATAC_State_Celltype <- ATAC_State_Celltype[order(ATAC_State_Celltype$Predicted_State),]
  ATAC_State_Celltype$Celltype_percent <- ATAC_State_Celltype$Count / rep(c(344, 30, 242, 3739, 313, 1885), each = 2)
  ATAC_State_Celltype$Predicted_State <- factor(ATAC_State_Celltype$Predicted_State,
                                                levels = c("5", "6", "3", "2", "1", "7"))
  pdf("Figure5_3_ATAC_predicted_state_celltype_percent.pdf", height = 4.5, width = 5.5)
  ggplot() +
    geom_bar(data = ATAC_State_Celltype, aes(x = Predicted_State, y = Celltype_percent, fill = Celltype),
             stat = "identity", position = "stack") +
    scale_fill_manual(values = celltype_color) +
    scale_y_continuous(labels = percent) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10), # angle = 45), # hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    labs(fill = "Cell-type", x = "Trajectory state", y = "Percent")
  dev.off()
  
  UMAP_Harmony <- getEmbedding(IM_ATAC_IB, embedding = "UMAP_Harmony")
  colnames(UMAP_Harmony) <- c("UMAP_1", "UMAP_2")
  ColData <- as.data.frame(getCellColData(IM_ATAC_IB))
  UMAP_Harmony <- as.data.frame(cbind(UMAP_Harmony,
                                      ColData))
  ggplot(data = UMAP_Harmony, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_State_Un)) +
    geom_point(size = 1) +
    scale_color_manual(values = state_color) + theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  temp <- as.data.frame(table(Meristem$State))
  colnames(temp) <- c("State", "RNA_num")
  table(IM_ATAC_IB$predictedGroup_State_Un)
  temp$ATAC_num_un <- c(344, 30, 242, 0, 3739, 313, 1885)
  temp$RNA_percent <- temp$RNA_num / sum(temp$RNA_num)
  temp$ATAC_percent_un <- temp$ATAC_num_un / sum(temp$ATAC_num_un)
  pdf("Figure5_3_state_percent_cor.pdf", width = 5.5, height = 5)
  ggplot(data = temp, aes(x = RNA_percent, y = ATAC_percent_un)) +
    geom_point(data = temp, aes(x = RNA_percent, y = ATAC_percent_un, color = State),
               size = 3) +
    geom_smooth(method = "lm", colour = "black") +
    ggpubr::stat_cor() +
    scale_color_manual(values = state_color) +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    labs(x = "Trajectory state percent", y = "Integration predicted state percent")
  dev.off()
  
  IM_ATAC_IB_GeneIntegrationMatrix <- getMatrixFromProject(IM_ATAC_IB,
                                                           useMatrix = "GeneIntegrationMatrix_State_Un")
  IM_ATAC_IB_GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix_State_Un"]]@Dimnames[[1]] <- IM_ATAC_IB_GeneIntegrationMatrix@elementMetadata@listData[["name"]]
  IM_ATAC_IB_GeneIntegrationMatrix <- IM_ATAC_IB_GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix_State_Un"]]
  
  IM_ATAC_IB_GeneIntegration <- CreateSeuratObject(IM_ATAC_IB_GeneIntegrationMatrix)
  IM_ATAC_IB_GeneIntegration$DataType <- "snATAC"
  IM_ATAC_IB_meta <- getCellColData(IM_ATAC_IB)
  IM_ATAC_IB_meta <- as.data.frame(IM_ATAC_IB_meta)
  IM_ATAC_IB_GeneIntegration$Celltype <- IM_ATAC_IB_meta[colnames(IM_ATAC_IB_GeneIntegration), "Celltype"]
  table(IM_ATAC_IB_GeneIntegration$Celltype)
  IM_ATAC_IB_GeneIntegration$DataType_Celltype <- paste0(IM_ATAC_IB_GeneIntegration$DataType,
                                                         "_",
                                                         IM_ATAC_IB_GeneIntegration$Celltype)
  IM_ATAC_IB_GeneIntegration <- NormalizeData(IM_ATAC_IB_GeneIntegration)
  Meristem$DataType <- "snRNA"
  Meristem$DataType_Celltype <- paste0(Meristem$DataType,
                                       "_",
                                       Meristem$Celltype)
  Meristem_ATAC_RNA_merge <- merge(IM_ATAC_IB_GeneIntegration,
                                   Meristem)
  Meristem_ATAC_RNA_merge <- NormalizeData(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- FindVariableFeatures(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- ScaleData(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- RunPCA(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- RunUMAP(Meristem_ATAC_RNA_merge,
                                     dims = 1:50)
  DimPlot(Meristem_ATAC_RNA_merge,
          group.by = "DataType")
  DimPlot(Meristem_ATAC_RNA_merge, 
          group.by = "Celltype")
  
  features <- SelectIntegrationFeatures(object.list = list(IM_ATAC_IB_GeneIntegration,
                                                           Meristem))
  anchors <- FindIntegrationAnchors(object.list = list(IM_ATAC_IB_GeneIntegration,
                                                       Meristem),
                                    anchor.features = features)
  Meristem_ATAC_RNA_integrate <- IntegrateData(anchorset = anchors)
  DefaultAssay(Meristem_ATAC_RNA_integrate) <- "integrated"
  Meristem_ATAC_RNA_integrate <- ScaleData(Meristem_ATAC_RNA_integrate, verbose = FALSE)
  Meristem_ATAC_RNA_integrate <- RunPCA(Meristem_ATAC_RNA_integrate, npcs = 50, verbose = FALSE)
  Meristem_ATAC_RNA_integrate <- RunUMAP(Meristem_ATAC_RNA_integrate, reduction = "pca", dims = 1:30)
  pdf("Figure5_3_Meristem_ATAC_RNA_integrate_celltype_UMAP.pdf", width = 6.7, height = 5)
  DimPlot(Meristem_ATAC_RNA_integrate,
          group.by = "Celltype", cols = celltype_color)
  dev.off()
  pdf("Figure5_3_Meristem_ATAC_RNA_integrate_datatype_UMAP.pdf", width = 5.5, height = 5)
  DimPlot(Meristem_ATAC_RNA_integrate,
          group.by = "DataType")
  dev.off()
  
  Meristem_ATAC_RNA_integrate@meta.data[IM_ATAC_IB$cellNames, "State"] <- IM_ATAC_IB$predictedGroup_State_Un
  table(Meristem_ATAC_RNA_integrate$DataType, Meristem_ATAC_RNA_integrate$State)
  pdf("Figure5_3_Meristem_ATAC_RNA_integrate_state_UMAP.pdf", width = 5, height = 5)
  DimPlot(Meristem_ATAC_RNA_integrate,
          group.by = "State", cols = state_color)
  dev.off()
  
  library(ggplot2)
  library(ggalluvial)
  IM_ATAC_IB$predictedCelltype_State_Un <- Meristem@meta.data[IM_ATAC_IB$predictedCell_State_Un, "Celltype"]
  table(IM_ATAC_IB$Celltype, IM_ATAC_IB$predictedCelltype_State_Un)
  (2554 + 3195) / (2554 + 3195 + 194 + 610)
  temp <- data.frame(ATAC = IM_ATAC_IB$Celltype, Predicted = IM_ATAC_IB$predictedCelltype_State_Un)
  df <- to_lodes_form(temp[,1:ncol(temp)],
                      axes = 1:ncol(temp),
                      id = "value")
  pdf("Figure5_3_ATAC_predicted_sankey.pdf", width = 3, height = 4.5)
  ggplot(df, aes(x = x, fill = stratum, label = stratum,
                 stratum = stratum, alluvium  = value)) + #数据
    geom_flow(width = 0.4,#连线宽度
              curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
              alpha = 0.5,#透明度
              color = 'white',#间隔颜色
              size = 0.1) + #间隔宽度
    geom_stratum(width = 0.35) + #图中方块的宽度
    geom_text(stat = 'stratum', size = 3.3, color = 'black') +
    scale_fill_manual(values = celltype_color[unique(IM_ATAC_IB$Celltype)]) + #自定义颜色
    theme_void() + #主题（无轴及网格线）
    theme(legend.position = 'none')
  dev.off()
}



# IM_Peaks_Enrich_2_TF_temp <- Annotation_genes[match(IM_Peaks_Enrich_2_TF,
#                                                     Annotation_genes$`CGSNL Gene Symbol`),]
# IM_Peaks_Enrich_2_TF_temp <- IM_Peaks_Enrich_2_TF_temp[which(IM_Peaks_Enrich_2_TF_temp$RAP_DB_Symbol %in% clustering$Genes),]
# 
# BM_Peaks_Enrich_2_TF_temp <- Annotation_genes[match(BM_Peaks_Enrich_2_TF,
#                                                     Annotation_genes$`CGSNL Gene Symbol`),]
# BM_Peaks_Enrich_2_TF_temp <- BM_Peaks_Enrich_2_TF_temp[which(BM_Peaks_Enrich_2_TF_temp$RAP_DB_Symbol %in% clustering$Genes),]
# 
# branched_heatmap_temp <- plot_genes_branched_heatmap(mycds[c(IM_Peaks_Enrich_2_TF_temp$RAP_DB_Symbol,
#                                                              BM_Peaks_Enrich_2_TF_temp$RAP_DB_Symbol),],
#                                                      branch_point = 3,
#                                                      num_clusters = 3, #这些基因被分成几个group
#                                                      cores = 1,
#                                                      branch_labels = c("Branch 1", "Branch 2"),
#                                                      hmcols = hmcols,
#                                                      branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
#                                                      use_gene_short_name = T,
#                                                      show_rownames = T,
#                                                      return_heatmap = T) #是否返回一些重要信息
# pdf("Figure5_2_monocle_branch_heatmap_IM_BM_enriched_TFs.pdf", height = 5, width = 4.5)
# branched_heatmap_temp$ph_res
# dev.off()





# markers
{
  # gene_id_map[which(gene_id_map$gene_id == "Os03g0753100"),]
  # FeaturePlot(Meristem, features = "PAP2", order = T)
  # DotPlot(Meristem, features = "PAP2", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os01g0848400"),]
  # FeaturePlot(Meristem, features = "qSH1", order = T)
  # DotPlot(Meristem, features = "qSH1", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os02g0743400"),]
  # FeaturePlot(Meristem, features = "PIN1A", order = T)
  # DotPlot(Meristem, features = "PIN1A", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os04g0411400"),]
  # FeaturePlot(Meristem, features = "Rcn4", order = T)
  # DotPlot(Meristem, features = "Rcn4", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os05g0455200"),]
  # DotPlot(Meristem, features = "Os05g0455200", group.by = "State")
  # DotPlot(Meristem, features = "Os05g0455200", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os06g0184500"),]
  # DotPlot(Meristem, features = "Os06g0184500", group.by = "State")
  # VlnPlot(Meristem, features = "Os06g0184500", group.by = "State")
  # DotPlot(Meristem, features = "Os06g0184500", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os06g0232300"),]
  # DotPlot(Meristem, features = "PIN1C", group.by = "State")
  # VlnPlot(Meristem, features = "PIN1C", group.by = "State")
  # DotPlot(Meristem, features = "PIN1C", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os06g0597500"),]
  # DotPlot(Meristem, features = "SPED1", group.by = "State")
  # VlnPlot(Meristem, features = "SPED1", group.by = "State")
  # DotPlot(Meristem, features = "SPED1", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os06g0665400"),]
  # DotPlot(Meristem, features = "APO1", group.by = "State")
  # VlnPlot(Meristem, features = "APO1", group.by = "State")
  # DotPlot(Meristem, features = "APO1", group.by = "Celltype")
  # 
  # gene_id_map[which(gene_id_map$gene_id == "Os07g0108900"),]
  # DotPlot(Meristem, features = "DEP", group.by = "State")
  # VlnPlot(Meristem, features = "DEP", group.by = "State")
  # DotPlot(Meristem, features = "DEP", group.by = "Celltype")
  # 
  # #
  # gene_id_map[which(gene_id_map$gene_id == "Os10g0478000"),]
  # DotPlot(Meristem, features = "G1L5", group.by = "State")
  # VlnPlot(Meristem, features = "G1L5", group.by = "State")
  # DotPlot(Meristem, features = "G1L5", group.by = "Celltype")
  # 
  # #
  # gene_id_map[which(gene_id_map$gene_id == "Os11g0152500"),]
  # DotPlot(Meristem, features = "RCN1", group.by = "State")
  # VlnPlot(Meristem, features = "RCN1", group.by = "State")
  # DotPlot(Meristem, features = "RCN1", group.by = "Celltype")
  # 
  # # 
  # gene_id_map[which(gene_id_map$gene_id == "Os12g0614600"),]
  # DotPlot(Meristem, features = "bif2", group.by = "State")
  # VlnPlot(Meristem, features = "bif2", group.by = "State")
  # DotPlot(Meristem, features = "bif2", group.by = "Celltype")
  # YES
  gene_id_map[which(gene_id_map$gene_id == "Os03g0752800"),]
  FeaturePlot(Meristem, features = "MADS14", order = T)
  DotPlot(Meristem, features = "MADS14", group.by = "Celltype")
  
  # YES
  gene_id_map[which(gene_id_map$gene_id == "Os07g0235800"),]
  DotPlot(Meristem, features = "SNB", group.by = "State")
  VlnPlot(Meristem, features = "SNB", group.by = "State")
  DotPlot(Meristem, features = "SNB", group.by = "Celltype")
  
  #  YES
  gene_id_map[which(gene_id_map$gene_id == "Os07g0605200"),]
  DotPlot(Meristem, features = "MADS18", group.by = "State")
  VlnPlot(Meristem, features = "MADS18", group.by = "State")
  DotPlot(Meristem, features = "MADS18", group.by = "Celltype")
  
  # YES
  gene_id_map[which(gene_id_map$gene_id == "Os02g0274100"),]
  DotPlot(Meristem, features = "AIM1", group.by = "State")
  VlnPlot(Meristem, features = "AIM1", group.by = "State")
  DotPlot(Meristem, features = "AIM1", group.by = "Celltype")
  
  gene_id_map[which(gene_id_map$gene_id == "Os01g0752200"),]
  DotPlot(Meristem, features = "NOG1", group.by = "State")
  VlnPlot(Meristem, features = "NOG1", group.by = "State")
  DotPlot(Meristem, features = "NOG1", group.by = "Celltype")
  
  mycds$PAP2 <- Meristem@assays$RNA@data["PAP2",colnames(mycds)]
  mycds$MADS14 <- Meristem@assays$RNA@data["MADS14",colnames(mycds)]
  mycds$SNB <- Meristem@assays$RNA@data["SNB",colnames(mycds)]
  mycds$MADS18 <- Meristem@assays$RNA@data["MADS18",colnames(mycds)]
  mycds$AIM1 <- Meristem@assays$RNA@data["AIM1",colnames(mycds)]
  mycds$OsIDD12 <- Meristem@assays$RNA@data["OsIDD12",colnames(mycds)]
  
  plot_cell_trajectory(mycds, color_by = "OsIDD12") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "gray", high = "red")
  
  DotPlot(Meristem, features = "Os03g0843700", group.by = "State")
}



plot_cell_trajectory(mycds, color_by = "MADS18") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "gray", high = "red")

rownames(Annotation_genes) <- Annotation_genes$RAP_DB_Symbol
sum(Annotation_genes$RAP_DB_Symbol %in% rownames(mycds))
temp <- Annotation_genes[which(Annotation_genes$RAP_DB_Symbol %in% rownames(mycds)),]
BEAM_res <- BEAM(mycds[temp$RAP_DB_Symbol,], branch_point = 3, cores = 1,
                 progenitor_method = "duplicate", verbose = TRUE)
# saveRDS(BEAM_res, file = "Figure_5_BEAM_res.rds")
BEAM_res <- BEAM_res[, c("gene_short_name","pval","qval")]
BEAM_res$mean_expression <- disp_table[match(BEAM_res$gene_short_name, disp_table$gene_id),"mean_expression"]
BEAM_res <- as.data.frame(na.omit(BEAM_res))
BEAM_res_filter <- BEAM_res[which(BEAM_res$qval < 0.05),]
BEAM_res_filter$Symbol <- Annotation_genes[match(BEAM_res_filter$gene_short_name, Annotation_genes$RAP_DB_Symbol),"CGSNL Gene Symbol"]
rownames(BEAM_res_filter) <- BEAM_res_filter$Symbol
# BEAM_res_filter_2 <- BEAM_res[which(BEAM_res$pval < 0.00001 & BEAM_res$mean_expression > 0.1),]

plotCol <- function(col, nrow=1, ncol=ceiling(length(col) / nrow), txt.col="black") {
  stopifnot(nrow >= 1, ncol >= 1)
  if(length(col) > nrow*ncol)
    warning("some colors will not be shown")
  grid.newpage()
  gl <- grid.layout(nrow, ncol)
  pushViewport(viewport(layout=gl))
  ic <- 1
  for(i in 1:nrow) {
    for(j in 1:ncol) {
      pushViewport(viewport(layout.pos.row=i, layout.pos.col=j))
      grid.rect(gp= gpar(fill=col[ic]))
      grid.text(col[ic], gp=gpar(col=txt.col))
      upViewport()
      ic <- ic+1
    }
  }
  upViewport()
  invisible(gl)
}
plotCol(colorRampPalette(colors = c("#0D25B9", "white", "#FD6585"))(62))
hmcols <- colorRampPalette(colors = c("#0000CD", "white", "#B22222"))(62)
# hmcols <- viridisLite::viridis(100, alpha = 1, begin = 0, end = 1, direction = 1, option = "B")
# hmcols <- gsub('.{2}$', '', hmcols)

TFs_temp <- c(BM_Peaks_Enrich_2_TF, IM_Peaks_Enrich_2_TF)
TFs_temp_RAP <- Annotation_genes[match(TFs_temp, Annotation_genes$`CGSNL Gene Symbol`), "RAP_DB_Symbol"]
sum(TFs_temp_RAP %in% BEAM_res_filter$gene_short_name)

rownames(temp) <- temp$RAP_DB_Symbol
mycds_temp <- mycds[temp$RAP_DB_Symbol,]
rownames(mycds_temp) <- temp$`CGSNL Gene Symbol`

BEAM_res_filter_2 <- BEAM_res_filter[c(BM_Peaks_Enrich_2_TF,
                                       IM_Peaks_Enrich_2_TF),]
BEAM_res_filter_2 <- as.data.frame(na.omit(BEAM_res_filter_2))

branched_heatmap <- plot_genes_branched_heatmap(mycds_temp[BEAM_res_filter_2$Symbol,],
                                                branch_point = 3,
                                                num_clusters = 3, #这些基因被分成几个group
                                                cores = 1,
                                                branch_labels = c("Cell fate 1", "Cell fate 2"),
                                                hmcols = hmcols,
                                                branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                                use_gene_short_name = T,
                                                show_rownames = T,
                                                return_heatmap = T) #是否返回一些重要信息
pdf("Figure5_2_monocle_branch_heatmap.pdf", height = 5, width = 4.5)
branched_heatmap$ph_res
dev.off()

clusters <- cutree(branched_heatmap$ph_res$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$Genes <- rownames(clustering)
table(clustering$Gene_Clusters)
# Annotation_genes <- Annotation_genes[Annotation_genes$DataSets_Symbol %in% rownames(mycds),]
heatmap_TF <- Annotation_genes[Annotation_genes$DataSets_Symbol %in% BEAM_res_filter_2$gene_short_name,]
Annotation_genes[Annotation_genes$DataSets_Symbol %in% BEAM_res_filter_2$gene_short_name,]
clustering[heatmap_TF$DataSets_Symbol,]



mycds$OsSPL10 <- Meristem@assays$RNA@data["OsSPL10",colnames(mycds)]
mycds$DOF25 <- Meristem@assays$RNA@data["DOF25",colnames(mycds)]
mycds$BZIP50 <- Meristem@assays$RNA@data["BZIP50",colnames(mycds)]
mycds$`OsMYB3R-2` <- Meristem@assays$RNA@data["OsMYB3R-2",colnames(mycds)]
mycds$GW8 <- Meristem@assays$RNA@data["GW8",colnames(mycds)]
mycds$DOF25 <- Meristem@assays$RNA@data["DOF25",colnames(mycds)]
mycds$`Roc5(t)` <- Meristem@assays$RNA@data["Roc5(t)",colnames(mycds)]
mycds$OsMADS2 <- Meristem@assays$RNA@data["OsMADS2",colnames(mycds)]
mycds$OsFBH1 <- Meristem@assays$RNA@data["OsFBH1",colnames(mycds)]
mycds$ARF8 <- Meristem@assays$RNA@data["ARF8",colnames(mycds)]
mycds$OsPHL3 <- Meristem@assays$RNA@data["OsPHL3",colnames(mycds)]
mycds$OsFBH1 <- Meristem@assays$RNA@data["OsFBH1",colnames(mycds)]
plot_cell_trajectory(mycds, color_by = "OsFBH1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "gray", high = "red")

##### TF deviation
{
  ##########
  ##### 整合
  ##########
  IM_ATAC_IB <- addGeneIntegrationMatrix(
      ArchRProj = IM_ATAC_IB,
      useMatrix = "GeneScoreMatrix",
      matrixName = "GeneIntegrationMatrix_State_Un",
      reducedDims = "IterativeLSI",
      seRNA = Meristem,
      addToArrow = TRUE,
      force= TRUE,
      groupRNA = "State",
      nameCell = "predictedCell_State_Un",
      nameGroup = "predictedGroup_State_Un",
      nameScore = "predictedScore_State_Un"
    )
  
  IM_ATAC_IB <- addGeneIntegrationMatrix(
    ArchRProj = IM_ATAC_IB,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix_Celltype_Un",
    reducedDims = "IterativeLSI",
    seRNA = Meristem,
    addToArrow = TRUE,
    force= TRUE,
    groupRNA = "Celltype",
    nameCell = "predictedCell_Celltype_Un",
    nameGroup = "predictedGroup_Celltype_Un",
    nameScore = "predictedScore_Celltype_Un"
  )
  
  getAvailableMatrices(IM_ATAC_IB)

  table(IM_ATAC_IB$predictedGroup_State_Un)
  table(IM_ATAC_IB$Celltype, IM_ATAC_IB$predictedGroup_Celltype_Un)
  
  UMAP_Harmony <- getEmbedding(IM_ATAC_IB, embedding = "UMAP_Harmony")
  colnames(UMAP_Harmony) <- c("UMAP_1", "UMAP_2")
  ColData <- as.data.frame(getCellColData(IM_ATAC_IB))
  UMAP_Harmony <- as.data.frame(cbind(UMAP_Harmony,
                                      ColData))
  pdf("Figure5_2_ATAC_UMAP_Celltype.pdf", height = 5, width = 6.8)
  ggplot(data = UMAP_Harmony, aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
    geom_point(size = 1) +
    scale_color_manual(values = celltype_color) + theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  dev.off()
  
  state_color <- c("#f8766d", "#c49a00", "#53b400", "#00c094", "#00b6eb", "#a58aff", "#fb61b7")
  names(state_color) <- 1:7
  ggplot(data = UMAP_Harmony, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_State_Un)) +
    geom_point(size = 1) +
    scale_color_manual(values = state_color) + theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(.~predictedGroup_State_Un)
  
  ### 限制性整合
  groupList <- SimpleList(
    "Branch meristems (BM)" = SimpleList(
      ATAC = IM_ATAC_IB$cellNames[IM_ATAC_IB$Celltype %in% "Branch meristems (BM)"],
      RNA = colnames(Meristem)[Meristem$Celltype %in% "Branch meristems (BM)"]
    ),
    "Inflorescence meristem (IM)" = SimpleList(
      ATAC = IM_ATAC_IB$cellNames[IM_ATAC_IB$Celltype %in% "Inflorescence meristem (IM)"],
      RNA = colnames(Meristem)[Meristem$Celltype %in% "Inflorescence meristem (IM)"]
    )    
  )
  
  IM_ATAC_IB <- addGeneIntegrationMatrix(
    ArchRProj = IM_ATAC_IB,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix_State_Co",
    reducedDims = "IterativeLSI",
    seRNA = Meristem,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "State",
    nameCell = "predictedCell_State_Co",
    nameGroup = "predictedGroup_State_Co",
    nameScore = "predictedScore_State_Co"
  )
  
  IM_ATAC_IB <- addGeneIntegrationMatrix(
    ArchRProj = IM_ATAC_IB,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix_Celltype_Co",
    reducedDims = "IterativeLSI",
    seRNA = Meristem,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "Celltype",
    nameCell = "predictedCell_Celltype_Co",
    nameGroup = "predictedGroup_Celltype_Co",
    nameScore = "predictedScore_Celltype_Co"
  )
  
  table(IM_ATAC_IB$predictedGroup_Celltype_Co, IM_ATAC_IB$Celltype)
  table(IM_ATAC_IB$predictedGroup_State_Co)
  
  UMAP_Harmony <- getEmbedding(IM_ATAC_IB, embedding = "UMAP_Harmony")
  colnames(UMAP_Harmony) <- c("UMAP_1", "UMAP_2")
  ColData <- as.data.frame(getCellColData(IM_ATAC_IB))
  UMAP_Harmony <- as.data.frame(cbind(UMAP_Harmony,
                                      ColData))
  ggplot(data = UMAP_Harmony, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_State_Co)) +
    geom_point(size = 1) +
    scale_color_manual(values = state_color) + theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  temp <- as.data.frame(table(Meristem$State))
  colnames(temp) <- c("State", "RNA_num")
  table(IM_ATAC_IB$predictedGroup_State_Co)
  temp$ATAC_num_co <- c(159, 41, 367, 0, 3341, 61, 2584)
  table(IM_ATAC_IB$predictedGroup_State_Un)
  temp$ATAC_num_un <- c(344, 30, 242, 0, 3739, 313, 1885)
  
  temp$RNA_percent <- temp$RNA_num / sum(temp$RNA_num)
  temp$ATAC_percent_co <- temp$ATAC_num_co / sum(temp$ATAC_num_co)
  temp$ATAC_percent_un <- temp$ATAC_num_un / sum(temp$ATAC_num_un)
  
  ggplot(data = temp, aes(x = RNA_percent, y = ATAC_percent_co)) +
    geom_point(data = temp, aes(x = RNA_percent, y = ATAC_percent_co, color = State)) +
    geom_smooth(method = "lm", colour = "black") +
    ggpubr::stat_cor() +
    scale_color_manual(values = state_color)
  
  pdf("Figure5_2_state_percent_cor.pdf", width = 5.5, height = 5)
  ggplot(data = temp, aes(x = RNA_percent, y = ATAC_percent_un)) +
    geom_point(data = temp, aes(x = RNA_percent, y = ATAC_percent_un, color = State),
               size = 3) +
    geom_smooth(method = "lm", colour = "black") +
    ggpubr::stat_cor() +
    scale_color_manual(values = state_color) +
    scale_y_continuous(labels = percent) +
    scale_x_continuous(labels = percent) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    labs(x = "Trajectory state percent", y = "Integration predicted state percent")
  dev.off()
  
  ##### 用非限制性整合的结果
  markersPeaks_State_IB <- getMarkerFeatures(
    ArchRProj = IM_ATAC_IB,
    useMatrix = "PeakMatrix",
    groupBy = "predictedGroup_State_Un",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  markerPeaksList_State_IB <- getMarkers(markersPeaks_State_IB, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  state_1_marker_peaks <- markerPeaksList_State_IB@listData[["1"]]
  state_1_marker_peaks <- as.data.frame(state_1_marker_peaks)
  state_1_marker_peaks$ID <- paste0(state_1_marker_peaks$seqnames, "_", state_1_marker_peaks$start, "_", state_1_marker_peaks$end)
  state_1_marker_peaks$nearestGene <- peaksets[state_1_marker_peaks$ID, "nearestGene"]
  state_1_marker_peaks$peakType <- peaksets[state_1_marker_peaks$ID, "peakType"]
  state_7_marker_peaks <- markerPeaksList_State_IB@listData[["7"]]
  state_7_marker_peaks <- as.data.frame(state_7_marker_peaks)
  state_7_marker_peaks$ID <- paste0(state_7_marker_peaks$seqnames, "_", state_7_marker_peaks$start, "_", state_7_marker_peaks$end)
  state_7_marker_peaks$nearestGene <- peaksets[state_7_marker_peaks$ID, "nearestGene"]
  state_7_marker_peaks$peakType <- peaksets[state_7_marker_peaks$ID, "peakType"]
  state_5_marker_peaks <- markerPeaksList_State_IB@listData[["5"]]
  state_5_marker_peaks <- as.data.frame(state_5_marker_peaks)
  state_5_marker_peaks$ID <- paste0(state_5_marker_peaks$seqnames, "_", state_5_marker_peaks$start, "_", state_5_marker_peaks$end)
  state_5_marker_peaks$nearestGene <- peaksets[state_5_marker_peaks$ID, "nearestGene"]
  state_5_marker_peaks$peakType <- peaksets[state_5_marker_peaks$ID, "peakType"]
  
  
  
  IM_ATAC_IB_GeneIntegrationMatrix <- getMatrixFromProject(IM_ATAC_IB,
                                                           useMatrix = "GeneIntegrationMatrix_State_Un")
  IM_ATAC_IB_GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix_State_Un"]]@Dimnames[[1]] <- IM_ATAC_IB_GeneIntegrationMatrix@elementMetadata@listData[["name"]]
  IM_ATAC_IB_GeneIntegrationMatrix <- IM_ATAC_IB_GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix_State_Un"]]
  
  RNA_pseudotime <- data.frame(cells = colnames(mycds),
                               pseudotime = mycds$Pseudotime,
                               row.names = colnames(mycds))
  IM_ATAC_IB$monocle2_pseudotime <- RNA_pseudotime[IM_ATAC_IB$predictedCell_State_Un, "pseudotime"]
  IM_ATAC_IB_meta <- getCellColData(IM_ATAC_IB)
  IM_ATAC_IB_meta <- as.data.frame(IM_ATAC_IB_meta)
  IM_ATAC_IB_meta <- IM_ATAC_IB_meta[colnames(IM_ATAC_IB_GeneIntegrationMatrix),]
  data <- as.matrix(IM_ATAC_IB_GeneIntegrationMatrix)
  which(rowSums(data) == 0)
  data <- data[-which(rowSums(data) == 0),]
  data <- round(data)
  data <- as(data, 'sparseMatrix')
  which(rowSums(data) == 0)
  data <- data[-which(rowSums(data) == 0),]
  temp <- CreateSeuratObject(data)
  temp <- NormalizeData(temp)
  temp <- FindVariableFeatures(temp)
  temp <- AddMetaData(temp, metadata = IM_ATAC_IB_meta)
  temp <- ScaleData(temp)
  temp <- RunPCA(temp)
  temp <- RunUMAP(temp, dims = 1:30)
  DimPlot(temp, group.by = "predictedGroup_State_Un", pt.size = 1)
  DimPlot(temp, group.by = "predictedGroup_State_Un", split.by = "predictedGroup_State_Un")
  FeaturePlot(temp, features = "monocle2_pseudotime", split.by = "Celltype")
  
  
  # 非限制性整合 monocle3
  library(monocle3)
  data <- GetAssayData(temp, assay = 'RNA', slot = 'counts')
  cell_metadata <- temp@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds_monocle3 <- new_cell_data_set(data,
                                    cell_metadata = cell_metadata,
                                    gene_metadata = gene_annotation)
  #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
  cds_monocle3 <- preprocess_cds(cds_monocle3, num_dim = 30)
  #umap降维
  cds_monocle3 <- reduce_dimension(cds_monocle3, preprocess_method = "PCA")
  plot_cells(cds_monocle3, reduction_method = "UMAP",
             color_cells_by = "predictedGroup_State_Un", cell_size = 1) +
    ggtitle('cds.umap')
    
  ##从seurat导入整合过的umap坐标
  cds_monocle3.embed <- cds_monocle3@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(temp, reduction = "umap")
  int.embed <- int.embed[rownames(cds_monocle3.embed),]
  cds_monocle3@int_colData$reducedDims$UMAP <- int.embed
  plot_cells(cds_monocle3, reduction_method = "UMAP",
             color_cells_by = "predictedGroup_State_Un", cell_size = 1) +
    ggtitle('cds.umap')
  plot_cells(cds_monocle3, reduction_method = "UMAP",
             color_cells_by = "Celltype", cell_size = 1) +
    ggtitle('cds.umap')
  cds_monocle3 <- cluster_cells(cds_monocle3, resolution = 0.005)
  plot_cells(cds_monocle3, show_trajectory_graph = FALSE, cell_size = 1,
             group_label_size = 5) + ggtitle("label by clusterID")
  plot_cells(cds_monocle3, color_cells_by = "partition",
             show_trajectory_graph = FALSE,
             cell_size = 1,
             group_label_size = 5) + 
    ggtitle("label by partitionID")
  cds_monocle3 <- learn_graph(cds_monocle3, close_loop = FALSE, use_partition = FALSE)
  plot_cells(cds_monocle3, label_groups_by_cluster = TRUE, label_leaves = FALSE, 
             label_branch_points = TRUE, cell_size = 1, group_label_size = 5)
  cds_monocle3 <- order_cells(cds_monocle3)
  plot_cells(cds_monocle3, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, label_leaves = FALSE, 
             label_branch_points = TRUE, cell_size = 1, group_label_size = 5)
  cds_monocle3_pseudotime <- pseudotime(cds_monocle3, reduction_method = 'UMAP')
  
  temp$monocle3_pseudotime <- cds_monocle3_pseudotime[colnames(temp)]
  IM_ATAC_IB$monocle3_pseudotime <- cds_monocle3_pseudotime[IM_ATAC_IB$cellNames]
  
  ggplot(data = IM_ATAC_IB_meta, aes(x = monocle2_pseudotime, y = monocle3_pseudotime)) +
    geom_point() +
    ggpubr::stat_cor(method = "spearman")
  
  Idents(temp) <- temp$predictedGroup_State_Un
  temp_state_markers <- FindAllMarkers(temp, only.pos = T)
  temp_state_markers_2 <- temp_state_markers[which(temp_state_markers$avg_log2FC > 0.5),]
  # temp_state_markers_2 <- temp_state_markers[which(temp_state_markers$cluster %in% c(1,7,5)),]
  # temp_state_markers_2 <- temp_state_markers_2[which(temp_state_markers_2$avg_log2FC > 0.5),]
  pd <- new('AnnotatedDataFrame', data = IM_ATAC_IB_meta)
  fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  identical(rownames(fd), rownames(data))
  # fd@data[["gene_short_name"]] <- rownames(data)
  # rownames(fd) <- rownames(data)
  ## 以下代码一律不得修改
  mycds_3 <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily = negbinomial.size()) # lowerDetectionLimit = 0.1
  
  mycds_3 <- estimateSizeFactors(mycds_3)
  mycds_3 <- estimateDispersions(mycds_3, cores = 4)
  disp_table_3 <- dispersionTable(mycds_3)
  rownames(disp_table_3) <- disp_table_3$gene_id
  disp_table_3[,"markers"] <- ifelse(disp_table_3[,"gene_id"] %in% temp_state_markers_2$gene,
                                     "YES", "NO")
  disp.genes_3 <- subset(disp_table_3, markers == "YES" & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  VariableFeatures(temp)
  # mycds_3 <- setOrderingFilter(mycds_3, unique(temp_state_markers_2$gene))
  mycds_3 <- setOrderingFilter(mycds_3, disp.genes_3)
  plot_ordering_genes(mycds_3)
  mycds_3 <- reduceDimension(mycds_3, max_components = 2, reduction_method = 'DDRTree')
  mycds_3 <- orderCells(mycds_3)
  plot_cell_trajectory(mycds_3, color_by = "State") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  mycds_3 <- orderCells(mycds_3, root_state = "5")
  plot_cell_trajectory(mycds_3, color_by = "Celltype") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  plot_cell_trajectory(mycds_3, color_by = "predictedGroup_State_Un") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  
  
  
  
  IM_ATAC_IB_GeneIntegration <- CreateSeuratObject(IM_ATAC_IB_GeneIntegrationMatrix)
  IM_ATAC_IB_GeneIntegration$DataType <- "snATAC"
  IM_ATAC_IB_meta <- getCellColData(IM_ATAC_IB)
  IM_ATAC_IB_meta <- as.data.frame(IM_ATAC_IB_meta)
  IM_ATAC_IB_GeneIntegration$Celltype <- IM_ATAC_IB_meta[colnames(IM_ATAC_IB_GeneIntegration), "Celltype"]
  table(IM_ATAC_IB_GeneIntegration$Celltype)
  IM_ATAC_IB_GeneIntegration$DataType_Celltype <- paste0(IM_ATAC_IB_GeneIntegration$DataType,
                                                         "_",
                                                         IM_ATAC_IB_GeneIntegration$Celltype)
  IM_ATAC_IB_GeneIntegration <- NormalizeData(IM_ATAC_IB_GeneIntegration)
  Meristem$DataType <- "snRNA"
  Meristem$DataType_Celltype <- paste0(Meristem$DataType,
                                       "_",
                                       Meristem$Celltype)
  Meristem_ATAC_RNA_merge <- merge(IM_ATAC_IB_GeneIntegration,
                                   Meristem)
  Meristem_ATAC_RNA_merge <- NormalizeData(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- FindVariableFeatures(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- ScaleData(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- RunPCA(Meristem_ATAC_RNA_merge)
  Meristem_ATAC_RNA_merge <- RunUMAP(Meristem_ATAC_RNA_merge,
                                     dims = 1:50)
  DimPlot(Meristem_ATAC_RNA_merge,
          group.by = "DataType")
  DimPlot(Meristem_ATAC_RNA_merge, 
          group.by = "Celltype")

  features <- SelectIntegrationFeatures(object.list = list(IM_ATAC_IB_GeneIntegration,
                                                           Meristem))
  anchors <- FindIntegrationAnchors(object.list = list(IM_ATAC_IB_GeneIntegration,
                                                       Meristem),
                                    anchor.features = features)
  Meristem_ATAC_RNA_integrate <- IntegrateData(anchorset = anchors)
  DefaultAssay(Meristem_ATAC_RNA_integrate) <- "integrated"
  Meristem_ATAC_RNA_integrate <- ScaleData(Meristem_ATAC_RNA_integrate, verbose = FALSE)
  Meristem_ATAC_RNA_integrate <- RunPCA(Meristem_ATAC_RNA_integrate, npcs = 50, verbose = FALSE)
  Meristem_ATAC_RNA_integrate <- RunUMAP(Meristem_ATAC_RNA_integrate, reduction = "pca", dims = 1:50)
  pdf("Figure5_2_Meristem_ATAC_RNA_integrate_celltype_UMAP.pdf", width = 6.7, height = 5)
  DimPlot(Meristem_ATAC_RNA_integrate,
          group.by = "Celltype", cols = celltype_color)
  dev.off()
  pdf("Figure5_2_Meristem_ATAC_RNA_integrate_datatype_UMAP.pdf", width = 5.5, height = 5)
  DimPlot(Meristem_ATAC_RNA_integrate,
          group.by = "DataType")
  dev.off()
  
  Meristem_ATAC_RNA_integrate@meta.data[IM_ATAC_IB$cellNames, "State"] <- IM_ATAC_IB$predictedGroup_State_Un
  table(Meristem_ATAC_RNA_integrate$DataType, Meristem_ATAC_RNA_integrate$State)
  pdf("Figure5_2_Meristem_ATAC_RNA_integrate_state_UMAP.pdf", width = 5, height = 5)
  DimPlot(Meristem_ATAC_RNA_integrate,
          group.by = "State", cols = state_color)
  dev.off()
  
  IM_ATAC_IB_meta <- getCellColData(IM_ATAC_IB)
  IM_ATAC_IB_meta <- as.data.frame(IM_ATAC_IB_meta)
  IM_ATAC_IB <- addBgdPeaks(IM_ATAC_IB, force = T)
  IM_ATAC_IB <- addDeviationsMatrix(
    ArchRProj = IM_ATAC_IB, 
    peakAnnotation = "TF-Motif",
    force = TRUE,
    matrixName = "TF-deviation"
  )
  IM_IB_TF_dev <- getMatrixFromProject(IM_ATAC_IB, useMatrix = "TF-deviation")
  IM_IB_TF_dev_matrix <- IM_IB_TF_dev@assays@data@listData[["z"]]
  IM_IB_TF_dev_matrix <- as.data.frame(IM_IB_TF_dev_matrix)
  IM_IB_TF_dev_matrix <- as.data.frame(t(IM_IB_TF_dev_matrix))
  IM_IB_TF_dev_matrix <- IM_IB_TF_dev_matrix[,clustering$Genes]
  IM_IB_TF_dev_matrix$State <- IM_ATAC_IB_meta[rownames(IM_IB_TF_dev_matrix), "predictedGroup_State_Un"]
  IM_IB_TF_dev_matrix$State <- paste0("State ", IM_IB_TF_dev_matrix$State)
  IM_IB_TF_dev_matrix$State <- factor(IM_IB_TF_dev_matrix$State,
                                      levels = c("State 5", "State 6", "State 3", "State 2", "State 1", "State 7"))
  State_branch <- data.frame(State = c("State 5", "State 6", "State 3", "State 2", "State 1", "State 7"),
                             row.names = c("State 5", "State 6", "State 3", "State 2", "State 1", "State 7"),
                             Branch = c("Pre-branch", "Pre-branch", "Pre-branch", "Pre-branch", "Branch 1", "Branch 2"))
  
  IM_IB_TF_dev_matrix$Branch <- State_branch[IM_IB_TF_dev_matrix$State, "Branch"]
  IM_IB_TF_dev_matrix$Branch <- factor(IM_IB_TF_dev_matrix$Branch,
                                       levels = c("Pre-branch", "Branch 1", "Branch 2"))
  gplot1_3 <- function(DATA, gene) {
    ggplot(DATA, aes(x = State, y = DATA[,gene], fill = Branch)) +
      geom_boxplot() +
      scale_fill_manual(values = c("gray", "#f05662", "#7990c8")) +
      theme_bw() +
      ggtitle(gene) +
      theme(axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 10, color = "black"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            plot.title = element_text(size = 12, hjust = 0.5, face = "bold")) +
      labs(x = "", y = "TF motif deviations")
  }
  
  gplot1_3(IM_IB_TF_dev_matrix, clustering$Genes[4])
  
  
  
}


# 2016-PJ
meta <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/2016-PJ/SraRunTable.txt",
                 sep = ",", header = T)
rownames(meta) <- meta$Run
PJ2016_count <- c()
for (i in meta$Run) {
  print(i)
  temp <- read.csv(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/2016-PJ/4.STAR_Count/",
                          i, "ReadsPerGene.out.tab"), sep = "\t", header = F)
  temp <- temp[-c(1:4),]
  temp <- temp[,c(1,2)]
  rownames(temp) <- temp$V1
  colnames(temp)[2] <- paste0(i, "_", meta[i, "Sample.Name"])
  PJ2016_count <- c(PJ2016_count, list(temp))
}

PJ2016_count <- reduce(PJ2016_count, merge)
PJ2016_count$V1 <- MSU_RAP[match(PJ2016_count$V1, MSU_RAP$MSU), "RAP"]
PJ2016_count <- as.data.frame(na.omit(PJ2016_count))
PJ2016_count <- PJ2016_count[!duplicated(PJ2016_count$V1),]
rownames(PJ2016_count) <- PJ2016_count$V1
PJ2016_count <- PJ2016_count[,-1]

library(DESeq2)

#创建配对比较的列表信息 
group <- unlist(lapply(colnames(PJ2016_count), function(x){
  unlist(strsplit(x, "_"))[2]
}))
names(group) <- colnames(PJ2016_count)
# 创建需要配对比较的列表函数,创建了三个分组。
createList <- function(group=NULL) {
  sample <- names(group)
  sampleList <- list()
  treatsamList <- list()
  treatnameList <- c()
  ctrlnameList <- c()
  #A-1: PBM vs 其他
  sampleList[[1]] <- sample
  treatsamList[[1]] <- intersect(sample, names(group[group=="PBM"])) # 亚型名称需要根据情况修改
  treatnameList[1] <- "PBM" # 该亚型的命名
  ctrlnameList[1] <- "Others" # 其他亚型的命名
  #A-2: ePBM/AM vs 其他
  sampleList[[2]] <- sample
  treatsamList[[2]] <- intersect(sample, names(group[group=="ePBM/AM"]))
  treatnameList[2] <- "ePBM_AM"
  ctrlnameList[2] <- "Others"
  #A-3: SM vs 其他
  sampleList[[3]] <- sample
  treatsamList[[3]] <- intersect(sample, names(group[group=="SM"]))
  treatnameList[3] <- "SM"
  ctrlnameList[3] <- "Others"
  #A-4: RM vs 其他
  sampleList[[4]] <- sample
  treatsamList[[4]] <- intersect(sample, names(group[group=="RM"]))
  treatnameList[4] <- "RM"
  ctrlnameList[4] <- "Others"
  #如果有更多类，按以上规律继续写
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}
complist <- createList(group = group)
# 配对DESeq2函数
twoclassDESeq2 <- function(res.path=NULL, countsTable=NULL, prefix=NULL, complist=NULL, overwt=FALSE) {
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  options(warn = 1)
  for (k in 1:length(sampleList)) { # 循环读取每一次比较的内容
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    compname <- paste(treatname, "_vs_", ctrlname, sep="") # 生成最终文件名
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(prefix, "_deseq2_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) { # 因为差异表达分析较慢，因此如果文件存在，在不覆盖的情况下（overwt=F）不再次计算差异表达
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    saminfo <- data.frame("Type"=tmp[samples],"SampleID"=samples,stringsAsFactors = F)
    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]
    # 差异表达过程，具体参数细节及输出结果解释，请参阅相关document
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = as.formula("~ Type")) # 设计矩阵仅包含亚型信息
    dds$Type <- relevel(dds$Type,ref = "control")
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("Type","treatment","control"))
    #将分析结果转化为数据框
    resData <- as.data.frame(res[order(res$padj),])
    #将行名作为id列
    resData$id <- rownames(resData)
    #提取想要的列数据
    resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    #修改列名
    colnames(resData) <- c("id","baseMean","log2FC","lfcSE","stat","PValue","FDR")
    #输出到文件
    write.table(resData, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
}

twoclassDESeq2(res.path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/2016-PJ/5.DESeq2", #所有配对差异表达结果都会输出在res.path路径下
               countsTable = PJ2016_count,
               prefix = "2016-PJ", #文件名以SKCM开头
               complist = complist,
               overwt = F)

ePBM_AM <- read.csv("./Ref/2016-PJ/5.DESeq2/2016-PJ_deseq2_test_result.ePBM_AM_vs_Others.txt",
                    sep = "\t", header = T)
ePBM_AM <- ePBM_AM[which(ePBM_AM$PValue < 0.05 & ePBM_AM$log2FC > 1),]

PBM <- read.csv("./Ref/2016-PJ/5.DESeq2/2016-PJ_deseq2_test_result.PBM_vs_Others.txt",
                sep = "\t", header = T)
PBM <- PBM[which(PBM$PValue < 0.05 & PBM$log2FC > 1),]

RM <- read.csv("./Ref/2016-PJ/5.DESeq2/2016-PJ_deseq2_test_result.RM_vs_Others.txt",
                sep = "\t", header = T)
RM <- RM[which(RM$PValue < 0.05 & RM$log2FC > 1),]

SM <- read.csv("./Ref/2016-PJ/5.DESeq2/2016-PJ_deseq2_test_result.SM_vs_Others.txt",
               sep = "\t", header = T)
SM <- SM[which(SM$PValue < 0.05 & SM$log2FC > 1),]

rownames(gene_id_map) <- gene_id_map$gene_id
ePBM_AM$Symbol <- gene_id_map[ePBM_AM$id, "symbol"]
ePBM_AM <- as.data.frame(na.omit(ePBM_AM))

PBM$Symbol <- gene_id_map[PBM$id, "symbol"]
PBM <- as.data.frame(na.omit(PBM))

RM$Symbol <- gene_id_map[RM$id, "symbol"]
RM <- as.data.frame(na.omit(RM))

SM$Symbol <- gene_id_map[SM$id, "symbol"]
SM <- as.data.frame(na.omit(SM))

rownames(Meristem)
Meristem <- AddModuleScore(Meristem, features = list(ePBM_AM = ePBM_AM$Symbol,
                                                     PBM = PBM$Symbol,
                                                     RM = RM$Symbol,
                                                     SM = SM$Symbol),
                           search = T)
colnames(Meristem@meta.data)[26:29] <- c("ePBM/AM", "PBM", "RM", "SM")

library(AUCell)
cells_rankings <- AUCell_buildRankings(Meristem@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(list(ePBM_AM = ePBM_AM$Symbol,
                                 PBM = PBM$Symbol,
                                 RM = RM$Symbol,
                                 SM = SM$Symbol),
                            cells_rankings, aucMaxRank = nrow(cells_rankings)*0.5)
cells_AUC <- cells_AUC@assays@data@listData[["AUC"]]
cells_AUC <- as.data.frame(t(cells_AUC))
colnames(cells_AUC)

mycds$`ePBM_AM` <- cells_AUC[colnames(mycds), "ePBM_AM"]
mycds$`PBM` <- cells_AUC[colnames(mycds), "PBM"]
mycds$`RM` <- cells_AUC[colnames(mycds), "RM"]
mycds$`SM` <- cells_AUC[colnames(mycds), "SM"]

plot_cell_trajectory(mycds, color_by = "RM") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "gray", high = "red")

Data <- data.frame(State = mycds$State,
                   ePBM_AM = mycds$ePBM_AM,
                   PBM = mycds$PBM,
                   RM = mycds$RM,
                   SM = mycds$SM)
Data$State <- as.character(Data$State)
Data$State[which(Data$State %in% c("5","6","4","3","2"))] <- "Pre-branch"
Data$State[which(Data$State %in% c("7"))] <- "Branch 2"
Data$State[which(Data$State %in% c("1"))] <- "Branch 1"
Data_melt <- melt(Data)
colnames(Data_melt) <- c("Branch", "Stage", "AUCell")
Data_melt$Stage <- factor(Data_melt$Stage,
                          levels = c("RM", "PBM", "ePBM_AM", "SM"))
AUCell_Stage <- ggplot(data = Data_melt, aes(x = Branch, y = AUCell)) +
  geom_violin(aes(fill = Branch, color = Branch),
              trim = FALSE) +
  geom_boxplot(width = 0.1, color = "white", fill = "black", outlier.shape = NA, linewidth = 0.8) +
  facet_wrap(. ~ Stage, scales = "free",) +
  ggpubr::stat_compare_means(comparisons = list(c("Pre-branch", "Branch 1"),
                                                c("Pre-branch", "Branch 2"),
                                                c("Branch 1", "Branch 2"))) +
  scale_color_manual(values = c("Pre-branch" = "#979797", "Branch 1" = "#f05662", "Branch 2" = "#7990c8")) +
  scale_fill_manual(values = c("Pre-branch" = "#979797", "Branch 1" = "#f05662", "Branch 2" = "#7990c8")) +
  labs(x = "", y = "AUCell Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
pdf("Figure5_3_Stage_AUCell_Branch.pdf", width = 6, height = 9)
AUCell_Stage
dev.off()

pdf("Figure5_3_Tissue_Batch.pdf", width = 8.5, height = 4.5)
plot_cell_trajectory(mycds, color_by = "Batch") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(.~tissues)
dev.off()

pdf("Figure5_3_nCount_RNA.pdf", width = 5.5, height = 4.5)
plot_cell_trajectory(mycds, color_by = "nCount_RNA") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figure5_3_nFeature_RNA.pdf", width = 5.5, height = 4.5)
plot_cell_trajectory(mycds, color_by = "nFeature_RNA") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figure5_3_percent.mt.pdf", width = 5.5, height = 4.5)
plot_cell_trajectory(mycds, color_by = "percent.mt") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figure5_3_percent.C.pdf", width = 5.5, height = 4.5)
plot_cell_trajectory(mycds, color_by = "percent.C") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdata <- as.data.frame(mycds@phenoData@data)
pdata <- pdata[,c("orig.ident", "Batch", "Pseudotime", "State")]
pdata <- pdata[order(pdata$Pseudotime, decreasing = F),]
nrow(pdata) / 100
Pseudobulk <- rep(1:100, each = 55)
Pseudobulk <- c(Pseudobulk, rep(55, 28))
pdata$Pseudobulk <- Pseudobulk

pdata$Pseudobulk <- ceiling(pdata$Pseudotime)
pdata$Pseudobulk[1] <- 1
pdata$Pseudobulk_label <- paste0(pdata$Pseudobulk - 1, "-", pdata$Pseudobulk)
Pseudobulk_temp <- data.frame(Pseudobulk = c("0-1", "1-2", "2-3", "3-4", "4-5",
                                             "5-6", "6-7", "7-8", "8-9", "9-10",
                                             "10-11", "11-12", "12-13", "13-14",
                                             "14-15", "15-16", "17-18", "18-19", "19-20"),
                              label = c("[0-1]", "(1-2]", "(2-3]", "(3-4]", "(4-5]",
                                        "(5-6]", "(6-7]", "(7-8]", "(8-9]", "(9-10]",
                                        "(10-11]", "(11-12]", "(12-13]", "(13-14]",
                                        "(14-15]", "(15-16]", "(17-18]", "(18-19]", "(19-20]"))
rownames(Pseudobulk_temp) <- Pseudobulk_temp$Pseudobulk
pdata$Pseudobulk_label <- Pseudobulk_temp[pdata$Pseudobulk_label, "label"]
pdata$Pseudobulk_label <- factor(pdata$Pseudobulk_label, levels = c("[0-1]", "(1-2]", "(2-3]", "(3-4]", "(4-5]",
                                                                    "(5-6]", "(6-7]", "(7-8]", "(8-9]", "(9-10]",
                                                                    "(10-11]", "(11-12]", "(12-13]", "(13-14]",
                                                                    "(14-15]", "(15-16]", "(17-18]", "(18-19]", "(19-20]"))
temp <- as.data.frame.array(table(pdata$Batch, pdata$Pseudobulk_label))
temp <- as.data.frame(t(temp))
temp$Label <- rownames(temp)
temp$Label <- factor(temp$Label, levels = temp$Label)

pdf("Figure5_3_Pseudotime_library_count_IM0.5cm.pdf", width = 5.5, height = 4.5)
ggplot() +
  geom_point(data = temp, aes(x = R1077, y = R1078, color = Label),
             size = 3.5) +
  theme_bw() +
  ggtitle("IM0.5cm") +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  geom_smooth(data = temp, aes(x = R1077, y = R1078),
              method = "lm", fill = "#0000CD") +
  ggpubr::stat_cor(data = temp, aes(x = R1077, y = R1078)) +
  labs(color = "Pseudotime", x = "Count of cells from R1077 library", y = "Count of cells from R1078 library") +
  scale_color_manual(values = colorRampPalette(colors = c("#00008B", "#0000CD", "#FF4500"))(20))
dev.off()

pdf("Figure5_3_Pseudotime_library_count_IM1cm.pdf", width = 5.5, height = 4.5)
ggplot() +
  geom_point(data = temp, aes(x = R1079, y = R1080, color = Label),
             size = 3.5) +
  theme_bw() +
  ggtitle("IM1cm") +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  geom_smooth(data = temp, aes(x = R1079, y = R1080),
              method = "lm", fill = "#0000CD") +
  ggpubr::stat_cor(data = temp, aes(x = R1079, y = R1080)) +
  labs(color = "Pseudotime", x = "Count of cells from R1079 library", y = "Count of cells from R1080 library") +
  scale_color_manual(values = colorRampPalette(colors = c("#00008B", "#0000CD", "#FF4500"))(20))
dev.off()


ggplot() +
  geom_density(data = pdata, aes(x = Pseudotime, group = Batch, color = Batch))


#### panicle genes
my_plot_cell_trajectory <- function(mycds, seuratObject, gene) {
  lib_info_with_pseudo <- mycds@phenoData@data
  reduced_dim_coords <- reducedDimK(mycds)
  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(mycds)
  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
  
  data_df <- t(monocle::reducedDimS(mycds)) %>%
    as.data.frame() %>%
    select_(data_dim_1 = 1, data_dim_2 = 2) %>%
    rownames_to_column("sample_name") %>%
    mutate(lib_info_with_pseudo$State) %>%
    left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")
  theta <- 0
  return_rotation_mat <- function(theta) {
    theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)
  mst_branch_nodes <- mycds@auxOrderingData[[mycds@dim_reduce_type]]$branch_points
  branch_point_df <- ica_space_df %>%
    slice(match(mst_branch_nodes, sample_name)) %>%
    mutate(branch_point_idx = seq_len(n()))
  rownames(data_df) <- data_df$sample_name
  data_df$exp <- seuratObject@assays$RNA@data[gene,rownames(data_df)]
  data_df <- data_df[order(data_df$exp),]
  g <- ggplot() +
    geom_segment(aes_string(x = "source_prin_graph_dim_1", y = "source_prin_graph_dim_2",
                            xend = "target_prin_graph_dim_1", yend = "target_prin_graph_dim_2"),
                 size = 0.75, linetype = "solid", na.rm = TRUE, data = edge_df) +
    geom_point(data = data_df, aes(x = data_dim_1, y = data_dim_2, color = exp)) +
    labs(x = "Component 1", y = "Component 2") 
  g <- g +
    geom_point(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"),
               size = 5, na.rm = TRUE, branch_point_df) +
    geom_text(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "branch_point_idx"),
              size = 4, color = "white", na.rm = TRUE, branch_point_df)
  
  return(g)
  
}

PanicleGenes <- xlsx::read.xlsx("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Panicle genes/PanicleGenes.xls",
                                sheetIndex = 1)
intersect(PanicleGenes$Locus.ID, state_1$Gene_id)
intersect(PanicleGenes$Locus.ID, state_7$Gene_id)

FeaturePlot(Meristem, features = gene_id_map["Os01g0140100",2], order = T)



pdf("Figure5_3_Os01g0140100 | Os01g0140100.pdf", width = 5.5, height = 4.5)
my_plot_cell_trajectory(mycds, Meristem, "Os01g0140100") +
  ggtitle("Os01g0140100 | Os01g0140100") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  scale_colour_gradientn(colours = c("gray", "#8FBC8F", "#FA8072", "#FF0000")) +
  labs(color = "Expression")
dev.off()

pdf("Figure5_3_Os09g0507100 | OsSPL18.pdf", width = 5.5, height = 4.5)
my_plot_cell_trajectory(mycds, Meristem, "OsSPL18") +
  ggtitle("Os09g0507100 | OsSPL18") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  scale_colour_gradientn(colours = c("gray", "#8FBC8F", "#FA8072", "#FF0000")) +
  labs(color = "Expression")
dev.off()
pdf("Figure5_3_Os09g0456100 | LP1.pdf", width = 5.5, height = 4.5)
my_plot_cell_trajectory(mycds, Meristem, "LP1") +
  ggtitle("Os09g0456100 | LP1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  scale_colour_gradientn(colours = c("gray", "#8FBC8F", "#FA8072", "#FF0000")) + # 20B2AA
  labs(color = "Expression")
dev.off()







PanicleGenes_2 <- xlsx::read.xlsx("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Panicle genes/Spikelet_branch_IM.xls",
                                  sheetIndex = 1)
intersect(PanicleGenes_2$Locus.ID, state_1$Gene_id)
intersect(PanicleGenes_2$Locus.ID, state_7$Gene_id)



mycds$branch <- ""
mycds$branch[which(mycds$State %in% c(2,3,4,5,6))] <- "Pre-branch"
mycds$branch[which(mycds$State %in% c(7))] <- "Branch 2"
mycds$branch[which(mycds$State %in% c(1))] <- "Branch 1"
table(mycds$branch)

mycds$Os01g0140100 <- Meristem@assays$RNA@data["Os01g0140100", colnames(mycds)]
mycds$OsLP1 <- Meristem@assays$RNA@data["LP1", colnames(mycds)]
mycds$OsSPL18 <- Meristem@assays$RNA@data["OsSPL18", colnames(mycds)]
mycds$RBPA <- Meristem@assays$RNA@data["RBPA", colnames(mycds)]

mycds$branch <- factor(mycds$branch, levels = c("Pre-branch", "Branch 1", "Branch 2"))

for (i in c("Os01g0140100", "OsLP1", "OsSPL18", "RBPA")) {
  temp <- data.frame(exp = mycds@phenoData@data[, i],
                     branch = mycds$branch)
  pdf(paste0("Figure5_C_", i, "_branch_expression.pdf"), width = 2, height = 3.5)
  p <- ggplot(temp, aes(x = branch, y = exp, fill = branch)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_violin() +
    ggpubr::stat_compare_means(aes(label = after_stat(p.signif)),
                               comparisons = list(c("Branch 1", "Branch 2"),
                                                  c("Branch 1", "Pre-branch"),
                                                  c("Branch 2", "Pre-branch"))) +
    theme_bw() +
    ggtitle(i) +
    scale_fill_manual(values = c("Branch 1" = "#f05662", "Branch 2" = "#7990c8", "Pre-branch" = "#979797")) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(size = 10, color = "black", hjust = 0.5)) +
    NoLegend() +
    labs(x = "", y = "Relative expression level")
  print(p)
  dev.off()
}

FeaturePlot(Meristem, features = gene_id_map["Os07g0605200",2], order = T)
pdf("Figure5_3_Os03g0333200 | DRUS1", width = 5.5, height = 4.5) # 细胞死亡，和GO符合
my_plot_cell_trajectory(mycds, Meristem, "DRUS1") +
  ggtitle("Os03g0333200 | DRUS1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  scale_colour_gradientn(colours = c("gray", "#8FBC8F", "#FA8072", "#FF0000")) +
  labs(color = "Expression")
dev.off()
pdf("Figure5_3_Os01g0769700 | DRUS2.pdf", width = 5.5, height = 4.5)
my_plot_cell_trajectory(mycds, Meristem, "DRUS2") +
  ggtitle("Os01g0769700 | DRUS2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  scale_colour_gradientn(colours = c("#20B2AA", "#8FBC8F", "#FA8072", "#FF0000")) +
  labs(color = "Expression")
dev.off()
pdf("Figure5_3_Os11g0637700 | RBPA.pdf", width = 5.5, height = 4.5)
my_plot_cell_trajectory(mycds, Meristem, "RBPA") +
  ggtitle("Os11g0637700 | RBPA") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black")) +
  scale_colour_gradientn(colours = c("gray", "#8FBC8F", "#FA8072", "#FF0000")) +
  labs(color = "Expression")
dev.off()



# D11
FeaturePlot(Meristem, features = gene_id_map["Os04g0469800",2], order = T)

# MFO1
FeaturePlot(Meristem, features = gene_id_map["Os02g0682200",2], order = T)

# SH5
FeaturePlot(Meristem, features = gene_id_map["Os05g0455200",2], order = T)

#OsLIS-L1
FeaturePlot(Meristem, features = gene_id_map["Os08g0162100",2], order = T)


######################################
######################################
######################################
######################################
######################################


TF_regulation <- readRDS("TF_regulation_filtered.rds")
TF_regulation_branch_1 <- TF_regulation[which(TF_regulation$Target %in% branch_1$Gene_id),]
TF_regulation_branch_1 <- TF_regulation_branch_1[which(TF_regulation_branch_1$TF %in% branch_1$Gene_id),]
TF_regulation_branch_1$TF_symbol <- Annotation_genes[TF_regulation_branch_1$TF, 2]
TF_regulation_branch_1$Target_symbol <- Annotation_genes[TF_regulation_branch_1$Target, 2]
write.table(TF_regulation_branch_1, "Figure5_TF_regulation_branch_1.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
TF_regulation_branch_1_node <- data.frame(node = c(unique(TF_regulation_branch_1$TF_symbol),
                                                   TF_regulation_branch_1$Target_symbol),
                                          type = c("TF", rep("Target", 7)))
write.table(TF_regulation_branch_1_node, "Figure5_TF_regulation_branch_1_node.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
branch_1[which(branch_1$Gene_id == "Os06g0194000"),]

TF_regulation_branch_2 <- TF_regulation[which(TF_regulation$Target %in% branch_2$Gene_id),]
TF_regulation_branch_2 <- TF_regulation_branch_2[which(TF_regulation_branch_2$TF %in% branch_2$Gene_id),]
TF_regulation_branch_2$TF_symbol <- Annotation_genes[TF_regulation_branch_2$TF, 2]
TF_regulation_branch_2$Target_symbol <- Annotation_genes[TF_regulation_branch_2$Target, 2]
write.table(TF_regulation_branch_2, "Figure5_TF_regulation_branch_2.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
TF_regulation_branch_2_node <- data.frame(node = c(unique(TF_regulation_branch_2$TF_symbol),
                                                   TF_regulation_branch_2$Target_symbol),
                                          type = c("TF", "TF", rep("Target", 17)))
write.table(TF_regulation_branch_2_node, "Figure5_TF_regulation_branch_2_node.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
branch_2[which(branch_2$Gene_id == "Os01g0229000"),]
BEAM_res["Os01g0229000",]
branch_2[which(branch_2$Gene_id == "Os01g0302500"),]


branch_TF <- Meristem@assays$RNA@data[c("ERF71",
                                        "Os01g0229000",
                                        "OSH6"),colnames(mycds)]
branch_TF <- as.matrix(branch_TF)
branch_TF <- t(branch_TF)[colnames(mycds),]
branch_TF_scaled <- scale(branch_TF)

mycds_meta <- mycds@phenoData@data
mycds_meta <- mycds_meta[,1:25]
mycds_meta <- as.data.frame(cbind(mycds_meta,
                                  branch_TF_scaled[rownames(mycds_meta),]))
pd <- new('AnnotatedDataFrame', data = mycds_meta)
mycds@phenoData <- pd

pdf("Figure5_cell_trajectory_OsERF71.pdf", width = 7.5, height = 5.5)
plot_cell_trajectory(mycds, color_by = "ERF71") +
  scale_color_viridis(option = "D") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  labs(color = "Log-Normalized Expression") +
  annotate("text", x = -4.3, y = 7.3, label = "OsERF71", size = 3.3)
dev.off()

pdf("Figure5_cell_trajectory_Os01g0229000.pdf", width = 7.5, height = 5.5)
plot_cell_trajectory(mycds, color_by = "Os01g0229000") +
  scale_color_viridis(option = "D") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  labs(color = "Log-Normalized Expression") +
  annotate("text", x = -4.3, y = 7.3, label = "OsMYB3R1L", size = 3.3)
dev.off()

pdf("Figure5_cell_trajectory_OSH6.pdf", width = 7.5, height = 5.5)
plot_cell_trajectory(mycds, color_by = "OSH6") +
  scale_color_viridis(option = "D") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  labs(color = "Log-Normalized Expression") +
  annotate("text", x = -4.3, y = 7.3, label = "OsOSH6", size = 3.3)
dev.off()

##
pdf("Figure5_genes_in_pseudotime_ERF71_branch_1.pdf", width = 5, height = 2.5)
plot_genes_in_pseudotime(cds_subset = mycds["ERF71",which(mycds$State %in% c("1","2"))],
                         color_by = "Celltype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")])
dev.off()

pdf("Figure5_genes_in_pseudotime_ERF71_branch_2.pdf", width = 5, height = 2.5)
plot_genes_in_pseudotime(cds_subset = mycds["ERF71",which(mycds$State %in% c("3","2"))],
                         color_by = "Celltype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")])
dev.off()

pdf("Figure5_genes_in_pseudotime_OSH6_branch_1.pdf", width = 5, height = 2.5)
plot_genes_in_pseudotime(cds_subset = mycds["OSH6",which(mycds$State %in% c("1","2"))],
                         color_by = "Celltype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")])
dev.off()

pdf("Figure5_genes_in_pseudotime_OSH6_branch_2.pdf", width = 5, height = 2.5)
plot_genes_in_pseudotime(cds_subset = mycds["OSH6",which(mycds$State %in% c("3","2"))],
                         color_by = "Celltype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")])
dev.off()

pdf("Figure5_genes_in_pseudotime_Os01g0229000_branch_1.pdf", width = 5, height = 2.5)
plot_genes_in_pseudotime(cds_subset = mycds["Os01g0229000",which(mycds$State %in% c("1","2"))],
                         color_by = "Celltype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")])
dev.off()

pdf("Figure5_genes_in_pseudotime_Os01g0229000_branch_2.pdf", width = 5, height = 2.5)
plot_genes_in_pseudotime(cds_subset = mycds["Os01g0229000",which(mycds$State %in% c("3","2"))],
                         color_by = "Celltype") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")])
dev.off()

##########

IM_ATAC_meristem <- IM_ATAC[which(IM_ATAC$Celltype %in% c("Branch meristems (BM)",
                                                          "Inflorescence meristem (IM)",
                                                          "Spikelet meristem (SM)")),]
p <- plotBrowserTrack(
  ArchRProj = IM_ATAC_meristem, 
  groupBy = "Celltype", 
  geneSymbol = c("OSH6","Os01g0229000","ERF71"), 
  upstream = 20000,
  downstream = 20000,
  pal = celltype_color[c("Branch meristems (BM)", "Inflorescence meristem (IM)", "Spikelet meristem (SM)")],
  baseSize = 10, facetbaseSize = 10
)
pdf("Figure5_IM_ATAC_meristem_OSH6.pdf", width = 6, height = 2.5)
grid::grid.newpage()
grid::grid.draw(p$OSH6)
dev.off()

pdf("Figure5_IM_ATAC_meristem_Os01g0229000.pdf", width = 6, height = 2.5)
grid::grid.newpage()
grid::grid.draw(p$Os01g0229000)
dev.off()

pdf("Figure5_IM_ATAC_meristem_ERF71.pdf", width = 6, height = 2.5)
grid::grid.newpage()
grid::grid.draw(p$ERF71)
dev.off()

#######################
#######################
#######################
colnames(Meristem@meta.data)
DimPlot(Meristem)
table(Meristem$State)
DimPlot(Meristem, group.by = "tissues")
FeaturePlot(Meristem, features = "tissues", order = T)
rownames(gene_id_map) <- gene_id_map$gene_id
rownames(Meristem)
FeaturePlot(Meristem, features = "OsUCL8", order = T)
DotPlot(Meristem, features = gene_id_map["Os11g0523700",2])
FeaturePlot(Meristem, features = gene_id_map["Os07g0669500",2], order = T)




#### Stage Genes
# StageGene1 <- xlsx::read.xlsx("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/IM_NewPhytologist/nph18008-sup-0005-tables4.xlsx",
#                               sheetName = "8.2006-PJ")
# 
# StageGene1[as.numeric(na.omit(match(state_1$Gene_id, StageGene1$RAP_ID))),]
# StageGene1[as.numeric(na.omit(match(state_7$Gene_id, StageGene1$RAP_ID))),]
# 
# StageGene2 <- xlsx::read.xlsx("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/IM_NewPhytologist/nph18008-sup-0005-tables4.xls",
#                               sheetName = "9.2016-PJ")
# StageGene2[as.numeric(na.omit(match(state_1$Gene_id, StageGene2$RapID))),]
# StageGene2[as.numeric(na.omit(match(state_7$Gene_id, StageGene2$RapID))),]
# table(StageGene2[as.numeric(na.omit(match(state_1$Gene_id, StageGene2$RapID))),"Name"])
# table(StageGene2[as.numeric(na.omit(match(state_7$Gene_id, StageGene2$RapID))),"Name"])

FeaturePlot(Meristem, features = gene_id_map["Os01g0140100",2], order = T)

FeaturePlot(Meristem, features = gene_id_map["Os09g0456100",2], order = T)

FeaturePlot(Meristem, features = gene_id_map["Os09g0507100",2], order = T)

FeaturePlot(Meristem, features = gene_id_map["Os07g0163500",2], order = T)

# D11
DotPlot(Meristem, features = gene_id_map["Os04g0469800",2])
FeaturePlot(Meristem, features = gene_id_map["Os04g0469800",2], order = T)

# MFO1
DotPlot(Meristem, features = gene_id_map["Os02g0682200",2])
FeaturePlot(Meristem, features = gene_id_map["Os02g0682200",2], order = T)

# SH5
DotPlot(Meristem, features = gene_id_map["Os05g0455200",2])
FeaturePlot(Meristem, features = gene_id_map["Os05g0455200",2], order = T)

#OsLIS-L1
DotPlot(Meristem, features = gene_id_map["Os08g0162100",2])
FeaturePlot(Meristem, features = gene_id_map["Os08g0162100",2], order = T)


DotPlot(Meristem, features = gene_id_map["Os01g0856500",2])
VlnPlot(Meristem, features = gene_id_map["Os01g0856500",2])
FeaturePlot(Meristem, features = gene_id_map["Os01g0856500",2], order = T)



