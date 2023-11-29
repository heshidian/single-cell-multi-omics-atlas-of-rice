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

proj_pass_filter <- readRDS("2.proj_pass_filter_annotated.rds")
unique(proj_pass_filter$Celltype)
proj_pass_filter$Celltype[which(proj_pass_filter$Celltype == "Cryptic bract/bract(cb/b)")] <- "Cryptic bract/bract (cb/b)"
proj_pass_filter$Celltype[which(proj_pass_filter$Celltype == "Inflorescence meristem(IM)")] <- "Inflorescence meristem (IM)"
proj_pass_filter$Celltype[which(proj_pass_filter$Celltype == "Branch meristems(BM)")] <- "Branch meristems (BM)"
proj_pass_filter$Celltype[which(proj_pass_filter$Celltype == "Spikelet meristem(SM)")] <- "Spikelet meristem (SM)"
unique(proj_pass_filter$Celltype)
table(proj_pass_filter$TissuesSub, proj_pass_filter$Celltype)
proj_pass_filter$TissuesSub_celltype <- paste0(proj_pass_filter$TissuesSub,
                                               "-",
                                               proj_pass_filter$Celltype)
table(proj_pass_filter$TissuesSub_celltype)
pathToMacs2 <- "/home/heshidian/mambaforge/envs/common/bin/macs2"
addArchRGenome("genome_annotation")

### TissuesSub_celltype
# table(proj_pass_filter$Celltype)
# proj_pass_filter <- addGroupCoverages(ArchRProj = proj_pass_filter, groupBy = "TissuesSub_celltype", force = TRUE)
# proj_pass_filter <- addReproduciblePeakSet( # each peak is 501 bp in length
#   ArchRProj = proj_pass_filter,
#   groupBy = "TissuesSub_celltype",
#   pathToMacs2 = pathToMacs2,
#   excludeChr = c("Mt", "Pt"),
#   force = TRUE,
#   geneAnnotation = gene_annotation,
#   genomeAnnotation = genome_annotation,
#   genomeSize = sum(genome_annotation@listData[["chromSizes"]]@ranges@width)
# )
# proj_pass_filter <- addPeakMatrix(proj_pass_filter)
# TissuesSub_celltype_peak_set <- getPeakSet(proj_pass_filter) # each peak is 501 bp in length
# getAvailableMatrices(proj_pass_filter)
# markersPeaks_TissuesSub_celltype <- getMarkerFeatures(
#   ArchRProj = proj_pass_filter,
#   useMatrix = "PeakMatrix",
#   groupBy = "TissuesSub_celltype",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# markerPeaksList_TissuesSub_celltype <- getMarkers(markersPeaks_TissuesSub_celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
# heatmapPeaks <- plotMarkerHeatmap(
#   seMarker = markersPeaks_TissuesSub_celltype,
#   cutOff = "FDR <= 0.05 & Log2FC >= 1",
#   transpose = F,
#   plotLog2FC = T,
#   labelMarkers = NULL,
#   labelRows = FALSE,
#   returnMatrix = F,
#   clusterCols = F,
#   nLabel = 1,
#   nPrint = 0
# )
# pdf("Figure2_TissuesSub_celltype_peaks_heatmap.pdf", width = 7, height = 10)
# draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# dev.off()



### celltype
set.seed(123)
proj_pass_filter$Celltype
proj_pass_filter <- addGroupCoverages(ArchRProj = proj_pass_filter, groupBy = "Celltype", force = TRUE)
GroupCoverages <- addGroupCoverages(ArchRProj = proj_pass_filter, groupBy = "Celltype", force = TRUE,
                                    returnGroups = TRUE)
proj_pass_filter <- addReproduciblePeakSet( # each peak is 501 bp in length
  ArchRProj = proj_pass_filter, 
  groupBy = "Celltype", 
  pathToMacs2 = pathToMacs2,
  excludeChr = c("Mt", "Pt"),
  force = TRUE,
  geneAnnotation = gene_annotation,
  genomeAnnotation = genome_annotation,
  genomeSize = sum(genome_annotation@listData[["chromSizes"]]@ranges@width)
)
proj_pass_filter <- addPeakMatrix(proj_pass_filter)
celltype_peak_set <- getPeakSet(proj_pass_filter) # each peak is 501 bp in length
celltype_peak_matrix <- getMatrixFromProject(proj_pass_filter, "PeakMatrix")
getAvailableMatrices(proj_pass_filter)
markersPeaks_Celltype <- getMarkerFeatures(
  ArchRProj = proj_pass_filter, 
  useMatrix = "PeakMatrix", 
  groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerPeaksList_Celltype <- getMarkers(markersPeaks_Celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks_Celltype, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = F,
  plotLog2FC = T,
  labelMarkers = NULL,
  labelRows = FALSE,
  returnMatrix = F,
  clusterCols = F,
  nLabel = 1,
  nPrint = 0
)
pdf("Figure2_celltype_peaks_heatmap.pdf", width = 7, height = 10)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

celltype_heatmapPeakMatrix <- plotMarkerHeatmap(
  seMarker = markersPeaks_Celltype, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = F,
  plotLog2FC = T,
  labelMarkers = NULL,
  labelRows = FALSE,
  clusterCols = F,
  returnMatrix = T
)

# Rice_TF_list <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/PlantTFDB_Oryza_sativa_TF_list/Osj_TF_list.txt",
#                          sep = "\t", header = T)
# rownames(Rice_TF_list) <- Rice_TF_list$Gene_ID
# TF_binding_motif_files <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs_individual")
# TF_MUS_ID <- unlist(lapply(TF_binding_motif_files, function(x){
#   unlist(strsplit(x, ".", fixed = T))[1]
# }))
# TF_ID <- Rice_TF_list[which(Rice_TF_list$Gene_ID %in% TF_MUS_ID),]
# 
# RAP_MUS <- read.csv("./Ref/RAP-MSU_2023-03-15.txt", sep = "\t", header = F)
# num <- c()
# MUS <- c()
# for (i in 1:nrow(RAP_MUS)) {
#   mus <- unlist(strsplit(RAP_MUS[i,2], ",", fixed = T))
#   mus <- unlist(lapply(mus, function(x){
#     unlist(strsplit(x, ".", fixed = T))[1]
#   }))
#   num <- c(num, rep(i, length(mus)))
#   MUS <- c(MUS, mus)
# }
# RAP_MUS <- data.frame(RAP = RAP_MUS$V1[num],
#                       MUS = MUS)
# RAP_MUS <- RAP_MUS[!duplicated(RAP_MUS),]
# RAP_MUS <- RAP_MUS[which(RAP_MUS$MUS != "None"),]
# TF_ID$TF_ID %in% RAP_MUS$MUS
# RAP_MUS <- RAP_MUS[which(RAP_MUS$MUS %in% TF_ID$TF_ID), ]
# rownames(RAP_MUS) <- RAP_MUS$MUS
# ppm2pwm <- function(ppm) {
#   ppm_matrix <- ppm@motif
#   bg <- ppm@bkg
#   pwm <- log2(ppm_matrix / bg)
#   pwm[is.infinite(pwm)] <- 0 # 频率为0的碱基取对数以后是负无穷，需要将其置为0
#   return(pwm)
# }
# pwm <- c()
# TF_binding_motif_files <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs_individual")
# for (m in TF_binding_motif_files) {
#   motif <- read_meme(file = paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs_individual/",
#                                    m))
#   motif <- convert_motifs(motif, class = "TFBSTools-PWMatrix")
#   motif <- PWMatrixList(motif)
#   # pwm_matrix <- ppm2pwm(motif)
#   temp_RAP_MUS <- RAP_MUS[which(RAP_MUS$MUS %in% motif@listData[[1]]@name),]
#   motif2 <- motif
#   for (i in 1:nrow(temp_RAP_MUS)) {
#     motif2@listData[[1]]@name <- temp_RAP_MUS[i,1]
#     pwm <- c(pwm, list(motif2))
#     names(pwm)[length(pwm)] <- temp_RAP_MUS[i,1]
#   }
# }

motif <- read_meme("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs.meme")
motif_all <- convert_motifs(motif, class = "TFBSTools-PWMatrix")
motif_all <- do.call(PWMatrixList, motif_all)
names(motif_all) <- ID(motif_all)
motif_all <- motif_all[!duplicated(names(motif_all))] # ID duplicated !!!

proj_pass_filter <- addMotifAnnotations(ArchRProj = proj_pass_filter, motifPWMs = motif_all,
                                        annoName = "TF-Motif", force = T)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks_Celltype,
  ArchRProj = proj_pass_filter,
  peakAnnotation = "TF-Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)


df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggUp

# pheatmap(enrichMotifs@assays@data@listData[["Enrichment"]],
#          scale = "row")
# pheatmap(enrichMotifs@assays@data@listData[["mlog10Padj"]],
#          scale = "row")
# plotEnrichHeatmap(enrichMotifs, transpose = TRUE)

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




Peak_sets <- getPeakSet(proj_pass_filter)
Peak_sets$Celltype <- Peak_sets@ranges@NAMES
Peak_sets$Chr <- as.character(Peak_sets@seqnames)
Peak_sets$Start <- Peak_sets@ranges@start
Peak_sets$ID <- paste0(Peak_sets$Celltype, " - ", Peak_sets$Chr, " - ", Peak_sets$Start)
Peak_sets <- as.data.frame(Peak_sets@elementMetadata)
rownames(Peak_sets) <- Peak_sets$ID

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
proj_pass_filter <- addPeakAnnotations(proj_pass_filter,
                                          list(CE),
                                          name = "CE",
                                          force = TRUE)
CE_peak_anno <- getPeakAnnotation(proj_pass_filter, name = "CE")
CE_peak_anno_match <- readRDS(CE_peak_anno[["Matches"]])
CE_peak_anno_match <- CE_peak_anno_match@assays@data@listData[["matches"]]
CE_peak_anno_match <- as.matrix(CE_peak_anno_match)
CE_peak_anno_match <- CE_peak_anno_match + 0
CE_peak_anno_match <- as.data.frame(CE_peak_anno_match)
rownames(CE_peak_anno_match) <- rownames(Peak_sets)

CE_peak_anno_match$Celltype <- Peak_sets[rownames(CE_peak_anno_match), "Celltype"]
CE_peak_anno_match$Celltype <- factor(CE_peak_anno_match$Celltype,
                                      levels = rev(sort(unique(CE_peak_anno_match$Celltype))))

CE_peak_anno_match$CE <- ifelse(CE_peak_anno_match$Region_1 == 0, "Not-matched CE", "Matched CE")
CE_peak_anno_match$CE <- factor(CE_peak_anno_match$CE,
                                levels = rev(c("Matched CE", "Not-matched CE")))

temp <- as.data.frame.array(table(CE_peak_anno_match$Celltype, CE_peak_anno_match$CE))
temp$`Matched CE` <- temp$`Matched CE` / rowSums(temp)
temp$x <- temp$`Matched CE` / 2
temp$Celltype <- rownames(temp)
temp$label <- round(temp$`Matched CE`, 3)

pdf("Figure2_conserved_element_percent.pdf", width = 7, height = 6)
ggplot() + 
  geom_bar(data = CE_peak_anno_match, aes(y = Celltype, fill = CE),
           stat = "count", position = "fill") +
  geom_text(data = temp, aes(x = x, y = Celltype, label = label), size = 3.3) +
  theme_bw(base_size = 10) +
  scale_x_continuous(labels = scales::percent) +
  # scale_fill_manual() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10)) +
  labs(x = "", y = "Percentage", fill = "Peak Type")
dev.off()



proj_pass_filter <- addPeakAnnotations(
  ArchRProj = proj_pass_filter,
  regions = list("H3K23ac_leaf" = H3K23ac_Osj_leaf_normal,
                 "H4K16ac_leaf" = H4K16ac_Osj_leaf_normal,
                 "H3K36me3_seedling" = H3K36me3_Osj_seedling_normal,
                 "H3K4me2_seedling" = H3K4me2_Osj_seedling_normal,
                 "H4K12ac_seedling" = H4K12ac_Osj_seedling_normal),
  name = "Histone_modifications",
  force = TRUE
)

Celltype_tissue <- data.frame(Celltype = proj_pass_filter$Celltype,
                              Tissue = proj_pass_filter$TissuesSub)
# colnames(Celltype_tissue) <- c("Celltype", "Tissue", "Num")

Peak_anno <- getPeakAnnotation(proj_pass_filter, name = "Histone_modifications")
Peak_anno_match <- readRDS(Peak_anno[["Matches"]])
Peak_anno_match <- Peak_anno_match@assays@data@listData[["matches"]]
Peak_anno_match <- as.matrix(Peak_anno_match)
Peak_anno_match <- Peak_anno_match + 0
Peak_anno_match <- as.data.frame(Peak_anno_match)
Peak_anno_match$ID <- rownames(Peak_sets)
Peak_anno_match <- reshape::melt(Peak_anno_match, id.vars = "ID")
Peak_anno_match$Celltype <- Peak_sets[Peak_anno_match$ID, "Celltype"]
Peak_anno_match$match <- Peak_anno_match$variable
Peak_anno_match$match <- as.character(Peak_anno_match$match)
Peak_anno_match$match[which(Peak_anno_match$value == 0)] <- "Not annotated"
ggplot(data = Peak_anno_match, aes(y = Celltype, fill = match)) +
  geom_bar(stat = "count", position = "fill")


Peak_anno_match <- Peak_anno_match[which(Peak_anno_match$value == 1),]
colnames(Peak_anno_match) <- c("Peak_id", "Histone_modifications", "value")
Peak_anno_match$Celltype <-  Peak_sets[Peak_anno_match$Peak_id, "Celltype"]

Peak_anno_match_unique <- Peak_anno_match[!duplicated(Peak_anno_match$ID),]
Peak_anno_match_unique <- as.data.frame(table(Peak_anno_match_unique$Celltype))
Peak_anno_match_unique$all <- table(Peak_sets$Celltype)[Peak_anno_match_unique$Var1]
Peak_anno_match_unique$percent <- Peak_anno_match_unique$Freq / Peak_anno_match_unique$all
Peak_anno_match_unique$Var1 <- factor(as.character(Peak_anno_match_unique$Var1),
                                      levels = rev(sort(as.character(Peak_anno_match_unique$Var1))))
Peak_anno_match_unique$Type <- "Annotated"
Peak_anno_match_unique2 <- Peak_anno_match_unique
Peak_anno_match_unique2$percent <- 1 - Peak_anno_match_unique2$percent
Peak_anno_match_unique2$Type <- "Not-annotated"

Peak_anno_match_unique <- as.data.frame(rbind(Peak_anno_match_unique,
                                              Peak_anno_match_unique2))
Peak_anno_match_unique$Type <- factor(Peak_anno_match_unique$Type,
                                      levels = c("Not-annotated", "Annotated"))
Peak_anno_match_unique_text <- Peak_anno_match_unique[which(Peak_anno_match_unique$Type == "Annotated"),]
Peak_anno_match_unique_text$percent2 <- Peak_anno_match_unique_text$percent
Peak_anno_match_unique_text$percent2 <- Peak_anno_match_unique_text$percent2 / 2
Peak_anno_match_unique_text$percent <- round(Peak_anno_match_unique_text$percent, 4)
Peak_anno_match_unique_text$percent <- Peak_anno_match_unique_text$percent * 100
Peak_anno_match_unique_text$percent <- paste0(Peak_anno_match_unique_text$percent,
                                              " %")
Peak_anno_match_percent <- ggplot() +
  geom_bar(data = Peak_anno_match_unique, aes(x = Var1, y = percent, fill = Type),
           stat = "identity", position = "fill") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  geom_text(data = Peak_anno_match_unique_text,
            aes(x = Var1, y = percent2, label = percent, hjust = 0.5)) +
  scale_fill_manual(values = c("#8FBC8F", "#F08080")) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(x = "", y = "Percentage", fill = "PeakType") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))


{
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
} # Multiple plot function

Celltype_tissue$Celltype <- factor(Celltype_tissue$Celltype,
                                   levels = rev(sort(unique(Celltype_tissue$Celltype))))
Celltype_tissue_bar <- ggplot(data = Celltype_tissue, aes(x = Celltype, fill = Tissue)) +
  geom_bar(stat = "count", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#488E44","#867C58","#618694","#A2A882","#eea29a","#f7786b")) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(x = "", y = "Percent") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

Peak_sets$Celltype <- factor(Peak_sets$Celltype,
                             levels = rev(sort(unique(Peak_sets$Celltype))))
Celltype_peakType <- as.data.frame(table(Peak_sets$Celltype, Peak_sets$peakType))
Celltype_peakType$Freq <- as.numeric(Celltype_peakType$Freq)
Celltype_peakType <- Celltype_peakType[order(Celltype_peakType$Var1),]
for (i in unique(as.character(Celltype_peakType$Var1))) {
  cell_num <- sum(Celltype_peakType[which(Celltype_peakType$Var1 == i),3])
  log2_num <- log2(cell_num)
  peaktype <- Celltype_peakType[which(Celltype_peakType$Var1 == i),3]
  peaktype <- peaktype / sum(peaktype)
  peaktype <- peaktype * log2_num
  Celltype_peakType[which(Celltype_peakType$Var1 == i),3] <- peaktype
}
Celltype_peakType_bar <- ggplot() +
  geom_bar(data = Celltype_peakType, aes(x = Var1, y = Freq, fill = Var2),
           stat = "identity", position = "stack") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  geom_text(data = as.data.frame(table(Peak_sets$Celltype)),
            aes(x = Var1, y = log2(Freq), label = Freq, hjust = 1.2)) +
  scale_fill_manual(values = c("#B18ED7","#ADC698","#A0C1D1","#D676A1")) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(x = "", y = "Log2(Count)", fill = "PeakType") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))



Peak_anno_match$Celltype <- factor(as.character(Peak_anno_match$Celltype),
                                   levels = rev(sort(unique(Peak_anno_match$Celltype))))
Peak_anno_match$Histone_modifications <- factor(as.character(Peak_anno_match$Histone_modifications),
                                                levels = c("H3K23ac_leaf", "H4K16ac_leaf", "H4K12ac_seedling", "H3K36me3_seedling", "H3K4me2_seedling"))
Celltype_peakMatch <- as.data.frame(table(Peak_anno_match$Celltype, Peak_anno_match$Histone_modifications))
Celltype_peakMatch$Freq <- as.numeric(Celltype_peakMatch$Freq)
Celltype_peakMatch <- Celltype_peakMatch[order(Celltype_peakMatch$Var1),]
Celltype_peakMatch$Freq2 <- Celltype_peakMatch$Freq
for (i in unique(as.character(Celltype_peakMatch$Var1))) {
  cell_num <- sum(Celltype_peakMatch[which(Celltype_peakMatch$Var1 == i),4])
  log2_num <- log2(cell_num)
  matchtype <- Celltype_peakMatch[which(Celltype_peakMatch$Var1 == i),4]
  matchtype <- matchtype / sum(matchtype)
  matchtype <- matchtype * log2_num
  Celltype_peakMatch[which(Celltype_peakMatch$Var1 == i),4] <- matchtype
}

Celltype_peakMatch$Freq3 <- Celltype_peakMatch$Freq2
for (cell in as.character(Celltype_peakMatch$Var1)) {
  temp_match <- Celltype_peakMatch[which(Celltype_peakMatch$Var1 == cell),]
  for (i in 1:nrow(temp_match)) {
    if (i == 1) {
      temp_match$Freq3[i] <- temp_match$Freq2[i] / 2
    } else {
      temp <- sum(temp_match$Freq2[1:(i-1)])
      temp_match$Freq3[i] <- sum(c(2*temp, temp_match$Freq2[i])) / 2
    }
  }
  Celltype_peakMatch[which(Celltype_peakMatch$Var1 == cell),] <- temp_match
}

Celltype_peakMatch$Var2 <- factor(as.character(Celltype_peakMatch$Var2),
                                  levels = c("H3K4me2_seedling",
                                             "H3K36me3_seedling",
                                             "H4K12ac_seedling",
                                             "H4K16ac_leaf",
                                             "H3K23ac_leaf"))
Celltype_peakMatch_bar <- ggplot() +
  geom_bar(data = Celltype_peakMatch, aes(x = Var1, y = Freq2, fill = Var2),
           stat = "identity", position = "stack") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  geom_text(data = Celltype_peakMatch,
            aes(x = Var1, y = Freq3, label = Freq), color = "black", alpha = 0.8) +
  scale_fill_manual(values = rev(c("#0ebeff","#4a9fec","#8780d9","#c361c6","#ff42b3"))) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(x = "", y = "Log2(Count)", fill = "Histone Modifications") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))




pdf("Figure2_peakSets.pdf", width = 24, height = 6)
multiplot(Celltype_tissue_bar, Celltype_peakType_bar, Celltype_peakMatch_bar, Peak_anno_match_percent, cols = 4)
dev.off()


##### marker peaks

Peak_sets$End <- Peak_sets$Start + 500
Peak_sets$Chr_start_end <- paste0(Peak_sets$Chr,":",Peak_sets$Start,"-",Peak_sets$End)
rownames(Peak_sets) <- Peak_sets$Chr_start_end
markerPeaksDF <- c()
for (i in names(markerPeaksList_Celltype@listData)) {
  temp <- as.data.frame(markerPeaksList_Celltype@listData[[i]])
  temp$Celltype <- i
  temp <- temp[order(temp$Log2FC, decreasing = T),]
  temp$Chr_start_end <- paste0(temp$seqnames,":",temp$start,"-",temp$end)
  temp$nearestGene <- Peak_sets[temp$Chr_start_end, "nearestGene"]
  markerPeaksDF <- as.data.frame(rbind(markerPeaksDF,
                                       temp))
}
sum(length(unique(markerPeaksDF$Chr_start_end)))
markerPeaksDF %>%
  group_by(Celltype) %>%
  top_n(n = 1, wt = Log2FC) -> top1
length(unique(top1$nearestGene))

celltype_heatmapPeakMatrix <- plotMarkerHeatmap(
  seMarker = markersPeaks_Celltype, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = F,
  plotLog2FC = T,
  labelMarkers = NULL,
  labelRows = FALSE,
  clusterCols = F,
  returnMatrix = T
)
sum(rownames(celltype_heatmapPeakMatrix) %in% Peak_sets$Chr_start_end)
celltype_heatmapPeakMatrix <- t(celltype_heatmapPeakMatrix)

peakmatrix <- getMatrixFromProject()
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
names(celltype_color) <- rownames(celltype_heatmapPeakMatrix)
row_anno <- rowAnnotation(df = data.frame(Celltype = rownames(celltype_heatmapPeakMatrix)), col = list(Celltype = celltype_color))
top1_peak <- HeatmapAnnotation(which = "column",
                               ano = anno_mark(at = match(top1$Chr_start_end,colnames(celltype_heatmapPeakMatrix)),
                                               labels = top1$nearestGene,
                                               labels_gp = gpar(fontsize = 12),
                                               link_width = unit(5, "mm"),
                                               extend = unit(0, "mm"),
                                               padding = unit(2, "mm")))
pdf("Figure2_celltype_peaks_heatmap.pdf", width = 14.5, height = 8)
Heatmap(celltype_heatmapPeakMatrix, col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        show_row_names = T, show_column_names = F, row_names_side = c("left"), name = paste0("Col Z-Scores","\n",ncol(celltype_heatmapPeakMatrix)," features","\n","Peak Matrix"),
        cluster_rows = F, cluster_columns = F, use_raster = T, raster_quality = 20,
        left_annotation = row_anno, top_annotation = top1_peak)
dev.off(Figure2_celltype_peaks_heatmap.pdf)

Celltype_peak_num <- as.data.frame(table(markerPeaksDF$Celltype))
Celltype_peak_num$Var1 <- factor(Celltype_peak_num$Var1,
                                 levels = rev(rownames(celltype_heatmapPeakMatrix)))
pdf("Figure2_Celltype_peak_num.pdf", height = 6, width = 5)
ggplot() +
  geom_bar(data = Celltype_peak_num, aes(x = Var1, y = log2(Freq), fill = Var1),
           stat = "identity") +
  scale_fill_manual(values = celltype_color) +
  geom_text(data = Celltype_peak_num, aes(x = Var1, y = log2(Freq), label = Freq),
            hjust = 1.2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "Log2(Count)") +
  coord_flip()
dev.off()


### 导入DAVID数据
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
}


###### GO
nrow(Peak_sets)
peakset <- getPeakSet(proj_pass_filter)
peakset@ranges@NAMES <- NULL
peakset <- as.data.frame(peakset)
rownames(peakset) <- paste0(peakset$seqnames, "_", peakset$start, "_", peakset$end)
{
  genes_list <- c()
  rownames(gene_id_map) <- gene_id_map$symbol
  for (i in names(markerPeaksList_Celltype)) {
    temp <- markerPeaksList_Celltype[[i]]
    temp <- as.data.frame(temp@listData)
    rownames(temp) <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
    temp <- peakset[rownames(temp), "nearestGene"]
    temp <- unique(temp)
    temp <- gene_id_map[temp, "gene_id"]
    genes_list <- c(genes_list, list(temp))
  }
  names(genes_list) <- names(markerPeaksList_Celltype)
  GO_result_celltype <- c()
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
                                                 Common_genes = paste(gene_count, collapse = "_"),
                                                 Gene_ratio = gene_ratio,
                                                 P_value = p))
      } else {
        next
      }
    } # GO term enrichment
    GO_result$Celltype <- i
    GO_result_celltype <- as.data.frame(rbind(GO_result_celltype, GO_result))
  }
  
  GO_result_celltype <- GO_result_celltype[which(GO_result_celltype$Ontology == "BP"),]
  GO_result_celltype <- GO_result_celltype[which(GO_result_celltype$P_value < 0.05),]
  
  GO_result_celltype_split <- split.data.frame(GO_result_celltype, f = list(GO_result_celltype$Celltype))
  for (i in names(GO_result_celltype_split)) {
    temp <- GO_result_celltype_split[[i]]
    temp <- temp[order(temp$P_value, decreasing = F),]
    GO_result_celltype_split[[i]] <- temp
    
  }
  
  names(GO_result_celltype_split)[1]
  temp <- GO_result_celltype_split[[1]]
  GO_BM <- c("developmental process", "anatomical structure development", "cell differentiation", "system development", "cell fate specification")
  GO_BM <- GO_result_celltype_split[[1]][which(GO_result_celltype_split[[1]]$Description %in% GO_BM), ]
  
  names(GO_result_celltype_split)[2]
  temp <- GO_result_celltype_split[[2]]
  GO_Cortex <- c("organic cyclic compound biosynthetic process", "regulation of RNA biosynthetic process", "cellular response to hormone stimulus", "regulation of biosynthetic process", "secondary metabolic process")
  GO_Cortex <- GO_result_celltype_split[[2]][which(GO_result_celltype_split[[2]]$Description %in% GO_Cortex), ]
  
  names(GO_result_celltype_split)[3]
  temp <- GO_result_celltype_split[[3]]
  GO_CBB <- c("developmental process", "anatomical structure development", "system development", "multicellular organism development", "developmental process involved in reproduction")
  GO_CBB <- GO_result_celltype_split[[3]][which(GO_result_celltype_split[[3]]$Description %in% GO_CBB), ]
  
  
  names(GO_result_celltype_split)[4]
  temp <- GO_result_celltype_split[[4]]
  GO_Endodermis <- c("cellular response to hormone stimulus", "cellular response to endogenous stimulus", "cell wall organization or biogenesis", "cellular response to organic substance", "regulation of macromolecule biosynthetic process")
  GO_Endodermis <- GO_result_celltype_split[[4]][which(GO_result_celltype_split[[4]]$Description %in% GO_Endodermis), ]
  
  
  names(GO_result_celltype_split)[5]
  temp <- GO_result_celltype_split[[5]]
  GO_Epi_ini <- c("response to light stimulus", "response to radiation", "regulation of stomatal movement", "cell differentiation", "cellular developmental process")
  GO_Epi_ini <- GO_result_celltype_split[[5]][which(GO_result_celltype_split[[5]]$Description %in% GO_Epi_ini), ]
  
  names(GO_result_celltype_split)[6]
  temp <- GO_result_celltype_split[[6]]
  GO_Epi <- c("ethylene-activated signaling pathway", "regulation of response to stimulus", "intracellular signal transduction", "cellular response to chemical stimulus", "signal transduction")
  GO_Epi <- GO_result_celltype_split[[6]][which(GO_result_celltype_split[[6]]$Description %in% GO_Epi), ]
  
  
  names(GO_result_celltype_split)[7]
  temp <- GO_result_celltype_split[[7]]
  GO_Fiber <- c("signal transduction", "response to stimulus", "response to auxin", "response to wounding", "cellular response to chemical stimulus")
  GO_Fiber <- GO_result_celltype_split[[7]][which(GO_result_celltype_split[[7]]$Description %in% GO_Fiber), ]
  
  
  names(GO_result_celltype_split)[8]
  temp <- GO_result_celltype_split[[8]]
  GO_IM <- c("developmental process", "system development", "cellular developmental process", "cell differentiation", "plant organ formation")
  GO_IM <- GO_result_celltype_split[[8]][which(GO_result_celltype_split[[8]]$Description %in% GO_IM), ]
  
  
  names(GO_result_celltype_split)[9]
  temp <- GO_result_celltype_split[[9]]
  GO_Lar_par_MO <- c("cell communication","response to water deprivation", "carbon fixation", "cellular response to light stimulus", "carbohydrate metabolic process")
  GO_Lar_par_MO <- GO_result_celltype_split[[9]][which(GO_result_celltype_split[[9]]$Description %in% GO_Lar_par_MO), ]
  
  
  names(GO_result_celltype_split)[10]
  temp <- GO_result_celltype_split[[10]]
  GO_Lar_par_PO <- c("cell communication","response to water deprivation", "carbon fixation", "response to light stimulus",
                     "carbohydrate metabolic process", "cell surface receptor signaling pathway", "response to hormone")
  GO_Lar_par_PO <- GO_result_celltype_split[[10]][which(GO_result_celltype_split[[10]]$Description %in% GO_Lar_par_PO), ]
  
  
  
  names(GO_result_celltype_split)[11]
  temp <- GO_result_celltype_split[[11]]
  GO_Mesophyll_MO <- c("response to red or far red light", "photosynthesis", "response to light stimulus",
                       "carbon fixation", "cellular response to light stimulus")
  GO_Mesophyll_MO <- GO_result_celltype_split[[11]][which(GO_result_celltype_split[[11]]$Description %in% GO_Mesophyll_MO), ]
  
  
  
  names(GO_result_celltype_split)[12]
  temp <- GO_result_celltype_split[[12]]
  GO_Mesophyll_PO <- c("photosynthesis", "response to light stimulus", "cellular carbohydrate biosynthetic process",
                       "photosynthesis, light harvesting", "cell fate commitment")
  GO_Mesophyll_PO <- GO_result_celltype_split[[12]][which(GO_result_celltype_split[[12]]$Description %in% GO_Mesophyll_PO), ]
  
  
  
  names(GO_result_celltype_split)[13]
  temp <- GO_result_celltype_split[[13]]
  GO_Mesophyll_init <- c("cell cycle", "mitotic cell cycle", "photosynthesis", "cell division", "developmental process")
  GO_Mesophyll_init <- GO_result_celltype_split[[13]][which(GO_result_celltype_split[[13]]$Description %in% GO_Mesophyll_init), ]
  
  
  
  names(GO_result_celltype_split)[14]
  temp <- GO_result_celltype_split[[14]]
  GO_Mesophyll_pre <- c("photosynthesis", "developmental process", "anatomical structure development", "cell morphogenesis", "multidimensional cell growth")
  GO_Mesophyll_pre <- GO_result_celltype_split[[14]][which(GO_result_celltype_split[[14]]$Description %in% GO_Mesophyll_pre), ]
  
  
  names(GO_result_celltype_split)[15]
  temp <- GO_result_celltype_split[[15]]
  GO_Mestome_sheath <- c("response to wounding", "response to water", "small molecule metabolic process",
                         "intracellular signal transduction", "signal transduction")
  GO_Mestome_sheath <- GO_result_celltype_split[[15]][which(GO_result_celltype_split[[15]]$Description %in% GO_Mestome_sheath), ]
  
  
  names(GO_result_celltype_split)[16]
  temp <- GO_result_celltype_split[[16]]
  GO_Pericycle <- c("plant organ morphogenesis", "plant organ formation", "asymmetric cell division",
                    "plant-type cell wall organization or biogenesis", "meristem determinacy")
  GO_Pericycle <- GO_result_celltype_split[[16]][which(GO_result_celltype_split[[16]]$Description %in% GO_Pericycle), ]
  
  
  names(GO_result_celltype_split)[17]
  temp <- GO_result_celltype_split[[17]]
  GO_Phloem <- c("anion transport", "amino acid transport", "carboxylic acid transport",
                 "amine transport", "regulation of transport")
  GO_Phloem <- GO_result_celltype_split[[17]][which(GO_result_celltype_split[[17]]$Description %in% GO_Phloem), ]
  
  
  names(GO_result_celltype_split)[18]
  temp <- GO_result_celltype_split[[18]]
  GO_Procambium <- c("organ growth", "cell communication", "response to hormone",
                     "response to auxin", "response to ethylene")
  GO_Procambium <- GO_result_celltype_split[[18]][which(GO_result_celltype_split[[18]]$Description %in% GO_Procambium), ]
  
  
  names(GO_result_celltype_split)[19]
  temp <- GO_result_celltype_split[[19]]
  GO_Rachis <- c("developmental process", "cell differentiation", "system development",
                     "cellular developmental process", "cell fate specification")
  GO_Rachis <- GO_result_celltype_split[[19]][which(GO_result_celltype_split[[19]]$Description %in% GO_Rachis), ]
  
  
  names(GO_result_celltype_split)[20]
  temp <- GO_result_celltype_split[[20]]
  GO_Root_hair <- c("response to ethylene", "ion transport", "response to acidic pH", "response to water", "organic acid transport")
  GO_Root_hair <- GO_result_celltype_split[[20]][which(GO_result_celltype_split[[20]]$Description %in% GO_Root_hair), ]
  
  
  names(GO_result_celltype_split)[21]
  temp <- GO_result_celltype_split[[21]]
  GO_SM <- c("cellular developmental process", "tissue development", "cell differentiation", "meristem initiation",
             "cell cycle process")
  GO_SM <- GO_result_celltype_split[[21]][which(GO_result_celltype_split[[21]]$Description %in% GO_SM), ]
  
  
  names(GO_result_celltype_split)[22]
  temp <- GO_result_celltype_split[[22]]
  GO_Stele <- c("lateral root development", "cell wall macromolecule metabolic process", "positive regulation of cell growth",
                "ammonium transport", "nitrate transport")
  GO_Stele <- GO_result_celltype_split[[22]][which(GO_result_celltype_split[[22]]$Description %in% GO_Stele), ]
  
  
  GO_select_celltype <- as.data.frame(rbind(GO_BM, GO_Cortex, GO_CBB, GO_Endodermis, GO_Epi_ini, GO_Epi, GO_Fiber, GO_IM, GO_Lar_par_MO,
                                            GO_Lar_par_PO, GO_Mesophyll_MO, GO_Mesophyll_PO, GO_Mesophyll_init, GO_Mesophyll_pre,
                                            GO_Mestome_sheath, GO_Pericycle, GO_Phloem, GO_Procambium, GO_Rachis, GO_Root_hair,
                                            GO_SM, GO_Stele))
  GO_select_celltype$Description <- Hmisc::capitalize(GO_select_celltype$Description)
  GO_select_celltype$Description <- paste0(GO_select_celltype$Celltype, "_", GO_select_celltype$Description)
  GO_select_celltype$Description <- factor(GO_select_celltype$Description,
                                           levels = rev(GO_select_celltype$Description))
  pdf("Figure2_Celltype_GO.pdf", width = 14, height = 15)
  ggplot(data = GO_select_celltype, aes(x = -log10(P_value), y = Description, fill = Celltype)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_vline(xintercept = -log10(0.05), colour = "red") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.ticks = element_line(color = "black")) +
    facet_wrap(.~Celltype, ncol = 2, scales = "free") +
    scale_fill_manual(values = celltype_color) +
    labs(y = "Enriched GO terms (Biological process)", y = "-log10(P-value)") +
    NoLegend()
  dev.off()
  
  {
    temp <- as.data.frame.array(table(GO_result_celltype$Description, GO_result_celltype$Celltype))
    GO_select <- names(which(rowSums(temp) == 1))
    temp <- as.data.frame(temp)
    temp <- temp[GO_select,]
    GO_result_celltype_temp <- GO_result_celltype[which(GO_result_celltype$Description %in% GO_select),]
    GO_result_celltype_temp <- split.data.frame(GO_result_celltype_temp, f = list(GO_result_celltype_temp$Celltype))
    
    GO_result_celltype_matrix <- matrix(NA, nrow = 19*5+3+4, ncol = 21)
    colnames(GO_result_celltype_matrix) <- names(GO_result_celltype_temp)
    GO_order <- c()
    for (i in colnames(GO_result_celltype_matrix)) {
      temp <- GO_result_celltype_temp[[i]]
      temp <- temp[order(temp$P_value),]
      if (nrow(temp) >= 5) {
        temp <- temp[1:5,]
      }
      GO_order <- c(GO_order, temp$Description)
    }
    rownames(GO_result_celltype_matrix) <- GO_order
    
    rownames(GO_result_celltype) <- paste0(GO_result_celltype$Celltype, "_", GO_result_celltype$Description)
    for (i in colnames(GO_result_celltype_matrix)) {
      for (j in rownames(GO_result_celltype_matrix)) {
        GO_result_celltype_matrix[j, i] <- GO_result_celltype[paste0(i, "_", j), "P_value"]
      }
    }
    
    GO_result_celltype_matrix <- GO_result_celltype_matrix[-6,]
    GO_result_celltype_matrix <- -log10(GO_result_celltype_matrix)
    rownames(GO_result_celltype_matrix) <- Hmisc::capitalize(rownames(GO_result_celltype_matrix))
    GO_result_celltype_matrix[which(GO_result_celltype_matrix > 3)] <- 3
    library(circlize)
    library(momr)
    col_fun <- colorRamp2(c(1.5, 2, 2.5, 3),
                          c("#F08080", "#FF4500", "#FF0000", "#8B0000"))
    col <- c(seq(1.5,3,0.05))
    color_bar <- col_fun(col)
    pdf("Figure2_GO_celltype.pdf", height = 20, width = 12)
    pheatmap(GO_result_celltype_matrix, cluster_rows = F, cluster_cols = F,
             color = color_bar, border_color = NA)
    dev.off()
  }
  
  
}







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

ATAC_common_tissues <- addGeneIntegrationMatrix(
  ArchRProj = ATAC_common_tissues, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = RNA_common_tissues,
  groupRNA = "Celltype",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  sampleCellsATAC = 15000,
  sampleCellsRNA = 15000,
  addToArrow = TRUE,
  force = TRUE
)
getAvailableMatrices(ATAC_common_tissues)
IntegrationData <- getMatrixFromProject(ATAC_common_tissues,
                                        useMatrix = "GeneIntegrationMatrix")
IntegrationData@assays@data@listData[["GeneIntegrationMatrix"]]@Dimnames[[1]] <- IntegrationData@elementMetadata@listData[["name"]]

common_genes <- intersect(rownames(RNA_common_tissues),
                          IntegrationData@elementMetadata@listData[["name"]])
ATAC <- CreateSeuratObject(IntegrationData@assays@data@listData[["GeneIntegrationMatrix"]],
                           meta.data = as.data.frame(IntegrationData@colData))
ATAC <- FindVariableFeatures(ATAC)
RNA <- RNA_common_tissues[common_genes,]
# table(ATAC_common_tissues$Celltype, ATAC_common_tissues$predictedGroup_Un)
RNA <- NormalizeData(RNA)
RNA <- FindVariableFeatures(RNA)

features <- SelectIntegrationFeatures(object.list = list(ATAC = ATAC, RNA = RNA))

RNA <- ScaleData(RNA, features = features, verbose = FALSE)
RNA <- RunPCA(RNA, features = features, verbose = FALSE)
ATAC <- ScaleData(ATAC, features = features, verbose = FALSE)
ATAC <- RunPCA(ATAC, features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(ATAC = ATAC, RNA = RNA),
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca")
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

integrated$Type <- "ATAC"
integrated@meta.data[colnames(RNA),"Type"] <- "RNA"
table(integrated$Type)

DimPlot(integrated, group.by = "Celltype", cols = celltype_color, split.by = "Type")
DimPlot(integrated, group.by = "TissuesSub", split.by = "Type")



########## change exon symbol
all_trans <- read.csv("./Ref/Oryza_sativa.IRGSP-1.0.52.chr.gtf",
                      sep = "\t", header = F)
trans <- unlist(lapply(all_trans$V9, function(x){
  temp <- unlist(strsplit(x, "transcript_id ", fixed = T))[2]
  unlist(strsplit(temp, "; ", fixed = T))[1]
}))
gene <- unlist(lapply(all_trans$V9, function(x){
  temp <- unlist(strsplit(x, "gene_id ", fixed = T))[2]
  unlist(strsplit(temp, ";", fixed = T))[1]
}))
all_trans <- data.frame(Trans_ID = trans,
                        Locus_ID = gene)
all_trans <- all_trans[!duplicated(all_trans),]
all_trans$Symbol <- gene_id_map[all_trans$Locus_ID, 2]
all_trans$Symbol[is.na(all_trans$Symbol)] <- all_trans$Locus_ID[is.na(all_trans$Symbol)]
rownames(all_trans) <- all_trans$Trans_ID
all_trans$Trans_ID_2 <- unlist(lapply(all_trans$Trans_ID, function(x){
  unlist(strsplit(x, "-", fixed = T))[1]
}))
all_trans <- all_trans[!duplicated(all_trans$Trans_ID_2),]
rownames(all_trans) <- all_trans$Trans_ID_2

exon_symbol <- all_trans[ATAC_common_tissues@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["symbol"]],"Symbol"]
exon_symbol[is.na(exon_symbol)] <- ATAC_common_tissues@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["symbol"]][is.na(exon_symbol)]
exon_id <- all_trans[ATAC_common_tissues@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]],"Locus_ID"]
exon_id[is.na(exon_id)] <- ATAC_common_tissues@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]][is.na(exon_id)]

ATAC_common_tissues@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]] <- exon_id
ATAC_common_tissues@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["symbol"]] <- exon_symbol

ATAC_common_tissues <- addPeak2GeneLinks(
  ArchRProj = ATAC_common_tissues,
  reducedDims = "IterativeLSI"
)
peak_to_gene <- getPeak2GeneLinks(
  ArchRProj = ATAC_common_tissues,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = F
)
idxATAC <- metadata(peak_to_gene)[[1]]
idxATAC <- as.data.frame(idxATAC)
idxRNA <- metadata(peak_to_gene)[[2]]
idxRNA <- as.data.frame(idxRNA)

peak_to_gene_DF <- as.data.frame(peak_to_gene)
sort(table(peak_to_gene_DF$idxRNA), decreasing = T)
peak_to_gene_DF[which(peak_to_gene_DF$idxRNA == "36271"),]$idxATAC

geneset <- as.data.frame(peak_to_gene@metadata[["geneSet"]]@elementMetadata@listData)
gene <- geneset[36271,]


p2g <- getPeak2GeneLinks(
  ArchRProj = ATAC_common_tissues_temp,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE
)

ATAC_common_tissues$Celltype
celltype_temp <- sort(unique(ATAC_common_tissues$Celltype))
celltype_temp <- celltype_temp[-13]
celltype_color_temp <- celltype_color[celltype_temp]

p <- plotBrowserTrack(
  ArchRProj = ATAC_common_tissues_temp, 
  groupBy = "Celltype", 
  geneSymbol = "Os09g0453100", 
  upstream = 250000,
  downstream = 250000,
  loops = p2g,
  pal = celltype_color_temp
)
pdf("Figure2_Os09g0453100_peak_to_gene.pdf", height = 6, width = 6)
grid::grid.newpage()
grid::grid.draw(p$Os09g0453100)
dev.off()

p <- plotPeak2GeneHeatmap(ArchRProj = ATAC_common_tissues, groupBy = "Celltype",
                          palGroup = celltype_color_temp, seed = 123, k = 5)
pdf("Figure2_peak_to_gene_linkage_2.pdf", width = 12, height = 5.5)
p
dev.off()

Os09g0453100_integrated_RNA <- data.frame(Os09g0453100 = IntegrationData@assays@data@listData[["GeneIntegrationMatrix"]]["Os09g0453100",],
                                          celltype = IntegrationData@colData@listData[["Celltype"]])

ggplot() +
  geom_violin(data = Os09g0453100_integrated_RNA,
              aes(x = celltype, y = Os09g0453100, fill = celltype)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

Peak2GeneHeatmapMatrix <- plotPeak2GeneHeatmap(ArchRProj = ATAC_common_tissues, groupBy = "Celltype",
                                               palGroup = celltype_color_temp, seed = 123, k = 5, returnMatrices = T)


# co-accessibility
all_trans <- read.csv("./Ref/Oryza_sativa.IRGSP-1.0.52.chr.gtf",
                      sep = "\t", header = F)
trans <- unlist(lapply(all_trans$V9, function(x){
  temp <- unlist(strsplit(x, "transcript_id ", fixed = T))[2]
  unlist(strsplit(temp, "; ", fixed = T))[1]
}))
gene <- unlist(lapply(all_trans$V9, function(x){
  temp <- unlist(strsplit(x, "gene_id ", fixed = T))[2]
  unlist(strsplit(temp, ";", fixed = T))[1]
}))
all_trans <- data.frame(Trans_ID = trans,
                        Locus_ID = gene)
all_trans <- all_trans[!duplicated(all_trans),]
all_trans$Symbol <- gene_id_map[all_trans$Locus_ID, 2]
all_trans$Symbol[is.na(all_trans$Symbol)] <- all_trans$Locus_ID[is.na(all_trans$Symbol)]
rownames(all_trans) <- all_trans$Trans_ID
all_trans$Trans_ID_2 <- unlist(lapply(all_trans$Trans_ID, function(x){
  unlist(strsplit(x, "-", fixed = T))[1]
}))
all_trans <- all_trans[!duplicated(all_trans$Trans_ID_2),]
rownames(all_trans) <- all_trans$Trans_ID_2

exon_symbol <- all_trans[proj_pass_filter@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["symbol"]],"Symbol"]
exon_symbol[is.na(exon_symbol)] <- proj_pass_filter@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["symbol"]][is.na(exon_symbol)]
exon_id <- all_trans[proj_pass_filter@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]],"Locus_ID"]
exon_id[is.na(exon_id)] <- proj_pass_filter@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]][is.na(exon_id)]

proj_pass_filter@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]] <- exon_id
proj_pass_filter@geneAnnotation@listData[["exons"]]@elementMetadata@listData[["symbol"]] <- exon_symbol


proj_pass_filter <- addCoAccessibility(
  ArchRProj = proj_pass_filter,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = proj_pass_filter,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA_temp1 <- data.frame(Start = cA@metadata[["peakSet"]]@ranges@start,
                       Celltype = cA@metadata[["peakSet"]]@ranges@NAMES,
                       Chr = as.character(cA@metadata[["peakSet"]]@seqnames))
cA_temp1$Chr_start_celltype <- paste0(cA_temp1$Chr, ":", cA_temp1$Start,
                                      "-", cA_temp1$Celltype)
sum(duplicated(cA_temp1$Chr_start))

cA_temp2 <- as.data.frame(cA@listData)
cA_temp2$queryHits_Chr_start_celltype <- cA_temp1[cA_temp2$queryHits,"Chr_start_celltype"]
cA_temp2$queryHits_celltype <- cA_temp1[cA_temp2$queryHits,"Celltype"]
cA_temp2$subjectHits_Chr_start_celltype <- cA_temp1[cA_temp2$subjectHits,"Chr_start_celltype"]
cA_temp2$subjectHits_celltype <- cA_temp1[cA_temp2$subjectHits,"Celltype"]

markerPeaksDF$Chr_start_celltype <- paste0(markerPeaksDF$seqnames,":",
                                           markerPeaksDF$start,"-",markerPeaksDF$Celltype)
rownames(markerPeaksDF) <- markerPeaksDF$Chr_start_celltype

Peak_sets$Chr_start_celltype <- paste0(Peak_sets$Chr,":",Peak_sets$Start,"-",Peak_sets$Celltype)
rownames(Peak_sets) <- Peak_sets$Chr_start_celltype

cA_temp2$queryHits_Type <- Peak_sets[cA_temp2$queryHits_Chr_start_celltype,"peakType"]
cA_temp2$subjectHits_Type <- Peak_sets[cA_temp2$subjectHits_Chr_start_celltype,"peakType"]
cA_temp2$queryHits_nearestGene <- Peak_sets[cA_temp2$queryHits_Chr_start_celltype,"nearestGene"]
cA_temp2$subjectHits_nearestGene <- Peak_sets[cA_temp2$subjectHits_Chr_start_celltype,"nearestGene"]

cA_temp2$queryHits_celltype_specific <- markerPeaksDF[cA_temp2$queryHits_Chr_start_celltype,"Celltype"]
cA_temp2$subjectHits_celltype_specific <- markerPeaksDF[cA_temp2$subjectHits_Chr_start_celltype,"Celltype"]

cA_temp2$queryHits_celltype_FC <- markerPeaksDF[cA_temp2$queryHits_Chr_start_celltype,"Log2FC"]
cA_temp2$subjectHits_celltype_FC <- markerPeaksDF[cA_temp2$subjectHits_Chr_start_celltype,"Log2FC"]

cA_temp2 <- as.data.frame(na.omit(cA_temp2))
keep_row <- (cA_temp2$queryHits_Type == "Promoter") + (cA_temp2$subjectHits_Type == "Promoter")
cA_temp2 <- cA_temp2[which(keep_row > 0),]

keep_row <- (cA_temp2$queryHits_celltype == cA_temp2$subjectHits_celltype) + 0
cA_temp2 <- cA_temp2[which(keep_row > 0),]


p <- plotBrowserTrack(
  ArchRProj = proj_pass_filter, 
  groupBy = "Celltype", 
  geneSymbol = "Os04g0438300", 
  upstream = 20000,
  downstream = 30000,
  loops = getCoAccessibility(proj_pass_filter),
  pal = celltype_color,
  baseSize = 7, facetbaseSize = 7
)

dev.off()
pdf("Figure2_co_access_linkage.pdf", width = 10, height = 7)
grid::grid.newpage()
grid::grid.draw(p$Os04g0438300)
dev.off()


### snRNA
RNA <- readRDS("./Rice-snRNA/all_data_annotated.rds")
Idents(RNA) <- RNA$Celltype
RNA_cell_type_marker <- FindAllMarkers(RNA, only.pos = T, logfc.threshold = 1.5)
table(RNA_cell_type_marker$cluster)

RNA_cell_type_marker$symbol <- gene_id_map[RNA_cell_type_marker$gene,"symbol"]
RNA_cell_type_marker_genes <- unique(RNA_cell_type_marker$symbol)
RNA_cell_type_marker_genes <- na.omit(RNA_cell_type_marker_genes)
RNA_cell_type_marker_genes_DF <- c()
for (i in RNA_cell_type_marker_genes) {
  celltype <- RNA_cell_type_marker[which(RNA_cell_type_marker$symbol == i),"cluster"]
  celltype <- as.character(celltype)
  pct1 <- RNA_cell_type_marker[which(RNA_cell_type_marker$symbol == i),"pct.1"]
  pct2 <- RNA_cell_type_marker[which(RNA_cell_type_marker$symbol == i),"pct.2"]
  avg_log2FC <- RNA_cell_type_marker[which(RNA_cell_type_marker$symbol == i),"avg_log2FC"]
  gene <- i
  temp <- data.frame(gene = gene,
                     pct1 = paste(pct1, collapse = " | "),
                     pct2 = paste(pct2, collapse = " | "),
                     avg_log2FC = paste(avg_log2FC, collapse = " | "),
                     celltype = paste(celltype, collapse = " | ")
                     )
  RNA_cell_type_marker_genes_DF <- as.data.frame(rbind(RNA_cell_type_marker_genes_DF,
                                                       temp))
}

rownames(RNA_cell_type_marker_genes_DF) <- RNA_cell_type_marker_genes_DF$gene

cA_temp2_RNA <- cA_temp2
#cA_temp2_RNA$subjectHits_RNA <- RNA_cell_type_marker_genes_DF[cA_temp2_RNA$subjectHits_nearestGene,"celltype"]
cA_temp2_RNA$queryHits_RNA <- RNA_cell_type_marker_genes_DF[cA_temp2_RNA$queryHits_nearestGene,"celltype"]
cA_temp2_RNA$queryHits_LogFC <- RNA_cell_type_marker_genes_DF[cA_temp2_RNA$queryHits_nearestGene,"avg_log2FC"]

cA_temp2_RNA <- cA_temp2_RNA[-which(is.na(cA_temp2_RNA$queryHits_RNA)),]
cA_temp2_RNA <- cA_temp2_RNA[which(cA_temp2_RNA$queryHits_celltype_specific == cA_temp2_RNA$queryHits_RNA),]
cA_temp2_RNA <- cA_temp2_RNA[which(cA_temp2_RNA$queryHits_nearestGene == cA_temp2_RNA$subjectHits_nearestGene),]

# OsLOX8

p <- plotBrowserTrack(
  ArchRProj = proj_pass_filter, 
  groupBy = "Celltype", 
  geneSymbol = "OsLOX8", 
  upstream = 20000,
  downstream = 20000,
  loops = getCoAccessibility(proj_pass_filter),
  pal = celltype_color
)

dev.off()
pdf("Figure2_co_access_linkage_2.pdf", width = 10, height = 7)
grid::grid.newpage()
grid::grid.draw(p$OsLOX8)
dev.off()

OsLOX8_integrated_RNA <- data.frame(OsLOX8 = IntegrationData@assays@data@listData[["GeneIntegrationMatrix"]]["OsLOX8",],
                                    celltype = IntegrationData@colData@listData[["Celltype"]])

ggplot() +
  geom_violin(data = OsLOX8_integrated_RNA,
              aes(x = celltype, y = OsLOX8, fill = celltype)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

id <- gene_id_map[which(gene_id_map$symbol == "OsLOX8"),][1,1]

RNA_OsLOX8 <- data.frame(OsLOX8 = RNA@assays$RNA@data[id,],
                         Celltype = RNA$Celltype)
DotPlot(RNA, group.by = "Celltype", features = id)


### integrated gene score
GeneScore_integrated <- CreateSeuratObject(IntegrationData@assays@data$GeneIntegrationMatrix)
GeneScore_integrated@meta.data <- as.data.frame(cbind(GeneScore_integrated@meta.data,
                                                      as.data.frame(IntegrationData@colData@listData)))
Idents(GeneScore_integrated) <- GeneScore_integrated$Celltype
GeneScore_integrated_cell_type_marker <- FindAllMarkers(GeneScore_integrated, only.pos = T, logfc.threshold = 1.5)
table(GeneScore_integrated_cell_type_marker$cluster)

GeneScore_integrated_cell_type_marker$symbol <- gene_id_map[GeneScore_integrated_cell_type_marker$gene,"symbol"]
GeneScore_integrated_cell_type_marker_genes <- unique(GeneScore_integrated_cell_type_marker$symbol)
GeneScore_integrated_cell_type_marker_genes <- na.omit(GeneScore_integrated_cell_type_marker_genes)
GeneScore_integrated_cell_type_marker_genes_DF <- c()
for (i in GeneScore_integrated_cell_type_marker_genes) {
  celltype <- GeneScore_integrated_cell_type_marker[which(GeneScore_integrated_cell_type_marker$symbol == i),"cluster"]
  celltype <- as.character(celltype)
  pct1 <- GeneScore_integrated_cell_type_marker[which(GeneScore_integrated_cell_type_marker$symbol == i),"pct.1"]
  pct2 <- GeneScore_integrated_cell_type_marker[which(GeneScore_integrated_cell_type_marker$symbol == i),"pct.2"]
  avg_log2FC <- GeneScore_integrated_cell_type_marker[which(GeneScore_integrated_cell_type_marker$symbol == i),"avg_log2FC"]
  gene <- i
  temp <- data.frame(gene = gene,
                     pct1 = paste(pct1, collapse = " | "),
                     pct2 = paste(pct2, collapse = " | "),
                     avg_log2FC = paste(avg_log2FC, collapse = " | "),
                     celltype = paste(celltype, collapse = " | ")
  )
  GeneScore_integrated_cell_type_marker_genes_DF <- as.data.frame(rbind(GeneScore_integrated_cell_type_marker_genes_DF,
                                                                        temp))
}

rownames(GeneScore_integrated_cell_type_marker_genes_DF) <- GeneScore_integrated_cell_type_marker_genes_DF$gene

cA_temp2_RNA <- cA_temp2
#cA_temp2_RNA$subjectHits_RNA <- RNA_cell_type_marker_genes_DF[cA_temp2_RNA$subjectHits_nearestGene,"celltype"]
cA_temp2_RNA$queryHits_RNA <- GeneScore_integrated_cell_type_marker_genes_DF[cA_temp2_RNA$queryHits_nearestGene,"celltype"]
cA_temp2_RNA$queryHits_LogFC <- GeneScore_integrated_cell_type_marker_genes_DF[cA_temp2_RNA$queryHits_nearestGene,"avg_log2FC"]
if ( sum(is.na(cA_temp2_RNA$queryHits_RNA)) > 0 ) {
  cA_temp2_RNA <- cA_temp2_RNA[-which(is.na(cA_temp2_RNA$queryHits_RNA)),]
}
cA_temp2_RNA <- cA_temp2_RNA[which(cA_temp2_RNA$queryHits_celltype_specific == cA_temp2_RNA$queryHits_RNA),]
cA_temp2_RNA <- cA_temp2_RNA[which(cA_temp2_RNA$queryHits_nearestGene == cA_temp2_RNA$subjectHits_nearestGene),]

# Os01g0507700
# Os12g0233100
RNA_ATAC_common_celltype <- intersect(unique(RNA$Celltype), unique(proj_pass_filter$Celltype))
p <- plotBrowserTrack(
  ArchRProj = proj_pass_filter, 
  groupBy = "Celltype", 
  # geneSymbol = "Os12g0233100",
  geneSymbol = "OSH6", 
  upstream = 20000,
  downstream = 20000,
  loops = getCoAccessibility(proj_pass_filter),
  groupPal = celltype_color
)
grid::grid.newpage()
grid::grid.draw(p$OSH6)

DotPlot(GeneScore_integrated, features = "OSH6", group.by = "Celltype")
RNA_sub <- subset(RNA, subset = Celltype %in% common_celltypes)
DotPlot(RNA_sub, features = "Os01g0302500", group.by = "Celltype")
pdf("Figure2_RNA_OSH6.pdf", width = 5, height = 8)
DotPlot(RNA, features = "Os01g0302500", group.by = "Celltype")
dev.off()

dev.off()
pdf("Figure2_co_access_linkage_OSH6.pdf", width = 10, height = 9)
grid::grid.newpage()
# grid::grid.draw(p$OsFCA)
grid::grid.draw(p$OSH6)
dev.off()

GeneScore_integrated$Celltype <- factor(GeneScore_integrated$Celltype,
                                        levels = sort(unique(GeneScore_integrated$Celltype), decreasing = T))
pdf("Figure2_Os12g0233100_expression_level.pdf")
DotPlot(GeneScore_integrated, features = "Os12g0233100", group.by = "Celltype")
dev.off()

id <- gene_id_map[which(gene_id_map$symbol == "Os12g0233100"),][1,1]

pdf("Figure2_Os12g0233100_expression_level.pdf")
DotPlot(GeneScore_integrated, group.by = "Celltype", features = id)
dev.off()


saveRDS(proj_pass_filter, "3.proj_pass_filter_withPeaks.rds")



#### random peaks
proj_pass_filter <- readRDS("3.proj_pass_filter_withPeaks.rds")
peaksets <- getPeakSet(proj_pass_filter)
peaksets$Celltype <- peaksets@ranges@NAMES
peaksets@ranges@NAMES <- NULL
peaksets_df <- as.data.frame(peaksets)
mean(table(peaksets_df$Celltype)) # 6000

set.seed(123)
library(GenomicRanges)
chr_size <- genome_annotation@listData[["chromSizes"]]@seqinfo@seqlengths[1:12]
names(chr_size) <- genome_annotation@listData[["chromSizes"]]@seqinfo@seqnames[1:12]

chr_random <- sample(c(1:12), 6000, replace = T,
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

random_range <- makeGRangesFromDataFrame(random_range)
random_range_temp <- as.data.frame(random_range)
random_peaks <- paste0(random_range_temp$seqnames, "_", random_range_temp$start, "_", random_range_temp$end)
rm(random_range_temp)
random_range@elementMetadata@listData$PeakID <- random_peaks
saveRDS(random_range, "Celltype_random_range.rds")

regions <- list("H3K23ac_leaf" = H3K23ac_Osj_leaf_normal,
                     "H4K16ac_leaf" = H4K16ac_Osj_leaf_normal,
                     "H3K36me3_seedling" = H3K36me3_Osj_seedling_normal,
                     "H3K4me2_seedling" = H3K4me2_Osj_seedling_normal,
                     "H4K12ac_seedling" = H4K12ac_Osj_seedling_normal,
                     "CE" = CE)
regionPositions <- lapply(seq_along(regions), function(x){
  regions[[x]]
}) %>% GRangesList

names(regionPositions) <- names(regions)
allPositions <- unlist(regionPositions, use.names=TRUE)
overlapRegions <- findOverlaps(random_range, allPositions, ignore.strand=TRUE)
regionMat <- Matrix::sparseMatrix(
  i = queryHits(overlapRegions),
  j = match(names(allPositions),names(regionPositions))[subjectHits(overlapRegions)],
  x = rep(TRUE, length(overlapRegions)),
  dims = c(length(random_range), length(regionPositions))
)
colnames(regionMat) <- names(regionPositions)
regionMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = regionMat), rowRanges = random_range)

nO <- Matrix::colSums(assay(regionMat))
rF <- names(which(nO == 0))
if(length(rF) > 0){
  regionPositions <- regionPositions[!(names(regionPositions) %in% rF)]
  regionMat <- regionMat[,names(regionPositions),drop=FALSE]
}

temp <- regionMat@assays@data@listData[["matches"]]
temp <- as.data.frame(temp)
sum(rowSums(temp[,1:5]) > 0) / 6000
sum(temp[,6] > 0) / 6000
colSums(temp) / 6000
