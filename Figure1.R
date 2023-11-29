setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas")
library(ArchR)
library(ggalt)
library(Seurat)
library(BSgenome.OSativa.NCBI.IRGSPv1.0)
library(parallel)
library(clustree)
library(dplyr)
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

sample_library <- read.table("样本情况.csv", sep = ",", header = T)
sample_library$Path <- paste0("./fragments/",sample_library$Library,".fragments.tsv.gz")
files <- sample_library$Path
names(files) <- sample_library$Library

arrowFiles <- createArrowFiles(inputFiles = files,
                               sampleNames = names(files),
                               outputNames = names(files),
                               minTSS = 0, minFrags = 0,
                               excludeChr = c("Pt", "Mt"),
                               addTileMat = T,
                               addGeneScoreMat = T,
                               subThreading = F,
                               geneAnnotation = gene_annotation,
                               genomeAnnotation = genome_annotation)

arrowFiles <- list.files(pattern = ".arrow", recursive = F, full.names = F)

proj_all <- ArchRProject(
  ArrowFiles = arrowFiles, 
  geneAnnotation = gene_annotation,
  genomeAnnotation = genome_annotation,
  outputDirectory = "Rice-snATAC",
  copyArrows = F #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveRDS(proj_all, "1.proj_all_before_any_filter.rds")

rownames(sample_library) <- sample_library$Library

proj_all$Tissues <- sample_library[proj_all$Sample, "Tissue"]
proj_all$TissuesSub <- sample_library[proj_all$Sample, "TissueSub"]

samples <- unique(proj_all$Sample)
Tissues <- unique(proj_all$Tissues)
TissueSub <- unique(proj_all$TissueSub)

tissues_pro <- c()
for (t in Tissues) {
  print(t)
  cellnames <- getCellNames(proj_all)[which(proj_all$Tissues == t)]
  temp <- subsetArchRProject(ArchRProj = proj_all,
                             cells = cellnames,
                             outputDirectory = paste0("./tissues_project/",t),
                             dropCells = TRUE)
  tissues_pro <- c(tissues_pro, list(temp))
}
rm(temp)
names(tissues_pro) <- Tissues

# TSS&Fragament QC Plots
{
  filterFrags <- 1500
  filterTSS <- 2
  tissue_samples <- data.frame(sample = samples,
                               tissue = c("10dRoot", "10dRoot", "10dLeaf", "10dLeaf",
                                          "10dShoot", "10dShoot", "60dRootTip", "60dRootTip",
                                          "IM0.5cm", "IM0.5cm", "IM1cm", "IM1cm"))
  names(arrowFiles) <- unlist(lapply(arrowFiles, function(x){
    unlist(strsplit(x, ".arrow", fixed = T))[1]
  }))
  QCDir <- "./tissues_project/1.QCplot/"
  for (i in unique(tissue_samples$tissue)) {
    print(i)
    Metadata_all <- c()
    for (j in tissue_samples$sample[which(tissue_samples$tissue == i)]) {
      fragSummary <- .fastFragmentInfo(
        ArrowFile = arrowFiles[j], 
        cellNames = .availableCells(arrowFiles[j]), 
        nucLength = 147,
        logFile = NULL
      )
      Metadata <- fragSummary[[1]]
      TSSParams <- list()
      TSSParams$TSS <- gene_annotation$TSS
      TSSParams$ArrowFile <- arrowFiles[j]
      TSSParams$cellNames <- Metadata$cellNames
      TSSOut <- do.call(.fastTSSEnrichment, TSSParams)
      Metadata$TSSEnrichment <- TSSOut$tssScores
      Metadata$ReadsInTSS <- TSSOut$tssReads
      Metadata$Keep <- 1*(Metadata$nFrags >= filterFrags & Metadata$TSSEnrichment >= filterTSS)
      Metadata <- as.data.frame(Metadata)
      Metadata_all <- as.data.frame(rbind(Metadata_all, Metadata))
    }
    
    ggtitle <- sprintf("%s\n%s\n%s",
                       paste0(i, "\nnCells Pass Filter = ", sum(Metadata_all$Keep)),
                       paste0("Median Frags = ", median(Metadata_all$nFrags[Metadata_all$Keep==1])),
                       paste0("Median TSS Enrichment = ", median(Metadata_all$TSSEnrichment[Metadata_all$Keep==1]))
    )
    gg <- ggPoint(
      x = pmin(log10(Metadata_all$nFrags), 5) + rnorm(length(Metadata_all$nFrags), sd = 0.00001),
      y = Metadata_all$TSSEnrichment + rnorm(length(Metadata_all$nFrags), sd = 0.00001), 
      colorDensity = TRUE,
      xlim = c(2.5, 5),
      ylim = c(0, max(Metadata_all$TSSEnrichment) * 1.05),
      baseSize = 6,
      continuousSet = "sambaNight",
      xlabel = "Log 10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      title = ggtitle,
      rastr = TRUE) + 
      geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
      geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)
    pdf(file.path(QCDir,paste0(i,"-TSS_by_Unique_Frags.pdf")),width=6,height=6,onefile=FALSE)
    .fixPlotSize(gg, plotWidth = 6, plotHeight = 6)
    dev.off()
  }
  
}

for (t in Tissues) {
  temp_pro <- tissues_pro[[t]]
  pass_filter <- temp_pro$cellNames[c(temp_pro$TSSEnrichment >= filterTSS & temp_pro$nFrags >= filterFrags)]
  temp_pro <- temp_pro[pass_filter,]
  tissues_pro[[t]] <- temp_pro
}


### Doublet Identification
samples <- unique(proj_all$Sample)
samples_pro <- list()
for (t in Tissues) {
  samples_t <- tissues_pro[[t]]@sampleColData@rownames
  temp_pro <- tissues_pro[[t]]
  for (s in samples_t) {
    cellnames <- temp_pro$cellNames[c(temp_pro$Sample == s)]
    temp_pro_sample <- temp_pro[cellnames,]
    temp_list <- list(temp_pro_sample)
    names(temp_list) <- s
    samples_pro <- c(samples_pro, temp_list)
  }
}

for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  temp_pro_sample <- addDoubletScores(
    input = temp_pro_sample,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    LSIMethod = 1
  )
  samples_pro[[s]] <- temp_pro_sample
}

samples_filter_doublet <- samples
for (s in samples_filter_doublet) {
  temp_pro_sample <- samples_pro[[s]]
  temp_pro_sample <- filterDoublets(temp_pro_sample)
  samples_pro[[s]] <- temp_pro_sample
}

####

####
cells_pass_filter <- c()
for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  cells_pass_filter <- c(cells_pass_filter, temp_pro_sample$cellNames)
}
proj_pass_filter <- proj_all[cells_pass_filter,]
proj_pass_filter <- addIterativeLSI(ArchRProj = proj_pass_filter, useMatrix = "TileMatrix", name = "IterativeLSI")
proj_pass_filter <- addHarmony(
  ArchRProj = proj_pass_filter,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

pathToMacs2 <- "/home/heshidian/mambaforge/envs/common/bin/macs2"
addArchRGenome("genome_annotation")

proj_pass_filter <- addGroupCoverages(ArchRProj = proj_pass_filter, groupBy = "TissuesSub", force = TRUE)
proj_pass_filter <- addReproduciblePeakSet( # each peak is 501 bp in length
  ArchRProj = proj_pass_filter, 
  groupBy = "TissuesSub", 
  pathToMacs2 = pathToMacs2,
  excludeChr = c("Mt", "Pt"),
  force = TRUE,
  geneAnnotation = gene_annotation,
  genomeAnnotation = genome_annotation,
  genomeSize = sum(genome_annotation@listData[["chromSizes"]]@ranges@width)
)
proj_pass_filter <- addPeakMatrix(proj_pass_filter)
Tissue_peak_set <- getPeakSet(proj_pass_filter) # each peak is 501 bp in length
getAvailableMatrices(proj_pass_filter)
markersPeaks_tissues <- getMarkerFeatures(
  ArchRProj = proj_pass_filter, 
  useMatrix = "PeakMatrix", 
  groupBy = "TissuesSub",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerPeaksList_tissues <- getMarkers(markersPeaks_tissues, cutOff = "FDR <= 0.05 & Log2FC >= 1")
rownames(Peak_sets) <- paste0(Peak_sets$Chr, "_", Peak_sets$Start, "_", Peak_sets$End)
for (i in names(markerPeaksList_tissues)) {
  temp <- markerPeaksList_tissues@listData[[i]]
  temp <- as.data.frame(temp)
  temp$peakid <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
  temp$peakgroup <- Peak_sets[temp$peakid, "peakType"]
  temp <- table(temp$peakgroup) / nrow(temp)
  print(sum(temp[which(names(temp) %in% c("Distal", "Promoter"))]))
}


# Tissue_peak_set_meta <- as.data.frame(Tissue_peak_set@elementMetadata@listData)
Tissue_peak_set_meta <- as.data.frame(Tissue_peak_set, row.names = 1:132008)
Tissue_peak_set_meta$peak_id <- paste0("Chr", Tissue_peak_set_meta$seqnames, "_",
                                       Tissue_peak_set_meta$idx)
rownames(Tissue_peak_set_meta) <- Tissue_peak_set_meta$peak_id

# each peak is 501 bp in length
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks_tissues, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = F,
  plotLog2FC = T,
  labelMarkers = NULL,
  labelRows = FALSE,
  returnMatrix = T,
  clusterCols = F,
  nLabel = 1,
  nPrint = 0
)
pal <- paletteContinuous(set = "solarExtra", n = 100)

Heatmap(heatmapPeaks, col = pal, name = paste0("Row Z-Scores\n", nrow(heatmapPeaks), " features\n"),
        show_row_names = F, show_column_names = T, show_column_dend = F,
        cluster_rows = F, cluster_columns = T, use_raster = T, raster_resize_mat = median)

pdf("Figure1_tissue_peaks_heatmap.pdf", width = 4.5, height = 8)
# heatmap(heatmapPeaks, Rowv = F, Colv = F, labRow = F, scale = "none")
# draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
Heatmap(heatmapPeaks, col = pal, name = paste0("Row Z-Scores\n", nrow(heatmapPeaks), " features\n"),
        show_row_names = F, show_column_names = T, show_column_dend = F,
        cluster_rows = F, cluster_columns = T, use_raster = T, raster_resize_mat = median)
dev.off()

tissues_sub <- names(markerPeaksList_tissues@listData)
gene_id_map <- as.data.frame(gene_annotation@listData[["genes"]]@elementMetadata)
rownames(gene_id_map) <- gene_id_map$symbol


proj_pass_filter <- addUMAP(
  ArchRProj = proj_pass_filter, 
  reducedDims = "IterativeLSI", 
  name = "UMAP_withBatchEffect", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
proj_pass_filter <- addUMAP(
  ArchRProj = proj_pass_filter, 
  reducedDims = "Harmony", 
  name = "UMAP_Harmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

proj_pass_filter <- addTSNE(
  ArchRProj = proj_pass_filter, 
  reducedDims = "IterativeLSI", 
  name = "TSNE_withBatchEffect", 
  perplexity = 30
)
proj_pass_filter <- addTSNE(
  ArchRProj = proj_pass_filter, 
  reducedDims = "Harmony", 
  name = "TSNE_Harmony", 
  perplexity = 30
)

TSNE_withBatchEffect <- getEmbedding(proj_pass_filter, embedding = "TSNE_withBatchEffect")
colnames(TSNE_withBatchEffect) <- c("tSNE_1", "tSNE_2")
ColData <- as.data.frame(getCellColData(proj_pass_filter))
TSNE_withBatchEffect <- as.data.frame(cbind(TSNE_withBatchEffect,
                                            ColData))
pdf("Figure1_ATAC_tissues_tSNE.pdf", width = 6.8, height = 6)
ggplot(data = TSNE_withBatchEffect[which(TSNE_withBatchEffect$TissuesSub == "10dLeaf"),], aes(x = tSNE_1, y = tSNE_2, color = TissuesSub)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#488E44","#867C58","#618694","#A2A882",
                                "#eea29a","#f7786b","#8A9FD1","#C06CAB",
                                "#bd5734")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
dev.off()

write.csv(gene_id_map, "gene_id_map.csv", quote = F)
### integration
{
  ### integration
  RNA <- readRDS("./Rice-snRNA/all_data.rds")
  colnames(colData(RNA))
  # Unconstrained integration
  cells_select <- as.data.frame(getCellColData(proj_pass_filter))
  cells_select <- cells_select[which(cells_select$TissuesSub %in% c("IM0.5cm",
                                                                    "IM1cm",
                                                                    "10dLeaf",
                                                                    "10dRoot",
                                                                    "60dRootTip")),]
  tissue_sub_main <- data.frame(sub = c("IM0.5cm","IM1cm","10dLeaf",
                                        "10dRoot","60dRootTip","120dFlagLeaf","90dFlagLeaf"),
                                main = c("IM","IM","Leaf","Root","RootTip","FlagLeaf","FlagLeaf"),
                                row.names = c("IM0.5cm","IM1cm","10dLeaf",
                                              "10dRoot","60dRootTip","120dFlagLeaf","90dFlagLeaf"))
  proj_pass_filter_select <- proj_pass_filter[rownames(cells_select),]
  proj_pass_filter_select$main <- tissue_sub_main[proj_pass_filter_select$TissuesSub, "main"]
  
  Idents(RNA) <- RNA$tissues
  RNA_select <- subset(RNA, idents = c("10dLeaf","10dRoot","60dRootTip", # "120dFlagLeaf", "90dFlagLeaf"
                                       "IM0.5cm","IM1cm"))
  RNA_select$main <- tissue_sub_main[RNA_select$tissues, "main"]
  RNA_select <- RNA_select[which(rownames(RNA_select) %in% gene_id_map$gene_id),]
  rownames(gene_id_map) <- gene_id_map$gene_id
  RNA_id_gene <- data.frame(id = rownames(RNA_select),
                            symbol = gene_id_map[rownames(RNA_select),2],
                            row.names = rownames(RNA_select))
  RNA_select@assays$RNA@counts@Dimnames[[1]] <- RNA_id_gene$symbol
  RNA_select@assays$RNA@data@Dimnames[[1]] <- RNA_id_gene$symbol
  rownames(RNA_select@assays$RNA@scale.data) <- RNA_id_gene[rownames(RNA_select@assays$RNA@scale.data),"symbol"]
  RNA_select <- CreateSeuratObject(counts = RNA_select@assays$RNA@counts,
                                   meta.data = RNA_select@meta.data)
  RNA_select = NormalizeData(object = RNA_select, verbose = FALSE)
  proj_pass_filter_select <- addGeneIntegrationMatrix(
    ArchRProj = proj_pass_filter_select, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = RNA_select,
    addToArrow = FALSE,
    groupRNA = "main",
    groupATAC = "main",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    sampleCellsATAC = 15000,
    sampleCellsRNA = 15000
  )
  table(proj_pass_filter_select$main, proj_pass_filter_select$predictedGroup_Un)
  (18394+16800+10540+4963+12216) / length(proj_pass_filter_select$main) # 0.8938
  library(ggalluvial)
  library(ggplot2)
  library(dplyr)
  library(networkD3)
  library(riverplot)
  proj_pass_filter_select$main2 <- proj_pass_filter_select$main
  proj_pass_filter_select$main2[which(proj_pass_filter_select$main %in% c("Root", "RootTip"))] <- "Root"
  proj_pass_filter_select$predictedGroup_Un2 <- proj_pass_filter_select$predictedGroup_Un
  proj_pass_filter_select$predictedGroup_Un2[which(proj_pass_filter_select$predictedGroup_Un %in% c("Root", "RootTip"))] <- "Root"
  meta <- as.data.frame(table(paste0("ATAC-",proj_pass_filter_select$main2), paste0("RNA-",proj_pass_filter_select$predictedGroup_Un2)))
  colnames(meta) <- c("source","target","value")
  df.nodes <- data.frame(
    name=c(as.character(meta$source), 
           as.character(meta$target)) %>% unique()
  )
  meta$IDsource <- match(meta$source, df.nodes$name)-1 
  meta$IDtarget <- match(meta$target, df.nodes$name)-1
  meta$group <- "flow"
  library(scales)
  my_color <- 'd3.scaleOrdinal() 
            .domain(["0","1","2","3","4","5","6","7","ATAC-IM","ATAC-Leaf","ATAC-Root","ATAC-RootTip","RNA-IM","RNA-Leaf","RNA-Root","RNA-RootTip","flow"]) 
            .range(["#F8766D","#E9842C","#D69100","#BC9D00","#9CA700","#6FB000","#00B813","#00BD61","#DE6C91","#488E44","#867C58","#A2A882","#DE6C91","#488E44","#867C58","#A2A882","#DEDEE9"])'
  
  Sankey.p <- sankeyNetwork(Links = meta, Nodes = df.nodes,
                            Source = "IDsource", Target = "IDtarget",
                            Value = "value", NodeID = "name", LinkGroup = "group", colourScale=my_color, 
                            sinksRight=T, nodeWidth=30, nodePadding=6, fontSize=16,width=300,height = 800)
  Sankey.p
  library(htmlwidgets)
  saveWidget(Sankey.p, file="Figure1_Sankey.html")
  library(webshot)
# Sys.setenv(OPENSSL_CONF="/dev/null")
  webshot("Figure1_Sankey.html", "Figure1_Sankey.pdf")
  # 
}

### correlation heatmap

# RNA_exp_data <- as.matrix(RNA@assays$RNA@data)
# RNA_exp_data <- as.data.frame(t(RNA_exp_data))
# RNA_exp_data <- data.frame(Tissue = RNA$tissues,
#                            RNA_exp_data)
# RNA_exp_data <- aggregate.data.frame(RNA_exp_data[,2:ncol(RNA_exp_data)],
#                                      by = list(RNA_exp_data$Tissue),
#                                      FUN = mean)

### Tissue Marker

RAPDB_list <- list.files("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Figure1_tissue_marker_genes",
                         recursive = F, full.names = F)
RAPDB_markers <- c()
for (i in RAPDB_list) {
  temp <- read.csv(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Figure1_tissue_marker_genes/",
                          i), sep = "\t", header = T)
  tissue <- unlist(strsplit(i, " "))[1]
  temp <- data.frame(tissue = tissue,
                     temp)
  RAPDB_markers <- as.data.frame(rbind(RAPDB_markers, temp))
}
RAPDB_markers <- RAPDB_markers[,c(1,2,4,10,11)]

rownames(gene_id_map) <- gene_id_map$symbol
ATAC_tissue_markers <- getMarkerFeatures(
  ArchRProj = proj_pass_filter, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "TissuesSub",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
ATAC_tissue_markers_list <- getMarkers(ATAC_tissue_markers, cutOff = "FDR <= 0.05 & Log2FC >= 1")

ATAC_tissue_markers_list_DF <- c()
for (i in names(ATAC_tissue_markers_list)) {
  temp <- do.call(cbind, ATAC_tissue_markers_list[[i]][,-1])
  temp <- as.data.frame(temp)
  temp <- temp[,c("name","Log2FC","FDR")]
  temp[,2] <- as.numeric(temp[,2])
  temp[,3] <- as.numeric(temp[,3])
  colnames(temp)[1] <- "Symbol"
  temp$tissue <- i
  ATAC_tissue_markers_list_DF <- as.data.frame(rbind(temp,
                                                     ATAC_tissue_markers_list_DF))
}
ATAC_tissue_markers_list_DF$gene_id <- gene_id_map[ATAC_tissue_markers_list_DF$Symbol, "gene_id"]
ATAC_tissue_markers_list_DF <- ATAC_tissue_markers_list_DF[which(ATAC_tissue_markers_list_DF$FDR < 0.05),]

ATAC_tissue_markers_list_DF_RAPDB <- merge(ATAC_tissue_markers_list_DF,
                                           RAPDB_markers, by.x = "gene_id", by.y = "Locus.ID",
                                           all = F)

ATAC_tissue_markers_unique_top <- c()
ATAC_tissue_markers_unique <- c()
for (i in names(ATAC_tissue_markers_list)) {
  temp1 <- ATAC_tissue_markers_list_DF[which(ATAC_tissue_markers_list_DF$tissue == i),]
  temp2 <- ATAC_tissue_markers_list_DF[which(ATAC_tissue_markers_list_DF$tissue != i),]
  temp <- setdiff(temp1$Symbol, temp2$Symbol)
  #temp1 <- temp1[which(temp1$Symbol %in% temp),]
  temp1 <- temp1[order(temp1$Log2FC, decreasing = T),]
  ATAC_tissue_markers_unique_top <- as.data.frame(rbind(ATAC_tissue_markers_unique_top,
                                                        temp1[1:5,]))
  ATAC_tissue_markers_unique <- as.data.frame(rbind(ATAC_tissue_markers_unique,
                                                    temp1))
}

gene_id_map <- do.call(cbind, gene_annotation@listData[["genes"]]@elementMetadata@listData)
gene_id_map <- as.data.frame(gene_id_map)
rownames(gene_id_map) <- gene_id_map$symbol
gene_anno_info <- read.csv("./Ref/IRGSP-1.0_representative_annotation_2021-11-11.tsv",
                           sep = "\t", header = T)
gene_anno_info <- gene_anno_info[!duplicated(gene_anno_info$Locus_ID),]
rownames(gene_anno_info) <- gene_anno_info$Locus_ID

ATAC_tissue_markers_unique_top <- data.frame(ATAC_tissue_markers_unique_top,
                                             gene_anno_info[ATAC_tissue_markers_unique_top$Symbol,])
rownames(ATAC_tissue_markers_unique_top) <- ATAC_tissue_markers_unique_top$Symbol
ATAC_tissue_markers_unique_top$id <- gene_id_map[ATAC_tissue_markers_unique_top$Symbol,"gene_id"]

proj_pass_filter <- addImputeWeights(proj_pass_filter)
Tissue_ImputeMatrix <- getImputeWeights(proj_pass_filter)
Tissue_GeneScoreMatrix <- getMatrixFromProject(
  ArchRProj = proj_pass_filter,
  useMatrix = "GeneScoreMatrix",
  verbose = TRUE,
  binarize = FALSE
)
rownames(Tissue_GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]) <- Tissue_GeneScoreMatrix@elementMetadata@listData[["name"]]
# MAGIC_data <- imputeMatrix(mat = Tissue_GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]][c(ATAC_tissue_markers_unique_top$Symbol,
#                                                                                                     "Os10g0113549","Prol-26"),],
#                            imputeWeights = Tissue_ImputeMatrix)
sum(ATAC_tissue_markers_list_DF_RAPDB$Symbol %in% rownames(Tissue_GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]))
MAGIC_data <- imputeMatrix(mat = Tissue_GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]][unique(ATAC_tissue_markers_list_DF_RAPDB$Symbol),],
                           imputeWeights = Tissue_ImputeMatrix, threads = 1)

MAGIC_data <- as(MAGIC_data, Class = "CsparseMatrix")
MAGIC_data_TSNE <- as.data.frame(MAGIC_data)
MAGIC_data_TSNE <- as.data.frame(t(MAGIC_data_TSNE))

TSNE_withBatchEffect <- getEmbedding(proj_pass_filter, embedding = "TSNE_withBatchEffect")
colnames(TSNE_withBatchEffect) <- c("tSNE_1", "tSNE_2")

MAGIC_data_TSNE <- as.data.frame(cbind(MAGIC_data_TSNE,
                                       TSNE_withBatchEffect))
### Tissue Markers featurePlot
library(circlize)
library(momr)
col_fun <- colorRamp2(c(0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.3),
                      c("#191970", "#7B68EE", "#6495ED", "#F08080", "#FF4500", "#FF0000", "#8B0000"))
col <- c(seq(0,0.004,0.002), seq(0.01, 0.3, 0.005))
color_bar <- col_fun(col)
plotCol(color_bar)
make_plot <- function(j, legend = F) {
  MAGIC_data_TSNE <- MAGIC_data_TSNE[order(MAGIC_data_TSNE[,j], decreasing = F),]
  p <- ggplot() +
    geom_point(data = MAGIC_data_TSNE, aes(x = tSNE_1, y = tSNE_2, color = MAGIC_data_TSNE[,j]),
               size = 0.3) +
    scale_color_gradientn(colours = color_bar) +
    ggtitle(paste0(j," | ",ATAC_tissue_markers_unique_top[j,"tissue"],
                   "\n","ATAC, Log2FC = ",round(ATAC_tissue_markers_unique_top[j,"Log2FC"],2))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
         # axis.text = element_text(size = 12),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
         # axis.title = element_text(size = 14, face = "bold"),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.line = element_blank(),
          panel.border = element_blank()) +
    labs(color = j)
  if (legend == F) {
    p <- p +
      theme(legend.position = "none")
  }
  return(p)
}
make_plot(j = "MIL1", legend = T) # 
make_plot(j = "NAC106") # 
make_plot(j = "OsbHLH128") # 
make_plot(j = "OsPIN1D") # 
make_plot(j = "Os04g0630300") # 
make_plot(j = "Os10g0122200") # 
# 
# make_plot(j = "Prol.26")
plist <- lapply(c("Os04g0520800","Os02g0193100",
                  "Os01g0262000","Os04g0418600",
                  "Os04g0630300","Os10g0122200"), make_plot)
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)

pdf("./Figure_1_ATAC_tissue_top_marker_featurePlot_(legend).pdf")
make_plot(j = "Os04g0520800", legend = T) # 10dLeaf
dev.off()

pdf("./Figure_1_ATAC_tissue_top_marker_featurePlot.pdf", width = 9, height = 7)
p
dev.off()

saveRDS(proj_pass_filter, "2.proj_pass_filter.rds")


###### Add cell type info
library(Hmisc)
library(momr)
Leaf_ATAC_annotated <- readRDS("3.Leaf_ATAC_annotated.rds")
IM_ATAC_annotated <- readRDS("3.IM_ATAC_annotated.rds")
Shoot_ATAC_annotated <- readRDS("3.Shoot_ATAC_annotated.rds")
Root_ATAC_annotated <- readRDS("3.Root_ATAC_annotated.rds")

Leaf_ATAC_meta <- getCellColData(Leaf_ATAC_annotated)
Leaf_ATAC_meta <- as.data.frame(Leaf_ATAC_meta)
as.data.frame.array(table(Leaf_ATAC_meta$TissuesSub, Leaf_ATAC_meta$Celltype))
IM_ATAC_meta <- getCellColData(IM_ATAC_annotated)
IM_ATAC_meta <- as.data.frame(IM_ATAC_meta)
IM_ATAC_meta$Celltype <- capitalize(IM_ATAC_meta$Celltype)
as.data.frame.array(table(IM_ATAC_meta$TissuesSub, IM_ATAC_meta$Celltype))
Shoot_ATAC_meta <- getCellColData(Shoot_ATAC_annotated)
Shoot_ATAC_meta <- as.data.frame(Shoot_ATAC_meta)
as.data.frame.array(table(Shoot_ATAC_meta$TissuesSub, Shoot_ATAC_meta$Celltype))
Root_ATAC_meta <- getCellColData(Root_ATAC_annotated)
Root_ATAC_meta <- as.data.frame(Root_ATAC_meta)
as.data.frame.array(table(Root_ATAC_meta$TissuesSub, Root_ATAC_meta$Celltype))


All_cell_meta <- as.data.frame(rbind(data.frame(Cells = rownames(Leaf_ATAC_meta),
                                                Celltype = Leaf_ATAC_meta$Celltype),
                                     data.frame(Cells = rownames(IM_ATAC_meta),
                                                Celltype = IM_ATAC_meta$Celltype),
                                     data.frame(Cells = rownames(Shoot_ATAC_meta),
                                                Celltype = Shoot_ATAC_meta$Celltype),
                                     data.frame(Cells = rownames(Root_ATAC_meta),
                                                Celltype = Root_ATAC_meta$Celltype)))
rownames(All_cell_meta) <- All_cell_meta$Cells

proj_pass_filter$Celltype <- All_cell_meta[proj_pass_filter@cellColData@rownames, "Celltype"]
saveRDS(proj_pass_filter, "2.proj_pass_filter_annotated.rds")
length(unique(proj_pass_filter$Celltype))
TSNE_withBatchEffect <- getEmbedding(proj_pass_filter, embedding = "TSNE_withBatchEffect")
colnames(TSNE_withBatchEffect) <- c("tSNE_1", "tSNE_2")
ColData <- as.data.frame(getCellColData(proj_pass_filter))
TSNE_withBatchEffect <- as.data.frame(cbind(TSNE_withBatchEffect,
                                            ColData))
celltype_color <- as.character(ArchRPalettes$kelly)
celltype_color <- c(celltype_color, "#0088A8", "#FF7744")
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
plotCol(celltype_color)
names(celltype_color) <- sort(unique(TSNE_withBatchEffect$Celltype))
pdf("Figure1_ATAC_celltype_tSNE.pdf", width = 9.5, height = 6)
ggplot(data = TSNE_withBatchEffect, aes(x = tSNE_1, y = tSNE_2, color = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = celltype_color) + theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
dev.off()





