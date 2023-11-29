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
library(Rmagic)
library(circlize)
library(momr)

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

gene_id_map <- readRDS("gene_id_map.rds")
rownames(gene_id_map) <- gene_id_map$gene_id

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

############# 10dRoot epidermis root hair
Root_ATAC <- proj_pass_filter[which(proj_pass_filter$Celltype %in% c("Epidermis", "Root hair")),]
Root_ATAC <- Root_ATAC[which(Root_ATAC$TissuesSub %in% c("10dRoot")),]

unique(Root_ATAC$Celltype)
table(Root_ATAC$Celltype)
Root_ATAC <- addHarmony(
  ArchRProj = Root_ATAC,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = T
)
Root_ATAC <- addUMAP(
  ArchRProj = Root_ATAC, 
  reducedDims = "IterativeLSI", 
  name = "UMAP_withBatchEffect", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
Root_ATAC <- addUMAP(
  ArchRProj = Root_ATAC, 
  reducedDims = "Harmony", 
  name = "UMAP_Harmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
Root_ATAC <- addClusters(
  input = Root_ATAC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,
  force = TRUE
)
Root_ATAC <- addImputeWeights(Root_ATAC)

plotEmbedding(ArchRProj = Root_ATAC, colorBy = "cellColData",
              name = "Clusters", embedding = "UMAP_Harmony",
              size = 0.3, baseSize = 20)
State_cells <- readRDS("Epidermis_Roothair.rds")
Root_ATAC$Branch <- State_cells[Root_ATAC$cellNames, "Branch"]
markersGS <- getMarkerFeatures(
  ArchRProj = Root_ATAC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Branch",
  bias = c("log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 2")
markerGene <- c()
for (i in names(markerList)) {
  markerGene <- c(markerGene,
                  markerList[[i]][,"name"])
}


plotEmbedding(ArchRProj = Root_ATAC, colorBy = "cellColData",
              name = "Clusters", embedding = "UMAP_Harmony")

p <- plotBrowserTrack(
  ArchRProj = Root_ATAC, 
  groupBy = "Clusters", 
  geneSymbol = "Os10g0452700", 
  upstream = 20000,
  downstream = 20000
)
grid::grid.newpage()
grid::grid.draw(p$Os10g0452700)

p <- plotEmbedding(
  ArchRProj = Root_ATAC, 
  colorBy = "GeneScoreMatrix", 
  name = "Os10g0452700", 
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(Root_ATAC)
)
p

p <- plotEmbedding(
  ArchRProj = Root_ATAC, 
  colorBy = "GeneScoreMatrix", 
  name = "Os10g0452700", 
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(Root_ATAC)
)
p

gene_id_map[gene_id_map$gene_id == "Os01g0248900",]
p <- plotEmbedding(
  ArchRProj = Root_ATAC, 
  colorBy = "GeneScoreMatrix", 
  name = "EXPA8", 
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(Root_ATAC)
)
p

colorMat <- .getMatrixValues(
  ArchRProj = Root_ATAC, 
  name = "Os10g0452700", 
  matrixName = "GeneScoreMatrix", 
  log2Norm = FALSE
)
colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                         imputeWeights = getImputeWeights(Root_ATAC))

Root_ATAC_GS <- getMatrixFromProject(Root_ATAC, useMatrix = "GeneScoreMatrix")
Root_ATAC_GS@assays@data@listData[["GeneScoreMatrix"]]@Dimnames[[1]] <- Root_ATAC_GS@elementMetadata@listData[["name"]]
Root_ATAC_GS <- Root_ATAC_GS@assays@data@listData[["GeneScoreMatrix"]]
Root_ATAC_GS <- Root_ATAC_GS[-which(rowSums(Root_ATAC_GS) == 0),]
sum(rowSums(Root_ATAC_GS) == 0)
Root_ATAC_imputated <- magic(Root_ATAC_GS)

cell_meta <- getCellColData(Root_ATAC)
cells <- cell_meta@rownames
cell_meta <- as.data.frame(cell_meta@listData)
rownames(cell_meta) <- cells
seu <- CreateSeuratObject(Root_ATAC_imputated[["result"]])
seu <- AddMetaData(seu, cell_meta)
# seu$Os10g0452700 <- colorMat[1,colnames(seu)]

seu %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  RunTSNE() -> seu


library(monocle)
exp <- seu@assays$RNA@data
sample_ann <-  seu@meta.data
gene_ann <- data.frame(gene_short_name = rownames(exp),
                       row.names =  rownames(exp))
pd <- new("AnnotatedDataFrame", data = sample_ann)
fd <- new("AnnotatedDataFrame", data = gene_ann)
sc_cds <- newCellDataSet(
  exp, 
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = uninormal()
)
# sc_cds <- detectGenes(sc_cds, min_expr = 0.1) 
# sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- setOrderingFilter(sc_cds, markerGene)
sc_cds@assayData$exprs@Dim
sc_cds <- reduceDimension(sc_cds, max_components = 2,
                          method = "DDRTree", norm_method = "none")
sc_cds <- orderCells(sc_cds)
plot_cell_trajectory(sc_cds, color_by = "State")
sc_cds <- orderCells(sc_cds, root_state = "2")

lib_info_with_pseudo <- pData(sc_cds)
reduced_dim_coords <- reducedDimK(sc_cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))
dp_mst <- minSpanningTree(sc_cds)
edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from", 
                                                          target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name", 
                                                                                                                source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                       by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name", 
                                                                                                                                             target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                    by = "target")
mst_branch_nodes <- sc_cds@auxOrderingData[[sc_cds@dim_reduce_type]]$branch_points
branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes, 
                                                sample_name)) %>% mutate(branch_point_idx = seq_len(n()))

plot_cell_trajectory(sc_cds, color_by = "Celltype")
plot_cell_trajectory(sc_cds, color_by = "Pseudotime")

monocle_dims <- as.data.frame(t(monocle::reducedDimS(sc_cds)))
colnames(monocle_dims) <- c("Component 1", "Component 2")

pseudo_data <- data.frame(Cells = colnames(sc_cds),
                          monocle_dims[colnames(sc_cds), ],
                          Pseudotime = sc_cds$Pseudotime,
                          State = sc_cds$State,
                          row.names = colnames(sc_cds),
                          check.names = F)
colorMat <- .getMatrixValues(
  ArchRProj = Root_ATAC, 
  name = "Os10g0452700", 
  matrixName = "GeneScoreMatrix", 
  log2Norm = FALSE
)
colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                         imputeWeights = getImputeWeights(Root_ATAC))
pseudo_data$Os10g0452700 <- colorMat[rownames(pseudo_data)]
meta <- data.frame(Cells = colnames(sc_cds),
                   Celltype = sc_cds$Celltype,
                   row.names = colnames(sc_cds))
pseudo_data$Celltype <- meta[pseudo_data$Cells, "Celltype"]

pseudo_data_state_1 <- pseudo_data[which(pseudo_data$State == 1),]
pseudo_data_state_1 <- pseudo_data_state_1[order(pseudo_data_state_1$Os10g0452700, decreasing = T),]
pseudo_data_state_2 <- pseudo_data[which(pseudo_data$State == 2),]
pseudo_data_state_2 <- pseudo_data_state_2[order(pseudo_data_state_2$Os10g0452700, decreasing = F),]
pseudo_data_state_3 <- pseudo_data[which(pseudo_data$State == 3),]
pseudo_data_state_3 <- pseudo_data_state_3[order(pseudo_data_state_3$Os10g0452700, decreasing = T),]

pseudo_data_2 <- as.data.frame(rbind(pseudo_data_state_1,
                                     pseudo_data_state_2,
                                     pseudo_data_state_3))

pseudo_data <- pseudo_data[order(pseudo_data$Celltype, decreasing = T),]
pdf("Epi_roothair_Celltype_monocle.pdf", width = 5.2, height = 4.5)
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = 0.75, 
               linetype = "solid", na.rm = TRUE, data = edge_df) +
  geom_point(data = pseudo_data, aes(x = `Component 1`, y = `Component 2`, color = Celltype),
             size = 1.5) +
  scale_color_manual(values = c("#FF8888", "#C85C6C")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10, colour = "black"),
        axis.line = element_blank(),
        panel.border = element_blank()) +
  labs(color = "Celltype") +
  geom_point(aes_string(x = "prin_graph_dim_1", 
                        y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
             branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                     y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                          size = 4, color = "white", na.rm = TRUE, branch_point_df)
dev.off()

pdf("Epi_roothair_pseudotime_monocle.pdf", width = 5.2, height = 4.5)
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = 0.75, 
               linetype = "solid", na.rm = TRUE, data = edge_df) +
  geom_point(data = pseudo_data, aes(x = `Component 1`, y = `Component 2`, color = Pseudotime),
             size = 1.5) +
  scale_color_gradientn(colours = viridis::viridis(n = 20, option = "C")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10, colour = "black"),
        axis.line = element_blank(),
        panel.border = element_blank()) +
  labs(color = "Pseudotime") +
  geom_point(aes_string(x = "prin_graph_dim_1", 
                        y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
             branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                     y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                          size = 4, color = "white", na.rm = TRUE, branch_point_df)
dev.off()

pdf("Epi_roothair_state_monocle.pdf", width = 5.2, height = 4.5)
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = 0.75, 
               linetype = "solid", na.rm = TRUE, data = edge_df) +
  geom_point(data = pseudo_data, aes(x = `Component 1`, y = `Component 2`, color = State),
             size = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10, colour = "black"),
        axis.line = element_blank(),
        panel.border = element_blank()) +
  labs(color = "State") +
  geom_point(aes_string(x = "prin_graph_dim_1", 
                        y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
             branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                     y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                          size = 4, color = "white", na.rm = TRUE, branch_point_df)
dev.off()


pdf("Epi_roothair_Os10g0452700_monocle.pdf", width = 5.8, height = 4.5)
ggplot() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = 0.75, 
               linetype = "solid", na.rm = TRUE, data = edge_df) +
  geom_point(data = pseudo_data_2, aes(x = `Component 1`, y = `Component 2`, color = Os10g0452700),
             size = 1.5) +
  scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10, colour = "black"),
        axis.line = element_blank(),
        panel.border = element_blank()) +
  labs(color = "Os10g0452700") +
  geom_point(aes_string(x = "prin_graph_dim_1", 
                        y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
             branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                     y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                          size = 4, color = "white", na.rm = TRUE, branch_point_df)
dev.off()


Root_ATAC$State <- pseudo_data[Root_ATAC$cellNames, "State"]
Root_ATAC$State <- as.character(Root_ATAC$State)
markersGS_State <- getMarkerFeatures(
  ArchRProj = Root_ATAC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "State",
  bias = c("log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList_State <- getMarkers(markersGS_State, cutOff = "FDR <= 0.05 & Log2FC >= 1")

# GO
{
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
  genes_list <- c()
  rownames(gene_id_map) <- gene_id_map$symbol
  for (i in names(markerList_State)) {
    temp <- markerList_State[[i]]
    temp <- as.data.frame(temp@listData)
    temp <- unique(temp$name)
    temp <- gene_id_map[temp, "gene_id"]
    genes_list <- c(genes_list, list(temp))
  }
  names(genes_list) <- names(markerList_State)
  GO_result_State <- c()
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
    GO_result$State <- i
    GO_result_State <- as.data.frame(rbind(GO_result_State, GO_result))
  }
  
  GO_result_State <- GO_result_State[which(GO_result_State$Ontology == "BP"),]
  GO_result_State$Branch <- ""
  GO_result_State$Branch[which(GO_result_State$State == "2")] <- "Pre-branch"
  GO_result_State$Branch[which(GO_result_State$State == "1")] <- "Branch 2"
  GO_result_State$Branch[which(GO_result_State$State == "3")] <- "Branch 1"
  GO_result_State_2 <- GO_result_State
  GO_result_State <- GO_result_State[which(GO_result_State$P_value < 0.05),]
  
  
  openxlsx::write.xlsx(GO_result_State, "Epi_root_hair_branch_GO.xlsx")
  
  GO_result_State_temp <- GO_result_State_2[which(GO_result_State_2$Description %in% c("mitotic cell cycle phase transition",
                                                                                       "meiotic nuclear division",
                                                                                       "cellular developmental process",
                                                                                       "cell wall organization",
                                                                                       "cell wall biogenesis",
                                                                                       "cell wall macromolecule metabolic process",
                                                                                       "response to stress",
                                                                                       "organic acid transport",
                                                                                       "ethylene-activated signaling pathway")),]
  GO_result_State_temp$Description <- factor(GO_result_State_temp$Description,
                                             levels = c("mitotic cell cycle phase transition",
                                                        "meiotic nuclear division",
                                                        "cellular developmental process",
                                                        "cell wall organization",
                                                        "cell wall biogenesis",
                                                        "cell wall macromolecule metabolic process",
                                                        "response to stress",
                                                        "organic acid transport",
                                                        "ethylene-activated signaling pathway"))
  GO_result_State_temp$P_value_log10 <- -log10(GO_result_State_temp$P_value)
  P <- matrix(NA, nrow = 9, ncol = 3)
  rownames(P) <- c("mitotic cell cycle phase transition",
                   "meiotic nuclear division",
                   "cellular developmental process",
                   "cell wall organization",
                   "cell wall biogenesis",
                   "cell wall macromolecule metabolic process",
                   "response to stress",
                   "organic acid transport",
                   "ethylene-activated signaling pathway")
  colnames(P) <- unique(GO_result_State_temp$Branch)
  for (i in 1:nrow(GO_result_State_temp)) {
    P[GO_result_State_temp[i,3], GO_result_State_temp[i, 9]] <- GO_result_State_temp[i, 7]
  }
  # P[P > 1.5] <- 1.5
  pdf("Epi_roothair_GO_heatmap.pdf", height = 4, width = 4.1)
  pheatmap(P, cluster_rows = F, cluster_cols = F)
  dev.off()
  
  ggplot(data = GO_result_State_temp, aes(x = P_value_log10, y = Description, fill = Branch)) +
    geom_bar(stat = "identity") +
    facet_wrap(.~Branch) +
    geom_vline(xintercept = -log10(0.05))
  
  
}

# marker gene
{
  # OsGSTU5
  gene_id_map[which(gene_id_map$gene_id == "Os09g0367700"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "OsGSTU5", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$OsGSTU5 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$OsGSTU5[which(pseudo_data$OsGSTU5 > 1.45)] <- 1.45
  {
    pdf("Epi_roothair_OsGSTU5_roothair.pdf", width = 5, height = 3)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$OsGSTU5),], aes(x = `Component 1`, y = `Component 2`, color = OsGSTU5),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsGSTU5") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  # OsCSLD1
  gene_id_map[which(gene_id_map$symbol == "CSLD1"),]
  gene_id_map[which(gene_id_map$gene_id == "Os10g0578200"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "CSLD1", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$CSLD1 <- colorMat[rownames(pseudo_data)]
  {
    pdf("Epi_roothair_CSLD1_epi.pdf", width = 5, height = 3)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$CSLD1),], aes(x = `Component 1`, y = `Component 2`, color = CSLD1),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsCSLD1") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  # LPL3
  gene_id_map[which(gene_id_map$symbol == "LPL3"),]
  gene_id_map[which(gene_id_map$gene_id == "Os08g0544500"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "LPL3", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$LPL3 <- colorMat[rownames(pseudo_data)]
  pseudo_data$LPL3[which(pseudo_data$LPL3 > 0.55)] <- 0.55
  {
    pdf("Epi_roothair_LPL3_epi.pdf", width = 5, height = 3)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$LPL3),], aes(x = `Component 1`, y = `Component 2`, color = LPL3),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsLPL3") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # LPL2
  gene_id_map[which(gene_id_map$symbol == "LPL2"),]
  gene_id_map[which(gene_id_map$gene_id == "Os03g0143800"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "LPL2", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$LPL2 <- colorMat[rownames(pseudo_data)]
  {
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$LPL2),], aes(x = `Component 1`, y = `Component 2`, color = LPL2),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "LPL2") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
  }
  
  
  
  # SRL1
  gene_id_map[which(gene_id_map$symbol == "SRL1"),]
  gene_id_map[which(gene_id_map$gene_id == "Os07g0102300"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "SRL1", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$SRL1 <- colorMat[rownames(pseudo_data)]
  {
    pdf("Epi_roothair_SRL1_epi.pdf", width = 5, height = 3)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$SRL1),], aes(x = `Component 1`, y = `Component 2`, color = SRL1),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "SRL1") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # OsbHLH126,上皮; OsGH9C1,上皮；OsEXPB2,根毛；OsFH1,上皮；OsRac3,上皮；OsbHLH127,上皮；OsbHLH128,根毛
  gene_id_map[which(gene_id_map$symbol == "OsEXPB2"),]
  gene_id_map[which(gene_id_map$gene_id == "Os10g0555700"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "OsEXPB2", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$OsEXPB2 <- colorMat[rownames(pseudo_data)]
  pseudo_data$OsEXPB2[which(pseudo_data$OsEXPB2 > 0.75)] <- 0.75
  {
    pdf("Epi_roothair_OsEXPB2_roothair.pdf", width = 5, height = 3)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$OsEXPB2),], aes(x = `Component 1`, y = `Component 2`, color = OsEXPB2),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsEXPB2") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  gene_id_map[which(gene_id_map$symbol == "OsbHLH128"),]
  gene_id_map[which(gene_id_map$gene_id == "Os07g0588400"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "OsbHLH128", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$OsbHLH128 <- colorMat[rownames(pseudo_data)]
  pseudo_data$OsbHLH128[which(pseudo_data$OsbHLH128 > 0.9)] <- 0.9
  {
    pdf("Epi_roothair_OsbHLH128_roothair.pdf", width = 5, height = 3)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$OsbHLH128),], aes(x = `Component 1`, y = `Component 2`, color = OsbHLH128),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsbHLH128") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  # epidermis, OsCesA4
  gene_id_map[which(gene_id_map$symbol == "OsCesA4"),]
  gene_id_map[which(gene_id_map$gene_id == "Os01g0750300"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "OsCesA4", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$OsCesA4 <- colorMat[rownames(pseudo_data)]
  pseudo_data$OsCesA4[which(pseudo_data$OsCesA4 > 0.4)] <- 0.4
  {
    pdf("Epi_roothair_sup_OsCesA4_epi.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$OsCesA4),], aes(x = `Component 1`, y = `Component 2`, color = OsCesA4),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsCesA4") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  # epidermis, OsUXS
  gene_id_map[which(gene_id_map$symbol == "UXS"),]
  gene_id_map[which(gene_id_map$gene_id == "Os03g0278000"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "UXS", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$UXS <- colorMat[rownames(pseudo_data)]
  pseudo_data$UXS[which(pseudo_data$UXS > 0.6)] <- 0.6
  {
    pdf("Epi_roothair_sup_OsUXS_epi.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$UXS),], aes(x = `Component 1`, y = `Component 2`, color = UXS),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "UXS") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # epidermis, OsCESA9
  gene_id_map[which(gene_id_map$symbol == "BC6"),]
  gene_id_map[which(gene_id_map$gene_id == "Os09g0422500"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "BC6", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$BC6 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$BC6[which(pseudo_data$BC6 > 0.6)] <- 0.6
  {
    pdf("Epi_roothair_sup_OsCESA9_epi.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$BC6),], aes(x = `Component 1`, y = `Component 2`, color = BC6),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsCESA9") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # epidermis, OsPRX39
  gene_id_map[which(gene_id_map$symbol == "PRX39"),]
  gene_id_map[which(gene_id_map$gene_id == "Os03g0234900"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "PRX39", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$PRX39 <- colorMat[rownames(pseudo_data)]
  pseudo_data$PRX39[which(pseudo_data$PRX39 > 0.75)] <- 0.75
  {
    pdf("Epi_roothair_sup_OsPRX39_epi.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$PRX39),], aes(x = `Component 1`, y = `Component 2`, color = PRX39),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsPRX39") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  # epidermis, TagI.4
  gene_id_map[which(gene_id_map$symbol == "TagI.4"),]
  gene_id_map[which(gene_id_map$gene_id == "Os08g0489300"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "TagI.4", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$TagI.4 <- colorMat[rownames(pseudo_data)]
  pseudo_data$TagI.4[which(pseudo_data$TagI.4 > 0.75)] <- 0.75
  {
    pdf("Epi_roothair_sup_TagI.4_epi.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$TagI.4),], aes(x = `Component 1`, y = `Component 2`, color = TagI.4),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "TagI.4") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # epidermis, prx71
  gene_id_map[which(gene_id_map$symbol == "prx71"),]
  gene_id_map[which(gene_id_map$gene_id == "Os05g0135500"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "prx71", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$prx71 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$prx71[which(pseudo_data$prx71 > 0.75)] <- 0.75
  {
    pdf("Epi_roothair_sup_prx71_epi.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$prx71),], aes(x = `Component 1`, y = `Component 2`, color = prx71),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "prx71") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # roothair, OsALDH2C1
  gene_id_map[which(gene_id_map$symbol == "OsALDH2C1"),]
  gene_id_map[which(gene_id_map$gene_id == "Os01g0591300"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "OsALDH2C1", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$OsALDH2C1 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$OsALDH2C1[which(pseudo_data$OsALDH2C1 > 0.75)] <- 0.75
  {
    pdf("Epi_roothair_sup_OsALDH2C1_roothair.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$OsALDH2C1),], aes(x = `Component 1`, y = `Component 2`, color = OsALDH2C1),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OsALDH2C1") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # roothair, OSINV2
  gene_id_map[which(gene_id_map$symbol == "OSINV2"),]
  gene_id_map[which(gene_id_map$gene_id == "Os04g0535600"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "OSINV2", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$OSINV2 <- colorMat[rownames(pseudo_data)]
  pseudo_data$OSINV2[which(pseudo_data$OSINV2 > 3)] <- 3
  {
    pdf("Epi_roothair_sup_OSINV2_roothair.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$OSINV2),], aes(x = `Component 1`, y = `Component 2`, color = OSINV2),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "OSINV2") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  
  # roothair, Os11g0149900
  gene_id_map[which(gene_id_map$symbol == "Os11g0149900"),]
  gene_id_map[which(gene_id_map$gene_id == "Os11g0149900"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "Os11g0149900", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$Os11g0149900 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$Os11g0149900[which(pseudo_data$Os11g0149900 > 3)] <- 3
  {
    pdf("Epi_roothair_sup_Os11g0149900_roothair.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$Os11g0149900),], aes(x = `Component 1`, y = `Component 2`, color = Os11g0149900),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "Os11g0149900") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # roothair, Os07g0659300
  gene_id_map[which(gene_id_map$symbol == "Os07g0659300"),]
  gene_id_map[which(gene_id_map$gene_id == "Os07g0659300"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "Os07g0659300", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$Os07g0659300 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$Os07g0659300[which(pseudo_data$Os07g0659300 > 3)] <- 3
  {
    pdf("Epi_roothair_sup_Os07g0659300_roothair.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$Os07g0659300),], aes(x = `Component 1`, y = `Component 2`, color = Os07g0659300),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "Os07g0659300") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  
  # roothair, HCT3
  gene_id_map[which(gene_id_map$symbol == "HCT3"),]
  gene_id_map[which(gene_id_map$gene_id == "Os09g0422000"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "HCT3", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$HCT3 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$HCT3[which(pseudo_data$HCT3 > 3)] <- 3
  {
    pdf("Epi_roothair_sup_HCT3_roothair.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$HCT3),], aes(x = `Component 1`, y = `Component 2`, color = HCT3),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "HCT3") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
  # roothair, Os11g0649400
  gene_id_map[which(gene_id_map$symbol == "Os11g0649400"),]
  gene_id_map[which(gene_id_map$gene_id == "Os11g0649400"),]
  colorMat <- .getMatrixValues(
    ArchRProj = Root_ATAC, 
    name = "Os11g0649400", 
    matrixName = "GeneScoreMatrix", 
    log2Norm = FALSE
  )
  colorMat <- imputeMatrix(mat = as.matrix(colorMat),
                           imputeWeights = getImputeWeights(Root_ATAC))
  pseudo_data$Os11g0649400 <- colorMat[rownames(pseudo_data)]
  # pseudo_data$Os11g0649400[which(pseudo_data$Os11g0649400 > 3)] <- 3
  {
    pdf("Epi_roothair_sup_Os11g0649400_roothair.pdf", width = 5, height = 4)
    ggplot() +
      geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                              y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                              yend = "target_prin_graph_dim_2"), size = 0.75, 
                   linetype = "solid", na.rm = TRUE, data = edge_df) +
      geom_point(data = pseudo_data[order(pseudo_data$Os11g0649400),], aes(x = `Component 1`, y = `Component 2`, color = Os11g0649400),
                 size = 1.5) +
      scale_color_gradientn(colours = c("#191970", "#7B68EE", "#6495ED", "#FFE4E1", "#FF4500", "#FF0000", "#8B0000")) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 10, colour = "black"),
            axis.line = element_blank(),
            panel.border = element_blank()) +
      labs(color = "Os11g0649400") +
      geom_point(aes_string(x = "prin_graph_dim_1", 
                            y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, 
                 branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                         y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                              size = 4, color = "white", na.rm = TRUE, branch_point_df)
    dev.off()
  }
  
  
}

#### peaks
Root_ATAC$Pseudotime <- pseudo_data[Root_ATAC$cellNames, "Pseudotime"]
# Root_ATAC_sub <- Root_ATAC[which(Root_ATAC$State %in% c("1", "2")),]
Root_ATAC$Branch <- ""
Root_ATAC$Branch[which(Root_ATAC$State == "2")] <- "Pre-branch"
Root_ATAC$Branch[which(Root_ATAC$State == "1")] <- "Branch 2"
Root_ATAC$Branch[which(Root_ATAC$State == "3")] <- "Branch 1"
# Root_ATAC_sub$Branch <- ifelse(Root_ATAC_sub$State == "2", "Pre-branch", "Branch 2")
markersGS_Branch <- getMarkerFeatures(
  ArchRProj = Root_ATAC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Branch",
  bias = c("log10(nFrags)", "TSSEnrichment"),
  testMethod = "wilcoxon"
)

markerGSList_Branch <- getMarkers(markersGS_Branch, cutOff = "FDR <= 0.05 & Log2FC >= 1")

markersPK_Branch <- getMarkerFeatures(
  ArchRProj = Root_ATAC, 
  useMatrix = "PeakMatrix", 
  groupBy = "Branch",
  bias = c("log10(nFrags)", "TSSEnrichment"),
  testMethod = "wilcoxon"
)

markerPKList_Branch <- getMarkers(markersPK_Branch, cutOff = "FDR <= 0.05 & Log2FC >= 1")

markerPKList_Branch_peak_id <- c()
for (i in names(markerPKList_Branch)) {
  temp <- markerPKList_Branch[[i]]
  temp <- as.data.frame(temp)
  temp$ID <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
  markerPKList_Branch_peak_id <- c(markerPKList_Branch_peak_id,
                                   list(temp))
}
names(markerPKList_Branch_peak_id) <- names(markerPKList_Branch)
marker_peak_id <- c(markerPKList_Branch_peak_id$`Branch 2`$ID,
                    markerPKList_Branch_peak_id$`Pre-branch`$ID,
                    markerPKList_Branch_peak_id$`Branch 1`$ID)
sum(duplicated(marker_peak_id))

duplicated_peaks <- marker_peak_id[duplicated(marker_peak_id)]
for (i in names(markerPKList_Branch_peak_id)) {
  temp <- markerPKList_Branch_peak_id[[i]]
  temp <- temp[-which(temp$ID %in% duplicated_peaks),]
  markerPKList_Branch_peak_id[[i]] <- temp
}
marker_peak_id <- c(markerPKList_Branch_peak_id$`Branch 2`$ID,
                    markerPKList_Branch_peak_id$`Pre-branch`$ID,
                    markerPKList_Branch_peak_id$`Branch 1`$ID)
sum(duplicated(marker_peak_id))

# peak matrix
Root_ATAC_peak_mat <- getMatrixFromProject(Root_ATAC,
                                           useMatrix = "PeakMatrix")
temp <- as.data.frame(Root_ATAC_peak_mat@rowRanges)
temp$ID <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
Root_ATAC_peak_mat <- Root_ATAC_peak_mat@assays@data@listData[["PeakMatrix"]]
rownames(Root_ATAC_peak_mat) <- temp$ID

branch_2 <- Root_ATAC[which(Root_ATAC$Branch %in% c("Branch 2")),]
min(branch_2$Pseudotime)
max(branch_2$Pseudotime)
branch_1 <- Root_ATAC[which(Root_ATAC$Branch %in% c("Branch 1")),]
min(branch_1$Pseudotime)
max(branch_1$Pseudotime)
pre_branch <- Root_ATAC[which(Root_ATAC$Branch %in% c("Pre-branch")),]

library(infotheo)
# pre branch
{
  Prebranch_mat <- Root_ATAC_peak_mat[c(markerPKList_Branch_peak_id$`Pre-branch`$ID,
                                        markerPKList_Branch_peak_id$`Branch 1`$ID,
                                        markerPKList_Branch_peak_id$`Branch 2`$ID), rownames(pre_branch)]
  Prebranch_mat <- t(Prebranch_mat)
  Prebranch_mat <- as.matrix(Prebranch_mat)
  Prebranch_mat <- as.data.frame(Prebranch_mat)
  Pseudotime <- data.frame(Cells = rownames(Prebranch_mat),
                           Pseudotime = pre_branch@cellColData[rownames(Prebranch_mat), c("Pseudotime")])
  nbins <- 60
  equal_width <- discretize(Pseudotime$Pseudotime, "equalwidth", nbins)
  Pseudotime$bin <- equal_width$X
  Pseudotime <- Pseudotime[order(Pseudotime$Pseudotime, Pseudotime$bin, decreasing = F),]
  Prebranch_mat_2 <- c()
  Prebranch_mat <- Prebranch_mat[Pseudotime$Cells, ]
  for (i in unique(Pseudotime$bin)) {
    print(i)
    temp <- Prebranch_mat[which(Pseudotime$bin == i), , drop = F]
    temp <- colMeans(temp)
    temp <- as.data.frame(temp)
    colnames(temp) <- i
    temp <- as.data.frame(t(temp))
    Prebranch_mat_2 <- as.data.frame(rbind(Prebranch_mat_2, temp))
  }
  Prebranch_mat_2 <- Prebranch_mat_2 + rnorm(n = length(Prebranch_mat_2), mean = 0.001, sd = 0.0001)
  Prebranch_mat_2 <- t(Prebranch_mat_2)
  Prebranch_mat_2 <- Prebranch_mat_2[, as.character(sort(as.numeric(colnames(Prebranch_mat_2))))]
  Prebranch_mat_2 <- sweep(Prebranch_mat_2 - rowMeans(Prebranch_mat_2), 1, matrixStats::rowSds(Prebranch_mat_2), `/`)
  
  Prebranch_mat_2_top1 <- Prebranch_mat_2[markerPKList_Branch_peak_id$`Pre-branch`$ID,]
  Prebranch_mat_2_top2 <- Prebranch_mat_2[markerPKList_Branch_peak_id$`Branch 1`$ID,]
  Prebranch_mat_2_top3 <- Prebranch_mat_2[markerPKList_Branch_peak_id$`Branch 2`$ID,]
  
  max(Prebranch_mat_2)
  min(Prebranch_mat_2)
  max(Prebranch_mat_2_top1)
  min(Prebranch_mat_2_top1)
  Prebranch_mat_2[Prebranch_mat_2 > 3] <- 3
  Prebranch_mat_2[Prebranch_mat_2 < -3] <- -3
  
  Prebranch_mat_2_top1 <- Prebranch_mat_2[markerPKList_Branch_peak_id$`Pre-branch`$ID,]
  Prebranch_mat_2_top2 <- Prebranch_mat_2[markerPKList_Branch_peak_id$`Branch 1`$ID,]
  Prebranch_mat_2_top3 <- Prebranch_mat_2[markerPKList_Branch_peak_id$`Branch 2`$ID,]
  
  idx <- order(apply(Prebranch_mat_2_top1, 1, which.max))
  Prebranch_mat_2_top1_ordered <- Prebranch_mat_2_top1[idx, ]
  Prebranch_mat_2_top1_ordered_ID <- rownames(Prebranch_mat_2_top1_ordered)
  
  Prebranch_mat_2_top1_right <- Prebranch_mat_2_top1
  colnames(Prebranch_mat_2_top1_right) <- paste0("prebranch_right_", colnames(Prebranch_mat_2_top1_right))
  Prebranch_mat_2_top1_left <- Prebranch_mat_2_top1
  Prebranch_mat_2_top1_left <- Prebranch_mat_2_top1_left[,ncol(Prebranch_mat_2_top1_left):1]
  colnames(Prebranch_mat_2_top1_left) <- paste0("prebranch_left_", colnames(Prebranch_mat_2_top1_left))
  
  Prebranch_mat_2_top2_right <- Prebranch_mat_2_top2
  colnames(Prebranch_mat_2_top2_right) <- paste0("prebranch_right_", colnames(Prebranch_mat_2_top2_right))
  Prebranch_mat_2_top2_left <- Prebranch_mat_2_top2
  Prebranch_mat_2_top2_left <- Prebranch_mat_2_top2_left[,ncol(Prebranch_mat_2_top2_left):1]
  colnames(Prebranch_mat_2_top2_left) <- paste0("prebranch_left_", colnames(Prebranch_mat_2_top2_left))
  
  Prebranch_mat_2_top3_right <- Prebranch_mat_2_top3
  colnames(Prebranch_mat_2_top3_right) <- paste0("prebranch_right_", colnames(Prebranch_mat_2_top3_right))
  Prebranch_mat_2_top3_left <- Prebranch_mat_2_top3
  Prebranch_mat_2_top3_left <- Prebranch_mat_2_top3_left[,ncol(Prebranch_mat_2_top3_left):1]
  colnames(Prebranch_mat_2_top3_left) <- paste0("prebranch_left_", colnames(Prebranch_mat_2_top3_left))
  
}


# branch 1
{
  Branch1_mat <- Root_ATAC_peak_mat[c(markerPKList_Branch_peak_id$`Pre-branch`$ID,
                                        markerPKList_Branch_peak_id$`Branch 1`$ID,
                                        markerPKList_Branch_peak_id$`Branch 2`$ID), rownames(branch_1)]
  Branch1_mat <- t(Branch1_mat)
  Branch1_mat <- as.matrix(Branch1_mat)
  Branch1_mat <- as.data.frame(Branch1_mat)
  Pseudotime <- data.frame(Cells = rownames(Branch1_mat),
                           Pseudotime = branch_1@cellColData[rownames(Branch1_mat), c("Pseudotime")])
  nbins <- 40
  equal_width <- discretize(Pseudotime$Pseudotime, "equalwidth", nbins)
  Pseudotime$bin <- equal_width$X
  Pseudotime <- Pseudotime[order(Pseudotime$Pseudotime, Pseudotime$bin, decreasing = F),]
  Branch1_mat_2 <- c()
  Branch1_mat <- Branch1_mat[Pseudotime$Cells, ]
  for (i in unique(Pseudotime$bin)) {
    print(i)
    temp <- Branch1_mat[which(Pseudotime$bin == i), , drop = F]
    temp <- colMeans(temp)
    temp <- as.data.frame(temp)
    colnames(temp) <- i
    temp <- as.data.frame(t(temp))
    Branch1_mat_2 <- as.data.frame(rbind(Branch1_mat_2, temp))
  }
  Branch1_mat_2 <- Branch1_mat_2 + rnorm(n = length(Branch1_mat_2), mean = 0.001, sd = 0.0001)
  Branch1_mat_2 <- t(Branch1_mat_2)
  Branch1_mat_2 <- Branch1_mat_2[, as.character(sort(as.numeric(colnames(Branch1_mat_2))))]
  Branch1_mat_2 <- sweep(Branch1_mat_2 - rowMeans(Branch1_mat_2), 1, matrixStats::rowSds(Branch1_mat_2), `/`)
  
  Branch1_mat_2_top1 <- Branch1_mat_2[markerPKList_Branch_peak_id$`Pre-branch`$ID,]
  Branch1_mat_2_top2 <- Branch1_mat_2[markerPKList_Branch_peak_id$`Branch 1`$ID,]
  Branch1_mat_2_top3 <- Branch1_mat_2[markerPKList_Branch_peak_id$`Branch 2`$ID,]
  
  max(Branch1_mat_2)
  min(Branch1_mat_2)
  max(Branch1_mat_2_top2)
  min(Branch1_mat_2_top2)
  Branch1_mat_2[Branch1_mat_2 > 3] <- 3
  Branch1_mat_2[Branch1_mat_2 < -3] <- -3
  
  Branch1_mat_2_top1 <- Branch1_mat_2[markerPKList_Branch_peak_id$`Pre-branch`$ID,]
  Branch1_mat_2_top2 <- Branch1_mat_2[markerPKList_Branch_peak_id$`Branch 1`$ID,]
  Branch1_mat_2_top3 <- Branch1_mat_2[markerPKList_Branch_peak_id$`Branch 2`$ID,]
  
  idx <- order(apply(Branch1_mat_2_top2, 1, which.max))
  Branch1_mat_2_top2_ordered <- Branch1_mat_2_top2[idx, ]
  Branch1_mat_2_top2_ordered_ID <- rownames(Branch1_mat_2_top2_ordered)
  
  Branch1_mat_2_top1_right <- Branch1_mat_2_top1
  colnames(Branch1_mat_2_top1_right) <- paste0("branch1_right_", colnames(Branch1_mat_2_top1_right))
  
  Branch1_mat_2_top2_right <- Branch1_mat_2_top2
  colnames(Branch1_mat_2_top2_right) <- paste0("branch1_right_", colnames(Branch1_mat_2_top2_right))
  
  Branch1_mat_2_top3_right <- Branch1_mat_2_top3
  colnames(Branch1_mat_2_top3_right) <- paste0("branch1_right_", colnames(Branch1_mat_2_top3_right))
  
}


# branch 2
{
  Branch2_mat <- Root_ATAC_peak_mat[c(markerPKList_Branch_peak_id$`Pre-branch`$ID,
                                      markerPKList_Branch_peak_id$`Branch 1`$ID,
                                      markerPKList_Branch_peak_id$`Branch 2`$ID), rownames(branch_2)]
  Branch2_mat <- t(Branch2_mat)
  Branch2_mat <- as.matrix(Branch2_mat)
  Branch2_mat <- as.data.frame(Branch2_mat)
  Pseudotime <- data.frame(Cells = rownames(Branch2_mat),
                           Pseudotime = branch_2@cellColData[rownames(Branch2_mat), c("Pseudotime")])
  nbins <- 40
  equal_width <- discretize(Pseudotime$Pseudotime, "equalwidth", nbins)
  Pseudotime$bin <- equal_width$X
  Pseudotime <- Pseudotime[order(Pseudotime$Pseudotime, Pseudotime$bin, decreasing = F),]
  Branch2_mat_2 <- c()
  Branch2_mat <- Branch2_mat[Pseudotime$Cells, ]
  for (i in unique(Pseudotime$bin)) {
    print(i)
    temp <- Branch2_mat[which(Pseudotime$bin == i), , drop = F]
    temp <- colMeans(temp)
    temp <- as.data.frame(temp)
    colnames(temp) <- i
    temp <- as.data.frame(t(temp))
    Branch2_mat_2 <- as.data.frame(rbind(Branch2_mat_2, temp))
  }
  Branch2_mat_2 <- Branch2_mat_2 + rnorm(n = length(Branch2_mat_2), mean = 0.001, sd = 0.0001)
  Branch2_mat_2 <- t(Branch2_mat_2)
  Branch2_mat_2 <- Branch2_mat_2[, as.character(sort(as.numeric(colnames(Branch2_mat_2))))]
  Branch2_mat_2 <- sweep(Branch2_mat_2 - rowMeans(Branch2_mat_2), 1, matrixStats::rowSds(Branch2_mat_2), `/`)
  
  Branch2_mat_2_top1 <- Branch2_mat_2[markerPKList_Branch_peak_id$`Pre-branch`$ID,]
  Branch2_mat_2_top2 <- Branch2_mat_2[markerPKList_Branch_peak_id$`Branch 1`$ID,]
  Branch2_mat_2_top3 <- Branch2_mat_2[markerPKList_Branch_peak_id$`Branch 2`$ID,]
  
  max(Branch2_mat_2)
  min(Branch2_mat_2)
  max(Branch2_mat_2_top3)
  min(Branch2_mat_2_top3)
  Branch2_mat_2[Branch2_mat_2 > 3] <- 3
  Branch2_mat_2[Branch2_mat_2 < -3] <- -3
  Branch2_mat_2 <- Branch2_mat_2[,ncol(Branch2_mat_2):1]
  
  Branch2_mat_2_top1 <- Branch2_mat_2[markerPKList_Branch_peak_id$`Pre-branch`$ID,]
  Branch2_mat_2_top2 <- Branch2_mat_2[markerPKList_Branch_peak_id$`Branch 1`$ID,]
  Branch2_mat_2_top3 <- Branch2_mat_2[markerPKList_Branch_peak_id$`Branch 2`$ID,]
  
  idx <- order(apply(Branch2_mat_2_top3, 1, which.max))
  Branch2_mat_2_top3_ordered <- Branch2_mat_2_top3[idx, ]
  Branch2_mat_2_top3_ordered_ID <- rownames(Branch2_mat_2_top3_ordered)
  
  Branch2_mat_2_top1_left <- Branch2_mat_2_top1
  colnames(Branch2_mat_2_top1_left) <- paste0("branch2_left_", colnames(Branch2_mat_2_top1_left))
  
  Branch2_mat_2_top2_left <- Branch2_mat_2_top2
  colnames(Branch2_mat_2_top2_left) <- paste0("branch2_left_", colnames(Branch2_mat_2_top2_left))
  
  Branch2_mat_2_top3_left <- Branch2_mat_2_top3
  colnames(Branch2_mat_2_top3_left) <- paste0("branch2_left_", colnames(Branch2_mat_2_top3_left))
  
}

peak_matrix_ordered <- as.data.frame(rbind(Prebranch_mat_2_top1_left[Prebranch_mat_2_top1_ordered_ID,],
                                           Prebranch_mat_2_top2_left[Branch1_mat_2_top2_ordered_ID,],
                                           Prebranch_mat_2_top3_left[rev(Branch2_mat_2_top3_ordered_ID),]))
Heatmap(as.matrix(peak_matrix_ordered), show_row_names = F,
        cluster_rows = F, use_raster = T, raster_resize_mat = F,
        show_column_names = F, cluster_columns = F,
        raster_device = "png", raster_quality = 5, raster_by_magick = F,
        col = paletteContinuous(set = "solarExtra", n = 100),
        name = paste0(nrow(peak_matrix_ordered), " Peaks"))
peak_matrix_ordered <- as.data.frame(cbind(peak_matrix_ordered,
                                           rbind(Prebranch_mat_2_top1_right[Prebranch_mat_2_top1_ordered_ID,],
                                                 Prebranch_mat_2_top2_right[Branch1_mat_2_top2_ordered_ID,],
                                                 Prebranch_mat_2_top3_right[rev(Branch2_mat_2_top3_ordered_ID),])))

peak_matrix_ordered <- as.data.frame(cbind(peak_matrix_ordered,
                                           rbind(Branch1_mat_2_top1_right[Prebranch_mat_2_top1_ordered_ID,],
                                                 Branch1_mat_2_top2_right[Branch1_mat_2_top2_ordered_ID,],
                                                 Branch1_mat_2_top3_right[rev(Branch2_mat_2_top3_ordered_ID),])))

peak_matrix_ordered <- as.data.frame(cbind(rbind(Branch2_mat_2_top1_left[Prebranch_mat_2_top1_ordered_ID,],
                                                 Branch2_mat_2_top2_left[Branch1_mat_2_top2_ordered_ID,],
                                                 Branch2_mat_2_top3_left[rev(Branch2_mat_2_top3_ordered_ID),]),
                                     peak_matrix_ordered))

pdf("Epi_roothair_branch_peaks_heatmap.pdf", height = 7.5, width = 6)
Heatmap(as.matrix(peak_matrix_ordered), show_row_names = F,
        cluster_rows = F, use_raster = T, raster_resize_mat = F,
        show_column_names = F, cluster_columns = F,
        raster_device = "png", raster_quality = 5, raster_by_magick = F,
        col = colorRampPalette(c("#000080","#3A5FCD","#87CEEB", "white", "#F08080","red","#B51F29"))(100),
        name = paste0(nrow(peak_matrix_ordered), " Peaks"),
        row_split = c(rep("Pre-branch", length(markerPKList_Branch_peak_id$`Pre-branch`$ID)),
                      rep("Branch 1", length(markerPKList_Branch_peak_id$`Branch 1`$ID)),
                      rep("Branch 2", length(markerPKList_Branch_peak_id$`Branch 2`$ID))))
 dev.off()



# TF motif
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
  
  markersPK_Branch <- getMarkerFeatures(
    ArchRProj = Root_ATAC, 
    useMatrix = "PeakMatrix", 
    groupBy = "Branch",
    bias = c("log10(nFrags)", "TSSEnrichment"),
    testMethod = "wilcoxon"
  )
  markerPKList_Branch <- getMarkers(markersPK_Branch, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  
  marker_peak_id <- c()
  for (i in names(markerPKList_Branch)) {
    temp <- markerPKList_Branch[[i]]
    temp$ID <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
    marker_peak_id <- c(marker_peak_id, temp$ID)
    markerPKList_Branch[[i]] <- temp
  }
  
  sum(duplicated(marker_peak_id))
  duplicated_peaks <- marker_peak_id[duplicated(marker_peak_id)]
  
  for (i in names(markerPKList_Branch)) {
    temp <- markerPKList_Branch[[i]]
    temp <- temp[-which(temp$ID %in% duplicated_peaks),]
    markerPKList_Branch[[i]] <- temp
  }
  
  peaksets <- getPeakSet(Root_ATAC)
  peaksets@ranges@NAMES <- NULL
  peaksets <- as.data.frame(peaksets)
  peaksets$ID <- paste0(peaksets$seqnames, "_", peaksets$start, "_", peaksets$end)
  rownames(peaksets) <- peaksets$ID
  
  for (i in names(markerPKList_Branch)) {
    temp <- markerPKList_Branch[[i]]
    temp$peakType <- peaksets[temp$ID, "peakType"]
    temp <- temp[which(temp$peakType %in% c("Promoter",
                                            "Distal")),]
    markerPKList_Branch[[i]] <- temp
  }
  
  Root_ATAC <- addMotifAnnotations(ArchRProj = Root_ATAC, motifPWMs = motif_all,
                                   annoName = "TF-Motif", force = T)
  # enrichMotifs <- peakAnnoEnrichment(
  #   seMarker = markerPKList_Branch,
  #   ArchRProj = Root_ATAC,
  #   peakAnnotation = "TF-Motif",
  #   cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
  # )
  
  Branch2_Peaks_PD <- markerPKList_Branch$`Branch 2`$ID
  Branch1_Peaks_PD <- markerPKList_Branch$`Branch 1`$ID
  Prebranch_Peaks_PD <- markerPKList_Branch$`Pre-branch`$ID
  
  matches <- getMatches(Root_ATAC, name = "TF-Motif")
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  rownames(matches) <- pr1
  
  peaksets_2_pd <- peaksets[c(Prebranch_Peaks_PD, Branch1_Peaks_PD, Branch2_Peaks_PD),]
  
  temp <- matrix(FALSE, nrow = nrow(peaksets_2_pd), ncol = 3)
  rownames(temp) <- peaksets_2_pd$ID
  colnames(temp) <- c("Pre-branch", "Branch 1", "Branch 2")
  temp[Prebranch_Peaks_PD, 1] <- TRUE
  temp[Branch1_Peaks_PD, 2] <- TRUE
  temp[Branch2_Peaks_PD, 3] <- TRUE
  
  matches <- matches[which(rownames(matches) %in% rownames(temp)),]
  matches <- matches[rownames(temp),]
  identical(rownames(temp), rownames(matches))
  
  Prebranch_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,1]), 1:nrow(matches))
  Prebranch_Peaks_Enrich$Enrichment_log <- log2(Prebranch_Peaks_Enrich$Enrichment)
  Prebranch_Peaks_Enrich$Enrichmented <- "NO"
  Prebranch_Peaks_Enrich$Enrichmented[which(Prebranch_Peaks_Enrich$Enrichment_log >= 0.25 & Prebranch_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Prebranch_Peaks_Enrich$Enrichmented)
  
  colnames(temp)
  Branch1_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,2]), 1:nrow(matches))
  Branch1_Peaks_Enrich$Enrichment_log <- log2(Branch1_Peaks_Enrich$Enrichment)
  Branch1_Peaks_Enrich$Enrichmented <- "NO"
  Branch1_Peaks_Enrich$Enrichmented[which(Branch1_Peaks_Enrich$Enrichment_log >= 0.25 & Branch1_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Branch1_Peaks_Enrich$Enrichmented)
  
  colnames(temp)
  Branch2_Peaks_Enrich <- .computeEnrichment(matches, which(temp[,3]), 1:nrow(matches))
  Branch2_Peaks_Enrich$Enrichment_log <- log2(Branch2_Peaks_Enrich$Enrichment)
  Branch2_Peaks_Enrich$Enrichmented <- "NO"
  Branch2_Peaks_Enrich$Enrichmented[which(Branch2_Peaks_Enrich$Enrichment_log >= 0.25 & Branch2_Peaks_Enrich$mlog10Padj >= 5)] <- "YES"
  table(Branch2_Peaks_Enrich$Enrichmented)
  
  gplot1 <- function(DATA) {
    DATA$Enrichmented <- factor(DATA$Enrichmented,
                                levels = c("YES", 'NO'))
    temp <- NULL
    if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") > 5) {
      temp <- DATA[which(DATA$Enrichmented == "YES"),]
      temp <- temp[order(temp$Enrichment_log, decreasing = T),]
      temp <- temp[1:5,]
    }
    if (sum(DATA$Enrichmented == "YES") > 0 & sum(DATA$Enrichmented == "YES") <= 5) {
      temp <- DATA[which(DATA$Enrichmented == "YES"),]
      temp <- temp[order(temp$Enrichment_log, decreasing = T),]
    }
    
    g <- ggplot() +
      geom_point(data = DATA, aes(x = Enrichment_log, y = mlog10Padj, color = Enrichmented)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
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
  
  Prebranch_Peaks_Enrich_ggplot <- gplot1(Prebranch_Peaks_Enrich)
  Prebranch_Peaks_Enrich_ggplot
  table(Prebranch_Peaks_Enrich$Enrichmented)
  Branch1_Peaks_Enrich_ggplot <- gplot1(Branch1_Peaks_Enrich)
  Branch1_Peaks_Enrich_ggplot
  table(Branch1_Peaks_Enrich$Enrichmented)
  Branch2_Peaks_Enrich_ggplot <- gplot1(Branch2_Peaks_Enrich)
  Branch2_Peaks_Enrich_ggplot
  table(Branch2_Peaks_Enrich$Enrichmented)
  
  pdf("Epi_roothair_TF_enrich.pdf", width = 10, height = 3)
  Prebranch_Peaks_Enrich_ggplot +
    Branch1_Peaks_Enrich_ggplot +
    Branch2_Peaks_Enrich_ggplot  +
    plot_layout(nrow = 1, byrow = FALSE)
  dev.off()
  
  Prebranch_TFs <- Prebranch_Peaks_Enrich[which(Prebranch_Peaks_Enrich$Enrichmented == "YES"),]
  Branch1_TFs <- Branch1_Peaks_Enrich[which(Branch1_Peaks_Enrich$Enrichmented == "YES"),]
  Branch2_TFs <- Branch2_Peaks_Enrich[which(Branch2_Peaks_Enrich$Enrichmented == "YES"),]
  
  TFs_specific <- c(Prebranch_TFs$feature, Branch1_TFs$feature, Branch2_TFs$feature)
  
  duplicated_TF <- TFs_specific[duplicated(TFs_specific)]
  
  # TFs_specific <- TFs_specific[!duplicated(TFs_specific)]
  
  nrow(Branch2_Peaks_Enrich)
  rownames(Prebranch_Peaks_Enrich) <- Prebranch_Peaks_Enrich$feature
  rownames(Branch1_Peaks_Enrich) <- Branch1_Peaks_Enrich$feature
  rownames(Branch2_Peaks_Enrich) <- Branch2_Peaks_Enrich$feature
  enrich_heatmap <- matrix(NA, nrow = 3, ncol = 263)
  colnames(enrich_heatmap) <- Branch2_Peaks_Enrich$feature
  rownames(enrich_heatmap) <- c("Pre-branch", "Branch 1", "Branch 2")
  
  enrich_heatmap <- as.data.frame(t(enrich_heatmap))
  enrich_heatmap$`Pre-branch` <- Prebranch_Peaks_Enrich[rownames(enrich_heatmap), "mlog10Padj"]
  enrich_heatmap$`Branch 1` <- Branch1_Peaks_Enrich[rownames(enrich_heatmap), "mlog10Padj"]
  enrich_heatmap$`Branch 2` <- Branch2_Peaks_Enrich[rownames(enrich_heatmap), "mlog10Padj"]
  duplicated_TF
  Prebranch_TFs <- Prebranch_TFs[which(!c(Prebranch_TFs$feature %in% duplicated_TF)),]
  Branch1_TFs <- Branch1_TFs[which(!c(Branch1_TFs$feature %in% duplicated_TF)),]
  Branch2_TFs <- Branch2_TFs[which(!c(Branch2_TFs$feature %in% duplicated_TF)),]
  TFs <- c(Prebranch_TFs$feature,
           Branch1_TFs$feature,
           Branch2_TFs$feature)
  enrich_heatmap <- enrich_heatmap[TFs,]
  pdf("Epi_roothair_TFs_enrich_heatmap.pdf", width = 11, height = 2)
  pheatmap(t(enrich_heatmap), scale = "none", cluster_rows = F,
           cluster_cols = F, border_color = F, angle_col = "90", fontsize = 7)
  dev.off()
}

# TF Family
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
  
  TF_Family_temp <- TF_Family
  rownames(TF_Family_temp) <- TF_Family_temp$RAP_ID
  
  rownames(Annotation_genes) <- Annotation_genes$`CGSNL Gene Symbol`
  Prebranch_TFs_ID <- Annotation_genes[Prebranch_TFs$feature, "Locus_ID"]
  Branch1_TFs_ID <- Annotation_genes[Branch1_TFs$feature, "Locus_ID"]
  Branch2_TFs_ID <- Annotation_genes[Branch2_TFs$feature, "Locus_ID"]
  
  TF_Family_temp <- TF_Family_temp[which(TF_Family_temp$RAP_ID %in% c(Prebranch_TFs_ID, Branch1_TFs_ID, Branch2_TFs_ID)),]
  TF_Family_temp$Stage <- ""
  TF_Family_temp$Stage[which(TF_Family_temp$RAP_ID %in% Prebranch_TFs_ID)] <- "Pre-branch"
  TF_Family_temp$Stage[which(TF_Family_temp$RAP_ID %in% Branch1_TFs_ID)] <- "Branch 1"
  TF_Family_temp$Stage[which(TF_Family_temp$RAP_ID %in% Branch2_TFs_ID)] <- "Branch 2"
  TF_Family_temp <- TF_Family_temp[order(TF_Family_temp$Stage),]
  table(TF_Family_temp$Stage)
  
  
  TF_family_stage <- as.data.frame.array(table(TF_Family_temp$TF_Family, TF_Family_temp$Stage))
  TF_family_stage_p <- matrix(0, nrow = nrow(TF_family_stage), ncol = ncol(TF_family_stage))
  rownames(TF_family_stage_p) <- rownames(TF_family_stage)
  colnames(TF_family_stage_p) <- colnames(TF_family_stage)
  TF_family_number <- as.data.frame(table(TF_Family$TF_Family))
  TF_total_number <- 1750
  for (i in colnames(TF_family_stage_p)) {
    for (j in rownames(TF_family_stage_p)) {
      common_num <- TF_family_stage[j,i]
      p <- phyper(common_num-1, sum(TF_family_stage[,i]), TF_total_number-sum(TF_family_stage[,i]), as.numeric(table(TF_Family$TF_Family)[j]), lower.tail=F)
      TF_family_stage_p[j,i] <- p
    }
  }
  
  TF_family_stage <- TF_family_stage[order(TF_family_stage[,3],
                                           TF_family_stage[,2],
                                           TF_family_stage[,1], decreasing = T),]
  TF_family_stage_p <- TF_family_stage_p[rownames(TF_family_stage),]
  TF_family_stage_p <- TF_family_stage_p[!c(TF_family_stage_p[,1] > 0.05 & TF_family_stage_p[,2] > 0.05 & TF_family_stage_p[,3] > 0.05),]
  TF_family_stage <- TF_family_stage[rownames(TF_family_stage_p),]
  
  TF_order <- c("MYB",
                "NAC", "WRKY", "HD-ZIP", "YABBY",
                "ERF", "TCP")
  stage_order <- colnames(TF_family_stage)
  TF_family_stage_melt <- melt(data.frame(TF = rownames(TF_family_stage),
                                          TF_family_stage), id.vars = "TF")
  TF_family_stage_p_melt <- melt(data.frame(TF = rownames(TF_family_stage_p),
                                            TF_family_stage_p), id.vars = "TF")
  TF_family_stage_p_melt$TF <- factor(TF_family_stage_p_melt$TF, levels = TF_order)
  
  TF_family_stage_p_melt$P_if <- ""
  for (i in 1:nrow(TF_family_stage_p_melt)) {
    if (TF_family_stage_p_melt$value[i] < 0.05) {
      TF_family_stage_p_melt$P_if[i] <- "*"
    }
    if (TF_family_stage_p_melt$value[i] < 0.01) {
      TF_family_stage_p_melt$P_if[i] <- "**"
    }
    if (TF_family_stage_p_melt$value[i] < 0.001) {
      TF_family_stage_p_melt$P_if[i] <- "***"
    }
    if (TF_family_stage_p_melt$value[i] < 0.0001) {
      TF_family_stage_p_melt$P_if[i] <- "****"
    }
  }
  
  
  TF_family_stage_melt$P <- TF_family_stage_p_melt$value
  # TF_family_stage_melt$label <- paste0(TF_family_stage_melt$value, " (", TF_family_stage_p_melt$P_if, ")")
  TF_family_stage_melt$label <- TF_family_stage_p_melt$P_if
  # TF_family_stage_melt <- TF_family_stage_melt[which(TF_family_stage_melt$value > 0),]
  TF_family_stage_melt$P[which(TF_family_stage_melt$P > 0.05)] <- NA
  TF_family_stage_melt$TF <- factor(TF_family_stage_melt$TF, levels = TF_order)
  pdf("Epi_roothair_TFs_family.pdf", width = 4, height = 4.5)
  ggplot() + 
    theme_bw() +
    geom_tile(data = TF_family_stage_melt, aes(x = variable, y = TF, fill = -log10(P))) +
    geom_text(data = TF_family_stage_melt, aes(x = variable, y = TF, label = label),
              size = 8) +
    scale_fill_gradient2(low = "gray40", mid = "gray80", high = "red") +
    labs(y = "TF family", x = "Trajectory stage", fill = "-Log10(P-value)") +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          panel.grid = element_blank()) +
    scale_x_discrete(expand = c(0,0))
  dev.off()
  
  TF_family_stage_melt_sig <- TF_family_stage_melt[which(TF_family_stage_melt$P < 0.05),]
  TF_family_stage_melt_notsig <- TF_family_stage_melt[which(TF_family_stage_melt$P > 0.05),]
  
  pdf("Figure4_sup_2_TF_family.pdf", width = 5, height = 7)
  ggplot() + 
    theme_bw() +
    geom_tile(data = TF_family_stage_p_melt, aes(x = variable, y = TF, fill = -log10(value))) +
    geom_text(data = TF_family_stage_melt_sig, aes(x = variable, y = TF, label = label),
              size = 4) +
    geom_text(data = TF_family_stage_melt_notsig, aes(x = variable, y = TF, label = value),
              size = 4) +
    scale_fill_gradient2(midpoint = -log10(0.05), low = "gray40", mid = "gray80", high = "red") +
    labs(y = "TF family", x = "Trajectory stage", fill = "-Log10(P-value)") +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          panel.grid = element_blank()) +
    scale_x_discrete(expand = c(0,0))
  dev.off()
  
  
}

# WRKY family
{
  table(Root_ATAC$Branch)
  Root_ATAC <- addGroupCoverages(ArchRProj = Root_ATAC,
                                 groupBy = "Branch", force = T)
  motifPositions <- getPositions(Root_ATAC, name = "TF-Motif")
  seFoot_Epi_root_hair <- getFootprints(
    ArchRProj = Root_ATAC, 
    positions = motifPositions, 
    groupBy = "Branch"
  )
  unique(Root_ATAC$Branch)
  plotFootprints(
    seFoot = seFoot_Epi_root_hair,
    ArchRProj = Root_ATAC, 
    normMethod = "Subtract",
    plotName = "Epi_root_hair_TF_branch",
    addDOC = FALSE,
    smoothWindow = 5,
    pal = c("#698EAF", "#F27C5F", "#DB7577")
  )
  
  Root_ATAC <- addBgdPeaks(Root_ATAC)
  Root_ATAC <- addDeviationsMatrix(
    ArchRProj = Root_ATAC, 
    peakAnnotation = "TF-Motif",
    force = TRUE
  )
  
  WRKY <- c("OsWRKY16", "OsWRKY71", "OsWRKY67", "OsWRKY51", "OsWRKY24", "OsWRKY96", "OsWRKY1",
            "OsWRKY37", "OsWRKY39", "OsWRKY72", "OsWRKY13", "OsWRKY3", "OsWRKY29", "OsWRKY23",
            "OsWRKY45", "OsWRKY43", "OsWRKY48", "OsWRKY87", "OsWRKY74", "OsWRKY28", "OsWRKY36",
            "OsWRKY30")
  WRKY <- data.frame(WRKY = WRKY,
                     Symbol = Annotation_genes[WRKY, "DataSets_Symbol"])
  
  
  markersGS_Branch <- getMarkerFeatures(
    ArchRProj = Root_ATAC, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Branch",
    bias = c("log10(nFrags)", "TSSEnrichment"),
    testMethod = "wilcoxon"
  )
  markerGSList_Branch <- getMarkers(markersGS_Branch, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  
  markerGSList_Branch_gene <- c()
  for (i in names(markerGSList_Branch)) {
    temp <- markerGSList_Branch[[i]]
    markerGSList_Branch_gene <- c(markerGSList_Branch_gene, temp@listData[["name"]])
  }
  sum(duplicated(markerGSList_Branch_gene))
  
  genes_duplicated <- markerGSList_Branch_gene[duplicated(markerGSList_Branch_gene)]
 
  # sum(WRKY$Symbol %in% rownames(Root_ATAC_GS))
  # Branch2_DEGs <- as.data.frame(markerGSList_Branch$`Branch 2`)
  # rownames(Branch2_DEGs) <- Branch2_DEGs$name
  # WRKY_filtered_2 <- Branch2_DEGs[WRKY$Symbol,]
  WRKY_filtered_2 <- as.data.frame(na.omit(WRKY_filtered_2))
  
  genes_meta <- as.data.frame(gene_annotation@listData[["genes"]])
  rownames(genes_meta) <- genes_meta$symbol
  genes_meta$matrix_id <- paste0(genes_meta$seqnames, ":", genes_meta$symbol)
  WRKY$matrix_id <- genes_meta[WRKY$Symbol, "matrix_id"]
  Root_ATAC_sub <- Root_ATAC[which(Root_ATAC$Branch %in% c("Pre-branch", "Branch 2")),]
  trajectory_GS <- getTrajectory(Root_ATAC_sub, name = "Pseudotime", smoothWindow = 5,
                                 useMatrix = "GeneScoreMatrix", Pseudotime_max = ceiling(max(Root_ATAC_sub$Pseudotime)))
  trajectory_GS_top_var <- plotTrajectoryHeatmap(trajectory_GS,  pal = paletteContinuous(set = "horizonExtra"), 
                                                 varCutOff = 0.9, returnMatrix = TRUE)
  
  WRKY_filtered <- WRKY[WRKY$matrix_id %in% rownames(trajectory_GS_top_var),]
  
  min(Root_ATAC$Pseudotime)

  ht1 <- plotTrajectoryHeatmap(trajectory_GS[WRKY_filtered$matrix_id,],  pal = paletteContinuous(set = "horizonExtra"), 
                               varCutOff = 0, idxLabel_custom = rownames(trajectory_GS), rowOrder = WRKY_filtered$matrix_id)
  TF_expression <- plotTrajectoryHeatmap(trajectory_GS[WRKY_filtered$matrix_id,],  pal = paletteContinuous(set = "horizonExtra"), 
                                         varCutOff = 0, idxLabel_custom = rownames(trajectory_GS), rowOrder = WRKY_filtered$matrix_id,
                                         returnMatrix = T)
  TF_expression <- colMeans(TF_expression)
  TF_expression <- ( TF_expression - min(TF_expression) ) / ( max(TF_expression) - min(TF_expression) )
  
  trajectory_MM <- getTrajectory(Root_ATAC_sub, name = "Pseudotime", smoothWindow = 10,
                                 useMatrix = "TF-MotifMatrix", Pseudotime_max = ceiling(max(Root_ATAC$Pseudotime)))
  ht2 <- plotTrajectoryHeatmap(trajectory_MM[paste0("z:", WRKY_filtered$WRKY),], pal = paletteContinuous(set = "solarExtra"),
                               varCutOff = 0, idxLabel_custom = rownames(trajectory_MM[paste0("z:", WRKY_filtered$WRKY),]),
                               rowOrder = paste0("z:", WRKY_filtered$WRKY), limits = c(-0.8, 0.215))
  ht2
  Motif_score <- plotTrajectoryHeatmap(trajectory_MM[paste0("z:", WRKY_filtered$WRKY),], pal = paletteContinuous(set = "solarExtra"),
                                       varCutOff = 0, idxLabel_custom = rownames(trajectory_MM[paste0("z:", WRKY_filtered$WRKY),]),
                                       rowOrder = paste0("z:", WRKY_filtered$WRKY), returnMatrix = T)
  Motif_score <- colMeans(Motif_score)
  Motif_score <- ( Motif_score - min(Motif_score) ) / ( max(Motif_score) - min(Motif_score) )
  
  
  pdf("Epi_roothair_WRKY_heatmap.pdf", width = 10, height = 5)
  ht1 + ht2
  dev.off()
  
  peaksets <- getPeakSet(Root_ATAC)
  motifPositions <- motifmatchr::matchMotifs(
    pwms = motif_all,
    subject = peaksets,
    genome = BSgenome.OSativa.NCBI.IRGSPv1.0, 
    out = "matches", 
    p.cutoff = 5e-05, 
    w = 7
  )
  peaksets@ranges@NAMES <- NULL
  peaksets_df <- as.data.frame(peaksets)
  peaksets_df$peak_id <- paste0(peaksets_df$seqnames, "_", peaksets_df$start, "_", peaksets_df$end)
  rownames(peaksets_df) <- peaksets_df$peak_id
  motifPositions@assays@data@listData[["motifMatches"]]@Dimnames[[2]] <- colnames(motifPositions)
  motifPositions@assays@data@listData[["motifMatches"]]@Dimnames[[1]] <- peaksets_df$peak_id
  motifPositions <- motifPositions@assays@data@listData[["motifMatches"]]
  motifPositions <- motifPositions + 0
  motifPositions_WRKY <- motifPositions[, WRKY_filtered$WRKY]
  motifPositions_WRKY <- motifPositions_WRKY[rowSums(motifPositions_WRKY) > 0,]
  
  peaks_WRKY <- peaksets_df[rownames(motifPositions_WRKY),]
  peaks_WRKY <- peaks_WRKY[which(peaks_WRKY$peakType == "Promoter"),]
  
  motifPositions_WRKY <- motifPositions_WRKY[peaks_WRKY$peak_id,]
  
  WRKY_list_gene <- c() 
  for (i in colnames(motifPositions_WRKY)) {
    temp <- rownames(motifPositions_WRKY)[which(motifPositions_WRKY[,i] > 0)]
    temp_genes <- peaks_WRKY[temp, "nearestGene"]
    WRKY_list_gene <- c(WRKY_list_gene, list(temp_genes))
  }
  names(WRKY_list_gene) <- colnames(motifPositions_WRKY)
  
  Epi_root_hair_GS <- getMatrixFromProject(Root_ATAC_sub, useMatrix = "GeneScoreMatrix")
  
  trajGS  <- getTrajectory(Root_ATAC_sub, name = "Pseudotime", smoothWindow = 5,
                           useMatrix = "GeneScoreMatrix", Pseudotime_max = ceiling(max(Root_ATAC_sub$Pseudotime)))
  smoothMat <- trajGS@assays@data@listData[["smoothMat"]]
  rownames(smoothMat) <- trajGS@elementMetadata@listData[["name"]]
  WRKY_target_genes_gs <- c()
  for (i in names(WRKY_list_gene)) {
    temp_genes <- WRKY_list_gene[[i]]
    temp <- smoothMat[temp_genes,]
    temp <- colSums(temp)
    temp <- temp / length(temp_genes)
    WRKY_target_genes_gs <- as.data.frame(cbind(WRKY_target_genes_gs,
                                                temp))
  }
  colnames(WRKY_target_genes_gs) <- names(WRKY_list_gene)
  pheatmap(t(WRKY_target_genes_gs), cluster_rows = F, cluster_cols = F)
  rowOrder = paste0("z:", WRKY_filtered$WRKY)
  WRKY_target_genes_gs <- t(WRKY_target_genes_gs)
  
  mat <- sweep(WRKY_target_genes_gs - rowMeans(WRKY_target_genes_gs), 1, matrixStats::rowSds(WRKY_target_genes_gs), `/`)
  mat[mat > 1.5] <- 1.5
  mat[mat < -1.5] <- -1.5
  
  pdf("Epi_roothair_WRKY_target_genes_gs.pdf", width = 5, height = 5)
  ComplexHeatmap::Heatmap(mat, cluster_rows = F, cluster_columns = F, show_column_names = F,
                          col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), border = NA,
                          )
  dev.off()
  
  Target_gene_mean <- colMeans(mat)
  Target_gene_mean <- ( Target_gene_mean - min(Target_gene_mean) ) / ( max(Target_gene_mean) - min(Target_gene_mean) )
  
  DF <- data.frame(Pseudotime = 1:length(TF_expression),
                   TF = TF_expression,
                   Motif = Motif_score,
                   Target = Target_gene_mean)
  DF$Pseudotime <- as.character(DF$Pseudotime)
  DF <- melt(DF, id.vars = "Pseudotime")
  DF$Pseudotime <- as.numeric(DF$Pseudotime)
  
  library(ggalt)
  pdf("Epi_roothair_TF_target_motif.pdf", width = 4, height = 3)
  ggplot(data = DF, aes(x = Pseudotime, y = value, group = variable, color = variable)) +
    # geom_point() +
    geom_xspline(spline_shape = 1.5) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw()
  dev.off()
  
  p <- plotBrowserTrack(
    ArchRProj = Root_ATAC, 
    groupBy = "Branch", 
    geneSymbol = c("OsWRKY72", "OsWRKY37"), 
    upstream = 10000,
    downstream = 10000,
    groupPal = c("#F27C5F", "#DB7577", "#698EAF"),
    baseSize = 8, facetbaseSize = 8
  )
  
  dev.off()
  pdf("Epi_roothair_genome_browse_OsWRKY72.pdf", width = 6, height = 2.5)
  grid::grid.newpage()
  grid::grid.draw(p$OsWRKY72)
  dev.off()
  
  pdf("Epi_roothair_genome_browse_OsWRKY37.pdf", width = 6, height = 2.5)
  grid::grid.newpage()
  grid::grid.draw(p$OsWRKY37)
  dev.off()
  
  dev.off()
  
  
  WRKY_list_gene
  
  
  
  genes_list <- WRKY_list_gene
  rownames(gene_id_map) <- gene_id_map$symbol
  for (i in names(genes_list)) {
    temp <- genes_list[[i]]
    temp <- gene_id_map[temp, "gene_id"]
    temp <- unique(temp)
    genes_list[[i]] <- temp
  }
  temp <- unlist(genes_list)
  temp <- sort(table(temp), decreasing = T)
  temp <- as.data.frame(temp)
  # temp <- temp[which(temp[,2] > 13),]
  genes_list <- list("WRKY" = temp$temp)
  GO_result_WRKY <- c()
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
    GO_result$TF <- i
    GO_result_WRKY <- as.data.frame(rbind(GO_result_WRKY, GO_result))
  }
  
  GO_result_WRKY <- GO_result_WRKY[which(GO_result_WRKY$Ontology == "BP"),]
  GO_result_WRKY_2 <- GO_result_WRKY
  GO_result_WRKY <- GO_result_WRKY[which(GO_result_WRKY$P_value < 0.05),]
  
  
  temp <- GO_result_WRKY[grep(pattern = "development|root|morphogenesis|differentiation|divition|organization", GO_result_WRKY$Description),]
  openxlsx::write.xlsx(GO_result_State, "Epi_root_hair_branch_GO.xlsx")
  
  
  
  
}

Pseudotime <- data.frame(Cells = rownames(Root_ATAC),
                         Pseudotime = Root_ATAC@cellColData[, c("Pseudotime")],
                         Branch = Root_ATAC$Branch,
                         Celltype = Root_ATAC$Celltype)
Pseudotime$Celltype <- as.character(Root_ATAC@cellColData[Pseudotime$Cells, "Celltype"])
pdf("Epi_roothair_density.pdf", height = 1.5, width = 5)
# ggplot() + 
#   geom_density(data = Pseudotime, aes(x = Pseudotime, color = Celltype),
#                alpha = 1, linewidth = 1, adjust = 0.5, position = "stack") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10)) +
#   scale_color_manual(values = c("#FF8888", "#c85c6c")) +
#   labs(x = "Pseudotime", y = "Density", fill = "Region") +
#   scale_x_continuous(expand = c(0,0))

ggplot() + 
  geom_density(data = Pseudotime, aes(x = Pseudotime, fill = Branch),
               alpha = 1, linewidth = 1, adjust = 0.5, position = "fill") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_fill_manual(values = c("#F27C5F", "#DB7577", "#698EAF")) +
  labs(x = "Pseudotime", y = "Density", fill = "Region") +
  scale_x_continuous(expand = c(0,0))

dev.off()


###############
###############
###############

if (FALSE) {
  
  table(Root_ATAC_sub$State)
  Root_ATAC_sub_peak_mat <- getMatrixFromProject(Root_ATAC_sub,
                                                 useMatrix = "PeakMatrix")
  temp <- as.data.frame(Root_ATAC_sub_peak_mat@rowRanges)
  temp$ID <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
  Root_ATAC_sub_peak_mat <- Root_ATAC_sub_peak_mat@assays@data@listData[["PeakMatrix"]]
  rownames(Root_ATAC_sub_peak_mat) <- temp$ID
  
  Prebranch_mat <- Root_ATAC_sub_peak_mat[markerPKList_Branch_peak_id$`Pre-branch`$ID, ]
  Prebranch_mat <- t(Prebranch_mat)
  Prebranch_mat <- as.matrix(Prebranch_mat)
  Prebranch_mat <- as.data.frame(Prebranch_mat)
  Pseudotime <- data.frame(Cells = rownames(Prebranch_mat),
                           Pseudotime = Root_ATAC@cellColData[rownames(Prebranch_mat), c("Pseudotime")])
  library(infotheo)
  nbins <- 100
  equal_width <- discretize(Pseudotime$Pseudotime, "equalwidth", nbins)
  Pseudotime$bin <- equal_width$X
  Pseudotime <- Pseudotime[order(Pseudotime$Pseudotime, Pseudotime$bin, decreasing = F),]
  Prebranch_mat_2 <- c()
  Prebranch_mat <- Prebranch_mat[Pseudotime$Cells, ]
  for (i in unique(Pseudotime$bin)) {
    print(i)
    temp <- Prebranch_mat[which(Pseudotime$bin == i), , drop = F]
    temp <- colMeans(temp)
    temp <- as.data.frame(temp)
    colnames(temp) <- i
    temp <- as.data.frame(t(temp))
    Prebranch_mat_2 <- as.data.frame(rbind(Prebranch_mat_2, temp))
  }
  Prebranch_mat_2 <- t(Prebranch_mat_2)
  Prebranch_mat_2 <- Prebranch_mat_2[, as.character(sort(as.numeric(colnames(Prebranch_mat_2))))]
  Prebranch_mat_2 <- sweep(Prebranch_mat_2 - rowMeans(Prebranch_mat_2), 1, matrixStats::rowSds(Prebranch_mat_2), `/`)
  Prebranch_mat_2[Prebranch_mat_2 > 3] <- 3
  Prebranch_mat_2[Prebranch_mat_2 < -2] <- -2
  idx <- order(apply(Prebranch_mat_2, 1, which.max))
  Prebranch_mat_2 <- Prebranch_mat_2[idx, ]
  pdf("Epi_roothair_peaks_prebranch_heatmap.pdf", width = 5, height = 6)
  Heatmap(Prebranch_mat_2, show_row_names = F,
          cluster_rows = F, use_raster = T, raster_resize_mat = F,
          show_column_names = F, cluster_columns = F,
          raster_device = "png", raster_quality = 5, raster_by_magick = F,
          col = paletteContinuous(set = "solarExtra", n = 100),
          name = paste0(nrow(Prebranch_mat_2), " Peaks"))
  dev.off()
  
  
  
  Branch2_mat <- Root_ATAC_sub_peak_mat[markerPKList_Branch_peak_id$`Branch 2`$ID, ]
  Branch2_mat <- t(Branch2_mat)
  Branch2_mat <- as.matrix(Branch2_mat)
  Branch2_mat <- as.data.frame(Branch2_mat)
  Pseudotime <- data.frame(Cells = rownames(Branch2_mat),
                           Pseudotime = Root_ATAC@cellColData[rownames(Branch2_mat), c("Pseudotime")])
  library(infotheo)
  nbins <- 100
  equal_width <- discretize(Pseudotime$Pseudotime, "equalwidth", nbins)
  Pseudotime$bin <- equal_width$X
  Pseudotime <- Pseudotime[order(Pseudotime$Pseudotime, Pseudotime$bin, decreasing = F),]
  Branch2_mat_2 <- c()
  Branch2_mat <- Branch2_mat[Pseudotime$Cells, ]
  for (i in unique(Pseudotime$bin)) {
    print(i)
    temp <- Branch2_mat[which(Pseudotime$bin == i), , drop = F]
    temp <- colMeans(temp)
    temp <- as.data.frame(temp)
    colnames(temp) <- i
    temp <- as.data.frame(t(temp))
    Branch2_mat_2 <- as.data.frame(rbind(Branch2_mat_2, temp))
  }
  Branch2_mat_2 <- t(Branch2_mat_2)
  Branch2_mat_2 <- Branch2_mat_2[, as.character(sort(as.numeric(colnames(Branch2_mat_2))))]
  Branch2_mat_2 <- sweep(Branch2_mat_2 - rowMeans(Branch2_mat_2), 1, matrixStats::rowSds(Branch2_mat_2), `/`)
  Branch2_mat_2[Branch2_mat_2 > 3] <- 3
  Branch2_mat_2[Branch2_mat_2 < -2] <- -2
  idx <- order(apply(Branch2_mat_2, 1, which.max))
  Branch2_mat_2 <- Branch2_mat_2[idx, ]
  pdf("Epi_roothair_peaks_branch2_heatmap.pdf", width = 5, height = 6)
  Heatmap(Branch2_mat_2, show_row_names = F,
          cluster_rows = F, use_raster = T, raster_resize_mat = F,
          show_column_names = F, cluster_columns = F,
          raster_device = "png", raster_quality = 5, raster_by_magick = F,
          col = paletteContinuous(set = "solarExtra", n = 100),
          name = paste0(nrow(Branch2_mat_2), " Peaks"))
  dev.off()
  
  Pseudotime <- data.frame(Cells = rownames(Prebranch_mat),
                           Pseudotime = Root_ATAC@cellColData[rownames(Prebranch_mat), c("Pseudotime")])
  library(infotheo)
  nbins <- 100
  equal_width <- discretize(Pseudotime$Pseudotime, "equalwidth", nbins)
  Pseudotime$bin <- equal_width$X
  Pseudotime <- as.data.frame(Pseudotime)
  Pseudotime$Celltype <- as.character(Root_ATAC@cellColData[Pseudotime$Cells, "Celltype"])
  pdf("Epi_roothair_density.pdf", height = 1.5, width = 5)
  ggplot() + 
    geom_density(data = Pseudotime, aes(x = Pseudotime, color = Celltype),
                 alpha = 1, linewidth = 1, adjust = 0.5, position = "stack") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    scale_color_manual(values = c("#FF8888", "#c85c6c")) +
    labs(x = "Pseudotime", y = "Density", fill = "Region") +
    scale_x_continuous(expand = c(0,0))
  dev.off()
  
  
  
  
  
  
  
  sc_cds$Os10g0452700 <- log2(seu@assays$RNA@counts["Os10g0452700", colnames(sc_cds)] + 1)
  plot_cell_trajectory(sc_cds, color_by = "Os10g0452700") +
    scale_color_gradient(low = "gray", high = "red")
  
  p <- plotEmbedding(
    ArchRProj = Root_ATAC, 
    colorBy = "GeneScoreMatrix", 
    name = c("Os10g0452700"), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(Root_ATAC)
  )
  p
  
  plotEmbedding(ArchRProj = Root_ATAC, colorBy = "cellColData",
                name = "Clusters", embedding = "UMAP_Harmony")
  
  
  getAvailableMatrices(Root_ATAC)
  Root_ATAC_GS <- getMatrixFromProject(project, useMatrix = "GeneScoreMatrix")
  Root_ATAC_GS@assays@data@listData[["GeneScoreMatrix"]]@Dimnames[[1]] <- Root_ATAC_GS@elementMetadata@listData[["name"]]
  Root_ATAC_GS <- Root_ATAC_GS@assays@data@listData[["GeneScoreMatrix"]]
  
  Archr_FeaturePlot <- function(project = Root_ATAC, matrix_use,
                                genes, reduction, dim_names) {
    embedding <- getEmbedding(Root_ATAC, embedding = reduction)
    colnames(embedding) <- dim_names
    exp <- matrix_use[genes, , drop = FALSE]
    exp <- as.matrix(t(exp), drop = FALSE)
    Data <- as.data.frame(cbind(embedding,
                                exp[rownames(embedding),,drop = FALSE]))
    ggplot(data = Data, aes(x = Data[,1], y = Data[,2], color = log2(Data[,3]))) +
      geom_point(size = 1) +
      scale_color_manual(values = celltype_color) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.line = element_line(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
  }
  
  UMAP_withBatchEffect <- getEmbedding(Root_ATAC, embedding = "UMAP_withBatchEffect")
  colnames(UMAP_withBatchEffect) <- c("UMAP_1", "UMAP_2")
  ColData <- as.data.frame(getCellColData(Root_ATAC))
  UMAP_withBatchEffect <- as.data.frame(cbind(UMAP_withBatchEffect,
                                              ColData))
  pdf("Figure4_补充_epi_roothair_Celltype_UMAP.pdf", width = 6.9, height = 6)
  ggplot(data = UMAP_withBatchEffect, aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
    geom_point(size = 1) +
    scale_color_manual(values = celltype_color) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggplot(data = UMAP_withBatchEffect, aes(x = UMAP_1, y = UMAP_2, color = TissuesSub)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggplot(data = UMAP_withBatchEffect, aes(x = UMAP_1, y = UMAP_2, color = Clusters)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  dev.off()
  
  
  
  
  
  getAvailableMatrices(Root_ATAC)
  
  Root_ATAC_GS <- getMatrixFromProject(Root_ATAC, useMatrix = "GeneScoreMatrix")
  Root_ATAC_GS@assays@data@listData[["GeneScoreMatrix"]]@Dimnames[[1]] <- Root_ATAC_GS@elementMetadata@listData[["name"]]
  Root_ATAC_GS_colData <- as.data.frame(Root_ATAC_GS@colData@listData)
  rownames(Root_ATAC_GS_colData) <- Root_ATAC_GS@colData@rownames
  Root_ATAC_GS <- Root_ATAC_GS@assays@data@listData[["GeneScoreMatrix"]]
  
  Root_ATAC_GS_seu <- CreateSeuratObject(Root_ATAC_GS, meta.data = Root_ATAC_GS_colData)
  Root_ATAC_GS_seu <- NormalizeData(Root_ATAC_GS_seu, scale.factor = 100000)
  DotPlot(Root_ATAC_GS_seu, features = c("Os10g0452700"), group.by = "Clusters")
  
  Root_ATAC_GS_seu %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = c("nFrags", "TSSEnrichment")) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    RunTSNE() -> Root_ATAC_GS_seu
  
  library(monocle)
  
  FeaturePlot(Root_ATAC_GS_seu, features = "Os10g0452700", order = T,
              pt.size = 1)
  
  DimPlot(Root_ATAC_GS_seu, group.by = "Celltype")
  
  UMAP_withBatchEffect <- getEmbedding(Root_ATAC, embedding = "UMAP_withBatchEffect")
  colnames(UMAP_withBatchEffect) <- c("UMAP_1", "UMAP_2")
  ColData <- as.data.frame(getCellColData(Root_ATAC))
  UMAP_withBatchEffect <- as.data.frame(cbind(UMAP_withBatchEffect,
                                              ColData))
  pdf("Figure4_补充_10dRoot_Celltype_UMAP.pdf", width = 6.9, height = 6)
  ggplot(data = UMAP_withBatchEffect, aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
    geom_point(size = 1) +
    scale_color_manual(values = celltype_color) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  dev.off()
  
  unique(UMAP_withBatchEffect$Celltype)
  trajectory <- c("Epidermis", "Root hair")
  Root_ATAC <- addTrajectory(
    ArchRProj = Root_ATAC, 
    name = "Epi_root_hair", 
    groupBy = "Celltype",
    trajectory = trajectory, 
    embedding = "UMAP_withBatchEffect", 
    force = TRUE
  )
  p <- plotTrajectory(Root_ATAC, trajectory = "Epi_root_hair", embedding = "UMAP_withBatchEffect",
                      colorBy = "cellColData", name = "Epi_root_hair", baseSize = 10, rastr = F,
                      continuousSet = "blueYellow")
  dev.off()
  pdf("Figure4_Epi_root_hair_trajectory.pdf", width = 6, height = 6)
  p[[1]]
  dev.off()
}

