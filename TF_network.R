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
RNA_celltype_markers <- FindAllMarkers(RNA, logfc.threshold = 0.25, only.pos = T)
RNA_celltype_markers_1 <- RNA_celltype_markers[which(RNA_celltype_markers$avg_log2FC > 1),]
RNA_celltype_markers_1 %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

DotPlot(RNA, features = c("Os04g0639100",
                          "Os05g0111300",
                          "Os01g0127600",
                          "Os10g0554800"), group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

Idents(RNA) <- RNA$tissues
RNA_tissue_markers <- FindAllMarkers(RNA, logfc.threshold = 0.25, only.pos = T)

RNA_imputated <- magic(RNA)
saveRDS(RNA_imputated, "RNA_imputated.rds")
RNA_imputated_matrix <- RNA_imputated@assays[["MAGIC_RNA"]]@data

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

TF_regulation <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/regulation_merged_Osj.txt",
                          sep = "\t", header = F)
table(TF_regulation$V5)
colnames(TF_regulation) <- c("TF_MSU", "Type", "Target_MSU", "Species", "Methods")
TF_regulation$TF_Target_MSU <- paste0(TF_regulation$TF_MSU ,"_" , TF_regulation$Target_MSU)
TF_regulation <- TF_regulation[!duplicated(TF_regulation$TF_Target_MSU),]

TF_regulation2 <- merge(x = TF_regulation, y = MSU_RAP, by.x = "TF_MSU", by.y = "MSU")
colnames(TF_regulation2)[ncol(TF_regulation2)] <- "TF_RAP"
TF_regulation2 <- merge(x = TF_regulation2, y = MSU_RAP, by.x = "Target_MSU", by.y = "MSU")
colnames(TF_regulation2)[ncol(TF_regulation2)] <- "Target_RAP"
TF_regulation <- as.data.frame(na.omit(TF_regulation2))
TF_regulation$Pairs <- paste0(TF_regulation$TF_RAP,"_",TF_regulation$Target_RAP)
sum(duplicated(TF_regulation$Pairs))
nrow(TF_regulation) - sum(duplicated(TF_regulation$Pairs))

TFs <- unique(TF_regulation$TF_RAP)
Targets <- setdiff(unique(TF_regulation$Target_RAP), unique(TF_regulation$TF_RAP))
TF_regulation_node_type <- data.frame(Nodes = c(TFs, Targets),
                                      Type = c(rep("TF",length(TFs)), rep("Target",length(Targets))))
rownames(TF_regulation_node_type) <- TF_regulation_node_type$Nodes

intersect(RNA_celltype_markers$gene, TF_regulation$TF_RAP)
intersect(RNA_celltype_markers$gene, TF_regulation$Target_RAP)

TF_regulation_list <- split.data.frame(TF_regulation, f = list(TF_regulation$TF_RAP))
Cor_result <- c()
for (i in names(TF_regulation_list)) {
  print(i)
  print(which(names(TF_regulation_list) == i))
  temp <- TF_regulation_list[[i]]
  temp_target <- temp$Target_RAP
  temp_target <- temp_target[which(temp_target %in% rownames(RNA_imputated_matrix))]
  temp_tf <- i
  if (temp_tf %in% rownames(RNA_imputated_matrix)) {
    temp_target <- t(RNA_imputated_matrix[temp_target,])
    temp_tf <- matrix(t(RNA_imputated_matrix[temp_tf,]))
    colnames(temp_tf) <- i
    rownames(temp_tf) <- colnames(RNA_imputated_matrix)
    temp_cor_result_dataframe <- c()
    for (j in 1:ncol(temp_target)) {
      print(j)
      temp_cor_test <- cor.test(temp_tf[,1], temp_target[,j])
      temp_cor <- temp_cor_test$estimate
      temp_cor_p <- temp_cor_test$p.value
      temp_cor_result <- data.frame(TF = i,
                                    Target = colnames(temp_target)[j],
                                    Pearson = temp_cor,
                                    P = temp_cor_p)
      temp_cor_result_dataframe <- as.data.frame(rbind(temp_cor_result_dataframe,
                                                       temp_cor_result))
    }
    Cor_result <- c(Cor_result, list(temp_cor_result_dataframe))
  }
}
Cor_result <- rbindlist(Cor_result)
Cor_result$FDR <- p.adjust(Cor_result$P, method = "fdr")
Cor_result$Pairs <- paste0(Cor_result$TF, "_", Cor_result$Target)
Cor_result2 <- Cor_result[!duplicated(Cor_result$Pairs),]

Cor_result_select <- Cor_result2[which(abs(Cor_result2$Pearson) > 0.65 & Cor_result2$FDR < 0.00001),]
which(Cor_result_select$Pearson < 0)
Cor_result_select <- as.data.frame(Cor_result_select)
rownames(Cor_result_select) <- Cor_result_select$Pairs
Cor_result_select <- Cor_result_select[which(Cor_result_select$Pearson > 0),]
saveRDS(Cor_result_select, "TF_regulation_filtered.rds")

RNA_celltype_markers_list <- split.data.frame(RNA_celltype_markers, f = factor(RNA_celltype_markers$cluster))
celltypes <- names(RNA_celltype_markers_list)
celltypes <- celltypes[-22]
celltypes <- sort(celltypes)
RNA_celltype_markers_list <- RNA_celltype_markers_list[celltypes]
RNA_celltype_markers_list_filtered <- RNA_celltype_markers_list
for (i in celltypes) {
  print(i)
  temp <- RNA_celltype_markers_list_filtered[[i]]
  temp <- temp[which(temp$gene %in% c(Cor_result_select$TF, Cor_result_select$Target)),]
  RNA_celltype_markers_list_filtered[[i]] <- temp
}
RNA_celltype_TF_network <- c()
for (i in celltypes) {
  print(i)
  temp <- RNA_celltype_markers_list_filtered[[i]]
  temp <- matrix(0, nrow = length(temp$gene), ncol = length(temp$gene))
  colnames(temp) <- RNA_celltype_markers_list_filtered[[i]]$gene
  rownames(temp) <- RNA_celltype_markers_list_filtered[[i]]$gene
  temp <- data.frame(gene = RNA_celltype_markers_list_filtered[[i]]$gene,
                     temp)
  temp <- reshape::melt.data.frame(temp, id.vars = "gene")
  temp$Pairs <- paste0(temp$gene, "_", temp$variable)
  rownames(temp) <- temp$Pairs
  temp <- temp[which(temp$Pairs %in% Cor_result_select$Pairs),]
  temp$Pearson <- Cor_result_select[temp$Pairs, "Pearson"]
  temp$FDR <- Cor_result_select[temp$Pairs, "FDR"]
  RNA_celltype_TF_network <- as.data.frame(rbind(RNA_celltype_TF_network,
                                                 data.frame(TF = temp$gene,
                                                            Target = temp$variable,
                                                            Pearson = temp$Pearson,
                                                            FDR = temp$FDR,
                                                            Pairs = temp$Pairs,
                                                            Celltype = i)))
}
sum(duplicated(RNA_celltype_TF_network$Pairs))
table(table(RNA_celltype_TF_network$Pairs))
RNA_celltype_TF_network_list <- split.data.frame(RNA_celltype_TF_network, f = factor(RNA_celltype_TF_network$Celltype))

RNA_celltype_TF <- as.data.frame.array(table(RNA_celltype_TF_network$TF, RNA_celltype_TF_network$Celltype))

##### 超几何检验
TF_target_num <- as.data.frame(table(Cor_result_select$TF))
colnames(TF_target_num) <- c("TF", "Target_num")
rownames(TF_target_num) <- TF_target_num$TF
RNA_celltype_TF_filtered <- RNA_celltype_TF
for (i in colnames(RNA_celltype_TF_filtered)) {
  for (j in rownames(RNA_celltype_TF_filtered)) {
    x <- RNA_celltype_TF_filtered[j,i]
    tf_num <- TF_target_num[j,2]
    k <- nrow(RNA_celltype_markers_list[[i]])
    p <- phyper(x-1, tf_num, 30000, k, lower.tail = F)
    RNA_celltype_TF_filtered[j,i] <- p
  }
}

RNA_celltype_TF_filtered_temp <- data.frame(TF = rownames(RNA_celltype_TF_filtered),
                                            RNA_celltype_TF_filtered,
                                            check.names = F)
RNA_celltype_TF_filtered_temp <- reshape::melt.data.frame(RNA_celltype_TF_filtered_temp,
                                                          id.vars = "TF")
colnames(RNA_celltype_TF_filtered_temp) <- c("TF", "Celltype", "P")

RNA_celltype_TF_temp <- data.frame(TF = rownames(RNA_celltype_TF),
                                   RNA_celltype_TF,
                                   check.names = F)
RNA_celltype_TF_temp <- reshape::melt.data.frame(RNA_celltype_TF_temp,
                                                 id.vars = "TF")
colnames(RNA_celltype_TF_temp) <- c("TF", "Celltype", "Target_num")

RNA_celltype_regulon <- data.frame(TF = RNA_celltype_TF_temp$TF,
                                   Celltype = RNA_celltype_TF_temp$Celltype,
                                   Target_num = RNA_celltype_TF_temp$Target_num,
                                   P = RNA_celltype_TF_filtered_temp$P)
RNA_celltype_regulon <- RNA_celltype_regulon[which(RNA_celltype_regulon$P < 0.05),]
table(RNA_celltype_regulon$Celltype)

RNA_celltype_regulon$Regulon <- paste0(RNA_celltype_regulon$TF,"(",
                                       RNA_celltype_regulon$Target_num,"g)")
sum(duplicated(RNA_celltype_regulon$Regulon))
Target_list <- c()
for (i in 1:nrow(RNA_celltype_regulon)) {
  temp <- RNA_celltype_TF_network_list[[RNA_celltype_regulon[i,"Celltype"]]]
  targets <- temp[which(temp$TF == RNA_celltype_regulon[i,"TF"]),"Target"]
  targets <- sort(as.character(targets))
  targets <- paste(targets, collapse = " | ")
  Target_list <- c(Target_list, targets)
}
sum(duplicated(Target_list))
RNA_celltype_regulon$Target_list <- Target_list
RNA_celltype_regulon$TF_Celltype <- paste0(RNA_celltype_regulon$TF, "_", RNA_celltype_regulon$Celltype)
rownames(RNA_celltype_regulon) <- RNA_celltype_regulon$TF_Celltype
RNA_celltype_regulon$Regulon_list <- paste0(RNA_celltype_regulon$Regulon, "_",
                                            RNA_celltype_regulon$Target_list)
for (i in unique(RNA_celltype_regulon$Regulon)) {
  temp <- RNA_celltype_regulon[which(RNA_celltype_regulon$Regulon == i),]
  if (nrow(temp) > 1) {
    if (length(unique(temp$Regulon_list)) == nrow(temp)) {
      temp$Regulon <- paste0(temp$Regulon, "_", 1:nrow(temp))
    }
    RNA_celltype_regulon[rownames(temp),"Regulon"] <- temp$Regulon
  }
}

for (i in names(which(table(RNA_celltype_regulon$Regulon) > 1))) {
  temp <- RNA_celltype_regulon[which(RNA_celltype_regulon$Regulon == i),]
  if (length(unique(temp$Target_list)) > 1) {
    temp2 <- as.data.frame(table(temp$Target_list))
    temp2$Regulon <- paste0(i, "_", temp2$Freq)
    rownames(temp2) <- temp2$Var1
    temp$Regulon <- temp2[temp$Target_list,"Regulon"]
    RNA_celltype_regulon[rownames(temp),"Regulon"] <- temp$Regulon
  }
}

table(RNA_celltype_regulon$Regulon)

RNA_celltype_regulon_pairs <- c()
for (i in 1:nrow(RNA_celltype_regulon)) {
  tf <- RNA_celltype_regulon[i, "TF"]
  celltype <- as.character(RNA_celltype_regulon[i, "Celltype"])
  regulon <- as.character(RNA_celltype_regulon[i, "Regulon"])
  targets <- unlist(strsplit(RNA_celltype_regulon[i, "Target_list"], " | ", fixed = T))
  temp <- data.frame(TF = tf,
                     Celltype = celltype,
                     Regulon = regulon,
                     Targets = targets)
  RNA_celltype_regulon_pairs <- as.data.frame(rbind(RNA_celltype_regulon_pairs,
                                                    temp))
}
length(unique(RNA_celltype_regulon_pairs$Targets))
Celltype_regulon <- as.data.frame.array(table(RNA_celltype_regulon_pairs$Regulon, RNA_celltype_regulon_pairs$Celltype))
apply(Celltype_regulon, 1, max)
Regulon_pairs <- RNA_celltype_regulon_pairs[,c("Regulon", "Targets")]
Regulon_pairs <- Regulon_pairs[!duplicated(Regulon_pairs),]
Regulon_pairs_list <- split.data.frame(Regulon_pairs, f = list(Regulon_pairs$Regulon))
Regulon_pairs_list <- lapply(Regulon_pairs_list, function(x) {
  x[,"Targets"]
})

# regulon中加入TF
for (i in names(Regulon_pairs_list)) {
  tf <- unlist(strsplit(i, "(", fixed = T))[1]
  Regulon_pairs_list[[i]] <- c(Regulon_pairs_list[[i]], tf)
}

# 有重复的regulon

openxlsx::write.xlsx(RNA_celltype_regulon_pairs,
                     "Supplementary_Table_2(RNA_celltype_regulon_pairs).xlsx")

### 对所有细胞进行 AUCell 打分
exprMatrix <- RNA@assays$RNA@data
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 30, plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(Regulon_pairs_list, cells_rankings)
cells_AUC_matrix <- cells_AUC@assays@data@listData[["AUC"]]
cells_AUC_matrix <- t(cells_AUC_matrix)
cells_AUC_celltype <- aggregate.data.frame(as.data.frame(cells_AUC_matrix), by = list(RNA$Celltype), FUN = mean)
cells_AUC_celltype <- cells_AUC_celltype[-which(cells_AUC_celltype$Group.1 == "Unidentified"),]
rownames(cells_AUC_celltype) <- cells_AUC_celltype$Group.1
cells_AUC_celltype <- cells_AUC_celltype[,-1]
cells_AUC_celltype <- cells_AUC_celltype[,]
cells_AUC_celltype <- cells_AUC_celltype[,RNA_celltype_regulon[!duplicated(RNA_celltype_regulon$Regulon),"Regulon"]]
pheatmap(cells_AUC_celltype, scale = "column", cluster_rows = T, cluster_cols = T,
         color = colorRampPalette(c("blue", "yellow", "red"))(100), cutree_cols = 4,
         border_color = "black")

ncol(cells_AUC_celltype)

# 导入motif
{
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
}

gene_id_map <- readRDS("gene_id_map.rds")
rownames(gene_id_map) <- gene_id_map$gene_id
TF_RAP <- unique(RNA_celltype_regulon$TF)
TF_Symbol <- gene_id_map[TF_RAP,2]
TF_id_map <- data.frame(TF_RAP = TF_RAP,
                        TF_Symbol = TF_Symbol,
                        row.names = TF_RAP)
rownames(Annotation_genes) <- Annotation_genes$DataSets_Symbol
TF_id_map$TF_Symbol %in% Annotation_genes$DataSets_Symbol
TF_id_map$names <- Annotation_genes[TF_id_map$TF_Symbol, "CGSNL Gene Symbol"]
TF_id_map$names[which(is.na(TF_id_map$names))] <- TF_id_map$TF_Symbol[which(is.na(TF_id_map$names))]
rownames(TF_id_map) <- TF_id_map$TF_RAP

RNA_celltype_regulon$TF_Symbol <- TF_id_map[RNA_celltype_regulon$TF, "names"]
# Regulon_pairs$TF <- unlist(lapply(Regulon_pairs$Regulon, function(x){
#   unlist(strsplit(x, "(", fixed = T))[1]
# }))
# Regulon_pairs$TF_Symbol <- TF_id_map[Regulon_pairs$TF, "TF_Symbol"]

RNA_celltype_regulon_pairs
RNA_celltype_regulon_pairs$TF_Symbol <- TF_id_map[RNA_celltype_regulon_pairs$TF, "names"]
Regulon_Symbol <- c()
for (i in 1:nrow(RNA_celltype_regulon_pairs)) {
  temp <- unlist(strsplit(RNA_celltype_regulon_pairs[i,"Regulon"],
                          "(", fixed = T))[2]
  Regulon_Symbol <- c(Regulon_Symbol,
                      paste0(RNA_celltype_regulon_pairs[i,"TF_Symbol"],"(",temp))
}
RNA_celltype_regulon_pairs$Regulon_Symbol <- Regulon_Symbol
Regulon_ID_Symbol <- RNA_celltype_regulon_pairs[,c("Regulon", "Regulon_Symbol")]
Regulon_ID_Symbol <- Regulon_ID_Symbol[!duplicated(Regulon_ID_Symbol),]
rownames(Regulon_ID_Symbol) <- Regulon_ID_Symbol$Regulon

colnames(cells_AUC_celltype) <- Regulon_ID_Symbol[colnames(cells_AUC_celltype), "Regulon_Symbol"]


RNA_celltype_color <- c("#faa818", # BM
                        "#62adac", # Columella initials
                        "#d5a0e3", # Cortex
                        "#fbdf72", # cb/b
                        "#f29f34", # Dividing inner
                        "#ff8250", # Dividing outer
                        "#367d7d", # Endodermis
                        "#69a2bd", # Endodermis-cortical initial
                        "#FF88C2", "#FF8888", # Epidermal
                        "#8ff0a4", # Fiber
                        "#5555FF", # Inflorescence meristem (IM)
                        "#ccc223", # Intermediate anther
                        "#B94FFF", # Large parenchyma (MO)
                        "#bd6259", # Late anther
                        "#6ec9bb", # Lemma (le)
                        "#5eb5c7", # Mesophyll (MO)
                        "#4c9fd4", # Mesophyll initial
                        "#3d8ede", # Mesophyll precursor
                        "#df9cbf", # Other Spikelet
                        "#5c7ada", # Phloem
                        "#62589d", # Procambium
                        "#74fad0", # Proliferating cell
                        "#A4A5EE", # Rachis
                        "#de87a6", # Shoot endodermis
                        "#c85c6c", # Shoot meristem initials
                        "#b4446c", # Shoot meristematic cell
                        "#A27BFF", # Spikelet meristem (SM)
                        "#FF8055", # Suspensor
                        "#B7C792", # Upper protoderm
                        "#B7D62D" # Vascular initial
)
names(RNA_celltype_color) <- rownames(cells_AUC_celltype)
row_anno <- rowAnnotation(df = data.frame(Celltype = rownames(cells_AUC_celltype)),
                          col = list(Celltype = RNA_celltype_color),
                          show_legend = F, show_annotation_name = F)
annotation_row <- data.frame(
  Celltype = rownames(cells_AUC_celltype)
)
rownames(annotation_row) <- rownames(cells_AUC_celltype)
ann_colors <- list(
  Celltype = RNA_celltype_color
  )

pdf("Figure3_改_celltype_regulon_AUCell_score_heatmap.pdf", width = 14, height = 5)
pt <- pheatmap(cells_AUC_celltype, scale = "column", cluster_rows = T, cluster_cols = T,
         # color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")))(100),
         color =  colorRampPalette(c("blue", "white", "red"))(100),
         cutree_cols = 4, show_colnames = F, annotation_colors = ann_colors, annotation_row = annotation_row,
         border_color = NULL)
dev.off()

tree_col <- pt$tree_col
cluster_tree_cut <- cutree(tree_col, k = 4)
table(cluster_tree_cut)
cut_tree_module <- data.frame(treecut = c(1, 2, 3, 4),
                              module = c("Module 3", "Module 2", "Module 4", "Module 1"))
rownames(cut_tree_module) <- cut_tree_module$treecut
cluster_tree_cut <- as.data.frame(cluster_tree_cut)
cluster_tree_cut$Module <- cut_tree_module[cluster_tree_cut$cluster_tree_cut, "module"]
cluster_tree_cut$regulon <- rownames(cluster_tree_cut)

genes_list <- c()
for (i in c("Module 1", "Module 2", "Module 3", "Module 4")) {
  temp <- cluster_tree_cut[which(cluster_tree_cut$Module == i),]
  temp_regulon_pairs <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Regulon_Symbol %in% temp$regulon),]
  temp_regulon_pairs <- temp_regulon_pairs[,c("Regulon", "Targets")]
  temp_regulon_pairs <- temp_regulon_pairs[!duplicated(temp_regulon_pairs),]
  temp_regulon_pairs_table <- as.data.frame.array(table(temp_regulon_pairs$Regulon, temp_regulon_pairs$Targets))
  temp_regulon_pairs_table <- colSums(temp_regulon_pairs_table)
  genes <- names(temp_regulon_pairs_table)[which(temp_regulon_pairs_table > 1)]
  
  temp <- cluster_tree_cut[which(cluster_tree_cut$Module == i),]
  temp_regulon_pairs <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Regulon_Symbol %in% temp$regulon),]
  tf_number <- length(unique(temp_regulon_pairs$TF))
  target_number <- length(unique(temp_regulon_pairs$Targets))
  regulon_number <- length(unique(temp_regulon_pairs$Regulon_Symbol))
  genes <- c(temp_regulon_pairs$TF, genes)
  genes <- unique(genes)
  
  temp <- list(module = i,
               tf_number = tf_number,
               target_number = target_number,
               regulon_number = regulon_number,
               gene_list = genes
               )
  genes_list <- c(genes_list, list(temp))
}
names(genes_list) <- c("Module 1", "Module 2", "Module 3", "Module 4")

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

### Gene Ontology
GO_result_module <- c()
for (i in names(genes_list)) {
  print(i)
  genes <- genes_list[[i]]$gene_list
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
  GO_result$Module <- i
  GO_result_module <- as.data.frame(rbind(GO_result_module, GO_result))
  
}

GO_result_module <- GO_result_module[which(GO_result_module$P_value < 0.05),]
GO_result_module_BP <- GO_result_module[which(GO_result_module$Ontology == "BP"),]
GO_result_module_BP <- GO_result_module_BP[order(GO_result_module_BP$P_value, decreasing = F),]

GO_result_module_BP_M1 <- GO_result_module_BP[which(GO_result_module_BP$Module == "Module 1"),]
GO_result_module_BP_M2 <- GO_result_module_BP[which(GO_result_module_BP$Module == "Module 2"),]
GO_result_module_BP_M3 <- GO_result_module_BP[which(GO_result_module_BP$Module == "Module 3"),]
GO_result_module_BP_M4 <- GO_result_module_BP[which(GO_result_module_BP$Module == "Module 4"),]

GO_result_module_BP_M1_select <- GO_result_module_BP_M1[c(2:4,7,16,20,21,24,35,43), ] # ,20,21,24
GO_result_module_BP_M2_select <- GO_result_module_BP_M2[c(1:9,13), ]
M3 <- c("cell cycle",
        "mitotic cell cycle", "cell cycle process", "nuclear division", "cell differentiation",
        "cellular developmental process", "multicellular organism development", "system development",
        "anatomical structure development", "positive regulation of transcription, DNA-templated")
GO_result_module_BP_M3_select <- GO_result_module_BP_M3[GO_result_module_BP_M3$Description %in% M3, ]
GO_result_module_BP_M4_select <- GO_result_module_BP_M4[c(1,3,4,7,11,13,14,15,18,19), ]

GO_result_module_BP_M1_select$log10P <- -log10(GO_result_module_BP_M1_select$P_value)
GO_result_module_BP_M1_select <- GO_result_module_BP_M1_select[order(GO_result_module_BP_M1_select$log10P,
                                                                     decreasing = T),]
GO_result_module_BP_M2_select$log10P <- -log10(GO_result_module_BP_M2_select$P_value)
GO_result_module_BP_M2_select <- GO_result_module_BP_M2_select[order(GO_result_module_BP_M2_select$log10P,
                                                                     decreasing = T),]
GO_result_module_BP_M3_select$log10P <- -log10(GO_result_module_BP_M3_select$P_value)
GO_result_module_BP_M3_select <- GO_result_module_BP_M3_select[order(GO_result_module_BP_M3_select$log10P,
                                                                     decreasing = T),]
GO_result_module_BP_M4_select$log10P <- -log10(GO_result_module_BP_M4_select$P_value)
GO_result_module_BP_M4_select <- GO_result_module_BP_M4_select[order(GO_result_module_BP_M4_select$log10P,
                                                                     decreasing = T),]
GO_result_module <- as.data.frame(rbind(GO_result_module_BP_M1_select,
                                        GO_result_module_BP_M2_select,
                                        GO_result_module_BP_M3_select,
                                        GO_result_module_BP_M4_select))
GO_result_module$Description <- paste0(GO_result_module$Description, "_", GO_result_module$Module)
GO_result_module$Description <- factor(GO_result_module$Description,
                                       levels = GO_result_module$Description)
pdf("Figure3_GO_result_module.pdf", width = 14, height = 7)
ggplot(GO_result_module, aes(x = Description, y = log10P, fill = Module)) +
  geom_bar(stat = "identity", width = 0.75) +
  facet_wrap(.~Module, nrow = 1, scales = "free") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

regulon_target <- RNA_celltype_regulon_pairs[,c("TF", "Regulon", "Targets", "TF_Symbol", "Regulon_Symbol")]
regulon_target <- regulon_target[!duplicated(regulon_target), ]
regulon_target$module <- cluster_tree_cut[regulon_target$Regulon_Symbol, "Module"]

RNA_celltype_regulon_pairs$Module <- cluster_tree_cut[RNA_celltype_regulon_pairs$Regulon_Symbol, "Module"]

# module中的靶基因被多少TF调节
{
  M1_regulon_target <- regulon_target[which(regulon_target$module == "Module 1"),]
  M1_TF_target_table <- as.data.frame.array(table(M1_regulon_target$Targets, M1_regulon_target$TF_Symbol))
  M1_TF_target_table <- as.matrix(M1_TF_target_table)
  M1_TF_target_table[which(M1_TF_target_table > 0)] <- 1
  table(rowSums(M1_TF_target_table))
  
  M2_regulon_target <- regulon_target[which(regulon_target$module == "Module 2"),]
  M2_TF_target_table <- as.data.frame.array(table(M2_regulon_target$Targets, M2_regulon_target$TF_Symbol))
  M2_TF_target_table <- as.matrix(M2_TF_target_table)
  M2_TF_target_table[which(M2_TF_target_table > 0)] <- 1
  table(rowSums(M2_TF_target_table))
  
  M3_regulon_target <- regulon_target[which(regulon_target$module == "Module 3"),]
  M3_TF_target_table <- as.data.frame.array(table(M3_regulon_target$Targets, M3_regulon_target$TF_Symbol))
  M3_TF_target_table <- as.matrix(M3_TF_target_table)
  M3_TF_target_table[which(M3_TF_target_table > 0)] <- 1
  table(rowSums(M3_TF_target_table))
  
  M4_regulon_target <- regulon_target[which(regulon_target$module == "Module 4"),]
  M4_TF_target_table <- as.data.frame.array(table(M4_regulon_target$Targets, M4_regulon_target$TF_Symbol))
  M4_TF_target_table <- as.matrix(M4_TF_target_table)
  M4_TF_target_table[which(M4_TF_target_table > 0)] <- 1
  table(rowSums(M4_TF_target_table))
  
  
}

# 在module中，同一个TF在不同细胞类型中调控不同的靶基因
{
  # M2举例
  temp <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Module == "Module 2"),]
  OsMADS2_M2 <- temp[which(temp$TF_Symbol == "OsMADS2"),]
  genes_list <- c()
  for (i in unique(OsMADS2_M2$Celltype)) {
    genes_list <- c(genes_list, list(OsMADS2_M2$Targets[which(OsMADS2_M2$Celltype == i)]))
  }
  names(genes_list) <- unique(OsMADS2_M2$Celltype)
  GO_result_OsMADS2_M2 <- c()
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
    GO_result_OsMADS2_M2 <- as.data.frame(rbind(GO_result_OsMADS2_M2, GO_result))
    
  }
  GO_result_OsMADS2_M2 <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$P_value < 0.05),]
  unique(GO_result_OsMADS2_M2$Celltype)
  GO_result_OsMADS2_M2 <- GO_result_OsMADS2_M2[order(GO_result_OsMADS2_M2$P_value, decreasing = F),]
  GO_result_OsMADS2_M2 <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Ontology == "BP"),]
  GO_result_OsMADS2_M2_CI <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Columella initials"),]
  GO_result_OsMADS2_M2_DO <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Dividing outer"),]
  GO_result_OsMADS2_M2_Epi <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Epidermis"),]
  GO_result_OsMADS2_M2_IA <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Intermediate anther"),]
  GO_result_OsMADS2_M2_LA <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Late anther"),]
  GO_result_OsMADS2_M2_MI <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Mesophyll initial"),]
  GO_result_OsMADS2_M2_MP <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Mesophyll precursor"),]
  GO_result_OsMADS2_M2_Sus <- GO_result_OsMADS2_M2[which(GO_result_OsMADS2_M2$Celltype == "Suspensor"),]
  
  GO_result_OsMADS2_M2_temp <- as.data.frame.array(table(GO_result_OsMADS2_M2$Description,
                                                         GO_result_OsMADS2_M2$Celltype))
  GO_result_OsMADS2_M2_temp <- rowSums(GO_result_OsMADS2_M2_temp)
  sort(GO_result_OsMADS2_M2_temp, decreasing = F)
  GO_select <- names(GO_result_OsMADS2_M2_temp)[which(GO_result_OsMADS2_M2_temp == 1)]
  
  
  GO_result_OsMADS2_M2_temp <- as.data.frame.array(table(GO_result_OsMADS2_M2$Description,
                                                         GO_result_OsMADS2_M2$Celltype))
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[GO_select,]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,1]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,2]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,3]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,4]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,5]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,6]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,7]),]
  GO_result_OsMADS2_M2_temp <- GO_result_OsMADS2_M2_temp[order(GO_result_OsMADS2_M2_temp[,8]),]
  rownames(GO_result_OsMADS2_M2) <- paste0(GO_result_OsMADS2_M2$Description, "_", GO_result_OsMADS2_M2$Celltype)
  for (i in 1:nrow(GO_result_OsMADS2_M2_temp)) {
    for (j in 1:ncol(GO_result_OsMADS2_M2_temp)) {
      if (GO_result_OsMADS2_M2_temp[i,j] == 1) {
        GO_result_OsMADS2_M2_temp[i,j] <- GO_result_OsMADS2_M2[paste0(rownames(GO_result_OsMADS2_M2_temp)[i], "_",
                                                                      colnames(GO_result_OsMADS2_M2_temp)[j]), "P_value"]
      }
    }
  }
  
  GO_result_OsMADS2_M2$`-log10(P)` <- -log10(GO_result_OsMADS2_M2$P_value)
  # GO_select <- c( #"anatomical structure formation involved in morphogenesis",
  # #                "macromolecule localization",
  # #                "vegetative to reproductive phase transition of meristem",
  # #                "organelle organization",
  # #                "cell wall organization",
  # #                "nucleobase-containing compound metabolic process"
  #                # "macromolecular complex assembly"
  #                # "seed maturation", "developmental process",
  #                # "reproductive structure development"
  #   "cellular developmental process",
  #   "establishment or maintenance of cell polarity",
  #   "post-embryonic plant organ morphogenesis",
  #   "specification of plant organ identity"
  #                )
  # GO_select <- names(GO_result_OsMADS2_M2_temp)[which(GO_result_OsMADS2_M2_temp == 1)]
  # GO_matrix <- matrix(NA, nrow = 8, ncol = length(GO_select))
  # rownames(GO_matrix) <- c("Mesophyll initial", "Mesophyll precursor", "Suspensor", "Dividing outer",
  #                          "Columella initials", "Intermediate anther", "Late anther", "Epidermis")
  # colnames(GO_matrix) <- GO_select
  # 
  # for (i in 1:nrow(GO_matrix)) {
  #   for (j in 1:ncol(GO_matrix)) {
  #     GO_matrix[i,j] <- GO_result_OsMADS2_M2[paste0(colnames(GO_matrix)[j], "_",
  #                                                   rownames(GO_matrix)[i]), "-log10(P)"]
  #   }
  # }
  GO_result_OsMADS2_M2_temp <- as.matrix(GO_result_OsMADS2_M2_temp)
  GO_result_OsMADS2_M2_temp[which(GO_result_OsMADS2_M2_temp == 0)] <- NA
  GO_result_OsMADS2_M2_temp <- -log10(GO_result_OsMADS2_M2_temp)
  pdf("Figure3_GO_matrix_1个细胞类型.pdf", width = 3.5, height = 5.5)
  pheatmap(GO_result_OsMADS2_M2_temp, cluster_rows = F, cluster_cols = F, fontsize = 20, border_color = NA,
           color =  colorRampPalette(c("#FFD700", "red"))(100), show_rownames = F)
  dev.off()

  
  # pdf("Figure3_GO_matrix_8个细胞类型.pdf", width = 5.5, height = 7.5)
  # pheatmap(GO_matrix, cluster_rows = F, cluster_cols = F, fontsize = 15, border_color = NA,
  #          color =  colorRampPalette(c("#FFD700", "red"))(100), show_rownames = T)
  # dev.off()
  
 
}


# 靶基因被多少TF 调控
{
  temp <- as.data.frame.array(table(RNA_celltype_regulon_pairs$TF, RNA_celltype_regulon_pairs$Targets))
  temp <- as.matrix(temp)
  temp[which(temp > 0)] <- 1
  table(colSums(temp))
}

# ATAC 验证
# 导入TF motif
{
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
}
{
  proj_pass_filter <- readRDS("3.proj_pass_filter_withPeaks.rds")
  pathToMacs2 <- "/home/heshidian/mambaforge/envs/common/bin/macs2"
  cells_select <- as.data.frame(getCellColData(proj_pass_filter))
  cells_select <- cells_select[which(cells_select$TissuesSub %in% c("IM0.5cm",
                                                                    "IM1cm",
                                                                    "10dLeaf",
                                                                    "10dRoot",
                                                                    "60dRootTip")),]
  # cells_select <- cells_select[which(cells_select$Celltype %in% common_cells_type),]
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
  common_cells_type <- intersect(unique(proj_pass_filter_select$Celltype),
                                 unique(RNA_select$Celltype))
  
  RNA_select <- subset(RNA_select, subset = Celltype %in% common_cells_type)
  proj_pass_filter_select <- proj_pass_filter_select[which(proj_pass_filter_select$Celltype %in% common_cells_type),]
  
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
  RNA_select <- NormalizeData(object = RNA_select, verbose = FALSE)
  common_cells_type
  common_cells_type_2 <- common_cells_type[-2]
  
  proj_pass_filter_select <- proj_pass_filter_select[which(proj_pass_filter_select$Celltype %in% common_cells_type_2),]
  RNA_select <- subset(RNA_select, subset = Celltype %in% common_cells_type_2)
  
  groupList <- SimpleList(
    Cortex = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Cortex"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Cortex")]
    ),
    # Epidermis = SimpleList( RNA中细胞数目太少
    #   ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Epidermis"],
    #   RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Epidermis")]
    # ),
    Endodermis = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Endodermis"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Endodermis")]
    ),
    Phloem = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Phloem"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Phloem")]
    ),
    Procambium = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Procambium"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Procambium")]
    ),
    Fiber = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Fiber"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Fiber")]
    ),
    "Mesophyll (MO)" = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Mesophyll (MO)"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Mesophyll (MO)")]
    ),
    "Mesophyll precursor" = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Mesophyll precursor"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Mesophyll precursor")]
    ),
    Rachis = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Rachis"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Rachis")]
    ),
    "Cryptic bract/bract (cb/b)" = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Cryptic bract/bract (cb/b)"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Cryptic bract/bract (cb/b)")]
    ),
    "Inflorescence meristem (IM)" = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Inflorescence meristem (IM)"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Inflorescence meristem (IM)")]
    ),
    "Branch meristems (BM)" = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Branch meristems (BM)"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Branch meristems (BM)")]
    ),
    "Spikelet meristem (SM)" = SimpleList(
      ATAC = proj_pass_filter_select$cellNames[proj_pass_filter_select$Celltype == "Spikelet meristem (SM)"],
      RNA = colnames(RNA_select)[which(RNA_select$Celltype == "Spikelet meristem (SM)")]
    )
  )
  
  addArchRThreads(threads = 1)
  proj_pass_filter_select <- addGeneIntegrationMatrix(
    ArchRProj = proj_pass_filter_select, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = RNA_select,
    addToArrow = FALSE,
    groupRNA = "Celltype",
    groupATAC = "Celltype",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co",
    sampleCellsATAC = 5000,
    sampleCellsRNA = 5000,
    groupList = groupList
  )
  
  GeneIntegrationMatrix <- getMatrixFromProject(proj_pass_filter_select,
                                                useMatrix = "GeneIntegrationMatrix")
  GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix"]]@Dimnames[[1]] <- GeneIntegrationMatrix@elementMetadata@listData[["name"]]
  
  GeneIntegrationMatrix <- CreateSeuratObject(GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix"]])
  library(Rmagic)
  RNA_imputated <-magic(GeneIntegrationMatrix)
  RNA_imputated <- RNA_imputated@assays[["MAGIC_RNA"]]
  
  
  gene_id_map <- readRDS("gene_id_map.rds")
  rownames(gene_id_map) <- gene_id_map$symbol
  GeneIntegrationMatrix_imputated <- RNA_imputated
  rm(RNA_imputated)
  gc()
  GeneIntegrationMatrix_imputated <- t(GeneIntegrationMatrix_imputated@data)
  sum(colnames(GeneIntegrationMatrix_imputated) %in% gene_id_map$symbol)
  dim(GeneIntegrationMatrix_imputated)
  colnames(GeneIntegrationMatrix_imputated) <- gene_id_map[colnames(GeneIntegrationMatrix_imputated), "gene_id"]
  
  sum(Regulon_pairs$TF %in% colnames(GeneIntegrationMatrix_imputated)) == nrow(Regulon_pairs)
  sum(Regulon_pairs$Targets %in% colnames(GeneIntegrationMatrix_imputated))
  sum(Regulon_pairs$Targets %in% colnames(GeneIntegrationMatrix_imputated)) == nrow(Regulon_pairs)
  Regulon_pairs$Pairs <- paste0(Regulon_pairs$TF,"_",Regulon_pairs$Targets)
  Regulon_pairs$cor <- Cor_result_select[Regulon_pairs$Pairs, "Pearson"]
  
  ATAC_cor <- c()
  for (i in 1:nrow(Regulon_pairs)) {
    tf <- Regulon_pairs$TF[i]
    target <- Regulon_pairs$Targets[i]
    if (target %in% colnames(GeneIntegrationMatrix_imputated)) {
      temp_cor <- cor(GeneIntegrationMatrix_imputated[,tf],
                      GeneIntegrationMatrix_imputated[,target])
      ATAC_cor <- c(ATAC_cor, temp_cor)
    } else {
      ATAC_cor <- c(ATAC_cor, NA)
    }
  }
  
  temp <- data.frame(Regulon_pairs,
                     ATAC_cor = ATAC_cor)
  temp <- na.omit(temp)
  
  ggplot(temp, aes(x = cor, y = ATAC_cor)) +
    geom_point() +
    geom_smooth() +
    ggpubr::stat_cor()
  
  
  
  proj_pass_filter_select$Celltype_2 <- proj_pass_filter_select$Celltype
  proj_pass_filter_select$Celltype_2[which(proj_pass_filter_select$Celltype_2 %in% c("Mesophyll precursor",
                                                                                     "Epidermis"))]
  
  Annotation_genes$`CGSNL Gene Symbol` == "OsMADS2"
  proj_pass_filter_select <- addGroupCoverages(proj_pass_filter_select, groupBy = "Celltype", force = TRUE)
  proj_pass_filter_select <- addBgdPeaks(proj_pass_filter_select)
  
  proj_pass_filter_select <- addMotifAnnotations(ArchRProj = proj_pass_filter_select, motifPWMs = motif_all,
                                                       annoName = "TF-Motif", force = T)
  
  set.seed(123)
  
  proj_pass_filter_select <- addMotifAnnotations(ArchRProj = proj_pass_filter_select, motifPWMs = motif_all,
                                                 annoName = "TF-Motif", force = T)
  motifPositions <- getPositions(proj_pass_filter_select, name = "TF-Motif")
  proj_pass_filter_select <- addDeviationsMatrix(
    ArchRProj = proj_pass_filter_select, 
    peakAnnotation = "TF-Motif",
    matrixName = "MotifMatrix",
    force = TRUE
  )
  plotVarDev <- getVarDeviations(proj_pass_filter_select, name = "MotifMatrix", plot = FALSE)
  motifdeviation_matrix <- getMatrixFromProject(proj_pass_filter_select, useMatrix = "MotifMatrix" )
  motifdeviation_matrix_dscore <- motifdeviation_matrix@assays@data@listData[["deviations"]]
  dim(motifdeviation_matrix_dscore)
  motifdeviation_matrix_dscore <- t(motifdeviation_matrix_dscore)
  motifdeviation_matrix_dscore <- aggregate.data.frame(motifdeviation_matrix_dscore, by = list(motifdeviation_matrix@colData[rownames(motifdeviation_matrix_dscore),"Celltype"]), FUN = median)
  rownames(motifdeviation_matrix_dscore) <- motifdeviation_matrix_dscore$Group.1
  motifdeviation_matrix_dscore <- motifdeviation_matrix_dscore[,-1]
  motifdeviation_matrix_dscore
  
  module_TF <- table(RNA_celltype_regulon_pairs$TF_Symbol, RNA_celltype_regulon_pairs$Module)
  module_TF <- as.data.frame.array(module_TF)
  module_TF <- module_TF[rownames(module_TF) %in% colnames(motifdeviation_matrix_dscore),]
  module_TF <- module_TF[-which(rownames(module_TF) == "OsCXC2"),]
  module1_TF <- rownames(module_TF)[which(module_TF$`Module 1` > 0)]
  module2_TF <- rownames(module_TF)[which(module_TF$`Module 2` > 0)]
  module3_TF <- rownames(module_TF)[which(module_TF$`Module 3` > 0)]
  module4_TF <- rownames(module_TF)[which(module_TF$`Module 4` > 0)]
  
  colanno <- data.frame(row.names = c(module1_TF,module2_TF,module3_TF,module4_TF),
                        Module = c(rep("module1", length(module1_TF)), rep("module2", length(module2_TF)),
                                   rep("module3", length(module3_TF)), rep("module4", length(module4_TF))))
  
  pheatmap(motifdeviation_matrix_dscore[,c(module1_TF,module2_TF,module3_TF,module4_TF)],
           annotation_col = colanno, cluster_cols = F)
  
  
  ### marker peak enrichment
  markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_pass_filter_select, 
    useMatrix = "PeakMatrix", 
    groupBy = "Celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_pass_filter_select,
    peakAnnotation = "TF-Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1"
  )
  
  mlog10Padj <- as.data.frame(enrichMotifs@assays@data@listData[["mlog10Padj"]])
  
  module_TF <- table(RNA_celltype_regulon_pairs$TF_Symbol, RNA_celltype_regulon_pairs$Module)
  module_TF <- as.data.frame.array(module_TF)
  module_TF <- module_TF[rownames(module_TF) %in% rownames(mlog10Padj),]
  module_TF <- module_TF[-which(rownames(module_TF) == "OsCXC2"),]
  module1_TF <- rownames(module_TF)[which(module_TF$`Module 1` > 0)]
  module2_TF <- rownames(module_TF)[which(module_TF$`Module 2` > 0)]
  module3_TF <- rownames(module_TF)[which(module_TF$`Module 3` > 0)]
  module4_TF <- rownames(module_TF)[which(module_TF$`Module 4` > 0)]
  pheatmap(mlog10Padj[c(module1_TF,
                        module2_TF,
                        module3_TF,
                        module4_TF),], cluster_rows = F, cluster_cols = F)

  
  ### 整合 RNA ATAC
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
  
  GeneIntegrationMatrix <- getMatrixFromProject(proj_pass_filter_select,
                                                useMatrix = "GeneIntegrationMatrix")
  GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix"]]@Dimnames[[1]] <- GeneIntegrationMatrix@elementMetadata@listData[["name"]]
  
  GeneIntegrationMatrix <- CreateSeuratObject(GeneIntegrationMatrix@assays@data@listData[["GeneIntegrationMatrix"]])
  saveRDS(GeneIntegrationMatrix, "Figure3_GeneIntegrationMatrix.rds")
  
  gene_id_map <- readRDS("gene_id_map.rds")
  rownames(gene_id_map) <- gene_id_map$symbol
  GeneIntegrationMatrix_imputated <- readRDS("Figure3_GeneIntegrationMatrix_imputated.rds")
  GeneIntegrationMatrix_imputated <- t(GeneIntegrationMatrix_imputated@data)
  sum(colnames(GeneIntegrationMatrix_imputated) %in% gene_id_map$symbol)
  dim(GeneIntegrationMatrix_imputated)
  colnames(GeneIntegrationMatrix_imputated) <- gene_id_map[colnames(GeneIntegrationMatrix_imputated), "gene_id"]
  
  sum(Regulon_pairs$TF %in% colnames(GeneIntegrationMatrix_imputated)) == nrow(Regulon_pairs)
  sum(Regulon_pairs$Targets %in% colnames(GeneIntegrationMatrix_imputated)) == nrow(Regulon_pairs)
  Regulon_pairs$Pairs <- paste0(Regulon_pairs$TF,"_",Regulon_pairs$Targets)
  Regulon_pairs$cor <- Cor_result_select[Regulon_pairs$Pairs, "Pearson"]
  
  ATAC_cor <- c()
  for (i in 1:nrow(Regulon_pairs)) {
    tf <- Regulon_pairs$TF[i]
    target <- Regulon_pairs$Targets[i]
    if (target %in% colnames(GeneIntegrationMatrix_imputated)) {
      temp_cor <- cor(GeneIntegrationMatrix_imputated[,tf],
                      GeneIntegrationMatrix_imputated[,target])
      ATAC_cor <- c(ATAC_cor, temp_cor)
    } else {
      ATAC_cor <- c(ATAC_cor, NA)
    }
  }
  
  temp <- data.frame(Regulon_pairs,
                     ATAC_cor = ATAC_cor)
  temp <- na.omit(temp)
  
  ggplot(temp, aes(x = cor, y = ATAC_cor)) +
    geom_point() +
    geom_smooth() +
    ggpubr::stat_cor()
  
  #### ATAC gene score AUCell
  
  
  
  
  
  
  
  
  
  
  
  seFoot <- getFootprints(
    ArchRProj = proj_pass_filter_common_cells, 
    positions = motifPositions, 
    groupBy = "Celltype"
  )
  
  keep <- names(seFoot@assays@data@listData) %in% RNA_celltype_regulon$TF
  keep <- names(seFoot@assays@data@listData)[keep]
  
  plotFootprints(
    seFoot = seFoot,
    names = keep,
    ArchRProj = proj_pass_filter_common_cells, 
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias_2",
    addDOC = FALSE,
    smoothWindow = 5,
    force = T,
    #pal = RNA_celltype_color[common_cells_type],
    pal = RNA_celltype_color["Spikelet meristem (SM)"],
    baseSize = 10
  )
}



heatmap(as.matrix(cells_AUC_celltype), col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")))(100),
        show_row_names = T, show_column_names = T, row_names_side = c("left"),
        name = paste0("Scaled AUCScore"),
        cluster_rows = T, cluster_columns = F, use_raster = F,
        show_heatmap_legend = T, row_names_gp = gpar(fontsize = 8),
        left_annotation = row_anno, show_row_dend = F,
        column_names_gp = gpar(fontsize = 8),
        top_annotation = column_names,
        scale = "column")



cells_AUC_celltype_scaled <- scale(cells_AUC_celltype)
cells_AUC_celltype_scaled[which(cells_AUC_celltype_scaled > 2.5)] <- 2.5
cells_AUC_celltype_scaled[which(cells_AUC_celltype_scaled < -2.5)] <- -2.5

col_order <- hclust(dist(t(cells_AUC_celltype_scaled)))

column_names <- HeatmapAnnotation(which = "column",
                               ano = anno_mark(at = 1:ncol(cells_AUC_celltype_scaled),
                                               labels = colnames(cells_AUC_celltype_scaled[,col_order[["order"]]]),
                                               labels_gp = gpar(fontsize = 4.5),
                                               link_width = unit(8, "mm"),
                                               extend = unit(0.5, "mm"),
                                               padding = unit(0.3, "mm")))

pdf("Figure3_celltype_regulon.pdf", width = 16, height = 6.5)
Heatmap(cells_AUC_celltype_scaled[,col_order[["order"]]], col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")))(100),
        show_row_names = T, show_column_names = F, row_names_side = c("left"), name = paste0("Scaled AUCScore"),
        cluster_rows = T, cluster_columns = F, use_raster = F, show_heatmap_legend = T, row_names_gp = gpar(fontsize = 8),
        left_annotation = row_anno, show_row_dend = F, column_names_gp = gpar(fontsize = 8), top_annotation = column_names)
dev.off()


RNA_celltype_regulon$Regulon_Symbol <- Regulon_ID_Symbol[RNA_celltype_regulon$Regulon,2]
temp <- RNA_celltype_regulon[,c("Regulon_Symbol", "Celltype")]
temp$Regulon_Symbol <- factor(temp$Regulon_Symbol,
                              levels = rev(colnames(cells_AUC_celltype_scaled[,col_order[["order"]]])))
pdf("Figure3_Regulon_Celltype_num.pdf", width = 15, height = 2.5)
ggplot() +
  geom_bar(data = temp, aes(x = Regulon_Symbol, fill = Celltype)) +
  scale_fill_manual(values = RNA_celltype_color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())
dev.off()

write.table(RNA_celltype_regulon_pairs, "Figure3_RNA_celltype_regulon_pairs.csv",
            sep = ",", row.names = F, col.names = T)
Node_Regulon <- unique(RNA_celltype_regulon_pairs$Regulon_Symbol)
Node_Target <- unique(RNA_celltype_regulon_pairs$Targets)
TF_Network_Node <- data.frame(Node = c(Node_Regulon, Node_Target),
                              Type = c(rep("Regulon",length(Node_Regulon)), rep("Target",length(Node_Target))))
write.table(TF_Network_Node, "Figure3_TF_Network_Node.csv",
            sep = ",", row.names = F, col.names = T)

Pie_infor <- as.data.frame.array(table(RNA_celltype_regulon$Regulon_Symbol,
                                       RNA_celltype_regulon$Celltype))
Pie_infor <- as.data.frame(cbind(Regulon = rownames(Pie_infor),
                                 Pie_infor))
write.table(Pie_infor, "Figure3_TF_Network_Node_Pie_infor.csv",
            sep = ",", row.names = F, col.names = T, quote = F)

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
TF_MSU_RAP <- TF_regulation[,c("TF_MSU", "TF_RAP")]
TF_MSU_RAP <- TF_MSU_RAP[!duplicated(TF_MSU_RAP),]
# unique(TF_MSU_RAP$TF_MSU) %in% name(motif_all)
# rownames(TF_MSU_RAP) <- TF_MSU_RAP$TF_RAP
# RNA_celltype_regulon$TF
MSU_ID <- unlist(lapply(motif_all@listData, function(x){
  x@name
}))
names(motif_all) <- MSU_ID
motif_all <- motif_all[unique(TF_MSU_RAP$TF_MSU)]
MSU_ID <- names(motif_all)
RAP_ID <- match(MSU_ID, TF_MSU_RAP$TF_MSU)
RAP_ID <- TF_MSU_RAP$TF_RAP[RAP_ID]
# 
motif_all <- motif_all[!duplicated(names(motif_all))] # ID duplicated !!!
RAP_ID <- MSU_RAP[match(names(motif_all), MSU_RAP$MSU), "RAP"]
names(motif_all) <- RAP_ID
#motif_all@listData <- motif_all@listData[-98]

proj_pass_filter <- readRDS("3.proj_pass_filter_withPeaks.rds")

pathToMacs2 <- "/home/heshidian/mambaforge/envs/common/bin/macs2"

common_cells_type <- intersect(unique(proj_pass_filter$Celltype),
                               unique(RNA$Celltype))
common_cells_type <- sort(common_cells_type)
proj_pass_filter_common_cells <- proj_pass_filter[which(proj_pass_filter$Celltype %in% common_cells_type),]
set.seed(123)

proj_pass_filter_common_cells <- addMotifAnnotations(ArchRProj = proj_pass_filter_common_cells, motifPWMs = motif_all,
                                                     annoName = "TF-Motif", force = T)
motifPositions <- getPositions(proj_pass_filter_common_cells, name = "TF-Motif")
# proj_pass_filter_common_cells <- addGroupCoverages(ArchRProj = proj_pass_filter_common_cells,
#                                                    groupBy = "Celltype", force = TRUE)
seFoot <- getFootprints(
  ArchRProj = proj_pass_filter_common_cells, 
  positions = motifPositions, 
  groupBy = "Celltype"
)

keep <- names(seFoot@assays@data@listData) %in% RNA_celltype_regulon$TF
keep <- names(seFoot@assays@data@listData)[keep]

plotFootprints(
  seFoot = seFoot,
  names = keep,
  ArchRProj = proj_pass_filter_common_cells, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias_2",
  addDOC = FALSE,
  smoothWindow = 5,
  force = T,
  #pal = RNA_celltype_color[common_cells_type],
  pal = RNA_celltype_color["Spikelet meristem (SM)"],
  baseSize = 10
)
unique(proj_pass_filter$Celltype)

RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os01g0229000"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os01g0195000"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os09g0475400"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os02g0635800"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os01g0883100"),]

RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os09g0475400"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os06g0164400"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os04g0671900"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os03g0782500"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os02g0766700"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os02g0635800"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os01g0924400"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os01g0883100"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os12g0233800"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF == "Os01g0924400"),]

RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "OsMYB2P-1"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "WRKY45"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "OSH1"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "OsTCP7"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "OsMYB3R1-L"),]
RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "OsbHLH032"),]

############
addMotifAnnotations <- function(
    ArchRProj = NULL,
    motifSet = "cisbp",
    annoName = "Motif",
    species = NULL,
    collection = "CORE",
    motifPWMs = NULL,
    cutOff = 5e-05, 
    width = 7,
    version = 2,
    force = FALSE,
    peak_selected = NULL,
    logFile = createLogFile("addMotifAnnotations"),
    ...
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = motifSet, name = "motifSet", valid = c("character", "null"))
  .validInput(input = annoName, name = "annoName", valid = c("character"))
  .validInput(input = species, name = "species", valid = c("character", "null"))
  .validInput(input = collection, name = "collection", valid = c("character", "null"))
  .validInput(input = cutOff, name = "cutOff", valid = c("numeric"))
  .validInput(input = width, name = "width", valid = c("integer"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  if(!is.null(motifPWMs)){
    if(!is(motifPWMs, "PWMatrixList")){
      stop("User Supplied motifPWMS must be a PWMatrixList!")
    }
    motifSet <- "Custom"
  }
  
  if(is.null(motifSet)){
    stop("Must provide motifSet or motifPWMs!")
  }
  
  .requirePackage("motifmatchr", installInfo='BiocManager::install("motifmatchr")')
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addMotifAnnotations Input-Parameters", logFile = logFile)
  
  if(annoName %in% names(ArchRProj@peakAnnotation)){
    if(force){
      message("peakAnnotation name already exists! Overriding.")
    }else{
      stop("peakAnnotation name already exists! set force = TRUE to override!")
    }
  }
  
  if(grepl("JASPAR|CISBP", motifSet, ignore.case = TRUE) & is.null(species)){
    if(grepl("hg19",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Homo sapiens"
    }
    if(grepl("hg38",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Homo sapiens"
    }
    if(grepl("mm9",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Mus musculus"
    }
    if(grepl("mm10",getGenomeAnnotation(ArchRProj)$genome, ignore.case = TRUE)){
      species <- "Mus musculus"
    }
  }
  
  #############################################################
  # Get PWM List adapted from chromVAR!
  #############################################################
  
  .logDiffTime(paste0("Gettting Motif Set, Species : ", species), t1 = tstart, verbose = TRUE, logFile = logFile)
  
  if(tolower(motifSet)=="jaspar2020"){
    
    .requirePackage("JASPAR2020",installInfo='BiocManager::install("JASPAR2020")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
    
  }else if(tolower(motifSet)=="jaspar2018"){
    
    .requirePackage("JASPAR2018",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
    
  }else if(tolower(motifSet)=="jaspar2016"){
    
    .requirePackage("JASPAR2016",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
    
  }else if(tolower(motifSet)=="cisbp"){
    
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    if(tolower(species) == "mus musculus"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("mouse_pwms_v1")
        motifs <- mouse_pwms_v1        
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("mouse_pwms_v2")
        motifs <- mouse_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else if(tolower(species) == "homo sapiens"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("human_pwms_v1")
        motifs <- human_pwms_v1        
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("human_pwms_v2")
        motifs <- human_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else{
      stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
    }
    
  }else if(tolower(motifSet)=="encode"){
    
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("encode_pwms")
    motifs <- encode_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
    
  }else if(tolower(motifSet)=="homer"){
    
    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("homer_pwms")
    motifs <- homer_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary
    
  }else if(tolower(motifSet)=="vierstra"){
    if(tolower(collection)=="individual"){
      url = "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Individual_Motifs.rds"
      message("Using Vierstra v1.0 motifs. See https://www.vierstra.org/resources/motif_clustering for more details.")
    } else if(tolower(collection == "archetype")){
      url = "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Archetype_Motifs_v2.1.rds"
      message("Using Vierstra v2.1beta motifs. See https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/ for more details.")
    } else {
      stop(paste0("Error! collection ", collection, " not recognized for motifSet ",motifSet,
                  ". Accepted values are 'individual' and 'archetype'"))
    }
    
    annoPath <- file.path(find.package("ArchR", NULL, quiet = TRUE), "data", "Annotations")
    dir.create(annoPath, showWarnings = FALSE)
    
    #Download
    if(!file.exists(file.path(annoPath, basename(url)))){
      message("Motif file ", basename(url)," does not exist! Downloading..")
      download.file(
        url = url, 
        destfile = file.path(annoPath, basename(url)),
        quiet = FALSE
      )
    }
    motifFile <- file.path(annoPath, basename(url))
    
    motifs <- readRDS(motifFile)
    obj <- NULL
    motifSummary <- NULL
    
  }else if(tolower(motifSet)=="custom"){
    
    obj <- NULL
    motifs <- motifPWMs
    motifSummary <- NULL
    
  }else{
    
    stop("Error MotifSet Not Recognized!")
    
  }
  
  .logThis(motifs, "motifs", logFile = logFile)
  .logThis(obj, "obj", logFile = logFile)
  .logThis(motifSummary, "motifSummary", logFile = logFile)
  
  #############################################################
  # Get BSgenome Information!
  #############################################################
  genome <- ArchRProj@genomeAnnotation$genome
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)
  
  #############################################################
  # Calculate Motif Positions
  #############################################################
  .logDiffTime("Finding Motif Positions with motifmatchr!", t1 = tstart, verbose = TRUE, logFile = logFile)
  peakSet <- ArchRProj@peakSet
  if(is.null(peakSet)){
    .logStop("peakSet is NULL. You need a peakset to run addMotifAnnotations! See addReproduciblePeakSet!", logFile = logFile)
  }
  if(!is.null(peak_selected)){
    peakSet <- peakSet[peak_selected,]
  }
  motifPositions <- motifmatchr::matchMotifs(
    pwms = motifs,
    subject = peakSet,
    genome = BSgenome, 
    out = "positions", 
    p.cutoff = cutOff, 
    w = width
  )
  
  #############################################################
  # Filter Motifs With No Matches
  #############################################################
  
  #Number of Overlaps
  nO <- lapply(motifPositions, length) %>% unlist
  mF <- names(which(nO == 0))
  
  if(all(nO == 0)){
    stop("No Overlaps Found! Please check your peakSet and genome!")
  }
  
  if(length(mF) > 0){
    .logDiffTime(paste0("Filtering Motif Annotations with 0 overlaps :\n\n ", paste(mF, collapse=", "), "\n\n"), t1 = tstart, verbose = TRUE, logFile = logFile)
    #Filter
    motifPositions <- motifPositions[nO > 0]
    motifSummary <- motifSummary[names(motifPositions),,drop=FALSE]
    motifs <- motifs[names(motifPositions)]
  }else{
    .logDiffTime(paste0("All Motifs Overlap at least 1 peak!"), t1 = tstart, verbose = TRUE, logFile = logFile)
  }  
  
  #############################################################
  # Motif Overlap Matrix
  #############################################################
  .logDiffTime("Creating Motif Overlap Matrix", t1 = tstart, verbose = TRUE, logFile = logFile)
  allPositions <- unlist(motifPositions)
  overlapMotifs <- findOverlaps(peakSet, allPositions, ignore.strand=TRUE)
  motifMat <- Matrix::sparseMatrix(
    i = queryHits(overlapMotifs),
    j = match(names(allPositions),names(motifPositions))[subjectHits(overlapMotifs)],
    x = rep(TRUE, length(overlapMotifs)),
    dims = c(length(peakSet), length(motifPositions))
  )
  colnames(motifMat) <- names(motifPositions)
  motifMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = motifMat), rowRanges = peakSet)
  .logDiffTime("Finished Getting Motif Info!", t1 = tstart, verbose = TRUE, logFile = logFile)
  
  out <- SimpleList(
    motifSummary = motifSummary,
    motifMatches = motifMat,
    motifPositions = motifPositions,
    motifList = motifs,
    date = Sys.Date()
  )
  
  dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), showWarnings=FALSE)
  savePositions <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(annoName,"-Positions-In-Peaks.rds"))
  saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(annoName,"-Matches-In-Peaks.rds"))
  
  ArchRProj@peakAnnotation[[annoName]]$Name <- annoName
  ArchRProj@peakAnnotation[[annoName]]$motifs <- motifs
  ArchRProj@peakAnnotation[[annoName]]$motifSummary <- motifSummary
  ArchRProj@peakAnnotation[[annoName]]$Positions <- savePositions
  ArchRProj@peakAnnotation[[annoName]]$Matches <- saveMatches
  
  .safeSaveRDS(out, file.path(getOutputDirectory(ArchRProj),  "Annotations", paste0(annoName,"-In-Peaks-Summary.rds")), compress = FALSE)
  .safeSaveRDS(out$motifPositions, savePositions, compress = FALSE)
  .safeSaveRDS(out$motifMatches, saveMatches, compress = FALSE)
  
  .endLogging(logFile = logFile)
  
  return(ArchRProj)
  
}



gene_id_map <- readRDS("gene_id_map.rds")
ArchR_temp <- subsetArchRProject(proj_pass_filter_common_cells,
                                 cells = getCellNames(proj_pass_filter_common_cells),
                                 outputDirectory = "TEMP"
)
set.seed(123)
ArchR_temp <- addGroupCoverages(ArchRProj = ArchR_temp,
                                groupBy = "Celltype", force = TRUE)
ArchR_temp_peakSets <- getPeakSet(ArchR_temp)
ArchR_temp_peakSets$Width <- ArchR_temp_peakSets@ranges@width
ArchR_temp_peakSets$Celltype <- ArchR_temp_peakSets@ranges@NAMES
ArchR_temp_peakSets$Chr <- as.character(ArchR_temp_peakSets@seqnames)
ArchR_temp_peakSets$Start <- ArchR_temp_peakSets@ranges@start
ArchR_temp_peakSets$ID <- paste0(ArchR_temp_peakSets$Celltype, " - ", ArchR_temp_peakSets$Chr, " - ", ArchR_temp_peakSets$Start)
ArchR_temp_peakSets$End <- ArchR_temp_peakSets$Start + ArchR_temp_peakSets$Width
ArchR_temp_peakSets <- as.data.frame(ArchR_temp_peakSets@elementMetadata)
rownames(ArchR_temp_peakSets) <- ArchR_temp_peakSets$ID

temp_regulon <- "Os12g0560900(35g)"
remp_regulon_list <- RNA_celltype_regulon[which(RNA_celltype_regulon$Regulon == temp_regulon), "Target_list"]
remp_regulon_list <- unlist(strsplit(remp_regulon_list, " | ", fixed = T))
remp_regulon_list <- gene_id_map[match(remp_regulon_list, gene_id_map$gene_id), "symbol"]
remp_regulon_list_peakSets <- ArchR_temp_peakSets[which(ArchR_temp_peakSets$nearestGene %in% remp_regulon_list),]
remp_regulon_list_peakSets <- remp_regulon_list_peakSets[which(remp_regulon_list_peakSets$peakType == "Promoter"),]
unique(remp_regulon_list_peakSets$nearestGene)
peak_selected <- match(remp_regulon_list_peakSets$ID, ArchR_temp_peakSets$ID)
# remp_regulon_list_peakSets <- remp_regulon_list_peakSets[order(remp_regulon_list_peakSets$Chr, remp_regulon_list_peakSets$Start),]
# remp_regulon_list_peakSets_GR <- GRanges(seqnames = remp_regulon_list_peakSets$Chr,
#                                          IRanges(start = remp_regulon_list_peakSets$Start,
#                                                  end = remp_regulon_list_peakSets$End),
#                                          gene = remp_regulon_list_peakSets$nearestGene)
# remp_regulon_list_peakSets_GR <- GRangesList(remp_regulon_list_peakSets_GR)


# ArchR_temp <- addPeakAnnotations(
#   ArchRProj = ArchR_temp,
#   regions = remp_regulon_list_peakSets_GR,
#   name = "Os12g0560900(35g)"
# )
# 
proj_pass_filter_common_cells <- addMotifAnnotations(ArchRProj = proj_pass_filter_common_cells, motifPWMs = motif_all["Os12g0560900"],
                                  annoName = "Os12g0560900(35g)", force = T, peak_selected = peak_selected)
motifPositions <- getPositions(proj_pass_filter_common_cells, name = "Os12g0560900(35g)")
seFoot <- getFootprints(
  ArchRProj = proj_pass_filter_common_cells, 
  positions = motifPositions, 
  groupBy = "Celltype"
)

keep <- names(seFoot@assays@data@listData) %in% RNA_celltype_regulon$TF
keep <- names(seFoot@assays@data@listData)[keep]

plotFootprints(
  seFoot = seFoot,
  names = keep,
  ArchRProj = proj_pass_filter_common_cells, 
  normMethod = "Subtract",
  addDOC = FALSE,
  smoothWindow = 5,
  force = T,
  pal = RNA_celltype_color[common_cells_type],
  baseSize = 10,
  plotName = keep
)

#### sequence logo
motif_pfm <- read_meme("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs.meme")
motif_pfm <- convert_motifs(motif_pfm, class = "TFBSTools-PFMatrix")
motif_pfm <- do.call(PFMatrixList, motif_pfm)
# TF_MSU_RAP <- c[,c("TF_MSU", "TF_RAP")]
# TF_MSU_RAP <- TF_MSU_RAP[!duplicated(TF_MSU_RAP),]
motif_pfm_MSU_ID <- unlist(lapply(motif_pfm@listData, function(x){
  x@name
}))
names(motif_pfm) <- motif_pfm_MSU_ID
motif_pfm_RAP_ID <- match(motif_pfm_MSU_ID, MSU_RAP$MSU)
motif_pfm_RAP_ID <- MSU_RAP$RAP[motif_pfm_RAP_ID]
motif_pfm <- motif_pfm[-which(is.na(motif_pfm_RAP_ID))]
motif_pfm_RAP_ID <- motif_pfm_RAP_ID[-which(is.na(motif_pfm_RAP_ID))]
names(motif_pfm) <- motif_pfm_RAP_ID
motif_pfm <- motif_pfm[!duplicated(names(motif_pfm))] # ID duplicated !!!


library(ggseqlogo)
pdf("Figure3_seqlogo.pdf", width = 10, height = 4)
ggseqlogo(list(Os01g0229000 = motif_pfm@listData$Os01g0229000@profileMatrix,
               Os01g0195000 = motif_pfm@listData$Os01g0195000@profileMatrix,
               Os09g0475400 = motif_pfm@listData$Os09g0475400@profileMatrix,
               Os02g0635800 = motif_pfm@listData$Os02g0635800@profileMatrix))
dev.off()



temp1 <- data.frame(start = RNA_celltype_regulon$TF,
                    end = RNA_celltype_regulon$Celltype)
temp2 <- data.frame(start = RNA_celltype_regulon$Celltype,
                    end = RNA_celltype_regulon$Regulon)
df <- as.data.frame(rbind(temp1,
                          temp2))


ggplot(data = RNA_celltype_regulon, aes(x = TF, y = Celltype, fill = log2(Target_num+1))) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black")) +
  scale_fill_gradientn(colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(20),)


length(unique(RNA_celltype_regulon$Regulon))
ggplot(data = RNA_celltype_regulon, aes(x = TF, y = Celltype, fill = Regulon)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black")) +
  scale_fill_manual(values = scales::hue_pal()(220)) +
  NoLegend()

pdf("Figure3_C.pdf", width = 16, height = 6)
ggplot(data = RNA_celltype_regulon, aes(x = TF_Symbol, y = Target_num, fill = Celltype)) +
  geom_point(size = 3, shape = 21, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black")) +
  scale_fill_manual(values = RNA_celltype_color) +
  labs(x = "TF", y = "Number of target genes of TF in cell-types") +
  NoLegend() +
  ggbreak::scale_y_break(c(115, 170), # 截断位置及范围
                space = 0.3, # 间距大小
                scales = 0.5) # 上下显示比例，大于1上面比例大，小于1下面比例大
dev.off()




library(ggalluvial)
RNA_celltype_regulon
df <- to_lodes_form(RNA_celltype_regulon[which(RNA_celltype_regulon$TF_Symbol == "OSH1"),c(9,2,10)],
                    axes = 1:3,
                    id = "value")
ggplot(df, aes(x = x, fill = stratum, label = stratum,
               stratum = stratum, alluvium  = value))+#数据
  geom_flow(width = 0.3,#连线宽度
            curve_type = "sine",#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
            alpha = 0.5,#透明度
            color = 'white',#间隔颜色
            size = 0.1)+#间隔宽度
  geom_stratum(width = 0.28)+#图中方块的宽度
  geom_text(stat = 'stratum', size = 5, color = 'black')+
  # scale_fill_manual(values = col)+#自定义颜色
  theme_void()+#主题（无轴及网格线）
  theme(legend.position = 'none')#去除图例



RNA_meta <- RNA@meta.data

M1 <- RNA_meta[which(RNA_meta$Celltype %in% c("Large parenchyma (MO)", "Mesophyll (MO)",
                                              "Procambium", "Fiber", "Phloem")),]
table(M1$tissues)
pdf("Figure3_M1_tissue_pie.pdf", width = 3.5, height = 3.5)
pie(c(8361,
      114,
      20431+580),
    labels = c("Leaf (28.3%)",
               "Root (0.4%)",
               "Flag leaf (71.3%)"),
    col = c("#488e44", "#867c58", "#A4738B"))
dev.off()

M2 <- RNA_meta[which(RNA_meta$Celltype %in% c("Intermediate anther", "Epidermis", "Late anther", "Mesophyll initial",
                                              "Mesophyll precursor", "Upper protoderm", "Suspensor", "Vascular initial",
                                              "Dividing inner", "Shoot meristematic cell", "Columella initials", "Dividing outer")),]
table(M2$tissues)
pdf("Figure3_M2_tissue_pie.pdf", width = 3.5, height = 3.5)
pie(c(199,
      8569,
      4306,
      1375,
      4003),
    labels = c("Leaf (1.1%)",
               "Flag leaf (46.4%)",
               "Mature pollen (23.3%)",
               "SAM (7.5%)",
               "Zygote (21.7%)"),
    col = c("#488e44", "#A4738B", "#4E80A2", "#B44562", "#BE853B"))
dev.off()

M3 <- RNA_meta[which(RNA_meta$Celltype %in% c("Epidermal cell", "Shoot endodermis", "Inflorescence meristem(IM)",
                                              "Proliferating cell", "Shoot meristematic cell", "Cryptic bract/bract(cb/b)",
                                              "Lemma(le)", "Rachis", "Spikelet meristem(SM)", "Branch meristems(BM)", "Other Spikelet")),]
table(M3$tissues)
pdf("Figure3_M3_tissue_pie.pdf", width = 3.5, height = 3.5)
pie(c(12445 + 14196,
      6797),
    labels = c("IM (79.7%)",
               "SAM (20.3%)"),
    col = c("#F7786B", "#B44562"))
dev.off()

M4 <- RNA_meta[which(RNA_meta$Celltype %in% c("Endodermis-cortical initial",
                                              "Cortex", "Endodermis")),]
table(M4$tissues)
pdf("Figure3_M4_tissue_pie.pdf", width = 3.5, height = 3.5)
pie(c(16130,
      12208),
    labels = c("Root (56.9%)",
               "Root tip (43.1%)"),
    col = c("#867c58", "#A2A882"))
dev.off()


################ 废弃代码
if (F) {
  # Mesophyll initial
  c("specification of floral organ identity", "floral meristem determinacy", "", "")
  
  temp <- as.data.frame.array(table(OsMADS2_M2$Celltype, OsMADS2_M2$Targets))
  sort(colSums(temp))
  library(UpSetR)
  OsMADS2_M2_celltype_target_list <- c()
  for (i in unique(OsMADS2_M2$Celltype)) {
    temp <- OsMADS2_M2[which(OsMADS2_M2$Celltype == i),]
    OsMADS2_M2_celltype_target_list <- c(OsMADS2_M2_celltype_target_list,
                                         list(temp$Targets))
  }
  names(OsMADS2_M2_celltype_target_list) <- unique(OsMADS2_M2$Celltype)
  
  OsMADS2_M2_celltype_target_list <- list_to_matrix(OsMADS2_M2_celltype_target_list)
  pdf("Figure3_OsMADS2_M2_celltype_target_UpSetPlot.pdf", height = 5, width = 12)
  upset(as.data.frame(OsMADS2_M2_celltype_target_list), nsets = 8, nintersects = NA, mb.ratio = c(0.5, 0.7),
        set_size.show = T, point.size = 3, text.scale = 2.5,
        order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
  dev.off()
  
  
  # M3举例
  temp <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Module == "Module 3"),]
  OsBZIP42_M3 <- temp[which(temp$TF_Symbol == "OsBZIP42"),]
  genes_list <- c()
  for (i in unique(OsBZIP42_M3$Celltype)) {
    genes_list <- c(genes_list, list(OsBZIP42_M3$Targets[which(OsBZIP42_M3$Celltype == i)]))
  }
  names(genes_list) <- unique(OsBZIP42_M3$Celltype)
  GO_result_OsBZIP42_M3 <- c()
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
    GO_result$Celltype <- i
    GO_result_OsBZIP42_M3 <- as.data.frame(rbind(GO_result_OsBZIP42_M3, GO_result))
    
  }
  GO_result_OsBZIP42_M3 <- GO_result_OsBZIP42_M3[which(GO_result_OsBZIP42_M3$P_value < 0.05),]
  unique(GO_result_OsBZIP42_M3$Celltype)
  GO_result_OsBZIP42_M3 <- GO_result_OsBZIP42_M3[order(GO_result_OsBZIP42_M3$P_value, decreasing = F),]
  GO_result_OsBZIP42_M3 <- GO_result_OsBZIP42_M3[which(GO_result_OsBZIP42_M3$Ontology == "BP"),]
  GO_result_OsBZIP42_M3_IM <- GO_result_OsBZIP42_M3[which(GO_result_OsBZIP42_M3$Celltype == "Inflorescence meristem (IM)"),]
  GO_result_OsBZIP42_M3_Epi <- GO_result_OsBZIP42_M3[which(GO_result_OsBZIP42_M3$Celltype == "Epidermal cell"),]
  GO_result_OsBZIP42_M3_Pro <- GO_result_OsBZIP42_M3[which(GO_result_OsBZIP42_M3$Celltype == "Proliferating cell"),]
  GO_result_OsBZIP42_M3_SM <- GO_result_OsBZIP42_M3[which(GO_result_OsBZIP42_M3$Celltype == "Shoot meristematic cell"),]
  
  
  
  
  
  
  temp <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Module == "Module 3"),]
  temp <- as.data.frame.array(table(temp$Celltype, temp$Targets))
  temp <- as.data.frame(t(temp))
  write.csv(temp, "Figure3_temp.csv")
  
  as.data.frame.array(table(temp$Celltype, temp$TF_Symbol))
  # OsMADS2, Os01g0883100
  Annotation_genes[which(Annotation_genes$`CGSNL Gene Symbol` == "OsBZIP42"),]
  
  {
    temp <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Module == "Module 3"),]
    OsBZIP42_M3 <- temp[which(temp$TF_Symbol == "OsBZIP42"),]
    temp <- as.data.frame.array(table(OsBZIP42_M3$Celltype, OsBZIP42_M3$Targets))
    sort(colSums(temp))
    library(UpSetR)
    OsBZIP42_M3_celltype_target_list <- c()
    for (i in unique(OsBZIP42_M3$Celltype)) {
      temp <- OsBZIP42_M3[which(OsBZIP42_M3$Celltype == i),]
      OsBZIP42_M3_celltype_target_list <- c(OsBZIP42_M3_celltype_target_list,
                                            list(temp$Targets))
    }
    names(OsBZIP42_M3_celltype_target_list) <- unique(OsBZIP42_M3$Celltype)
    
    OsBZIP42_M3_celltype_target_list <- list_to_matrix(OsBZIP42_M3_celltype_target_list)
    # pdf("Figure1_tissue_peakGene_UpSetPlot.pdf", height = 5, width = 20)
    upset(as.data.frame(OsBZIP42_M3_celltype_target_list), nsets = 8, nintersects = NA, mb.ratio = c(0.5, 0.5), set_size.show = T, point.size = 5, text.scale = 2.5,
          order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
    temp <- as.data.frame.array(table(OsBZIP42_M3$Targets, OsBZIP42_M3$Celltype))
    temp$target_symbol <- gene_id_map[rownames(temp), "symbol"]
    write.csv(temp, "Figure3_temp.csv")
  }
  
  # 
  
  
  
  OsMADS2_M2 <- temp[which(temp$TF_Symbol == "OsMADS2"),]
  temp <- as.data.frame.array(table(OsMADS2_M2$Celltype, OsMADS2_M2$Targets))
  sort(colSums(temp))
  library(UpSetR)
  OsMADS2_M2_celltype_target_list <- c()
  for (i in unique(OsMADS2_M2$Celltype)) {
    temp <- OsMADS2_M2[which(OsMADS2_M2$Celltype == i),]
    OsMADS2_M2_celltype_target_list <- c(OsMADS2_M2_celltype_target_list,
                                         list(temp$Targets))
  }
  names(OsMADS2_M2_celltype_target_list) <- unique(OsMADS2_M2$Celltype)
  
  OsMADS2_M2_celltype_target_list <- list_to_matrix(OsMADS2_M2_celltype_target_list)
  # pdf("Figure1_tissue_peakGene_UpSetPlot.pdf", height = 5, width = 20)
  upset(as.data.frame(OsMADS2_M2_celltype_target_list), nsets = 8, nintersects = NA, mb.ratio = c(0.5, 0.5), set_size.show = T, point.size = 5, text.scale = 2.5,
        order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
  dev.off()
  
  
  RNA_mean <- matrix(0, nrow = length(unique(rownames(RNA))), ncol = length(unique(RNA$Celltype)))
  rownames(RNA_mean) <- unique(rownames(RNA))
  colnames(RNA_mean) <- unique(RNA$Celltype)
  for (i in unique(RNA$Celltype)) {
    temp <- subset(RNA, subset = Celltype == i)
    temp_sum <- rowSums(temp@assays$RNA@data)
    temp_mean <- temp_sum / ncol(temp)
    RNA_mean[,i] <- temp_mean[rownames(RNA_mean)]
  }
  temp <- RNA_celltype_regulon_pairs[which(RNA_celltype_regulon_pairs$Module == "Module 2"),]
  OsMADS2_M2 <- temp[which(temp$TF_Symbol == "OsMADS2"),]
  OsMADS2_Os05g0557600 <- data.frame(Celltype = unique(OsMADS2_M2$Celltype))
  target_expression_mean <- c()
  for (i in 1:nrow(OsMADS2_Os05g0557600)) {
    target_expression_mean <- c(target_expression_mean,
                                RNA_mean["Os05g0557600", OsMADS2_Os05g0557600$Celltype[i]])
  }
  tf_expression_mean <- c()
  for (i in 1:nrow(OsMADS2_Os05g0557600)) {
    tf_expression_mean <- c(tf_expression_mean,
                            RNA_mean["Os01g0883100", OsMADS2_Os05g0557600$Celltype[i]])
  }
  
  OsMADS2_Os05g0557600$tf_expression_mean <- tf_expression_mean
  OsMADS2_Os05g0557600$target_expression_mean <- target_expression_mean
  OsMADS2_Os05g0557600$tf_target <- OsMADS2_Os05g0557600$tf_expression_mean + OsMADS2_Os05g0557600$target_expression_mean
  
  temp <- data.frame(word = OsMADS2_Os05g0557600$Celltype,
                     freq = OsMADS2_Os05g0557600$tf_target)
  library(jiebaR)
  library(wordcloud2)
  my_graph <- wordcloud2(
    data = temp, size = 0.2, shape = 'star'
  )
  my_graph
  
  temp <- as.data.frame.array(table(RNA$tissues, RNA$Celltype))
  temp <- temp[,unique(OsMADS2_M2$Celltype)]
  
  OsMADS2_M2$Target_symbol <- gene_id_map[OsMADS2_M2$Targets, "symbol"]
  temp <- as.data.frame.array(table(OsMADS2_M2$Celltype, OsMADS2_M2$Target_symbol))
  temp <- as.data.frame(t(temp))
  write.csv(temp, "Figure3_temp.csv")
  temp <- data.frame(tissue = c("Zygote", "Zygote", "MaturePollen", "MaturePollen", "MaturePollen", "FlagLeaf", "FlagLeaf", "Zygote"),
                     temp)
  temp <- aggregate.data.frame(temp[,-1], by = list(temp$tissue), FUN = sum)
  rownames(temp) <- temp$Group.1
  temp <- temp[,-1]
  temp <- as.data.frame(t(temp))
  temp$gene <- rownames(temp)
  write.csv(temp, "Figure3_temp.csv")
  gene_id_map[which(gene_id_map$symbol == "NIP3"),]
  gene_id_map[which(gene_id_map$symbol == "OsCW.ZF6"),]
}
