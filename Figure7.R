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
library(ggExtra)
library(VennDiagram)
library(ggpubr)
library(qqman)
library(CMplot)
library(goftest)

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

# load("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Figure6_2.RData")
enrichment_p <- function(pathway_genes, diff_genes) {
  x <- intersect(pathway_genes, diff_genes)
  x <- length(x) - 1
  m <- length(unique(pathway_genes))
  n <- 37960 - m
  k <- length(unique(diff_genes))
  phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE)
}

pheno <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/QC/pheno.txt", sep = " ", header = T)
rownames(pheno) <- pheno$FID
sample_location <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/sample_location.csv",
                            sep = ",", header = T, check.names = F)
table(sample_location$Subpopulation)
rownames(sample_location) <- sample_location$`Cultivar ID`
LAMBDA <- c()

meta_name <- colnames(pheno)[-c(1,2)]
# Heading_date
m <- 1
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 25)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 60 & temp < 130)]
  hist(temp)
  dev.off()
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
                    )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
                    )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 12308486 rice4k_sample_filtered_SNP_filtered.map
  # 145 rice4k_sample_filtered_SNP_filtered.ped
  # 12308631 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于1e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 37194 variants and 145 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  dev.off()
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                                   sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 0.00005)
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Plant_height
m <- 2
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 50)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp <- temp[which(temp > 50)]
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 12310202 rice4k_sample_filtered_SNP_filtered.map
  # 158 rice4k_sample_filtered_SNP_filtered.ped
  # 12310360 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48520 variants and 157 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Num_panicles
m <- 3
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 0 & temp < 30)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 12310202 rice4k_sample_filtered_SNP_filtered.map
  # 158 rice4k_sample_filtered_SNP_filtered.ped
  # 12310360 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48520 variants and 157 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Num_effective_panicles
m <- 4
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 0 & temp < 30)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 12310202 rice4k_sample_filtered_SNP_filtered.map
  # 158 rice4k_sample_filtered_SNP_filtered.ped
  # 12310360 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 38498 variants and 158 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Yield
m <- 5
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 0 & temp < 60)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 12315310 rice4k_sample_filtered_SNP_filtered.map
  # 157 rice4k_sample_filtered_SNP_filtered.ped
  # 12315467 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48373 variants and 156 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  dev.off()
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Grain_weight
m <- 6
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 10 & temp < 35)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 14270639 rice4k_sample_filtered_SNP_filtered.map
  # 150 rice4k_sample_filtered_SNP_filtered.ped
  # 14270789 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 46751 variants and 156 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  dev.off()
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Spikelet_length
m <- 7
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 10)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 14190792 rice4k_sample_filtered_SNP_filtered.map
  # 157 rice4k_sample_filtered_SNP_filtered.ped
  # 14190949 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48520 variants and 157 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  dev.off()
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Grain_length
m <- 8
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 5)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 14190792 rice4k_sample_filtered_SNP_filtered.map
  # 157 rice4k_sample_filtered_SNP_filtered.ped
  # 14190949 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48520 variants and 157 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  dev.off()
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Grain_width
m <- 9
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 0)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 14190792 rice4k_sample_filtered_SNP_filtered.map
  # 157 rice4k_sample_filtered_SNP_filtered.ped
  # 14190949 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48520 variants and 157 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  dev.off()
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}
# Grain_thickness
m <- 10
{
  # 表型正态分布检验
  temp <- pheno[,meta_name[m]]
  names(temp) <- pheno$FID
  hist(temp, breaks = 50)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  quantile(temp)
  temp <- temp[which(temp > 0)]
  hist(temp)
  ks.test(temp, "pnorm", mean = mean(temp), sd = sqrt(var(temp)))
  cvm.test(temp)
  temp_pheno <- pheno[names(temp),]
  temp_sample <- temp_pheno[,c(1,2)]
  dir.create(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  path <- paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m], "/")
  write.table(temp_sample, paste0(path, "sample_keep.txt"),
              sep = " ", quote = F, row.names = F, col.names = F)
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # 提取样本
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_onlySNP ",
                    "--keep ", "sample_keep.txt ",
                    "--recode --out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 格式转换
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--out rice4k_sample_filtered_2"
  )
  system(command)
  
  # 查看个体缺失的位点数，每个SNP缺失的个体数目
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--make-bed ",
                    "--missing"
  )
  system(command)
  
  # 删除SNP缺失率>0.1的个体
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_2 ",
                    "--mind 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_0"
  )
  system(command)
  
  # 删除缺失率>0.1的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_0 ",
                    "--geno 0.1 ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--recode --out rice4k_sample_filtered_SNP_filtered "
  )
  system(command)
  
  # 统计样本数目和SNP数量
  command <- paste0("wc -l rice4k_sample_filtered_SNP_filtered.map rice4k_sample_filtered_SNP_filtered.ped"
  )
  system(command)
  # 14190792 rice4k_sample_filtered_SNP_filtered.map
  # 157 rice4k_sample_filtered_SNP_filtered.ped
  # 14190949 总计
  
  # 计算SNP位点基因频率
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--freq --out MAF_check"
  )
  system(command)
  
  # 去除MAF小于0.05的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered ",
                    "--maf 0.05 --make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered"
  )
  system(command)
  
  # 去除相互关联的SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--indep-pairwise 500 50 0.5 --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 基于LD信息过滤SNP
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered ",
                    "--extract rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep.prune.in  --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 格式转化
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "-make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep"
  )
  system(command)
  
  # 计算HWE平衡p值
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--hardy"
  )
  system(command)
  
  # 选择HWE p值小于5e-5的SNP
  command <- paste0("awk '{ if ($9 < 0.00005) print $0 }' plink.hwe > plink_HWE_filtered.hwe"
  )
  system(command)
  command <- paste0("awk '{ print $2 }' plink_HWE_filtered.hwe > plink_HWE_filtered.hwe.SNP"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep ",
                    "--extract plink_HWE_filtered.hwe.SNP --recode --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--make-bed --out rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final"
  )
  system(command)
  # 48520 variants and 157 people pass filters and QC
  
  # 样本PCA
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--pca 10"
  )
  system(command)
  
  # 选择PC数量
  pca <- read.csv(paste0(path, "plink.eigenvec"), sep = " ", header = F)
  temp_location <- sample_location[pca$V1,]
  temp_pca_location <- data.frame(pca[,-c(1,2)],
                                  population = temp_location$Subpopulation)
  ggplot(data = temp_pca_location, aes(x = V3, y = V4, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V3, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V4, y = V5, color = population)) +
    geom_point()
  ggplot(data = temp_pca_location, aes(x = V5, y = V6, color = population)) +
    geom_point()
  
  pca_temp <- data.frame(1,
                         pca[,c("V3","V4")])
  write.table(pca_temp, paste0(path, "cov.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pca_temp <- data.frame(pca[,c("V1","V2","V3","V4")])
  colnames(pca_temp) <- c("FID", "IID", "COV1", "COV2")
  write.table(pca_temp, paste0(path, "cov_plink.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,meta_name[m]], paste0(path, "pheno.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  pheno_temp <- pheno[pca$V1,]
  write.table(pheno_temp[,c("FID", "IID", meta_name[m])], paste0(path, "pheno_plink.txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  setwd(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/", meta_name[m]))
  
  # Plink LM线性模型
  command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                    "--file ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "--allow-no-sex --adjust --linear --pheno pheno_plink.txt --all-pheno --covar cov_plink.txt --out plink_result"
  )
  system(command)
  
  
  # GEMMA LMM线性混合模型
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-gk 2 -p pheno.txt"
  )
  system(command)
  
  command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
                    "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                    "-k ./output/result.sXX.txt -lmm 1 -p pheno.txt -c cov.txt"
  )
  system(command)
  
  # command <- paste0("/media/heshidian/RAID5_42TB/3.Software/gemma-0.98.1-linux-static ",
  #                   "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
  #                   "-lm 1 -p pheno.txt -c cov.txt"
  # )
  # system(command)
  
  # 曼哈顿图、qqplot、膨胀系数
  # GEMMA
  temp <- data.table::fread(paste0(path, "output/result.assoc.txt"),
                            sep = "\t", header = T)
  temp <- as.data.frame(temp)
  temp2 <- temp[,c("rs", "chr", "ps", "p_wald")]
  colnames(temp2) <- c("SNP", "CHR", "BP", "P")
  sum(temp2$P < 5e-5)
  dev.off()
  pdf(paste0(path, "manhattan_plot.pdf"), width = 8, height = 6.5)
  CMplot(temp2, plot.type = "m", threshold = c(0.00005),
         amplify = T, signal.cex = c(1,1), signal.pch = c(20,20),
         signal.col = c("red","blue"), multracks = F, file.output = F, file = "pdf")
  dev.off()
  CMplot(temp2, plot.type = "q", threshold = 0.05, main = meta_name[m], main.cex = 2, file = "pdf")
  p_value <- temp2$P
  chisq <- qchisq(1-p_value, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  LAMBDA <- c(LAMBDA, lambda)
}

names(LAMBDA) <- meta_name


Chr_length <- as.data.frame(genome_annotation@listData[["chromSizes"]])
rownames(Chr_length) <- Chr_length$seqnames

### GWAS atlas
GWAS_atlas_vcf_files <- list.files("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs/GWAS_atlas/vcf",
                                  full.names = F, recursive = F)
GWAS_atlas_vcf <- c()
for (i in GWAS_atlas_vcf_files) {
  print(i)
  temp_vcf <- data.table::fread(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs/GWAS_atlas/vcf/", i),
                                sep = "\t")
  GWAS_atlas_vcf <- as.data.frame(rbind(GWAS_atlas_vcf,
                                        temp_vcf))
}
rm(temp_vcf)
table(GWAS_atlas_vcf$REF)
table(GWAS_atlas_vcf$ALT)
table(nchar(GWAS_atlas_vcf$ALT))
nrow(GWAS_atlas_vcf)
GWAS_atlas_vcf <- GWAS_atlas_vcf[which(nchar(GWAS_atlas_vcf$ALT) == 1),]
unique(GWAS_atlas_vcf$`#CHROM`)
chr_name <- data.frame(original = unique(GWAS_atlas_vcf$`#CHROM`),
                       rename = c(1,10,11,12,2,3,4,5,6,7,8,9),
                       row.names = unique(GWAS_atlas_vcf$`#CHROM`))
GWAS_atlas_vcf$`#CHROM` <- chr_name[GWAS_atlas_vcf$`#CHROM`, "rename"]
GWAS_atlas_vcf <- GWAS_atlas_vcf[, c(1,2,4,5)]
colnames(GWAS_atlas_vcf)[1] <- "CHR"
GWAS_atlas_vcf$SNP_CHR_BP <- paste0(GWAS_atlas_vcf$CHR, "_", GWAS_atlas_vcf$POS)
rownames(GWAS_atlas_vcf) <- GWAS_atlas_vcf$SNP_CHR_BP

GWAS_atlas <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs/GWAS_atlas/gwas_association_result_for_rice.txt",
                        sep = "\t", header = T, check.names = F)
chr_rename <- data.frame(ch = paste0("chr", 1:12),
                         re = 1:12,
                         row.names = paste0("chr", 1:12))
GWAS_atlas$Chr <- chr_rename[GWAS_atlas$Chr, "re"]
table(GWAS_atlas$`Penotype type`)
table(GWAS_atlas$`Ref ver`)
table(GWAS_atlas$Environment)
table(GWAS_atlas$`Species/Population`)
unique(GWAS_atlas$`Species/Population`)
GWAS_atlas$Chr_length <- Chr_length[GWAS_atlas$Chr, "width"]
GWAS_atlas <- GWAS_atlas[-which(GWAS_atlas$Pos > GWAS_atlas$Chr_length),]
GWAS_atlas$`P-value` <- as.numeric(GWAS_atlas$`P-value`)
GWAS_atlas <- GWAS_atlas[!is.na(GWAS_atlas$`P-value`),]
GWAS_atlas <- GWAS_atlas[which(GWAS_atlas$`P-value` < 0.00001),]
GWAS_atlas$SNP_CHR_BP <- paste0(GWAS_atlas$Chr, "_", GWAS_atlas$Pos)
GWAS_atlas <- GWAS_atlas[which(GWAS_atlas$SNP_CHR_BP %in% GWAS_atlas_vcf$SNP_CHR_BP),]

unique(GWAS_atlas$Trait)
unique(GWAS_atlas$Chr)
GWAS_atlas$Trait <- unlist(lapply(GWAS_atlas$Trait, function(x){
  gsub(pattern = " ", replacement = "_", x)
}))

GWAS_atlas_SNP <- data.frame(SNP = paste0("GWAS_atlas_", 1:nrow(GWAS_atlas)),
                             CHR = GWAS_atlas$Chr,
                             BP = GWAS_atlas$Pos,
                             TRAIT = GWAS_atlas$Trait,
                             Negative_Log10P = -log10(GWAS_atlas$`P-value`),
                             Chr_length = Chr_length[GWAS_atlas$Chr, "width"])
GWAS_atlas_SNP$SNP_CHR_BP <- paste0(GWAS_atlas_SNP$CHR, "_", GWAS_atlas_SNP$BP)
GWAS_atlas_SNP$REF <- GWAS_atlas_vcf[GWAS_atlas_SNP$SNP_CHR_BP, "REF"]
GWAS_atlas_SNP$ALT <- GWAS_atlas_vcf[GWAS_atlas_SNP$SNP_CHR_BP, "ALT"]
GWAS_atlas_SNP$SNP <- paste0("GWAS_", GWAS_atlas_SNP$SNP_CHR_BP)
table(GWAS_atlas_SNP$ALT)
table(GWAS_atlas_SNP$REF)

### SNPSeek
files <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs/SNPSeek/",
                    recursive = F, full.names = F, pattern = ".csv$")
trait <- unlist(lapply(files, function(x){
  unlist(strsplit(x, ".csv", fixed = T))[1]
}))
SNPSeek <- data.frame(row.names = files,
                      trait = trait,
                      file = files)
SNPSeek_SNP <- c()
for (i in SNPSeek$file) {
  temp <- read.csv(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs/SNPSeek/", i),
                   sep = ",", header = T)
  temp <- data.frame(# SNP = paste0("GWAS_atlas_", 1:nrow(GWAS_atlas)),
                     CHR = temp$Chr,
                     BP = temp$Position,
                     TRAIT = SNPSeek[i,"trait"],
                     Negative_Log10P = temp$Value,
                     Chr_length = Chr_length[temp$Chr, "width"])
  SNPSeek_SNP <- as.data.frame(rbind(SNPSeek_SNP,
                                     temp))
}
SNPSeek_SNP$SNP_CHR_BP <- paste0(SNPSeek_SNP$CHR, "_", SNPSeek_SNP$BP)
SNPSeek_SNP <- data.frame(SNP = paste0("SNPSeek_", SNPSeek_SNP$SNP_CHR_BP),
                          SNPSeek_SNP)
sum(as.numeric(SNPSeek_SNP$BP) < as.numeric(SNPSeek_SNP$Chr_length)) == nrow(SNPSeek_SNP) # TRUE
SNPSeek_vcf <- data.table::fread("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs/SNPSeek/NB_final_snp.bim.gz",
                                 sep = "\t", header = F)
SNPSeek_vcf <- data.frame(CHR = SNPSeek_vcf$V1,
                          POS = SNPSeek_vcf$V4,
                          REF = SNPSeek_vcf$V6,
                          ALT = SNPSeek_vcf$V5)
SNPSeek_vcf$SNP_CHR_BP <- paste0(SNPSeek_vcf$CHR, "_", SNPSeek_vcf$POS)
rownames(SNPSeek_vcf) <- SNPSeek_vcf$SNP_CHR_BP
nrow(SNPSeek_SNP)
SNPSeek_SNP <- SNPSeek_SNP[which(SNPSeek_SNP$SNP_CHR_BP %in% SNPSeek_vcf$SNP_CHR_BP),]
nrow(SNPSeek_SNP)
SNPSeek_SNP$REF <- SNPSeek_vcf[SNPSeek_SNP$SNP_CHR_BP, "REF"]
SNPSeek_SNP$ALT <- SNPSeek_vcf[SNPSeek_SNP$SNP_CHR_BP, "ALT"]
table(SNPSeek_SNP$REF)
table(SNPSeek_SNP$ALT)

#### RiceVarMap2.0
RiceVarMap_phone <- data.frame(original = c("Heading_date", "Grain_thickness", "Grain_width", "Grain_length", "Spikelet_length",
                                            "Grain_weight", "Yield", "Num_effective_panicles", "Num_panicles", "Plant_height"),
                               after = c("heading_date", "grain_thickness", "grain_width", "grain_length", "spikelet_length",
                                         "grain_weight", "yield", "effective_panicles_number", "panicles_number", "plant_height"),
                               row.names = c("Heading_date", "Grain_thickness", "Grain_width", "Grain_length", "Spikelet_length",
                                             "Grain_weight", "Yield", "Num_effective_panicles", "Num_panicles", "Plant_height"))
RiceVarMap2.0_SNP <- c()
for (i in RiceVarMap_phone$original) {
  temp <- data.table::fread(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/",i,"/", "plink_result.",i,".assoc.linear"),
                            sep = " ", header = T)
  if (length(unique(temp$TEST)) != 3) {
    message(paste0(i, " 出错了"))
  }
  temp <- temp[which(temp$P < 0.00001),]
  temp <- temp[which(temp$TEST == "ADD"),]
  temp <- data.frame(CHR = temp$CHR,
                     BP = temp$BP,
                     TRAIT = RiceVarMap_phone[i, "after"],
                     Negative_Log10P = -log10(temp$P),
                     Chr_length = Chr_length[temp$CHR, "width"])
  RiceVarMap2.0_SNP <- as.data.frame(rbind(RiceVarMap2.0_SNP,
                                           temp))
}
RiceVarMap2.0_SNP$SNP_CHR_BP <- paste0(RiceVarMap2.0_SNP$CHR, "_", RiceVarMap2.0_SNP$BP)
RiceVarMap2.0_SNP <- data.frame(SNP = paste0("RiceVarMap2.0_", RiceVarMap2.0_SNP$SNP_CHR_BP),
                                RiceVarMap2.0_SNP)
sum(RiceVarMap2.0_SNP$BP < RiceVarMap2.0_SNP$Chr_length) == nrow(RiceVarMap2.0_SNP)

RiceVarMap2.0_vcf <- data.table::fread("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/rice4k_geno_no_del.bim",
                                       sep = "\t", header = F)
RiceVarMap2.0_vcf <- data.frame(CHR = RiceVarMap2.0_vcf$V1,
                                POS = RiceVarMap2.0_vcf$V4,
                                REF = RiceVarMap2.0_vcf$V6,
                                ALT = RiceVarMap2.0_vcf$V5)
RiceVarMap2.0_vcf$SNP_CHR_BP <- paste0(RiceVarMap2.0_vcf$CHR, "_", RiceVarMap2.0_vcf$POS)
rownames(RiceVarMap2.0_vcf) <- RiceVarMap2.0_vcf$SNP_CHR_BP
RiceVarMap2.0_vcf <- RiceVarMap2.0_vcf[which(RiceVarMap2.0_vcf$REF %in% c("A","T","C","G")),]
RiceVarMap2.0_vcf <- RiceVarMap2.0_vcf[which(RiceVarMap2.0_vcf$ALT %in% c("A","T","C","G")),]
table(RiceVarMap2.0_vcf$REF)
table(RiceVarMap2.0_vcf$ALT)
RiceVarMap2.0_SNP$REF <- RiceVarMap2.0_vcf[RiceVarMap2.0_SNP$SNP_CHR_BP, "REF"]
RiceVarMap2.0_SNP$ALT <- RiceVarMap2.0_vcf[RiceVarMap2.0_SNP$SNP_CHR_BP, "ALT"]
table(RiceVarMap2.0_SNP$REF)
table(RiceVarMap2.0_SNP$ALT)

###
GWAS_atlas_SNP$Database <- "GWAS_atlas"
SNPSeek_SNP$Database <- "SNPSeek"
RiceVarMap2.0_SNP$Database <- "RiceVarMap2.0"

All_SNP <- as.data.frame(rbind(GWAS_atlas_SNP,
                               SNPSeek_SNP,
                               RiceVarMap2.0_SNP))

# 重新整理Trait
# traits <- sort(unique(All_SNP$TRAIT))
setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/SNPs")
# traits <- data.frame(Original_label = traits)
# openxlsx::write.xlsx(traits, "traits_2.xlsx")
traits <- openxlsx::read.xlsx("traits_2.xlsx")
length(unique(traits$Changed_label)) # 348
rownames(traits) <- traits$Original_label
sum(All_SNP$TRAIT %in% traits$Original_label) == nrow(All_SNP)
All_SNP$TRAIT_Renamed <- traits[All_SNP$TRAIT, "Changed_label"]
sum(is.na(All_SNP$TRAIT_Renamed))
All_SNP$SNP_name_database <- paste0(All_SNP$CHR,
                                    "_",
                                    All_SNP$BP,
                                    "_",
                                    All_SNP$TRAIT_Renamed,
                                    "_",
                                    All_SNP$Database)
sum(duplicated(All_SNP$SNP_name_database)) # 25249 association duplicated
All_SNP <- All_SNP[!duplicated(All_SNP$SNP_name_database),] # 55164
# All_SNP$SNP_CHR_BP <- paste0(All_SNP$CHR, "_", All_SNP$BP)
# All_SNP[which(All_SNP$BP > All_SNP$Chr_length),]
sum(duplicated(All_SNP$SNP_CHR_BP)) # 8993
unique(All_SNP$Database)
length(unique(All_SNP$SNP_CHR_BP))
# 55164 associations, 46171 SNPs
length(unique(All_SNP$TRAIT_Renamed))
# 312 traits

All_SNP_position <- data.frame(chr = All_SNP$CHR,
                               start = All_SNP$BP,
                               end = All_SNP$BP+1)
All_SNP_position$Database <- All_SNP$Database
All_SNP_position$TRAIT_Renamed <- All_SNP$TRAIT_Renamed
All_SNP_position$SNP_CHR_BP <- All_SNP$SNP_CHR_BP
All_SNP_position$Negative_Log10P <- All_SNP$Negative_Log10P
All_SNP_position$Chr_length <- All_SNP$Chr_length
All_SNP_position$REF <- All_SNP$REF
All_SNP_position$ALT <- All_SNP$ALT

All_SNP_position_filtered <- All_SNP_position[which(All_SNP_position$Negative_Log10P > 5),]

## 暂时不过滤掉主等位与参考基因组不一致的SNP
if (FALSE) {
  fasta_sequences <- readDNAStringSet("./Ref/Oryza_sativa.IRGSP-1.0.fa")
  names(fasta_sequences)
  names(fasta_sequences) <- c("1", "Mt", "Pt", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10", "11", "12")
  SNP_range <- data.frame(chr = All_SNP_position_filtered$chr,
                          start = All_SNP_position_filtered$start - 24,
                          end = All_SNP_position_filtered$end + 25,
                          strand = "*")
  sub_sequence <- c()
  for (i in 1:nrow(SNP_range)) {
    print(i)
    temp_seq <- subseq(fasta_sequences[[as.character(SNP_range[i,1])]], start = as.numeric(SNP_range[i,2]), end = as.numeric(SNP_range[i,3]))
    sub_sequence <- c(sub_sequence, list(temp_seq))
  }
  sub_sequence <- unlist(lapply(sub_sequence, as.character))
  sub_sequence_ref <- sub_sequence
  ref_sum <- 0
  not_ref <- c()
  for (i in 1:length(sub_sequence_ref)) {
    if (substr(sub_sequence_ref[i], 25, 25) == All_SNP_position_filtered[i,"REF"]) {
      ref_sum <- ref_sum + 1
    } else {
      not_ref <- c(not_ref, i)
    }
  }
  not_ref_df <- All_SNP_position_filtered_unique[not_ref,]
}

All_SNP_position_filtered_ref <- All_SNP_position_filtered
All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$SNP_CHR_BP == "1_31283533"),]

GWAS_atlas <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$Database == "GWAS_atlas"), "SNP_CHR_BP"]
GWAS_atlas <- unique(GWAS_atlas) # 24690
SNPSeek <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$Database == "SNPSeek"), "SNP_CHR_BP"]
SNPSeek <- unique(SNPSeek) # 21264
RiceVarMap2.0 <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$Database == "RiceVarMap2.0"), "SNP_CHR_BP"]
RiceVarMap2.0 <- unique(RiceVarMap2.0) # 645 (437)

library(gplots)
library(VennDiagram)
setwd("../../")
pdf("Figure7_3SNP_database_venn.pdf", width = 5.2, height = 5)
Tvenn <- venn.diagram(list("GWAS Atlas" =GWAS_atlas,
                           "SNPSeek" = SNPSeek,
                           "RiceVarMap V2.0" = RiceVarMap2.0),
                      filename = NULL,
                      lwd = 1, lty = 2, alpha = 0.8, cex = 1,
                      col = c("#1998D3", "#41A592", "#D05146"),
                      fill = c("#1998D3", "#41A592", "#D05146"),
                      cat.col = c("#1998D3", "#41A592", "#D05146"),
                      reverse = TRUE)
grid.draw(Tvenn)
dev.off()

All_SNP_position_filtered_ref$SNP_name <- paste0(All_SNP_position_filtered_ref$chr,
                                             "_",
                                             All_SNP_position_filtered_ref$start,
                                             "_",
                                             All_SNP_position_filtered_ref$TRAIT_Renamed)
All_SNP_position_filtered_ref$SNP_name_database <- paste0(All_SNP_position_filtered_ref$SNP_name,
                                                      "_",
                                                      All_SNP_position_filtered_ref$Database)
sum(duplicated(All_SNP_position_filtered_ref$SNP_name)) # 94
sum(duplicated(All_SNP_position_filtered_ref$SNP_name_database)) # 0

Unique_SNP <- All_SNP_position_filtered_ref[!duplicated(All_SNP_position_filtered_ref$SNP_CHR_BP),]
Unique_SNP_GR <- makeGRangesFromDataFrame(Unique_SNP,
                                          keep.extra.columns = TRUE,
                                          seqnames.field = "chr",
                                          start.field = "start",
                                          end.field = "end",
                                          ignore.strand = TRUE)
Unique_SNP_GR_nearstgene <- .fastAnnoPeaks(peaks = Unique_SNP_GR,
                                           BSgenome = BSgenome.OSativa.NCBI.IRGSPv1.0,
                                           geneAnnotation = gene_annotation,
                                           promoterRegion = c(2000, 100))
Unique_SNP_GR_nearstgene_DF <- as.data.frame(Unique_SNP_GR_nearstgene)
rownames(Unique_SNP_GR_nearstgene_DF) <- Unique_SNP_GR_nearstgene_DF$SNP_CHR_BP
colnames(Unique_SNP_GR_nearstgene_DF)
All_SNP_position_filtered_ref$distToGeneStart <- Unique_SNP_GR_nearstgene_DF[All_SNP_position_filtered_ref$SNP_CHR_BP, "distToGeneStart"]
All_SNP_position_filtered_ref$nearestGene <- Unique_SNP_GR_nearstgene_DF[All_SNP_position_filtered_ref$SNP_CHR_BP, "nearestGene"]
All_SNP_position_filtered_ref$peakType <- Unique_SNP_GR_nearstgene_DF[All_SNP_position_filtered_ref$SNP_CHR_BP, "peakType"]
All_SNP_position_filtered_ref$nearestTSS <- Unique_SNP_GR_nearstgene_DF[All_SNP_position_filtered_ref$SNP_CHR_BP, "nearestTSS"]
All_SNP_position_filtered_ref$distToTSS <- Unique_SNP_GR_nearstgene_DF[All_SNP_position_filtered_ref$SNP_CHR_BP, "distToTSS"]
table(All_SNP_position_filtered_ref$peakType)
All_SNP_position_filtered_ref$SNP_region <- All_SNP_position_filtered_ref$peakType
All_SNP_position_filtered_ref$SNP_region[which(All_SNP_position_filtered_ref$SNP_region %in% c("Exonic", "Intronic"))] <- "Intragenic"
All_SNP_position_filtered_ref$SNP_region[which(All_SNP_position_filtered_ref$SNP_region %in% c("Promoter", "Distal"))] <- "Intergenic"
table(All_SNP_position_filtered_ref$SNP_region)

library(ggplot2)
library(ggforce)
Intergenic_temp <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$SNP_region == "Intergenic"),]
Intergenic_temp <- Intergenic_temp[!duplicated(Intergenic_temp$SNP_CHR_BP),]
Intragenic_temp <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$SNP_region == "Intragenic"),]
Intragenic_temp <- Intragenic_temp[!duplicated(Intragenic_temp$SNP_CHR_BP),]
temp <- as.data.frame(rbind(Intergenic_temp,
                            Intragenic_temp))
temp <- as.data.frame(table(temp$SNP_region))
colnames(temp) <- c("PeakType", "Count")
temp$Ratio <- temp$Count / sum(temp$Count)
temp$Label <- c("Intergenic(n = 37,649)", "Intragenic(n = 8,522)")
pdf("Figure7_PeakType_Percentage.pdf", width = 5, height = 5)
pie(temp$Count, labels = temp$Label,
    radius = 1.0, clockwise = T,
    col = c("#E59CC4", "#57C3F3"),
    main = "PeakType Percentage")
dev.off()

GWAS_atlas <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$Database == "GWAS_atlas"), ]
GWAS_atlas <- GWAS_atlas[!duplicated(GWAS_atlas$SNP_CHR_BP),]
SNPSeek <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$Database == "SNPSeek"),]
SNPSeek <- SNPSeek[!duplicated(SNPSeek$SNP_CHR_BP),]
RiceVarMap2.0 <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$Database == "RiceVarMap2.0"),]
RiceVarMap2.0 <- RiceVarMap2.0[!duplicated(RiceVarMap2.0$SNP_CHR_BP),]
temp <- as.data.frame(rbind(GWAS_atlas, SNPSeek, RiceVarMap2.0))
pdf("Figure7_Database_SNP_region_percent.pdf", height = 5.7, width = 5)
ggplot(data = temp,
       aes(x = Database, fill = SNP_region)) +
  geom_bar(stat = "count", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#E59CC4", "#57C3F3")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black")) +
  labs(x = "", y = "Percentage")
dev.off()


### 
gene_id_map <- readRDS("gene_id_map.rds")
rownames(gene_id_map) <- gene_id_map$symbol

######## Intergenic
Intergenic <- All_SNP_position_filtered_ref[which(All_SNP_position_filtered_ref$SNP_region == "Intergenic"),]
length(unique(Intergenic$TRAIT_Renamed)) # 299

sum(duplicated(Intergenic$SNP_CHR_BP)) # 7491
Intergenic_unique_SNP <- Intergenic[!duplicated(Intergenic$SNP_CHR_BP),] # 提取SNP，做peak注释
nrow(Intergenic_unique_SNP) # 37649 SNPs

# SNP annotated to peaks
Chr_length_vector <- Chr_length$width
names(Chr_length_vector) <- Chr_length$seqnames
Intergenic_GR <- makeGRangesFromDataFrame(Intergenic_unique_SNP, keep.extra.columns = TRUE,
                                          seqnames.field = "chr",
                                          start.field = "start",
                                          end.field = "end",
                                          ignore.strand = TRUE,
                                          seqinfo = Chr_length_vector)
proj_pass_filter <- readRDS("3.proj_pass_filter_withPeaks.rds")
peaksets <- getPeakSet(proj_pass_filter)
peaksets@ranges@NAMES <- NULL
peaksets_df <- as.data.frame(peaksets)
peaksets_df$peak_id <- paste0(peaksets_df$seqnames, "_", peaksets_df$start, "_", peaksets_df$end)
rownames(peaksets_df) <- peaksets_df$peak_id

overlapRegions <- findOverlaps(peaksets, Intergenic_GR, ignore.strand = TRUE, minoverlap = 2)
# 8454 intergenic SNPs overlaped with peaks

Intergenic_GR_DF <- as.data.frame(Intergenic_GR)

overlapRegions_DF <- data.frame(Intergenic_SNP = paste0(Intergenic_GR_DF$seqnames, "_", Intergenic_GR_DF$start)[subjectHits(overlapRegions)],
                                Peakset = peaksets_df$peak_id[queryHits(overlapRegions)])
rownames(overlapRegions_DF) <- overlapRegions_DF$Intergenic_SNP
Intergenic_overlap_with_peaks <- Intergenic[which(Intergenic$SNP_CHR_BP %in% overlapRegions_DF$Intergenic_SNP),]
nrow(Intergenic_overlap_with_peaks)
# 10254 associations
length(unique(Intergenic_overlap_with_peaks$SNP_CHR_BP)) # 8454
length(unique(Intergenic_overlap_with_peaks$TRAIT_Renamed)) # 219

Intergenic$OnPeak <- ifelse(Intergenic$SNP_CHR_BP %in% Intergenic_overlap_with_peaks$SNP_CHR_BP,
                            TRUE, FALSE)
table(Intergenic$OnPeak)
temp <- Intergenic[!duplicated(Intergenic$SNP_CHR_BP),]
table(temp$OnPeak)
temp <- as.data.frame(table(temp$OnPeak))
colnames(temp) <- c("SNP_region", "Count")
temp$Ratio <- temp$Count / sum(temp$Count)
temp$Label <- c("Not located in peak region(n = 29,195)", "Located in peak region(n = 8,454)")
pdf("Figure7_37649_SNP_peak_percent.pdf", width = 5, height = 5)
pie(temp$Count, labels = temp$Label,
    radius = 1.0, clockwise = T,
    col = c("#E6E6FA", "#6495ED"),
    main = "")
dev.off()

Intergenic_overlap_with_peaks$peak_id <- overlapRegions_DF[Intergenic_overlap_with_peaks$SNP_CHR_BP, "Peakset"]
temp <- as.data.frame.array(table(Intergenic_overlap_with_peaks$TRAIT_Renamed,
                                  Intergenic_overlap_with_peaks$peak_id))
temp <- as.matrix(temp)
temp[which(temp > 1)] <- 1
trait_peaks <- rowSums(temp)
temp <- as.data.frame(trait_peaks)
temp$index <- ""
temp <- temp[order(as.numeric(temp$trait_peaks), decreasing = F),]
temp$index <- 1:nrow(temp)
temp$trait_peaks <- as.numeric(temp$trait_peaks)
pdf("Figure7_trait_peak_5.pdf", width = 6, height = 5)
ggplot(data = temp, aes(x = index, y = trait_peaks, color = trait_peaks)) +
  geom_point(size = 1, alpha = 0.6) +
  theme_bw() +
  geom_hline(yintercept = 5, color = "red") +
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black")) +
  labs(x = "Trait index (Ordered by peak number)", y = "Peaks number", color = "Count") +
  scale_color_viridis(option = "H")
dev.off()

trait_peaks <- trait_peaks[which(trait_peaks >= 5)] # 89 trait passed QC
Intergenic_overlap_with_peaks_filtered <- Intergenic_overlap_with_peaks[which(Intergenic_overlap_with_peaks$TRAIT_Renamed %in% names(trait_peaks)),]
length(unique(Intergenic_overlap_with_peaks_filtered$SNP_CHR_BP)) # 8240
nrow(Intergenic_overlap_with_peaks_filtered) # 10003

length(unique(Intergenic_overlap_with_peaks_filtered$peak_id))

Intergenic_overlap_with_peaks_filtered_unique <- Intergenic_overlap_with_peaks_filtered[!duplicated(Intergenic_overlap_with_peaks_filtered$SNP_CHR_BP),]
nrow(Intergenic_overlap_with_peaks_filtered_unique)
Celltype_marker_peaks <- c()
for (i in 1:length(markerPeaksList_Celltype)) {
  temp <- markerPeaksList_Celltype[[i]]
  temp <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
  Celltype_marker_peaks <- c(Celltype_marker_peaks, temp)
}
Tissue_marker_peaks <- c()
for (i in 1:length(markerPeaksList_tissue)) {
  temp <- markerPeaksList_tissue[[i]]
  temp <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
  Tissue_marker_peaks <- c(Tissue_marker_peaks, temp)
}

Intergenic_overlap_with_peaks_filtered_unique$Celltype_Marker_Peaks <- Intergenic_overlap_with_peaks_filtered_unique$peak_id %in% Celltype_marker_peaks
table(Intergenic_overlap_with_peaks_filtered_unique$Celltype_Marker_Peaks)

Intergenic_overlap_with_peaks_filtered_unique$Tissue_Marker_Peaks <- Intergenic_overlap_with_peaks_filtered_unique$peak_id %in% Tissue_marker_peaks
table(Intergenic_overlap_with_peaks_filtered_unique$Tissue_Marker_Peaks)

table(Intergenic_overlap_with_peaks_filtered_unique$Celltype_Marker_Peaks,
      Intergenic_overlap_with_peaks_filtered_unique$Tissue_Marker_Peaks)

nrow(Intergenic_overlap_with_peaks_filtered_unique) # 8240

temp <- Intergenic_overlap_with_peaks_filtered_unique[!duplicated(Intergenic_overlap_with_peaks_filtered_unique$peak_id),]
table(temp$Celltype_Marker_Peaks,
      temp$Tissue_Marker_Peaks)
table(temp$Tissue_Marker_Peaks)
temp <- as.data.frame(table(temp$Tissue_Marker_Peaks))
colnames(temp) <- c("PeakType", "Count")
temp$Ratio <- temp$Count / sum(temp$Count)
temp$Label <- c("Not tissue marker peaks(n = 3,533)", "Tissue marker peaks(n = 1,031)")
pdf("Figure7_4564_peaks_tissue_marker.pdf", width = 5, height = 5)
pie(temp$Count, labels = temp$Label,
    radius = 1.0, clockwise = T,
    col = c("#FAF0E6", "#FF6347"),
    main = "")
dev.off()

temp <- Intergenic_overlap_with_peaks_filtered_unique[!duplicated(Intergenic_overlap_with_peaks_filtered_unique$peak_id),]
temp <- temp[which(temp$Tissue_Marker_Peaks == FALSE),]
table(temp$Celltype_Marker_Peaks)
temp <- as.data.frame(table(temp$Celltype_Marker_Peaks))
colnames(temp) <- c("PeakType", "Count")
temp$Ratio <- temp$Count / sum(temp$Count)
temp$Label <- c("Not cell-type marker peaks(n = 2,407)", "Cell-type marker peaks(n = 1,126)")
pdf("Figure7_3533_peaks_celltype_marker.pdf", width = 5, height = 5)
pie(temp$Count, labels = temp$Label,
    radius = 1.0, clockwise = T,
    col = c("#FAF0E6", "#FF6347"),
    main = "")
dev.off()

gene_trait <- Intergenic_overlap_with_peaks_filtered[,c("nearestGene", "TRAIT_Renamed")]
gene_trait <- gene_trait[!duplicated(gene_trait),]
gene_trait_num <- as.data.frame(sort(table(gene_trait$nearestGene), decreasing = T))

Intergenic_trait_peaks <- split(Intergenic_overlap_with_peaks_filtered, Intergenic_overlap_with_peaks_filtered$TRAIT_Renamed)
length(unique(Intergenic_overlap_with_peaks_filtered$SNP_CHR_BP))
# 10003 association, 8240 SNPs
# celltype enrichment
{
  markersPeaks_Celltype <- getMarkerFeatures(
    ArchRProj = proj_pass_filter, 
    useMatrix = "PeakMatrix", 
    groupBy = "Celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerPeaksList_Celltype <- getMarkers(markersPeaks_Celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  markerPeaksList_Celltype <- markerPeaksList_Celltype@listData
  celltype_enrichment_matrix <- matrix(1, nrow = length(markerPeaksList_Celltype),
                                       ncol = length(names(Intergenic_trait_peaks)))
  rownames(celltype_enrichment_matrix) <- names(markerPeaksList_Celltype)
  colnames(celltype_enrichment_matrix) <- names(Intergenic_trait_peaks)
  celltype_enrichment_matrix_binary <-  matrix(0, nrow = length(markerPeaksList_Celltype),
                                               ncol = length(names(Intergenic_trait_peaks)))
  rownames(celltype_enrichment_matrix_binary) <- names(markerPeaksList_Celltype)
  colnames(celltype_enrichment_matrix_binary) <- names(Intergenic_trait_peaks)
  for (cell in rownames(celltype_enrichment_matrix)) {
    print(cell)
    for (trait in colnames(celltype_enrichment_matrix)) {
      cell_markers <- markerPeaksList_Celltype[[cell]]
      cell_markers <- as.data.frame(cell_markers)
      cell_markers <- paste0(cell_markers$seqnames, "_", cell_markers$start, "_", cell_markers$end)
      trait_peaks <- Intergenic_trait_peaks[[trait]]
      trait_peaks <- unique(trait_peaks$peak_id)
      common <- intersect(cell_markers, trait_peaks)
      p <- phyper((length(common) - 1), length(trait_peaks),
                  (132008-length(trait_peaks)),
                  length(cell_markers), lower.tail = F)
      celltype_enrichment_matrix[cell, trait] <- p
      if (p < 0.05) {
        celltype_enrichment_matrix_binary[cell, trait] <- 1
      }
    }
  }
  celltype_enrichment_matrix2 <- celltype_enrichment_matrix
  trait_sig <- colSums(celltype_enrichment_matrix_binary)
  trait_sig <- trait_sig[which(trait_sig > 0)]
  celltype_enrichment_matrix <- celltype_enrichment_matrix[,which(colnames(celltype_enrichment_matrix) %in% names(trait_sig))]
  celltype_sig <- rowSums(celltype_enrichment_matrix_binary)
  celltype_sig <- celltype_sig[which(celltype_sig > 0)]
  celltype_enrichment_matrix <- celltype_enrichment_matrix[which(rownames(celltype_enrichment_matrix) %in% names(celltype_sig)),]
  trait_order <- sort(names(trait_sig))
  trait_order <- data.frame(trait = trait_order)
  # Physiological
  # Morphological
  # 
  trait_order$trait_group <- c("Physiological traits", "Physiological traits", "Physiological traits", "Physiological traits", "Morphological traits",
                               "Morphological traits", "Morphological traits", "Morphological traits", "Physiological traits", "Physiological traits",
                               "Physiological traits", "Morphological traits", "Physiological traits", "Physiological traits", "Yield-related traits",
                               "Yield-related traits", "Yield-related traits", "Yield-related traits", "Yield-related traits", "Yield-related traits",
                               "Yield-related traits", "Morphological traits", "Morphological traits", "Physiological traits", "Yield-related traits",
                               "Morphological traits", "Morphological traits", "Yield-related traits", "Physiological traits", "Stress-response traits",
                               "Ecological traits", "Stress-response traits", "Physiological traits", "Yield-related traits", "Morphological traits",
                               "Morphological traits", "Morphological traits", "Yield-related traits"
                               )
  trait_order$trait_group <- factor(trait_order$trait_group,
                                    levels = c("Yield-related traits", "Stress-response traits",
                                               "Physiological traits", "Morphological traits",
                                               "Ecological traits"))
  trait_order <- trait_order[order(trait_order$trait_group, trait_order$trait),]
  celltype_order <- c("Branch meristems (BM)", "Cryptic bract/bract (cb/b)", "Inflorescence meristem (IM)", "Rachis", "Spikelet meristem (SM)",
                      "Mestome sheath", "Procambium", "Large parenchyma (PO)", "Fiber", "Large parenchyma (MO)", "Mesophyll (MO)",
                      "Mesophyll precursor", "Mesophyll (PO)", "Mesophyll initial", "Phloem",
                      "Cortex", "Endodermis", "Pericycle", "Root hair"
                      )
  
  temp <- data.frame(Celltype = proj_pass_filter$Celltype,
                     Tissue = proj_pass_filter$Tissues)
  temp <- temp[which(temp$Celltype %in% celltype_order),]
  temp$Celltype <- factor(temp$Celltype, levels = celltype_order)
  pdf("Figure7_Celltype_enrich_x_tissue_percent.pdf", width = 7, height = 3)
  ggplot(data = temp, aes(x = Celltype, fill = Tissue)) +
    geom_bar(stat = "count", position = "fill") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
    scale_y_continuous(labels = scales::percent)
  dev.off()
  
  celltype_enrichment_matrix <- celltype_enrichment_matrix[celltype_order, trait_order$trait]
  if (!is.null(celltype_enrichment_matrix)) {
    pmt <- celltype_enrichment_matrix
    pmt <- as.data.frame(pmt)
    ssmt <- celltype_enrichment_matrix <= 0.0001
    pmt[ssmt] <- "****"
    smt <- celltype_enrichment_matrix > 0.0001 & celltype_enrichment_matrix <= 0.001
    pmt[smt] <- "***"
    smt2 <- celltype_enrichment_matrix > 0.001 & celltype_enrichment_matrix <= 0.01
    pmt[smt2] <- "**"
    smt3 <- celltype_enrichment_matrix > 0.01 & celltype_enrichment_matrix <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  
  celltype_enrichment_matrix_log <- -log10(celltype_enrichment_matrix)
  celltype_enrichment_matrix_log[which(celltype_enrichment_matrix_log == -log10(1))] <- NA
  pdf("Figure7_celltype_enrichment_matrix_peak.pdf", height = 10, width = 12)
  pheatmap(t(celltype_enrichment_matrix_log),
           color = colorRampPalette(c("yellow", "red"))(20),
           scale = "none", cluster_row = F,
           cluster_col = F,
           treeheight_col = 0,
           display_numbers = t(pmt),
           fontsize_number = 14,
           number_color = "black",
           cellwidth = 15, cellheight = 10)
  dev.off()
  
  if (!is.null(celltype_enrichment_matrix)) {
    pmt <- celltype_enrichment_matrix
    pmt <- as.data.frame(pmt)
    ssmt <- celltype_enrichment_matrix <= 0.0001
    pmt[ssmt] <- 1
    smt <- celltype_enrichment_matrix > 0.0001 & celltype_enrichment_matrix <= 0.001
    pmt[smt] <- 1
    smt2 <- celltype_enrichment_matrix > 0.001 & celltype_enrichment_matrix <= 0.01
    pmt[smt2] <- 1
    smt3 <- celltype_enrichment_matrix > 0.01 & celltype_enrichment_matrix <= 0.05
    pmt[smt3] <- 1
    pmt[!ssmt&!smt&!smt2&!smt3] <- 0
  }
  
  temp <- data.frame(rowSums(pmt))
  temp$Celltype <- rownames(temp)
  colnames(temp)[1] <- "Count"
  temp <- temp[order(temp$Count, decreasing = T),]
  temp$Celltype <- factor(temp$Celltype, levels = temp$Celltype)
  
  pdf("Figure7.celltype_enriched_trait_count_peak.pdf", height = 4, width = 6)
  ggplot(data = temp, aes(x = Celltype, y = Count)) +
    geom_bar(stat = "identity", width = 0.6, fill = "#476D87") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "", y = "Enriched trait count")
  dev.off()
  peak_celltype_traits <- names(trait_sig)
}
# tissue enrichment
{
  markersPeaks_tissue <- getMarkerFeatures(
    ArchRProj = proj_pass_filter, 
    useMatrix = "PeakMatrix", 
    groupBy = "Tissues",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerPeaksList_tissue <- getMarkers(markersPeaks_tissue, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  markerPeaksList_tissue <- markerPeaksList_tissue@listData
  tissue_enrichment_matrix <- matrix(1, nrow = length(markerPeaksList_tissue),
                                     ncol = length(names(Intergenic_trait_peaks)))
  rownames(tissue_enrichment_matrix) <- names(markerPeaksList_tissue)
  colnames(tissue_enrichment_matrix) <- names(Intergenic_trait_peaks)
  tissue_enrichment_matrix_binary <- matrix(0, nrow = length(markerPeaksList_tissue),
                                     ncol = length(names(Intergenic_trait_peaks)))
  rownames(tissue_enrichment_matrix_binary) <- names(markerPeaksList_tissue)
  colnames(tissue_enrichment_matrix_binary) <- names(Intergenic_trait_peaks)
  
  for (tissue in rownames(tissue_enrichment_matrix)) {
    print(tissue)
    for (trait in colnames(tissue_enrichment_matrix)) {
      tissue_markers <- markerPeaksList_tissue[[tissue]]
      tissue_markers <- as.data.frame(tissue_markers)
      tissue_markers <- paste0(tissue_markers$seqnames, "_", tissue_markers$start, "_", tissue_markers$end)
      trait_peaks <- Intergenic_trait_peaks[[trait]]
      trait_peaks <- unique(trait_peaks$peak_id)
      common <- intersect(tissue_markers, trait_peaks)
      p <- phyper((length(common) - 1), length(trait_peaks),
                  (132008-length(trait_peaks)),
                  length(tissue_markers), lower.tail = F)
      tissue_enrichment_matrix[tissue, trait] <- p
      if (p < 0.05) {
        tissue_enrichment_matrix_binary[tissue, trait] <- 1
      }
    }
  }
  tissue_enrichment_matrix2 <- tissue_enrichment_matrix
  
  trait_sig <- colSums(tissue_enrichment_matrix_binary)
  trait_sig <- trait_sig[which(trait_sig > 0)]
  trait_order <- sort(names(trait_sig))
  trait_order <- data.frame(trait = trait_order)
  trait_order$trait_group <- c("Morphological traits", "Morphological traits", "Morphological traits", "Morphological traits",
                               "Physiological traits", "Yield-related traits", "Yield-related traits", "Yield-related traits",
                               "Morphological traits", "Morphological traits", "Morphological traits", "Physiological traits",
                               "Stress-response traits", "Morphological traits"
                               )
  trait_order$trait_group <- factor(trait_order$trait_group,
                                    levels = c("Yield-related traits", "Stress-response traits",
                                               "Physiological traits", "Morphological traits",
                                               "Ecological traits"))
  trait_order <- trait_order[order(trait_order$trait_group, trait_order$trait),]
  tissue_enrichment_matrix <- tissue_enrichment_matrix[c("IM", "Shoot", "Leaf", "Root"), trait_order$trait]
  
  if (!is.null(tissue_enrichment_matrix)) {
    pmt <- tissue_enrichment_matrix
    pmt <- as.data.frame(pmt)
    ssmt <- tissue_enrichment_matrix <= 0.0001
    pmt[ssmt] <- "****"
    smt <- tissue_enrichment_matrix > 0.0001 & tissue_enrichment_matrix <= 0.001
    pmt[smt] <- "***"
    smt2 <- tissue_enrichment_matrix > 0.001 & tissue_enrichment_matrix <= 0.01
    pmt[smt2] <- "**"
    smt3 <- tissue_enrichment_matrix > 0.01 & tissue_enrichment_matrix <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  # 
  # trait_order2 <- sort(names(trait_sig))
  # trait_order2 <- c(trait_order2[c(6:8,1:5,9:14)])
  # 
  tissue_enrichment_matrix_log <- -log10(tissue_enrichment_matrix)
  tissue_enrichment_matrix_log[which(tissue_enrichment_matrix_log == -log10(1))] <- NA
  # tissue_enrichment_matrix_log <- tissue_enrichment_matrix_log[,trait_order2]
  # pmt <- pmt[,trait_order2]
  pdf("Figure7_tissue_enrichment_matrix_peak.pdf", height = 6, width = 6)
  pheatmap(t(tissue_enrichment_matrix_log),
           color = colorRampPalette(c("yellow", "red"))(20),
           scale = "none", cluster_row = F,
           cluster_col = F,      
           treeheight_col = 0,
           display_numbers = t(pmt),
           fontsize_number = 12,
           number_color = "black",
           cellwidth = 15, cellheight = 10)
  dev.off()
  
  
  if (!is.null(tissue_enrichment_matrix)) {
    pmt <- tissue_enrichment_matrix
    pmt <- as.data.frame(pmt)
    ssmt <- tissue_enrichment_matrix <= 0.0001
    pmt[ssmt] <- 1
    smt <- tissue_enrichment_matrix > 0.0001 & tissue_enrichment_matrix <= 0.001
    pmt[smt] <- 1
    smt2 <- tissue_enrichment_matrix > 0.001 & tissue_enrichment_matrix <= 0.01
    pmt[smt2] <- 1
    smt3 <- tissue_enrichment_matrix > 0.01 & tissue_enrichment_matrix <= 0.05
    pmt[smt3] <- 1
    pmt[!ssmt&!smt&!smt2&!smt3] <- 0
  }
  
  temp <- data.frame(rowSums(pmt))
  temp$tissue <- rownames(temp)
  colnames(temp)[1] <- "Count"
  temp <- temp[order(temp$Count, decreasing = T),]
  temp$tissue <- factor(temp$tissue, levels = temp$tissue)
  
  pdf("Figure7.tissue_enriched_trait_count_peak.pdf", height = 4, width = 1.85)
  ggplot(data = temp, aes(x = tissue, y = Count)) +
    geom_bar(stat = "identity", width = 0.6, fill = "#476D87") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "", y = "Enriched trait count")
  dev.off()
  
  peak_tissue_traits <- names(trait_sig)
}


###### 修改后的汉明距离
{
  saveRDS(markerPeaksList_Celltype, "./random_peaks/markerPeaksList_Celltype.rds")
  saveRDS(peaksets_df, "./random_peaks/peaksets_df.rds")
  saveRDS(markerPeaksList_tissue, "./random_peaks/markerPeaksList_tissue.rds")
  
  markerPeaksList_Celltype
  
  Celltypes <- names(markerPeaksList_Celltype)
  for (i in Celltypes) {
    command <- paste0("/usr/bin/Rscript ", "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/random_peaks/Celltypes_peaks_random.R ", '"', i,'"', "&")
    system(command)
  }
  Celltype_marker_peaks_random_list <- c()
  for (i in Celltypes) {
    print(i)
    temp <- readRDS(paste0("./random_peaks/", gsub(pattern = "/", replacement = "_", fixed = T, x = i), ".rds"))
    temp <- list(temp)
    names(temp) <- i
    Celltype_marker_peaks_random_list <- c(Celltype_marker_peaks_random_list, temp)
    rm(temp)
    Sys.sleep(15)
  }
  
  Tissues <- names(markerPeaksList_tissue)
  Tissue_marker_peaks_random_list <- c()
  for (i in Tissues) {
    print(i)
    temp <- readRDS(paste0("./random_peaks/", gsub(pattern = "/", replacement = "_", fixed = T, x = i), ".rds"))
    temp <- list(temp)
    names(temp) <- i
    Tissue_marker_peaks_random_list <- c(Tissue_marker_peaks_random_list, temp)
    rm(temp)
    Sys.sleep(15)
  }
  
  Intergenic_trait_peaks
  length(Intergenic_trait_peaks)
  Trait_peak_list <- c()
  for (i in names(Intergenic_trait_peaks)) {
    temp <- Intergenic_trait_peaks[[i]]
    temp <- temp$peak_id
    temp_list <- rep(0, 132008)
    names(temp_list) <- peaksets_df$peak_id
    temp_list[temp] <- 1
    names(temp_list) <- NULL
    Trait_peak_list <- c(Trait_peak_list, list(temp_list))
  }
  names(Trait_peak_list) <- names(Intergenic_trait_peaks)
  
  Celltypes_marker_peaks_matrix <- matrix(0, nrow = 132008, ncol = length(Celltypes))
  rownames(Celltypes_marker_peaks_matrix) <- peaksets_df$peak_id
  colnames(Celltypes_marker_peaks_matrix) <- Celltypes
  for (i in Celltypes) {
    temp <- markerPeaksList_Celltype[[i]]
    temp <- paste0(temp$seqnames,"_",temp$start,"_",temp$end)
    print(length(temp) == sum(temp %in% rownames(Celltypes_marker_peaks_matrix)))
    Celltypes_marker_peaks_matrix[which(rownames(Celltypes_marker_peaks_matrix) %in% temp), i] <- 1
  }
  
  Tissues <- names(markerPeaksList_tissue)
  Tissues_marker_peaks_matrix <- matrix(0, nrow = 132008, ncol = length(Tissues))
  rownames(Tissues_marker_peaks_matrix) <- paste0(peaksets_df$peak_id)
  colnames(Tissues_marker_peaks_matrix) <- Tissues
  for (i in Tissues) {
    temp <- markerPeaksList_tissue[[i]]
    temp <- paste0(temp$seqnames,"_",temp$start,"_",temp$end)
    print(length(temp) == sum(temp %in% rownames(Tissues_marker_peaks_matrix)))
    Tissues_marker_peaks_matrix[which(rownames(Tissues_marker_peaks_matrix) %in% temp), i] <- 1
  }
  
  
  
  Celltype_trait_enrich_matrix_ham <- matrix(0,
                                             nrow = length(Trait_peak_list),
                                             ncol = length(Celltype_marker_peaks_random_list))
  rownames(Celltype_trait_enrich_matrix_ham) <- names(Trait_peak_list)
  colnames(Celltype_trait_enrich_matrix_ham) <- names(Celltype_marker_peaks_random_list)
  for (j in names(Celltype_marker_peaks_random_list)) {
    start <- Sys.time()
    for (i in names(Trait_peak_list)) {
      print(paste0(i, "  ", j))
      Random_common <- lapply(1:10000, function(x){
        temp <- Celltype_marker_peaks_random_list[[j]][[x]]
        temp <- temp + Trait_peak_list[[i]]
        sum(temp == 2)
      })
      Random_common <- unlist(Random_common)
      temp <- Celltypes_marker_peaks_matrix[,j]
      temp <- temp + Trait_peak_list[[i]]
      Real_common <- sum(temp == 2)
      Celltype_trait_enrich_matrix_ham[i,j] <- sum(Random_common >= Real_common) / 10000
    }
    end <- Sys.time()
  }
  
  Celltype_trait_enrich_matrix_ham_binary <- matrix(0,
                                                    nrow = length(Trait_peak_list),
                                                    ncol = length(Celltype_marker_peaks_random_list))
  rownames(Celltype_trait_enrich_matrix_ham_binary) <- names(Trait_peak_list)
  colnames(Celltype_trait_enrich_matrix_ham_binary) <- names(Celltype_marker_peaks_random_list)
  for (i in 1:nrow(Celltype_trait_enrich_matrix_ham_binary)) {
    for (j in 1:ncol(Celltype_trait_enrich_matrix_ham_binary)) {
      p <- Celltype_trait_enrich_matrix_ham[i,j]
      if (p < 0.05) {
        Celltype_trait_enrich_matrix_ham_binary[i, j] <- 1
      }
    }
  }
  trait_sig_ham <- rowSums(Celltype_trait_enrich_matrix_ham_binary)
  trait_sig_ham <- trait_sig_ham[which(trait_sig_ham > 0)]
  trait_order <- sort(names(trait_sig_ham))
  trait_order <- data.frame(trait = trait_order)
  trait_order$trait_group <- c("Physiological traits", "Physiological traits", "Physiological traits", "Physiological traits", "Morphological traits",
                               "Morphological traits", "Morphological traits", "Morphological traits", "Physiological traits", "Physiological traits",
                               "Physiological traits", "Morphological traits", "Physiological traits", "Physiological traits", "Yield-related traits",
                               "Yield-related traits", "Yield-related traits", "Yield-related traits", "Yield-related traits", "Yield-related traits",
                               "Yield-related traits", "Morphological traits", "Morphological traits", "Physiological traits", "Yield-related traits",
                               "Morphological traits", "Morphological traits", "Yield-related traits", "Physiological traits", "Stress-response traits",
                               "Ecological traits", "Stress-response traits", "Physiological traits", "Yield-related traits", "Morphological traits",
                               "Morphological traits", "Morphological traits", "Yield-related traits"
  )
  trait_order$trait_group <- factor(trait_order$trait_group,
                                    levels = c("Yield-related traits", "Stress-response traits",
                                               "Physiological traits", "Morphological traits",
                                               "Ecological traits"))
  trait_order <- trait_order[order(trait_order$trait_group, trait_order$trait),]
  celltype_order <- c("Branch meristems (BM)", "Cryptic bract/bract (cb/b)", "Inflorescence meristem (IM)", "Rachis", "Spikelet meristem (SM)",
                      "Mestome sheath", "Procambium", "Large parenchyma (PO)", "Fiber", "Large parenchyma (MO)", "Mesophyll (MO)",
                      "Mesophyll precursor", "Mesophyll (PO)", "Mesophyll initial", "Phloem",
                      "Cortex", "Endodermis", "Pericycle", "Root hair"
  )
  
  Celltype_trait_enrich_matrix_ham_2 <- Celltype_trait_enrich_matrix_ham[trait_order$trait, celltype_order]
  if (!is.null(Celltype_trait_enrich_matrix_ham_2)) {
    pmt <- Celltype_trait_enrich_matrix_ham_2
    pmt <- as.data.frame(pmt)
    ssmt <- Celltype_trait_enrich_matrix_ham_2 <= 0.0001
    pmt[ssmt] <- "****"
    smt <- Celltype_trait_enrich_matrix_ham_2 > 0.0001 & Celltype_trait_enrich_matrix_ham_2 <= 0.001
    pmt[smt] <- "***"
    smt2 <- Celltype_trait_enrich_matrix_ham_2 > 0.001 & Celltype_trait_enrich_matrix_ham_2 <= 0.01
    pmt[smt2] <- "**"
    smt3 <- Celltype_trait_enrich_matrix_ham_2 > 0.01 & Celltype_trait_enrich_matrix_ham_2 <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  
  Celltype_trait_enrich_matrix_ham_log <- -log10(Celltype_trait_enrich_matrix_ham_2)
  Celltype_trait_enrich_matrix_ham_log[which(Celltype_trait_enrich_matrix_ham_log == -log10(1))] <- NA
  pdf("Figure7_celltype_enrichment_matrix_peak_hamming.pdf", height = 10, width = 12)
  pheatmap(Celltype_trait_enrich_matrix_ham_log,
           color = colorRampPalette(c("yellow", "red"))(20),
           scale = "none", cluster_row = F,
           cluster_col = F,
           treeheight_col = 0,
           display_numbers = pmt,
           fontsize_number = 14,
           number_color = "black",
           cellwidth = 15, cellheight = 10)
  dev.off()
  saveRDS(Celltype_trait_enrich_matrix_ham_2, "Figure7_Celltype_trait_enrich_matrix_ham.rds")
  
  if (!is.null(Celltype_trait_enrich_matrix_ham_2)) {
    pmt <- Celltype_trait_enrich_matrix_ham_2
    pmt <- as.data.frame(pmt)
    ssmt <- Celltype_trait_enrich_matrix_ham_2 <= 0.0001
    pmt[ssmt] <- 1
    smt <- Celltype_trait_enrich_matrix_ham_2 > 0.0001 & Celltype_trait_enrich_matrix_ham_2 <= 0.001
    pmt[smt] <- 1
    smt2 <- Celltype_trait_enrich_matrix_ham_2 > 0.001 & Celltype_trait_enrich_matrix_ham_2 <= 0.01
    pmt[smt2] <- 1
    smt3 <- Celltype_trait_enrich_matrix_ham_2 > 0.01 & Celltype_trait_enrich_matrix_ham_2 <= 0.05
    pmt[smt3] <- 1
    pmt[!ssmt&!smt&!smt2&!smt3] <- 0
  }
  
  temp <- data.frame(colSums(pmt))
  temp$Celltype <- rownames(temp)
  colnames(temp)[1] <- "Count"
  temp <- temp[order(temp$Count, decreasing = T),]
  temp$Celltype <- factor(temp$Celltype, levels = temp$Celltype)
  
  pdf("Figure7.celltype_enriched_trait_count_peak_hamming.pdf", height = 4, width = 6.5)
  ggplot(data = temp, aes(x = Celltype, y = Count)) +
    geom_bar(stat = "identity", width = 0.6, fill = "#476D87") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "", y = "Enriched trait count")
  dev.off()
  
  Celltype_trait_enrich_matrix_ham_binary_df <- reshape2::melt(Celltype_trait_enrich_matrix_ham)
  # Celltype_trait_enrich_matrix_ham_binary_df <- Celltype_trait_enrich_matrix_ham_binary_df[which(Celltype_trait_enrich_matrix_ham_binary_df$value < 0.05),]
  Celltype_trait_enrich_matrix_ham_binary_df$pair <- paste0(Celltype_trait_enrich_matrix_ham_binary_df$Var1,
                                                            "_",
                                                            Celltype_trait_enrich_matrix_ham_binary_df$Var2)
  celltype_enrichment_matrix_binary_df <- reshape2::melt(celltype_enrichment_matrix2)
  # celltype_enrichment_matrix_binary_df <- celltype_enrichment_matrix_binary_df[which(celltype_enrichment_matrix_binary_df$value < 0.05),]
  celltype_enrichment_matrix_binary_df$pair <- paste0(celltype_enrichment_matrix_binary_df$Var2,
                                                      "_",
                                                      celltype_enrichment_matrix_binary_df$Var1)
  
  rownames(Celltype_trait_enrich_matrix_ham_binary_df) <- Celltype_trait_enrich_matrix_ham_binary_df$pair
  rownames(celltype_enrichment_matrix_binary_df) <- celltype_enrichment_matrix_binary_df$pair
  celltype_enrichment_matrix_binary_df <- celltype_enrichment_matrix_binary_df[Celltype_trait_enrich_matrix_ham_binary_df$pair,]
  
  temp <- data.frame(Enrichment_P = celltype_enrichment_matrix_binary_df$value,
                     Overlap_Imputation_P = Celltype_trait_enrich_matrix_ham_binary_df$value)
  saveRDS(temp, "Figure7_hamming_enrich_P.rds")
  pdf("Figure7_Cell-type-Trait_Association_P_point.pdf", width = 4.5, height = 4.6)
  ggplot(data = temp, aes(x = -log10(Enrichment_P), y = -log10(Overlap_Imputation_P))) +
    geom_point(alpha = 0.6) +
    geom_smooth(colour = "red") +
    theme_bw() +
    ggtitle("Cell-type-Trait Association") +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)
          ) +
    ggpubr::stat_cor(size = 4)
  dev.off()
  
  library(VennDiagram)
  Tvenn <- venn.diagram(list("hamming" = Celltype_trait_enrich_matrix_ham_binary_df$pair,
                             "enrichment" = celltype_enrichment_matrix_binary_df$pair),
                        filename = NULL,
                        lwd = 1, lty = 2,
                        col = c("#00BFC4", "#F8766D") ,
                        fill = c("#00BFC4", "#F8766D"),
                        cat.col = c("black", "black"),
                        main = "Celltype", main.cex = 2,
                        reverse = TRUE)
  grid.draw(Tvenn)
  
  Tissue_trait_enrich_matrix_ham <- matrix(0,
                                           nrow = length(Trait_peak_list),
                                           ncol = length(Tissue_marker_peaks_random_list))
  rownames(Tissue_trait_enrich_matrix_ham) <- names(Trait_peak_list)
  colnames(Tissue_trait_enrich_matrix_ham) <- names(Tissue_marker_peaks_random_list)
  for (j in names(Tissue_marker_peaks_random_list)) {
    start <- Sys.time()
    for (i in names(Trait_peak_list)) {
      print(paste0(i, "  ", j))
      Random_common <- lapply(1:10000, function(x){
        temp <- Tissue_marker_peaks_random_list[[j]][[x]]
        temp <- temp + Trait_peak_list[[i]]
        sum(temp == 2)
      })
      Random_common <- unlist(Random_common)
      temp <- Tissues_marker_peaks_matrix[,j]
      temp <- temp + Trait_peak_list[[i]]
      Real_common <- sum(temp == 2)
      Tissue_trait_enrich_matrix_ham[i,j] <- sum(Random_common >= Real_common) / 10000
    }
    end <- Sys.time()
  }
  
  Tissue_trait_enrich_matrix_ham_binary <- matrix(0,
                                                  nrow = length(Trait_peak_list),
                                                  ncol = length(Tissue_marker_peaks_random_list))
  rownames(Tissue_trait_enrich_matrix_ham_binary) <- names(Trait_peak_list)
  colnames(Tissue_trait_enrich_matrix_ham_binary) <- names(Tissue_marker_peaks_random_list)
  for (i in 1:nrow(Tissue_trait_enrich_matrix_ham_binary)) {
    for (j in 1:ncol(Tissue_trait_enrich_matrix_ham_binary)) {
      p <- Tissue_trait_enrich_matrix_ham[i,j]
      if (p < 0.05) {
        Tissue_trait_enrich_matrix_ham_binary[i, j] <- 1
      }
    }
  }
  trait_sig_ham <- rowSums(Tissue_trait_enrich_matrix_ham_binary)
  trait_sig_ham <- trait_sig_ham[which(trait_sig_ham > 0)]
  
  trait_order <- sort(names(trait_sig_ham))
  trait_order <- data.frame(trait = trait_order)
  trait_order$trait_group <- c("Morphological traits", "Morphological traits", "Morphological traits", "Morphological traits",
                               "Physiological traits", "Yield-related traits", "Yield-related traits", "Yield-related traits",
                               "Morphological traits", "Morphological traits", "Morphological traits", "Physiological traits",
                               "Stress-response traits", "Morphological traits"
  )
  trait_order$trait_group <- factor(trait_order$trait_group,
                                    levels = c("Yield-related traits", "Stress-response traits",
                                               "Physiological traits", "Morphological traits",
                                               "Ecological traits"))
  trait_order <- trait_order[order(trait_order$trait_group, trait_order$trait),]
  Tissue_trait_enrich_matrix_ham_2 <- Tissue_trait_enrich_matrix_ham[trait_order$trait, c("IM", "Shoot", "Leaf", "Root")]
  
  if (!is.null(Tissue_trait_enrich_matrix_ham_2)) {
    pmt <- Tissue_trait_enrich_matrix_ham_2
    pmt <- as.data.frame(pmt)
    ssmt <- Tissue_trait_enrich_matrix_ham_2 <= 0.0001
    pmt[ssmt] <- "****"
    smt <- Tissue_trait_enrich_matrix_ham_2 > 0.0001 & Tissue_trait_enrich_matrix_ham_2 <= 0.001
    pmt[smt] <- "***"
    smt2 <- Tissue_trait_enrich_matrix_ham_2 > 0.001 & Tissue_trait_enrich_matrix_ham_2 <= 0.01
    pmt[smt2] <- "**"
    smt3 <- Tissue_trait_enrich_matrix_ham_2 > 0.01 & Tissue_trait_enrich_matrix_ham_2 <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  Tissue_trait_enrich_matrix_ham_2_log <- -log10(Tissue_trait_enrich_matrix_ham_2)
  Tissue_trait_enrich_matrix_ham_2_log[which(Tissue_trait_enrich_matrix_ham_2_log == -log10(1))] <- NA
  pdf("Figure7_tissue_enrichment_matrix_peak_hamming.pdf", height = 6, width = 6)
  pheatmap(Tissue_trait_enrich_matrix_ham_2_log,
           color = colorRampPalette(c("yellow", "red"))(20),
           scale = "none", cluster_row = F,
           cluster_col = F,      
           treeheight_col = 0,
           display_numbers = pmt,
           fontsize_number = 12,
           number_color = "black",
           cellwidth = 15, cellheight = 10)
  dev.off()
  
  
  Tissue_trait_enrich_matrix_ham_2 <- Tissue_trait_enrich_matrix_ham[which(rownames(Tissue_trait_enrich_matrix_ham) %in% names(trait_sig_ham)),]
  Tissue_trait_enrich_matrix_ham_2 <- Tissue_trait_enrich_matrix_ham_2[trait_order2,]
  if (!is.null(Tissue_trait_enrich_matrix_ham_2)) {
    pmt <- Tissue_trait_enrich_matrix_ham_2
    pmt <- as.data.frame(pmt)
    ssmt <- Tissue_trait_enrich_matrix_ham_2 <= 0.0001
    pmt[ssmt] <- "****"
    smt <- Tissue_trait_enrich_matrix_ham_2 > 0.0001 & Tissue_trait_enrich_matrix_ham_2 <= 0.001
    pmt[smt] <- "***"
    smt2 <- Tissue_trait_enrich_matrix_ham_2 > 0.001 & Tissue_trait_enrich_matrix_ham_2 <= 0.01
    pmt[smt2] <- "**"
    smt3 <- Tissue_trait_enrich_matrix_ham_2 > 0.01 & Tissue_trait_enrich_matrix_ham_2 <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  
  pdf("Figure7_tissue_enrichment_matrix_peak_hamming.pdf", height = 10, width = 12)
  pheatmap(-log10(Tissue_trait_enrich_matrix_ham_2),
           color = colorRampPalette(c("yellow", "red"))(20),
           scale = "none", cluster_row = F,
           cluster_col = T, treeheight_col = 0,             
           display_numbers = pmt,
           fontsize_number = 12,
           number_color = "black",
           cellwidth = 15, cellheight = 10)
  dev.off()
  saveRDS(Tissue_trait_enrich_matrix_ham_2,
          "Figure7_Tissue_trait_enrich_matrix_ham.rds")
  
  Tissue_trait_enrich_matrix_ham_binary_df <- reshape2::melt(Tissue_trait_enrich_matrix_ham)
  Tissue_trait_enrich_matrix_ham_binary_df$pair <- paste0(Tissue_trait_enrich_matrix_ham_binary_df$Var1,
                                                          "_",
                                                          Tissue_trait_enrich_matrix_ham_binary_df$Var2)
  tissue_enrichment_matrix_binary_df <- reshape2::melt(tissue_enrichment_matrix2)
  tissue_enrichment_matrix_binary_df$pair <- paste0(tissue_enrichment_matrix_binary_df$Var2,
                                                    "_",
                                                    tissue_enrichment_matrix_binary_df$Var1)
  
  rownames(Tissue_trait_enrich_matrix_ham_binary_df) <- Tissue_trait_enrich_matrix_ham_binary_df$pair
  rownames(tissue_enrichment_matrix_binary_df) <- tissue_enrichment_matrix_binary_df$pair
  tissue_enrichment_matrix_binary_df <- tissue_enrichment_matrix_binary_df[Tissue_trait_enrich_matrix_ham_binary_df$pair,]
  
  temp <- data.frame(Enrichment_P = tissue_enrichment_matrix_binary_df$value,
                     Overlap_Imputation_P = Tissue_trait_enrich_matrix_ham_binary_df$value)
  saveRDS(temp,
          "Figure7_Tissue_hamming_enrich_P.rds")
  pdf("Figure7_Tissue-Trait_Association_P_point.pdf", width = 4.5, height = 4.6)
  ggplot(data = temp, aes(x = -log10(Enrichment_P), y = -log10(Overlap_Imputation_P))) +
    geom_point(alpha = 0.6) +
    geom_smooth(colour = "red") +
    theme_bw() +
    ggtitle("Tissue-Trait Association") +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    ggpubr::stat_cor(size = 4)
  dev.off()
  
  
  library(VennDiagram)
  Tvenn <- venn.diagram(list("hamming" = Tissue_trait_enrich_matrix_ham_binary_df$pair,
                             "enrichment" = tissue_enrichment_matrix_binary_df$pair),
                        filename = NULL,
                        lwd = 1, lty = 2,
                        col = c("#00BFC4", "#F8766D") ,
                        fill = c("#00BFC4", "#F8766D"),
                        cat.col = c("black", "black"),
                        main = "Celltype", main.cex = 2,
                        reverse = TRUE)
  grid.draw(Tvenn)
  
}


####
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

if (FALSE) {
  RNA <- readRDS("./Rice-snRNA/all_data_annotated.rds")
  RNA$Celltype[which(RNA$Celltype == "Branch meristems(BM)")] <- "Branch meristems (BM)"
  RNA$Celltype[which(RNA$Celltype == "Inflorescence meristem(IM)")] <- "Inflorescence meristem (IM)"
  RNA$Celltype[which(RNA$Celltype == "Spikelet meristem(SM)")] <- "Spikelet meristem (SM)"
  RNA$Celltype[which(RNA$Celltype == "Lemma(le)")] <- "Lemma (le)"
  RNA$Celltype[which(RNA$Celltype == "Cryptic bract/bract(cb/b)")] <- "Cryptic bract/bract (cb/b)"
  
  common_celltypes <- intersect(unique(RNA$Celltype),
                                unique(proj_pass_filter$Celltype))
  
  # RNA <- subset(RNA, subset = Celltype %in% common_celltypes)
  Idents(RNA) <- RNA$Celltype
  RNA_celltype_marker <- FindAllMarkers(RNA, logfc.threshold = 1, only.pos = T)
  RNA_celltype_marker_TF_filter <- RNA_celltype_marker[which(RNA_celltype_marker$gene %in% Annotation_genes$Locus_ID),]
  temp <- as.data.frame(table(RNA_celltype_marker_TF_filter$cluster))
  temp <- temp[order(temp$Freq, decreasing = T),]
  temp <- temp[which(temp$Var1 %in% common_celltypes),]
  temp$Var1 <- factor(temp$Var1, levels = temp$Var1)
  pdf("Figure7_RNA_细胞类型差异表达的TF数目.pdf", width = 7, height = 6)
  ggplot(data = temp, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(data = temp, aes(x = Var1, y = Freq + 0.5, label = Freq)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(y = "TF marker count", x = "")
  dev.off()
}

#### 导入TF-Target
{
  TF_target_filtered <- readRDS("./TF_regulation_filtered.rds")
  TF_regulation <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/regulation_merged_Osj.txt",
                            sep = "\t", header = F)
  table(TF_regulation$V5)
  colnames(TF_regulation) <- c("TF_MSU", "Type", "Target_MSU", "Species", "Methods")
  TF_regulation$TF_Target_MSU <- paste0(TF_regulation$TF_MSU ,"_" , TF_regulation$Target_MSU)
  TF_regulation <- TF_regulation[!duplicated(TF_regulation$TF_Target_MSU),]
  
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
  
  TF_regulation2 <- merge(x = TF_regulation, y = MSU_RAP, by.x = "TF_MSU", by.y = "MSU")
  colnames(TF_regulation2)[ncol(TF_regulation2)] <- "TF_RAP"
  TF_regulation2 <- merge(x = TF_regulation2, y = MSU_RAP, by.x = "Target_MSU", by.y = "MSU")
  colnames(TF_regulation2)[ncol(TF_regulation2)] <- "Target_RAP"
  TF_regulation <- as.data.frame(na.omit(TF_regulation2))
  TF_regulation$Pairs <- paste0(TF_regulation$TF_RAP,"_",TF_regulation$Target_RAP)
  sum(duplicated(TF_regulation$Pairs))
  TF_regulation_all <- TF_regulation
  TF_regulation_all <- TF_regulation_all[!duplicated(TF_regulation_all$Pairs),]
}

###### 产生SNP版本的基因组
{
  fasta_sequences <- readDNAStringSet("./Ref/Oryza_sativa.IRGSP-1.0.fa")
  names(fasta_sequences)
  names(fasta_sequences) <- c("1", "Mt", "Pt", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10", "11", "12")
  for (i in c("1", "Mt", "Pt", "2", "3", "4", "5",
              "6", "7", "8", "9", "10", "11", "12")) {
    temp <- getSeq(BSgenome.OSativa.NCBI.IRGSPv1.0, names = i, as.character = TRUE)
    temp2 <- fasta_sequences[i]
    temp2 <- as.character(temp2)
    print(temp == temp2)
  }
  
  SNP_ref_site <- c()
  for (i in 1:nrow(Unique_SNP)) {
    print(i)
    temp <- Unique_SNP[i,]
    SNP_site <- fasta_sequences[as.character(temp$chr)]
    SNP_site <- subseq(SNP_site, start = as.numeric(temp$start), end = as.numeric(temp$end))
    SNP_site <- as.character(SNP_site)
    SNP_site <- substr(SNP_site, 1, 1)
    SNP_ref_site <- c(SNP_ref_site, SNP_site)
  }
  Unique_SNP$SNP_ref_site <- SNP_ref_site
  sum(Unique_SNP$REF == Unique_SNP$SNP_ref_site) # 45970
  sum(Unique_SNP$ALT == Unique_SNP$SNP_ref_site) # 201
  temp <- Unique_SNP[which(Unique_SNP$ALT == Unique_SNP$SNP_ref_site),]
  
  fasta_sequences_without_mt_pt <- fasta_sequences[c("1", "2", "3", "4", "5", "6",
                                                     "7", "8", "9", "10", "11", "12")]
  fasta_sequences_without_mt_pt_ref <- fasta_sequences_without_mt_pt
  fasta_sequences_without_mt_pt_alt <- fasta_sequences_without_mt_pt
  for (i in 1:nrow(Unique_SNP)) {
    print(i)
    subseq(fasta_sequences_without_mt_pt_ref[as.character(Unique_SNP$chr[i])], start = as.numeric(Unique_SNP$start[i]), end = as.numeric(Unique_SNP$start[i])) <- DNAString(Unique_SNP$REF[i])
    subseq(fasta_sequences_without_mt_pt_alt[as.character(Unique_SNP$chr[i])], start = as.numeric(Unique_SNP$start[i]), end = as.numeric(Unique_SNP$start[i])) <- DNAString(Unique_SNP$ALT[i])
  }
  
  writeXStringSet(fasta_sequences_without_mt_pt_alt, filepath = paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/REF_ALT_BSgenome/fasta_sequences_without_mt_pt_alt.fasta"), format = "fasta")
  writeXStringSet(fasta_sequences_without_mt_pt_ref, filepath = paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/REF_ALT_BSgenome/fasta_sequences_without_mt_pt_ref.fasta"), format = "fasta")
}

# 产生BSgenome
{
  library(BSgenome)
  setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/REF_ALT_BSgenome")
  forgeBSgenomeDataPkg("./BSgenome.Oryza.sativa.IRGSP.1.0.REF-seed",
                       seqs_srcdir = getwd(), destdir = getwd(),
                       verbose = T)
  system(command = "R CMD build BSgenome.Oryza.sativa.IRGSP.1.0.REF")
  system(command = "R CMD check --no-manual BSgenome.Oryza.sativa.IRGSP.1.0.REF_1.0.0.tar.gz")
  system(command = "R CMD INSTALL BSgenome.Oryza.sativa.IRGSP.1.0.REF_1.0.0.tar.gz")
  
  forgeBSgenomeDataPkg("./BSgenome.Oryza.sativa.IRGSP.1.0.ALT-seed",
                       seqs_srcdir = getwd(), destdir = getwd(),
                       verbose = T)
  system(command = "R CMD build BSgenome.Oryza.sativa.IRGSP.1.0.ALT")
  system(command = "R CMD check --no-manual BSgenome.Oryza.sativa.IRGSP.1.0.ALT_1.0.0.tar.gz")
  system(command = "R CMD INSTALL BSgenome.Oryza.sativa.IRGSP.1.0.ALT_1.0.0.tar.gz")
}

library(BSgenome.Oryza.sativa.IRGSP.1.0.REF)
library(BSgenome.Oryza.sativa.IRGSP.1.0.ALT)

SNP_8240_peak <- data.frame(SNP_chr = Intergenic_overlap_with_peaks_filtered$chr,
                            SNP_loci = Intergenic_overlap_with_peaks_filtered$start,
                            SNP = Intergenic_overlap_with_peaks_filtered$SNP_CHR_BP,
                            SNP_REF = Intergenic_overlap_with_peaks_filtered$REF,
                            SNP_ALT = Intergenic_overlap_with_peaks_filtered$ALT,
                            SNP_nearestGene = Intergenic_overlap_with_peaks_filtered$nearestGene,
                            SNP_peak_id = Intergenic_overlap_with_peaks_filtered$peak_id,
                            SNP_trait = Intergenic_overlap_with_peaks_filtered$TRAIT_Renamed)

#### SNP range
if (FALSE) {
  SNP_range <- data.frame(chr = SNP_8240_peak$SNP_chr,
                          start = SNP_8240_peak$SNP_loci,
                          end = SNP_8240_peak$SNP_loci + 1,
                          strand = "*")
  SNP_range <- makeGRangesFromDataFrame(SNP_range,
                                        seqinfo = SeqinfoForBSGenome(BSgenome.Oryza.sativa.IRGSP.1.0.REF))
  SNP_range <- extendGR(SNP_range, downstream = 49, upstream = 50)
  SNP_range$SNP_REF <- SNP_8240_peak$SNP_REF
  SNP_range$SNP_ALT <- SNP_8240_peak$SNP_ALT
  SNP_range$SNP_nearestGene <- SNP_8240_peak$SNP_nearestGene
  SNP_range$SNP_peak_id <- SNP_8240_peak$SNP_peak_id
  SNP_range$SNP_trait <- SNP_8240_peak$SNP_trait
  SNP_range$SNP <- SNP_8240_peak$SNP
  SNP_range$SNP_chr <- SNP_8240_peak$SNP_chr
  SNP_range$SNP_loci <- SNP_8240_peak$SNP_loci
  
  #### REF BSgenome 与 Motif 比对
  {
    motifPositions <- motifmatchr::matchMotifs(
      pwms = motif_all,
      subject = SNP_range,
      genome = BSgenome.Oryza.sativa.IRGSP.1.0.REF, 
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
    overlapMotifs <- findOverlaps(SNP_range, allPositions, ignore.strand = TRUE)
    overlapMotifs <- as.data.frame(overlapMotifs)
    peak_id <- SNP_range$SNP_peak_id
    allPositions$TF <- allPositions@ranges@NAMES
    allPositions@ranges@NAMES <- NULL
    allPositions_df <- as.data.frame(allPositions)
    ref_peak_TF <- data.frame(SNP_trait = SNP_range$SNP_trait[overlapMotifs$queryHits],
                              SNP_chr = SNP_range$SNP_chr[overlapMotifs$queryHits],
                              SNP_loci = SNP_range$SNP_loci[overlapMotifs$queryHits],
                              SNP = SNP_range$SNP[overlapMotifs$queryHits],
                              SNP_REF = SNP_range$SNP_REF[overlapMotifs$queryHits],
                              SNP_ALT = SNP_range$SNP_ALT[overlapMotifs$queryHits],
                              SNP_nearestGene = SNP_range$SNP_nearestGene[overlapMotifs$queryHits],
                              ref_snp_motif_chr = allPositions_df$seqnames[overlapMotifs$subjectHits],
                              ref_snp_motif_start = allPositions_df$start[overlapMotifs$subjectHits],
                              ref_snp_motif_end = allPositions_df$end[overlapMotifs$subjectHits],
                              ref_snp_TF = allPositions_df$TF[overlapMotifs$subjectHits],
                              ref_snp_motif_strand = allPositions_df$strand[overlapMotifs$subjectHits],
                              ref_snp_motif_score = allPositions_df$score[overlapMotifs$subjectHits],
                              ref_snp_motif_width = allPositions_df$width[overlapMotifs$subjectHits],
                              ref_snp_motif_peakid = peak_id[overlapMotifs$queryHits],
                              ref_snp_gene = SNP_range$SNP_nearestGene[overlapMotifs$queryHits]
    )
    length(unique(ref_peak_TF$ref_snp_motif_peakid)) # ref snp 3963
    gene_annotation_df <- as.data.frame(gene_annotation@listData[["genes"]])
    rownames(gene_annotation_df) <- gene_annotation_df$symbol
    ref_peak_TF$ref_snp_gene_strand <- as.character(gene_annotation_df[ref_peak_TF$SNP_nearestGene, "strand"])
    ref_peak_TF <- ref_peak_TF[which(as.character(ref_peak_TF$ref_snp_motif_strand) == ref_peak_TF$ref_snp_gene_strand),]
    ref_peak_TF_SNP <- ref_peak_TF
    ref_peak_TF_SNP <- ref_peak_TF_SNP[which(ref_peak_TF_SNP$SNP_loci >= ref_peak_TF_SNP$ref_snp_motif_start &
                                             ref_peak_TF_SNP$SNP_loci <= ref_peak_TF_SNP$ref_snp_motif_end),]
  }
  
  #### ALT BSgenome 与 peak Motif 比对
  {
    motifPositions <- motifmatchr::matchMotifs(
      pwms = motif_all,
      subject = SNP_range,
      genome = BSgenome.Oryza.sativa.IRGSP.1.0.ALT, 
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
    overlapMotifs <- findOverlaps(SNP_range, allPositions, ignore.strand = TRUE)
    overlapMotifs <- as.data.frame(overlapMotifs)
    peak_id <- SNP_range$SNP_peak_id
    allPositions$TF <- allPositions@ranges@NAMES
    allPositions@ranges@NAMES <- NULL
    allPositions_df <- as.data.frame(allPositions)
    alt_peak_TF <- data.frame(SNP_trait = SNP_range$SNP_trait[overlapMotifs$queryHits],
                              SNP_chr = SNP_range$SNP_chr[overlapMotifs$queryHits],
                              SNP_loci = SNP_range$SNP_loci[overlapMotifs$queryHits],
                              SNP = SNP_range$SNP[overlapMotifs$queryHits],
                              SNP_REF = SNP_range$SNP_REF[overlapMotifs$queryHits],
                              SNP_ALT = SNP_range$SNP_ALT[overlapMotifs$queryHits],
                              SNP_nearestGene = SNP_range$SNP_nearestGene[overlapMotifs$queryHits],
                              alt_snp_motif_chr = allPositions_df$seqnames[overlapMotifs$subjectHits],
                              alt_snp_motif_start = allPositions_df$start[overlapMotifs$subjectHits],
                              alt_snp_motif_end = allPositions_df$end[overlapMotifs$subjectHits],
                              alt_snp_TF = allPositions_df$TF[overlapMotifs$subjectHits],
                              alt_snp_motif_strand = allPositions_df$strand[overlapMotifs$subjectHits],
                              alt_snp_motif_score = allPositions_df$score[overlapMotifs$subjectHits],
                              alt_snp_motif_width = allPositions_df$width[overlapMotifs$subjectHits],
                              alt_snp_motif_peakid = peak_id[overlapMotifs$queryHits]
    )
    length(unique(alt_peak_TF$alt_snp_motif_peakid)) # alt snp 3964
    gene_annotation_df <- as.data.frame(gene_annotation@listData[["genes"]])
    rownames(gene_annotation_df) <- gene_annotation_df$symbol
    alt_peak_TF$alt_snp_gene_strand <- as.character(gene_annotation_df[alt_peak_TF$SNP_nearestGene, "strand"])
    alt_peak_TF <- alt_peak_TF[which(as.character(alt_peak_TF$alt_snp_motif_strand) == alt_peak_TF$alt_snp_gene_strand),]
    alt_peak_TF_SNP <- alt_peak_TF
    alt_peak_TF_SNP <- alt_peak_TF_SNP[which(alt_peak_TF_SNP$SNP_loci >= alt_peak_TF_SNP$alt_snp_motif_start &
                                             alt_peak_TF_SNP$SNP_loci <= alt_peak_TF_SNP$alt_snp_motif_end),]
  }
  
  ### 整合 ref 和 alt 的结合结果
  {
    ref_peak_TF_SNP_2 <- ref_peak_TF_SNP
    ref_peak_TF_SNP_2 <- split.data.frame(ref_peak_TF_SNP_2, f = as.factor(paste0(ref_peak_TF_SNP_2$SNP,
                                                                                  "_",
                                                                                  ref_peak_TF_SNP_2$ref_snp_TF)))
    ref_peak_TF_SNP_2 <- lapply(ref_peak_TF_SNP_2, function(x){
      x[which.max(x$ref_snp_motif_score),]
    })
    ref_peak_TF_SNP_2 <- rbindlist(ref_peak_TF_SNP_2)
    ref_peak_TF_SNP_2 <- apply(ref_peak_TF_SNP_2, 2, as.character)
    ref_peak_TF_SNP_2 <- as.data.frame(ref_peak_TF_SNP_2)
    
    
    alt_peak_TF_SNP_2 <- alt_peak_TF_SNP
    alt_peak_TF_SNP_2 <- split.data.frame(alt_peak_TF_SNP_2, f = as.factor(paste0(alt_peak_TF_SNP_2$SNP,
                                                                                  "_",
                                                                                  alt_peak_TF_SNP_2$alt_snp_TF)))
    alt_peak_TF_SNP_2 <- lapply(alt_peak_TF_SNP_2, function(x){
      x[which.max(x$alt_snp_motif_score),]
    })
    alt_peak_TF_SNP_2 <- rbindlist(alt_peak_TF_SNP_2)
    alt_peak_TF_SNP_2 <- apply(alt_peak_TF_SNP_2, 2, as.character)
    alt_peak_TF_SNP_2 <- as.data.frame(alt_peak_TF_SNP_2)
    
    length(unique(ref_peak_TF_SNP_2$SNP))
    length(unique(alt_peak_TF_SNP_2$SNP))
    
    ref_peak_TF_SNP_2$TF <- ref_peak_TF_SNP_2$ref_snp_TF
    alt_peak_TF_SNP_2$TF <- alt_peak_TF_SNP_2$alt_snp_TF
    ref_alt_peak_TF_SNP <- merge(x = ref_peak_TF_SNP_2, y = alt_peak_TF_SNP_2,
                                 by = c("SNP_chr", "SNP_loci", "SNP", "SNP_REF", "SNP_ALT", "SNP_nearestGene", "SNP_trait", "TF"), all = TRUE)
    ref_alt_peak_TF_SNP[is.na(ref_alt_peak_TF_SNP)] <- ""
  }
  
  #### 对于每一个细胞类型和表型对，分别统计
  celltype_enrichment_matrix_binary_df <- reshape2::melt(celltype_enrichment_matrix_binary)
  celltype_enrichment_matrix_binary_df <- celltype_enrichment_matrix_binary_df[which(celltype_enrichment_matrix_binary_df$value == 1),]
  colnames(celltype_enrichment_matrix_binary_df) <- c("Celltype", "Trait", "SigIf")
  Celltype_Stat <- c()
  Celltype_Trait_SNP_remove <- c()
  Celltype_Trait_SNP_overlap <- c() # 和细胞类型有overlap
  Celltype_Trait_SNP_all <- c() # 显著富集的trait的所有SNP
  for (i in 1:nrow(celltype_enrichment_matrix_binary_df)) {
    print(i)
    trait <- celltype_enrichment_matrix_binary_df[i, 2]
    trait <- as.character(trait)
    trait_peaks <- Intergenic_trait_peaks[[trait]]
    trait_peaks <- unique(trait_peaks$peak_id)
    trait_SNP <- Intergenic_trait_peaks[[trait]]
    Celltype_Trait_SNP_all <- c(Celltype_Trait_SNP_all, unique(trait_SNP$SNP_CHR_BP))
    cell <- celltype_enrichment_matrix_binary_df[i, 1]
    cell <- as.character(cell)
    cell_markers <- markerPeaksList_Celltype[[cell]]
    cell_markers <- as.data.frame(cell_markers)
    cell_markers <- paste0(cell_markers$seqnames, "_", cell_markers$start, "_", cell_markers$end)
    common_peaks <- intersect(cell_markers, trait_peaks)
    for (peak in common_peaks) {
      trait_SNP <- Intergenic_trait_peaks[[trait]]
      trait_SNP <- trait_SNP[which(trait_SNP$peak_id == peak),]
      trait_SNP <- trait_SNP[!duplicated(trait_SNP$SNP_CHR_BP),]
      rownames(trait_SNP) <- trait_SNP$SNP_CHR_BP
      for (snp in trait_SNP$SNP_CHR_BP) { # 对每一个SNP循环
        temp_trait_SNP <- trait_SNP[snp,]
        Celltype_Trait_SNP_overlap <- as.data.frame(rbind(Celltype_Trait_SNP_overlap,
                                                          data.frame(Celltype = cell,
                                                                     Trait = trait,
                                                                     SNP = temp_trait_SNP$SNP_CHR_BP,
                                                                     Peaks = temp_trait_SNP$peak_id)))
        temp_ref_alt_peak_TF_SNP <- ref_alt_peak_TF_SNP[which(ref_alt_peak_TF_SNP$SNP %in% snp),]
        if (nrow(temp_ref_alt_peak_TF_SNP) == 0) {
          Celltype_Trait_SNP_remove <- as.data.frame(rbind(Celltype_Trait_SNP_remove,
                                                           data.frame(Celltype = cell,
                                                                      Trait = trait,
                                                                      SNP = temp_trait_SNP$SNP_CHR_BP,
                                                                      Peaks = temp_trait_SNP$peak_id)))
        } else {
          print("Yes!")
          Celltype_Stat_temp <- data.frame(Celltype = cell,
                                           Trait = trait,
                                           temp_ref_alt_peak_TF_SNP)
          
          Celltype_Stat <- as.data.frame(rbind(Celltype_Stat,
                                               Celltype_Stat_temp))
          
        }
      }
    }
  }
  
  Celltype_Stat$Peak <- overlapRegions_DF[Celltype_Stat$SNP, 2]
  for (i in 1:nrow(Celltype_Stat)) {
    if (Celltype_Stat$ref_snp_TF[i] == "" & Celltype_Stat$alt_snp_TF[i] != "") {
      Celltype_Stat$TF[i] <- Celltype_Stat$alt_snp_TF[i]
    }
    if (Celltype_Stat$ref_snp_TF[i] != "" & Celltype_Stat$alt_snp_TF[i] == "") {
      Celltype_Stat$TF[i] <- Celltype_Stat$ref_snp_TF[i]
    }
    if (Celltype_Stat$ref_snp_TF[i] != "" & Celltype_Stat$alt_snp_TF[i] != "") {
      if (Celltype_Stat$ref_snp_TF[i] == Celltype_Stat$alt_snp_TF[i]) {
        Celltype_Stat$TF[i] <- Celltype_Stat$ref_snp_TF[i]
      } else {
        Celltype_Stat$TF[i] <- ""
      }
    }
  }
  
  length(unique(Celltype_Stat$SNP)) # 272 SNPs
  length(unique(Celltype_Stat$Celltype)) # 17 Celltypes
  length(unique(Celltype_Stat$Trait)) # 26 Traits
  length(unique(Celltype_Stat$TF)) # 189 TFs
  
  heatmap_value <- c()
  for (i in 1:nrow(Celltype_Stat)) {
    temp_ref <- Celltype_Stat$ref_snp_motif_strand[i] == ""
    temp_alt <- Celltype_Stat$alt_snp_motif_strand[i] == ""
    if (temp_ref == TRUE & temp_alt == FALSE) { # 突变导致结合
      heatmap_value <- c(heatmap_value, "1")
    }
    if (temp_ref == FALSE & temp_alt == FALSE) { # 突变不影响结合
      heatmap_value <- c(heatmap_value, "2")
    }
    if (temp_ref == FALSE & temp_alt == TRUE) { # 突变导致不能结合
      heatmap_value <- c(heatmap_value, "3")
    }
  }
  
  table(Celltype_Stat$ref_snp_motif_strand, Celltype_Stat$alt_snp_motif_strand)
  table(heatmap_value)
  
  openxlsx::write.xlsx(Celltype_Stat, "Figure7_Celltype_Stat_100bp.xlsx")
  
  ### ATAC diff

  # proj_pass_filter$Celltype_2 <- ifelse(proj_pass_filter$Celltype %in% c("Branch meristems (BM)",
  #                                                                        "Cryptic bract/bract (cb/b)",
  #                                                                        "Rachis",
  #                                                                        "Spikelet meristem (SM)"),
  #                                       "YES",)
  markersPeaks_Celltype <- getMarkerFeatures(
    ArchRProj = proj_pass_filter, 
    useMatrix = "PeakMatrix", 
    groupBy = "Celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerPeaksList_Celltype <- getMarkers(markersPeaks_Celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  markerPeaksList_Celltype <- markerPeaksList_Celltype@listData
  
  markerGenesList_celltype <- getMarkers(markersGenes_celltype, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
  markerGenesList_celltype <- markerGenesList_celltype@listData
  
  
  GS_matrix <- GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]
  GS_matrix <- GS_matrix[, which(GeneScoreMatrix@colData@listData[["Celltype"]] %in% rownames(celltype_enrichment_matrix_log))]
  PK_matrix <- getMatrixFromProject(proj_pass_filter, useMatrix = "PeakMatrix")
  PK_df <- as.data.frame(PK_matrix@rowRanges)
  PK_matrix@assays@data@listData[["PeakMatrix"]]@Dimnames[[1]] <- paste0(PK_df$seqnames, "_", PK_df$start, "_", PK_df$end)
  PK_matrix <- PK_matrix@assays@data@listData[["PeakMatrix"]]
  
  Peak_diff <- c()
  Gene_diff <- c()
  TF_diff <- c()
  Ref_alt_fc <- c()
  Gene_peak_cor <- c()
  TF_gene_cor <- c()
  
  cellData <- as.data.frame(proj_pass_filter@cellColData)
  Sys.time()
  for (i in 1:nrow(Celltype_Stat)) {
    GS_matrix_temp <- GS_matrix[Celltype_Stat[i, "SNP_nearestGene"], ]
    GS_matrix_temp <- as.data.frame(GS_matrix_temp)
    GS_matrix_temp$Celltype <- cellData[rownames(GS_matrix_temp), "Celltype"]
    GS_matrix_temp_mean <- aggregate.data.frame(GS_matrix_temp[,1],
                                                by = list(GS_matrix_temp$Celltype),
                                                FUN = mean)
    colnames(GS_matrix_temp_mean)[1] <- "Celltype"
    rownames(GS_matrix_temp_mean) <- GS_matrix_temp_mean$Celltype
    
    PK_matrix_temp <- as.data.frame(PK_matrix[Celltype_Stat[i, "Peak"],])
    PK_matrix_temp$Celltype <- cellData[rownames(PK_matrix_temp), "Celltype"]
    PK_matrix_temp <- PK_matrix_temp[which(PK_matrix_temp$Celltype %in% rownames(celltype_enrichment_matrix_log)),]
    PK_matrix_temp_mean <- aggregate.data.frame(PK_matrix_temp[,1],
                                                by = list(PK_matrix_temp$Celltype),
                                                FUN = mean)
    colnames(PK_matrix_temp_mean)[1] <- "Celltype"
    rownames(PK_matrix_temp_mean) <- PK_matrix_temp_mean$Celltype
    
    GS_PK_temp <- data.frame(GS_matrix_temp_mean,
                             PK_matrix_temp_mean[GS_matrix_temp_mean$Celltype,c(2)])
    cor_temp <- cor(GS_PK_temp[,2], GS_PK_temp[,3])
    Gene_peak_cor <- c(Gene_peak_cor, cor_temp)
    
    print(i)
    cell_markers <- markerPeaksList_Celltype[[Celltype_Stat$Celltype[i]]]
    cell_markers <- as.data.frame(cell_markers)
    cell_markers <- paste0(cell_markers$seqnames, "_", cell_markers$start, "_", cell_markers$end)
    
    tf <- Celltype_Stat[i,"TF"]
    tf_id <- Annotation_genes[tf, "Locus_ID"]
    tf_symbol <- gene_id_map[which(gene_id_map$gene_id == tf_id), "symbol"]
    TF_diff <- c(TF_diff, ifelse(tf_symbol %in% markerGenesList_celltype[[Celltype_Stat$Celltype[i]]][,"name"],
                                 1, 0))
    
    GS_matrix_temp_TF <- GS_matrix[tf_symbol, ]
    GS_matrix_temp_TF <- as.data.frame(GS_matrix_temp_TF)
    GS_matrix_temp_TF$Celltype <- cellData[rownames(GS_matrix_temp_TF), "Celltype"]
    GS_matrix_temp_TF_mean <- aggregate.data.frame(GS_matrix_temp_TF[,1],
                                                by = list(GS_matrix_temp_TF$Celltype),
                                                FUN = mean)
    colnames(GS_matrix_temp_TF_mean)[1] <- "Celltype"
    rownames(GS_matrix_temp_TF_mean) <- GS_matrix_temp_TF_mean$Celltype
    
    gene_TF_temp <- data.frame(GS_matrix_temp_mean,
                               GS_matrix_temp_TF_mean[GS_matrix_temp_mean$Celltype,c(2)])
    cor_temp <- cor(gene_TF_temp[,2], gene_TF_temp[,3])
    TF_gene_cor <- c(Gene_peak_cor, cor_temp)
    
    gene <- Celltype_Stat[i, "SNP_nearestGene"]
    gene_symbol <- gene
    Gene_diff <- c(Gene_diff, ifelse(gene_symbol %in% markerGenesList_celltype[[Celltype_Stat$Celltype[i]]][,"name"],
                                     1, 0))
    Peak_diff <- c(Peak_diff,
                   ifelse(Celltype_Stat$Peak[i] %in% cell_markers, 1, 0))
    
    fc_temp <- as.numeric(Celltype_Stat$alt_peak_motif_score[i]) / as.numeric(Celltype_Stat$ref_peak_motif_score[i])
    Ref_alt_fc <- c(Ref_alt_fc, fc_temp)
  }
  Sys.time()
  
  
  Celltype_Stat$Peak_diff <- Peak_diff
  Celltype_Stat$Gene_diff <- Gene_diff
  Celltype_Stat$TF_diff <- TF_diff
  Celltype_Stat$Ref_alt_fc <- Ref_alt_fc
  Celltype_Stat$Gene_peak_cor <- Gene_peak_cor
  # Celltype_Stat$TF_gene_cor <- TF_gene_cor
  
  Celltype_Stat_filtered <- Celltype_Stat[which(Celltype_Stat$Ref_alt_fc > 1),]
  Celltype_Stat_filtered <- Celltype_Stat_filtered[which(Celltype_Stat_filtered$Celltype %in% c("Cryptic bract/bract (cb/b)",
                                                                                                "Inflorescence meristem (IM)",
                                                                                                "Branch meristems (BM)",
                                                                                                "Rachis")),]
  Celltype_Stat_filtered <- Celltype_Stat_filtered[which(Celltype_Stat_filtered$Gene_peak_cor > 0.2),]
    
  Celltype_Stat_filtered <- Celltype_Stat[which(Celltype_Stat$Gene_diff == 1 & Celltype_Stat$TF_diff == 1),]
  
  openxlsx::write.xlsx(Celltype_Stat_filtered, "Figure7_Celltype_Stat_100bp_filtered.xlsx")
  
  
  
  
  
  
}

#### Peak range
{
  #### REF BSgenome 与 peak Motif 比对
  {
    motifPositions <- motifmatchr::matchMotifs(
      pwms = motif_all,
      subject = peaksets,
      genome = BSgenome.Oryza.sativa.IRGSP.1.0.REF, 
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
    overlapMotifs <- findOverlaps(peaksets, allPositions, ignore.strand = TRUE)
    overlapMotifs <- as.data.frame(overlapMotifs)
    peak_id <- paste0(peaksets_df$peak_id)
    allPositions$TF <- allPositions@ranges@NAMES
    allPositions@ranges@NAMES <- NULL
    allPositions_df <- as.data.frame(allPositions)
    ref_peak_TF <- data.frame(ref_peak_motif_chr = allPositions_df$seqnames[overlapMotifs$subjectHits],
                              ref_peak_motif_start = allPositions_df$start[overlapMotifs$subjectHits],
                              ref_peak_motif_end = allPositions_df$end[overlapMotifs$subjectHits],
                              ref_peak_TF = allPositions_df$TF[overlapMotifs$subjectHits],
                              ref_peak_motif_strand = allPositions_df$strand[overlapMotifs$subjectHits],
                              ref_peak_motif_score = allPositions_df$score[overlapMotifs$subjectHits],
                              ref_peak_motif_width = allPositions_df$width[overlapMotifs$subjectHits],
                              ref_peak_motif_peakid = peak_id[overlapMotifs$queryHits]
    )
    length(unique(ref_peak_TF$ref_peak_motif_peakid)) # ref peak 131982
    peaksets_df
    ref_peak_TF$ref_peak_gene <- peaksets_df[ref_peak_TF$ref_peak_motif_peakid, "nearestGene"]
    gene_annotation_df <- as.data.frame(gene_annotation@listData[["genes"]])
    rownames(gene_annotation_df) <- gene_annotation_df$symbol
    ref_peak_TF$ref_peak_gene_strand <- as.character(gene_annotation_df[ref_peak_TF$ref_peak_gene, "strand"])
    ref_peak_TF <- ref_peak_TF[which(as.character(ref_peak_TF$ref_peak_motif_strand) == ref_peak_TF$ref_peak_gene_strand),]
    
    ref_peak_TF_SNP <- merge(x = ref_peak_TF, y = SNP_8240_peak, by.x = "ref_peak_motif_peakid", by.y = "SNP_peak_id")
    ref_peak_TF_SNP <- ref_peak_TF_SNP[which(ref_peak_TF_SNP$SNP_loci >= ref_peak_TF_SNP$ref_peak_motif_start &
                                               ref_peak_TF_SNP$SNP_loci <= ref_peak_TF_SNP$ref_peak_motif_end),]
  }
  
  #### ALT BSgenome 与 peak Motif 比对
  {
    motifPositions <- motifmatchr::matchMotifs(
      pwms = motif_all,
      subject = peaksets,
      genome = BSgenome.Oryza.sativa.IRGSP.1.0.ALT, 
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
    overlapMotifs <- findOverlaps(peaksets, allPositions, ignore.strand = TRUE)
    overlapMotifs <- as.data.frame(overlapMotifs)
    peak_id <- paste0(peaksets_df$peak_id)
    allPositions$TF <- allPositions@ranges@NAMES
    allPositions@ranges@NAMES <- NULL
    allPositions_df <- as.data.frame(allPositions)
    alt_peak_TF <- data.frame(alt_peak_motif_chr = allPositions_df$seqnames[overlapMotifs$subjectHits],
                              alt_peak_motif_start = allPositions_df$start[overlapMotifs$subjectHits],
                              alt_peak_motif_end = allPositions_df$end[overlapMotifs$subjectHits],
                              alt_peak_TF = allPositions_df$TF[overlapMotifs$subjectHits],
                              alt_peak_motif_strand = allPositions_df$strand[overlapMotifs$subjectHits],
                              alt_peak_motif_score = allPositions_df$score[overlapMotifs$subjectHits],
                              alt_peak_motif_width = allPositions_df$width[overlapMotifs$subjectHits],
                              alt_peak_motif_peakid = peak_id[overlapMotifs$queryHits]
    )
    length(unique(alt_peak_TF$alt_peak_motif_peakid)) # alt peak 131983
    peaksets_df
    alt_peak_TF$alt_peak_gene <- peaksets_df[alt_peak_TF$alt_peak_motif_peakid, "nearestGene"]
    gene_annotation_df <- as.data.frame(gene_annotation@listData[["genes"]])
    rownames(gene_annotation_df) <- gene_annotation_df$symbol
    alt_peak_TF$alt_peak_gene_strand <- as.character(gene_annotation_df[alt_peak_TF$alt_peak_gene, "strand"])
    alt_peak_TF <- alt_peak_TF[which(as.character(alt_peak_TF$alt_peak_motif_strand) == alt_peak_TF$alt_peak_gene_strand),]
    
    alt_peak_TF_SNP <- merge(x = alt_peak_TF, y = SNP_8240_peak, by.x = "alt_peak_motif_peakid", by.y = "SNP_peak_id")
    alt_peak_TF_SNP <- alt_peak_TF_SNP[which(alt_peak_TF_SNP$SNP_loci >= alt_peak_TF_SNP$alt_peak_motif_start &
                                               alt_peak_TF_SNP$SNP_loci <= alt_peak_TF_SNP$alt_peak_motif_end),]
  }
  
  ### 整合 ref 和 alt 的结合结果
  {
    ref_peak_TF_SNP_2 <- ref_peak_TF_SNP
    ref_peak_TF_SNP_2 <- split.data.frame(ref_peak_TF_SNP_2, f = as.factor(paste0(ref_peak_TF_SNP_2$SNP,
                                                                                  "_",
                                                                                  ref_peak_TF_SNP_2$ref_peak_TF)))
    ref_peak_TF_SNP_2 <- lapply(ref_peak_TF_SNP_2, function(x){
      x[which.max(x$ref_peak_motif_score)[1],]
    })
    ref_peak_TF_SNP_2 <- rbindlist(ref_peak_TF_SNP_2)
    ref_peak_TF_SNP_2 <- apply(ref_peak_TF_SNP_2, 2, as.character)
    ref_peak_TF_SNP_2 <- as.data.frame(ref_peak_TF_SNP_2)
    ref_peak_TF_SNP_2$SNP_TF <- paste0(ref_peak_TF_SNP_2$SNP, "_", ref_peak_TF_SNP_2$ref_peak_TF)
    sum(duplicated(ref_peak_TF_SNP_2$SNP_TF)) # 0
    
    alt_peak_TF_SNP_2 <- alt_peak_TF_SNP
    alt_peak_TF_SNP_2 <- split.data.frame(alt_peak_TF_SNP_2, f = as.factor(paste0(alt_peak_TF_SNP_2$SNP,
                                                                                  "_",
                                                                                  alt_peak_TF_SNP_2$alt_peak_TF)))
    alt_peak_TF_SNP_2 <- lapply(alt_peak_TF_SNP_2, function(x){
      x[which.max(x$alt_peak_motif_score)[1],]
    })
    alt_peak_TF_SNP_2 <- rbindlist(alt_peak_TF_SNP_2)
    alt_peak_TF_SNP_2 <- apply(alt_peak_TF_SNP_2, 2, as.character)
    alt_peak_TF_SNP_2 <- as.data.frame(alt_peak_TF_SNP_2)
    alt_peak_TF_SNP_2$SNP_TF <- paste0(alt_peak_TF_SNP_2$SNP, "_", alt_peak_TF_SNP_2$alt_peak_TF)
    sum(duplicated(alt_peak_TF_SNP_2$SNP_TF)) # 0
    
    length(unique(ref_peak_TF_SNP_2$SNP))
    length(unique(alt_peak_TF_SNP_2$SNP))
    
    ref_peak_TF_SNP_2$TF <- ref_peak_TF_SNP_2$ref_peak_TF
    alt_peak_TF_SNP_2$TF <- alt_peak_TF_SNP_2$alt_peak_TF
    ref_peak_TF_SNP_2 <- ref_peak_TF_SNP_2[,-which(colnames(ref_peak_TF_SNP_2) == "SNP_trait")]
    alt_peak_TF_SNP_2 <- alt_peak_TF_SNP_2[,-which(colnames(alt_peak_TF_SNP_2) == "SNP_trait")]
    ref_alt_peak_TF_SNP <- merge(x = ref_peak_TF_SNP_2, y = alt_peak_TF_SNP_2,
                                 by = c("SNP_chr", "SNP_loci", "SNP", "SNP_REF", "SNP_ALT", "SNP_nearestGene", "TF"), all = TRUE)
    ref_alt_peak_TF_SNP[is.na(ref_alt_peak_TF_SNP)] <- ""
  }
  
  #### 对于每一个细胞类型和表型对，分别统计
  celltype_enrichment_matrix_binary_df <- reshape2::melt(celltype_enrichment_matrix_binary)
  celltype_enrichment_matrix_binary_df <- celltype_enrichment_matrix_binary_df[which(celltype_enrichment_matrix_binary_df$value == 1),]
  colnames(celltype_enrichment_matrix_binary_df) <- c("Celltype", "Trait", "SigIf")
  Celltype_Stat <- c()
  Celltype_Trait_SNP_remove <- c()
  Celltype_Trait_SNP_overlap <- c() # 和细胞类型有overlap
  Celltype_Trait_SNP_all <- c() # 显著富集的trait的所有SNP
  for (i in 1:nrow(celltype_enrichment_matrix_binary_df)) {
    print(i)
    trait <- celltype_enrichment_matrix_binary_df[i, 2]
    trait <- as.character(trait)
    trait_peaks <- Intergenic_trait_peaks[[trait]]
    trait_peaks <- unique(trait_peaks$peak_id)
    trait_SNP <- Intergenic_trait_peaks[[trait]]
    Celltype_Trait_SNP_all <- c(Celltype_Trait_SNP_all, unique(trait_SNP$SNP_CHR_BP))
    cell <- celltype_enrichment_matrix_binary_df[i, 1]
    cell <- as.character(cell)
    cell_markers <- markerPeaksList_Celltype[[cell]]
    cell_markers <- as.data.frame(cell_markers)
    cell_markers <- paste0(cell_markers$seqnames, "_", cell_markers$start, "_", cell_markers$end)
    common_peaks <- intersect(cell_markers, trait_peaks)
    for (peak in common_peaks) {
      trait_SNP <- Intergenic_trait_peaks[[trait]]
      trait_SNP <- trait_SNP[which(trait_SNP$peak_id == peak),]
      trait_SNP <- trait_SNP[!duplicated(trait_SNP$SNP_CHR_BP),]
      rownames(trait_SNP) <- trait_SNP$SNP_CHR_BP
      for (snp in trait_SNP$SNP_CHR_BP) { # 对每一个SNP循环
        temp_trait_SNP <- trait_SNP[snp,]
        Celltype_Trait_SNP_overlap <- as.data.frame(rbind(Celltype_Trait_SNP_overlap,
                                                          data.frame(Celltype = cell,
                                                                     Trait = trait,
                                                                     SNP = temp_trait_SNP$SNP_CHR_BP,
                                                                     Peaks = temp_trait_SNP$peak_id)))
        temp_ref_alt_peak_TF_SNP <- ref_alt_peak_TF_SNP[which(ref_alt_peak_TF_SNP$SNP %in% snp),]
        if (nrow(temp_ref_alt_peak_TF_SNP) == 0) {
          Celltype_Trait_SNP_remove <- as.data.frame(rbind(Celltype_Trait_SNP_remove,
                                                           data.frame(Celltype = cell,
                                                                      Trait = trait,
                                                                      SNP = temp_trait_SNP$SNP_CHR_BP,
                                                                      Peaks = temp_trait_SNP$peak_id)))
        } else {
          print("Yes!")
          Celltype_Stat_temp <- data.frame(Celltype = cell,
                                           Trait = trait,
                                           temp_ref_alt_peak_TF_SNP)
          
          Celltype_Stat <- as.data.frame(rbind(Celltype_Stat,
                                               Celltype_Stat_temp))
          
        }
      }
    }
  }
  Celltype_Stat$Peak <- overlapRegions_DF[Celltype_Stat$SNP, 2]
  for (i in 1:nrow(Celltype_Stat)) {
    if (Celltype_Stat$ref_peak_TF[i] == "" & Celltype_Stat$alt_peak_TF[i] != "") {
      Celltype_Stat$TF[i] <- Celltype_Stat$alt_peak_TF[i]
    }
    if (Celltype_Stat$ref_peak_TF[i] != "" & Celltype_Stat$alt_peak_TF[i] == "") {
      Celltype_Stat$TF[i] <- Celltype_Stat$ref_peak_TF[i]
    }
    if (Celltype_Stat$ref_peak_TF[i] != "" & Celltype_Stat$alt_peak_TF[i] != "") {
      if (Celltype_Stat$ref_peak_TF[i] == Celltype_Stat$alt_peak_TF[i]) {
        Celltype_Stat$TF[i] <- Celltype_Stat$ref_peak_TF[i]
      } else {
        Celltype_Stat$TF[i] <- ""
      }
      
    }
  }
  
  Celltype_Trait_SNP_overlap_IM <- Celltype_Trait_SNP_overlap[which(Celltype_Trait_SNP_overlap$Celltype %in% 
                                                                      c("Branch meristems (BM)", "Cryptic bract/bract (cb/b)",
                                                                        "Inflorescence meristem (IM)", "Rachis",
                                                                        "Spikelet meristem (SM)")),]
  Celltype_Trait_SNP_overlap_IM <- Celltype_Trait_SNP_overlap_IM[which(Celltype_Trait_SNP_overlap_IM$Trait %in% colnames(celltype_enrichment_matrix_log)[1:11]),]
  Celltype_Trait_SNP_overlap_IM$genes <- peaksets_df[Celltype_Trait_SNP_overlap_IM$Peaks, "nearestGene"]
  Celltype_Stat_ordered_IM_gene_list <- unique(Celltype_Trait_SNP_overlap_IM$genes)
  
  nrow(Celltype_Stat) # 746
  length(unique(Celltype_Stat$SNP)) # 279 SNPs
  length(unique(Celltype_Stat$Celltype)) # 17 Celltypes
  length(unique(Celltype_Stat$Trait)) # 26 Traits
  length(unique(Celltype_Stat$TF)) # 190 TFs

  
  SNP_TF_group <- c()
  for (i in 1:nrow(Celltype_Stat)) {
    temp_ref <- Celltype_Stat$ref_peak_motif_strand[i] == ""
    temp_alt <- Celltype_Stat$alt_peak_motif_strand[i] == ""
    if (temp_ref == TRUE & temp_alt == FALSE) { # 突变导致结合
      SNP_TF_group <- c(SNP_TF_group, "1")
    }
    if (temp_ref == FALSE & temp_alt == FALSE) { # 突变不影响结合
      SNP_TF_group <- c(SNP_TF_group, "2")
    }
    if (temp_ref == FALSE & temp_alt == TRUE) { # 突变导致不能结合
      SNP_TF_group <- c(SNP_TF_group, "3")
    }
  }
  table(SNP_TF_group)
  Celltype_Stat$SNP_TF_group <- SNP_TF_group
  Celltype_Stat$ref_peak_motif_score <- as.numeric(Celltype_Stat$ref_peak_motif_score)
  Celltype_Stat$alt_peak_motif_score <- as.numeric(Celltype_Stat$alt_peak_motif_score)
  which(Celltype_Stat$ref_peak_motif_score == 0)
  which(Celltype_Stat$alt_peak_motif_score == 0)
  Celltype_Stat$ref_peak_motif_score[which(is.na(Celltype_Stat$ref_peak_motif_score))] <- 0
  Celltype_Stat$alt_peak_motif_score[which(is.na(Celltype_Stat$alt_peak_motif_score))] <- 0
  
  Celltype_Stat$ref_alt_score_fc <- Celltype_Stat$ref_peak_motif_score / Celltype_Stat$alt_peak_motif_score
  
  Celltype_Stat_ordered <- c()
  
  temp <- Celltype_Stat[which(Celltype_Stat$SNP_TF_group == "1"),]
  temp <- temp[order(temp$alt_peak_motif_score, decreasing = T),]
  pheatmap(temp[,c("ref_peak_motif_score", "alt_peak_motif_score")],
           cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
  Celltype_Stat_ordered <- as.data.frame(rbind(Celltype_Stat_ordered,
                                               temp))
  
  temp <- Celltype_Stat[which(Celltype_Stat$SNP_TF_group == "2"),]
  temp <- temp[order(temp$ref_alt_score_fc, decreasing = F),]
  pheatmap(temp[,c("ref_peak_motif_score", "alt_peak_motif_score")],
           cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
  Celltype_Stat_ordered <- as.data.frame(rbind(Celltype_Stat_ordered,
                                               temp))
  
  temp <- Celltype_Stat[which(Celltype_Stat$SNP_TF_group == "3"),]
  temp <- temp[order(temp$ref_peak_motif_score, decreasing = T),]
  pheatmap(temp[,c("ref_peak_motif_score", "alt_peak_motif_score")],
           cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
  Celltype_Stat_ordered <- as.data.frame(rbind(Celltype_Stat_ordered,
                                               temp))
  
  Celltype_Stat_ordered$ref_peak_motif_score[which(Celltype_Stat_ordered$ref_peak_motif_score == 0)] <- NA
  Celltype_Stat_ordered$alt_peak_motif_score[which(Celltype_Stat_ordered$alt_peak_motif_score == 0)] <- NA
  Celltype_Stat_ordered$rownames <- paste0(Celltype_Stat_ordered$Celltype,
                                            "_",
                                            Celltype_Stat_ordered$Trait,
                                            "_",
                                            Celltype_Stat_ordered$SNP,
                                            "_",
                                            Celltype_Stat_ordered$TF)
  rownames(Celltype_Stat_ordered) <- Celltype_Stat_ordered$rownames
  temp <- Celltype_Stat_ordered[which(Celltype_Stat_ordered$rownames == "Branch meristems (BM)_grain_length_3_16728550_Os10g0115200"),]
  
  annotation_row <- Celltype_Stat_ordered[, c("Celltype","SNP_TF_group")]
  colnames(annotation_row) <- c("Celltype","SNP_TF_group")
  rownames(annotation_row) <- Celltype_Stat_ordered$rownames
  annotation_row$Celltype <- factor(annotation_row$Celltype,
                                    levels = sort(unique(annotation_row$Celltype)))
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
  length(celltype_color)
  names(celltype_color) <- sort(unique(proj_pass_filter$Celltype))
  ann_colors <- list(
    Celltype = celltype_color[levels(annotation_row$Celltype)],
    SNP_TF_group = c("1" = "slateblue3", "2" = "red2", "3" = "black")
  )
  
  paste0("Branch meristems (BM)",
         "_",
         "grain_length",
         "_",
         "3_16735437",
         "_",
         "OsOSH1")
  
  pdf("Figure7_Celltype_enriched_SNP_TF.pdf", width = 7, height = 10)
  pheatmap(Celltype_Stat_ordered[,c("ref_peak_motif_score", "alt_peak_motif_score")],
           cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F,
           annotation_row = annotation_row, annotation_colors = ann_colors,
           use_raster = TRUE, raster_resize_mat = TRUE, labels_row = c(12),
           color = rev(colorRampPalette(viridisLite::viridis(n = 20, option = "B"))(20)))
  Heatmap(Celltype_Stat_ordered[,c("ref_peak_motif_score", "alt_peak_motif_score")],
          cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
          use_raster = TRUE, raster_resize_mat = FALSE, col = rev(colorRampPalette(viridisLite::viridis(n = 20, option = "B"))(20))) +
    rowAnnotation(link = anno_mark(at = which(rownames(Celltype_Stat_ordered) %in% c("Branch meristems (BM)_grain_length_3_16735437_OsOSH1",
                                                                                     "Branch meristems (BM)_grain_length_3_16729065_OsRFL")), 
                                   labels = c("Branch meristems (BM)_grain_length_3_16735437_OsOSH1",
                                              "Branch meristems (BM)_grain_length_3_16729065_OsRFL"), labels_gp = gpar(fontsize = 10)))
  dev.off()
  
  
  
  ## IM
  Celltype_Stat_ordered_IM <- Celltype_Stat_ordered[which(Celltype_Stat_ordered$Celltype %in% c("Branch meristems (BM)", "Cryptic bract/bract (cb/b)",
                                                                                                "Inflorescence meristem (IM)", "Rachis",
                                                                                                "Spikelet meristem (SM)")),]
  Celltype_Stat_ordered_IM <- Celltype_Stat_ordered_IM[which(Celltype_Stat_ordered_IM$Trait %in% colnames(celltype_enrichment_matrix_log)[1:11]),]
  Celltype_Stat_ordered_IM_gene_list <- unique(Celltype_Stat_ordered_IM$SNP_nearestGene)
  rownames(gene_id_map) <- gene_id_map$symbol
  Celltype_Stat_ordered_IM_gene_list <- gene_id_map[Celltype_Stat_ordered_IM_gene_list, "gene_id"]
  genes_list <- list(IM = Celltype_Stat_ordered_IM_gene_list)
  # Trait Ontology
  TO_gene_list <- readRDS("TO_gene_list.rds")
  TO_gene_num <- table(TO_gene_list$TO)
  TO_gene_num <- TO_gene_num[which(TO_gene_num >= 5)]
  TO_gene_list <- TO_gene_list[which(TO_gene_list$TO %in% names(TO_gene_num)),]
  TO_gene_list[which(TO_gene_list$TO == "TO:0000190"), "Description"] <- "seed coat color"
  TO_result_Celltype_Stat_ordered_IM <- c()
  for (i in names(genes_list)) {
    print(i)
    genes <- genes_list[[i]]
    TO_result <- c()
    for (z in unique(TO_gene_list$TO)) {
      temp <- TO_gene_list[which(TO_gene_list$TO==z),]
      temp <- temp[!duplicated(temp$Gene),]
      gene_count <- intersect(genes,
                              temp$Gene)
      if (length(gene_count) > 0) {
        gene_ratio <- length(gene_count) / nrow(temp)
        ontology <- "Trait Ontology"
        description <- unique(temp$Description)
        x <- length(gene_count) - 1
        m <- length(unique(temp$Gene))
        n <- 37000 - m
        k <- length(genes)
        p <- phyper(x, m, n, k, lower.tail = F)
        fc <- (length(gene_count) / length(gene_count)) / (nrow(temp) / 37000)
        TO_result <- rbind.data.frame(TO_result,
                                      data.frame(TO_term = z,
                                                 Ontology = ontology,
                                                 TO_term_gene_count = nrow(temp),
                                                 Description = description,
                                                 Gene_count = length(gene_count),
                                                 Gene_ratio = gene_ratio,
                                                 P_value = p,
                                                 FC = fc,
                                                 Celltype_gene_count = length(genes)))
      } else {
        next
      }
    } # TO term enrichment
    TO_result_Celltype_Stat_ordered_IM <- as.data.frame(rbind(TO_result_Celltype_Stat_ordered_IM,
                                                        TO_result))
    
  }
  TO_result_Celltype_Stat_ordered_IM <- TO_result_Celltype_Stat_ordered_IM[which(TO_result_Celltype_Stat_ordered_IM$P_value < 0.05),]
  
  ### Gene Ontology
  GO_result <- c()
  for (i in names(genes_list)) {
    print(i)
    genes <- genes_list[[i]]
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
  }
  
  GO_result_Sig <- GO_result[which(GO_result$P_value < 0.05),]
  GO_result_Sig <- GO_result_Sig[which(GO_result_Sig$Ontology == "BP"),]
  
  annotation_row_IM <- Celltype_Stat_ordered_IM[, c("Celltype","SNP_TF_group")]
  colnames(annotation_row_IM) <- c("Celltype","SNP_TF_group")
  rownames(annotation_row_IM) <- Celltype_Stat_ordered_IM$rownames
  annotation_row_IM$Celltype <- factor(annotation_row_IM$Celltype,
                                       levels = sort(unique(annotation_row_IM$Celltype)))
  ann_colors_IM <- list(
    Celltype = celltype_color[levels(annotation_row_IM$Celltype)],
    SNP_TF_group = c("1" = "slateblue3", "2" = "red2", "3" = "black")
  )
  pdf("Figure7_Celltype_enriched_SNP_TF_IM.pdf", width = 7, height = 10)
  pheatmap(Celltype_Stat_ordered_IM[,c("ref_peak_motif_score", "alt_peak_motif_score")],
           cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F,
           annotation_row = annotation_row_IM, annotation_colors = ann_colors_IM,
           use_raster = TRUE, raster_resize_mat = TRUE, labels_row = c(12),
           color = rev(colorRampPalette(viridisLite::viridis(n = 20, option = "B"))(20)))
  Heatmap(Celltype_Stat_ordered_IM[,c("ref_peak_motif_score", "alt_peak_motif_score")],
          cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
          use_raster = TRUE, raster_resize_mat = FALSE, col = rev(colorRampPalette(viridisLite::viridis(n = 20, option = "B"))(20))) +
    rowAnnotation(link = anno_mark(at = which(rownames(Celltype_Stat_ordered_IM) %in% c("Branch meristems (BM)_grain_length_3_16735437_OsOSH1",
                                                                                     "Branch meristems (BM)_grain_length_3_16729065_OsRFL")), 
                                   labels = c("Branch meristems (BM)_grain_length_3_16735437_OsOSH1",
                                              "Branch meristems (BM)_grain_length_3_16729065_OsRFL"), labels_gp = gpar(fontsize = 10)))
  dev.off()
  
  markerGenesList_celltype
  table(Celltype_Stat$ref_peak_motif_strand, Celltype_Stat$alt_peak_motif_strand)
  
  openxlsx::write.xlsx(Celltype_Stat, "Figure7_Celltype_Stat_Peak.xlsx")
  
  ### ATAC diff
  Peak_diff <- c()
  Gene_diff <- c()
  TF_diff <- c()
  for (i in 1:nrow(Celltype_Stat)) {
    print(i)
    cell_markers <- markerPeaksList_Celltype[[Celltype_Stat$Celltype[i]]]
    cell_markers <- as.data.frame(cell_markers)
    cell_markers <- paste0(cell_markers$seqnames, "_", cell_markers$start, "_", cell_markers$end)
    
    tf <- Celltype_Stat[i,"TF"]
    tf_id <- Annotation_genes[tf, "Locus_ID"]
    tf_symbol <- gene_id_map[which(gene_id_map$gene_id == tf_id), "symbol"]
    TF_diff <- c(TF_diff, ifelse(tf_symbol %in% markerGenesList_celltype[[Celltype_Stat$Celltype[i]]][,"name"],
                                 1, 0))
    
    gene <- Celltype_Stat[i, "SNP_nearestGene"]
    gene_symbol <- gene
    Gene_diff <- c(Gene_diff, ifelse(gene_symbol %in% markerGenesList_celltype[[Celltype_Stat$Celltype[i]]][,"name"],
                                     1, 0))
    
    Peak_diff <- c(Peak_diff,
                   ifelse(Celltype_Stat$Peak[i] %in% cell_markers, 1, 0))
  }
  
  Celltype_Stat$Peak_diff <- Peak_diff
  Celltype_Stat$Gene_diff <- Gene_diff
  Celltype_Stat$TF_diff <- TF_diff
  
  Celltype_Stat_filtered <- Celltype_Stat[which(Celltype_Stat$Gene_diff == 1 & Celltype_Stat$TF_diff == 1),]
  Celltype_Stat_filtered <- Celltype_Stat_filtered[which(Celltype_Stat_filtered$Celltype == "Mesophyll initial"),]
  
  
  openxlsx::write.xlsx(Celltype_Stat_filtered, "Figure7_Celltype_Stat_Peak_filtered.xlsx")
  
  #### 细胞类型富集分析
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
  
  rownames(gene_id_map) <- gene_id_map$symbol
  # genes_list <- c()
  # for (i in unique(Celltype_Stat_ordered$Celltype)) {
  #   temp <- Celltype_Stat_ordered[which(Celltype_Stat_ordered$Celltype == i),]
  #   temp <- temp$SNP_nearestGene
  #   temp <- gene_id_map[temp, "gene_id"]
  #   genes_list <- c(genes_list, list(temp))
  # }
  # names(genes_list) <- unique(Celltype_Stat_ordered$Celltype)
  # genes_list <- genes_list[-17]
  genes_list <- c()
  markerGenesList_celltype <- getMarkers(markersGenes_celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  markerGenesList_celltype <- markerGenesList_celltype@listData
  marker_genes_celltypes <- c()
  for (i in names(markerGenesList_celltype)) {
    temp <- unique(markerGenesList_celltype[[i]]$name)
    temp <- gene_id_map[temp, "gene_id"]
    marker_genes_celltypes <- c(marker_genes_celltypes, list(temp))
  }
  names(marker_genes_celltypes) <- names(markerGenesList_celltype)
  marker_genes_celltypes <- marker_genes_celltypes[rownames(celltype_enrichment_matrix_log)]
  genes_list <- marker_genes_celltypes
  
  ### Gene Ontology
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
                  paste0("./Figure7_CellStat_",gsub(pattern = "/", replacement = "_", i, fixed = T),"_enrich_GO_result.tsv"))
      GO_result_Sig <- paste0("There is no GO pathways were significantly enriched.") 
      write.table(GO_result_Sig, quote = F, row.names = F, col.names = F, sep = "\t",
                  paste0("./Figure7_CellStat_",gsub(pattern = "/", replacement = "_", i, fixed = T),"_enrich_GO_result_Sig.tsv"))
    } else {
      GO_result$FDR <- p.adjust(GO_result$P_value, method = "fdr")
      GO_result <- GO_result[order(GO_result$Ontology, GO_result$P_value),]
      GO_result_Sig <- GO_result[which(GO_result$P_value<0.05),]
      if (nrow(GO_result_Sig) < 1) {
        write.table(GO_result, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure7_CellStat_",gsub(pattern = "/", replacement = "_", i, fixed = T),"_enrich_GO_result.tsv"))
        GO_result_Sig <- paste0("There is no GO pathways were significantly enriched.")
        write.table(GO_result_Sig, quote = F, row.names = F, col.names = F, sep = "\t",
                    paste0("./Figure7_CellStat_",gsub(pattern = "/", replacement = "_", i, fixed = T),"_enrich_GO_result_Sig.tsv"))
      } else {
        write.table(GO_result, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure7_CellStat_",gsub(pattern = "/", replacement = "_", i, fixed = T),"_enrich_GO_result.tsv"))
        write.table(GO_result_Sig, quote = F, row.names = F, col.names = T, sep = "\t",
                    paste0("./Figure7_CellStat_",gsub(pattern = "/", replacement = "_", i, fixed = T),"_enrich_GO_result_Sig.tsv"))
      }
    }
  }
  
  ### Trait Ontology
  TO_gene_list <- readRDS("TO_gene_list.rds")
  TO_gene_num <- table(TO_gene_list$TO)
  TO_gene_num <- TO_gene_num[which(TO_gene_num >= 5)]
  TO_gene_list <- TO_gene_list[which(TO_gene_list$TO %in% names(TO_gene_num)),]
  TO_gene_list[which(TO_gene_list$TO == "TO:0000190"), "Description"] <- "seed coat color"
  TO_result_celltype <- c()
  for (i in names(genes_list)) {
    print(i)
    genes <- genes_list[[i]]
    TO_result <- c()
    for (z in unique(TO_gene_list$TO)) {
      temp <- TO_gene_list[which(TO_gene_list$TO==z),]
      temp <- temp[!duplicated(temp$Gene),]
      gene_count <- intersect(genes,
                              temp$Gene)
      if (length(gene_count) > 0) {
        gene_ratio <- length(gene_count) / nrow(temp)
        ontology <- "Trait Ontology"
        description <- unique(temp$Description)
        x <- length(gene_count) - 1
        m <- length(unique(temp$Gene))
        n <- 37000 - m
        k <- length(genes)
        p <- phyper(x, m, n, k, lower.tail = F)
        fc <- (length(gene_count) / length(gene_count)) / (nrow(temp) / 37000)
        TO_result <- rbind.data.frame(TO_result,
                                      data.frame(TO_term = z,
                                                 Ontology = ontology,
                                                 TO_term_gene_count = nrow(temp),
                                                 Description = description,
                                                 Gene_count = length(gene_count),
                                                 Gene_ratio = gene_ratio,
                                                 P_value = p,
                                                 FC = fc,
                                                 Celltype_gene_count = length(genes)))
      } else {
        next
      }
    } # TO term enrichment
    
    TO_result$Celltype <- i
    TO_result_celltype <- as.data.frame(rbind(TO_result_celltype,
                                              TO_result))
    
  }
  TO_BM_Sig <- TO_result_celltype[which(TO_result_celltype$Celltype == "Branch meristems (BM)"),]
  TO_BM_Sig <- TO_BM_Sig[which(TO_BM_Sig$P_value < 0.05),]
  TO_BM_Sig <- TO_BM_Sig[order(TO_BM_Sig$Description),]
  traits_selected <- TO_BM_Sig$Description[24:31]
  TO_BM_Sig <- TO_BM_Sig[which(TO_BM_Sig$Description %in% traits_selected),]
  TO_BM_Sig <- TO_BM_Sig[order(TO_BM_Sig$P_value),]
  
  TO_result_celltype_selected <- TO_result_celltype[which(TO_result_celltype$Description %in% traits_selected),]
  TO_result_celltype_selected_matrx <- matrix(1, nrow = length(celltype_order), ncol = length(traits_selected))
  colnames(TO_result_celltype_selected_matrx) <- traits_selected
  rownames(TO_result_celltype_selected_matrx) <- celltype_order
  for (i in 1:nrow(TO_result_celltype_selected)) {
    TO_result_celltype_selected_matrx[TO_result_celltype_selected[i,"Celltype"],
                                      TO_result_celltype_selected[i,"Description"]] <- TO_result_celltype_selected[i,"P_value"]
  }
  
  if (!is.null(TO_result_celltype_selected_matrx)) {
    pmt <- TO_result_celltype_selected_matrx
    pmt <- as.data.frame(pmt)
    ssmt <- TO_result_celltype_selected_matrx <= 0.0001
    pmt[ssmt] <- "****"
    smt <- TO_result_celltype_selected_matrx > 0.0001 & TO_result_celltype_selected_matrx <= 0.001
    pmt[smt] <- "***"
    smt2 <- TO_result_celltype_selected_matrx > 0.001 & TO_result_celltype_selected_matrx <= 0.01
    pmt[smt2] <- "**"
    smt3 <- TO_result_celltype_selected_matrx > 0.01 & TO_result_celltype_selected_matrx <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  
  TO_result_celltype_selected_matrx_log <- -log10(TO_result_celltype_selected_matrx)
  TO_result_celltype_selected_matrx_log[which(TO_result_celltype_selected_matrx_log == -log10(1))] <- NA
  pdf("Figure7_TO_result_celltype_selected_matrx_log.pdf", height = 10, width = 12)
  pheatmap(TO_result_celltype_selected_matrx_log,
           color = colorRampPalette(c("#7AC5CD", "yellow", "#FA8072", "red"))(20),
           scale = "none", cluster_row = F,
           cluster_col = F,
           treeheight_col = 0,
           display_numbers = pmt,
           fontsize_number = 14,
           number_color = "black",
           cellwidth = 30, cellheight = 20)
  dev.off()
  
  
  TO_result_celltype_selected$Description <- factor(TO_result_celltype_selected$Description,
                                                    levels = rev(TO_BM_Sig$Description))
  TO_result_celltype_selected$Celltype <- factor(TO_result_celltype_selected$Celltype,
                                                 levels = celltype_order)
  TO_result_celltype_selected$P_value_log10 <- -log10(TO_result_celltype_selected$P_value)
  pdf("Figure7_TO_result_celltype_selected.pdf", width = 10, height = 8)
  ggplot(data = TO_result_celltype_selected, aes(x = P_value_log10, y = Description, fill = Celltype)) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(.~Celltype, nrow = 4) +
    geom_vline(xintercept = -log10(0.05), color = "red") +
    scale_fill_manual(values = celltype_color) +
    labs(x = "-Log10(P-Value)", y = "TO terms") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black"),
          panel.grid = element_blank()) +
    NoLegend()
  dev.off()
}




##########
##########
##########

#### 可视化TF 和 SNP
{
  library(BSgenome.Oryza.sativa.IRGSP.1.0.ALT)
  library(BSgenome.Oryza.sativa.IRGSP.1.0.REF)
  motif_pfm <- read_meme("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/Osj_TF_binding_motifs.meme")
  motif_pfm <- convert_motifs(motif_pfm, class = "TFBSTools-PFMatrix")
  motif_pfm <- do.call(PFMatrixList, motif_pfm)
  motif_pfm_MSU_ID <- unlist(lapply(motif_pfm@listData, function(x){
    x@name
  }))
  names(motif_pfm) <- motif_pfm_MSU_ID
  motif_pfm_RAP_ID <- match(motif_pfm_MSU_ID, MSU_RAP$MSU)
  motif_pfm_RAP_ID <- MSU_RAP$RAP[motif_pfm_RAP_ID]
  sum(is.na(motif_pfm_RAP_ID))
  motif_pfm <- motif_pfm[-which(is.na(motif_pfm_RAP_ID))]
  motif_pfm_RAP_ID <- motif_pfm_RAP_ID[-which(is.na(motif_pfm_RAP_ID))]
  names(motif_pfm) <- motif_pfm_RAP_ID
  motif_pfm <- motif_pfm[!duplicated(names(motif_pfm))] # ID duplicated !!!
  rownames(gene_id_map) <- gene_id_map$gene_id
  motif_pfm <- motif_pfm[names(motif_pfm) %in% gene_id_map$gene_id]
  Symbol <- match(names(motif_pfm), gene_id_map$gene_id)
  Symbol <- gene_id_map[Symbol,]
  rownames(Annotation_genes) <- Annotation_genes$Locus_ID
  names(motif_pfm) <- Annotation_genes[names(motif_pfm), "CGSNL Gene Symbol"]
  sum(is.na(names(motif_pfm)))
  
  pdf("Figure7_ggseqlogo_OsRFL.pdf", width = 6, height = 4)
  ggseqlogo::ggseqlogo(list(OsRFL = motif_all_pfm@listData$OsRFL@profileMatrix))
  dev.off()
  pdf("Figure3_ggseqlogo_OsMADS2.pdf", width = 6, height = 4)
  ggseqlogo::ggseqlogo(list(OsMADS2 = motif_all_pfm@listData$OsMADS2@profileMatrix))
  dev.off()
  
  ggseqlogo::ggseqlogo(list(OsERF53 = motif_all_pfm@listData$OsERF53@profileMatrix))
  getSeq(BSgenome.Oryza.sativa.IRGSP.1.0.ALT, names = "10", start = 19768190, end = 19768210, strand = "-")
  getSeq(BSgenome.Oryza.sativa.IRGSP.1.0.REF, names = "10", start = 19768190, end = 19768210, strand = "-")
  
  peakSet <- getPeakSet(proj_pass_filter)
  Chr <- as.data.frame(seqnames(peakSet))
  Range <- as.data.frame(peakSet@ranges)
  peakSet$Peak_id <- paste0(Chr$value, "_", Range$start, "_", Range$end)
  
  peak_group_selected <- peak_group[which(peak_group$nearestGene %in% genes_DF$gene),]
  temp <- peakSet[which(peakSet$Peak_id == "3_16735106_16735606")]
  
  featureList <- SimpleList("3_16735106_16735606" = temp)
  
  
  p_Maincluster <- plotBrowserTrack(ArchRProj = proj_pass_filter,
                                    groupBy = "Celltype",
                                    geneSymbol = "GS3",
                                    upstream = 15000, downstream = 15000,
                                    features = featureList,
                                    # peakPal = peak_type_color,
                                    # groupPal = Maincluster_color,
                                    facetbaseSize = 10)
  pdf("Figure7_GS3_peak_3_16735106_16735606.pdf")
  grid::grid.newpage()
  grid::grid.draw(p_Maincluster$GS3)
  dev.off()
  
  
  pdf("Figure7_ggseqlogo_OsOSH1.pdf", width = 6, height = 4)
  ggseqlogo::ggseqlogo(list(OsOSH1 = motif_all_pfm@listData$OsOSH1@profileMatrix))
  dev.off()
  
  pdf("Figure7_ggseqlogo_OsRFL.pdf", width = 6, height = 4)
  ggseqlogo::ggseqlogo(list(OsRFL = motif_all_pfm@listData$OsRFL@profileMatrix))
  dev.off()
  
  pdf("Figure7_叶肉起始细胞.pdf")
  ggseqlogo::ggseqlogo(list(OsERF54 = motif_all_pfm@listData$OsERF54@profileMatrix))
  dev.off()
  
  getSeq(BSgenome.Oryza.sativa.IRGSP.1.0.REF, names = "10", start = 19768190, end = 19768210, strand = "-")
  getSeq(BSgenome.Oryza.sativa.IRGSP.1.0.ALT, names = "10", start = 19768190, end = 19768210, strand = "-")
  
  
  temp <- peakSet[which(peakSet$Peak_id == "3_16728824_16729324")]
  featureList <- SimpleList("3_16728824_16729324" = temp)
  p_Maincluster <- plotBrowserTrack(ArchRProj = proj_pass_filter,
                                    groupBy = "Celltype",
                                    geneSymbol = "GS3",
                                    upstream = 15000, downstream = 15000,
                                    features = featureList,
                                    # peakPal = peak_type_color,
                                    # groupPal = Maincluster_color,
                                    facetbaseSize = 10)
  pdf("Figure7_GS3_peak_3_16728824_16729324.pdf")
  grid::grid.newpage()
  grid::grid.draw(p_Maincluster$GS3)
  dev.off()
  
  getSeq(BSgenome.Oryza.sativa.IRGSP.1.0.REF, names = "3", start = 16729061, end = 16729079, strand = "-")
  getSeq(BSgenome.Oryza.sativa.IRGSP.1.0.ALT, names = "3", start = 16729061, end = 16729079, strand = "-")
  
  
  GS_matrix <- GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]
  GS_matrix <- GS_matrix[, which(GeneScoreMatrix@colData@listData[["Celltype"]] %in% rownames(celltype_enrichment_matrix_log))]
  c(Annotation_genes["OsERF53","DataSets_Symbol"],
    Annotation_genes["OsERF54","DataSets_Symbol"],
    "OsRFPHC-10") %in% rownames(GS_matrix)
  GS_matrix_temp <- GS_matrix[c(Annotation_genes["OsERF53","DataSets_Symbol"],
                                Annotation_genes["OsERF54","DataSets_Symbol"], # OsERF54, OsERF.3
                                "OsRFPHC-10"), ]
  GS_matrix_temp <- as.data.frame(t(GS_matrix_temp))
  cellData <- as.data.frame(proj_pass_filter@cellColData)
  GS_matrix_temp$Celltype <- cellData[rownames(GS_matrix_temp), "Celltype"]
  length(unique(GS_matrix_temp$Celltype))
  GS_matrix_temp_mean <- aggregate.data.frame(GS_matrix_temp[,1:3],
                                              by = list(GS_matrix_temp$Celltype),
                                              FUN = mean)
  colnames(GS_matrix_temp_mean)[1] <- "Celltype"
  rownames(GS_matrix_temp_mean) <- GS_matrix_temp_mean$Celltype
  ggplot() +
    geom_point(data = GS_matrix_temp_mean, aes(y = `OsRFPHC-10`, x = OsERF.3, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_matrix_temp_mean, aes(y = `OsRFPHC-10`, x = OsERF.3), method = "pearson") +
    geom_smooth(data = GS_matrix_temp_mean, aes(y = `OsRFPHC-10`, x = OsERF.3), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) +
    labs(x = "Average accessibility of OsERF54", y = "Average gene score\nof OsSIP13") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  
  
  
  PK_matrix <- getMatrixFromProject(proj_pass_filter, useMatrix = "PeakMatrix")
  PK_df <- as.data.frame(PK_matrix@rowRanges)
  PK_matrix@assays@data@listData[["PeakMatrix"]]@Dimnames[[1]] <- paste0(PK_df$seqnames, "_", PK_df$start, "_", PK_df$end)
  PK_matrix <- PK_matrix@assays@data@listData[["PeakMatrix"]]
  PK_matrix_temp <- as.data.frame(t(PK_matrix[c("3_16728824_16729324", "3_16735106_16735606","10_19767946_19768446"),]))
  PK_matrix_temp$Celltype <- cellData[rownames(PK_matrix_temp), "Celltype"]
  PK_matrix_temp <- PK_matrix_temp[which(PK_matrix_temp$Celltype %in% rownames(celltype_enrichment_matrix_log)),]
  PK_matrix_temp_mean <- aggregate.data.frame(PK_matrix_temp[,1:3],
                                              by = list(PK_matrix_temp$Celltype),
                                              FUN = mean)
  colnames(PK_matrix_temp_mean)[1] <- "Celltype"
  rownames(PK_matrix_temp_mean) <- PK_matrix_temp_mean$Celltype
  
  GS_PK_temp <- data.frame(GS_matrix_temp_mean,
                           PK_matrix_temp_mean[GS_matrix_temp_mean$Celltype,c(2,3,4)])
  
  pdf("Figure7_GS3_3_16728824_16729324_point_plot.pdf", width = 3.2, height = 2.2)
  ggplot() +
    geom_point(data = GS_PK_temp, aes(y = GS3, x = X3_16728824_16729324, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_PK_temp, aes(y = GS3, x = X3_16728824_16729324), ) +
    geom_smooth(data = GS_PK_temp, aes(y = GS3, x = X3_16728824_16729324), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) +
    labs(x = "Average accessibility of peak\n(chr3:16728824-16729324)", y = "Average gene score\nof GS3") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  dev.off()
  
  pdf("Figure7_GS3_3_16735106_16735606_point_plot.pdf", width = 3.2, height = 2.2)
  ggplot() +
    geom_point(data = GS_PK_temp, aes(y = GS3, x = X3_16735106_16735606, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_PK_temp, aes(y = GS3, x = X3_16735106_16735606), ) +
    geom_smooth(data = GS_PK_temp, aes(y = GS3, x = X3_16735106_16735606), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) +
    labs(x = "Average accessibility of peak\n(chr3:16735106-16735606)", y = "Average gene score\nof GS3") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  dev.off()
  
  pdf("Figure7_GS3_Os10g0115200_point_plot.pdf", width = 3.2, height = 2.1)
  ggplot() +
    geom_point(data = GS_PK_temp, aes(y = GS3, x = Os10g0115200, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_PK_temp, aes(y = GS3, x = Os10g0115200), method = "pearson") +
    geom_smooth(data = GS_PK_temp, aes(y = GS3, x = Os10g0115200), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) +
    labs(x = "Average accessibility of Os10g0115200", y = "Average gene score\nof GS3") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  dev.off()
  
  pdf("Figure7_GS3_RFL_point_plot.pdf", width = 3.2, height = 2.1)
  ggplot() +
    geom_point(data = GS_PK_temp, aes(y = GS3, x = RFL, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_PK_temp, aes(y = GS3, x = RFL), ) +
    geom_smooth(data = GS_PK_temp, aes(y = GS3, x = RFL), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) +
    labs(x = "Average accessibility of RFL", y = "Average gene score\nof GS3") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  dev.off()
  
  
  pdf("Figure7_OsSIP13_OsERF54_point_plot_叶肉.pdf", width = 3.2, height = 2.1)
  ggplot() +
    geom_point(data = GS_PK_temp, aes(y = OsRFPHC.10, x = OsERF.3, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_PK_temp, aes(y = OsRFPHC.10, x = OsERF.3), ) +
    geom_smooth(data = GS_PK_temp, aes(y = OsRFPHC.10, x = OsERF.3), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) +
    labs(x = "Average expression of OsERF54\n", y = "Average expression \nof OsSIP13") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  dev.off()
  
  
  pdf("Figure7_OsSIP13_10_19767946_19768446_point_plot_叶肉.pdf", width = 3.2, height = 2.1)
  ggplot() +
    geom_point(data = GS_PK_temp, aes(y = OsRFPHC.10, x = X10_19767946_19768446, color = Celltype), size = 1.5) +
    theme_bw() +
    stat_cor(data = GS_PK_temp, aes(y = OsRFPHC.10, x = X10_19767946_19768446), ) +
    geom_smooth(data = GS_PK_temp, aes(y = OsRFPHC.10, x = X10_19767946_19768446), method = "lm") +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.ticks = element_line(color = "black")) + 
    labs(x = "Average accessibility of peak\n(chr10:19767946-19768446)", y = "Average gene score\nof OsSIP13") +
    scale_color_manual(values = celltype_color) +
    NoLegend()
  dev.off()
}


#################
sample_plink <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_selected_SNP_2.ped",
                         sep = " ", header = F)
sample_plink <- read_plink("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/529samples/rice4k_sample_filtered_selected_SNP")
sample_plink <- sample_plink[["X"]]
sample_plink <- as.data.frame(t(sample_plink))
sample_pheno <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/phenos.csv",
                         sep = ",", header = T)
rownames(sample_pheno) <- sample_pheno$FID
sample_plink$Grain_weight <- sample_pheno[rownames(sample_plink), "Grain_weight"]
sample_plink$Grain_length <- sample_pheno[rownames(sample_plink), "Grain_length"]
sample_plink <- sample_plink[which(sample_plink$Grain_length != -9),]
sample_plink <- sample_plink[which(sample_plink$Grain_weight != -9),]

sample_plink$vg0316729065 <- as.character(sample_plink$vg0316729065)
sample_plink$vg0316735437 <- as.character(sample_plink$vg0316735437)
sample_plink$vg1019768203 <- as.character(sample_plink$vg1019768203)

ggplot(sample_plink, aes(x = vg0316735437, y = Grain_length)) +
  geom_boxplot() +
  ggpubr::stat_compare_means()

ggplot(sample_plink, aes(x = vg0316729065, y = Grain_length)) +
  geom_boxplot() +
  ggpubr::stat_compare_means()

ggplot(sample_plink, aes(x = vg1019768203, y = Grain_weight)) +
  geom_boxplot() +
  ggpubr::stat_compare_means()

RiceVarMap2.0_SNP[which(RiceVarMap2.0_SNP$SNP == "RiceVarMap2.0_10_19768203"),]

#########################
#########################
#########################
#########################
#########################
message("以下为废弃代码！")
message("以下为废弃代码！")
message("以下为废弃代码！")
message("以下为废弃代码！")
message("以下为废弃代码！")
message("以下为废弃代码！")
message("以下为废弃代码！")

if (FALSE) {
  # SNPs of trait annotated more than 5 genes
  
  
  temp <- as.data.frame(table(All_SNP_position_GR_nearstgene_DF$TRAIT_Renamed,
                              All_SNP_position_GR_nearstgene_DF$nearestGene))
  colnames(temp) <- c("Trait", "SNP_Count")
  pdf("")
  ggplot(data = temp, aes(x = log2(SNP_Count))) +
    geom_density() +
    scale_x_continuous(limits = c(0, 100)) +
    geom_vline(xintercept = log2(5), color = "red")
  dev.off()
  
  # gene annotation
  {
    gene_annotation_df <- as.data.frame(gene_annotation)
    gene_annotation_df <- gene_annotation_df[gene_annotation_df$group_name == "genes",]
    gene_annotation_MAGMA <- data.frame(id = gene_annotation_df$gene_id,
                                        chr = gene_annotation_df$seqnames,
                                        start = gene_annotation_df$start,
                                        end = gene_annotation_df$end,
                                        strand = gene_annotation_df$strand,
                                        gene = gene_annotation_df$symbol)
    write.table(gene_annotation_MAGMA,
                "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/MAGMA/rice.gene.loc", sep = "\t", quote = F,
                row.names = F, col.names = F)
  }
  # gene level analysis
  {
    MAGMA_path <- "/media/heshidian/RAID5_42TB/3.Software/magma_v1.10/magma"
    setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA")
    for (i in meta_name) {
      print(i)
      path <- "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/"
      path <- paste0(path, i)
      setwd(path)
      # SNP annotation
      command <- paste0(MAGMA_path, " --annotate window=1,1 ",
                        "--snp-loc rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final.bim ",
                        "--gene-loc /media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/MAGMA/rice.gene.loc ",
                        "--out MAGMA_SNP_anno")
      system(command = command)
      # 给fam文件更新表型
      command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                        "-bfile ", "rice4k_sample_filtered_SNP_filtered_MAF_filtered_indep_final ",
                        "--pheno pheno_plink.txt --allow-no-sex --make-bed --out input_with_pheno "
      )
      system(command)
      command <- paste0("/home/heshidian/mambaforge/envs/common/bin/plink ",
                        "-bfile ", "input_with_pheno ",
                        "--recode --allow-no-sex --out input_with_pheno"
      )
      system(command)
      
      # gene analysis
      SNP <- data.table::fread(paste0(path, "/output/result.assoc.txt"),
                               sep = "\t", header = T)
      SNP <- SNP[,c("rs", "p_wald")]
      write.table(SNP, paste0(path, "/SNP_P.txt"),
                  sep = " ", row.names = F, col.names = F, quote = F)
      sample_number <- read.csv("pheno_plink.txt", sep = "\t", header = T)
      sample_number <- nrow(sample_number)
      command <- paste0(MAGMA_path, " --bfile input_with_pheno ",
                        "--pval SNP_P.txt N=", sample_number, " ",
                        "--gene-annot MAGMA_SNP_anno.genes.annot ",
                        "--out MAGMA_gene_analysis")
      system(command)
      
    }
    
    RNA <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Rice-snRNA/all_data_annotated.rds")
    Idents(RNA) <- RNA$Celltype
    RNA_celltype_marker <- FindAllMarkers(RNA, only.pos = T)
    gene_anno <- as.data.frame(gene_annotation)
    gene_anno <- gene_anno[which(gene_anno$group_name == "genes"),]
    rownames(gene_anno) <- gene_anno$gene_id
    chr_length <- as.data.frame(genome_annotation@listData[["chromSizes"]])
    rownames(chr_length) <- chr_length$seqnames
    
    for (cell in unique(as.character(RNA_celltype_marker$cluster))) {
      temp_celltype_marker <- RNA_celltype_marker[which(RNA_celltype_marker$cluster == cell),]
      temp_celltype_marker <- temp_celltype_marker[which(temp_celltype_marker$avg_log2FC >= 1),]
      temp_celltype_marker <- temp_celltype_marker[which(temp_celltype_marker$p_val_adj < 0.05),]
      temp_celltype_marker <- as.data.frame(na.omit(temp_celltype_marker))
      temp_bed <- gene_anno[temp_celltype_marker$gene,]
      temp_bed <- temp_bed[, c("seqnames", "start", "end", "gene_id", "strand")]
      temp_bed <- as.data.frame(na.omit(temp_bed))
      celltype <- gsub(pattern = " ", replacement = "_", cell)
      celltype <- gsub(pattern = "/", replacement = "_", celltype, fixed = T)
      temp_bed$celltype <- celltype
      temp_bed$start <- temp_bed$start - 100000
      temp_bed$start[which(temp_bed$start <= 0)] <- 1
      temp_bed$end <- temp_bed$end + 100000
      max_length <- chr_length[temp_bed$seqnames, "width"]
      temp_bed$end[temp_bed$end >= max_length] <- max_length[temp_bed$end >= max_length]
      write.table(temp_bed, paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/LDSC/Celltype_marker_bed/", celltype, "_marker_bed.txt"),
                  sep = "\t", quote = F, row.names = F, col.names = F)
    }
    
    Idents(RNA) <- RNA$tissues
    RNA_tissue_marker <- FindAllMarkers(RNA, only.pos = T)
    
    celltype_enrichment_matrix <- matrix(1, nrow = length(unique(RNA_celltype_marker$cluster)),
                                         ncol = length(meta_name))
    rownames(celltype_enrichment_matrix) <- unique(RNA_celltype_marker$cluster)
    colnames(celltype_enrichment_matrix) <- meta_name
    
    tissue_enrichment_matrix <- matrix(1, nrow = length(unique(RNA_tissue_marker$cluster)),
                                       ncol = length(meta_name))
    rownames(tissue_enrichment_matrix) <- unique(RNA_tissue_marker$cluster)
    colnames(tissue_enrichment_matrix) <- meta_name
    
    for (i in meta_name) {
      path <- "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/"
      path <- paste0(path, i)
      MAGMA <- read.csv(paste0(path, "/MAGMA_gene_analysis.genes.out"), sep = "\t", header = T)
      MAGMA2 <- apply(MAGMA, 1, function(x){
        temp <- unlist(strsplit(x, split = " "))
        temp <- temp[-which(temp == "")]
        temp
      })
      MAGMA2 <- as.data.frame(t(MAGMA2))
      colnames(MAGMA2) <- c("GENE", "CHR", "START", "STOP",
                            "NSNPS", 'NPARAM', "N", "ZSTAT", "P")
      MAGMA2$ZSTAT <- as.numeric(MAGMA2$ZSTAT)
      MAGMA2$P <- as.numeric(MAGMA2$P)
      rownames(MAGMA2) <- MAGMA2$GENE
      Trait_genes <- MAGMA2$GENE[which(MAGMA2$P < 0.05)]
      for (cell in unique(as.character(RNA_celltype_marker$cluster))) {
        temp_celltype_marker <- RNA_celltype_marker[which(RNA_celltype_marker$cluster == cell),]
        temp_celltype_marker <- temp_celltype_marker[which(temp_celltype_marker$avg_log2FC >= 1),]
        p <- enrichment_p(pathway_genes = Trait_genes, diff_genes = temp_celltype_marker$gene)
        celltype_enrichment_matrix[cell, i] <- p
      }
      for (tissue in unique(as.character(RNA_tissue_marker$cluster))) {
        temp_tissue_marker <- RNA_tissue_marker[which(RNA_tissue_marker$cluster == tissue),]
        temp_tissue_marker <- temp_tissue_marker[which(temp_tissue_marker$avg_log2FC >= 1),]
        p <- enrichment_p(pathway_genes = Trait_genes, diff_genes = temp_tissue_marker$gene)
        tissue_enrichment_matrix[tissue, i] <- p
      }
    }
    
    library(ggplot2)
    library(tidyverse)
    library(reshape2)
    library(ggtree) # 聚类
    library(aplot) # 拼图
    p_sig_star <- function(x) {
      star <- c()
      for (i in x) {
        if (i <= 0.0001) {
          star <- c(star, "****")
        }
        if (i > 0.0001 & i <= 0.001) {
          star <- c(star, "***")
        }
        if (i > 0.001 & i <= 0.01) {
          star <- c(star, "**")
        }
        if (i > 0.01 & i <= 0.05) {
          star <- c(star, "*")
        }
        if (i > 0.05) {
          star <- c(star, "")
        }
      }
      star
    }
    
    celltype_enrichment_matrix_2 <- celltype_enrichment_matrix[-which(rownames(celltype_enrichment_matrix) == "Unidentified"),]
    if (!is.null(celltype_enrichment_matrix_2)) {
      pmt <- celltype_enrichment_matrix_2
      pmt <- as.data.frame(pmt)
      ssmt <- celltype_enrichment_matrix_2 <= 0.0001
      pmt[ssmt] <- "****"
      smt <- celltype_enrichment_matrix_2 > 0.0001 & celltype_enrichment_matrix_2 <= 0.001
      pmt[smt] <- "***"
      smt2 <- celltype_enrichment_matrix_2 > 0.001 & celltype_enrichment_matrix_2 <= 0.01
      pmt[smt2] <- "**"
      smt3 <- celltype_enrichment_matrix_2 > 0.01 & celltype_enrichment_matrix_2 <= 0.05
      pmt[smt3] <- "*"
      pmt[!ssmt&!smt&!smt2&!smt3] <- ''
    }
    
    pdf("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Figure7_celltype_enrichment.pdf", width = 10, height = 10)
    pheatmap(celltype_enrichment_matrix_2,
             color = colorRampPalette(c("#B83D3D", "white", "#1A5592"))(100),
             scale = "none", cluster_row = T,
             cluster_col = T,                
             display_numbers = pmt,
             fontsize_number = 12,
             number_color = "black",
             cellwidth = 20, cellheight = 16)
    dev.off()
    
    if (!is.null(tissue_enrichment_matrix)) {
      pmt <- tissue_enrichment_matrix
      pmt <- as.data.frame(pmt)
      ssmt <- tissue_enrichment_matrix <= 0.0001
      pmt[ssmt] <- "****"
      smt <- tissue_enrichment_matrix > 0.0001 & tissue_enrichment_matrix <= 0.001
      pmt[smt] <- "***"
      smt2 <- tissue_enrichment_matrix > 0.001 & tissue_enrichment_matrix <= 0.01
      pmt[smt2] <- "**"
      smt3 <- tissue_enrichment_matrix > 0.01 & tissue_enrichment_matrix <= 0.05
      pmt[smt3] <- "*"
      pmt[!ssmt&!smt&!smt2&!smt3] <- ''
    }
    
    pdf("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Figure7_tissue_enrichment.pdf", width = 10, height = 10)
    pheatmap(tissue_enrichment_matrix,
             color = colorRampPalette(c("#B83D3D", "white", "#1A5592"))(100),
             scale = "none", cluster_row = T,
             cluster_col = T,                
             display_numbers = pmt,
             fontsize_number = 12,
             number_color = "black",
             cellwidth = 20, cellheight = 16)
    dev.off()
    
    pheno <- "Heading_date Grain_thickness Grain_width Grain_length Spikelet_length Grain_weight Yield Num_effective_panicles Num_panicles"
    Celltypes <- "Branch_meristems(BM) Columella_initials Cortex Cryptic_bract_bract(cb_b) Dividing_inner Dividing_outer Endodermis-cortical_initial Endodermis Epidermal_cell Epidermis Fiber Inflorescence_meristem(IM) Intermediate_anther Large_parenchyma_(MO) Late_anther Lemma(le) Mesophyll_initial Mesophyll_(MO) Mesophyll_precursor Other_Spikelet Phloem Procambium Proliferating_cell Rachis Shoot_endodermis Shoot_meristematic_cell Shoot_meristem_initials Spikelet_meristem(SM) Suspensor Upper_protoderm Vascular_initial"
    
    pheno <- unlist(strsplit(pheno, " "))
    Celltypes <- unlist(strsplit(Celltypes, " "))
    celltype_matrix <- matrix(1, nrow = length(pheno), ncol = length(Celltypes))
    rownames(celltype_matrix) <- pheno
    colnames(celltype_matrix) <- Celltypes
    enrich_matrix <- matrix(1, nrow = length(pheno), ncol = length(Celltypes))
    rownames(enrich_matrix) <- pheno
    colnames(enrich_matrix) <- Celltypes
    for (cell in Celltypes) {
      for (phe in pheno) {
        temp <- read.csv(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/LDSC/LDSC_result/",
                                phe, "_", cell, "_baselineLD.results"), sep = "\t", header = T)
        enrich_matrix[phe, cell] <- temp$Enrichment
        celltype_matrix[phe, cell] <- temp$Enrichment_p
      }
    }
    
    
    library(clusterProfiler)
    P_list <- c()
    NES_list <- c()
    for (i in meta_name) {
      path <- "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/RiceVarMap2.0/GEMMA/"
      path <- paste0(path, i)
      MAGMA <- read.csv(paste0(path, "/MAGMA_gene_analysis.genes.out"), sep = "\t", header = T)
      MAGMA2 <- apply(MAGMA, 1, function(x){
        temp <- unlist(strsplit(x, split = " "))
        temp <- temp[-which(temp == "")]
        temp
      })
      MAGMA2 <- as.data.frame(t(MAGMA2))
      colnames(MAGMA2) <- c("GENE", "CHR", "START", "STOP",
                            "NSNPS", 'NPARAM', "N", "ZSTAT", "P")
      MAGMA2$ZSTAT <- as.numeric(MAGMA2$ZSTAT)
      MAGMA2$P <- as.numeric(MAGMA2$P)
      rownames(MAGMA2) <- MAGMA2$GENE
      # RAPDB <- openxlsx::read.xlsx(paste0(path, "/", "RAPDB_heading.xlsx"), colNames = FALSE)
      # RAPDB <- data.frame(MAGMA_P = MAGMA2[RAPDB$X1,c("ZSTAT","P")],
      #                     RAPDB)
      # 
      # ggplot()
      MAGMA_gene <- MAGMA2$ZSTAT
      names(MAGMA_gene) <- MAGMA2$GENE
      MAGMA_gene <- sort(MAGMA_gene, decreasing = T)
      MAGMA_gene <- MAGMA_gene + rnorm(length(MAGMA_gene), mean = 0.0001, sd = 0.002)
      MAGMA_gene <- sort(MAGMA_gene, decreasing = T)
      RNA_celltype_marker_2 <- RNA_celltype_marker[which(RNA_celltype_marker$gene %in% names(MAGMA_gene)),]
      RNA_celltype_marker_2 <- RNA_celltype_marker_2[which(RNA_celltype_marker_2$avg_log2FC >= 1),]
      celltypes <- unique(RNA_celltype_marker_2$cluster)
      celltypes <- as.character(celltypes)
      P <- c()
      NES <- c()
      for (cell in celltypes) {
        temp_cell_markers <- RNA_celltype_marker_2[which(RNA_celltype_marker_2$cluster == cell),]
        # temp_cell_markers <- temp_cell_markers[which(temp_cell_markers$avg_log2FC >= 1),]
        gsea_res <- GSEA(MAGMA_gene, TERM2GENE = temp_cell_markers[,c("cluster", "gene")], pvalueCutoff = 1, pAdjustMethod = "BH")
        if (nrow(gsea_res@result) == 0) {
          p <- NaN
          nes <- NaN
        } else {
          p <- gsea_res@result$pvalue
          nes <- gsea_res@result$NES
        }
        P <- c(P, p)
        NES <- c(NES, nes)
      }
      names(P) <- celltypes
      names(NES) <- celltypes
      P_list <- c(P_list, list(P))
      NES_list <- c(NES_list, list(NES))
      
    }
    
    
    SNP <- data.table::fread(paste0(path, "/output/result.assoc.txt"),
                             sep = "\t", header = T)
    SNP <- as.data.frame(SNP)
    SNP_position <- data.frame(chr = SNP$chr,
                               start = SNP$ps,
                               end = SNP$ps+1)
    SNP_position$ID <- SNP$rs
    SNP_position$allele1 <- SNP$allele1
    SNP_position$allele0 <- SNP$allele0
    SNP_position$p_wald <- SNP$p_wald
    SNP_position_GR <- makeGRangesFromDataFrame(SNP_position,
                                                keep.extra.columns = TRUE,
                                                seqnames.field = "chr",
                                                start.field = "start",
                                                end.field = "end",
                                                ignore.strand = TRUE)
    SNP_position_GR_nearstgene <- .fastAnnoPeaks(peaks = SNP_position_GR,
                                                 BSgenome = BSgenome.OSativa.NCBI.IRGSPv1.0,
                                                 geneAnnotation = gene_annotation,
                                                 promoterRegion = c(2000, 100))
    SNP_position_GR_nearstgene_DF <- as.data.frame(SNP_position_GR_nearstgene)
    table(SNP_position_GR_nearstgene_DF$peakType)
    
  }
  
  ###### QTL/QTG
  Yield_components <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 1)
  Heading_date <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 2)
  Plant_Architecture <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 3)
  Seed_Morphology <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 4)
  Taste_quality <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 5)
  Secondary_metabolism <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 6)
  Biotic_Stress <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 7)
  Abiotic_Stress <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 8)
  Others <- openxlsx::read.xlsx(xlsxFile = "./Ref/RiceNavi/riceQTL.xlsx", sheet = 9)
  
  QTL <- as.data.frame(rbind(Yield_components, Heading_date, Plant_Architecture,
                             Seed_Morphology, Taste_quality, Secondary_metabolism,
                             Biotic_Stress, Abiotic_Stress, Others))
  QTL <- QTL[-which(QTL$Chr == "mitochondrial"),]
  QTL$Pos_7.0 <- as.numeric(QTL$Pos_7.0)
  QTL <- QTL[!is.na(QTL$Pos_7.0),]
  unique(QTL$Chr)
  # QTL <- QTL[-which(is.na(QTL$RAP.ID)),]
  QTL <- QTL[!is.na(QTL$RAP.ID),]
  QTL <- QTL[!is.na(QTL$Alt_Allele_Function),]
  
  unique(QTL$Chr)
  CHR <- data.frame(Chr = paste0("Chr", 1:12),
                    Number = c(1:12),
                    row.names = paste0("Chr", 1:12))
  QTL_position <- data.frame(chr = QTL$Chr,
                             start = QTL$Pos_7.0,
                             end = QTL$Pos_7.0+1)
  QTL_position$gene <- QTL$Gene
  QTL_position$rapid <- QTL$RAP.ID
  QTL_position$chr_num <- CHR[QTL_position$chr, "Number"]
  # QTL_position$Chr <- QTL$Chr
  QTL_position$position <- QTL$Pos_7.0
  
  QTL_position_GR <- makeGRangesFromDataFrame(QTL_position,
                                              keep.extra.columns = TRUE,
                                              seqnames.field = "chr_num",
                                              start.field = "start",
                                              end.field = "end",
                                              ignore.strand = TRUE)
  
  gene_id_map <- readRDS("gene_id_map.rds")
  rownames(gene_id_map) <- gene_id_map$gene_id
  QTL_position_GR$symbol <- gene_id_map[QTL_position_GR$rapid,"symbol"]
  
  QTL_position_GR_nearstgene <- .fastAnnoPeaks(peaks = QTL_position_GR,
                                               BSgenome = BSgenome.OSativa.NCBI.IRGSPv1.0,
                                               geneAnnotation = gene_annotation,
                                               promoterRegion = c(2000, 100))
  
}
if (F) {
  ######## Intragenic
  Intragenic <- All_SNP_position_filtered[which(All_SNP_position_filtered$SNP_region == "Intragenic"),]
  length(unique(Intragenic$TRAIT_Renamed)) # 257
  temp <- as.data.frame.array(table(Intragenic$TRAIT_Renamed, Intragenic$nearestGene))
  temp <- as.matrix(temp)
  temp[which(temp > 1)] <- 1
  temp <- as.data.frame(temp)
  trait_gene <- rowSums(temp)
  trait_gene <- trait_gene[which(trait_gene >= 5)]
  Intragenic_filtered <- Intragenic[which(Intragenic$TRAIT_Renamed %in% names(trait_gene)),] # 113 trait passed QC
  Intragenic_filtered$gene_id <- gene_id_map[Intragenic_filtered$nearestGene, "gene_id"]
  sum(is.na(Intragenic_filtered$gene_id)) # 0
  Intragenic_trait_genes <- split(Intragenic_filtered, Intragenic_filtered$TRAIT_Renamed)
  nrow(Intragenic_filtered) # 11459
  length(unique(Intragenic_filtered$SNP_CHR_BP)) # 9758
  length(unique(Intragenic_filtered$TRAIT_Renamed)) # 113
  
  # celltype enrichment
  {
    RNA_celltype_marker <- readRDS("RNA_celltype_marker.rds")
    RNA_celltype_marker <- RNA_celltype_marker[which(RNA_celltype_marker$avg_log2FC >= 1),]
    RNA_celltype_marker <- RNA_celltype_marker[which(RNA_celltype_marker$p_val_adj < 0.05),]
    unique(RNA_celltype_marker$cluster)
    RNA_celltype_marker <- RNA_celltype_marker[-which(RNA_celltype_marker$cluster == "Unidentified"),]
    table(RNA_celltype_marker$cluster)
    celltype_enrichment_matrix <- matrix(1, nrow = length(unique(as.character(RNA_celltype_marker$cluster))),
                                         ncol = length(names(Intragenic_trait_genes)))
    rownames(celltype_enrichment_matrix) <- unique(as.character(RNA_celltype_marker$cluster))
    colnames(celltype_enrichment_matrix) <- names(Intragenic_trait_genes)
    celltype_enrichment_matrix_binary <- matrix(0, nrow = length(unique(as.character(RNA_celltype_marker$cluster))),
                                                ncol = length(names(Intragenic_trait_genes)))
    rownames(celltype_enrichment_matrix_binary) <- unique(as.character(RNA_celltype_marker$cluster))
    colnames(celltype_enrichment_matrix_binary) <- names(Intragenic_trait_genes)
    # test_result <- c()
    for (cell in rownames(celltype_enrichment_matrix)) {
      print(cell)
      for (trait in colnames(celltype_enrichment_matrix)) {
        cell_markers <- RNA_celltype_marker[which(RNA_celltype_marker$cluster == cell),]
        cell_markers <- cell_markers$gene
        cell_markers <- cell_markers[which(cell_markers %in% gene_id_map$gene_id)]
        trait_genes <- Intragenic_trait_genes[[trait]]
        trait_genes <- unique(trait_genes$gene_id)
        common <- intersect(cell_markers, trait_genes)
        p <- phyper((length(common) - 1), length(trait_genes),
                    (nrow(gene_id_map)-length(trait_genes)),
                    length(cell_markers), lower.tail = F)
        celltype_enrichment_matrix[cell, trait] <- p
        if (p < 0.05) {
          celltype_enrichment_matrix_binary[cell, trait] <- 1
        }
        # hit_value <- length(common) / length(cell_markers) # 共同marker / 细胞类型marker
        # bg_value <- length(trait_genes) / nrow(gene_id_map) # trait的marker / 总基因数目
        # fe <- hit_value / bg_value
        # test_result <- as.data.frame(rbind(test_result,
        #                                    data.frame(cell = cell,
        #                                           trait = trait,
        #                                               p = p,
        #                                               fold_enrichment = fe)))
      }
    }
    # ggplot(data = test_result, aes(x = p, y = fold_enrichment)) +
    #   geom_point()
    # min(test_result[which(test_result$p < 0.05), "fold_enrichment"])
    trait_sig <- colSums(celltype_enrichment_matrix_binary)
    trait_sig <- trait_sig[which(trait_sig > 0)]
    length(trait_sig)
    celltype_enrichment_matrix <- celltype_enrichment_matrix[,which(colnames(celltype_enrichment_matrix) %in% names(trait_sig))]
    
    
    if (!is.null(celltype_enrichment_matrix)) {
      pmt <- celltype_enrichment_matrix
      pmt <- as.data.frame(pmt)
      ssmt <- celltype_enrichment_matrix <= 0.0001
      pmt[ssmt] <- "++++"
      smt <- celltype_enrichment_matrix > 0.0001 & celltype_enrichment_matrix <= 0.001
      pmt[smt] <- "+++"
      smt2 <- celltype_enrichment_matrix > 0.001 & celltype_enrichment_matrix <= 0.01
      pmt[smt2] <- "++"
      smt3 <- celltype_enrichment_matrix > 0.01 & celltype_enrichment_matrix <= 0.05
      pmt[smt3] <- "+"
      pmt[!ssmt&!smt&!smt2&!smt3] <- ''
    }
    
    pdf("Figure7_celltype_enrichment_matrix_gene.pdf", height = 10, width = 11.5)
    pheatmap(t(celltype_enrichment_matrix),
             color = colorRampPalette(c("#B83D3D", "white", "#1A5592"))(100),
             scale = "none", cluster_row = T,
             cluster_col = T,                
             display_numbers = t(pmt),
             fontsize_number = 12,
             number_color = "black",
             cellwidth = 15, cellheight = 10)
    dev.off()
    
    if (!is.null(celltype_enrichment_matrix)) {
      pmt <- celltype_enrichment_matrix
      pmt <- as.data.frame(pmt)
      ssmt <- celltype_enrichment_matrix <= 0.0001
      pmt[ssmt] <- 1
      smt <- celltype_enrichment_matrix > 0.0001 & celltype_enrichment_matrix <= 0.001
      pmt[smt] <- 1
      smt2 <- celltype_enrichment_matrix > 0.001 & celltype_enrichment_matrix <= 0.01
      pmt[smt2] <- 1
      smt3 <- celltype_enrichment_matrix > 0.01 & celltype_enrichment_matrix <= 0.05
      pmt[smt3] <- 1
      pmt[!ssmt&!smt&!smt2&!smt3] <- 0
    }
    
    temp <- data.frame(rowSums(pmt))
    temp$Celltype <- rownames(temp)
    colnames(temp)[1] <- "Count"
    temp <- temp[order(temp$Count, decreasing = T),]
    temp$Celltype <- factor(temp$Celltype, levels = temp$Celltype)
    
    pdf("Figure7.celltype_enriched_trait_count_gene.pdf", height = 4, width = 10.5)
    ggplot(data = temp, aes(x = Celltype, y = Count)) +
      geom_bar(stat = "identity", width = 0.6, fill = "#476D87") +
      theme_bw() +
      theme(axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10, color = "black"),
            legend.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10, color = "black")) +
      labs(x = "", y = "Enriched trait count")
    dev.off()
    
    gene_celltype_traits <- names(trait_sig)
  }
  # tissue enrichment
  RNA <- readRDS("./Rice-snRNA/all_data_annotated.rds")
  table(RNA$tissues)
  tissues_db <- data.frame(tissue = c("10dLeaf", "10dRoot", "120dFlagLeaf", "60dRootTip",
                                      "90dFlagLeaf", "IM0.5cm", "IM1cm", "MaturePollen", "SAM", "Zygote"),
                           main = c("Leaf", "Root", "FlagLeaf", "Root",
                                    "FlagLeaf", "IM", "IM", "MaturePollen", "SAM", "Zygote"),
                           row.names = c("10dLeaf", "10dRoot", "120dFlagLeaf", "60dRootTip",
                                         "90dFlagLeaf", "IM0.5cm", "IM1cm", "MaturePollen", "SAM", "Zygote"))
  RNA$tissues_main <- tissues_db[RNA$tissues, "main"]
  table(RNA$tissues_main)
  Idents(RNA) <- RNA$tissues_main
  RNA_tissue_marker <- FindAllMarkers(RNA, only.pos = T, logfc.threshold = 1)
  RNA_tissue_marker <- RNA_tissue_marker[which(RNA_tissue_marker$p_val_adj < 0.05),]
  {
    unique(RNA_tissue_marker$cluster)
    table(RNA_tissue_marker$cluster)
    tissue_enrichment_matrix <- matrix(1, nrow = length(unique(as.character(RNA_tissue_marker$cluster))),
                                       ncol = length(names(Intragenic_trait_genes)))
    rownames(tissue_enrichment_matrix) <- unique(as.character(RNA_tissue_marker$cluster))
    colnames(tissue_enrichment_matrix) <- names(Intragenic_trait_genes)
    tissue_enrichment_matrix_binary <- matrix(0, nrow = length(unique(as.character(RNA_tissue_marker$cluster))),
                                              ncol = length(names(Intragenic_trait_genes)))
    rownames(tissue_enrichment_matrix_binary) <- unique(as.character(RNA_tissue_marker$cluster))
    colnames(tissue_enrichment_matrix_binary) <- names(Intragenic_trait_genes)
    for (tissue in rownames(tissue_enrichment_matrix)) {
      print(tissue)
      for (trait in colnames(tissue_enrichment_matrix)) {
        tissue_markers <- RNA_tissue_marker[which(RNA_tissue_marker$cluster == tissue),]
        tissue_markers <- tissue_markers$gene
        tissue_markers <- tissue_markers[which(tissue_markers %in% gene_id_map$gene_id)]
        trait_genes <- Intragenic_trait_genes[[trait]]
        trait_genes <- unique(trait_genes$gene_id)
        common <- intersect(tissue_markers, trait_genes)
        p <- phyper((length(common) - 1), length(trait_genes),
                    (nrow(gene_id_map)-length(trait_genes)),
                    length(tissue_markers), lower.tail = F)
        tissue_enrichment_matrix[tissue, trait] <- p
        if (p < 0.05) {
          tissue_enrichment_matrix_binary[tissue, trait] <- 1
        }
      }
    }
    
    trait_sig <- colSums(tissue_enrichment_matrix_binary)
    trait_sig <- trait_sig[which(trait_sig > 0)]
    tissue_enrichment_matrix <- tissue_enrichment_matrix[,which(colnames(tissue_enrichment_matrix) %in% names(trait_sig))]
    
    
    if (!is.null(tissue_enrichment_matrix)) {
      pmt <- tissue_enrichment_matrix
      pmt <- as.data.frame(pmt)
      ssmt <- tissue_enrichment_matrix <= 0.0001
      pmt[ssmt] <- "++++"
      smt <- tissue_enrichment_matrix > 0.0001 & tissue_enrichment_matrix <= 0.001
      pmt[smt] <- "+++"
      smt2 <- tissue_enrichment_matrix > 0.001 & tissue_enrichment_matrix <= 0.01
      pmt[smt2] <- "++"
      smt3 <- tissue_enrichment_matrix > 0.01 & tissue_enrichment_matrix <= 0.05
      pmt[smt3] <- "+"
      pmt[!ssmt&!smt&!smt2&!smt3] <- ''
    }
    
    pdf("Figure7_tissue_enrichment_matrix_gene.pdf", height = 6, width = 6)
    pheatmap(t(tissue_enrichment_matrix),
             color = colorRampPalette(c("#B83D3D", "white", "#1A5592"))(100),
             scale = "none", cluster_row = T,
             cluster_col = T,                
             display_numbers = t(pmt),
             fontsize_number = 12,
             number_color = "black",
             cellwidth = 15, cellheight = 10)
    dev.off()
    
    if (!is.null(tissue_enrichment_matrix)) {
      pmt <- tissue_enrichment_matrix
      pmt <- as.data.frame(pmt)
      ssmt <- tissue_enrichment_matrix <= 0.0001
      pmt[ssmt] <- 1
      smt <- tissue_enrichment_matrix > 0.0001 & tissue_enrichment_matrix <= 0.001
      pmt[smt] <- 1
      smt2 <- tissue_enrichment_matrix > 0.001 & tissue_enrichment_matrix <= 0.01
      pmt[smt2] <- 1
      smt3 <- tissue_enrichment_matrix > 0.01 & tissue_enrichment_matrix <= 0.05
      pmt[smt3] <- 1
      pmt[!ssmt&!smt&!smt2&!smt3] <- 0
    }
    
    temp <- data.frame(rowSums(pmt))
    temp$Tissue <- rownames(temp)
    colnames(temp)[1] <- "Count"
    temp <- temp[order(temp$Count, decreasing = T),]
    temp$Tissue <- factor(temp$Tissue, levels = temp$Tissue)
    
    pdf("Figure7.tissue_enriched_trait_count_gene.pdf", height = 4, width = 3)
    ggplot(data = temp, aes(x = Tissue, y = Count)) +
      geom_bar(stat = "identity", width = 0.6, fill = "#476D87") +
      theme_bw() +
      theme(axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10, color = "black"),
            legend.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10, color = "black")) +
      labs(x = "", y = "Enriched trait count")
    dev.off()
    
    gene_tissue_traits <- names(trait_sig)
  }
  
  ## RNA & ATAC
  library(VennDiagram)
  Tvenn <- venn.diagram(list("snRNA enriched traits" = gene_celltype_traits,
                             "ATAC enriched traits" = peak_celltype_traits),
                        filename = NULL,
                        lwd = 1, lty = 2,
                        col = c("#00BFC4", "#F8766D") ,
                        fill = c("#00BFC4", "#F8766D"),
                        cat.col = c("black", "black"),
                        main = "Celltype", main.cex = 2,
                        reverse = TRUE)
  pdf("Figure7_celltype_enriched_traits_ATAC_RNA.pdf", height = 3.5, width = 4)
  grid.draw(Tvenn)
  dev.off()
  
  Tvenn <- venn.diagram(list("snRNA enriched traits" = gene_tissue_traits,
                             "ATAC enriched traits" = peak_tissue_traits),
                        filename = NULL,
                        lwd = 1, lty = 2,
                        col = c("#00BFC4", "#F8766D") ,
                        fill = c("#00BFC4", "#F8766D"),
                        cat.col = c("black", "black"),
                        main = "Tissue", main.cex = 2,
                        reverse = TRUE)
  pdf("Figure7_tissue_enriched_traits_ATAC_RNA.pdf", height = 3.5, width = 4)
  grid.draw(Tvenn)
  dev.off()
  
  ## RNA & ATAC celltype sankey plot
  {
    celltype_enrichment_matrix_peak <- matrix(1, nrow = length(markerPeaksList_Celltype),
                                              ncol = length(names(Intergenic_trait_peaks)))
    rownames(celltype_enrichment_matrix_peak) <- names(markerPeaksList_Celltype)
    colnames(celltype_enrichment_matrix_peak) <- names(Intergenic_trait_peaks)
    celltype_enrichment_matrix_peak_binary_peak <-  matrix(0, nrow = length(markerPeaksList_Celltype),
                                                           ncol = length(names(Intergenic_trait_peaks)))
    rownames(celltype_enrichment_matrix_peak_binary_peak) <- names(markerPeaksList_Celltype)
    colnames(celltype_enrichment_matrix_peak_binary_peak) <- names(Intergenic_trait_peaks)
    for (cell in rownames(celltype_enrichment_matrix_peak)) {
      print(cell)
      for (trait in colnames(celltype_enrichment_matrix_peak)) {
        cell_markers <- markerPeaksList_Celltype[[cell]]
        cell_markers <- as.data.frame(cell_markers)
        cell_markers <- paste0(cell_markers$seqnames, "_", cell_markers$start, "_", cell_markers$end)
        trait_peaks <- Intergenic_trait_peaks[[trait]]
        trait_peaks <- unique(trait_peaks$peak_id)
        common <- intersect(cell_markers, trait_peaks)
        p <- phyper((length(common) - 1), length(trait_peaks),
                    (132008-length(trait_peaks)),
                    length(cell_markers), lower.tail = F)
        celltype_enrichment_matrix_peak[cell, trait] <- p
        if (p < 0.05) {
          celltype_enrichment_matrix_peak_binary_peak[cell, trait] <- 1
        }
      }
    }
    celltype_enrichment_matrix_peak_binary_peak <- reshape2::melt(celltype_enrichment_matrix_peak_binary_peak)
    celltype_enrichment_matrix_peak_binary_peak <- celltype_enrichment_matrix_peak_binary_peak[which(celltype_enrichment_matrix_peak_binary_peak$value == 1),]
    
    celltype_enrichment_matrix_gene <- matrix(1, nrow = length(unique(as.character(RNA_celltype_marker$cluster))),
                                              ncol = length(names(Intragenic_trait_genes)))
    rownames(celltype_enrichment_matrix_gene) <- unique(as.character(RNA_celltype_marker$cluster))
    colnames(celltype_enrichment_matrix_gene) <- names(Intragenic_trait_genes)
    celltype_enrichment_matrix_gene_binary_gene <- matrix(0, nrow = length(unique(as.character(RNA_celltype_marker$cluster))),
                                                          ncol = length(names(Intragenic_trait_genes)))
    rownames(celltype_enrichment_matrix_gene_binary_gene) <- unique(as.character(RNA_celltype_marker$cluster))
    colnames(celltype_enrichment_matrix_gene_binary_gene) <- names(Intragenic_trait_genes)
    # test_result <- c()
    for (cell in rownames(celltype_enrichment_matrix_gene)) {
      print(cell)
      for (trait in colnames(celltype_enrichment_matrix_gene)) {
        cell_markers <- RNA_celltype_marker[which(RNA_celltype_marker$cluster == cell),]
        cell_markers <- cell_markers$gene
        cell_markers <- cell_markers[which(cell_markers %in% gene_id_map$gene_id)]
        trait_genes <- Intragenic_trait_genes[[trait]]
        trait_genes <- unique(trait_genes$gene_id)
        common <- intersect(cell_markers, trait_genes)
        p <- phyper((length(common) - 1), length(trait_genes),
                    (nrow(gene_id_map)-length(trait_genes)),
                    length(cell_markers), lower.tail = F)
        celltype_enrichment_matrix_gene[cell, trait] <- p
        if (p < 0.05) {
          celltype_enrichment_matrix_gene_binary_gene[cell, trait] <- 1
        }
      }
    }
    celltype_enrichment_matrix_gene_binary_gene <- reshape2::melt(celltype_enrichment_matrix_gene_binary_gene)
    celltype_enrichment_matrix_gene_binary_gene <- celltype_enrichment_matrix_gene_binary_gene[which(celltype_enrichment_matrix_gene_binary_gene$value == 1),]
    
    colnames(celltype_enrichment_matrix_peak_binary_peak) <- c("from", "to", "omics_type")
    celltype_enrichment_matrix_peak_binary_peak$omics_type <- "ATAC"
    celltype_enrichment_matrix_peak_binary_peak$from <- as.character(celltype_enrichment_matrix_peak_binary_peak$from)
    celltype_enrichment_matrix_peak_binary_peak$to <- as.character(celltype_enrichment_matrix_peak_binary_peak$to)
    celltype_enrichment_matrix_peak_binary_peak$omics_type <- as.character(celltype_enrichment_matrix_peak_binary_peak$omics_type)
    
    colnames(celltype_enrichment_matrix_gene_binary_gene) <- c("to", "from", "omics_type")
    celltype_enrichment_matrix_gene_binary_gene$omics_type <- "RNA"
    celltype_enrichment_matrix_gene_binary_gene$from <- as.character(celltype_enrichment_matrix_gene_binary_gene$from)
    celltype_enrichment_matrix_gene_binary_gene$to <- as.character(celltype_enrichment_matrix_gene_binary_gene$to)
    celltype_enrichment_matrix_gene_binary_gene$omics_type <- as.character(celltype_enrichment_matrix_gene_binary_gene$omics_type)
    
    celltype_enrichment_matrix_gene_binary_gene <- celltype_enrichment_matrix_gene_binary_gene[,c("from", "to", "omics_type")]
    
    unique(celltype_enrichment_matrix_peak_binary_peak$from)
    unique(celltype_enrichment_matrix_gene_binary_gene$to)
    celltype_enrichment_matrix_gene_binary_gene$to[which(celltype_enrichment_matrix_gene_binary_gene$to == "Spikelet meristem(SM)")] <- "Spikelet meristem (SM)"
    celltype_enrichment_matrix_gene_binary_gene$to[which(celltype_enrichment_matrix_gene_binary_gene$to == "Cryptic bract/bract(cb/b)")] <- "Cryptic bract/bract (cb/b)"
    celltype_enrichment_matrix_gene_binary_gene$to[which(celltype_enrichment_matrix_gene_binary_gene$to == "Branch meristems(BM)")] <- "Branch meristems (BM)"
    celltype_enrichment_matrix_gene_binary_gene$to[which(celltype_enrichment_matrix_gene_binary_gene$to == "Lemma(le)")] <- "Lemma (le)"
    celltype_enrichment_matrix_gene_binary_gene$to[which(celltype_enrichment_matrix_gene_binary_gene$to == "Inflorescence meristem(IM)")] <- "Inflorescence meristem (IM)"
    unique(celltype_enrichment_matrix_gene_binary_gene$to)
    
    unique(celltype_enrichment_matrix_peak_binary_peak$from)
    celltype_enrichment_matrix_peak_binary_peak$from <- paste0("ATAC-", celltype_enrichment_matrix_peak_binary_peak$from)
    unique(celltype_enrichment_matrix_gene_binary_gene$to)
    celltype_enrichment_matrix_gene_binary_gene$to <- paste0("RNA-", celltype_enrichment_matrix_gene_binary_gene$to)
    
    common_celltype_traits <- intersect(gene_celltype_traits,
                                        peak_celltype_traits)
    celltype_enrichment_matrix_peak_binary_peak <- celltype_enrichment_matrix_peak_binary_peak[which(celltype_enrichment_matrix_peak_binary_peak$to %in% common_celltype_traits),]
    celltype_enrichment_matrix_gene_binary_gene <- celltype_enrichment_matrix_gene_binary_gene[which(celltype_enrichment_matrix_gene_binary_gene$from %in% common_celltype_traits),]
    
    celltype_trait_enrich <- as.data.frame(rbind(celltype_enrichment_matrix_peak_binary_peak,
                                                 celltype_enrichment_matrix_gene_binary_gene))
    unique(celltype_trait_enrich$from)
    unique(celltype_trait_enrich$to)
    
    library(highcharter)
    sankey_data <- data.table(from = celltype_trait_enrich$from,
                              to = celltype_trait_enrich$to,
                              weight = rep(1, nrow(celltype_trait_enrich)))
    highchart() %>%
      hc_title(text = "桑基图") %>%
      hc_add_series(data = sankey_data, type = "sankey",
                    hcaes(from = from,to = to,weight = weight)) %>%
      hc_add_theme(hc_theme_google()) %>%
      hc_chart(events = list(load = JS("function() {
      var chart = this;
      chart.update({
        plotOptions: {
          series: {
            dataLabels: {
              style: {
                fontSize: '20px'
              }
            }
          }
        }
      })
    }")))
  }
  
}
#### GCP LCP
if (FALSE) {
  peak_group <- readRDS("Figure6_peak_group.rds")
  peak_type_color <- c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD")
  names(peak_type_color) <- c("Gene-linked peaks",
                              "Lineage constitutive peaks",
                              "Global constitutive peaks",
                              "Other peaks")
  ##### 1
  length(unique(Celltype_Trait_SNP_overlap$SNP))
  Celltype_Trait_SNP_overlap$peak_group <- peak_group[Celltype_Trait_SNP_overlap$Peaks, "peakType"]
  sum(is.na(Celltype_Trait_SNP_overlap$peak_group))
  # NA: 7
  Celltype_Trait_SNP_overlap$peak_group[is.na(Celltype_Trait_SNP_overlap$peak_group)] <- "Other peaks"
  temp <- as.data.frame(table(Celltype_Trait_SNP_overlap$peak_group))
  colnames(temp) <- c("Peak_Type", "Count")
  pdf("Figure7_893_SNP_Peak_type_count.pdf", width = 5.5, height = 5.5)
  pie(c(temp$Count[1],
        temp$Count[3],
        temp$Count[2],
        temp$Count[4]),
      labels = c("Gene-linked peaks (n = 372)",
                 "Lineage constitutive peaks (n = 419)",
                 "Global constitutive peaks (n = 14)",
                 "Other peaks (n = 424)"),
      col = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD"))
  dev.off()
  
  
  
  Intergenic_overlap_with_peaks_filtered_unique$peak_group <- peak_group[Intergenic_overlap_with_peaks_filtered_unique$peak_id, "peakType"]
  sum(is.na(Intergenic_overlap_with_peaks_filtered_unique$peak_group))
  # NA: 106
  Intergenic_overlap_with_peaks_filtered_unique$peak_group[is.na(Intergenic_overlap_with_peaks_filtered_unique$peak_group)] <- "Other peaks"
  temp <- as.data.frame(table(Intergenic_overlap_with_peaks_filtered_unique$peak_group))
  colnames(temp) <- c("Peak_Type", "Count")
  
  pdf("Figure7_8240_SNP_Peak_type_count.pdf", width = 5.5, height = 5.5)
  pie(c(temp$Count[1],
        temp$Count[3],
        temp$Count[2],
        temp$Count[4]),
      labels = c("Gene-linked peaks (n = 1,067)",
                 "Lineage constitutive peaks (n = 1,717)",
                 "Global constitutive peaks (n = 299)",
                 "Other peaks (n = 5,157)"),
      col = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD"))
  dev.off()
  
  # X-squared = 12, df = 9, p-value = 0.2133
  chisq.test(temp$Count, as.data.frame(table(peak_group$peakType))[,2], correct = FALSE)
  
  ##### 2
  temp <- Intergenic_overlap_with_peaks_filtered_unique[which(Intergenic_overlap_with_peaks_filtered_unique$Marker_Peaks == TRUE),]
  nrow(temp)
  table(temp$peak_group)
  temp2 <- as.data.frame(table(temp$peak_group))
  colnames(temp2) <- c("Peak_Type", "Count")
  pdf("Figure7_3442_SNP_Peak_type_count.pdf", width = 5.5, height = 5.5)
  pie(c(temp2$Count[1],
        temp2$Count[3],
        temp2$Count[2],
        temp2$Count[4]),
      labels = c("Gene-linked peaks (n = 689)",
                 "Lineage constitutive peaks (n = 985)",
                 "Global constitutive peaks (n = 56)",
                 "Other peaks (n = 1,712)"),
      col = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD"))
  dev.off()
  
  # X-squared = 12, df = 9, p-value = 0.2133
  chisq.test(temp2$Count, as.data.frame(table(peak_group$peakType))[,2], correct = FALSE)
  
  pdf("Figure7_SNP_peak_group_marker_peak_percent.pdf", height = 6, width = 4)
  ggplot(data = Intergenic_overlap_with_peaks_filtered_unique,
         aes(x = peak_group, fill = Marker_Peaks)) +
    geom_bar(stat = "count", position = "fill") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 10),
          axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "", y = "Percentage")
  dev.off()
  
  ##### 3, 4798 SNPs
  temp <- Intergenic_overlap_with_peaks_filtered_unique[which(Intergenic_overlap_with_peaks_filtered_unique$Marker_Peaks == FALSE),]
  table(temp$peak_group)
  temp2 <- as.data.frame(table(temp$peak_group))
  colnames(temp2) <- c("Peak_Type", "Count")
  pdf("Figure7_4798_SNP_Peak_type_count.pdf", width = 5.5, height = 5.5)
  pie(c(temp2$Count[1],
        temp2$Count[3],
        temp2$Count[2],
        temp2$Count[4]),
      labels = c("Gene-linked peaks (n = 378)",
                 "Lineage constitutive peaks (n = 732)",
                 "Global constitutive peaks (n = 243)",
                 "Other peaks (n = 3,445)"),
      col = c("#6A5ACD", "#FD8F52", "#FD5F52", "#96CDCD"))
  dev.off()
  
  temp_GCPs <- temp[which(temp$peak_group == "Global constitutive peaks"),]
  temp_LCPs <- temp[which(temp$peak_group == "Lineage constitutive peaks"),]
  temp_GCPs_enrich_traits <- c()
  
  temp_GCPs_LCPs_enrich_matrix <- matrix(1, nrow = 1,
                                         ncol = length(names(Intergenic_trait_peaks)))
  rownames(temp_GCPs_LCPs_enrich_matrix) <- c("Con_peaks")
  colnames(temp_GCPs_LCPs_enrich_matrix) <- names(Intergenic_trait_peaks)
  temp_GCPs_LCPs_enrich_matrix_binary <-  matrix(0, nrow = 1,
                                                 ncol = length(names(Intergenic_trait_peaks)))
  rownames(temp_GCPs_LCPs_enrich_matrix_binary) <- c("Con_peaks")
  colnames(temp_GCPs_LCPs_enrich_matrix_binary) <- names(Intergenic_trait_peaks)
  temp_GCPs_LCPs_enrich_matrix_OR <- matrix(1, nrow = 1,
                                            ncol = length(names(Intergenic_trait_peaks)))
  rownames(temp_GCPs_LCPs_enrich_matrix_OR) <- c("Con_peaks")
  colnames(temp_GCPs_LCPs_enrich_matrix_OR) <- names(Intergenic_trait_peaks)
  
  # temp_peak_list <- list(GCPs = unique(temp_GCPs$peak_id),
  #                        LCPs = unique(temp_LCPs$peak_id))
  
  temp_peak_list <- list(GCPs = unique(peak_group$peakID[which(peak_group$peakType == "Global constitutive peaks")]),
                         LCPs = unique(peak_group$peakID[which(peak_group$peakType == "Lineage constitutive peaks")]))
  peaks_c <- c(temp_peak_list$GCPs, temp_peak_list$LCPs)
  peaks_c <- unique(peaks_c)
  length(peaks_c) # 30499
  
  # LCP+GCP中有多少是trait相关的peak
  {
    Intergenic_overlap_with_peaks_filtered
    Intergenic_overlap_with_peaks_filtered
    LCP_GCP_trait <- data.frame(LCP_GCP = peaks_c,
                                Trait_related = ifelse(peaks_c %in% Intergenic_overlap_with_peaks_filtered$peak_id,
                                                       1, 0))
    table(LCP_GCP_trait$Trait_related) # 29344, 1155, 0.03787009
    
    Intergenic_overlap_with_peaks_filtered$Constitutive_peaks <- ifelse(Intergenic_overlap_with_peaks_filtered$peak_id %in% peaks_c,
                                                                        1, 0)
    temp <- Intergenic_overlap_with_peaks_filtered[!duplicated(Intergenic_overlap_with_peaks_filtered[,c("TRAIT_Renamed","Constitutive_peaks")]),]
    temp <- temp[which(temp$Constitutive_peaks == 1),]
    nrow(temp) / 89
  }
  
  # 细胞类型的marker peaks有多少是trait相关的
  {
    Celltype_marker_trait_percent <- c()
    for (i in names(markerPeaksList_Celltype)) {
      temp <- markerPeaksList_Celltype[[i]]
      temp <- as.data.frame(temp)
      temp_celltype_peaks <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
      temp_percent <- sum(temp_celltype_peaks %in% Intergenic_overlap_with_peaks_filtered$peak_id) / length(temp_celltype_peaks)
      Celltype_marker_trait_percent <- c(Celltype_marker_trait_percent, temp_percent)
    }
    Celltype_marker_trait_percent <- data.frame(Celltype = names(markerPeaksList_Celltype),
                                                Trait_percent = Celltype_marker_trait_percent)
  }
  
  # Trait 角度
  {
    # 与细胞类型差异peak比较
    {
      markerPeaksList_Celltype
      Trait_celltype_marker_percent <- c()
      for (i in names(markerPeaksList_Celltype)) {
        temp <- markerPeaksList_Celltype[[i]]
        temp <- as.data.frame(temp)
        temp_celltype_peaks <- paste0(temp$seqnames, "_", temp$start, "_", temp$end)
        temp <- Intergenic_overlap_with_peaks_filtered
        temp$Celltype_marker_peaks <- ifelse(temp$peak_id %in% temp_celltype_peaks,
                                             1, 0)
        temp <- temp[!duplicated(temp[,c("TRAIT_Renamed","Celltype_marker_peaks")]),]
        temp <- temp[which(temp$Celltype_marker_peaks == 1),]
        Trait_celltype_marker_percent <- c(Trait_celltype_marker_percent, nrow(temp) / 89)
      }
      
      Trait_celltype_marker_percent <- data.frame(Celltype = names(markerPeaksList_Celltype),
                                                  Trait_percent = Trait_celltype_marker_percent)
      
      Intergenic_overlap_with_peaks_filtered$Celltype_marker_peaks <- ifelse(Intergenic_overlap_with_peaks_filtered$peak_id %in% Celltype_marker_peaks,
                                                                             1, 0)
      temp <- Intergenic_overlap_with_peaks_filtered[!duplicated(Intergenic_overlap_with_peaks_filtered[,c("TRAIT_Renamed","Celltype_marker_peaks")]),]
      temp <- temp[which(temp$Celltype_marker_peaks == 1),]
      nrow(temp) / 89 # 98.87%
      
    }
    # 与组成型peaks比较
    {
      Intergenic_overlap_with_peaks_filtered$Constitutive_peaks <- ifelse(Intergenic_overlap_with_peaks_filtered$peak_id %in% peaks_c,
                                                                          1, 0)
      temp <- Intergenic_overlap_with_peaks_filtered[!duplicated(Intergenic_overlap_with_peaks_filtered[,c("TRAIT_Renamed","Constitutive_peaks")]),]
      temp <- temp[which(temp$Constitutive_peaks == 1),]
      nrow(temp) / 89 # 0.9550562 的trait中含有组合型peak
      
      Trait_Con_Number <- c()
      for (i in names(Intergenic_trait_peaks)) {
        print(i)
        temp_trait <- Intergenic_trait_peaks[[i]]
        temp_trait <- unique(temp_trait$peak_id)
        temp_trait_con <- length(intersect(temp_trait, peaks_c))
        temp_trait_con_percent <- temp_trait_con / length(temp_trait)
        temp_sample <- lapply(1:100, function(x){
          temp_c <- sample(setdiff(peaksets_df$peak_id, peaks_c), length(unique(peaks_c)), replace = F)
          # temp_c <- sample(peaksets_df$peak_id, length(unique(peaks_c)), replace = F)
          sample_trait_con <- length(intersect(temp_trait, temp_c))
          sample_trait_con_percent <- sample_trait_con / length(temp_trait)
          c(sample_trait_con, sample_trait_con_percent)
        })
        temp_sample <- unlist(temp_sample)
        temp_sample_trait_con <- temp_sample[seq(from = 1, to = 200, by = 2)]
        temp_sample_trait_con_percent <- temp_sample[seq(from = 2, to = 200, by = 2)]
        sample_trait <- data.frame(temp_sample_trait_con = temp_sample_trait_con,
                                   temp_sample_trait_con_percent = temp_sample_trait_con_percent)
        temp <- list(trait = i, trait_con_number = temp_trait_con, trait_con_percent = temp_trait_con_percent,
                     random_trait_con_number = sample_trait$temp_sample_trait_con,
                     random_trait_con_percent = sample_trait$temp_sample_trait_con_percent,
                     fc = temp_trait_con_percent / (30499 / 132008))
        Trait_Con_Number <- c(Trait_Con_Number, list(temp))
      }
      
      Trait_Con_Number_df <- c()
      for (i in 1:length(Trait_Con_Number)) {
        temp <- Trait_Con_Number[[i]]
        percent <- temp$trait_con_percent
        p <- sum(temp$random_trait_con_percent >= temp$trait_con_percent) / 100
        percent_SE <- plotrix::std.error(temp$random_trait_con_percent)
        percent_SD <- sd(temp$random_trait_con_percent)
        temp <- data.frame(trait = temp$trait,
                           percent = percent,
                           p = p,
                           percent_SE = p_SE,
                           percent_SD = p_SD,
                           percent_MEAN = mean(temp$random_trait_con_percent),
                           fc = temp$fc)
        Trait_Con_Number_df <- as.data.frame(rbind(Trait_Con_Number_df,
                                                   temp))
      }
      
      
      
      
      # traits中含有非组成型peaks的比例
      Intergenic_overlap_with_peaks_filtered$Not_constitutive_peaks <- ifelse(Intergenic_overlap_with_peaks_filtered$peak_id %in% setdiff(peaksets_df$peak_id, peaks_c),
                                                                              1, 0)
      temp <- Intergenic_overlap_with_peaks_filtered[!duplicated(Intergenic_overlap_with_peaks_filtered[,c("TRAIT_Renamed","Not_constitutive_peaks")]),]
      temp <- temp[which(temp$Not_constitutive_peaks == 1),]
      nrow(temp) / 89 # 1 
      
      # 对组成型peaks随机
      random_c_10 <- lapply(1:10, function(x){
        # sample(peaksets_df$peak_id, length(unique(peaks_c)), replace = F)
        sample(setdiff(peaksets_df$peak_id, peaks_c), length(unique(peaks_c)), replace = F) # 从非组成型peak里面随机取
      })
      random_c_percent <- c()
      for (i in 1:10) {
        print(i)
        temp_c_peaks <- random_c_10[[i]]
        temp <- Intergenic_overlap_with_peaks_filtered
        temp$Constitutive_peaks <- ifelse(temp$peak_id %in% temp_c_peaks,
                                          1, 0)
        temp <- temp[!duplicated(temp[,c("TRAIT_Renamed","Constitutive_peaks")]),]
        temp <- temp[which(temp$Constitutive_peaks == 1),]
        random_c_percent <- c(random_c_percent, nrow(temp) / 89)
      }
      mean(random_c_percent)
      sd(random_c_percent)
      
      
      
      # 对trait的peaks随机
      trait_peak_num <- lapply(Intergenic_trait_peaks, function(x){
        length(unique(x$peak_id))
      })
      trait_peak_num <- unlist(trait_peak_num)
      
      random_c_percent <- c()
      for (i in 1:100) {
        print(i)
        traits <- names(Intergenic_trait_peaks)
        traits_random <- c()
        traits2 <- c()
        for (t in traits) {
          traits_random <- c(traits_random, sample(peaksets_df$peak_id, trait_peak_num[t], replace = F))
          traits2 <- c(traits2, rep(t, trait_peak_num[t]))
        }
        temp <- data.frame("Peaks" = traits_random,
                           "Traits" =  traits2)
        
        temp$Constitutive_peaks <- ifelse(temp$Peaks %in% peaks_c,
                                          1, 0)
        temp <- temp[!duplicated(temp[,c("Traits","Constitutive_peaks")]),]
        temp <- temp[which(temp$Constitutive_peaks == 1),]
        random_c_percent <- c(random_c_percent, nrow(temp) / 89)
      }
      mean(random_c_percent)
      sd(random_c_percent)
      
    }
  }
  
  
  
  
  
  random_c <- lapply(1:10000, function(x){
    sample(peaksets_df$peak_id, length(peaks_c), replace = F)
  })
  trait_c_common <- c()
  trait_c_common_p <- c()
  for (i in names(Intergenic_trait_peaks)) {
    print(i)
    common <- length(unique(intersect(peaks_c, Intergenic_trait_peaks[[i]]$peak_id)))
    random_common <- lapply(1:10000, function(x){
      length(unique(intersect(random_c[[x]], Intergenic_trait_peaks[[i]]$peak_id)))
    })
    random_common <- unlist(random_common)
    trait_c_common <- c(trait_c_common, length(unique(intersect(peaks_c, Intergenic_trait_peaks[[i]]$peak_id))))
    trait_c_common_p <- c(trait_c_common_p, sum(random_common >= common) / 10000)
  }
  
  trait_c_df <- data.frame(trait = names(Intergenic_trait_peaks),
                           common_peaks_number = trait_c_common,
                           P = trait_c_common_p)
  sum(trait_c_df$P < 0.05)
  
  
  # for (p in c("GCPs", "LCPs")) {
  #   print(p)
  #   for (trait in colnames(temp_GCPs_LCPs_enrich_matrix)) {
  #     temp_c <- temp_peak_list[[p]]
  #     not_temp_c <- setdiff(peaksets_df$peak_id, temp_c)
  #     temp_trait <- Intergenic_trait_peaks[[trait]]
  #     temp_trait <- unique(temp_trait$peak_id)
  #     not_temp_trait <- setdiff(peaksets_df$peak_id, temp_trait)
  #     temp_matrix <- matrix(c(length(intersect(temp_c, temp_trait)),
  #                             length(intersect(not_temp_c, temp_trait)),
  #                             length(intersect(temp_c, not_temp_trait)),
  #                             length(intersect(not_temp_c, not_temp_trait))
  #     ), nrow = 2)
  #     if (sum(temp_matrix <= 5) >= 1) {
  #       P <- chisq.test(temp_matrix, correct = T)
  #     } else {
  #       P <- chisq.test(temp_matrix, correct = F)
  #     }
  #     odds_ratio <- (length(intersect(temp_c, temp_trait)) / length(temp_c)) / (length(intersect(not_temp_c, temp_trait)) / length(not_temp_c))
  #     temp_GCPs_LCPs_enrich_matrix[p, trait] <- P$p.value
  #     temp_GCPs_LCPs_enrich_matrix_OR[p, trait] <- odds_ratio
  #     if (P$p.value < 0.05 & odds_ratio > 1) {
  #       temp_GCPs_LCPs_enrich_matrix_binary[p, trait] <- 1
  #     }
  #   }
  # }
  # 
  
  for (p in c("Con_peaks")) {
    print(p)
    for (trait in colnames(temp_GCPs_LCPs_enrich_matrix)) {
      peak_markers <- peaks_c
      trait_peaks <- Intergenic_trait_peaks[[trait]]
      trait_peaks <- unique(trait_peaks$peak_id)
      common <- intersect(peak_markers, trait_peaks)
      P <- phyper((length(common) - 1), length(trait_peaks),
                  (132008-length(trait_peaks)),
                  length(peak_markers), lower.tail = F)
      temp_GCPs_LCPs_enrich_matrix[p, trait] <- P
      if (P < 0.05) {
        temp_GCPs_LCPs_enrich_matrix_binary[p, trait] <- 1
      }
      # temp_c <- temp_peak_list[[p]]
      # not_temp_c <- setdiff(peaksets_df$peak_id, temp_c)
      # temp_trait <- Intergenic_trait_peaks[[trait]]
      # temp_trait <- unique(temp_trait$peak_id)
      # not_temp_trait <- setdiff(peaksets_df$peak_id, temp_trait)
      # temp_matrix <- matrix(c(length(intersect(temp_c, temp_trait)),
      #                         length(intersect(not_temp_c, temp_trait)),
      #                         length(intersect(temp_c, not_temp_trait)),
      #                         length(intersect(not_temp_c, not_temp_trait))
      #                         ), nrow = 2)
      # if (sum(temp_matrix <= 5) >= 1) {
      #   P <- chisq.test(temp_matrix, correct = T)
      # } else {
      #   P <- chisq.test(temp_matrix, correct = F)
      # }
      # odds_ratio <- (length(intersect(temp_c, temp_trait)) / length(temp_c)) / (length(intersect(not_temp_c, temp_trait)) / length(not_temp_c))
      # temp_GCPs_LCPs_enrich_matrix[p, trait] <- P$p.value
      # temp_GCPs_LCPs_enrich_matrix_OR[p, trait] <- odds_ratio
      # if (P$p.value < 0.05 & odds_ratio > 1) {
      #   temp_GCPs_LCPs_enrich_matrix_binary[p, trait] <- 1
      # }
    }
  }
  trait_sig <- colSums(temp_GCPs_LCPs_enrich_matrix_binary)
  trait_sig <- trait_sig[which(trait_sig > 0)]
  temp_GCPs_LCPs_enrich_matrix <- temp_GCPs_LCPs_enrich_matrix[,which(colnames(temp_GCPs_LCPs_enrich_matrix) %in% names(trait_sig))]
  
  if (!is.null(temp_GCPs_LCPs_enrich_matrix)) {
    pmt <- temp_GCPs_LCPs_enrich_matrix
    pmt <- as.data.frame(pmt)
    ssmt <- temp_GCPs_LCPs_enrich_matrix <= 0.0001
    pmt[ssmt] <- "****"
    smt <- temp_GCPs_LCPs_enrich_matrix > 0.0001 & temp_GCPs_LCPs_enrich_matrix <= 0.001
    pmt[smt] <- "***"
    smt2 <- temp_GCPs_LCPs_enrich_matrix > 0.001 & temp_GCPs_LCPs_enrich_matrix <= 0.01
    pmt[smt2] <- "**"
    smt3 <- temp_GCPs_LCPs_enrich_matrix > 0.01 & temp_GCPs_LCPs_enrich_matrix <= 0.05
    pmt[smt3] <- "*"
    pmt[!ssmt&!smt&!smt2&!smt3] <- ''
  }
  
  pdf("Figure7_4798_GCPs_LCPs_enrich_matrix.pdf", height = 5, width = 20)
  pheatmap(-log10(temp_GCPs_LCPs_enrich_matrix),
           color = colorRampPalette(c("yellow", "red"))(20),
           scale = "none", cluster_row = T,
           cluster_col = T,                
           display_numbers = pmt,
           fontsize_number = 6,
           number_color = "black",
           cellwidth = 15, cellheight = 20)
  dev.off()
  
  ## 挨个统计
  #### 对于每一个细胞类型和表型对，分别统计
  GCP_LCP_enrichment_matrix_binary_df <-reshape::melt(temp_GCPs_LCPs_enrich_matrix_binary)
  GCP_LCP_enrichment_matrix_binary_df <- GCP_LCP_enrichment_matrix_binary_df[which(GCP_LCP_enrichment_matrix_binary_df$value == 1),]
  GCP_LCP_Stat <- c()
  GCP_LCP_Trait_SNP_remove <- c()
  GCP_LCP_Trait_SNP_overlap <- c() # 和GCP、LCP有overlap
  GCP_LCP_Trait_SNP_all <- c() # 显著富集的trait的所有SNP
  for (i in 1:nrow(GCP_LCP_enrichment_matrix_binary_df)) {
    print(i)
    trait <- GCP_LCP_enrichment_matrix_binary_df[i, 2]
    trait <- as.character(trait)
    trait_peaks <- Intergenic_trait_peaks[[trait]]
    trait_peaks <- unique(trait_peaks$peak_id)
    trait_SNP <- Intergenic_trait_peaks[[trait]]
    GCP_LCP_Trait_SNP_all <- c(GCP_LCP_Trait_SNP_all, unique(trait_SNP$SNP_CHR_BP))
    p <- GCP_LCP_enrichment_matrix_binary_df[i, 1]
    p <- as.character(p)
    peak_markers <- temp_peak_list[[p]]
    common_peaks <- intersect(peak_markers, trait_peaks)
    for (peak in common_peaks) {
      trait_SNP <- Intergenic_trait_peaks[[trait]]
      trait_SNP <- trait_SNP[which(trait_SNP$peak_id == peak),]
      trait_SNP <- trait_SNP[!duplicated(trait_SNP$SNP_CHR_BP),]
      rownames(trait_SNP) <- trait_SNP$SNP_CHR_BP
      for (snp in trait_SNP$SNP_CHR_BP) { # 对每一个SNP循环
        temp_trait_SNP <- trait_SNP[snp,]
        GCP_LCP_Trait_SNP_overlap <- as.data.frame(rbind(GCP_LCP_Trait_SNP_overlap,
                                                         data.frame(Peaks_group = p,
                                                                    Trait = trait,
                                                                    SNP = temp_trait_SNP$SNP_CHR_BP,
                                                                    Peaks = temp_trait_SNP$peak_id)))
        temp_ref_alt_peak_TF_SNP <- ref_alt_peak_TF_SNP[which(ref_alt_peak_TF_SNP$SNP %in% snp),]
        if (nrow(temp_ref_alt_peak_TF_SNP) == 0) {
          GCP_LCP_Trait_SNP_remove <- as.data.frame(rbind(GCP_LCP_Trait_SNP_remove,
                                                          data.frame(Peaks_group = p,
                                                                     Trait = trait,
                                                                     SNP = temp_trait_SNP$SNP_CHR_BP,
                                                                     Peaks = temp_trait_SNP$peak_id)))
        } else {
          print("Yes!")
          GCP_LCP_Stat_temp <- data.frame(Peaks_group = p,
                                          Trait = trait,
                                          temp_ref_alt_peak_TF_SNP)
          
          GCP_LCP_Stat <- as.data.frame(rbind(GCP_LCP_Stat,
                                              GCP_LCP_Stat_temp))
        }
      }
    }
  }
  
  length(unique(GCP_LCP_Trait_SNP_overlap$SNP)) # 967
  length(unique(GCP_LCP_Trait_SNP_overlap$Trait)) # 76
  
  GCP_LCP_Stat$Peak <- overlapRegions_DF[GCP_LCP_Stat$SNP, 2]
  for (i in 1:nrow(GCP_LCP_Stat)) {
    if (GCP_LCP_Stat$ref_peak_TF[i] == "" & GCP_LCP_Stat$alt_peak_TF[i] != "") {
      GCP_LCP_Stat$TF[i] <- GCP_LCP_Stat$alt_peak_TF[i]
    }
    if (GCP_LCP_Stat$ref_peak_TF[i] != "" & GCP_LCP_Stat$alt_peak_TF[i] == "") {
      GCP_LCP_Stat$TF[i] <- GCP_LCP_Stat$ref_peak_TF[i]
    }
    if (GCP_LCP_Stat$ref_peak_TF[i] != "" & GCP_LCP_Stat$alt_peak_TF[i] != "") {
      if (GCP_LCP_Stat$ref_peak_TF[i] == GCP_LCP_Stat$alt_peak_TF[i]) {
        GCP_LCP_Stat$TF[i] <- GCP_LCP_Stat$ref_peak_TF[i]
      } else {
        GCP_LCP_Stat$TF[i] <- ""
      }
    }
  }
  
  length(unique(GCP_LCP_Stat$SNP)) # 316 SNPs
  length(unique(GCP_LCP_Stat$Peaks_group)) # 2 Celltypes
  length(unique(GCP_LCP_Stat$Trait)) # 49 Traits
  length(unique(GCP_LCP_Stat$TF)) # 211 TFs
  nrow(GCP_LCP_Stat)
  
  heatmap_value <- c()
  for (i in 1:nrow(GCP_LCP_Stat)) {
    temp_ref <- GCP_LCP_Stat$ref_peak_motif_strand[i] == ""
    temp_alt <- GCP_LCP_Stat$alt_peak_motif_strand[i] == ""
    if (temp_ref == TRUE & temp_alt == FALSE) { # 突变导致结合
      heatmap_value <- c(heatmap_value, "1")
    }
    if (temp_ref == FALSE & temp_alt == FALSE) { # 突变不影响结合
      heatmap_value <- c(heatmap_value, "2")
    }
    if (temp_ref == FALSE & temp_alt == TRUE) { # 突变导致不能结合
      heatmap_value <- c(heatmap_value, "3")
    }
  }
  table(heatmap_value)
  
  
  #### 组成型开放的基因 ATAC
  getAvailableMatrices(proj_pass_filter)
  GeneScoreMatrix <- getMatrixFromProject(proj_pass_filter, useMatrix = "GeneScoreMatrix")
  Genes <- GeneScoreMatrix@elementMetadata@listData[["name"]]
  GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]@Dimnames[[1]] <- Genes
  Celltypes
  temp_matrix <- matrix(0, nrow = length(Genes), ncol = length(Celltypes))
  rownames(temp_matrix) <- Genes
  colnames(temp_matrix) <- Celltypes
  for (i in Celltypes) {
    print(i)
    temp <- GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]][,which(GeneScoreMatrix@colData@listData[["Celltype"]] == i)]
    temp <- rowSums(temp) / ncol(temp)
    threshold <- quantile(temp, probs = 0.5)
    temp <- ifelse(temp > threshold, 1, 0)
    temp_matrix[names(temp), i] <- temp
  }
  temp_matrix_2 <- rowSums(temp_matrix)
  constitutive_genes <- names(temp_matrix_2)[which(temp_matrix_2 == length(Celltypes))]
  
  markerGenesList_celltype <- getMarkers(markersGenes_celltype, cutOff = "FDR <= 0.05 & Log2FC >= 1")
  markerGenesList_celltype <- markerGenesList_celltype@listData
  marker_genes_celltypes <- c()
  for (i in Celltypes) {
    temp <- unique(markerGenesList_celltype[[i]]$name)
    marker_genes_celltypes <- c(marker_genes_celltypes, temp)
  }
  marker_genes_celltypes <- unique(marker_genes_celltypes)
  constitutive_genes <- constitutive_genes[-which(constitutive_genes %in% marker_genes_celltypes)]
  
  Gene_diff <- c()
  TF_diff <- c()
  for (i in 1:nrow(GCP_LCP_Stat)) {
    tf <- GCP_LCP_Stat$TF[i]
    tf_id <- Annotation_genes[tf, "Locus_ID"]
    tf_symbol <- gene_id_map[which(gene_id_map$gene_id == tf_id), "symbol"]
    
    gene <- GCP_LCP_Stat$SNP_nearestGene[i]
    
    Gene_diff <- c(Gene_diff, ifelse(gene %in% constitutive_genes, 1, 0))
    TF_diff <- c(TF_diff, ifelse(tf %in% constitutive_genes, 1, 0))
  }
  
  GCP_LCP_Stat$Gene_diff <- Gene_diff
  GCP_LCP_Stat$TF_diff <- TF_diff
  
  GCP_LCP_Stat_filtered <- GCP_LCP_Stat[which(GCP_LCP_Stat$Gene_diff == 1 & GCP_LCP_Stat$TF_diff == 1),]
  openxlsx::write.xlsx(GCP_LCP_Stat_filtered, "Figure7_GCP_LCP_Stat_filtered.xlsx")
}