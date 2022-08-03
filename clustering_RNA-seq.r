library(ClustVarLV)
library(ClusterR)
library(pheatmap)
library(grid)
library(amap)
library(RColorBrewer)

setwd("C:/project-hsf/clustering_RNA-seq")
# week2 dataset
df.all <- read.table("TPM_all.txt",
                     header = T,
                     row.names = 1,
                     sep = "\t")
grep("Col",colnames(df.all),value = T)

wt <- df.all[,c(grep("Col",colnames(df.all),value = T))]

wt.mean <- data.frame(Col_0_22C = (wt$Col_0_22C_Rep1_217R + wt$Col_0_22C_Rep3_217R + wt$Col_0_22C_218R + wt$Col_0_22C_rep2_219R)/4,
                      Col_0_37C5min = (wt$Col_0_37C5min_218R + wt$Col_0_37C5min_rep2_219R)/2,
                      Col_0_37C15min = (wt$Col_0_37C15min_218R + wt$Col_0_37C15min_rep2_219R)/2,
                      Col_0_37C30min = (wt$Col_0_37C_Rep1_217R + wt$Col_0_37C_Rep3_217R)/2)
rownames(wt.mean) <- rownames(wt)

wt.mean <- log2(wt.mean + 1)
# wt.mean <- wt.mean[apply(wt.mean,1,max)>2, ]

DEG_list <- read.table("C:/project-hsf/edgeR/DEGs_37Cvs22C_atLeastInOnePairwise.q0.01_lfc1.geneID")
wt.mean <- wt.mean[as.vector(DEG_list$V1),]
wt.mean <- wt.mean[grep("AT",rownames(wt.mean)),]

wt.hs_vs_rt <- data.frame(WT_37C5min_vs_WT_22C = wt.mean$Col_0_37C5min - wt.mean$Col_0_22C,
                          WT_37C15min_vs_WT_22C = wt.mean$Col_0_37C15min - wt.mean$Col_0_22C,
                          WT_37C30min_vs_WT_22C = wt.mean$Col_0_37C30min - wt.mean$Col_0_22C)

rownames(wt.hs_vs_rt) <- rownames(wt.mean)

wt.hs_vs_rt[wt.hs_vs_rt > 2] <- 2
wt.hs_vs_rt[wt.hs_vs_rt < -2] <- -2
data <- wt.hs_vs_rt
opt_gmm <- Optimal_Clusters_GMM(data,
                                max_clusters = 8,
                                criterion = "BIC",
                                dist_mode = "eucl_dist",
                                seed_mode = "static_spread",
                                km_iter = 30,
                                em_iter = 30,
                                var_floor = 1e-10,
                                plot_data = T)
clusterN <- 8
gmm <- GMM(data,
           gaussian_comps = clusterN,
           dist_mode = "eucl_dist",
           seed_mode = "static_spread",
           km_iter = 30,
           em_iter = 30,
           verbose = F)
pr <- predict_GMM(data,
                  gmm$centroids,
                  gmm$covariance_matrices,
                  gmm$weights)
names(pr$cluster_labels) <- rownames(data)
clustersTree <- sort(pr$cluster_labels)
gaps <- which((clustersTree[-1] - clustersTree[-length(clustersTree)]) != 0)

gene.cluster <- as.data.frame(clustersTree)
df.cluster <- cbind(gene.cluster,data[rownames(gene.cluster),])
data_clusters <- cbind(gene.cluster,data[rownames(gene.cluster),])
write.table(data_clusters,
            file = paste("WT_HS_vs_WT_RT.DEG_list.cluster",clusterN,".gmm.plotOrder.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

my_gene_clusters <- data.frame(cluster = df.cluster[,1])
rownames(my_gene_clusters) <- rownames(df.cluster)
colorsVec <- c(brewer.pal(12, "Set3"),
               brewer.pal(8, "Accent"))
names(colorsVec) <- c(1:clusterN)
my_colors <- list(cluster = colorsVec)

pdf(paste("WT_HS_vs_WT_RT.DEG_list.cluster",clusterN,".gmm.plotOrder.pdf",sep = ""), width = 2, height = 8)
phm <- pheatmap(df.cluster[,-1],
                color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                border_color = "snow",
                # scale = "row",
                cluster_rows = F,
                cluster_cols = F,                
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                labels_col = colnames(df.cluster[,-1]),
                main = paste("n = ",length(rownames(df.cluster)), sep = ""),
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)
dev.off()

# WT clusters as ref, other samples follow the order of WT clusters:
# qK_37C_0min/Col_0_22C_0min
# qK_37C_5min/Col_0_22C_5min
# qK_37C_15min/Col_0_22C_15min

# prd <- df.all[,c(grep("PrD",colnames(df.all),value = T))]
qk <- df.all[,c(grep("217R|218R|219R",grep("QK",colnames(df.all),value = T),value = T))]

qk.mean <- data.frame(
  qk_22C = (qk$QK_22C_218R + qk$QK_22C_rep2_219R)/2,
  qk_37C5min = (qk$QK_37C5min_218R + qk$QK_37C5min_rep2_219R)/2,
  qk_37C15min = (qk$QK_37C15min_218R + qk$QK_37C15min_rep2_219R)/2
)

rownames(qk.mean) <- rownames(qk)
qk.mean <- qk.mean[as.vector(DEG_list$V1),]
qk.mean <- qk.mean[grep("AT",rownames(df.cluster),value = T),]
qk.mean <- log2(qk.mean + 1)

wt.mean.clustered <- wt.mean[rownames(df.cluster),]

df.show.qk <- data.frame(qk_37C5min_vs_qk_22C = qk.mean$qk_37C5min - qk.mean$qk_22C,
                         qk_37C15min_vs_qk_22C = qk.mean$qk_37C15min - qk.mean$qk_22C)

head(df.show.qk)
head(df.cluster)

rownames(df.show.qk) <- rownames(qk.mean)

df.cluster.allsamples <- cbind(df.cluster,
                               df.show.qk)

my_gene_clusters <- data.frame(cluster = df.cluster[,1])
rownames(my_gene_clusters) <- rownames(df.cluster)

#apply(df.cluster.allsamples,1,)>2
df.cluster.allsamples[df.cluster.allsamples > 2] <- 2
df.cluster.allsamples[df.cluster.allsamples < -2] <- -2

colnames(df.cluster.allsamples[,-1])

my_gaps_col <- c(3)

pdf(paste("Col0_37C_vs_Col0_22C_as_ref.qk_37C_vs_qk_22C",clusterN,".gmm.pdf",sep = ""),
    height = 8,
    width = 3)
phm <- pheatmap(df.cluster.allsamples[,-1],
                color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                border_color = "snow",
                # scale = "row",
                cluster_rows = F,
                cluster_cols = F,                
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                gaps_col = my_gaps_col,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                labels_col = colnames(df.cluster.allsamples[,-1]),
                main = paste("n = ",length(rownames(df.cluster.allsamples)), sep = ""),
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)
dev.off()


# mark out A1 class bound genes
HsfA1ClassBoundGenes <- read.table("C:/project-hsf/clustering_RNA-seq/genesBoundByatLeast2HSFs/HsfA1_class_bound_confident_genes.geneID", header = F, sep = "\t")
rownames(HsfA1ClassBoundGenes) <- HsfA1ClassBoundGenes$V1

head(HsfA1ClassBoundGenes)
head(df.cluster)

colorsVec <- brewer.pal(8, "Set3")
names(colorsVec) <- c(1:8)
my_colors <- list(cluster = colorsVec)

df.cluster.allsamples$HSFA1Class_bound[rownames(df.cluster.allsamples) %in% rownames(HsfA1ClassBoundGenes)] <- 1
df.cluster.allsamples$HSFA1Class_bound[!(rownames(df.cluster.allsamples) %in% rownames(HsfA1ClassBoundGenes))] <- 0

df.cluster.allsamples$clustersTree <- df.cluster$clustersTree

write.table(df.cluster.allsamples,
            file = paste("WT_HS_vs_WT_RT.DEG_list.cluster8.gmm.plotOrder.HSFA1Class_bound.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

df.cluster.allsamples.sort <- df.cluster.allsamples[order(df.cluster.allsamples$clustersTree,df.cluster.allsamples$HSFA1Class_bound,decreasing = c(F,F)),]

write.table(df.cluster.allsamples.sort,
            file = paste("WT_HS_vs_WT_RT.DEG_list.cluster8.gmm.plotOrder.sort_by_HSFA1Class_bound.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

my_gaps_col <- c(3)

my_gene_clusters <- data.frame(cluster = df.cluster.allsamples.sort[,1],
                               bound = df.cluster.allsamples.sort[,7])
rownames(my_gene_clusters) <- rownames(df.cluster.allsamples.sort)


# clusterColors <- c(brewer.pal(12, "Set3"),
#                    brewer.pal(8, "Accent"))
clusterColors <- c(rep(c(brewer.pal(12, "Set3")[9]),8))
names(clusterColors) <- c(1:clusterN)
# boundColors <- c("grey60","indianred2")
boundColors <- c(brewer.pal(12, "Set3")[c(9,4)])
names(boundColors) <- c(0,1)
my_colors <- list(cluster = clusterColors,
                  bound = boundColors)

pdf(paste("Col0_37C_vs_Col0_22C_as_ref.bound_by_at_least_2_HSFA1.qk_37C_vs_qk_22C",clusterN,".gmm.pdf",sep = ""),
    height = 8,
    width = 3)
phm <- pheatmap(df.cluster.allsamples.sort[,c(2:6)],
                color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                border_color = "snow",
                # scale = "row",
                cluster_rows = F,
                cluster_cols = F,                
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                gaps_col = my_gaps_col,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                # labels_col = colnames(df.cluster.allsamples.sort[,c(7,2:6)]),
                labels_col = c("5min","15min","30min",
                               "5min","15min"),
                main = paste("n = ",length(rownames(df.cluster.allsamples.sort)), sep = ""),
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)
dev.off()