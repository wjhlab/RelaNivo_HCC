#Myeloid sub clusters=======
data_mye <- data_full_scaled[data_full_scaled$cluster1m=="Myeloid",]
submarkers <- c("CD68", "PDPN", "VIM", "CK", "CD138", "HEPPAR", "CD31", "CD206", "CD163",
                "PDL1", "FGL1", "CD15", "PNAd", "DCLAMP", "CD45", "KI67", "HLADR", "CD86")


data_mye_expr <- data.matrix(data_mye[,submarkers])


#flowSOM
som_submye <- ReadInput(data_mye_expr, transform = FALSE, scale = F) %>% BuildSOM(colsToUse = submarkers)

## Get the cell clustering into 100 SOM codes
cell_clustering_som_submye <- som_submye$map$mapping[,1]

numclusters=50
plottitle='consensus_plots'
## Metaclustering into numclusters with ConsensusClusterPlus
codes_submye <- som_submye$map$codes
mc_submye <- ConsensusClusterPlus(t(codes_submye), maxK = numclusters, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plottitle, 
                           plot = "png", clusterAlg = "hc", 
                           innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering_submye <- mc_submye[[numclusters]]$consensusClass#metaclusters consensus
cell_clustering_submye <- code_clustering_submye[cell_clustering_som_submye]#cell clustering from som



if(length(cell_clustering_submye) != nrow(data_mye)) {
  stop("Mismatch: clustering length doesn't match data rows!")
}

#Add back in metacluster assignment
data_mye$myeloid_cluster <- cell_clustering_submye

#SAVE RDS file of all relevant dataframes and matrices
flowSOM_clusterfcs_submye <- list(data_mye_expr, codes_submye, mc_submye, code_clustering_submye, cell_clustering_submye, cell_clustering_som_submye)
saveRDS(flowSOM_clusterfcs_submye, "flowSOM_clusterfcs_myeloidsubcluster.RDS")

# Get counts for clusters
counts <- as.data.frame(table(data_mye$myeloid_cluster))
colnames(counts) <- c("myeloid_cluster", "n_count")
counts$myeloid_cluster <- as.integer(as.character(counts$myeloid_cluster))

ordered_counts <- counts
counts <- counts[order(counts$n_count, decreasing = TRUE), ]
print(counts)

write.xlsx(ordered_counts, file = "flowSOM_myeloidsubcluster_counts.xlsx", rowNames = FALSE)



# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_mye)[, submarkers],
  by = list(cluster = as.integer(data_mye$myeloid_cluster)),
  FUN = mean,
  na.rm = TRUE
) %>%
  arrange(cluster) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

cluster_scaled<-t(scale(t(metacluster_means)))
rownames(cluster_scaled)<-rownames(metacluster_means)

# Add metacluster size annotation (SAME SORT ORDER)
metacluster_size <- as.numeric(
  table(data_mye$myeloid_cluster)[
    as.character(sort(as.integer(unique(data_mye$myeloid_cluster))))
  ]
)
# Check alignment
print(data.frame(
  cluster = as.integer(rownames(metacluster_means)),
  size = metacluster_size
))

# ===== ANNOTATION SETUP (matching original code exactly) =====
annotation_row <- data.frame(Cluster = factor(rownames(metacluster_means)))
rownames(annotation_row) <- rownames(metacluster_means)


# Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned <- kovesi.rainbow_bgyrm_35_85_c69(length(unique(rownames(metacluster_means))))
names(colorassigned) <- sort(unique(as.integer(rownames(metacluster_means))))
color_list = list(clusters = colorassigned)
color_list_byoriginal = colorassigned[match((rownames(metacluster_means)), names(colorassigned))]


#Cluster annotation
rAbar <- rowAnnotation(
  clusters = factor(rownames(metacluster_means), levels = sort(unique(as.integer(rownames(metacluster_means))))),
  col = color_list,
  gp = gpar(col = "white", lwd = 0.5),
  counts = anno_barplot(
    metacluster_size,  
    gp = gpar(fill = colorassigned),
    border = F,
    bar_width = 0.75,
    width = unit(2, "cm")
  )
)

# ===== CREATE HEATMAP =====
pdf("MyeloidSubcluster_UnannotatedFlowSOM_Heatmap.pdf", width = 10, height = 14)
Heatmap(cluster_scaled,
        column_title = "Myeloid Subcluster FlowSOM Metaclusters",
        name = "scaled",
        col = rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,  
        border = NA,
        rect_gp = gpar(col = "white", lwd = 0.5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(at = seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled) * unit(4, "mm"),
        height = nrow(cluster_scaled) * unit(4, "mm"))
dev.off()


#save as RDS file
mye_data <- list(data_mye, data_mye_expr)
saveRDS(mye_data, "Myeloid_Subcluster.RDS")
mye_data <- readRDS("Myeloid_Subcluster.RDS")
data_mye <- data.frame(mye_data[1])
data_mye_expr <- data.frame(mye_data[2])
data_mye_expr01 <- data.frame(mye_data[3])




# Generate density box plots
counts_table_m<-table(data_mye$myeloid_cluster, data_mye$sample_id)
props_table_m <- t(t(counts_table_m) / colSums(counts_table_m)) * 100
counts_m <- as.data.frame.matrix(counts_table_m)
props_m <-as.data.frame.matrix(props_table_m)

areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities_m <- t(t(counts_m)/areas$TotalArea)

write.xlsx(densities_m, file = "myeloidsubcluster_densities.xlsx", rowNames = FALSE)
#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_m), densities_m, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$roi <- factor(md$ROI[match(ggdfd$sample_id,md$sample_id)], levels=ROIlevels)
ggdfd$case <- factor(md$Case[match(ggdfd$sample_id,md$sample_id)], levels=caselevels)
ggdfd$tx <- factor(md$Tx[match(ggdfd$sample_id,md$sample_id)], levels=txlevels)
ggdfd$response <- factor(md$Response[match(ggdfd$sample_id,md$sample_id)], levels=resplevels)
ggdfd$presponse <- factor(md$pResponse[match(ggdfd$sample_id,md$sample_id)], levels=presplevels)


#plot box plots
# CELLS
ggp2<-ggplot(ggdfd,aes(x=case,y=densities,fill=case))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("Density (# of total area)")+
  ggtitle("J20121 Case - flowSOM Myeloid Subclusters")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=4, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11)
  )
pdf('MyeloidSubcluster_UnannotatedFlowSOM_Density_box_case.pdf',width=20,height=9)
ggp2
dev.off()


# Generate density box plots
counts_table_m<-table(data_full_scaled$cluster, data_full_scaled$sample_id)
props_table_m <- t(t(counts_table_m) / colSums(counts_table_m)) * 100
counts_m <- as.data.frame.matrix(counts_table_m)
props_m <-as.data.frame.matrix(props_table_m)

areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities_m <- t(t(counts_m)/areas$TotalArea)


#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_m), densities_m, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$roi <- factor(md$ROI[match(ggdfd$sample_id,md$sample_id)], levels=ROIlevels)
ggdfd$case <- factor(md$Case[match(ggdfd$sample_id,md$sample_id)], levels=caselevels)
ggdfd$tx <- factor(md$Tx[match(ggdfd$sample_id,md$sample_id)], levels=txlevels)
ggdfd$response <- factor(md$Response[match(ggdfd$sample_id,md$sample_id)], levels=resplevels)
ggdfd$presponse <- factor(md$pResponse[match(ggdfd$sample_id,md$sample_id)], levels=presplevels)


#plot box plots
# CELLS
ggp2<-ggplot(ggdfd,aes(x=case,y=densities,fill=case, shape=roi))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("Density (% of total area)")+
  ggtitle("J20121 Case - Unannotated Direct flowSOM Clusters")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=4, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11)
  )
pdf('MyeloidSubcluster_UnannotatedFlowSOM_Density_box_CaseByROI.pdf',width=20,height=9)
ggp2
dev.off()



#######CLUSTER ANNOTATIONS
#annotate each cluster based on expression, add in annotated clusters + respective cluster number in merge.xlsx ===========
clusterMergeFile = paste0(work,"/Config/merge_myeloid.xlsx") #create a column of clusters numbered + another column of respective annotation
cluster_merging <- read_excel(clusterMergeFile)


#Annotations based on cluster heatmap
clusterlevels=c("DC",
                "Gran",
                "Mac_I",
                "Mac_II",
                "Mac_III",
                "Mac_IV",
                "Mac_V",
                "Mac_VI",
                "Mac_VII",
                "Prolif_Mac",
                "Stroma",
                "Hep",
                "UA")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_mye$myeloid_cluster, cluster_merging$original_cluster)
data_mye$myeloid_cluster1m <- cluster_merging$new_cluster[mm1]


# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_mye)[, submarkers],
  by = list(cluster = data_mye$myeloid_cluster1m),
  FUN = mean,
  na.rm = TRUE
) %>%
  arrange(cluster) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

metacluster_means <- metacluster_means[clusterlevels, ]

cluster_scaled<-t(scale(t(metacluster_means)))
rownames(cluster_scaled)<-rownames(metacluster_means)

# Add metacluster size annotation (SAME SORT ORDER)
metacluster_size <- as.numeric(
  table(data_mye$myeloid_cluster1m)[clusterlevels]
)

# Check alignment
print(data.frame(
  cluster = rownames(metacluster_means),
  size = metacluster_size
))

# ===== ANNOTATION SETUP (matching original code exactly) =====

legend_breaks = seq(from = 0, to = 1, by = 0.2)
clusternames<-clusterlevels
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))
names(colorassigned)<-clusternames
color_list = list(clusters=colorassigned)


#Cluster annotation
rAbar <- rowAnnotation(
  clusters = factor(rownames(metacluster_means), levels = clusterlevels),
  col = color_list,
  gp = gpar(col = "white", lwd = 0.5),
  counts = anno_barplot(
    metacluster_size,  
    gp = gpar(fill = colorassigned),
    border = F,
    bar_width = 0.75,
    width = unit(2, "cm")
  )
)

# ===== CREATE HEATMAP =====
pdf("MyeloidSubcluster_MergedFlowSOM_Heatmap.pdf", width = 10, height = 8)
Heatmap(cluster_scaled,
        column_title = "FlowSOM Myeloid Merged Metaclusters",
        name = "scaled",
        col = rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,  
        border = NA,
        rect_gp = gpar(col = "white", lwd = 0.5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(at = seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled) * unit(4, "mm"),
        height = nrow(cluster_scaled) * unit(4, "mm"))
dev.off()

#save as RDS file
mye_data <- list(data_mye, data_mye_expr)
saveRDS(mye_data, "Myeloid_Subcluster.RDS")



#read as RDS file
mye_data <- readRDS("Myeloid_Subcluster.RDS")
data_mye <- data.frame(mye_data[1])
data_mye_expr <- data.frame(mye_data[2])






#Generate density plots for broad mac clusters
#Group by Broad Annotations
# Create mapping for cluster1m to cluster_broad
cluster_broad_mapping <- c("DC" = "DC",
                "Gran" = "Gran",
                "Mac_I" = "Mac",
                "Mac_II" = "Mac",
                "Mac_III" = "Mac",
                "Mac_IV" = "Mac",
                "Mac_V" = "Mac",
                "Mac_VI" = "Mac",
                "Mac_VII" = "Mac",
                "Prolif_Mac" = "Mac",
                "Stroma" = "Stroma",
                "Hep" = "Hep",
                "UA" = "UA")

# Add cluster_broad todata_mye
data_mye$cluster_broad <- cluster_broad_mapping[data_mye$myeloid_cluster1m]

#Broad Annotations
counts_table_b<-table(data_mye$cluster_broad, data_mye$sample_id)
props_table_b <- t(t(counts_table_b) / colSums(counts_table_b)) * 100
counts_b <- as.data.frame.matrix(counts_table_b)
props_b <-as.data.frame.matrix(props_table_b)

raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities_broad <- t(t(counts_b)/areas$TotalArea)


#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_broad), densities_broad, check.names = FALSE),
              id.vars = "cluster", value.name = "densities_b", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$roi <- factor(md$ROI[match(ggdfd$sample_id,md$sample_id)], levels=ROIlevels)
ggdfd$case <- factor(md$Case[match(ggdfd$sample_id,md$sample_id)], levels=caselevels)
ggdfd$tx <- factor(md$Tx[match(ggdfd$sample_id,md$sample_id)], levels=txlevels)
ggdfd$response <- factor(md$Response[match(ggdfd$sample_id,md$sample_id)], levels=resplevels)
ggdfd$presponse <- factor(md$pResponse[match(ggdfd$sample_id,md$sample_id)], levels=presplevels)


ggp2<-ggplot(ggdfd,aes(x=case,y=densities_b,fill=case))+
  geom_boxplot(outlier.shape=NA, lwd=0.25)+
  geom_jitter(width=0.2)+
  facet_wrap(~cluster,ncol=3,scales="free")+
  ylab("Density (# of cells)")+
  ggtitle("J20121 Density Merged")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=3, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 10)
  )
pdf('MyeloidSubcluster_MergedflowSOM_MacBroadClusters_Density_box_case.pdf',width=14,height=8)
ggp2
dev.off()



#Detailed Annotations
counts_table_m <-table(data_mye$myeloid_cluster1m, data_mye$sample_id)
props_table_m <- t(t(counts_table_m) / colSums(counts_table_m)) * 100
counts_m <- as.data.frame.matrix(counts_table_m)
props_m <-as.data.frame.matrix(props_table_m)

raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities_m <- t(t(counts_m)/areas$TotalArea)


#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_m), densities_m, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$roi <- factor(md$ROI[match(ggdfd$sample_id,md$sample_id)], levels=ROIlevels)
ggdfd$case <- factor(md$Case[match(ggdfd$sample_id,md$sample_id)], levels=caselevels)
ggdfd$tx <- factor(md$Tx[match(ggdfd$sample_id,md$sample_id)], levels=txlevels)
ggdfd$response <- factor(md$Response[match(ggdfd$sample_id,md$sample_id)], levels=resplevels)
ggdfd$presponse <- factor(md$pResponse[match(ggdfd$sample_id,md$sample_id)], levels=presplevels)


ggp2<-ggplot(ggdfd,aes(x=case,y=densities,fill=case))+
  geom_boxplot(outlier.shape=NA, lwd=0.25)+
  geom_jitter(width=0.2)+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("Density (# of cells)")+
  ggtitle("J20121 Density Merged")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=3, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 10)
  )
pdf('MyeloidSubcluster_MergedflowSOM_Density_box_case.pdf',width=14,height=10)
ggp2
dev.off()
