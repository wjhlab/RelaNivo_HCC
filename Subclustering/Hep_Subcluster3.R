#Hep sub clusters=======
data_hep <- data_full_scaled[data_full_scaled$cluster1m=="Hep",]
submarkers <- c("SMA","PDPN", "VIM", "CK", "CD138", "HEPPAR", "CD31",
                "PDL1", "FGL1", "PNAd", "AID", "KI67", "HLADR")


data_hep_expr <- data.matrix(data_hep[,submarkers])


#flowSOM
som_subhep <- ReadInput(data_hep_expr, transform = FALSE, scale = F) %>% BuildSOM(colsToUse = submarkers)

## Get the cell clustering into 100 SOM codes
cell_clustering_som_subhep <- som_subhep$map$mapping[,1]

numclusters=20
plottitle='consensus_plots'
## Metaclustering into numclusters with ConsensusClusterPlus
codes_subhep <- som_subhep$map$codes
mc_subhep <- ConsensusClusterPlus(t(codes_subhep), maxK = numclusters, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plottitle, 
                           plot = "png", clusterAlg = "hc", 
                           innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering_subhep <- mc_subhep[[numclusters]]$consensusClass#metaclusters consensus
cell_clustering_subhep <- code_clustering_subhep[cell_clustering_som_subhep]#cell clustering from som



if(length(cell_clustering_subhep) != nrow(data_hep)) {
  stop("Mismatch: clustering length doesn't match data rows!")
}

#Add back in metacluster assignment
data_hep$hep_cluster <- cell_clustering_subhep

#SAVE RDS file of all relevant dataframes and matrices
flowSOM_clusterfcs_subhep <- list(data_hep_expr, codes_subhep, mc_subhep, code_clustering_subhep, cell_clustering_subhep, cell_clustering_som_subhep)
saveRDS(flowSOM_clusterfcs_subhep, "flowSOM_clusterfcs_hepsubcluster.RDS")

# Get counts for clusters
counts <- as.data.frame(table(data_hep$hep_cluster))
colnames(counts) <- c("hep_cluster", "n_count")
counts$hep_cluster <- as.integer(as.character(counts$hep_cluster))

ordered_counts <- counts
counts <- counts[order(counts$n_count, decreasing = TRUE), ]
print(counts)

write.xlsx(ordered_counts, file = "flowSOM_hepsubcluster_counts.xlsx", rowNames = FALSE)



# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_hep)[, submarkers],
  by = list(cluster = as.integer(data_hep$hep_cluster)),
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
  table(data_hep$hep_cluster)[
    as.character(sort(as.integer(unique(data_hep$hep_cluster))))
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
pdf("HepSubcluster_UnannotatedFlowSOM_Heatmap.pdf", width = 10, height = 14)
Heatmap(cluster_scaled,
        column_title = "Hep Subcluster FlowSOM Metaclusters",
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
hep_data <- list(data_hep, data_hep_expr)
saveRDS(hep_data, "Hep_Subcluster.RDS")
hep_data <- readRDS("Hep_Subcluster.RDS")
data_hep <- data.frame(hep_data[1])
data_hep_expr <- data.frame(hep_data[2])
data_hep_expr01 <- data.frame(hep_data[3])




# Generate density box plots
counts_table_m<-table(data_hep$hep_cluster, data_hep$sample_id)
props_table_m <- t(t(counts_table_m) / colSums(counts_table_m)) * 100
counts_m <- as.data.frame.matrix(counts_table_m)
props_m <-as.data.frame.matrix(props_table_m)

areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities_m <- t(t(counts_m)/areas$TotalArea)

write.xlsx(densities_m, file = "hepsubcluster_densities.xlsx", rowNames = FALSE)
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
  ggtitle("J20121 Case - flowSOM Hep Subclusters")+
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
pdf('HepSubcluster_UnannotatedFlowSOM_Density_box_case.pdf',width=20,height=9)
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


#######CLUSTER ANNOTATIONS
#annotate each cluster based on expression, add in annotated clusters + respective cluster number in merge.xlsx ===========
clusterMergeFile = paste0(work,"/Config/merge_hep.xlsx") #create a column of clusters numbered + another column of respective annotation
cluster_merging <- read_excel(clusterMergeFile)


#Annotations based on cluster heatmap
clusterlevels=c("Hep_I",
                "Hep_II",
                "Hep_III",
                "Hep_IV",
                "Hep_V",
                "Hep_VI",
                "Prolif_Hep",
                "Stroma")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_hep$hep_cluster, cluster_merging$original_cluster)
data_hep$hep_cluster1m <- cluster_merging$new_cluster[mm1]


# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_hep)[, submarkers],
  by = list(cluster = data_hep$hep_cluster1m),
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
  table(data_hep$hep_cluster1m)[clusterlevels]
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
pdf("HepSubcluster_MergedFlowSOM_Heatmap.pdf", width = 10, height = 8)
Heatmap(cluster_scaled,
        column_title = "FlowSOM Hep Merged Metaclusters",
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
data_hep <- list(data_hep, data_hep_expr)
saveRDS(data_hep, "Hep_Subcluster.RDS")



#Read as RDS file
hep_data <- readRDS("Hep_Subcluster.RDS")
data_hep <- data.frame(hep_data[1])
data_hep_expr <- data.frame(hep_data[2])


#Annotations
counts_table_m <-table(data_hep$hep_cluster1m, data_hep$sample_id)
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
  facet_wrap(~cluster,ncol=7,scales="free")+
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
pdf('HepSubcluster_MergedFlowSOM_Density_box_case.pdf',width=20,height=7.5)
ggp2
dev.off()
