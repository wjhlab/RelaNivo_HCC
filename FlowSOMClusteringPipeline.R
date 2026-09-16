# WJH Lab
# Sarah Shin
# Version R 4.5.3

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

dataDirectory = paste0(work,"/Data")


## QC Check - Make sure all files in metadata are present in datadirectory
if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])){
  print(paste('ERR: not all filenames in metadata present in data folder - missing',
              md$file_name[!which(data$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),
                                                                                   pattern = '.csv')])],'Subsetting...'))
  md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])),]
}


## Read csv into csv_raw =========
file_paths <- paste0(dataDirectory, "/", md$file_name)
csv_list <- lapply(seq_along(file_paths), function(i) {
  df <- read.csv(file_paths[i], stringsAsFactors = FALSE)
  df$source_file <- basename(file_paths[i])   # record which csv this row came from
  df$sample_id_from_metadata <- as.character(md$sample_id[i]) # attach sample id from metadata
  df
})

csv_raw_full <- plyr::rbind.fill(csv_list)


#export raw channel names to clean/rename panel
#Ignore 'Error in colSums(csv_raw_full): 'x' must be numeric
rawcolnames <- c()
rawcolnames$name <- colnames(csv_raw_full)
rawcolnames$sum <- colSums(csv_raw_full)
#write.csv(rawcolnames, 'rawpanel.csv')

#clean panel by creating new 'cleanpanel.xlsx' file and adding new columns for cleaned up marker names + selecting which markers to use for clustering/subclustering
cleanpanel <- read_xlsx('cleanpanel.xlsx')
colnames(csv_raw_full) <- cleanpanel$clean_names
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]
csv_full <- csv_raw_full[,colnames(csv_raw_full) %in% panel]
csv_full <- csv_full %>% select(sample_id, everything())

#Once cleanpanel.xlsx is generated and cleaned, run this line
panelDataFile = paste0(work,"/cleanpanel.xlsx")

#sort panels into different categories
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]


#Add in case assignment for scaling by case======
data_full <- as.data.frame(csv_full)
data_full <- merge(data_full, md[, c("sample_id", "Case")], by = "sample_id", all.x = TRUE)
data_full_scaled <- data_full

# QC: Check unique row numbers per patient
patient_rows <- split(rownames(data_full), data_full$Case)
cat("Rows per patient:\n")
print(sapply(patient_rows, length))

# Save patient-row mapping
patient_row_mapping <- data.frame(
  patient_id = names(patient_rows),
  n_rows = sapply(patient_rows, length),
  row_indices = sapply(patient_rows, function(x) paste(x, collapse = ","))
)

# ===== SCALING BY PATIENT =====
for (patient in unique(data_full$Case)) {
  patient_idx <- which(data_full$Case == patient)
  
  # Extract patient data
  patient_data <- data.matrix(data_full[patient_idx, fullmarkerlist])
  patient_data <- asinh(patient_data / 0.8)
  
  # Calculate quantiles for this patient only
  rng <- colQuantiles(patient_data, probs = c(0.01, 0.99))
  
  # Scale to 0-1 for this patient
  data01_patient <- t((t(patient_data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  data01_patient[data01_patient < 0] <- 0
  data01_patient[data01_patient > 1] <- 1
  
  # Update scaled data in dataframe
  data_full_scaled[patient_idx, fullmarkerlist] <- data01_patient
}

#saveRDS(data_full_scaled, "data_full_scaled.RDS")

# ===== FLOWSOM ANALYSIS =====
set.seed(1234)
marker_cols <- c("CD68", "PDPN", "VIM", "CK", "CD138", "HEPPAR", "CD31",  
                 "FGL1", "CD3", "CD4", "CD8", "CD20", "CD15", "CD57", "PNAd",
                 "DCLAMP", "CD45")
expr <- as.matrix(data_full_scaled[marker_cols])
som <- ReadInput(expr, transform = FALSE, scale = T) %>% BuildSOM(colsToUse = marker_cols)

## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]

numclusters=50
plottitle='consensus_plots'
## Metaclustering into numclusters with ConsensusClusterPlus
codes <- som$map$codes
mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plottitle, 
                           plot = "png", clusterAlg = "hc", 
                           innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som



if(length(cell_clustering) != nrow(data_full_scaled)) {
  stop("Mismatch: clustering length doesn't match data rows!")
}

#Add back in metacluster assignment
data_full_scaled$cluster <- cell_clustering

# Get counts for clusters
counts <- as.data.frame(table(data_full_scaled$cluster))
colnames(counts) <- c("cluster", "n_count")
counts$cluster <- as.integer(as.character(counts$cluster))

ordered_counts <- counts
counts <- counts[order(counts$n_count, decreasing = TRUE), ]
print(counts)

write.xlsx(ordered_counts, file = "UnannotatedflowSOMscaled_counts.xlsx", rowNames = FALSE)

data_full_scaled<-as.data.frame(data_full_scaled)

# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_full_scaled)[, marker_cols],
  by = list(cluster = as.integer(data_full_scaled$cluster)),
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
  table(data_full_scaled$cluster)[
    as.character(sort(as.integer(unique(data_full_scaled$cluster))))
  ]
)
# Check alignment
print(data.frame(
  cluster = as.integer(rownames(metacluster_means)),
  size = metacluster_size
))

# ===== ANNOTATION SETUP =====
annotation_row <- data.frame(Cluster = factor(rownames(metacluster_means)))
rownames(annotation_row) <- rownames(metacluster_means)

color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

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
pdf("UnannotatedFlowSOMScaled_Heatmap.pdf", width = 10, height = 14)
Heatmap(cluster_scaled,
        column_title = "Original FlowSOM Scaled Metaclusters",
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


#update saved RDS file - this RDS has all subtyping markers re-clustered by flowSOM
global_data <- data_full_scaled
saveRDS(global_data, "data_full.RDS")


#######CLUSTER ANNOTATIONS
#annotate each cluster based on expression, add in annotated clusters + respective cluster number in merge.xlsx ===========
clusterMergeFile = paste0(work,"/Config/merge.xlsx") #create a column of clusters numbered + another column of respective annotation
cluster_merging <- read_excel(clusterMergeFile)


#Annotate first as broad cell classifications on cluster heatmap
clusterlevels=c("Lymphoid",
                "Myeloid",
                "Hep",
                "Endoth",
                "Stroma",
                "UA")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
names(colorassigned)<-clusterlevels
mm1 <- match(data_full_scaled$cluster, cluster_merging$original_cluster)
data_full_scaled$cluster1m <- cluster_merging$new_cluster[mm1]
 

# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_full_scaled)[, marker_cols],
  by = list(cluster = data_full_scaled$cluster1m),
  FUN = mean,
  na.rm = TRUE
) %>%
  arrange(cluster) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

broadcluster_scaled<-t(scale(t(metacluster_means)))
rownames(broadcluster_scaled)<-rownames(metacluster_means)

# Add metacluster size annotation (SAME SORT ORDER)
metacluster_size <- as.numeric(
  table(data_full_scaled$cluster1m)[clusterlevels]
)
  
# Check alignment
print(data.frame(
  cluster = rownames(metacluster_means),
  size = metacluster_size
))

# ===== ANNOTATION SETUP =====
annotation_row <- data.frame(Cluster = factor(rownames(metacluster_means)))
rownames(annotation_row) <- rownames(metacluster_means)

if(!is.null(clusterMergeFile)){
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
  annotation_row$Merged <- cluster_merging$new_cluster
  color_clusters2 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Merged))
  names(color_clusters2) <- levels(cluster_merging$new_cluster)
  annotation_colors$Merged <- color_clusters2
}

legend_breaks = seq(from = 0, to = 1, by = 0.2)
clusternames<-clusterlevels
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))
names(colorassigned)<-clusternames
color_list = list(clusters=colorassigned)


metacluster_means <- metacluster_means[clusterlevels, ]
broadcluster_scaled <- broadcluster_scaled[clusterlevels,]
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
pdf("InitialBroadClustersFlowSOMScaled_Heatmap.pdf", width = 10, height = 8)
Heatmap(broadcluster_scaled,
        column_title = "FlowSOM Broad Metaclusters",
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
        heatmap_legend_param = list(at = seq(from = round(min(broadcluster_scaled)), to = round(max(broadcluster_scaled)))),
        width = ncol(broadcluster_scaled) * unit(4, "mm"),
        height = nrow(broadcluster_scaled) * unit(4, "mm"))
dev.off()



# Go to subclustering for broad cluster population --> Lymphoid_Subcluster.R --> Myeloid_Subcluster2.R
# --> Hep_Subcluster3.R --> Stroma_Subcluster4.R --> Endoth_Subcluster5.R
# After first round of final clusters - and differential plots + final merged heatmap - we reviewed cluster names
# and updated the naming of certain cluster names (updatedclusternamesuponreview.xlsx) to cluster1m (in updateclusternames.R) 
# and saved updated dataframe as data_full.RDS <- this is now the MOST UPDATED dataframe

data_full <- readRDS("./Subclustering/data_full.RDS")
unique(data_full$cluster1m)

# Calculate mean expression per metacluster 
metacluster_means <- aggregate(
  as.data.frame(data_full)[, fullmarkerlist],
  by = list(cluster = data_full$cluster1m),
  FUN = mean,
  na.rm = TRUE
) %>%
  arrange(cluster) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

clusterlevels=c("B_I","B_II","B_III","Plasma",
                "Naive_Tc","Tc_I","Act_Tc_I","Act_Tc_II","Act_Tc_III","Act_Tc_IV","Act_Prolif_Tc","Mem_Tc_I","Mem_Tc_II",
                "Mem_Th_I","Mem_Th_II","Mem_Th_III","Mem_Prolif_Th","Treg",
                "Act_NK","NKT",
                "DC","Gran",
                "Mac_I","Mac_II","Mac_III","Mac_IV","Mac_V","Mac_VI","Mac_VII","Prolif_Mac",
                "Endoth", "HEV",
                "Stroma_I","Stroma_II","Stroma_III",
                "Hep_I","Hep_II", "Hep_III","Hep_IV","Hep_V","Hep_VI","Prolif_Hep",
                "UA")
                
metacluster_means <- metacluster_means[clusterlevels, ]

cluster_scaled<-t(scale(t(metacluster_means)))
rownames(cluster_scaled)<-rownames(metacluster_means)

# Add metacluster size annotation (SAME SORT ORDER)
metacluster_size <- as.numeric(
  table(data_full$cluster1m)[clusterlevels]
)

# Check alignment
print(data.frame(
  cluster = rownames(metacluster_means),
  size = metacluster_size
))

# ===== ANNOTATION SETUP (matching original code exactly) =====
#metacluster_means <- metacluster_means[clusterlevels, ]
annotation_row <- data.frame(Cluster = factor(rownames(metacluster_means)))
rownames(annotation_row) <- rownames(metacluster_means)

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

##Figure 3g =====
# ===== CREATE HEATMAP =====
pdf("0FINAL_AllSubclusters_flowSOM_heatmap_MERGED.pdf", width = 10, height = 12)
Heatmap(cluster_scaled,
        column_title = "Final FlowSOM Metaclusters",
        name = "scaled",
        col = rev(brewer.rdbu(100)),
        cluster_columns = F,
        cluster_rows = F,  
        border = NA,
        rect_gp = gpar(col = "white", lwd = 0.5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        heatmap_legend_param = list(at = seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled) * unit(4, "mm"),
        height = nrow(cluster_scaled) * unit(4, "mm"))
dev.off()


# Create Broad cateogries for all updated clusters

cluster_broad_mapping <- c("B_I" = "B",
                           "B_II" = "B",
                           "B_III" = "B",
                           "Plasma" = "Plasma",
                           "Naive_Tc" = "Tc",
                           "Tc_I" = "Tc",
                           "Act_Tc_I" = "Tc",
                           "Act_Tc_II" = "Tc",
                           "Act_Tc_III" = "Tc",
                           "Act_Tc_IV" = "Tc",
                           "Act_Prolif_Tc" = "Tc",
                           "Mem_Tc_I" = "Tc",
                           "Mem_Tc_II" = "Tc",
                           "Mem_Th_I" = "Th",
                           "Mem_Th_II" = "Th",
                           "Mem_Th_III" = "Th",
                           "Mem_Prolif_Th" = "Th",
                           "Treg" = "Treg",
                           "Act_NK" = "NK",
                           "NKT" = "NKT",
                           "DC" = "DC",
                           "Gran" = "Gran",
                           "Mac_I" = "Mac",
                           "Mac_II" = "Mac",
                           "Mac_III" = "Mac",
                           "Mac_IV" = "Mac",
                           "Mac_V" = "Mac",
                           "Mac_VI" = "Mac",
                           "Mac_VII" = "Mac",
                           "Prolif_Mac" = "Mac",
                           "Endoth" = "Endoth", 
                           "HEV" = "Endoth",
                           "Hep_I" = "Hep",
                           "Hep_II" = "Hep",
                           "Hep_III" = "Hep",
                           "Hep_IV" = "Hep",
                           "Hep_V" = "Hep",
                           "Hep_VI" = "Hep",
                           "Prolif_Hep" = "Hep",
                           "Stroma_I" = "Stroma",
                           "Stroma_II" = "Stroma",
                           "Stroma_III" = "Stroma",
                           "UA" = "UA")

broadclusterlevels=c("B", "Plasma", "Th", "Tc", "Treg",
                     "NKT", "NK", "DC", "Gran", "Mac",
                     "HEV", "Endoth", "Stroma", "Hep", "UA")


# Add cluster_broad to data_full_scaled
data_full$cluster_broad <- cluster_broad_mapping[data_full$cluster1m]
data_full$tx <- factor(md$Tx[match(data_full$sample_id,md$sample_id)], levels=txlevels)
data_full$response <- factor(md$Response[match(data_full$sample_id,md$sample_id)], levels=resplevels)
data_full$respgroup <- factor(md$Response_Grouped[match(data_full$sample_id,md$sample_id)], levels=respgrouplevels)
data_full$presponse <- factor(md$pResponse[match(data_full$sample_id,md$sample_id)], levels=presplevels)
data_full$prespgroup <- factor(md$pResponse_Grouped[match(data_full$sample_id,md$sample_id)], levels=prespgrouplevels)
data_full$presp_by_tx <- factor(md$pRespByTx[match(data_full$sample_id,md$sample_id)])

#saveRDS(data_full, "df_output.rds")

