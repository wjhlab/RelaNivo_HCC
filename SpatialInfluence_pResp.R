# Load function & data
df_output <- readRDS("./df_output.rds")
spatwt <- readRDS("./spatwt.RDS")

## CD8+T TME analysis====
# Tc
idx <- which(grepl("Tc", df_output$cluster1m))

spatwt_df_cd8T_2m<- spatwt[idx, ]
spatwt_df_cd8T_2m$sample_id <- factor(spatwt_df_cd8T_2m$sample_id, levels=samplevels)
spatwt_df_cd8T_2m$Case <- factor(md$Case[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_cd8T_2m$pRespByTx <- factor(md$pRespByTx[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)])

funcMarkers<- c("BCL6", "CD25", "TOX", "PD1", "GZMB", "KI67", "LAG3", "HLADR")

expr_cd8T<- as.data.frame(df_output[idx, funcMarkers])
expr_cd8T$sample_id<- spatwt_df_cd8T_2m$sample_id
expr_cd8T$Case<- spatwt_df_cd8T_2m$Case
expr_cd8T$pRespByTx <- spatwt_df_cd8T_2m$pRespByTx

do_CN_analysis(spatwt_df = spatwt_df_cd8T_2m,
               expr0 = expr_cd8T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "pRespByTx", 
               grouplevels = c("Nivo_NR", "Nivo_R", "Nivo_Rela_NR", "Nivo_Rela_R"),
               comparisons = list(c("Nivo_R", "Nivo_Rela_R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_Tc_15clus_pRespByTx_")


spatwt_df_cd8T_2m$prespgroup <- factor(md$pResponse_Grouped[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)], levels=prespgrouplevels)
expr_cd8T$prespgroup <- spatwt_df_cd8T_2m$prespgroup
do_CN_analysis(spatwt_df = spatwt_df_cd8T_2m,
               expr0 = expr_cd8T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "prespgroup", 
               grouplevels = c("NR", "R"),
               comparisons = list(c("NR", "R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_Tc_15clus_pRespGroup_")




## CD4+T TME analysis====
# Th
idx <- which(grepl("Th", df_output$cluster1m))

spatwt_df_cd4T_2m<- spatwt[idx, ]
spatwt_df_cd4T_2m$sample_id <- factor(spatwt_df_cd4T_2m$sample_id, levels=samplevels)
spatwt_df_cd4T_2m$Case <- factor(md$Case[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_cd4T_2m$pRespByTx <- factor(md$pRespByTx[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)])

funcMarkers<- c("BCL6", "CD25", "TOX", "PD1", "FOXP3", "KI67", "LAG3", "HLADR")
expr_cd4T<- as.data.frame(df_output[idx, funcMarkers])
expr_cd4T$sample_id<- spatwt_df_cd4T_2m$sample_id
expr_cd4T$Case<- spatwt_df_cd4T_2m$Case
expr_cd4T$pRespByTx <- spatwt_df_cd4T_2m$pRespByTx

do_CN_analysis(spatwt_df = spatwt_df_cd4T_2m,
               expr0 = expr_cd4T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "pRespByTx", 
               grouplevels = c("Nivo_NR", "Nivo_R", "Nivo_Rela_NR", "Nivo_Rela_R"),
               comparisons = list(c("Nivo_R", "Nivo_Rela_R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_Th_15clus_pRespByTx_")



spatwt_df_cd4T_2m$prespgroup <- factor(md$pResponse_Grouped[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)], levels=prespgrouplevels)
expr_cd4T$prespgroup <- spatwt_df_cd4T_2m$prespgroup
do_CN_analysis(spatwt_df = spatwt_df_cd4T_2m,
               expr0 = expr_cd4T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "prespgroup", 
               grouplevels = c("NR", "R"),
               comparisons = list(c("NR", "R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_Th_15clus_pRespGroup_")

## B TME analysis====
# B
idx <- which(grepl("B_", df_output$cluster1m))

spatwt_df_B_2m<- spatwt[idx, ]
spatwt_df_B_2m$sample_id <- factor(spatwt_df_B_2m$sample_id, levels=samplevels)
spatwt_df_B_2m$Case <- factor(md$Case[match(spatwt_df_B_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_B_2m$pRespByTx <- factor(md$pRespByTx[match(spatwt_df_B_2m$sample_id,md$sample_id)])

funcMarkers<- c("PNAd", "AID", "CD21", "CD23", "KI67", "HLADR")
expr_B<- as.data.frame(df_output[idx, funcMarkers])
expr_B$sample_id<- spatwt_df_B_2m$sample_id
expr_B$Case<- spatwt_df_B_2m$Case
expr_B$pRespByTx <- spatwt_df_B_2m$pRespByTx


do_CN_analysis(spatwt_df = spatwt_df_B_2m,
               expr0 = expr_B,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "pRespByTx", 
               grouplevels = c("Nivo_NR", "Nivo_R", "Nivo_Rela_NR", "Nivo_Rela_R"),
               comparisons = list(c("Nivo_R", "Nivo_Rela_R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_B_15clus_pRespByTx_")

spatwt_df_B_2m$prespgroup <- factor(md$pResponse_Grouped[match(spatwt_df_B_2m$sample_id,md$sample_id)], levels=prespgrouplevels)
expr_B$prespgroup <- spatwt_df_B_2m$prespgroup
do_CN_analysis(spatwt_df = spatwt_df_B_2m,
               expr0 = expr_B,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "prespgroup", 
               grouplevels = c("NR", "R"),
               comparisons = list(c("NR", "R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_B_15clus_pRespGroup_")

## Hep TME analysis====
# Hep
idx <- which(grepl("Hep", df_output$cluster1m))

spatwt_df_Hep_2m<- spatwt[idx, ]
spatwt_df_Hep_2m$sample_id <- factor(spatwt_df_Hep_2m$sample_id, levels=samplevels)
spatwt_df_Hep_2m$Case <- factor(md$Case[match(spatwt_df_Hep_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_Hep_2m$pRespByTx <- factor(md$pRespByTx[match(spatwt_df_Hep_2m$sample_id,md$sample_id)])

funcMarkers<- c("PDL1", "HLADR", "KI67", "CD86", "FGL1")
expr_Hep<- as.data.frame(df_output[idx, funcMarkers])
expr_Hep$sample_id<- spatwt_df_Hep_2m$sample_id
expr_Hep$Case<- spatwt_df_Hep_2m$Case
expr_Hep$pRespByTx <- spatwt_df_Hep_2m$pRespByTx


do_CN_analysis(spatwt_df = spatwt_df_Hep_2m,
               expr0 = expr_Hep,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 20, 
               scaleOption = F,
               groupvar = "pRespByTx", 
               grouplevels = c("Nivo_NR", "Nivo_R", "Nivo_Rela_NR", "Nivo_Rela_R"),
               comparisons = list(c("Nivo_R", "Nivo_Rela_R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_Hep_pRespByTx_")

spatwt_df_Hep_2m$prespgroup <- factor(md$pResponse_Grouped[match(spatwt_df_Hep_2m$sample_id,md$sample_id)], levels=prespgrouplevels)
expr_Hep$prespgroup <- spatwt_df_Hep_2m$prespgroup
do_CN_analysis(spatwt_df = spatwt_df_Hep_2m,
               expr0 = expr_Hep,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 20, 
               scaleOption = F,
               groupvar = "prespgroup", 
               grouplevels = c("NR", "R"),
               comparisons = list(c("NR", "R")),
               group_cols = NULL, 
               filenameprefix = "FunCN_Hep_pRespGroup_")