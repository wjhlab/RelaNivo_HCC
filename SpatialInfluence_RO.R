df_output <- readRDS("./df_output.rds")
spatwt <- readRDS("./spatwt.RDS")


##CD4_LAG3_RO =====
df_output$CD4_LAG3_RO <- factor(md$CD4_LAG3_RO[match(df_output$sample_id,md$sample_id)])

# Tc
idx <- which(grepl("Tc", df_output$cluster1m) & df_output$CD4_LAG3_RO != "NA")

spatwt_df_cd8T_2m<- spatwt[idx, ]
spatwt_df_cd8T_2m$sample_id <- factor(spatwt_df_cd8T_2m$sample_id, levels=samplevels)
spatwt_df_cd8T_2m$Case <- factor(md$Case[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_cd8T_2m$CD4_LAG3_RO <- factor(md$CD4_LAG3_RO[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)])


funcMarkers<- c("BCL6", "CD25", "TOX", "PD1", "GZMB", "KI67", "LAG3", "HLADR")

expr_cd8T<- as.data.frame(df_output[idx, funcMarkers])
expr_cd8T$sample_id<- spatwt_df_cd8T_2m$sample_id
expr_cd8T$Case<- spatwt_df_cd8T_2m$Case
expr_cd8T$CD4_LAG3_RO <- spatwt_df_cd8T_2m$CD4_LAG3_RO


do_CN_analysis(spatwt_df = spatwt_df_cd8T_2m,
               expr0 = expr_cd8T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "CD4_LAG3_RO", 
               grouplevels = c("Low", "High"),
               comparisons = list(c("Low","High")),
               group_cols = c("Low" = "#009E73", "High" = "#EC407A"), 
               filenameprefix = "FunCN_Tc_15clus_CD4LAG3RO_")


#Th
idx <- which(grepl("Th", df_output$cluster1m) & df_output$CD4_LAG3_RO != "NA")

spatwt_df_cd4T_2m<- spatwt[idx, ]
spatwt_df_cd4T_2m$sample_id <- factor(spatwt_df_cd4T_2m$sample_id, levels=samplevels)
spatwt_df_cd4T_2m$Case <- factor(md$Case[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_cd4T_2m$CD4_LAG3_RO <- factor(md$CD4_LAG3_RO[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)])

funcMarkers<- c("BCL6", "CD25", "TOX", "PD1", "FOXP3", "KI67", "LAG3", "HLADR")
expr_cd4T<- as.data.frame(df_output[idx, funcMarkers])
expr_cd4T$sample_id<- spatwt_df_cd4T_2m$sample_id
expr_cd4T$Case<- spatwt_df_cd4T_2m$Case
expr_cd4T$CD4_LAG3_RO <- spatwt_df_cd4T_2m$CD4_LAG3_RO


do_CN_analysis(spatwt_df = spatwt_df_cd4T_2m,
               expr0 = expr_cd4T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "CD4_LAG3_RO", 
               grouplevels = c("Low", "High"),
               comparisons = list(c("Low","High")),
               group_cols = c("Low" = "#009E73", "High" = "#EC407A"), 
               filenameprefix = "FunCN_Th_15clus_CD4LAG3RO_")


##CD8_LAG3_RO =====
df_output$CD8_LAG3_RO <- factor(md$CD8_LAG3_RO[match(df_output$sample_id,md$sample_id)])

# Tc
idx <- grep("Tc", df_output$cluster1m)
idx <- idx[df_output$CD8_LAG3_RO[idx] != "NA"]

spatwt_df_cd8T_2m<- spatwt[idx, ]
spatwt_df_cd8T_2m$sample_id <- factor(spatwt_df_cd8T_2m$sample_id, levels=samplevels)
spatwt_df_cd8T_2m$Case <- factor(md$Case[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_cd8T_2m$CD8_LAG3_RO <- factor(md$CD8_LAG3_RO[match(spatwt_df_cd8T_2m$sample_id,md$sample_id)])

funcMarkers<- c("BCL6", "CD25", "TOX", "PD1", "GZMB", "KI67", "LAG3", "HLADR")
expr_cd8T<- as.data.frame(df_output[idx, funcMarkers])
expr_cd8T$sample_id<- spatwt_df_cd8T_2m$sample_id
expr_cd8T$Case<- spatwt_df_cd8T_2m$Case
expr_cd8T$CD8_LAG3_RO <- spatwt_df_cd8T_2m$CD8_LAG3_RO


res <- do_CN_analysis(spatwt_df = spatwt_df_cd8T_2m,
                      expr0 = expr_cd8T,
                      funcMarkers = funcMarkers,
                      clusterby = clusterlevels,
                      numClusters = 15, 
                      scaleOption = F,
                      groupvar = "CD8_LAG3_RO", 
                      grouplevels = c("Low", "High"),
                      comparisons = list(c("Low","High")),
                      group_cols = c("Low" = "#009E73", "High" = "#EC407A"), 
                      filenameprefix = "FunCN_Tc_15clus_CD8LAG3RO_")

#Ext. Figure 4e
markerExpr_CD8T <- expr_cd8T[, c(funcMarkers, "CD8_LAG3_RO", "sample_id")]
markerExpr_CD8T$CN <- as.factor(paste0("CN",res))

# CN2 Tc 
CN2_Tc <- markerExpr_CD8T[markerExpr_CD8T$CN == "CN2", ] 
CN2_Tc$CD8_LAG3_RO <- factor(CN2_Tc$CD8_LAG3_RO, levels=c("Low", "High"))

# GZMB expr CD8 LAG3 RO low vs. high
figS4e<- ggplot(CN2_Tc, aes(x = CD8_LAG3_RO, y = GZMB, fill = CD8_LAG3_RO)) +
  geom_violin(alpha=0.6)+
  geom_boxplot(width=0.3, fill = NA, color = "black", alpha = 0.5)+
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Low", "High")),
    hide.ns = FALSE,
    label="p.format")+ 
  scale_fill_manual(values = c("Low" = "#009E73", "High" = "#EC407A")) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+ 
  theme_bw()+ 
  xlab("") + ylab("GZMB Expression")+ ggtitle("CN2 Tc")+
  theme(strip.text = element_text(size=12, color="black"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=12, color="black"))
pdf("./output/FigS4e_GZMB_expr_CD8LAG3RO.pdf", width=3, height=3)
print(figS4e)
dev.off()


#Th
idx <- grep("Th", df_output$cluster1m)
idx <- idx[df_output$CD8_LAG3_RO[idx] != "NA"]

spatwt_df_cd4T_2m<- spatwt[idx, ]
spatwt_df_cd4T_2m$sample_id <- factor(spatwt_df_cd4T_2m$sample_id, levels=samplevels)
spatwt_df_cd4T_2m$Case <- factor(md$Case[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)], levels=caselevels)
spatwt_df_cd4T_2m$CD8_LAG3_RO <- factor(md$CD8_LAG3_RO[match(spatwt_df_cd4T_2m$sample_id,md$sample_id)])

funcMarkers<- c("BCL6", "CD25", "TOX", "PD1", "FOXP3", "KI67", "LAG3", "HLADR")
expr_cd4T<- as.data.frame(df_output[idx, funcMarkers])
expr_cd4T$sample_id<- spatwt_df_cd4T_2m$sample_id
expr_cd4T$Case<- spatwt_df_cd4T_2m$Case
expr_cd4T$CD8_LAG3_RO <- spatwt_df_cd4T_2m$CD8_LAG3_RO


do_CN_analysis(spatwt_df = spatwt_df_cd4T_2m,
               expr0 = expr_cd4T,
               funcMarkers = funcMarkers,
               clusterby = clusterlevels,
               numClusters = 15, 
               scaleOption = F,
               groupvar = "CD8_LAG3_RO", 
               grouplevels = c("Low", "High"),
               comparisons = list(c("Low","High")),
               group_cols = c("Low" = "#009E73", "High" = "#EC407A"), 
               filenameprefix = "FunCN_Th_15clus_CD8LAG3RO_")

