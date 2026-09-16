#--------------------------------#
#          CN analysis           #
#--------------------------------#
# For CD4 LAG3 RO
source("./functions/clusterfcs_dim.R")

do_CN_analysis <- function(spatwt_df = spatwt_df, #cell of interest filtered, should contain Case, sample_id
                           expr0 = expr0, 
                           funcMarkers = funcMarkers, 
                           clusterby = clusterlevels, #should match the independent variables of spatwt
                           numClusters = 10, 
                           scaleOption = F, 
                           groupvar = "CD4_LAG3_RO", 
                           grouplevels = c("Low", "High"),
                           comparisons = list(c("Low","High")),
                           group_cols = NULL,
                           filenameprefix = "FunCN_Tc_15clus_CD4LAG3RO_")
{
  require(reshape2); require(ggplot2); require(pheatmap);
  require(dplyr); require(Hmisc);require(grid); require(RColorBrewer); 
  require(dittoSeq); require(tidyverse); require(ggpubr);
  
  # used clusterfunc using median metric 
  somClustering <- clusterfcs_dim(fcs=as.matrix(spatwt_df[, clusterby]),
                                  cluster_by = clusterby,
                                  xdims=10, ydims=10, 
                                  numclusters=numClusters, 
                                  scaleoption = scaleOption)
  
  spatwt_df$CN <- somClustering$cell_clustering
  #print(table(spatwt_df$CN, spatwt_df$sample_id))
  
  CNcolors<- dittoColors(reps = 1, get.names = FALSE)[1:numClusters]
  
  ### Dot plot for CN clusters 
  chTME<-c()
  # median of neighborhood weights 
  for (i in 1:numClusters){
    tmp<- apply(spatwt_df[spatwt_df$CN==i,clusterby],2, median)
    chTME<- cbind(chTME, tmp)
    colnames(chTME)[i]<-paste0('CN', i)
  }
  scaled_chTME <- t(apply(chTME, 1, scale))
  colnames(scaled_chTME)<- colnames(chTME)
  
  df1<- reshape2::melt(scaled_chTME)
  names(df1)<-c("celltype",'tme', "zscore")
  df2<- reshape2::melt(chTME)
  names(df2)<-c("celltype",'tme', "median")
  
  df_joined <- left_join(df1, df2, by = c("celltype", "tme"))
  
  df_joined$tme <- factor(df_joined$tme, levels = rev(levels(df_joined$tme)))
  
  # plot only Median weights > 0.01 
  df_filtered <- df_joined %>% 
  dplyr::filter(!is.na(zscore))%>% 
  dplyr::filter(median>0.01) # greater than 1% influence 
  
  ## dot plot for CN clustering ===
  p2<-ggplot(df_filtered, aes(x=celltype, y=tme)) + 
    geom_point(pch = 21, stroke = 0.5, col="black", aes(size = median, fill = zscore)) + 
    theme_classic() +
    labs(size = "Median\nCellular Influences", fill="z-score")+ 
    ylab("") + 
    xlab("")+
    theme(axis.line = element_blank(), 
          axis.text.x= element_text(colour ="black", size = 12),
          axis.text.y= element_text(colour ="black", size = 12),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.7),  
          panel.grid.major = element_line(size = 0.1, linetype = 'solid'),  
          panel.grid.minor = element_line(size = 0.1, linetype = 'solid'), 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7), 
          legend.key.size = unit(0.2, "cm"),
          legend.position="right") + 
    scale_color_discrete("black")+
    scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0)+
    scale_y_discrete(position = "left") + 
    scale_x_discrete(
      guide = guide_axis(angle = 45),position = "top") + 
    scale_size_continuous(breaks = seq(0, max(df_joined$median, na.rm=TRUE), by = 0.1))
  
  legend <- get_legend(p2)
  
  
  pdf(paste0('./output/', filenameprefix, 'CNdotplot.pdf'), width=6, height=4)
  print(p2+theme(legend.position = "none"))
  dev.off()
  
  # CN functional characteristic (CN relative)
  markerExpr <- expr0[, c(funcMarkers, groupvar)]
  markerExpr$CN <- as.factor(paste0("CN", spatwt_df$CN))
  data_for_plot <- reshape2::melt(markerExpr)

  summary_df <- data_for_plot %>%
    dplyr::group_by(variable, CN, .data[[groupvar]])%>% 
    dplyr::summarise(
      median_value = median(value), .groups = "drop")%>%
    dplyr::group_by(variable)%>% 
    dplyr::mutate(scaled_median = as.numeric(scale(median_value)))


  summary_df$CN <- factor(summary_df$CN, levels = paste0("CN", 1:numClusters))

  if (is.null(group_cols)) {
  group_cols <- setNames(scales::hue_pal()(length(grouplevels)), grouplevels)
  }

  # scaled Marker expression (across CN)
  p<- ggplot(summary_df, aes(variable, scaled_median, fill = .data[[groupvar]])) +
    geom_bar(width = 1, stat = "identity", color = "black",  alpha = 0.7)  + 
    scale_y_continuous(breaks = scales::breaks_width(1))+  
    scale_fill_manual(values = group_cols) + 
    geom_hline(yintercept = 0, color = "blue", size = 0.7) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 7),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",  # Explicitly show legend
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_polar() +
    facet_wrap(~CN, ncol=5)

  pdf(paste0("./output/", filenameprefix, "functionalcharac.pdf"), width=8, height = 6)
  print(p)
  dev.off()

  #CN comparison
  markerExpr <- expr0[, c(funcMarkers, groupvar, "sample_id")]
  markerExpr$CN <- as.factor(paste0("CN", spatwt_df$CN))

  CNprop <- as.data.frame(table(markerExpr$CN,markerExpr$sample_id))
  names(CNprop)<-c("CN", "sample_id", "Freq")
  CNprop <- left_join(CNprop, unique(expr0[, c("sample_id", groupvar)]), by="sample_id")
  CNprop$CN <- factor(CNprop$CN, levels = paste0("CN", 1:numClusters))

  # Compare Tumor CN composition between LAG3RO
  CNprop_summarized <- CNprop %>%
    dplyr::filter(!is.na(.data[[groupvar]])) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(
      total_counts = sum(Freq),
      prop = Freq / total_counts * 100
    ) %>%
    dplyr::ungroup()


  CNprop_summarized[[groupvar]] <- factor(CNprop_summarized[[groupvar]], levels = grouplevels)

  pdf(paste0("./output/", filenameprefix, "CNcomparison_violinplots_pformat.pdf"), width=8, height = 6)
  CNcomp<- ggplot(CNprop_summarized, aes(x = .data[[groupvar]], y = prop, fill = .data[[groupvar]])) +
    geom_violin(alpha=0.6)+
    geom_boxplot(width=0.3, fill = NA, color = "black", alpha = 0.5)+
    facet_wrap(~CN, scales = "free_y", ncol=5) + 
    stat_compare_means(
      method = "wilcox.test",
      comparisons = comparisons,
      hide.ns = FALSE,
      label="p.format")+ 
    scale_fill_manual(values = group_cols) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+ 
    theme_bw()+ 
    xlab("") + ylab("Proportion of CN")+
    theme(strip.text = element_text(size=12, color="black"), 
          axis.text.x = element_text(angle=45,  hjust =1), 
          axis.text = element_text(size=12, color="black"))
  #legend.position = "none")
  print(CNcomp)
  dev.off()
  
  return(somClustering$cell_clustering)
}