#==========================================#
#                                          #
#            Abundanace analysis           #
#                                          #
#==========================================#

# load full annotated data
df_output <- readRDS("./df_output.rds")

####Abundance plots======
counts_table<-table(df_output$cluster1m, df_output$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <-as.data.frame.matrix(props_table)

#Densities
densities <- t(t(counts)/areas$TotalArea)

densities_df <- as.data.frame(densities)
densities_df <- data.frame(Cluster = rownames(densities_df), densities_df)
write.xlsx(densities_df, file = "FinalMergedCluster_densities.xlsx", rowNames = FALSE)

cols_to_add <- c("sample_id", "ROI", "Case", "Tx", "Response", "Response_Grouped", 
                 "pResponse", "pResponse_Grouped", "pRespByTx", "CD4_LAG3_RO", "CD8_LAG3_RO")


colorassigned <- c(
  "#9F0712", "#E7180B", "#FB2C36", "#FF6467",
  "#283593", "#1447E6", "#2B7FFF", "#34A6F4", "#42D3F2", "#1C69A8","#5DADE2", "#8EC5FF", "#B8E6FE", "#C4B4FF", "#C27AFF", "#8E51FF", "#8207DB", "#7F22FE",
  "#FF6F00", "#F39C12",
  "#FDD835",
  "#C81CDE",
  "#116045", "#178236", "#1F7A55", "#2AA63E", "#31C950", "#31D492", "#7BF1A8", "#A4F4CF",
  "#C71D36", "#FF2056",
  "#FA058C", "#FC5FB5", "#FDA5D5",
  "#1F2A35", "#2F3E4E", "#445A6A", "#5F7A8A", "#7F9AA8", "#A3B8C2", "#C7D6DD",
  "#000000"
)  # creates 43 distinct colors

# Colors (assigned to clusters)
colorassigned <- setNames(colorassigned, clusterlevels)
#Re-order stackedclusterlevels to put Hep at the bottom
# Convert props to data frame and reshape
props_df <- as.data.frame(props)
props_df$cluster <- rownames(props_df)

props_long <- props_df %>%
  pivot_longer(cols = -cluster, names_to = "sample_id", values_to = "pct") %>%
  dplyr::select(sample_id, everything())  # sample_id first, then everything else

# Set order
# Order samples by increasing Hep proportion
props_long$cluster <- factor(props_long$cluster, levels = clusterlevels)
props_long$response <- factor(md$Response[match(props_long$sample_id,md$sample_id)], levels=resplevels)
props_long$presponse <- factor(md$pResponse_Grouped[match(props_long$sample_id,md$sample_id)], levels=prespgrouplevels)
props_long$case <- factor(md$Case[match(props_long$sample_id,md$sample_id)], levels=caselevels)
props_long$tx <- factor(md$Tx[match(props_long$sample_id,md$sample_id)], levels=txlevels)
props_long$CD4_LAG3_RO <- factor(md$CD4_LAG3_RO[match(props_long$sample_id,md$sample_id)], levels=ROlevels)
props_long$CD8_LAG3_RO <- factor(md$CD8_LAG3_RO[match(props_long$sample_id,md$sample_id)], levels=ROlevels)


#Hep order
props_long_ordered <- props_long %>%
  dplyr::filter(grepl("Hep", cluster)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(total_hep_pct = sum(pct)) %>%
  ungroup()

# Create the ordering vector based on response first, then hep_order
hep_order_vec <- props_long_ordered %>%
  distinct(sample_id, presponse, total_hep_pct, case, tx) %>%
  arrange(presponse, total_hep_pct) %>%  # order by response first, then hep_order
  pull(sample_id)

pdf("FinalMerged_StackedBarProp_SampleID.pdf", width = 21, height = 7.5)

##Figure 3f =====
# Main plot WITHOUT legend
p_main <- ggplot(props_long, aes(x = factor(sample_id, levels = hep_order_vec), y = pct, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.25),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.ticks = element_line(size = 0.25),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(5, 5, 0, 5)
  ) +
  ylab("% of Total Cells") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Cell Composition") +
  scale_fill_manual(values = colorassigned)

# Get unique metadata for annotation bars
anno_data <- props_long %>%
  distinct(sample_id, presponse, case, tx)

# Response annotation bar WITHOUT legend
p_presponse <- ggplot(anno_data, aes(x = factor(sample_id, levels = hep_order_vec), y = 1, fill = presponse)) +
  geom_tile() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  ) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name = "Response")

# Case annotation bar WITHOUT legend
p_case <- ggplot(anno_data, aes(x = factor(sample_id, levels = hep_order_vec), y = 1, fill = case)) +
  geom_tile() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  ) +
  scale_fill_viridis_d(name = "Case", option = "turbo")

# Treatment annotation bar WITHOUT legend
p_tx <- ggplot(anno_data, aes(x = factor(sample_id, levels = hep_order_vec), y = 1, fill = tx)) +
  geom_tile() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 5, 5)
  ) +
  scale_fill_manual(values = c("#66FD8E", "#9052F3"), name = "Treatment")

# Combine plots WITHOUT legends
p_combined <- p_main / p_presponse / p_case / p_tx + 
  plot_layout(heights = c(20, 1, 1, 1))

print(p_combined)

dev.off()


pdf("FinalMerged_StackedBarProp_SampleIDLegend.pdf", width = 4, height = 18)

# Create versions WITH legends to extract them
p_main_legend <- p_main + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "lines")  # smaller key size for 43 items
  ) +
  guides(fill = guide_legend(ncol = 1, title = "Cluster"))

p_presponse_legend <- p_presponse + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  )

p_case_legend <- p_case + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  )

p_tx_legend <- p_tx + 
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# Extract legends
legend_main <- get_legend(p_main_legend)
legend_presponse <- get_legend(p_presponse_legend)
legend_case <- get_legend(p_case_legend)
legend_tx <- get_legend(p_tx_legend)

# Combine all legends with proportional heights
# 43 clusters : 18 cases : 2 treatments ≈ 43:18:2 ratio
plot_grid(
  legend_main,
  legend_presponse,
  legend_case,
  legend_tx,
  ncol = 1,
  align = "v",
  rel_heights = c(43, 2, 18, 2)  # proportional to number of items
)

dev.off()

#Props dataframe
ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf <- left_join(ggdf, md[, cols_to_add], by = "sample_id")

colnames(ggdf) <- gsub("pResponse_Grouped", "prespgroup", colnames(ggdf))
colnames(ggdf) <- gsub("pRespByTx", "presp_by_tx", colnames(ggdf))
colnames(ggdf) <- gsub("Tx", "tx", colnames(ggdf))

ggdf$CD4_LAG3_RO <- factor(ggdf$CD4_LAG3_RO, levels = ROlevels)
ggdf$CD8_LAG3_RO <- factor(ggdf$CD8_LAG3_RO, levels = ROlevels)

#Density dataframe
ggdfd <- reshape2::melt(data.frame(cluster = rownames(densities), densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$cluster <- factor(ggdfd$cluster, levels=clusterlevels)
ggdfd <- left_join(ggdfd, md[, cols_to_add], by = "sample_id")

colnames(ggdfd) <- gsub("pResponse_Grouped", "prespgroup", colnames(ggdfd))
colnames(ggdfd) <- gsub("pRespByTx", "presp_by_tx", colnames(ggdfd))
colnames(ggdfd) <- gsub("Tx", "tx", colnames(ggdfd))

ggdfd$CD4_LAG3_RO <- factor(ggdfd$CD4_LAG3_RO, levels = ROlevels)
ggdfd$CD8_LAG3_RO <- factor(ggdfd$CD8_LAG3_RO, levels = ROlevels)

# ggplot theme 
my_theme <- theme(
  axis.text.x       = element_blank(),
  axis.text.y       = element_text(size = 10, color = "black"),
  axis.title.x      = element_blank(),
  axis.line.x       = element_line(size = 0.25, color = "black"),
  axis.line.y       = element_line(size = 0.25, color = "black"),
  axis.ticks.x      = element_blank(),
  axis.title.y      = element_text(size = 12, color = "black"),
  strip.background  = element_rect(fill = NA),
  strip.text        = element_text(size = 4,  color = "black"),
  panel.background  = element_rect(fill = "white"),
  legend.title      = element_blank(),
  legend.key.size   = unit(2, "lines"),
  legend.text       = element_text(size = 8),
  legend.key        = element_rect(fill = "white"),
  strip.text.x      = element_text(size = 11)
)

#plot box plots - use ggdf for proportion/ggdfd for densities
##Figure 4a (R vs. NR) =====
ggp2<-ggplot(ggdfd,aes(x=prespgroup,y=densities,fill=prespgroup))+
  geom_boxplot(outlier.shape=NA, lwd=0.1)+
  geom_jitter(width=0.1)+
  facet_wrap(~cluster,ncol=5,scales="free")+
  ylab("Density (# of cells)")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 1.5,
                     na.rm = TRUE,
                     ref.group = "NR") +
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
        strip.text.x = element_text(size = 6)
  )
pdf('Density_box_FINALmerged_pRespGrouped_pformat.pdf',width=8,height=10)
ggp2
dev.off()

##Figure 4c (Response stratified by Tx) =====
ggp2<-ggplot(ggdfd,aes(x=presp_by_tx,y=densities,fill=presp_by_tx))+
  geom_boxplot(outlier.shape=NA, lwd=0.1)+
  geom_jitter(width=0.1)+
  facet_wrap(~cluster,ncol=5,scales="free")+
  ylab("Density (# of cells)")+
  ggtitle("J20121 Density Merged")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 1.5,
                     na.rm = TRUE,
                     ref.group = "Nivo_R") +
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
        strip.text.x = element_text(size = 6)
  )
pdf('Density_box_FINALmerged_pRespByTx_pformat.pdf',width=9,height=4)
ggp2
dev.off()

##Figure 4c (Nivo vs. Nivo_Rela) =====
ggp2<-ggplot(ggdfd,aes(x=tx,y=densities,fill=tx))+
  geom_boxplot(outlier.shape=NA, lwd=0.1)+
  geom_jitter(width=0.1)+
  facet_wrap(~cluster,ncol=5,scales="free")+
  scale_fill_manual(values = c("#66FD8E", "#9052F3"))+
  ylab("Density (# of cells)")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 1.5,
                     na.rm = TRUE,
                     ref.group = "Nivo") +
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
        strip.text.x = element_text(size = 6)
  )
pdf('Density_box_FINALmerged_Tx.pdf',width=8,height=6)
ggp2
dev.off()

##Figure 3j =====
# CD4 LAG3 RO + CD8 LAG3 RO BoxPlots - can edit script below to adapt to proportion (ggdf) or density (ggdfd)
ggdf_CD4LAG3RO <- ggdf %>% filter(CD4_LAG3_RO != "NA")
ggdf_CD8LAG3RO <- ggdf %>% filter(CD8_LAG3_RO != "NA")
ggdfd_CD4LAG3RO <- ggdfd %>% filter(CD4_LAG3_RO != "NA")
ggdfd_CD8LAG3RO <- ggdfd %>% filter(CD8_LAG3_RO != "NA")

ggdfd_CD4LAG3RO$CD4_LAG3_RO <- factor(ggdfd_CD4LAG3RO$CD4_LAG3_RO, levels = ROlevels)
ggdfd_CD8LAG3RO$CD8_LAG3_RO <- factor(ggdfd_CD8LAG3RO$CD8_LAG3_RO, levels = ROlevels)

# Density (CD4 LAG3RO)
ggp2<-ggplot(ggdfd_CD4LAG3RO,aes(x=CD4_LAG3_RO,y=densities,fill=CD4_LAG3_RO))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  scale_fill_manual(values = c("Low" = "#009E73", "High" = "#EC407A")) + 
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("Density")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 1.5,
                     na.rm = TRUE,
                     ref.group = "Low") +
  ggtitle("J20121 CD4 LAG3 RO - Merged flowSOM Clusters Densities")+
  my_theme
pdf('0FINAL_AllSubclusters_Density_CD4LAG3RO_pformat.pdf',width=12,height=14)
ggp2
dev.off()

# Density (CD8 LAG3RO)
ggp2<-ggplot(ggdfd_CD8LAG3RO,aes(x=CD8_LAG3_RO,y=densities,fill=CD8_LAG3_RO))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  scale_fill_manual(values = c("Low" = "#009E73", "High" = "#EC407A")) + 
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("Density")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 1.5,
                     na.rm = TRUE,
                     ref.group = "Low") +
  ggtitle("J20121 CD8 LAG3 RO - Merged flowSOM Clusters Densities")+
  my_theme
pdf('0FINAL_AllSubclusters_Density_CD8LAG3RO_pformat.pdf',width=12,height=14)
ggp2
dev.off()


##Figure 3h, 3i =====
# CD4 LAG3 RO
Hep_subtypes <- c("Hep_I", "Hep_II", "Hep_III", "Hep_IV", "Hep_V", "Hep_VI", "Prolif_Hep")

data_hep_subset <- df_output[df_output$cluster1m%in%Hep_subtypes,]
data_hep_subset$CD4_LAG3_RO <- factor(md$CD4_LAG3_RO[match(data_hep_subset$sample_id,md$sample_id)], levels=ROlevels)
data_hep_subset <- data_hep_subset %>% filter(CD4_LAG3_RO != "NA")
data_hep_subset <- data_hep_subset %>%
  mutate(FGL1_group = case_when(
    cluster1m %in% c("Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep") ~ "FGL1low",
    cluster1m %in% c("Hep_II", "Hep_III", "Hep_IV") ~ "FGL1high",
    TRUE ~ NA_character_
  ))

# Get counts and merge Hep clusters as CD4_LAG3_RO high vs low
counts_hep <- table(data_hep_subset$cluster1m, data_hep_subset$CD4_LAG3_RO)
counts_hep <- counts_hep[Hep_subtypes, colnames(counts_hep) != "NA"]

# Get proportion table
props_hep <- t(t(counts_hep) / colSums(counts_hep)) * 100

# Convert to long format for ggplot
props_hep_long <- data.frame(
  cluster = rep(rownames(props_hep), ncol(props_hep)),
  CD4_LAG3_RO = rep(colnames(props_hep), each = nrow(props_hep)),
  proportion = as.vector(props_hep)
)

props_hep_long$cluster <- factor(
  props_hep_long$cluster,
  c("Hep_II", "Hep_III", "Hep_IV", "Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep")
)

props_hep_long$FGL1_group <- ifelse(
  props_hep_long$cluster %in% c("Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep"),
  "FGL1low Hep",
  "FGL1high Hep"
)

# Create color palette:  shades for High,  shades for Low
FGL1_high <- c("Hep_II", "Hep_III", "Hep_IV")
FGL1_low <-  c("Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep")
orange_colors <- c("#FACF8F", "#F8BC63", "#F5A327" )
#hcl.colors(length(high_cases), palette = "Reds", rev = TRUE)[1:length(high_cases)]
purple_colors <- c("#E8BBFC", "#DA8FFA", "#CB63F8", "#B727F5")
#hcl.colors(length(low_cases), palette = "Blues", rev = TRUE)[1:length(low_cases)]
HepColors <- c(setNames(orange_colors, FGL1_high), setNames(purple_colors, FGL1_low))
FGL1Colors <- c(setNames("#F4970B", "FGL1high Hep"), setNames("#BA0BF4", "FGL1low Hep"))

props_hep_long$CD4_LAG3_RO <- factor(props_hep_long$CD4_LAG3_RO, levels = c("Low", "High"))

#Stacked Bar Plot for Hep - with counts + proportion calculation embedded
ggp_bar_prop <- ggplot(props_hep_long, aes(x = CD4_LAG3_RO, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = HepColors) +
  ylab("Proportion (%)") +
  xlab("") +
  ggtitle("All Hep Proportion \nCD4_LAG3_RO High vs. Low") +
  scale_y_continuous(limits = c(0, 100))+
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 18))

ggsave("FinalMerged_HepClusters_stacked_bar_prop_CD4LAG3RO.pdf", ggp_bar_prop, width = 3.25, height = 5)

#Stacked Bar Plot for Hep - with counts + proportion calculation embedded
ggp_bar_prop <- ggplot(props_hep_long, aes(x = CD4_LAG3_RO, y = proportion, fill = FGL1_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  #scale_fill_manual() +
  ylab("Proportion (%)") +
  xlab("") +
  scale_fill_manual(values = FGL1Colors) +
  ggtitle("FGL1 Hep Proportion \nCD4_LAG3_RO High vs. Low") +
  scale_y_continuous(limits = c(0, 100))+
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 18))

ggsave("FinalMerged_HepClustersByFGL1_stacked_bar_prop_CD4LAG3RO.pdf", ggp_bar_prop, width = 3.5, height = 5)


#Generate FGL1high vs. FGL1low Hep clusters to compare CD4 LAG3RO Proportion statistically

# Create grouped categories
FGL1high_Hep <- colSums(counts_hep[c("Hep_II", "Hep_III", "Hep_IV"), ])
FGL1low_Hep <- colSums(counts_hep[c("Hep_I", "Hep_IV", "Hep_V", "Hep_VI", "Prolif_Hep"), ])

# Create contingency table for Fisher's exact test
hep_table <- rbind(FGL1high_Hep, FGL1low_Hep)
colnames(hep_table) <- c("Low CD4 LAG3 RO", "High CD4 LAG3 RO")

print(hep_table)

# Run Fisher's exact test
fisher.test(hep_table)

"Fisher's Exact Test for Count Data

data:  hep_table
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.3080655 0.3281291
sample estimates:
odds ratio 
 0.3179507" 

# CD8 LAG3 RO
# Get data_full_scale and add in CD4LAG3RO metadata 
data_hep_subset <- df_output[df_output$cluster1m%in%Hep_subtypes,]
data_hep_subset$CD8_LAG3_RO <- factor(md$CD8_LAG3_RO[match(data_hep_subset$sample_id,md$sample_id)], levels=ROlevels)
data_hep_subset <- data_hep_subset %>% dplyr::filter(CD8_LAG3_RO != "NA")
data_hep_subset <- data_hep_subset %>%
  mutate(FGL1_group = case_when(
    cluster1m %in% c("Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep") ~ "FGL1low",
    cluster1m %in% c("Hep_II", "Hep_III", "Hep_IV") ~ "FGL1high",
    TRUE ~ NA_character_
  ))


# Get counts and merge Hep clusters as CD4_LAG3_RO high vs low
counts_hep <- table(data_hep_subset$cluster1m, data_hep_subset$CD8_LAG3_RO)
counts_hep <- counts_hep[Hep_subtypes, colnames(counts_hep) != "NA"]

# Get proportion table
props_hep <- t(t(counts_hep) / colSums(counts_hep)) * 100

# Convert to long format for ggplot
props_hep_long <- data.frame(
  cluster = rep(rownames(props_hep), ncol(props_hep)),
  CD8_LAG3_RO = rep(colnames(props_hep), each = nrow(props_hep)),
  proportion = as.vector(props_hep)
)

props_hep_long$cluster <- factor(
  props_hep_long$cluster,
  c("Hep_II", "Hep_III", "Hep_IV", "Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep")
)

props_hep_long$FGL1_group <- ifelse(
  props_hep_long$cluster %in% c("Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep"),
  "FGL1low Hep",
  "FGL1high Hep"
)

# Create color palette:  shades for High,  shades for Low
FGL1_high <- c("Hep_II", "Hep_III", "Hep_IV")
FGL1_low <-  c("Hep_I", "Hep_V", "Hep_VI", "Prolif_Hep")
orange_colors <- c("#FACF8F", "#F8BC63", "#F5A327" )
#hcl.colors(length(high_cases), palette = "Reds", rev = TRUE)[1:length(high_cases)]
purple_colors <- c("#E8BBFC", "#DA8FFA", "#CB63F8", "#B727F5")
#hcl.colors(length(low_cases), palette = "Blues", rev = TRUE)[1:length(low_cases)]
HepColors <- c(setNames(orange_colors, FGL1_high), setNames(purple_colors, FGL1_low))
FGL1Colors <- c(setNames("#F4970B", "FGL1high Hep"), setNames("#BA0BF4", "FGL1low Hep"))

props_hep_long$CD8_LAG3_RO <- factor(props_hep_long$CD8_LAG3_RO, levels = ROlevels)


#Stacked Bar Plot for Hep - with counts + proportion calculation embedded
ggp_bar_prop <- ggplot(props_hep_long, aes(x = CD8_LAG3_RO, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = HepColors) +
  ylab("Proportion (%)") +
  xlab("") +
  ggtitle("All Hep Proportion \nCD8_LAG3_RO High vs. Low") +
  scale_y_continuous(limits = c(0, 100))+
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 18))

ggsave("FinalMerged_HepClusters_stacked_bar_prop_CD8LAG3RO.pdf", ggp_bar_prop, width = 3.25, height = 5)

#Stacked Bar Plot for Hep - with counts + proportion calculation embedded
ggp_bar_prop <- ggplot(props_hep_long, aes(x = CD8_LAG3_RO, y = proportion, fill = FGL1_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  #scale_fill_manual() +
  ylab("Proportion (%)") +
  xlab("") +
  scale_fill_manual(values = FGL1Colors) +
  ggtitle("FGL1 Hep Proportion \nCD8_LAG3_RO High vs. Low") +
  scale_y_continuous(limits = c(0, 100))+
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 18))

ggsave("FinalMerged_HepClustersByFGL1_stacked_bar_prop_CD8LAG3RO.pdf", ggp_bar_prop, width = 3.5, height = 5)

#Generate FGL1high vs. FGL1low Hep clusters to compare CD8 LAG3RO Proportion statistically

FGL1high_Hep <- colSums(counts_hep[c("Hep_II", "Hep_III", "Hep_IV"), ])
FGL1low_Hep <- colSums(counts_hep[c("Hep_I", "Hep_IV", "Hep_V", "Hep_VI", "Prolif_Hep"), ])

# Create contingency table for Fisher's exact test
hep_table <- rbind(FGL1high_Hep, FGL1low_Hep)
colnames(hep_table) <- c("Low CD8 LAG3 RO", "High CD8 LAG3 RO")

print(hep_table)

# Run Fisher's exact test
fisher.test(hep_table)

"Fisher's Exact Test for Count Data

data:  hep_table
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.6391388 0.6769800
sample estimates:
odds ratio 
 0.6577917 " 



##Figure 3i Wilcox test =====
# HEP SUBSET - Hep Subtype Proportion of CD4 and CD8 LAG3 RO
data_hep_subset <- df_output[df_output$cluster1m%in%Hep_subtypes,]
hep_counts_table<-table(data_hep_subset$cluster1m, data_hep_subset$sample_id)
hep_props_table <- t(t(hep_counts_table) / colSums(hep_counts_table)) * 100
hep_counts <- as.data.frame.matrix(hep_counts_table)
hep_props <-as.data.frame.matrix(hep_props_table)

hep_props <- hep_props[Hep_subtypes, ]

hep_ggdf <- melt(data.frame(cluster = rownames(hep_props), hep_props, check.names = FALSE),
                 id.vars = "cluster", value.name = "proportion", 
                 variable.name = "sample_id")
hep_ggdf$sample_id <- factor(hep_ggdf$sample_id, levels=samplevels)
hep_ggdf$cluster <- factor(hep_ggdf$cluster, levels=clusterlevels)
hep_ggdf <- left_join(hep_ggdf, md[, cols_to_add], by = "sample_id")
hep_ggdf$CD4_LAG3_RO <- factor(hep_ggdf$CD4_LAG3_RO, levels = ROlevels)
hep_ggdf$CD8_LAG3_RO <- factor(hep_ggdf$CD8_LAG3_RO, levels = ROlevels)

# CD4 LAG3 RO + CD8 LAG3 RO
hep_ggdf_CD4LAG3RO <- hep_ggdf %>% filter(CD4_LAG3_RO != "NA")
hep_ggdf_CD8LAG3RO <- hep_ggdf %>% filter(CD8_LAG3_RO != "NA")


# CD4 LAG3 RO
ggp2<-ggplot(hep_ggdf_CD4LAG3RO,aes(x=CD4_LAG3_RO,y=proportion,fill=CD4_LAG3_RO))+
  geom_boxplot(outlier.shape=NA, lwd=0.2)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  scale_fill_manual(values = c("Low" = "#009E73", "High" = "#EC407A")) + 
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("Proportion")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 3,
                     na.rm = TRUE,
                     ref.group = "Low") +
  ggtitle("J20121 CD4 LAG3 RO - Hep Clusters")+
  my_theme
pdf('0FINAL_OnlyHepSubclusters_Proportion_CD4LAG3RO_pformat.pdf',width=10,height=4)
ggp2
dev.off()


# CD8 LAG3 RO
ggp2<-ggplot(hep_ggdf_CD8LAG3RO,aes(x=CD8_LAG3_RO,y=proportion,fill=CD8_LAG3_RO))+
  geom_boxplot(outlier.shape=NA, lwd=0.2)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  scale_fill_manual(values = c("Low" = "#009E73", "High" = "#EC407A")) + 
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("Proportion")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = "top",
                     size = 3,
                     na.rm = TRUE,
                     ref.group = "Low") +
  ggtitle("J20121 CD8 LAG3 RO - Hep Clusters")+
  my_theme
pdf('0FINAL_OnlyHepSubclusters_Proportion_CD8LAG3RO_pformat.pdf',width=10,height=5)
ggp2
dev.off()