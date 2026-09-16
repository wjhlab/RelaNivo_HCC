#==========================================#
#                                          #
#        IMC J20121 - Spatial Influence    #
#                                          #
#==========================================#

# Having df_output.rds, start from here ====================================

## Reformat df_output to run spatial code ====
# Reformat data_full_scaled to have detailed clusters as cluster, large cluster categories as cluster1m, take out expr matrix
# Only run for the first time, afterwards simply load saved RDS file "df_output.RDS"
df_output <- data.frame(readRDS(paste0(work, "/df_output.rds")))

# check factor levels
str(df_output)

# set levels for celltypes
df_output$cluster1m <- factor(df_output$cluster1m, levels=clusterlevels)

# Spatwt
## Apply FunCN & compare TME ====
# Dot plot for pairwise Cell-cell interaction in pResponse_Grouped
spatwt <- compute_FunCN(expr = df_output,  # should contain sample_id, X, Y position, celltype in the df
                        x_col = "X_coord",
                        y_col = "Y_coord",
                        celltype_col = "cluster1m",
                        multisample = TRUE,
                        sample_col = "sample_id", 
                        kernels = "gaussian",
                        sigma = 10) 
spatwt$sample_id<- factor(spatwt$sample_id, levels=samplevels)
# save spatwt 
saveRDS(spatwt, "spatwt.rds")
