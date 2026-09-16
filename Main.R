#Set directory ======
work<-getwd()

# create output folder if running for the first time 
output_folder <- file.path(work, "output")

# Check if the folder exists
if (!dir.exists(output_folder)) {
  # Create the folder if it does not exist
  dir.create(output_folder)
  cat("output folder is created.\n")
} else {
  cat("output folder already exists.\n")
}

#Load libraries ======

## Check and install missing packages
required_packages <- unique(c(
  "readxl", "harmony", "pals", "flowCore", "ggridges", "plyr","dplyr", "ggplot2", 
  "tibble", "reshape2", "cowplot", "RColorBrewer", "FNN", "stringr", "fields", 
  "tidyr", "ComplexHeatmap", "pheatmap", "circlize", "ggpubr", "tidyverse", 
  "ggsignif", "broom", "TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
  "survival", "survminer","ggh4x", "BiocManager","spatstat","reshape","Hmisc","grid",
  "dittoSeq","ConsensusClusterPlus","tiff","raster", "scales", "DataEditR", 
  "FlowSOM","terra", "spatstat", "matrixStats", "Spectre", "openxlsx", "patchwork"
))


# Install missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load all required packages
lapply(required_packages, require, character.only = TRUE)

#Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load functions ====
source('./functions/clusterfcs_dim.R')
source('./functions/do_CN_analysis.R')
source('./functions/compute_FunCN.R')


#Set up levels ======
samplevels=c("J20121_P2_ROI1","J20121_P2_ROI2","J20121_P2_ROI3","J20121_P3_ROI1","J20121_P3_ROI2","J20121_P3_ROI3","J20121_P3_ROI4",
             "J20121_P4_ROI1","J20121_P4_ROI2","J20121_P4_ROI3","J20121_P4_ROI4","J20121_P6_ROI1","J20121_P6_ROI2","J20121_P6_ROI3","J20121_P6_ROI4",
             "J20121_P7_ROI1","J20121_P7_ROI2","J20121_P14_ROI1","J20121_P14_ROI2","J20121_P14_ROI3","J20121_P14_ROI4","J20121_P15_ROI1","J20121_P15_ROI2","J20121_P15_ROI3",
             "J20121_P17_ROI1","J20121_P17_ROI2","J20121_P17_ROI3","J20121_P17_ROI4","J20121_P21_ROI1","J20121_P21_ROI2","J20121_P21_ROI3",
             "J20121_P22_ROI1","J20121_P22_ROI2","J20121_P22_ROI3","J20121_P22_ROI4","J20121_P23_ROI1","J20121_P23_ROI2","J20121_P23_ROI3","J20121_P23_ROI4",
             "J20121_P24_ROI1","J20121_P24_ROI2","J20121_P24_ROI3","J20121_P24_ROI4","J20121_P24_ROI5","J20121_P25_ROI1","J20121_P25_ROI2","J20121_P25_ROI3","J20121_P25_ROI4",
             "J20121_P27_ROI1","J20121_P27_ROI2","J20121_P27_ROI3","J20121_P27_ROI4","J20121_P28_ROI1","J20121_P28_ROI2","J20121_P28_ROI3","J20121_P28_ROI4",
             "J20121_P29_ROI1","J20121_P29_ROI2","J20121_P29_ROI3","J20121_P29_ROI4","J20121_P29_ROI5","J20121_P30_ROI1","J20121_P30_ROI2","J20121_P30_ROI3",
             "J20121_P31_ROI1","J20121_P31_ROI2","J20121_P31_ROI3")

caselevels=c("J20121_P2","J20121_P3","J20121_P4","J20121_P6","J20121_P7",
             "J20121_P14","J20121_P15","J20121_P17","J20121_P21","J20121_P22",
             "J20121_P23","J20121_P24","J20121_P25","J20121_P27","J20121_P28",
             "J20121_P29","J20121_P30","J20121_P31")

ROIlevels=c("ROI1", "ROI2", "ROI3", "ROI4", "ROI5")

txlevels=c("Nivo", "Nivo_Rela")

resplevels=c("PD", "SD", "PR")

respgrouplevels=c("NR", "R")

presplevels=c("No", "Minor", "Major", "Complete")

prespgrouplevels=c("NR", "R")

ROlevels=c("Low", "High", "NA")

prespbytxlevels=c("Nivo_NR", "Nivo_R", "Nivo_Rela_NR", "Nivo_Rela_R")

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


## Metadata =======
areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
metaDataFile = paste0(work,"/Config/metadata.xlsx")

ifelse(grepl(metaDataFile,pattern='.xlsx'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
md$file_name <- factor(md$file_name)
md$sample_id <- factor(md$sample_id, levels=samplevels)
md$Case <- factor(md$Case)
md$ROI <- factor(md$ROI)
md$Tx <- factor(md$Tx)
md$Response <- factor(md$Response)
md$Response_Grouped <- factor(md$Response_Grouped)
md$pResponse <- factor(md$pResponse)
md$pResponse_Grouped <- factor(md$pResponse_Grouped)
md$pRespByTx <- factor(md$pRespByTx)
md$CD4_LAG3_RO <- factor(md$CD4_LAG3_RO)
md$CD8_LAG3_RO <- factor(md$CD8_LAG3_RO)


fullmarkerlist <- c("CD45", "CD20", "CD21", "CD23","AID", "CD138", 
                    "CD3", "CD4", "CD8", "CD45RA", "CD45RO", "BCL6", "CD25", "GZMB", "PD1", "TOX", "LAG3", "FOXP3",
                    "CD57", "DCLAMP", "CD15",
                    "CD68", "HLADR", "CD86", "CD163", "CD206", "PDL1",
                    "PNAd", "CD31",
                    "CK", "HEPPAR", "FGL1",
                    "SMA", "PDPN", "VIM", "KI67")
                     