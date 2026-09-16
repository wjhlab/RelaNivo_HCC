# RelaNivo_HCC

All necessary scripts for analyzing the IMC dataset are in this repository. The "Config" folder contains relevant metadata and clustering annotation text. 
All associated RDS files are in Zenodo at: 

## 1. Cell-type clustering and annotation
FlowSOMClusteringPipeline.R takes all single-cell csv files, performs patient-level normalization, and generates a cluster heatmap for annotation of broad cell types.

Scripts within the "Subclustering" folder are then used to further resolve each broadly annotated cell type into more detailed subclusters and cell phenotypes.

Subclustering should be performed sequentially in the following order:
Lymphoid --> Myeloid --> Hep --> Stroma --> Endoth

Because each subclustering step generates an updated dataframe, subclustering must be run sequentially in this order to obtain the final 'df_output' dataset containing all detailed cell annotations.


## 2. Spatial analysis setup
Once a fully annotated dataset 'df_output.RDS' has been generated, spatial analyses can be performed using Main.R

Main.R reformats the annotated dataset for spatial analysis. Compute_FunCN is applied to cell coordinates and cell annotations to generate the 'spatwt.RDS' spatial-weight object used in downstream spatial analyses. 

The "functions" folder includes necessary code to execute spatial analyses and functions: compute_FunCN.R, clusterfcs_dim.R, and do_CN_analysis


## 3. Spatial influence analysis
After generating 'spatwt.RDS' the spatial influence analysis scripts can be run.

SpatialInfluenceDataframe.R contains the code for running the main spatial influence analysis pipeline and code for generating the associated plots. 

SpatialInfluence_pResp.R performs spatial influence clustering and downstream analyses stratified by response.

SpatialInfluence_RO.R performs spatial influence clustering and downstream analyses stratified by LAG3 receptor occupancy.
