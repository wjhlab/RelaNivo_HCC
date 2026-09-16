#Merge lymphoid subclusters
# Load all files
global_data <- readRDS("data_full.RDS")
parent_df <- data.frame(global_data)

# Load and combine subclusters

data_lym <- data.frame(readRDS("Lymphoid_Subcluster.rds")[[1]])

# Create mapping and update
mapping_df <- data_lym %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  select(cell_key, cluster1m_new = lymphoid_cluster1m)

parent_df_updated <- parent_df %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  left_join(mapping_df, by = "cell_key") %>%
  mutate(cluster1m = ifelse(!is.na(cluster1m_new), cluster1m_new, cluster1m)) %>%
  select(-cell_key, -cluster1m_new)

# Save
saveRDS(parent_df_updated, "parent_dataframe_updatedwithlymphoidsubcluster.rds")
data_full_scaled <- parent_df_updated




#Merge myeloid subclusters
# Load all files
global_data <- readRDS("parent_dataframe_updatedwithlymphoidsubcluster.rds")
parent_df <- data.frame(global_data)

# Load and combine subclusters

data_mye <- data.frame(readRDS("Myeloid_Subcluster.rds")[[1]])

# Create mapping and update
mapping_df <- data_mye %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  select(cell_key, cluster1m_new = myeloid_cluster1m)

parent_df_updated <- parent_df %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  left_join(mapping_df, by = "cell_key") %>%
  mutate(cluster1m = ifelse(!is.na(cluster1m_new), cluster1m_new, cluster1m)) %>%
  select(-cell_key, -cluster1m_new)

# Save
saveRDS(parent_df_updated, "parent_dataframe_updatedwithlymphoid_myeloidsubcluster.rds")
data_full_scaled <- parent_df_updated




#Merge Hep subclusters
# Load all files
global_data <- readRDS("parent_dataframe_updatedwithlymphoid_myeloidsubcluster.rds")
parent_df <- data.frame(global_data)

# Load and combine subclusters

data_hep <- data.frame(readRDS("Hep_Subcluster.rds")[[1]])

# Create mapping and update
mapping_df <- data_hep %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  select(cell_key, cluster1m_new = hep_cluster1m)

parent_df_updated <- parent_df %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  left_join(mapping_df, by = "cell_key") %>%
  mutate(cluster1m = ifelse(!is.na(cluster1m_new), cluster1m_new, cluster1m)) %>%
  select(-cell_key, -cluster1m_new)

# Save
saveRDS(parent_df_updated, "parent_dataframe_updatedwithlymphoid_myeloid_hepsubcluster.rds")
data_full_scaled <- parent_df_updated



#Merge Stroma subclusters
# Load all files
global_data <- readRDS("parent_dataframe_updatedwithlymphoid_myeloid_hepsubcluster.rds")
parent_df <- data.frame(global_data)

# Load and combine subclusters

data_str <- data.frame(readRDS("Stroma_Subcluster.rds")[[1]])

# Create mapping and update
mapping_df <- data_str %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  select(cell_key, cluster1m_new = str_cluster1m)

parent_df_updated <- parent_df %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  left_join(mapping_df, by = "cell_key") %>%
  mutate(cluster1m = ifelse(!is.na(cluster1m_new), cluster1m_new, cluster1m)) %>%
  select(-cell_key, -cluster1m_new)

# Save
saveRDS(parent_df_updated, "parent_dataframe_updatedwithlymphoid_myeloid_hep_stromasubcluster.rds")
data_full_scaled <- parent_df_updated






#Merge Endoth subclusters
# Load all files
global_data <- readRDS("parent_dataframe_updatedwithlymphoid_myeloid_hep_stromasubcluster.rds")
parent_df <- data.frame(global_data)

# Load and combine subclusters

data_endoth <- data.frame(readRDS("Endoth_Subcluster.rds")[[1]])

# Create mapping and update
mapping_df <- data_endoth %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  select(cell_key, cluster1m_new = endoth_cluster1m)

parent_df_updated <- parent_df %>%
  mutate(cell_key = paste(sample_id, CellId, sep = "_")) %>%
  left_join(mapping_df, by = "cell_key") %>%
  mutate(cluster1m = ifelse(!is.na(cluster1m_new), cluster1m_new, cluster1m)) %>%
  select(-cell_key, -cluster1m_new)

# Save
saveRDS(parent_df_updated, "parent_dataframe_updatedwithlymphoid_myeloid_hep_stroma_endothsubcluster.rds")
data_full_scaled <- parent_df_updated

