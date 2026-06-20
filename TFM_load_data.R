#### Carga de datos

# 1. Get a list of all CSV files in the folder   
library(tidyverse)
files <- list.files(pattern = "^merged.*\\.csv", full.names = TRUE)

# 2. Read them all and bind them into one data frame
df <- files %>%
  map_df(~read_csv(.))
# Assume 'df' is your long-format data: Lyric, LLM, Index1...Index
index_cols <- c('Ratio_VD_S', 'Ratio_Nom', 'Ratio_LW_S', 
                'av_nominal_deps_NN', 'cl_av_deps', 'nn_all_nominal_deps_NN_struct', 
                'volume', 'abundance', 'mattr_50', 'evenness_50', 'disparity', 
                'dispersion_per_100', 'QWE_Score', 'BERTScore_F1')
df <- df %>% mutate(ID = factor(ID), LLM = factor(LLM))

df_clean <- as_tibble(df) %>% dplyr::select(-volume)

df_clean$disparity_adj <- df_clean$disparity - 1
df$disparity_adj <- df$disparity - 1
