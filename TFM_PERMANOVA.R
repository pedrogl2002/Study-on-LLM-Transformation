# Load required package for PERMANOVA
# install.packages("vegan")
library(vegan)

# 1. Define the final selected variables (retaining 'abundance', dropping 'volume')
selected_indices <- c("Ratio_VD_S", "Ratio_Nom", "Ratio_LW_S", "av_nominal_deps_NN", 
                      "cl_av_deps", "nn_all_nominal_deps_NN_struct", "abundance", 
                      "mattr_50", "evenness_50", "disparity", 
                      "dispersion_per_100", "QWE_Score", "BERTScore_F1")
lyrics_data <- df_clean
# 2. Subset the data
# Ensure 'ID' is a factor for the permutations to work correctly
lyrics_data <- lyrics_data %>%
  mutate(ID = as.factor(ID),
         LLM = as.factor(LLM))

# Extract only the selected variables
perm_data <- lyrics_data %>% select(all_of(selected_indices))

# 3. Scale the data
# This step is critical so that variables with larger ranges (like abundance) 
# do not completely dominate the Euclidean distance calculation.
perm_data_scaled <- scale(perm_data)

# 4. Run PERMANOVA
# We use Euclidean distance on the scaled data.
# 'strata = lyrics_data$ID' restricts permutations to within each song block.
set.seed(42) # For reproducibility of permutations
permanova_result <- adonis2(perm_data_scaled ~ LLM, 
                            data = lyrics_data, 
                            method = "euclidean", 
                            strata = lyrics_data$ID,
                            permutations = 999)

# Print the results
print(permanova_result)

# Calculate multivariate dispersions
dispersion_test <- betadisper(vegdist(perm_data_scaled, method = "euclidean"), df_clean$LLM)

# Test for significance (ANOVA on the distances to group centroids)
anova(dispersion_test)
