# Load required libraries
library(glmmTMB)
library(car)      # For Type II Wald chi-square omnibus tests
library(emmeans)  # For Estimated Marginal Means and pairwise contrasts

# 1. Define the Master Dictionary of Families
final_families <- list(
  "Ratio_VD_S"                    = ordbeta(link = "logit"),
  "Ratio_Nom"                    = Gamma(link = "log"),
  "Ratio_LW_S"                    = tweedie(link = "log"),
  "av_nominal_deps_NN"            = gaussian(link = "identity"),
  "cl_av_deps"                    = gaussian(link = "identity"),
  "nn_all_nominal_deps_NN_struct" = tweedie(link = "log"),
  "abundance"                     = nbinom2(link = "log"),
  "mattr_50"                      = beta_family(link = "logit"),
  "evenness_50"                   = beta_family(link = "logit"),
  "disparity"                     = tweedie(link = "log"), 
  "dispersion_per_100"            = gaussian(link = "identity"),
  "QWE_Score"                     = beta_family(link = "logit"),
  "BERTScore_F1"                  = gaussian(link = "identity")
)

# 2. Initialize a master list to safely store all outputs
final_analysis_results <- list()

# 3. Execute the Master Loop
for (index in names(final_families)) {
  
  message("\n==================================================")
  message("Processing Index: ", index)
  
  fam <- final_families[[index]]
  
  # Handle the disparity shift dynamically
  if (index == "disparity") {
    response_var <- "disparity_adj"
    message("Note: Using shifted variable 'disparity_adj'")
  } else {
    response_var <- index
  }
  
  # Construct the uniform random intercept formula
  form <- as.formula(paste(response_var, "~ LLM + (1 | ID)"))
  
  # Wrap in tryCatch to ensure one problematic index doesn't crash the whole loop
  tryCatch({
    
    # --- A. Fit the Model ---
    mod <- glmmTMB(form, data = df_clean, family = fam)
    
    # --- B. Omnibus Test (Type II Wald Chi-Square) ---
    # Tests the main question: Does 'LLM' have a significant overall effect?
    omnibus_test <- car::Anova(mod, type = "II")
    p_value <- omnibus_test$`Pr(>Chisq)`[1]
    
    message("Omnibus Test p-value for LLM: ", format.pval(p_value, eps = 0.001))
    
    # --- C. Estimated Marginal Means & Pairwise Contrasts ---
    # Computes the means for each LLM and tests every pairwise combination (Tukey adjusted)
    # Note: type = "response" ensures EMMs are back-transformed to the original data scale
    emms <- emmeans(mod, pairwise ~ LLM, adjust = "tukey", type = "response")
    
    # --- D. Store Results ---
    final_analysis_results[[index]] <- list(
      model        = mod,
      omnibus      = omnibus_test,
      emmeans_obj  = emms$emmeans,
      contrasts    = emms$contrasts
    )
    
    message("Successfully fitted model and extracted EMMs.")
    
  }, error = function(e) {
    message(" -> [Error] processing ", index, ": ", e$message)
  })
}

message("\nAnalysis complete! All models, omnibus tests, and contrasts are safely stored in 'final_analysis_results'.")

# Load required libraries
library(ggplot2)
library(dplyr)
library(purrr)

# ---------------------------------------------------------
# 1. Extract and Bind Data from All Indices
# ---------------------------------------------------------
# Initialize an empty list to collect the dataframes
master_emms_list <- list()

# Loop through the final_analysis_results list
for (index in names(final_analysis_results)) {
  
  # Check if the emmeans object exists for this index to avoid errors
  if (!is.null(final_analysis_results[[index]]$emmeans_obj)) {
    
    # Extract and convert to dataframe
    df_temp <- as.data.frame(final_analysis_results[[index]]$emmeans_obj)
    
    # Standardize the Mean column name (catches 'response' and changes it to 'emmean')
    if ("response" %in% colnames(df_temp)) {
      colnames(df_temp)[colnames(df_temp) == "response"] <- "emmean"
    }
    
    # Standardize Confidence Interval column names dynamically
    colnames(df_temp)[grepl("LCL|lower.CL", colnames(df_temp))] <- "LCL"
    colnames(df_temp)[grepl("UCL|upper.CL", colnames(df_temp))] <- "UCL"
    
    # Add a column specifying which index this data belongs to
    df_temp$Index <- index
    
    # Store in our collection list
    master_emms_list[[index]] <- df_temp
  }
}

# Combine all individual dataframes into one massive dataframe
df_master <- bind_rows(master_emms_list)

# ---------------------------------------------------------
# 2. Data Preparation: Alphabetical Y-Axis
# ---------------------------------------------------------
# Force the LLM column to be a factor ordered alphabetically
df_master$LLM <- factor(df_master$LLM, levels = sort(unique(as.character(df_master$LLM))))

# (Optional) Clean up the Index names so they look better in the plot headers
# df_master$Index <- gsub("_", " ", df_master$Index)

# ---------------------------------------------------------
# 3. Generate the Master Faceted Forest Plot
# ---------------------------------------------------------
p_master <- ggplot(df_master, aes(x = emmean, y = LLM, color = LLM)) +
  
  # Add error bars first
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), 
                 height = 0.3, 
                 linewidth = 0.7, 
                 show.legend = FALSE) + # Hide legend for error bars
  
  # Add the means
  geom_point(size = 2.5, show.legend = FALSE) + 
  
  # Create a separate panel for each index. 
  # 'scales = "free_x"' is CRITICAL here so each metric gets its own proper x-axis scale
  facet_wrap(~ Index, scales = "free_x", ncol = 3) + 
  
  # Color palette (Optional: makes distinguishing the LLMs easier across panels)
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  
  # Labels and Title
  labs(title = "Estimated Marginal Means Across All Linguistic Indices",
       subtitle = "Error bars represent 95% Confidence Intervals",
       x = "Estimated Value (Scale varies by index)",
       y = "Large Language Model") +
  
  # Clean Academic Theme
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    strip.background = element_rect(fill = "#f0f0f0", color = "grey50"),
    strip.text = element_text(face = "bold", size = 10), # Makes the Index names pop
    axis.text.y = element_text(face = "bold", size = 9, color = "black"),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines") # Adds breathing room between the panels
  )

# Display the master plot
print(p_master)

# Save the master plot (Recommend a larger canvas to fit all 13 panels comfortably)
ggsave("TFM-Master_Forest_Plot_Global.png", plot = p_master, width = 12, height = 10, dpi = 300)

for (index in names(final_families)) {
  message("### ---------  ###")
  message(" -> Processing ", index)
  print(final_analysis_results[[index]]$omnibus)
  print(final_analysis_results[[index]]$contrasts)
  print(final_analysis_results[[index]]$emmeans_obj)
}
