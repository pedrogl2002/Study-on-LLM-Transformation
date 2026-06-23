# Load required libraries
library(glmmTMB)
library(car)         
library(emmeans)     
library(performance) 
library(multcomp)    
library(dplyr)       

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

# 2. Initialize storage
final_analysis_results <- list()
summary_table <- data.frame(
  Index = character(),
  Wald_Test = character(),
  R2_Marginal = character(),   # Changed to character to allow "NA" text
  R2_Conditional = character(),# Changed to character to allow "NA" text
  EMM_Range = character(),
  Hierarchy_CLD = character(),
  stringsAsFactors = FALSE
)

# 3. Execute the Master Loop
for (index in names(final_families)) {
  
  message("\nProcessing Index: ", index, "...")
  fam <- final_families[[index]]
  
  response_var <- ifelse(index == "disparity", "disparity_adj", index)
  form <- as.formula(paste(response_var, "~ LLM + (1 | ID)"))
  #form <- as.formula(paste(response_var, "~ LLM + (1 + LLM | ID)"))  #maximal model
  tryCatch({
    # --- A. Fit Model ---
    mod <- glmmTMB(form, data = df_clean, family = fam)
    
    # --- B. Omnibus Test ---
    omnibus <- car::Anova(mod, type = "II")
    chisq_val <- round(omnibus$`Chisq`[1], 2)
    df_val <- omnibus$`Df`[1]
    p_val <- omnibus$`Pr(>Chisq)`[1]
    
    p_formatted <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
    wald_string <- sprintf("χ²(%d) = %.2f, %s", df_val, chisq_val, p_formatted)
    
    # --- C. Safely Extract R-Squared ---
    # If this fails for complex families, it will just assign "NA" and move on
    r2_m <- "NA"
    r2_c <- "NA"
    tryCatch({
      suppressWarnings({
        r2_res <- r2_nakagawa(mod)
        if (!is.null(r2_res)) {
          r2_m <- as.character(round(r2_res$R2_marginal, 3))
          r2_c <- as.character(round(r2_res$R2_conditional, 3))
        }
      })
    }, error = function(e) {
      message("  -> [Note]: R-squared unsupported for this family. Returning NA.")
    })
    
    # --- D. Estimated Marginal Means & Range ---
    emms <- emmeans(mod, ~ LLM, type = "response")
    emms_df <- as.data.frame(emms)
    
    # Bulletproof column renaming: The mean estimate is ALWAYS the 2nd column
    colnames(emms_df)[2] <- "emmean"
    
    min_row <- emms_df[which.min(emms_df$emmean), ]
    max_row <- emms_df[which.max(emms_df$emmean), ]
    emm_range_str <- sprintf("%.3f (%s) - %.3f (%s)", 
                             min_row$emmean, min_row$LLM, 
                             max_row$emmean, max_row$LLM)
    
    # --- E. Safely Generate Performance Hierarchy (CLD) ---
    cld_res <- cld(emms, Letters = letters, adjust = "sidak")
    
    # NEW FIX: Standardize the mean column name inside cld_res before sorting
    if ("response" %in% colnames(cld_res)) {
      colnames(cld_res)[colnames(cld_res) == "response"] <- "emmean"
    }
    
    # Now it will safely sort, regardless of what family it is
    cld_res <- cld_res[order(cld_res$emmean), ]
    
    # Dynamically find the grouping column (handles '.group', '.groups', 'group', etc.)
    group_col <- grep("group", colnames(cld_res), value = TRUE, ignore.case = TRUE)
    if(length(group_col) > 0) {
      group_letters <- trimws(cld_res[[group_col[1]]])
    } else {
      group_letters <- rep("?", nrow(cld_res)) # Fallback if cld fails completely
    }
    
    hierarchy_str <- paste(paste0(cld_res$LLM, " (", group_letters, ")"), collapse = " < ")
    
    # --- F. Store Results ---
    final_analysis_results[[index]] <- list(
      model = mod,
      omnibus = omnibus,
      emmeans_obj = emms
    )
    
    summary_table <- rbind(summary_table, data.frame(
      Index = index,
      Wald_Test = wald_string,
      R2_Marginal = r2_m,
      R2_Conditional = r2_c,
      EMM_Range = emm_range_str,
      Hierarchy_CLD = hierarchy_str,
      stringsAsFactors = FALSE
    ))
    
    message("  -> Success!")
    
  }, error = function(e) {
    message(" -> [CRITICAL Error] processing ", index, ": ", e$message)
  })
}

message("\nAnalysis complete! Displaying final summary table:\n")
print(summary_table)

write.csv(summary_table, "TFM_LLM_Linguistic_Summary_Table.csv", row.names = FALSE)
