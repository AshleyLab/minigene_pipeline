library(matrixStats)
library(tidyverse)
library(dplyr)
library(stringr)
library(readxl)
library(this.path)
library(yaml)

# number_obs_control: BNP_low, MYBPC3_high, UPS_nonGFP
# number_obs_experiment: BNP_high, MYBPC3_low, UPS_GFP

config_file <- file.path(this.dir(), "vesper_config.yaml")
config <- yaml::read_yaml(config_file)

protein <- config$protein
control_population <- config$control_population
experiment_population <- config$experimental_population

vesperii_output_dir = config$vesperii_output_dir
data_SNV_pairs <- read.csv(file.path(vesperii_output_dir, "mm2_hisat2_paired_counts", sprintf("%s_co_travelling_snv_pairs.csv", protein)))
data_snvs <- read.csv(file.path(vesperii_output_dir, "combined_t1t2_hisat2_count_files", sprintf("%s_T1T2_MAF_filtered.csv", protein)))
all_MYBPC3_vars_consequence <- read.csv(config$all_MYBPC3_vars_consequence)

# get positional values of Variant1 and Variant2 from data_SNV_pairs
data_SNV_pairs <- data_SNV_pairs %>% mutate(
  Variant1_pos = str_extract(Variant1, "-?\\d+"),
  Variant2_pos = str_extract(Variant2, "-?\\d+")
)
data_SNV_pairs$Variant1_pos <- as.numeric(data_SNV_pairs$Variant1_pos)
data_SNV_pairs$Variant2_pos <- as.numeric(data_SNV_pairs$Variant2_pos)

# Add REF and observed_nucleotide columns
data_snvs$REF <- sub(".*([A-Z])>.*", "\\1", data_snvs$Variant)
data_snvs$observed_nucleotide <- sub(".*>([A-Z])", "\\1", data_snvs$Variant)

# Remove reference allele rows and keep variant columns
all_consq_clean <- all_MYBPC3_vars_consequence %>%
  filter(Consequence != "REFERENCE ALLELE") %>%
  distinct()

## Map AA_consequence / Consequence / Pathogenic into data_snvs by Variant
data_snvs <- data_snvs %>%
  left_join(all_consq_clean %>% select(Variant, AA_consequence, Consequence, Pathogenic, PLP_BLB_VUS),
    by = "Variant"
  )


colnames(data_snvs) <- trimws(colnames(data_snvs))

colnames(data_snvs)[colnames(data_snvs) == "REF"] <- "reference_nucleotide"
colnames(data_snvs)[colnames(data_snvs) == "AA_consequence"] <- "Var_p"
colnames(data_snvs)[colnames(data_snvs) == "Consequence"] <- "Var_type"
colnames(data_snvs)[colnames(data_snvs) == "Pathogenic"] <- "P_LP"
colnames(data_snvs)[colnames(data_snvs) == "PLP_BLB_VUS"] <- "PLP_BLB_VUS"

proptest1_dir <- file.path(vesperii_output_dir, "R_Proptest_Results", "Proptest1")
# Create the proptest1 directory if it doesn't already exist
if (!dir.exists(proptest1_dir)) {
  dir.create(proptest1_dir, recursive = TRUE)
}

proptest2_dir <- file.path(vesperii_output_dir, "R_Proptest_Results", "Proptest2")

# Create the proptest2 directory if it doesn't already exist
if (!dir.exists(proptest2_dir)) {
  dir.create(proptest2_dir, recursive = TRUE)
}

merged_delta_proptest_dir <- file.path(vesperii_output_dir, "R_Proptest_Results", "Merged_Delta_Proptest")

# Create the merged_delta_proptest_dir directory if it doesn't already exist
if (!dir.exists(merged_delta_proptest_dir)) {
  dir.create(merged_delta_proptest_dir, recursive = TRUE)
}


# proptest control vs experiment
pt1_varname <- sprintf("proptest_results_%s_%s_vs_%s_%s", protein, control_population, protein, experiment_population)
assign(pt1_varname, data.frame(matrix(ncol = 10, nrow = 0)))

# remove rows where SNV is not found in hisat2 count files
data_snvs <- data_snvs %>% filter(!is.na(PLP_BLB_VUS))

# Build dynamic column names
count_con_col <- sprintf("count_%s_%s_T1T2", protein, control_population)
count_exp_col <- sprintf("count_%s_%s_T1T2", protein, experiment_population)

depth_con_col <- sprintf("depth_%s_%s_T1T2", protein, control_population)
depth_exp_col <- sprintf("depth_%s_%s_T1T2", protein, experiment_population)

# proptest1 control vs experiment (BNP_low vs BNP_high, MYBPC3_high vs MYBPC3_low, UPS_nonGFP vs UPS_GFP)
proptest1_df <- get(pt1_varname)

for (i in 1:nrow(data_snvs)) {
  enrichment <- matrix()
  data_snvs$number_obs_exp<- data_snvs[[count_exp_col]]
  data_snvs$number_obs_con <- data_snvs[[count_con_col]]
  data_snvs$depth_exp <- data_snvs[[depth_exp_col]]
  data_snvs$depth_con <- data_snvs[[depth_con_col]]
  results <- prop.test(x=c(data_snvs$number_obs_exp[i],data_snvs$number_obs_con[i]),n=c(data_snvs$depth_exp[i],data_snvs$depth_con[i])) #apply prop.test to row
  single_result <- c(data_snvs$Variant[i],data_snvs$reference_nucleotide[i], data_snvs$observed_nucleotide[i], data_snvs$Var_type[i], data_snvs$PLP_BLB_VUS[i], results$estimate, results$conf.int, results$p.value) #extract test results
  single_result_row <- t(single_result)
  proptest1_df <- rbind(proptest1_df,single_result_row) #add to proptest results df
}           

colnames(proptest1_df) <- c("Variant", "reference", "observed_nucleotide", "Var_type", "PLP_BLB_VUS", "proportion_exp", "proportion_con", "CI95_lower", "CI95_upper", "p_value")

proptest1_df$proportion_exp <- as.numeric(proptest1_df$proportion_exp)
proptest1_df$proportion_con <- as.numeric(proptest1_df$proportion_con)
proptest1_df$CI95_upper <- as.numeric(proptest1_df$CI95_upper)
proptest1_df$CI95_lower <- as.numeric(proptest1_df$CI95_lower)
proptest1_df$p_value <- as.numeric(proptest1_df$p_value)

proptest1_df$exp_enrichment <- (proptest1_df$proportion_exp/proptest1_df$proportion_con)
proptest1_df$neglogP <- -(log10(proptest1_df$p_value))
proptest1_df$Var_p <-data_snvs$Var_p
proptest1_df$Var_type <-data_snvs$Var_type
proptest1_df$P_LP <-data_snvs$P_LP
# proptest1_df$PLP_Stopgain_Syn_VUS <- data_snvs$PLP_Stopgain_Syn_VUS

output_filename <- sprintf("proptest1_results_%s_%s_vs_%s_%s.csv", protein, control_population, protein, experiment_population)
output_filepath <- file.path(proptest1_dir, output_filename)
write.csv(proptest1_df, file = output_filepath, row.names = FALSE)

proptest1_var_type_plot <- ggplot(proptest1_df, aes(x = Var_type, y=exp_enrichment, fill=Var_type)) +
  geom_violin(alpha=0.6, color=NA) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
#   stat_summary(fun = median, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
#   stat_summary(fun = median, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values=c("firebrick3", "palevioletred1", "plum3", "lightsteelblue2", "#f258ad")) +
  xlab("Variant Type") +
  ylab("Experiment Enrichment Score") +
#   scale_y_continuous(breaks = seq(0, max(proptest1_df$exp_enrichment, na.rm = TRUE)+1, by = 1)) +
  theme_classic() +
  theme(
    # panel.grid.major = element_line(color = "grey90", size = 0.2),
    # panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(-9, 0, 0, 0),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  ) +
  ggtitle(sprintf("Experiment Enrichment Score of %s_%s vs %s_%s variants", protein, control_population, protein, experiment_population))
proptest1_var_type_plot_path <- file.path(proptest1_dir, sprintf("Exp_enrichment_CX_SNVs_%s_%s_vs_%s_Var_type.pdf", protein, control_population, protein, experiment_population))
ggsave(proptest1_var_type_plot_path, plot=proptest1_var_type_plot, width=8, height=6)

proptest1_plp_plot <- ggplot(proptest1_df, aes(x = PLP_BLB_VUS, y=exp_enrichment, fill=PLP_BLB_VUS)) +
  geom_violin(alpha=0.6, color=NA) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
#   stat_summary(fun = median, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
#   stat_summary(fun = median, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values=c("firebrick3", "palevioletred1", "plum3")) +
  xlab("Variant Type") +
  ylab("Experiment Enrichment Score") +
#   scale_y_continuous(breaks = seq(0, max(proptest1_df$exp_enrichment, na.rm = TRUE)+1, by = 1)) +
  theme_classic() +
  theme(
    # panel.grid.major = element_line(color = "grey90", size = 0.2),
    # panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(-9, 0, 0, 0),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  ) +
  ggtitle(sprintf("Experiment Enrichment Score of %s_%s vs %s_%s variants", protein, control_population, protein, experiment_population))
proptest1_plp_plot_path <- file.path(proptest1_dir, sprintf("Exp_enrichment_CX_SNVs_%s_%s_vs_%s_PLP.pdf", protein, control_population, protein, experiment_population))
ggsave(proptest1_plp_plot_path, plot=proptest1_plp_plot, width=8, height=6)





# #######end proptests for different phenotypes, now individually subtract each pair 
# # and re-do prop-test and visualize the spread of new values
# # take the statistically significant enrichment value closest to 1 to correct for any effect of other variants
data_SNV_pairs[[sprintf("SNVcount_Var1_Count_%s_%s", protein, experiment_population)]] <- ""
data_SNV_pairs[[sprintf("SNVcount_Var1_Count_%s_%s", protein, control_population)]] <- ""

data_SNV_pairs[[sprintf("SNVdepth_Var1_Depth_%s_%s", protein, experiment_population)]] <- ""
data_SNV_pairs[[sprintf("SNVdepth_Var1_Depth_%s_%s", protein, control_population)]] <- ""

# Loop through rows and fill in values
for (i in 1:nrow(data_SNV_pairs)) {
  df <- data_snvs[which(data_snvs$Variant == data_SNV_pairs$Variant1[i]), ]
  if (nrow(df) == 0) {
    next
  }

  var1_count_exp_col <- sprintf("count_%s_%s_T1T2", protein, experiment_population)
  var1_count_con_col <- sprintf("count_%s_%s_T1T2", protein, control_population)
  var1_depth_exp_col <- sprintf("depth_%s_%s_T1T2", protein, experiment_population)
  var1_depth_con_col <- sprintf("depth_%s_%s_T1T2", protein, control_population)

  data_SNV_pairs[[sprintf("SNVcount_Var1_Count_%s_%s", protein, experiment_population)]][i] <- df[[var1_count_exp_col]]
  data_SNV_pairs[[sprintf("SNVdepth_Var1_Depth_%s_%s", protein, experiment_population)]][i] <- df[[var1_depth_exp_col]]
  data_SNV_pairs[[sprintf("SNVcount_Var1_Count_%s_%s", protein, control_population)]][i] <- df[[var1_count_con_col]]
  data_SNV_pairs[[sprintf("SNVdepth_Var1_Depth_%s_%s", protein, control_population)]][i] <- df[[var1_depth_con_col]]
}



### do the same for variant 2 of variant pair ###
data_SNV_pairs[[sprintf("SNVcount_Var2_Count_%s_%s", protein, experiment_population)]] <- ""
data_SNV_pairs[[sprintf("SNVcount_Var2_Count_%s_%s", protein, control_population)]] <- ""

data_SNV_pairs[[sprintf("SNVdepth_Var2_Depth_%s_%s", protein, experiment_population)]] <- ""
data_SNV_pairs[[sprintf("SNVdepth_Var2_Depth_%s_%s", protein, control_population)]] <- ""

# Loop through rows and fill in values
for (i in 1:nrow(data_SNV_pairs)) {
  df <- data_snvs[which(data_snvs$Variant == data_SNV_pairs$Variant2[i]), ]
  if (nrow(df) == 0) {
    next
  }

  var2_count_exp_col <- sprintf("count_%s_%s_T1T2", protein, experiment_population)
  var2_count_con_col <- sprintf("count_%s_%s_T1T2", protein, control_population)
  var2_depth_exp_col <- sprintf("depth_%s_%s_T1T2", protein, experiment_population)
  var2_depth_con_col <- sprintf("depth_%s_%s_T1T2", protein, control_population)

  data_SNV_pairs[[sprintf("SNVcount_Var2_Count_%s_%s", protein, experiment_population)]][i] <- df[[var2_count_exp_col]]
  data_SNV_pairs[[sprintf("SNVdepth_Var2_Depth_%s_%s", protein, experiment_population)]][i] <- df[[var2_depth_exp_col]]
  data_SNV_pairs[[sprintf("SNVcount_Var2_Count_%s_%s", protein, control_population)]][i] <- df[[var2_count_con_col]]
  data_SNV_pairs[[sprintf("SNVdepth_Var2_Depth_%s_%s", protein, control_population)]][i] <- df[[var2_depth_con_col]]
}

##Now recalculate the counts in each group per variant when counts of reads with variant pairs are subtracted
### get delta count for variant1 experimental condition ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, experiment_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVcount_Var1_Count_%s_%s", protein, experiment_population)]]) - 
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant1_mpileup_count_%s_%s", protein, experiment_population)]])
### if delta count is 0, change to 1 to avoid errors in prop.test ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, experiment_population)]] <- ifelse(data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, experiment_population)]] == 0,1,
                                                                                        data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, experiment_population)]])
### get delta depth for variant 1 experimental condition ###
data_SNV_pairs[[sprintf("deltaSNVdepth_Var1_%s_%s", protein, experiment_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVdepth_Var1_Depth_%s_%s", protein, experiment_population)]]) -
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant1_mpileup_count_%s_%s", protein, experiment_population)]])


### do the same for variant 1 control condition ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, control_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVcount_Var1_Count_%s_%s", protein, control_population)]]) - 
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant1_mpileup_count_%s_%s", protein, control_population)]])
### if delta count is 0, change to 1 to avoid errors in prop.test ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, control_population)]] <- ifelse(data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, control_population)]] == 0,1,
                                                                                        data_SNV_pairs[[sprintf("deltaSNVcount_Var1_%s_%s", protein, control_population)]])
### get delta depth for variant 1 control condition ###
data_SNV_pairs[[sprintf("deltaSNVdepth_Var1_%s_%s", protein, control_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVdepth_Var1_Depth_%s_%s", protein, control_population)]]) -
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant1_mpileup_count_%s_%s", protein, control_population)]])


### do the same for variant 2 experimental condition ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, experiment_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVcount_Var2_Count_%s_%s", protein, experiment_population)]]) - 
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant2_mpileup_count_%s_%s", protein, experiment_population)]])
### if delta count is 0, change to 1 to avoid errors in prop.test ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, experiment_population)]] <- ifelse(data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, experiment_population)]] == 0,1,
                                                                                        data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, experiment_population)]])
### get delta depth for variant 2 condition 1 ###
data_SNV_pairs[[sprintf("deltaSNVdepth_Var2_%s_%s", protein, experiment_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVdepth_Var2_Depth_%s_%s", protein, experiment_population)]]) -
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant2_mpileup_count_%s_%s", protein, experiment_population)]])


### do the same for variant 2 control condition ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, control_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVcount_Var2_Count_%s_%s", protein, control_population)]]) - 
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant2_mpileup_count_%s_%s", protein, control_population)]])
### if delta count is 0, change to 1 to avoid errors in prop.test ###
data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, control_population)]] <- ifelse(data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, control_population)]] == 0,1,
                                                                                        data_SNV_pairs[[sprintf("deltaSNVcount_Var2_%s_%s", protein, control_population)]])
### get delta depth for variant 2 control condition ###
data_SNV_pairs[[sprintf("deltaSNVdepth_Var2_%s_%s", protein, control_population)]] <- (as.numeric(data_SNV_pairs[[sprintf("SNVdepth_Var2_Depth_%s_%s", protein, control_population)]]) -
                                                                                              data_SNV_pairs[[sprintf("calculated_Variant2_mpileup_count_%s_%s", protein, control_population)]])

### save updated data_SNV_pairs with delta counts and depths in merged_delta_proptest_dir ###
# output_filepath_delta <- file.path(merged_delta_proptest_dir, sprintf("delta_counts_depths_SNV_pairs_%s_T1T2.csv", protein, tile))
# write.csv(data_SNV_pairs, file = output_filepath_delta, row.names = FALSE)

##now do proptests across variants with delta counts for each phenotype

### proptest for variant 1 control vs experiment ###
proptest2_variant1_df <- sprintf("proptest2_variant1_%s_%s_vs_%s_%s", protein, control_population, protein, experiment_population)
assign(proptest2_variant1_df, data.frame(matrix(ncol = 14, nrow = 0)))

for (i in 1:nrow(data_SNV_pairs)) {
  proptest2_var1_count_control_col <- sprintf("deltaSNVcount_Var1_%s_%s", protein, control_population)
  proptest2_var1_count_exp_col <- sprintf("deltaSNVcount_Var1_%s_%s", protein, experiment_population)
  proptest2_var1_depth_control_col <- sprintf("deltaSNVdepth_Var1_%s_%s", protein, control_population)
  proptest2_var1_depth_exp_col <- sprintf("deltaSNVdepth_Var1_%s_%s", protein, experiment_population)

  if (is.na(data_SNV_pairs[[proptest2_var1_count_control_col]][i]) | is.na(data_SNV_pairs[[proptest2_var1_count_exp_col]][i])) {
    next
  }

  results <- prop.test(
    x = c(data_SNV_pairs[[proptest2_var1_count_exp_col]][i], data_SNV_pairs[[proptest2_var1_count_control_col]][i]),
    n = c(data_SNV_pairs[[proptest2_var1_depth_exp_col]][i], data_SNV_pairs[[proptest2_var1_depth_control_col]][i])
  )
  single_result <- c(data_SNV_pairs$Variant1[i], data_SNV_pairs$Variant2[i], results$estimate, results$conf.int, results$p.value) #extract test results
  single_result_row <- t(single_result)

  assign(proptest2_variant1_df, rbind(get(proptest2_variant1_df), single_result_row))

}

# Set column names
tmp_proptest2_variant1_df <- get(proptest2_variant1_df)
colnames(tmp_proptest2_variant1_df) <- c(
  "Variant1",
  "Variant2",
  "proportion_exp",
  "proportion_con",
  "CI95_lower",
  "CI95_upper",
  "p_value"
)
assign(proptest2_variant1_df, tmp_proptest2_variant1_df)


# Convert to numeric and add enrichment calculations
pt2_v1_df <- get(proptest2_variant1_df)
pt2_v1_df$proportion_exp <- as.numeric(pt2_v1_df$proportion_exp)
pt2_v1_df$proportion_con <- as.numeric(pt2_v1_df$proportion_con)
pt2_v1_df$exp_enrichment <- pt2_v1_df$proportion_exp / pt2_v1_df$proportion_con
pt2_v1_df$log2_exp_enrichment <- log2(pt2_v1_df$exp_enrichment)
pt2_v1_df$CI95_upper <- as.numeric(pt2_v1_df$CI95_upper)
pt2_v1_df$CI95_lower <- as.numeric(pt2_v1_df$CI95_lower)
pt2_v1_df$p_value <- as.numeric(pt2_v1_df$p_value)
pt2_v1_df$neglogP <- -log10(pt2_v1_df$p_value)

# Save the updated df back to original name
assign(proptest2_variant1_df, pt2_v1_df)

# Write to CSV
output_file <- sprintf("proptest2_variant1_%s_%s_vs_%s_%s.csv", protein, control_population, protein, experiment_population)
write.csv(pt2_v1_df, file.path(proptest2_dir, output_file), row.names = FALSE)


### proptest for variant 2 neg vs pos ###
proptest2_variant2_df <- sprintf("proptest2_variant2_%s_%s_vs_%s_%s", protein, control_population, protein, experiment_population)
assign(proptest2_variant2_df, data.frame(matrix(ncol = 14, nrow = 0)))

for (i in 1:nrow(data_SNV_pairs)) {
  proptest2_var2_count_control_col <- sprintf("deltaSNVcount_Var2_%s_%s", protein, control_population)
  proptest2_var2_count_exp_col <- sprintf("deltaSNVcount_Var2_%s_%s", protein, experiment_population)
  proptest2_var2_depth_control_col <- sprintf("deltaSNVdepth_Var2_%s_%s", protein, control_population)
  proptest2_var2_depth_exp_col <- sprintf("deltaSNVdepth_Var2_%s_%s", protein, experiment_population)

  if (is.na(data_SNV_pairs[[proptest2_var2_count_control_col]][i]) | is.na(data_SNV_pairs[[proptest2_var2_count_exp_col]][i])) {
    next
  }
  results <- prop.test(
    x = c(data_SNV_pairs[[proptest2_var2_count_exp_col]][i], data_SNV_pairs[[proptest2_var2_count_control_col]][i]),
    n = c(data_SNV_pairs[[proptest2_var2_depth_exp_col]][i], data_SNV_pairs[[proptest2_var2_depth_control_col]][i])
  )
  single_result <- c(data_SNV_pairs$Variant2[i], data_SNV_pairs$Variant1[i], results$estimate, results$conf.int, results$p.value) #extract test results
  single_result_row <- t(single_result)

  assign(proptest2_variant2_df, rbind(get(proptest2_variant2_df), single_result_row))
}

# Set column names
tmp_proptest2_variant2_df <- get(proptest2_variant2_df)
colnames(tmp_proptest2_variant2_df) <- c(
  "Variant1",
  "Variant2",
  "proportion_exp",
  "proportion_con",
  "CI95_lower",
  "CI95_upper",
  "p_value"
)
assign(proptest2_variant2_df, tmp_proptest2_variant2_df)

# Convert to numeric and add enrichment calculations
pt2_v2_df <- get(proptest2_variant2_df)
pt2_v2_df$proportion_exp <- as.numeric(pt2_v2_df$proportion_exp)
pt2_v2_df$proportion_con <- as.numeric(pt2_v2_df$proportion_con)
pt2_v2_df$exp_enrichment <- pt2_v2_df$proportion_exp / pt2_v2_df$proportion_con
pt2_v2_df$log2_exp_enrichment <- log2(pt2_v2_df$exp_enrichment)
pt2_v2_df$CI95_upper <- as.numeric(pt2_v2_df$CI95_upper)
pt2_v2_df$CI95_lower <- as.numeric(pt2_v2_df$CI95_lower)
pt2_v2_df$p_value <- as.numeric(pt2_v2_df$p_value)
pt2_v2_df$neglogP <- -log10(pt2_v2_df$p_value)

# Save the updated df back to original name
assign(proptest2_variant2_df, pt2_v2_df)

# Write to CSV
output_file <- sprintf("proptest2_variant2_%s_%s_vs_%s_%s.csv", protein, control_population, protein, experiment_population)
write.csv(pt2_v2_df, file.path(proptest2_dir, output_file), row.names = FALSE)

### Merge proptest2 results for Var1 and Var2 to get a combined view ###
delta_proptest_df <- sprintf(
  "delta_proptests_combined_%s_%s_vs_%s_%s",
  protein, control_population, protein, experiment_population
)

pt2_v1_slim <- pt2_v1_df[, c("Variant1", "Variant2", "exp_enrichment", "p_value")]
pt2_v2_slim <- pt2_v2_df[, c("Variant1", "Variant2", "exp_enrichment", "p_value")]

# Now rbind is happy because names match
delta_proptest_df <- rbind(pt2_v1_slim, pt2_v2_slim)

# Keep only rows where either Variant1 or Variant2 exists in data_snvs
delta_proptest_df <- delta_proptest_df[
  delta_proptest_df$Variant1 %in% data_snvs$Variant |
  delta_proptest_df$Variant2 %in% data_snvs$Variants, 
]


merged_delta_proptest_df <- 
  delta_proptest_df %>%
  group_by(Variant1)  %>%
  summarize(delta_enrichment_min = min(exp_enrichment),
            delta_enrichment_min_p_value = p_value[which.min(exp_enrichment)],
            delta_enrichment_max = max(exp_enrichment),
            delta_enrichment_max_p_value = p_value[which.max(exp_enrichment)],
            delta_enrichment_mean = mean(exp_enrichment),
            delta_enrichment_median = median(exp_enrichment),
            delta_enrichment_SD = sd(exp_enrichment),
            delta_enrichment_closest_to_1 = exp_enrichment[which.min(abs(exp_enrichment-1))],
            delta_enrichment_closest_to_1_p_value = p_value[which.min(abs(exp_enrichment-1))]
  ) 

#annotate by adding back to data_snvs
merged_delta_proptest_df$reference_nucleotide <- ""
merged_delta_proptest_df$observed_nucleotide <- ""
merged_delta_proptest_df$Var_p <- ""
merged_delta_proptest_df$Var_type <- "" 
merged_delta_proptest_df$PLP_BLB_VUS <-""
merged_delta_proptest_df$P_LP <- ""

for (i in 1:nrow(merged_delta_proptest_df)) {
  df <- proptest1_df[which(proptest1_df$Variant==merged_delta_proptest_df$Variant1[i]),]
  if (nrow(df)==0) {
    next }
  merged_delta_proptest_df$reference_nucleotide[i] <- df$reference
  merged_delta_proptest_df$observed_nucleotide[i] <- df$observed_nucleotide
  merged_delta_proptest_df$Var_p[i] <- df$Var_p
  merged_delta_proptest_df$Var_type[i] <- df$Var_type
  merged_delta_proptest_df$PLP_BLB_VUS[i] <-df$PLP_BLB_VUS
  merged_delta_proptest_df$P_LP[i] <- df$P_LP
}


# #_#_#_#_#_#_#_#

# # #calculate FDR p_adj
merged_delta_proptest_df$p_adj_ct1 <- p.adjust(merged_delta_proptest_df$delta_enrichment_closest_to_1_p_value,"fdr")
merged_delta_proptest_df$p_adj_min <- p.adjust(merged_delta_proptest_df$delta_enrichment_min_p_value,"fdr")
merged_delta_proptest_df$p_adj_max <- p.adjust(merged_delta_proptest_df$delta_enrichment_max_p_value,"fdr")
write.csv(merged_delta_proptest_df, file.path(merged_delta_proptest_dir, sprintf("Merged_delta_proptest_by_SNV_%s_%s_vs_%s_%s.csv", protein, control_population, protein, experiment_population)), row.names = FALSE)

merged_delta_var_type_proptest_df_plot <- ggplot(merged_delta_proptest_df, aes(x = Var_type, y=delta_enrichment_closest_to_1, fill=Var_type)) +
  geom_violin(alpha=0.6, color=NA) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
#   stat_summary(fun = median, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
#   stat_summary(fun = median, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values=c("firebrick3", "palevioletred1", "plum3", "lightsteelblue2", "#f258ad")) +
  xlab("Variant Type") +
  ylab("Delta Enrichment CT1") +
#   scale_y_continuous(breaks = seq(0, max(merged_delta_proptest_df$delta_enrichment_closest_to_1, na.rm = TRUE)+1, by = 1)) +
  theme_classic() +
  theme(
    # panel.grid.major = element_line(color = "grey90", size = 0.2),
    # panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(-9, 0, 0, 0),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  ) +
  ggtitle(sprintf("Delta Enrichment CT1 of SNVs in %s_%s vs %s_%s", protein, control_population, protein, experiment_population))    
merged_delta_var_type_proptest_df_plot_path <- file.path(merged_delta_proptest_dir, sprintf("Merged_delta_proptest_by_SNV_%s_%s_vs_%s_%s_ct1_var_type.pdf", protein, control_population, protein, experiment_population))
ggsave(merged_delta_var_type_proptest_df_plot_path, plot=merged_delta_var_type_proptest_df_plot, width=8, height=6) 


merged_delta_plp_proptest_df_plot <- ggplot(merged_delta_proptest_df, aes(x = PLP_BLB_VUS, y=delta_enrichment_closest_to_1, fill=PLP_BLB_VUS)) +
  geom_violin(alpha=0.6, color=NA) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
#   stat_summary(fun = median, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
#   stat_summary(fun = median, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values=c("firebrick3", "palevioletred1", "plum3")) +
  xlab("Variant Type") +
  ylab("Delta Enrichment CT1") +
#   scale_y_continuous(breaks = seq(0, max(merged_delta_proptest_df$delta_enrichment_closest_to_1, na.rm = TRUE)+1, by = 1)) +
  theme_classic() +
  theme(
    # panel.grid.major = element_line(color = "grey90", size = 0.2),
    # panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(-9, 0, 0, 0),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  ) +
  ggtitle(sprintf("Delta Enrichment CT1 of SNVs in %s_%s vs %s_%s", protein, control_population, protein, experiment_population))    
merged_delta_plp_proptest_df_plot_path <- file.path(merged_delta_proptest_dir, sprintf("Merged_delta_proptest_by_SNV_%s_%s_vs_%s_%s_ct1_PLP.pdf", protein, control_population, protein, experiment_population))
ggsave(merged_delta_plp_proptest_df_plot_path, plot=merged_delta_plp_proptest_df_plot, width=8, height=6) 

###### median version plots ##########
median_merged_delta_var_type_proptest_df_plot <- ggplot(merged_delta_proptest_df, aes(x = Var_type, y=delta_enrichment_median, fill=Var_type)) +
  geom_violin(alpha=0.6, color=NA) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
#   stat_summary(fun = median, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
#   stat_summary(fun = median, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values=c("firebrick3", "palevioletred1", "plum3", "lightsteelblue2", "#f258ad")) +
  xlab("Variant Type") +
  ylab("Delta Enrichment Median") +
#   scale_y_continuous(breaks = seq(0, max(merged_delta_proptest_df$delta_enrichment_median, na.rm = TRUE)+1, by = 1)) +
  theme_classic() +
  theme(
    # panel.grid.major = element_line(color = "grey90", size = 0.2),
    # panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(-9, 0, 0, 0),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  ) +
  ggtitle(sprintf("Delta Enrichment Median of SNVs in %s_%s vs %s_%s", protein, control_population, protein, experiment_population))    
median_merged_delta_var_type_proptest_df_plot_path <- file.path(merged_delta_proptest_dir, sprintf("Median_Merged_delta_proptest_by_SNV_%s_%s_vs_%s_%s_ct1_var_type.pdf", protein, control_population, protein, experiment_population))
ggsave(median_merged_delta_var_type_proptest_df_plot_path, plot=median_merged_delta_var_type_proptest_df_plot, width=8, height=6) 


median_merged_delta_plp_proptest_df_plot <- ggplot(merged_delta_proptest_df, aes(x = PLP_BLB_VUS, y=delta_enrichment_median, fill=PLP_BLB_VUS)) +
  geom_violin(alpha=0.6, color=NA) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
#   stat_summary(fun = median, geom = "point", shape = 9, size = 2.5, color = "violetred4", show.legend = FALSE) +
#   stat_summary(fun = median, geom = "text", vjust = 0, hjust = -0.7, aes(label = round(..y.., 2)), color = "violetred4", size = 2.5) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values=c("firebrick3", "palevioletred1", "plum3")) +
  xlab("Variant Type") +
  ylab("Delta Enrichment Median") +
#   scale_y_continuous(breaks = seq(0, max(merged_delta_proptest_df$delta_enrichment_median, na.rm = TRUE)+1, by = 1)) +
  theme_classic() +
  theme(
    # panel.grid.major = element_line(color = "grey90", size = 0.2),
    # panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(-9, 0, 0, 0),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  ) +
  ggtitle(sprintf("Delta Enrichment Median of SNVs in %s_%s vs %s_%s", protein, control_population, protein, experiment_population))    
median_merged_delta_plp_proptest_df_plot_path <- file.path(merged_delta_proptest_dir, sprintf("Median_Merged_delta_proptest_by_SNV_%s_%s_vs_%s_%s_ct1_PLP.pdf", protein, control_population, protein, experiment_population))
ggsave(median_merged_delta_plp_proptest_df_plot_path, plot=median_merged_delta_plp_proptest_df_plot, width=8, height=6) 
