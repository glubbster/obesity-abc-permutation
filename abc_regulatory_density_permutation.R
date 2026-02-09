# ABC Regulatory Density Permutation Test
# Compares regulatory density (genes/SNPs) between female and male obesity variants
# Author: [Your Name]
# Date: February 2026

# Load data ---------------------------------------------------------------
# Input files should have columns: CellType, Tissue, gene_name, variant_id, sex
# Replace these paths with your actual file paths
female_data <- read.csv("female_abc_mappings.csv")
male_data <- read.csv("male_abc_mappings.csv")

# Combine datasets
female_data$sex <- "Female"
male_data$sex <- "Male"
combined <- rbind(female_data, male_data)

# Calculate observed regulatory density -----------------------------------
tissues <- unique(combined$Tissue)
observed_results <- data.frame()

for (tissue in tissues) {
  tissue_data <- combined[combined$Tissue == tissue, ]
  
  # Female metrics
  female_subset <- tissue_data[tissue_data$sex == "Female", ]
  female_genes <- length(unique(female_subset$gene_name))
  female_snps <- length(unique(female_subset$variant_id))
  female_ratio <- female_genes / female_snps
  
  # Male metrics
  male_subset <- tissue_data[tissue_data$sex == "Male", ]
  male_genes <- length(unique(male_subset$gene_name))
  male_snps <- length(unique(male_subset$variant_id))
  male_ratio <- male_genes / male_snps
  
  # Store results
  observed_results <- rbind(observed_results, data.frame(
    Tissue = tissue,
    Set1genes = female_genes,
    Set1snps = female_snps,
    Set1ratio = female_ratio,
    Set2genes = male_genes,
    Set2snps = male_snps,
    Set2ratio = male_ratio,
    Realdifference = female_ratio - male_ratio
  ))
}

# Permutation test --------------------------------------------------------
set.seed(12345)  # For reproducibility
n_permutations <- 10000

permuted_diffs <- matrix(NA, nrow = length(tissues), ncol = n_permutations)
rownames(permuted_diffs) <- tissues

for (i in 1:n_permutations) {
  # Shuffle sex labels
  combined$sex_permuted <- sample(combined$sex)
  
  for (j in 1:length(tissues)) {
    tissue <- tissues[j]
    tissue_data <- combined[combined$Tissue == tissue, ]
    
    # Recalculate ratios with permuted labels
    perm_female <- tissue_data[tissue_data$sex_permuted == "Female", ]
    perm_male <- tissue_data[tissue_data$sex_permuted == "Male", ]
    
    female_ratio <- length(unique(perm_female$gene_name)) / length(unique(perm_female$variant_id))
    male_ratio <- length(unique(perm_male$gene_name)) / length(unique(perm_male$variant_id))
    
    permuted_diffs[j, i] <- female_ratio - male_ratio
  }
}

# Calculate empirical p-values --------------------------------------------
observed_results$Pempirical <- NA

for (i in 1:nrow(observed_results)) {
  tissue <- observed_results$Tissue[i]
  obs_diff <- observed_results$Realdifference[i]
  perm_diffs <- permuted_diffs[tissue, ]
  
  # Two-sided test
  observed_results$Pempirical[i] <- sum(abs(perm_diffs) >= abs(obs_diff)) / n_permutations
}

# FDR correction ----------------------------------------------------------
observed_results$Padj <- p.adjust(observed_results$Pempirical, method = "BH")

# Calculate fold change ---------------------------------------------------
observed_results$Foldchange <- observed_results$Set1ratio / observed_results$Set2ratio

# Export results ----------------------------------------------------------
observed_results$Comparison <- "FemaleObesity vs MaleObesity"
observed_results <- observed_results[, c("Comparison", "Tissue", "Set1genes", "Set1snps", 
                                          "Set1ratio", "Set2genes", "Set2snps", "Set2ratio", 
                                          "Realdifference", "Pempirical", "Padj", "Foldchange")]

write.csv(observed_results, "regulatory_density_permutation_results.csv", row.names = FALSE)

cat("Analysis complete. Results saved to regulatory_density_permutation_results.csv\n")
