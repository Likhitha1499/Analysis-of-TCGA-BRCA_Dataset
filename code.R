# Install required packages
install.packages(c("ggpubr", "survminer", "ggplot2", "survival", "dplyr", "factoextra"))

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(factoextra)

# Load gene expression data
gene_expression_data <- read.table("C:/Users/likhi/Downloads/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load clinical data
clinical_data <- read.table("C:/Users/likhi/Downloads/brca_tcga_pan_can_atlas_2018_clinical_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Prepare gene expression data
# Remove rows with missing Hugo_Symbol
gene_expression_data <- gene_expression_data %>%
  filter(!is.na(Hugo_Symbol))

# Handle missing values by replacing them with 0
gene_expression_data[is.na(gene_expression_data)] <- 0

# Convert data types to numeric (excluding Hugo_Symbol and Entrez_Gene_Id)
gene_expression_data_numeric <- gene_expression_data %>%
  select(-Hugo_Symbol, -Entrez_Gene_Id) %>%
  mutate(across(everything(), as.numeric))

# Add Hugo_Symbol back to the numeric data
gene_expression_data_numeric <- cbind(Hugo_Symbol = gene_expression_data$Hugo_Symbol, gene_expression_data_numeric)

# Sum the values across the columns (samples) for each gene (row)
gene_expression_data_numeric <- gene_expression_data_numeric %>%
  mutate(sum_expression = rowSums(select(., -Hugo_Symbol)))

# Filter out genes with low expression (sum of expression values < 10)
gene_expression_data_filtered <- gene_expression_data_numeric %>%
  filter(sum_expression >= 10)

# Remove the sum_expression column as it is no longer needed
gene_expression_data_filtered <- gene_expression_data_filtered %>%
  select(-sum_expression)

# Remove genes whose expression is zero in more than 50% of the samples
threshold <- (ncol(gene_expression_data_filtered) - 1) / 2  # subtract 1 for the Hugo_Symbol column
gene_expression_data_filtered <- gene_expression_data_filtered %>%
  filter(rowSums(select(., -Hugo_Symbol) == 0) <= threshold)

# Log transformation of the gene expression values (log2(x + 1))
gene_expression_data_log <- gene_expression_data_filtered %>%
  mutate(across(-Hugo_Symbol, ~ log2(.x + 1)))

# Z-score normalization
gene_expression_data_zscore <- gene_expression_data_log %>%
  mutate(across(-Hugo_Symbol, ~ ( .x - mean(.x) ) / sd(.x)))

# Transpose gene expression data
gene_expression_data_t <- t(gene_expression_data_zscore[-1])
colnames(gene_expression_data_t) <- gene_expression_data_zscore$Hugo_Symbol
gene_expression_data_t <- as.data.frame(gene_expression_data_t)
gene_expression_data_t$patient_id <- gsub("\\.", "-", colnames(gene_expression_data_zscore)[-1])

# Prepare clinical data
# Extract relevant columns for survival analysis
survival_data <- clinical_data %>%
  select(Patient.ID, Overall.Survival..Months., Overall.Survival.Status) %>%
  rename(patient_id = Patient.ID, time_to_event = Overall.Survival..Months., event_status = Overall.Survival.Status) %>%
  mutate(event_status = ifelse(event_status == "1:DECEASED", 1, 0))

# Clean patient IDs
# Ensure patient IDs are strings and consistent
survival_data$patient_id <- trimws(tolower(as.character(survival_data$patient_id)))
gene_expression_data_t$patient_id <- trimws(tolower(as.character(gene_expression_data_t$patient_id)))

# Remove the suffix "-01" from gene expression data patient IDs
gene_expression_data_t$patient_id <- gsub("-01$", "", gene_expression_data_t$patient_id)

# Merge gene expression and clinical data
merged_data <- merge(survival_data, gene_expression_data_t, by = "patient_id")

# Save the merged data for further analysis
write.table(merged_data, file = "merged_clinical_gene_expression_data.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# For feature selection
# Prepare data for Cox Proportional Hazards Model
survival_data_for_analysis <- merged_data %>%
  select(-patient_id) %>%
  as.data.frame()

# Fit Cox Proportional Hazards Model for each gene
p_values <- sapply(colnames(survival_data_for_analysis)[-(1:2)], function(gene) {
  tryCatch({
    cox_model <- coxph(Surv(time_to_event, event_status) ~ survival_data_for_analysis[[gene]], data = survival_data_for_analysis)
    summary(cox_model)$coefficients[, "Pr(>|z|)"]
  }, warning = function(w) {
    NA  # Assign NA if there's a warning
  }, error = function(e) {
    NA  # Assign NA if there's an error
  })
})

# Convert p-values to a data frame and sort by p-value
p_values_df <- data.frame(Gene = names(p_values), PValue = unlist(p_values))
top_genes <- p_values_df %>%
  arrange(PValue) %>%
  head(100)  # Top 100 genes with the smallest p-values

# Save top 100 genes
write.table(top_genes, file = "top_100_genes.csv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Top 100 genes saved.\n")
print(head(top_genes))

# Clustering
# Load top 100 genes
top_genes <- read.table("top_100_genes.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
selected_gene_names <- top_genes$Gene
gene_expression_top100 <- merged_data %>%
  select(patient_id, all_of(selected_gene_names))

# Determine the optimal number of clusters
gene_expression_matrix <- as.matrix(gene_expression_top100 %>% select(-patient_id))

# Elbow method to determine the optimal number of clusters
fviz_nbclust(gene_expression_matrix, kmeans, method = "wss") + 
  geom_vline(xintercept = 3, linetype = 2) + 
  labs(subtitle = "Elbow method")

# Silhouette method to determine the optimal number of clusters
fviz_nbclust(gene_expression_matrix, kmeans, method = "silhouette") + 
  labs(subtitle = "Silhouette method")

# Perform k-means clustering and assign each patient to a cluster
set.seed(123)
k <- 3
kmeans_result <- kmeans(gene_expression_matrix, centers = k, nstart = 25)

# Add the cluster assignments to the gene expression data
gene_expression_top100$cluster <- kmeans_result$cluster

# Save the clustering result
write.table(gene_expression_top100, file = "patient_clusters.csv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Clustering completed and saved.\n")
print(head(gene_expression_top100))

# Randomly select one gene from the top 100 genes
set.seed(123)
selected_gene <- selected_gene_names[66] #select 66th gene from the list of top 100
cat("Selected Gene:", selected_gene, "\n")

# Survival analysis in clusters
# Fit Cox Proportional Hazards Model within each cluster
coxph_results <- list()

for (cluster in unique(gene_expression_top100$cluster)) {
  cluster_data <- merged_data %>%
    filter(cluster == cluster)
  
  if (nrow(cluster_data) == 0 || sum(cluster_data$event_status) == 0) {
    cat("No events observed in Cluster", cluster, ". Skipping Cox Proportional Hazards Model.\n")
    next
  }
  
  # Fit the cox model
  cox_model <- coxph(as.formula(paste("Surv(time_to_event, event_status) ~", selected_gene)), data = cluster_data)
  coxph_results[[paste0("Cluster_", cluster)]] <- summary(cox_model)
}

# Print Cox Proportional Hazards Model results
for (cluster in names(coxph_results)) {
  cat(cluster, "Cox Proportional Hazards Model:\n")
  print(coxph_results[[cluster]])
}

# Kaplan-Meier survival analysis within each cluster
for (cluster in unique(gene_expression_top100$cluster)) {
  cluster_data <- merged_data %>%
    filter(cluster == cluster)
  
  if (nrow(cluster_data) == 0 || sum(cluster_data$event_status) == 0) {
    cat("No events observed in Cluster", cluster, ". Skipping Kaplan-Meier analysis.\n")
    next
  }
  
  # Split patients into high and low expression groups based on the median
  median_expression <- median(cluster_data[[selected_gene]], na.rm = TRUE)
  cluster_data$expression_group <- ifelse(cluster_data[[selected_gene]] > median_expression, "High", "Low")
  
  # Fit Kaplan-Meier survival curves
  km_fit <- survfit(Surv(time_to_event, event_status) ~ expression_group, data = cluster_data)
  
  # Plot Kaplan-Meier survival curves
  plot_title <- paste("Kaplan-Meier Survival Curve for Cluster", cluster)
  ggsurvplot(km_fit, data = cluster_data, pval = TRUE, conf.int = TRUE, 
             title = plot_title,
             xlab = "Time (Months)", ylab = "Survival Probability")
}

