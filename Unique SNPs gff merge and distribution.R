## step1:Load Required Packages
#dplyr, tidyr, ggplot2 already installed
# Install Bioconductor package manager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install rtracklayer from Bioconductor
BiocManager::install("rtracklayer")

# Load the package
library(rtracklayer)

install.packages(c("R.utils"), dependencies = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(rtracklayer)
library(R.utils)

## step2: Unzip and Load the GFF File
# Path to the GFF file (Compressed)
gff_compressed <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/GCF_015228065.2_NSTDA_Pmon_1_genomic.gff.gz"

# Unzip the GFF file if not already unzipped
gff_unzipped <- sub(".gz$", "", gff_compressed)  # Remove .gz extension
if (!file.exists(gff_unzipped)) {
  gunzip(gff_compressed, overwrite = FALSE)
  message("✅ GFF file unzipped successfully.")
}

# Load GFF file
gff_data <- import(gff_unzipped)

# Convert to DataFrame
gff_df <- as.data.frame(gff_data)

## step3: Extract Relevant Gene Information

# Keep only gene-related data
genes_data <- gff_df %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, ID)  # Selecting Chromosome, Start, End, and Gene ID

# Rename columns for clarity
colnames(genes_data) <- c("CHROM_Gff", "Gene_Start", "Gene_End", "Gene_ID")

# Display extracted gene data
head(genes_data)

## step4: Merge SNP Data with Gene Annotations

# Merge SNPs from Control Group with Gene Data
merged_control <- control_only %>%
  rowwise() %>%
  mutate(Gene_ID = list(genes_data$Gene_ID[which(genes_data$Gene_Start <= POS & 
                                                   genes_data$Gene_End >= POS & 
                                                   genes_data$CHROM_Gff == CHROM)])) %>%
  unnest(cols = c(Gene_ID)) %>%
  filter(!is.na(Gene_ID))  # Remove unmatched SNPs

# Merge SNPs from Survived Group with Gene Data
merged_survived <- survived_only %>%
  rowwise() %>%
  mutate(Gene_ID = list(genes_data$Gene_ID[which(genes_data$Gene_Start <= POS & 
                                                   genes_data$Gene_End >= POS & 
                                                   genes_data$CHROM_Gff == CHROM)])) %>%
  unnest(cols = c(Gene_ID)) %>%
  filter(!is.na(Gene_ID))  # Remove unmatched SNPs

# Save the merged SNP-Gene data as CSV files
write.csv(merged_control, "Unique_Control_SNPs.csv", row.names = FALSE)
write.csv(merged_survived, "Unique_Survived_SNPs.csv", row.names = FALSE)
message("✅ SNP tables saved as CSV files.")

# =============================
# 5️⃣ Create SNP Distribution Table
# =============================

# SNP Distribution for Control Group
snp_distribution_control <- merged_control %>%
  group_by(Gene_ID) %>%
  summarise(SNP_count = n()) %>%
  mutate(SNP_bin = cut(SNP_count,
                       breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, Inf),
                       labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                  "11-15", "16-20", "21-25", "26-30", ">30"),
                       right = FALSE)) %>%
  group_by(SNP_bin) %>%
  summarise(Number_of_genes = n()) %>%
  mutate(Genotype = "Control")

# SNP Distribution for Survived Group
snp_distribution_survived <- merged_survived %>%
  group_by(Gene_ID) %>%
  summarise(SNP_count = n()) %>%
  mutate(SNP_bin = cut(SNP_count,
                       breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, Inf),
                       labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                  "11-15", "16-20", "21-25", "26-30", ">30"),
                       right = FALSE)) %>%
  group_by(SNP_bin) %>%
  summarise(Number_of_genes = n()) %>%
  mutate(Genotype = "Survived")

# Combine Both Groups for Final Plot
snp_distribution_combined <- bind_rows(snp_distribution_control, snp_distribution_survived)

# =============================
# 6️⃣ Plot SNP Distribution Per Gene (Publication-Ready)
# =============================

# Create the bar plot
snp_plot <- ggplot(snp_distribution_combined, aes(x = SNP_bin, y = Number_of_genes, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Border for clarity
  labs(title = "SNP Distribution per Gene (Control vs Survived)",
       x = "Number of SNPs per Gene",
       y = "Number of Genes") +
  theme_minimal(base_size = 15) +  # Clean publication theme
  scale_fill_manual(values = c("Control" = "#5E81AC", "Survived" = "#BF616A")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels

# Save plot as PNG
ggsave("SNPs_per_Gene_Comparison.png", plot = snp_plot, width = 8, height = 6, dpi = 300)

# Display success message
message("✅ Plot saved successfully: SNPs_per_Gene_Comparison.png")

# Display the plot in RStudio
print(snp_plot)
