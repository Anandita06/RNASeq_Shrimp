## Step1: Extract Unique SNPs and plot them
# Define required packages
packages <- c("vcfR", "dplyr", "tidyr", "ggplot2")

# Install any missing packages
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)

# Load packages
lapply(packages, library, character.only = TRUE)

# Function to read and extract data from VCF files
read_vcf_data <- function(file_path) {
  vcf <- read.vcfR(file_path)
  data <- data.frame(
    CHROM = vcf@fix[,1],  # Chromosome/Scaffold
    POS = vcf@fix[,2],    # Position
    REF = vcf@fix[,4],    # Reference Allele
    ALT = vcf@fix[,5],    # Alternate Allele
    ANN = vcf@fix[,8]     # Annotation
  )
  return(data)
}

# File paths
control_vcf_path <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/NEW_DATA-20250215T055250Z-001/NEW_DATA/CONTROL_36_SNP_annotation/36_SNP.ann.vcf"
survived_vcf_path <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/NEW_DATA-20250215T055250Z-001/NEW_DATA/SURVIVED_35_SNP_annotation/35_SNP.ann.vcf"

# Read and extract data
control_data <- read_vcf_data(control_vcf_path)
survived_data <- read_vcf_data(survived_vcf_path)

# Extract unique variants while keeping ANN
control_variants_new <- control_data %>%
  group_by(CHROM, POS, REF, ALT) %>%
  summarise(ANN = paste(unique(ANN), collapse = ";"), .groups = "drop")  # Preserve ANN

survived_variants_new <- survived_data %>%
  group_by(CHROM, POS, REF, ALT) %>%
  summarise(ANN = paste(unique(ANN), collapse = ";"), .groups = "drop")  # Preserve ANN

# Identify unique and common variants while keeping ANN
control_only_new <- anti_join(control_variants_new, survived_variants_new, by = c("CHROM", "POS", "REF", "ALT")) 
survived_only_new <- anti_join(survived_variants_new, control_variants_new, by = c("CHROM", "POS", "REF", "ALT")) 
common_variants_new <- inner_join(control_variants_new, survived_variants_new, by = c("CHROM", "POS", "REF", "ALT"), suffix = c("_control", "_survived"))

# Count variants in each category
variant_counts_new <- data.frame(
  Category = c("Control Only", "Survived Only", "Common"),
  Count = c(nrow(control_only_new), nrow(survived_only_new), nrow(common_variants_new))
)

# Plot with labels and high-quality colors
SNP_counts_barplot_new <- ggplot(variant_counts_new, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  # Subtle border for clarity
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +      # Add count labels
  theme_minimal(base_size = 15) +                              # Adjust font size
  labs(title = "Unique and Common Variants Count", x = "", y = "Variant Count") +
  scale_fill_manual(values = c("Control Only" = "#5E81AC",   # Soft blue
                               "Survived Only" = "#BF616A",  # Muted red
                               "Common" = "#A3BE8C")) +      # Subtle green
  theme(legend.position = "none")  # Hide legend (optional)

# Save plot as PNG
ggsave("SNP_counts_barplot_new.png", plot = SNP_counts_barplot, width = 8, height = 6, dpi = 300)


# Display success message
message("Plot saved successfully as: SNP_counts_barplot_new.png")

# Display the plot in RStudio
print(SNP_counts_barplot_new)

## For Step2 (to be followed) the packages: dplyr, tidyr, ggplot2 already installed
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

## Step2: Unzip and Load the GFF File
# Path to the GFF file (Compressed)
gff_compressed <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/GCF_015228065.2_NSTDA_Pmon_1_genomic.gff.gz"

# Unzip the GFF file if not already unzipped
gff_unzipped <- sub(".gz$", "", gff_compressed)  # Remove .gz extension
if (!file.exists(gff_unzipped)) {
  gunzip(gff_compressed, overwrite = FALSE)
  message("GFF file unzipped successfully.")
}

# Load GFF file
gff_data <- import(gff_unzipped)

# Convert to DataFrame
gff_df <- as.data.frame(gff_data)

## Step3: Extract Relevant Gene Information

# Keep only gene-related data
genes_data <- gff_df %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, ID)  # Selecting Chromosome, Start, End, and Gene ID

# Rename columns for clarity
colnames(genes_data) <- c("CHROM_Gff", "Gene_Start", "Gene_End", "Gene_ID")

# Display extracted gene data
head(genes_data)

## Step4: Merge unique SNPs data with gff annotations
# Merge SNPs from Control Group with Gene Data from gff file (keeping ANN)
merged_control_new <- control_only_new %>%
  rowwise() %>%
  mutate(Gene_ID = list(genes_data$Gene_ID[which(genes_data$Gene_Start <= POS & 
                                                   genes_data$Gene_End >= POS & 
                                                   genes_data$CHROM_Gff == CHROM)])) %>%
  unnest(cols = c(Gene_ID)) %>%
  filter(!is.na(Gene_ID)) %>%
  select(CHROM, POS, REF, ALT, ANN, Gene_ID)  # Keep ANN

# Merge SNPs from Survived Group with Gene Data (keeping ANN)
merged_survived_new <- survived_only_new %>%
  rowwise() %>%
  mutate(Gene_ID = list(genes_data$Gene_ID[which(genes_data$Gene_Start <= POS & 
                                                   genes_data$Gene_End >= POS & 
                                                   genes_data$CHROM_Gff == CHROM)])) %>%
  unnest(cols = c(Gene_ID)) %>%
  filter(!is.na(Gene_ID)) %>%
  select(CHROM, POS, REF, ALT, ANN, Gene_ID)  # Keep ANN

# save merged control and survived files in csv format
write.csv(merged_control_new, "Unique_Control_SNPs.csv", row.names = FALSE)
write.csv(merged_survived_new, "Unique_Survived_SNPs.csv", row.names = FALSE)
message("SNP tables saved as CSV files.")

# Step5: Create SNP Distribution Table
# SNP Distribution for Control Group (counting SNPs per gene)
snp_distribution_control_new <- merged_control_new %>%
  group_by(Gene_ID) %>%
  summarise(SNP_count = n(),
            Variant_Effects = paste(unique(ANN), collapse = ";")) %>%  # Preserve ANN information
  mutate(SNP_bin = cut(SNP_count,
                       breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, Inf),
                       labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                  "11-15", "16-20", "21-25", "26-30", ">30"),
                       right = FALSE)) %>%
  group_by(SNP_bin) %>%
  summarise(Number_of_genes = n()) %>%
  mutate(Genotype = "Control")

# SNP Distribution for Survived Group (counting SNPs per gene)
snp_distribution_survived_new <- merged_survived_new %>%
  group_by(Gene_ID) %>%
  summarise(SNP_count = n(),
            Variant_Effects = paste(unique(ANN), collapse = ";")) %>%  # Preserve ANN information
  mutate(SNP_bin = cut(SNP_count,
                       breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, Inf),
                       labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                  "11-15", "16-20", "21-25", "26-30", ">30"),
                       right = FALSE)) %>%
  group_by(SNP_bin) %>%
  summarise(Number_of_genes = n()) %>%
  mutate(Genotype = "Survived")

# Combine Both Groups for Final Plot
snp_distribution_combined_new <- bind_rows(snp_distribution_control_new, snp_distribution_survived_new)

## Step6: Plot SNP Distribution Per Gene
# Create the bar plot
snp_distribution_plot <- ggplot(snp_distribution_combined_new, aes(x = SNP_bin, y = Number_of_genes, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Border for clarity
  labs(title = "SNP Distribution per Gene (Control vs Survived)",
       x = "Number of SNPs per Gene",
       y = "Number of Genes") +
  theme_minimal(base_size = 15) +  # Clean publication theme
  scale_fill_manual(values = c("Control" = "#5E81AC", "Survived" = "#BF616A")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels

# Save plot as PNG
ggsave("SNPs_per_Gene_distribution.png", plot = snp_distribution_plot, width = 8, height = 6, dpi = 300)

# Display success message
message("Plot saved successfully: SNPs_per_Gene_distributionp.png")

# Display the plot in RStudio
print(snp_distribution_plot)

## Step7: Impact plots
# 7.1: For Control
# Count High, Moderate, Low impact variants
impact_counts <- merged_control_new %>%
  mutate(Impact = case_when(
    grepl("HIGH", ANN) ~ "High Impact",
    grepl("MODERATE", ANN) ~ "Moderate Impact",
    grepl("LOW", ANN) ~ "Low Impact",
    TRUE ~ "Modifier"
  )) %>%
  group_by(Impact) %>%
  summarise(Count = n())

# Plot impact distribution
variant_imapct_distribution_plot <- ggplot(impact_counts, aes(x = Impact, y = Count, fill = Impact)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Variant Impact Distribution", x = "Impact Level", y = "Number of Variants")

# Save plot as PNG
ggsave("Control_unique_SNPs__impact_distribution.png", plot = variant_imapct_distribution_plot, width = 8, height = 6, dpi = 300)

# Display success message
message("Plot saved successfully: Scontrol_NPs_pimpact _istributionp.png")

# Display the plot in RStudio
print(variant_imapct_distribution_plot)
