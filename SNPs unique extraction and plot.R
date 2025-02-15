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


# Extract unique variants based on CHROM, POS, REF, and ALT
control_variants <- control_data %>%
  distinct(CHROM, POS, REF, ALT)  # Unique SNPs in Control

survived_variants <- survived_data %>%
  distinct(CHROM, POS, REF, ALT)  # Unique SNPs in Survived

# Identify unique and common variants
control_only <- anti_join(control_variants, survived_variants, by = c("CHROM", "POS", "REF", "ALT"))
survived_only <- anti_join(survived_variants, control_variants, by = c("CHROM", "POS", "REF", "ALT"))
common_variants <- inner_join(control_variants, survived_variants, by = c("CHROM", "POS", "REF", "ALT"))

# Count variants in each category
variant_counts <- data.frame(
  Category = c("Control Only", "Survived Only", "Common"),
  Count = c(nrow(control_only), nrow(survived_only), nrow(common_variants))
)

# Plot with labels and improved colors
SNP_counts_barplot <- ggplot(variant_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  # Subtle border for clarity
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +      # Add count labels
  theme_minimal(base_size = 15) +                              # Adjust font size
  labs(title = "Unique and Common Variants Count", x = "", y = "Variant Count") +
  scale_fill_manual(values = c("Control Only" = "#5E81AC",   # Soft blue
                               "Survived Only" = "#BF616A",  # Muted red
                               "Common" = "#A3BE8C")) +      # Subtle green
  theme(legend.position = "none")  # Hide legend (optional)

# Save plot as PNG
ggsave("SNP_counts_barplot.png", plot = SNP_counts_barplot, width = 8, height = 6, dpi = 300)


# Display success message
message("Plot saved successfully as: SNP_counts_barplot.png")

# Display the plot in RStudio
print(SNP_counts_barplot)
