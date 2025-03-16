

## Step1: Extract Unique Insertions and Deletions (Indels) and plot them

# Define required packages
packages <- c("vcfR", "dplyr", "tidyr", "ggplot2")

# Install any missing packages
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)

# Load packages
lapply(packages, library, character.only = TRUE)

# Function to read and extract data from VCF files
read_vcf_data_indels <- function(file_path) {
  vcf_indels <- read.vcfR(file_path)
  data_indels <- data.frame(
    CHROM = vcf_indels@fix[,1],  # Chromosome/Scaffold
    POS = vcf_indels@fix[,2],    # Position
    REF = vcf_indels@fix[,4],    # Reference Allele
    ALT = vcf_indels@fix[,5],    # Alternate Allele
    ANN = vcf_indels@fix[,8]     # Annotation
  )
  return(data_indels)
}

# File paths
control_vcf_path_indels <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/NEW_DATA-20250215T055250Z-001/NEW_DATA/CONTROL_36_INDEL_annotation/36_INDEL.ann.vcf"
survived_vcf_path_indels <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/NEW_DATA-20250215T055250Z-001/NEW_DATA/SURVIVED_35_INDEL_annotation/35_INDEL.ann.vcf"

# Read and extract data
control_data_indels <- read_vcf_data_indels(control_vcf_path_indels)
survived_data_indels <- read_vcf_data_indels(survived_vcf_path_indels)

# Create copies for each group
control_indels_copy <- control_data_indels
survived_indels_copy <- survived_data_indels


# Function to classify variants
classify_variant <- function(ref, alt) {
  ref_alleles <- strsplit(ref, ",")[[1]]
  alt_alleles <- strsplit(alt, ",")[[1]]
  
  insertion <- any(nchar(alt_alleles) > nchar(ref_alleles))
  deletion <- any(nchar(ref_alleles) > nchar(alt_alleles))
  
  if (insertion & deletion) {
    return("Mixed")
  } else if (insertion) {
    return("Insertion")
  } else if (deletion) {
    return("Deletion")
  } else {
    return("Other")
  }
}

# Apply classification
control_indels_copy$Type <- mapply(classify_variant, control_indels_copy$REF, control_indels_copy$ALT)
survived_indels_copy$Type <- mapply(classify_variant, survived_indels_copy$REF, survived_indels_copy$ALT)

# Save classified copies
write.csv(control_indels_copy, gzfile("control_indels_copy.csv.gz"), row.names = FALSE)
write.csv(survived_indels_copy, gzfile("survived_indels_copy.csv.gz"), row.names = FALSE)

# Identify and Count multi-allelic variants
control_multiallelic <- control_indels_copy %>% filter(grepl(",", REF) | grepl(",", ALT))
survived_multiallelic <- survived_indels_copy %>% filter(grepl(",", REF) | grepl(",", ALT))

# Check if they are grouped as one or separated
head(control_multiallelic)
head(survived_multiallelic)

control_multi_allelic_count <- nrow(control_multiallelic)
survived_multi_allelic_count <- nrow(survived_multiallelic)
message("Multi-allelic variants in control: ", control_multi_allelic_count)
message("Multi-allelic variants in survived: ", survived_multi_allelic_count)

# Save multi-allelic variants
write.csv(control_multiallelic, gzfile("control_multiallelic.csv.gz"), row.names = FALSE)
write.csv(survived_multiallelic, gzfile("survived_multiallelic.csv.gz"), row.names = FALSE)

# Count types
control_each_indeltype_counts <- table(control_multiallelic$Type)
survived_each_indeltype_counts <- table(survived_multiallelic$Type)

message("Control Multi-allelic Counts (Insertions, Deletions, Mixed): ", toString(control_each_indeltype_counts))
message("Survived Multi-allelic Counts (Insertions, Deletions, Mixed): ", toString(survived_each_indeltype_counts))

# Filter classified variants
control_insertions <- filter(control_indels_copy, Type == "Insertion")
control_deletions <- filter(control_indels_copy, Type == "Deletion")
control_mixed <- filter(control_indels_copy, Type == "Mixed")

survived_insertions <- filter(survived_indels_copy, Type == "Insertion")
survived_deletions <- filter(survived_indels_copy, Type == "Deletion")
survived_mixed <- filter(survived_indels_copy, Type == "Mixed")

# Function to process variants
generate_variant_counts <- function(control, survived, title, filename) {
  control_variants <- control %>%
    group_by(CHROM, POS, REF, ALT) %>%
    summarise(ANN = paste(unique(ANN), collapse = ";"), .groups = "drop")
  
  survived_variants <- survived %>%
    group_by(CHROM, POS, REF, ALT) %>%
    summarise(ANN = paste(unique(ANN), collapse = ";"), .groups = "drop")
  
  control_only <- anti_join(control_variants, survived_variants, by = c("CHROM", "POS", "REF", "ALT")) 
  survived_only <- anti_join(survived_variants, control_variants, by = c("CHROM", "POS", "REF", "ALT")) 
  common_variants <- inner_join(control_variants, survived_variants, by = c("CHROM", "POS", "REF", "ALT"), suffix = c("_control", "_survived"))
  
  variant_counts <- data.frame(
    Category = c("Control Only", "Survived Only", "Common"),
    Count = c(nrow(control_only), nrow(survived_only), nrow(common_variants))
  )
  
  plot <- ggplot(variant_counts, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(aes(label = Count), vjust = -0.5, size = 5) +
    theme_minimal(base_size = 15) +
    labs(title = title, x = "", y = "Variant Count") +
    scale_fill_manual(values = c("Control Only" = "#5E81AC", "Survived Only" = "#BF616A", "Common" = "#A3BE8C")) +
    theme(legend.position = "none")
  
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300)
  message("Plot saved successfully as: ", filename)
  print(plot)
}

# Generate plots for Insertions, Deletions, and Mixed Variants
insertions_count_plot <-generate_variant_counts(control_insertions, survived_insertions, "Unique and Common Insertions Count", "Insertions_counts_barplot.png")
deletions_count_plot <- generate_variant_counts(control_deletions, survived_deletions, "Unique and Common Deletions Count", "Deletions_counts_barplot.png")
mixed_count_plot <- generate_variant_counts(control_mixed, survived_mixed, "Unique and Common Mixed Variants Count", "Mixed_counts_barplot.png")

# Save plot as PNG
ggsave("Insertions_counts_barplot.png", plot = insertions_count_plot, width = 8, height = 6, dpi = 300)
ggsave("Deletions_counts_barplot.png", plot = deletions_count_plot, width = 8, height = 6, dpi = 300)
ggsave("Mixed_counts_barplot.png", plot = mixed_count_plot, width = 8, height = 6, dpi = 300)

