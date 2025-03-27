# Load required libraries
library(dplyr)
library(readr)

# Define the list of files
files <- c("C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/RNASeq_Shrimp/unique_control_indels_insertions.csv.gz",
           "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/RNASeq_Shrimp/unique_control_indels_deletions.csv.gz",
           "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/RNASeq_Shrimp/unique_control_indels_mixed.csv.gz",
           "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/RNASeq_Shrimp/unique_survived_indels_insertions.csv.gz",
           "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/RNASeq_Shrimp/unique_survived_indels_deletions.csv.gz",
           "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/RNASeq_Shrimp/unique_survived_indels_mixed.csv.gz")

# Function to read, process, and merge files
process_and_merge_indels <- function(file_list) {
  processed_files <- lapply(file_list, function(file) {
    # Extract group and indel type from filename
    parts <- strsplit(basename(file), "_")[[1]]
    group <- parts[2]      # "control" or "survived"
    indel_type <- gsub(".csv.gz", "", parts[4]) # "insertions", "deletions", or "mixed"
    
    # Read the file
    data <- read_csv(file) 
    
    # Add columns
    data <- data %>%
      mutate(Group = group, Indel_Type = indel_type)
    
    return(data)
  })
  
  # Combine all processed files
  data_combined <- bind_rows(processed_files)
  
  # Split into control and survived datasets
  control_data <- filter(data_combined, Group == "control")
  survived_data <- filter(data_combined, Group == "survived")
  
  # Save the merged datasets
  write_csv(control_data, "Merged_Control_Indels.csv")
  write_csv(survived_data, "Merged_Survived_Indels.csv")
  
  message("Merged files saved successfully.")
}

# Run the function with your file list
process_and_merge_indels(files)

**
# Load required libraries
library(ggplot2)
library(tidyr)

# Step 1: Create a copy of CSV files as CSV.GFF
file.copy("Merged_Control_Indels.csv", "Merged_Control_Indels.csv.gz")
file.copy("Merged_Survived_Indels.csv", "Merged_Survived_Indels.csv.gz")

# Step 2: Read the GFF file (assuming tab-delimited format)
gff_file <- "C:/Users/hp/Downloads/Shrimp/RNASeq_Shrimp/GCF_015228065.2_NSTDA_Pmon_1_genomic.gff"  # Update with your actual GFF file path

# Assuming GFF format: Chromosome (col 1), Start (col 4), End (col 5), Gene_ID (from attributes in col 9)
gff_data <- read.delim(gff_file, header = FALSE, comment.char = "#", sep = "\t")

# Extract relevant columns and Gene_ID from attributes
gff_data <- gff_data %>%
  filter(V3 == "gene") %>%  # Keep only gene entries
  mutate(Gene_ID = sub(".*ID=([^;]+).*", "\\1", V9)) %>%
  select(CHROM_Gff = V1, Gene_Start = V4, Gene_End = V5, Gene_ID)

# ---- Step 2: Read variant data ----
control_data <- read_csv("Merged_Control_Indels.csv")
survived_data <- read_csv("Merged_Survived_Indels.csv")

# ---- Step 3: Function to merge variants with GFF annotations ----
merge_with_gff <- function(variant_data, genes_data) {
  merged_data <- variant_data %>%
    left_join(genes_data, by = c("CHROM" = "CHROM_Gff")) %>%
    filter(Gene_Start <= POS & Gene_End >= POS) %>%  # Keep only variants within gene regions
    select(CHROM, POS, REF, ALT, ANN, Group, Indel_Type, Gene_ID) %>%
    distinct()  # Remove duplicates if any

  return(merged_data)
}

# Apply function to datasets
gff_merged_control_indels <- merge_with_gff(control_data, gff_data)
gff_merged_survived_indels <- merge_with_gff(survived_data, gff_data)
# Save the merged datasets
  write_csv(gff_merged_control_indels, "gff_merged_control_indels.csv")
  write_csv(gff_merged_survived_indels, "gff_merged_survived_indels.csv")

# ---- Step 4: Function to calculate Indel distribution per gene ----
calculate_distribution <- function(merged_data, genotype) {
  indel_distribution <- merged_data %>%
    group_by(Gene_ID, Indel_Type) %>%
    summarise(Indel_count = n(),
              Variant_Effects = paste(unique(ANN), collapse = ";"),
              .groups = "drop") %>%
    mutate(Indel_bin = factor(cut(Indel_count,
                                  breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, Inf),
                                  labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                             "11-15", "16-20", "21-25", "26-30", ">30"),
                                  right = FALSE),
                             levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                        "11-15", "16-20", "21-25", "26-30", ">30"))) %>%
    group_by(Indel_bin, Indel_Type) %>%
    summarise(Number_of_genes = n(), .groups = "drop") %>%
    mutate(Genotype = genotype)
  
  return(indel_distribution)
}

# ---- Step 5: Calculate Indel Distribution ----
indel_distribution_control <- calculate_distribution(gff_merged_control_indels, "Control")
indel_distribution_survived <- calculate_distribution(gff_merged_survived_indels, "Survived")

# revised part
#third
library(ggplot2)
library(dplyr)

# ---- Step 6: Define Colors ----
color_palette <- c("Insertion_Control" = "#1F78B4",   # Dark Blue (Control)
                   "Insertion_Survived" = "#A6CEE3",  # Light Blue (Survived)
                   "Deletion_Control" = "#E31A1C",   # Dark Red (Control)
                   "Deletion_Survived" = "#FB9A99",  # Light Red (Survived)
                   "Mixed_Control" = "#33A02C",      # Dark Green (Control)
                   "Mixed_Survived" = "#B2DF8A")     # Light Green (Survived)

# ---- Step 7: Modify Data for Custom Fill Mapping ----
indel_distribution_combined <- indel_distribution_combined %>%
  mutate(Fill_Group = paste(Indel_Type, Genotype, sep = "_"))  # Create unique fill categories

# ---- Step 8: Plot Stacked Bar Chart with Colors & Counts ----
indel_distribution_plot <- ggplot(indel_distribution_combined, 
                                  aes(x = Indel_bin, y = Number_of_genes, fill = Fill_Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Stacked bars with border
  geom_text(aes(label = Number_of_genes), 
            position = position_stack(vjust = 0.5),  # Text placed at center of each stack
            size = 5, color = "black") +  
  scale_fill_manual(values = color_palette) +  # Apply color mapping
  facet_wrap(~ Genotype, ncol = 1) +  # Separate Control and Survived
  labs(title = "Indel Distribution per Gene (Control vs Survived)",
       x = "Number of Indels per Gene",
       y = "Number of Genes") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save plot
ggsave("Indel_Variant_Distribution_Per_Gene.png", plot = indel_distribution_plot, width = 12, height = 6, dpi = 300)

# Show the plot
print(indel_distribution_plot)