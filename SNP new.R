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

# Plot with labels and improved colors
SNP_counts_barplot_new <- ggplot(variant_counts_new, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  # Subtle border for clarity
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +      # Add count labels
  theme_minimal(base_size = 15) +                              # Adjust font size
  labs(title = "Unique and Common Variants Count", x = "", y = "Variant Count") +
  scale_fill_manual(values = c("Control Only" = "#5E81AC",   # Soft blue
                               "Survived Only" = "#BF616A",  # Muted red
                               "Common" = "#A3BE8C")) +      # Subtle green
  theme(legend.position = "none")  # Hide legend (optional)