# Calculate 97% representative species from MetaPhlAn4 data
# This script processes the metaphlan4.Rds file to identify the top 97% 
# most abundant taxa and exports them to CSV

# Load required libraries
library(phyloseq)

# Load the MetaPhlAn4 data
cat("Loading MetaPhlAn4 data...\n")
data <- readRDS('data/metaphlan4.Rds')

# Get OTU table as matrix
otu_data <- as.matrix(otu_table(data))

# Calculate total abundance for each taxon across all samples
total_abundance <- rowSums(otu_data)

# Calculate relative abundance (percentage of total)
total_sum <- sum(total_abundance)
relative_abundance <- (total_abundance / total_sum) * 100

# Create a data frame with taxon names, abundance, and relative abundance
abundance_df <- data.frame(
  taxon = names(total_abundance),
  abundance = total_abundance,
  relative_abundance = relative_abundance,
  stringsAsFactors = FALSE
)

# Sort by relative abundance in descending order
abundance_df <- abundance_df[order(abundance_df$relative_abundance, decreasing = TRUE), ]

# Calculate cumulative relative abundance
abundance_df$cumulative_rel_abundance <- cumsum(abundance_df$relative_abundance)

# Find taxa that make up the top 97% cumulative relative abundance
top_97_percent <- abundance_df[abundance_df$cumulative_rel_abundance <= 97, ]

# Add the next taxon that crosses the 97% threshold (if any)
if(nrow(top_97_percent) < nrow(abundance_df)) {
  next_row_idx <- nrow(top_97_percent) + 1
  if(next_row_idx <= nrow(abundance_df)) {
    top_97_percent <- rbind(top_97_percent, abundance_df[next_row_idx, ])
  }
}

# Remove the cumulative column for final output
final_df <- top_97_percent[, c('taxon', 'abundance', 'relative_abundance')]

# Print summary information
cat('=== Summary ===\n')
cat('Total number of taxa in dataset:', nrow(abundance_df), '\n')
cat('Number of taxa representing 97% abundance:', nrow(final_df), '\n')
cat('Cumulative relative abundance of selected taxa:', round(sum(final_df$relative_abundance), 2), '%\n')
cat('Top 10 most abundant taxa:\n')
print(head(final_df, 10))

# Save to CSV
write.csv(final_df, 'data/97_representative_species.csv', row.names = FALSE)
cat('\nCSV file saved as: data/97_representative_species.csv\n')