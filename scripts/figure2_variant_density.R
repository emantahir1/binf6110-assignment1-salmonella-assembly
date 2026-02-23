# Figure 2: Variant Density Comparison
# Bar chart showing difference in variant density

# Create output directory if needed
dir.create("figures", showWarnings = FALSE)

# Data
contigs <- c("Chromosome", "Plasmid")
density <- c(0.9, 78.52)
colors <- c("#4A90E2", "#E74C3C")  # Blue for chromosome, red for plasmid

# Create PNG
png("figures/figure2_variant_density.png", width = 2800, height = 2000, res = 300)

# Set margins
par(mar = c(5, 5, 3, 2))

# Create bar plot
bp <- barplot(density,
              names.arg = contigs,
              col = colors,
              border = NA,
              ylim = c(0, 85),
              ylab = "Variant Density (variants/kb)",
              main = "Variant Density by Contig",
              cex.names = 1.3,
              cex.axis = 1.2,
              cex.lab = 1.3,
              cex.main = 1.5,
              las = 1)

# Add value labels on bars
text(bp, density + 3, 
     labels = c("0.9", "78.52"),
     cex = 1.2,
     font = 2)

# Add annotation for fold-difference
text(bp[2], 45, 
     labels = "87-fold higher\nthan chromosome",
     cex = 1.1,
     font = 2)

# Add grid for readability
abline(h = seq(0, 80, 20), col = "gray90", lty = 2)

dev.off()

cat("Figure 2 created: figures/figure2_variant_density.png\n")
