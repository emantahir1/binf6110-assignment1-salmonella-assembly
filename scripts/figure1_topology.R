# Figure 1: Assembly Topology
# Simple, clean pie chart showing contig topology

# Create output directory if needed
dir.create("figures", showWarnings = FALSE)

# Data
topology <- c(2, 1)  # 2 linear, 1 circular
labels <- c("Linear contigs\n(Chromosome fragments)", 
            "Circular contig\n(Plasmid)")
colors <- c("#4A90E2", "#50C878")  # Professional blue and green

# Create PNG
png("figures/figure1_topology.png", width = 2400, height = 2000, res = 300)

# Set margins
par(mar = c(2, 2, 3, 2))

# Create pie chart
pie(topology, 
    labels = labels,
    col = colors,
    main = "Assembly Topology",
    cex.main = 1.5,
    cex = 1.2,
    border = "white",
    lwd = 2)

# Add legend
legend("bottomright", 
       legend = c("Linear (n=2)", "Circular (n=1)"),
       fill = colors,
       border = "white",
       cex = 1.1,
       bty = "n")

dev.off()

cat("Figure 1 created: figures/figure1_topology.png\n")
