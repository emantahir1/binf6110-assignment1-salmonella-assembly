library(tidyverse)
library(ggplot2)
library(patchwork)

# Read coverage data
coverage <- read.table("detailed_coverage.txt", header=TRUE, comment.char="")
colnames(coverage)[1] <- "contig"

# Create cleaner labels
coverage <- coverage %>%
  mutate(label = ifelse(grepl("003197", contig), "Chromosome", "Plasmid"))

# Plot 1: Reference Coverage %
p1 <- ggplot(coverage, aes(x=label, y=coverage)) +
  geom_col(aes(fill=label), width=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("Chromosome"="#3498DB", "Plasmid"="#E74C3C")) +
  geom_hline(yintercept=95, linetype="dashed", color="darkgreen", size=1) +
  geom_text(aes(label=paste0(round(coverage, 1), "%")), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="A. Reference Genome Coverage",
       x="Contig",
       y="Coverage (%)") +
  ylim(0, 105) +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold")) +
  annotate("text", x=1.5, y=98, label="95% threshold", 
           color="darkgreen", size=3.5, fontface="italic")

# Plot 2: Sequencing Depth
p2 <- ggplot(coverage, aes(x=label, y=meandepth)) +
  geom_col(aes(fill=label), width=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("Chromosome"="#3498DB", "Plasmid"="#E74C3C")) +
  geom_text(aes(label=paste0(round(meandepth), "×")), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="B. Mean Sequencing Depth",
       x="Contig",
       y="Depth (×)") +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold"))

# Plot 3: Mapping Quality
p3 <- ggplot(coverage, aes(x=label, y=meanmapq)) +
  geom_col(aes(fill=label), width=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("Chromosome"="#3498DB", "Plasmid"="#E74C3C")) +
  geom_hline(yintercept=60, linetype="dashed", color="darkgreen", size=1) +
  geom_text(aes(label=round(meanmapq)), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="C. Mean Mapping Quality (MAPQ)",
       x="Contig",
       y="MAPQ Score") +
  ylim(0, 65) +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold")) +
  annotate("text", x=1.5, y=62, label="High quality (60)", 
           color="darkgreen", size=3.5, fontface="italic")

# Plot 4: Number of aligned reads
p4 <- ggplot(coverage, aes(x=label, y=numreads)) +
  geom_col(aes(fill=label), width=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("Chromosome"="#3498DB", "Plasmid"="#E74C3C")) +
  geom_text(aes(label=format(numreads, big.mark=",")), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="D. Number of Aligned Reads",
       x="Contig",
       y="Read Count") +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold"))

# Combine
final <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title="Figure 2: Reference Alignment Quality Assessment",
    subtitle="Raw Reads Aligned to Salmonella enterica Reference Genome",
    theme=theme(plot.title=element_text(size=17, face="bold", hjust=0.5),
                plot.subtitle=element_text(size=13, hjust=0.5, face="italic"))
  )

ggsave("figure2_alignment_final.png", final, width=14, height=11, dpi=300)

cat("\n✅ SUCCESS!\n")
cat("Saved: figure2_alignment_final.png\n")
cat("\nShows:\n")
cat("- 97.8% chromosome coverage vs 45.7% plasmid\n")
cat("- 151× chromosome depth vs 98× plasmid\n")
cat("- High mapping quality (MAPQ 60 and 52)\n")
cat("- 183,082 reads aligned to chromosome\n")
