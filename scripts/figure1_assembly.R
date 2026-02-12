library(tidyverse)
library(ggplot2)
library(patchwork)

# Read assembly info from Flye
assembly <- read.table("flye_output/assembly_info.txt", 
                       header=TRUE, sep="\t", comment.char="")
colnames(assembly)[1] <- "contig"

# Plot 1: Contig lengths
p1 <- ggplot(assembly, aes(x=reorder(contig, -length), y=length/1e6)) +
  geom_col(fill="#2E86AB", width=0.6) +
  geom_text(aes(label=sprintf("%.2f Mb", length/1e6)), 
            vjust=-0.5, size=4.5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="A. Assembled Contig Lengths",
       x="Contig",
       y="Length (Mb)") +
  theme(axis.text.x=element_text(size=12, face="bold"),
        plot.title=element_text(face="bold", size=14)) +
  ylim(0, max(assembly$length/1e6)*1.15)

# Plot 2: Coverage per contig
p2 <- ggplot(assembly, aes(x=reorder(contig, -length), y=cov.)) +
  geom_col(fill="#A23B72", width=0.6) +
  geom_text(aes(label=paste0(round(cov.), "×")), 
            vjust=-0.5, size=4.5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="B. Sequencing Coverage per Contig",
       x="Contig",
       y="Coverage (×)") +
  theme(axis.text.x=element_text(size=12, face="bold"),
        plot.title=element_text(face="bold", size=14)) +
  ylim(0, max(assembly$cov.)*1.15)

# Plot 3: Circularity status
circ_data <- assembly %>%
  mutate(status = ifelse(circ. == "Y", "Circular", "Linear")) %>%
  count(status)

p3 <- ggplot(circ_data, aes(x="", y=n, fill=status)) +
  geom_col(width=1, color="white", size=2) +
  coord_polar(theta="y") +
  scale_fill_manual(values=c("Circular"="#27AE60", "Linear"="#3498DB")) +
  theme_void(base_size=14) +
  labs(title="C. Contig Topology",
       fill="Status") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14),
        legend.text=element_text(size=12, face="bold")) +
  geom_text(aes(label=paste0(n, "\n", status)), 
            position=position_stack(vjust=0.5), 
            size=5, fontface="bold", color="white")

# Plot 4: Assembly statistics as text plot
stats_text <- paste0(
  "Total Assembly Length: 5.10 Mb\n",
  "Number of Contigs: 3\n",
  "Largest Contig: 3.32 Mb\n",
  "Contig N50: 3.32 Mb\n",
  "Mean Coverage: 160×\n",
  "Circular Contigs: 1 (plasmid)"
)

p4 <- ggplot() +
  annotate("text", x=0.5, y=0.5, label=stats_text, 
           size=5, hjust=0, vjust=0.5, fontface="bold", lineheight=1.3) +
  theme_void() +
  labs(title="D. Assembly Summary") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14)) +
  xlim(0, 1) + ylim(0, 1)

# Combine
final <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title="Figure 1: Genome Assembly Quality Assessment",
    subtitle="Salmonella enterica - Oxford Nanopore Long-Read Sequencing",
    theme=theme(plot.title=element_text(size=17, face="bold", hjust=0.5),
                plot.subtitle=element_text(size=13, hjust=0.5, face="italic"))
  )

ggsave("figure1_assembly_final.png", final, width=14, height=11, dpi=300)

cat("\n✅ SUCCESS!\n")
cat("Saved: figure1_assembly_final.png\n")
cat("\nShows:\n")
cat("- Contig lengths (3 contigs: 3.32 Mb, 1.68 Mb, 0.11 Mb)\n")
cat("- Coverage per contig (153×, 169×, 245×)\n")
cat("- Topology (1 circular, 2 linear)\n")
cat("- Summary statistics\n")
