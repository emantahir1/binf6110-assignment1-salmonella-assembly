library(tidyverse)
library(ggplot2)
library(patchwork)

# Read variant data
variants <- read.table("variant_data.txt", 
                       col.names=c("contig", "pos", "ref_len", "alt_len", "qual"))

# Classify variants
variants <- variants %>%
  mutate(type = case_when(
    ref_len == 1 & alt_len == 1 ~ "SNP",
    ref_len < alt_len ~ "Insertion",
    ref_len > alt_len ~ "Deletion",
    TRUE ~ "Other"
  )) %>%
  mutate(contig_label = ifelse(grepl("003197", contig), "Chromosome", "Plasmid"))

# Calculate variant density
variant_summary <- variants %>%
  group_by(contig_label) %>%
  summarise(
    total = n(),
    length_mb = ifelse(first(contig_label) == "Chromosome", 4.86, 0.09),
    density = n() / length_mb / 1000,
    .groups = "drop"
  )

# Variant counts by type and contig
type_by_contig <- variants %>%
  count(contig_label, type)

# Overall type distribution
type_summary <- variants %>%
  count(type) %>%
  mutate(percentage = n / sum(n) * 100)

# Plot 1: Variant counts by contig
p1 <- ggplot(variant_summary, aes(x=contig_label, y=total)) +
  geom_col(aes(fill=contig_label), width=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("Chromosome"="#3498DB", "Plasmid"="#E74C3C")) +
  geom_text(aes(label=format(total, big.mark=",")), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="A. Total Variants by Contig",
       x="Contig",
       y="Number of Variants") +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold"))

# Plot 2: Variant density (variants per kb)
p2 <- ggplot(variant_summary, aes(x=contig_label, y=density)) +
  geom_col(aes(fill=contig_label), width=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("Chromosome"="#3498DB", "Plasmid"="#E74C3C")) +
  geom_text(aes(label=round(density, 2)), 
            vjust=-0.5, size=5, fontface="bold") +
  theme_minimal(base_size=14) +
  labs(title="B. Variant Density (variants/kb)",
       x="Contig",
       y="Variants per kb") +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold")) +
  annotate("text", x=1.5, y=max(variant_summary$density)*0.5, 
           label=sprintf("%.0f× higher on plasmid", 
                        max(variant_summary$density)/min(variant_summary$density)),
           size=4, fontface="italic", color="red")

# Plot 3: Variant type composition (pie chart)
p3 <- ggplot(type_summary, aes(x="", y=percentage, fill=type)) +
  geom_col(width=1, color="white", size=1.5) +
  coord_polar(theta="y") +
  scale_fill_manual(values=c("SNP"="#2E86AB", "Insertion"="#A23B72", 
                              "Deletion"="#F18F01", "Other"="#95A5A6")) +
  theme_void(base_size=14) +
  labs(title="C. Genome-Wide Variant Composition",
       fill="Type") +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=14),
        legend.text=element_text(size=11)) +
  geom_text(aes(label=paste0(round(percentage, 1), "%")), 
            position=position_stack(vjust=0.5), 
            size=4.5, fontface="bold", color="white")

# Plot 4: Variant types by contig (stacked bar)
p4 <- ggplot(type_by_contig, aes(x=contig_label, y=n, fill=type)) +
  geom_col(position="stack", width=0.6) +
  scale_fill_manual(values=c("SNP"="#2E86AB", "Insertion"="#A23B72", 
                              "Deletion"="#F18F01", "Other"="#95A5A6")) +
  theme_minimal(base_size=14) +
  labs(title="D. Variant Type Distribution by Contig",
       x="Contig",
       y="Number of Variants",
       fill="Variant Type") +
  theme(plot.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=12, face="bold"),
        legend.position="right")

# Combine
final <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title="Figure 3: Variant Landscape Analysis",
    subtitle="Genome-Wide Variant Distribution and Composition",
    theme=theme(plot.title=element_text(size=17, face="bold", hjust=0.5),
                plot.subtitle=element_text(size=13, hjust=0.5, face="italic"))
  )

ggsave("figure3_variants_final.png", final, width=14, height=11, dpi=300)

cat("\n✅ SUCCESS!\n")
cat("Saved: figure3_variants_final.png\n")
cat("\nKey Findings:\n")
cat("- Total variants: ", sum(variant_summary$total), "\n", sep="")
cat("- Chromosome: ", variant_summary$
