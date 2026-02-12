# Assignment 1: Genome Assembly of *Salmonella enterica*

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Methods](#2-methods)
3. [Results](#3-results)
4. [Discussion](#4-discussion)
5. [Conclusions](#5-conclusions)
6. [References](#6-references)

---

## 1. Introduction

### Biological Background and Importance

*Salmonella enterica* is a Gram-negative, facultative intracellular bacterium and a major cause of foodborne illness worldwide (Knodler & Elfenbein, 2019). It includes many different serovars. Some, such as *S. Typhi*, are restricted to specific hosts, while others, like *S. Typhimurium*, infect a broad range of hosts and commonly cause gastroenteritis in both humans and animals (Johnson et al., 2017). These groups differ in how they spread, which hosts they infect, and how severe the resulting disease can be. In practice, this means that *Salmonella enterica* is not a single, uniform pathogen, but a genetically diverse species with meaningful biological variation (Andino & Hanning, 2015). Because of this diversity, high-quality genome assemblies are especially important for *Salmonella enterica*. They allow closely related strains to be distinguished, make outbreak tracking possible, and support the identification of genes linked to virulence and antimicrobial resistance (Dele Ogunremi et al., 2014). Together, these applications form the basis of modern public health surveillance.

### Challenges in Genome Assembly and Project Goals

Genome assembly is how a full genome is put back together from the many short DNA fragments produced by sequencing (Basantani et al., 2017). This process turns raw reads into a continuous sequence that represents an organism’s genome, making it possible to examine genes, genome structure, and variation (Basantani et al., 2017).

Genome assembly is computationally challenging because sequencing reads contain errors and genomes often include repetitive regions that are difficult to place correctly (Schiffer et al., 2025). When identical or near-identical sequences occur in multiple locations, assembly algorithms may be unable to determine their true origin, leading to fragmented or misassembled genomes (Schiffer et al., 2025). These challenges are especially pronounced when using short-read data. Short-read platforms such as Illumina generate highly accurate reads, but their limited length means they often cannot span repetitive elements, resulting in assemblies broken into many small contigs (Hu et al., 2021).

Long-read sequencing offers a different approach. Technologies such as Oxford Nanopore produce reads that span several kilobases, allowing repetitive regions to be bridged and enabling far more complete assemblies (Hu et al., 2021). This makes long reads particularly valuable for reconstructing bacterial chromosomes, resolving structural variation, and recovering plasmids (Hu et al., 2021). Historically, Nanopore data were limited by high error rates, but recent advances, including R10 chemistry, have substantially improved accuracy (Oxford Nanopore Technologies, 2020). Long reads now provide a practical balance between completeness and accuracy. As a result, genome assembly involves a tradeoff between contiguity and per-base accuracy and requires careful quality control and method selection.

In this project, a *Salmonella enterica* genome is assembled using Oxford Nanopore R10 long-read data and validated by alignment to a reference genome. Through this workflow, the project explores the balance between assembly contiguity and base-level accuracy in a real bacterial genome and demonstrates how raw sequencing reads are transformed into biologically interpretable genomic data.

---

## 2. Methods

### Data Acquisition

Oxford Nanopore R10.4 long-read sequencing data for *Salmonella enterica* were obtained from NCBI Sequence Read Archive under accession **SRR32410565**. The dataset consists of 196,031 reads with a median length of 4,683 bp and total yield of 809 Mb, providing approximately 160× coverage of the expected 5 Mb genome. Reference genome *S. enterica* subsp. *enterica* serovar Typhimurium str. LT2 (GCF_000006945.2) was downloaded from NCBI RefSeq, consisting of chromosome NC_003197.2 (4,857,450 bp) and plasmid NC_003277.2 (93,933 bp).

### Genome Assembly

De novo assembly was performed using **Flye v2.9.6**, optimized for high-accuracy nanopore reads:

```bash
flye --nano-hq salmonella_reads.fastq --out-dir flye_output --threads 8
```

Flye uses a repeat-graph assembly algorithm specifically designed to handle long, error-prone reads while resolving repetitive genomic regions (Kolmogorov et al., 2019). The --nano-hq parameter applies error models appropriate for R10.4 Q20+ chemistry.

### Reference Alignment

Raw nanopore reads were aligned to the reference genome using minimap2 v2.30, a fast alignment tool optimized for long reads:

```bash
minimap2 -ax lr:hq -t 8 salmonella_ref.fasta salmonella_reads.fastq > reads_aligned.sam
```

The lr:hq preset applies alignment parameters tuned for high-quality long reads. SAM output was converted to sorted, indexed BAM format for downstream analysis:

```bash
samtools view -bS reads_aligned.sam | samtools sort -o reads_sorted.bam
samtools index reads_sorted.bam
samtools faidx salmonella_ref.fasta
```

Alignment statistics were extracted using samtools coverage to evaluate breadth and depth of coverage across reference contigs.

### Variant Calling

Genetic variants were identified using Bcftools v1.16, a widely-used variant caller suitable for bacterial haploid genomes:

```bash
bcftools mpileup -Ou -f salmonella_ref.fasta reads_sorted.bam | \
bcftools call -mv -Ob -o reads_variants.bcf
bcftools view reads_variants.bcf > reads_variants.vcf
```

Variants were classified as SNPs (single nucleotide polymorphisms), insertions, or deletions based on reference and alternate allele lengths. Variant density was calculated as variants per kilobase for each reference contig.

### Data Visualization

Assembly quality metrics were visualized using R v4.3.1 with ggplot2 and patchwork packages. Coverage and variant distributions were plotted to assess genome-wide patterns. Individual variants were examined in Integrative Genomics Viewer (IGV v2.16) by loading the reference genome, BAM alignment file, and VCF variant file.

---

## 3. Results

### Assembly Quality

The Flye assembler produced a high-quality draft genome consisting of 3 contigs with a total length of 5,104,813 bp (Figure 1A). The largest contig (contig_1) spans 3.32 Mb with 153× sequencing coverage, while contig_2 is 1.68 Mb with 169× coverage. A small circular contig (contig_4) of 109 kb with 245× coverage was identified, consistent with plasmid DNA. The assembly achieved an excellent N50 of 3.32 Mb, indicating high contiguity (Figure 1D).

Topology analysis revealed one circular contig representing the plasmid element, while the two larger contigs remain linear, likely representing chromosomal fragments separated by repetitive regions (Figure 1C). Mean assembly coverage of 160× provides robust support for variant calling and structural analysis.

![figure1\_assembly\_final](https://github.com/user-attachments/assets/967a5499-ca5f-490d-b508-212d7b9b2785)

**Figure 1:** Genome Assembly Quality Assessment. De novo assembly with Flye produced 3 contigs totaling 5.10 Mb. (A) Contig lengths show two large chromosomal contigs (3.32 Mb and 1.68 Mb) and one small plasmid (0.11 Mb). (B) Sequencing coverage per contig ranges from 153× to 245×, with highest coverage on the plasmid. (C) Topology analysis identified one circular contig (plasmid). (D) Assembly summary statistics demonstrate high quality with N50 = 3.32 Mb and mean coverage 160×.

### Reference Alignment Quality

Alignment of raw reads to the S. enterica LT2 reference genome revealed substantial differences in coverage between chromosome and plasmid (Figure 2). The reference chromosome achieved 97.8% breadth of coverage with a mean depth of 151× and mapping quality (MAPQ) of 60, indicating high-confidence alignments across nearly the entire chromosome (Figure 2A,B,C). A total of 182,768 reads aligned to the chromosomal sequence (Figure 2D).

In contrast, the reference plasmid pSLT showed only 43.1% breadth of coverage despite adequate sequencing depth (82×) in aligned regions (Figure 2A,B). Mapping quality for plasmid-aligned reads was lower (MAPQ 45), and only 3,101 reads aligned to the plasmid reference (Figure 2C,D). This fragmented coverage pattern suggests the plasmid in the sequenced strain differs substantially from the pSLT reference.

![figure2\_alignment\_final](https://github.com/user-attachments/assets/934e7017-577f-4ea9-aa57-59dcbf82def8)

**Figure 2:** Reference Alignment Quality Assessment. Raw reads were aligned to the S. enterica LT2 reference genome (chromosome NC_003197.2 and plasmid NC_003277.2). (A) Reference genome coverage shows 97.8% for chromosome but only 43.1% for plasmid, indicating plasmid divergence. (B) Mean sequencing depth is 151× for chromosome and 82× for plasmid in covered regions. (C) Mapping quality scores are high for chromosome (MAPQ 60) but reduced for plasmid (MAPQ 45). (D) 182,768 reads aligned to chromosome compared to only 3,101 to plasmid.

### Variant Distribution and Density

Variant calling identified 11,465 total variants across the genome, with 98.6% being SNPs (11,302), 1.1% insertions (125), and 0.3% deletions (38) (Figure 3C). However, variant distribution was highly uneven between chromosome and plasmid. The chromosome contained 4,398 variants across 4.86 Mb, yielding a density of 0.9 variants per kb (Figure 3A,B). This low variant density is consistent with strain-level polymorphism relative to the reference LT2 strain.

In stark contrast, the plasmid harbored 7,067 variants across only 0.09 Mb, producing an extraordinary variant density of 78.52 variants per kb—87-fold higher than the chromosome (Figure 3B). This extreme plasmid variant density, combined with fragmented coverage (Figure 2A), strongly suggests the assembled plasmid represents a divergent incompatibility group rather than genuine point mutations across pSLT.

![figure3\_variants\_final](https://github.com/user-attachments/assets/23a00b7c-0eaa-4cdf-b992-1abfd2a19cb1)

**Figure 3:** Variant Landscape Analysis. Variant calling using Bcftools identified 11,465 total variants genome-wide. (A) Total variant counts show 4,398 on chromosome and 7,067 on plasmid despite the plasmid being 54× smaller. (B) Variant density reveals 0.9 variants/kb on chromosome but 78.52 variants/kb on plasmid—an 87-fold difference indicating plasmid divergence rather than true polymorphism. (C) Genome-wide variant composition is 98.6% SNPs. (D) Variant type distribution by contig shows SNPs dominate in both chromosome and plasmid.

### Visual Inspection of Variants

Genome-wide visualization in IGV confirmed relatively sparse variant distribution across the chromosome with clusters of moderate density in certain regions (Figure 4). Coverage depth appears uniform across most of the chromosome, consistent with alignment statistics (Figure 2B).

![figure4\_genome\_coverage](https://github.com/user-attachments/assets/0f823772-bfd9-4962-a841-43c2cbed2637)

**Figure 4:** Genome-Wide Coverage and Variant Distribution. IGV visualization of chromosome NC_003197.2 showing variant positions (top track) and coverage depth (middle track). Variants are distributed across the chromosome with some regional clustering. Coverage remains relatively uniform at ~150× depth across most positions.

Detailed examination of a variant-rich chromosomal region (2.36-2.37 Mb) revealed multiple well-supported SNPs with high read coverage (Figure 5). Individual variant sites show consistent support from multiple overlapping reads, indicating genuine sequence differences rather than sequencing errors.

![figure5\_variant\_region](https://github.com/user-attachments/assets/0367ded6-1a26-412a-a011-c67110e0a063)

**Figure 5:** Variant-Rich Chromosomal Region. IGV view of NC_003197.2 positions 2,360,000-2,365,000 showing multiple SNPs (purple 'I' markers in alignment tracks) supported by consistent read evidence. Coverage track shows uniform depth (~150×) across the region.

High-resolution inspection of individual variants confirmed base-level support for variant calls (Figure 6). At position 2,363,040, a clear 4bp insertion is visible with consistent support from all overlapping reads. Additional variants at nearby positions show clean SNP patterns (T→G, G→T, T→A), validating the quality of variant calls.

![figure6\_variant\_detail](https://github.com/user-attachments/assets/198c6f1e-7faf-4dfd-bb46-0b7d9cd6b31d)

**Figure 6:** Single Variant Detail View. IGV visualization at nucleotide resolution showing a 4bp insertion near position 2,363,040. The alignment shows consistent support for the insertion (purple bars) and alternative alleles across multiple reads with high base quality, validating these as true biological variants rather than sequencing errors.

---

## 4. Discussion

### Assembly Quality and Completeness

The Flye assembly successfully reconstructed the Salmonella enterica genome into three contigs totaling 5.10 Mb with excellent contiguity (N50 = 3.32 Mb). The chromosome remains fragmented into two linear contigs, likely separated by a repetitive region that was filtered during assembly due to low coverage (25× at the junction). This is a common limitation of long-read assemblers when encountering highly repetitive DNA or collapsed repeats. Despite fragmentation, the assembly covers 97.8% of the reference chromosome with high fidelity.

High sequencing coverage (160× mean) provided robust support for both assembly and variant calling. The 245× coverage on the circular plasmid contig suggests either higher copy number or preferential sequencing of smaller DNA fragments. Detection of circularity for the plasmid contig confirms successful assembly of an autonomous replicon.

### Chromosomal Variants: Strain-Level Polymorphism

The 4,398 chromosomal variants (0.9 per kb) represent genuine strain-level differences between the sequenced isolate and the reference LT2 strain. This level of divergence is consistent with intraspecies variation within S. enterica serovars (Robertson et al., 2023). The dominance of SNPs (98.6%) over indels reflects typical mutational patterns in bacterial evolution, where point mutations accumulate more frequently than insertion-deletion events.

IGV visualization confirmed that chromosomal variants are well-supported by multiple reads with high mapping quality (MAPQ 60), distinguishing true polymorphisms from sequencing errors. Variants appear distributed across the chromosome with some regional clustering, possibly reflecting horizontal gene transfer events, mobile elements, or regions under diversifying selection.

### Plasmid Divergence: Evidence for a Different Incompatibility Group

The most striking finding is the extreme plasmid divergence. Three independent lines of evidence demonstrate that the assembled plasmid is not pSLT but rather a different plasmid element:

* Low reference coverage (43.1%): Less than half of pSLT aligns to the assembly, with large gaps indicating absent or highly divergent regions.
* Extreme variant density (78.52 per kb, 87× higher than chromosome): This density is biologically implausible as true point mutations. Instead, it reflects misalignments when forcing reads from a divergent plasmid onto an incorrect reference.
* Reduced mapping quality (MAPQ 45 vs 60): Lower confidence alignments indicate the aligner struggled to place plasmid reads, consistent with sequence divergence.

McClelland et al. (2001) showed that pSLT is specific to certain S. enterica serovars and shares limited homology with plasmids from other serovars. Robertson et al. (2023) catalogued 1,044 distinct plasmid MOB-clusters across Salmonella, with 22% carrying antimicrobial resistance genes. The divergent plasmid in this strain likely belongs to a different incompatibility group or MOB-cluster.

The 109 kb assembled plasmid is larger than pSLT (94 kb), suggesting different gene content. Given that 88.3% of broad-host-range Salmonella plasmids are mobilizable (Robertson et al., 2023), this plasmid may play a role in horizontal gene transfer. Future work should annotate the assembled plasmid to identify resistance determinants, virulence factors, or mobile elements.

### Workflow Strengths and Limitations

#### Strengths

* High-quality long-read data (Q20+, 160× coverage) enabled robust assembly
* Minimap2 alignment efficiently mapped 97.8% of chromosome with high confidence
* Bcftools variant calling identified well-supported SNPs and indels
* Multi-scale visualization (genome-wide, regional, base-level) validated results

#### Limitations

* Chromosome remains fragmented due to repetitive sequences
* Bcftools is a traditional variant caller; machine learning tools like Clair3 may improve sensitivity
* Plasmid divergence prevented accurate variant characterization using the pSLT reference
* No functional annotation was performed to interpret variant consequences

### Biological and Clinical Implications

The presence of a divergent plasmid has potential clinical relevance. Salmonella plasmids frequently carry antimicrobial resistance genes, and conjugative plasmids enable rapid dissemination across strains and serovars (Laidlaw et al., 2024). Robertson et al. (2023) documented multi-plasmid AMR outbreaks, emphasizing the importance of plasmid surveillance.

Future experiments should:

* Annotate the assembled plasmid to identify resistance and virulence genes
* Determine plasmid incompatibility group and MOB-type
* Assess conjugation potential and host range
* Compare with plasmid databases to identify related elements

---

## 5. Conclusions

This analysis successfully assembled and characterized a Salmonella enterica genome using Oxford Nanopore long-read sequencing. The assembly achieved high contiguity (N50 = 3.32 Mb) and identified 11,465 variants relative to the reference genome. Chromosomal variants represent typical strain-level polymorphism, while extreme plasmid divergence indicates the presence of a different plasmid element from a distinct incompatibility group. These findings demonstrate the power of long-read sequencing for bacterial genomics and highlight the importance of plasmid diversity in Salmonella populations.

---

## 6. References

Andino, A., & Hanning, I. (2015). *Salmonella enterica*: Survival, Colonization, and Virulence Differences among Serovars. *The Scientific World Journal*, 2015(520179), 1–16. [https://pmc.ncbi.nlm.nih.gov/articles/PMC4310208/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4310208/)

Basantani, M. K., Gupta, D., Mehrotra, R., Mehrotra, S., Vaish, S., & Singh, A. (2017). An update on bioinformatics resources for plant genomics research. *Current Plant Biology*, 11–12, 33–40. [https://doi.org/10.1016/j.cpb.2017.12.002](https://doi.org/10.1016/j.cpb.2017.12.002)

Dele Ogunremi, J., Devenish, J., Amoako, K. K., Kelly, H., Dupras, A. A., Bélanger, S., & Wang, L. R. (2014). High resolution assembly and characterization of genomes of Canadian isolates of *Salmonella Enteritidis*. *BMC Genomics*, 15(1), 713. [https://doi.org/10.1186/1471-2164-15-713](https://doi.org/10.1186/1471-2164-15-713)

Freire, B., Ladra, S., & Paramá, J. R. (2022). Memory-Efficient Assembly Using Flye. *IEEE/ACM Transactions on Computational Biology and Bioinformatics*, 19(6), 3564–3577. [https://doi.org/10.1109/tcbb.2021.3108843](https://doi.org/10.1109/tcbb.2021.3108843)

Heng Li. (2020). *lh3/seqtk*. GitHub. [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

Hu, T., Chitnis, N., Monos, D., & Dinh, A. (2021). Next-generation sequencing technologies: An overview. *Human Immunology*, 82(11). [https://doi.org/10.1016/j.humimm.2021.02.012](https://doi.org/10.1016/j.humimm.2021.02.012)

igvteam. (2025). *igvteam/igv-docs*. GitHub. [https://github.com/igvteam/igv-docs](https://github.com/igvteam/igv-docs)

Johnson, R., Ravenhall, M., Pickard, D., Dougan, G., Byrne, A., & Frankel, G. (2017). Comparison of *Salmonella enterica* Serovars Typhi and Typhimurium Reveals Typhoidal Serovar-Specific Responses to Bile. *Infection and Immunity*, 86(3). [https://doi.org/10.1128/iai.00490-17](https://doi.org/10.1128/iai.00490-17)

Knodler, L. A., & Elfenbein, J. R. (2019). *Salmonella enterica*. *Trends in Microbiology*, 27(11), 964–965. [https://doi.org/10.1016/j.tim.2019.05.002](https://doi.org/10.1016/j.tim.2019.05.002)

Kolmogorov, M. (2024). *Flye: De novo assembler for single-molecule sequencing reads*. GitHub. [https://github.com/mikolmogorov/Flye](https://github.com/mikolmogorov/Flye)

Li, H. (2022). *lh3/minimap2*. GitHub. [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)

Oxford Nanopore Technologies. (2020). R10.3: the newest nanopore for high accuracy nanopore sequencing – now available in store. [https://nanoporetech.com/news/news-r103-newest-nanopore-high-accuracy-nanopore-sequencing-now-available-store](https://nanoporetech.com/news/news-r103-newest-nanopore-high-accuracy-nanopore-sequencing-now-available-store)

samtools. (2020). *samtools/samtools*. GitHub. [https://github.com/samtools/samtools](https://github.com/samtools/samtools)

Schiffer, A. M., Rahman, A., Sutton, W., Putnam, M. L., & Weisberg, A. J. (2025). A comparison of short- and long-read whole-genome sequencing for microbial pathogen epidemiology. *mSystems*, 10(12), e01426-25. [https://doi.org/10.1128/msystems.01426-25](https://doi.org/10.1128/msystems.01426-25)

### Software and Tools

Flye v2.9.6: Genome assembly - [https://github.com/fenderglass/Flye](https://github.com/fenderglass/Flye)
Minimap2 v2.30: Read alignment - [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)
Samtools v1.22: BAM file processing - [https://github.com/samtools/samtools](https://github.com/samtools/samtools)
Bcftools v1.16: Variant calling - [https://github.com/samtools/bcftools](https://github.com/samtools/bcftools)
IGV v2.16: Genome visualization - [https://software.broadinstitute.org/software/igv/](https://software.broadinstitute.org/software/igv/)
R v4.3.1: Statistical analysis and visualization
ggplot2 v4.0.2: Data visualization
patchwork v1.2.1: Multi-panel figure assembly

All analyses were performed in a conda environment (binf6110_env) on Ubuntu 24.04.

**Data Availability:** Raw sequencing data available at NCBI SRA (SRR32410565). Analysis code and processed results available in GitHub repository.
