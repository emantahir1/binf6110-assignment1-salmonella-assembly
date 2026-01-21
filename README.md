# Assignment 1 – Part 1  
## Genome Assembly of *Salmonella enterica*

## 1. Introduction

### Biological Background and Importance  
*Salmonella enterica* is a Gram-negative, facultative intracellular bacterium and a major cause of foodborne illness worldwide (Knodler & Elfenbein, 2019). It includes many different serovars. Some, such as *S. Typhi*, are restricted to specific hosts, while others, like *S. Typhimurium*, infect a broad range of hosts and commonly cause gastroenteritis in both humans and animals (Johnson et al., 2017). These groups differ in how they spread, which hosts they infect, and how severe the resulting disease can be. In practice, this means that *Salmonella enterica* is not a single, uniform pathogen, but a genetically diverse species with meaningful biological variation (Andino & Hanning, 2015). Because of this diversity, high-quality genome assemblies are especially important for *Salmonella enterica*. They allow closely related strains to be distinguished, make outbreak tracking possible, and support the identification of genes linked to virulence and antimicrobial resistance (Dele Ogunremi et al., 2014). Together, these applications form the basis of modern public health surveillance.

### Challenges in Genome Assembly and Project Goals  
Genome assembly is how a full genome is put back together from the many short DNA fragments produced by sequencing (Basantani et al., 2017). This process turns raw reads into a continuous sequence that represents an organism’s genome, making it possible to examine genes, genome structure, and variation (Basantani et al., 2017).

Genome assembly is computationally challenging because sequencing reads contain errors and genomes often include repetitive regions that are difficult to place correctly (Schiffer et al., 2025). When identical or near-identical sequences occur in multiple locations, assembly algorithms may be unable to determine their true origin, leading to fragmented or misassembled genomes (Schiffer et al., 2025). These challenges are especially pronounced when using short-read data. Short-read platforms such as Illumina generate highly accurate reads, but their limited length means they often cannot span repetitive elements, resulting in assemblies broken into many small contigs (Hu et al., 2021).

Long-read sequencing offers a different approach. Technologies such as Oxford Nanopore produce reads that span several kilobases, allowing repetitive regions to be bridged and enabling far more complete assemblies (Hu et al., 2021). This makes long reads particularly valuable for reconstructing bacterial chromosomes, resolving structural variation, and recovering plasmids (Hu et al., 2021). Historically, Nanopore data were limited by high error rates, but recent advances, including R10 chemistry, have substantially improved accuracy (Oxford Nanopore Technologies, 2020). Long reads now provide a practical balance between completeness and accuracy. As a result, genome assembly involves a tradeoff between contiguity and per-base accuracy and requires careful quality control and method selection.

In this project, a *Salmonella enterica* genome is assembled using Oxford Nanopore R10 long-read data and validated by alignment to a reference genome. Through this workflow, the project explores the balance between assembly contiguity and base-level accuracy in a real bacterial genome and demonstrates how raw sequencing reads are transformed into biologically interpretable genomic data.

---

## 2. Proposed Methods

### Data Acquisition and Quality Control  
Raw sequencing data for *Salmonella enterica* will be retrieved from the NCBI Sequence Read Archive. The dataset consists of long reads generated using Oxford Nanopore R10 chemistry. Initial quality control will be performed using **seqtk (v1.5)** to confirm basic FASTQ integrity and summarize read characteristics (Heng Li, 2020). Specifically, seqtk will be used to verify that the FASTQ file can be parsed correctly and contains the expected number of reads, summarize read length distributions, and summarize per-read quality score distributions. Read length summary statistics will include estimation of the read N50, and quality summaries will be used to confirm that the dataset meets the expected Q20+ standard for high-quality long-read assembly. These QC checks ensure the input data are suitable for downstream assembly.

### Genome Assembly  
De novo genome assembly will be performed using **Flye (v2.9)**, a long-read assembler designed for error-prone long-read sequencing data and commonly used for bacterial genome assembly (Kolmogorov, 2024). Flye will be run using the `--nano-hq` option to optimize assembly for high-accuracy Nanopore reads, and multithreading will be enabled using `--threads 8` to improve computational efficiency (Kolmogorov, 2024). All other parameters will be kept at their default settings unless QC results suggest that filtering or parameter adjustments are necessary.

### Reference Alignment and Visualization  
To evaluate the assembly, the resulting contigs will be aligned to a standard *Salmonella enterica* reference genome obtained from NCBI RefSeq. Alignment will be performed using **Minimap2 (v2.25)**, which is optimized for fast and accurate alignment of long sequences (Li, 2022). The alignment will be generated in SAM format using the Nanopore-optimized preset (`-ax map-ont`) (Li, 2022).

Alignment files will then be processed using **Samtools (v1.23)** (samtools, 2020). The SAM file will be converted to BAM, sorted by genomic coordinates, and indexed to support efficient visualization. The finalized BAM alignment will be visualized using the **Integrative Genomics Viewer (IGV v2.19.7)** to inspect overall coverage patterns, evaluate structural consistency between the assembly and the reference, and examine representative regions for potential discrepancies (igvteam, 2025).

Together, these steps form a complete long-read genome assembly and validation workflow, transforming raw Oxford Nanopore sequencing reads into a coherent *Salmonella enterica* genome assembly that can be evaluated against a trusted reference sequence.

---

## References

Andino, A., & Hanning, I. (2015). *Salmonella enterica*: Survival, Colonization, and Virulence Differences among Serovars. *The Scientific World Journal*, 2015(520179), 1–16. https://pmc.ncbi.nlm.nih.gov/articles/PMC4310208/

Basantani, M. K., Gupta, D., Mehrotra, R., Mehrotra, S., Vaish, S., & Singh, A. (2017). An update on bioinformatics resources for plant genomics research. *Current Plant Biology*, 11–12, 33–40. https://doi.org/10.1016/j.cpb.2017.12.002

Dele Ogunremi, J., Devenish, J., Amoako, K. K., Kelly, H., Dupras, A. A., Bélanger, S., & Wang, L. R. (2014). High resolution assembly and characterization of genomes of Canadian isolates of *Salmonella Enteritidis*. *BMC Genomics*, 15(1), 713. https://doi.org/10.1186/1471-2164-15-713

Freire, B., Ladra, S., & Paramá, J. R. (2022). Memory-Efficient Assembly Using Flye. *IEEE/ACM Transactions on Computational Biology and Bioinformatics*, 19(6), 3564–3577. https://doi.org/10.1109/tcbb.2021.3108843

Heng Li. (2020). *lh3/seqtk*. GitHub. https://github.com/lh3/seqtk

Hu, T., Chitnis, N., Monos, D., & Dinh, A. (2021). Next-generation sequencing technologies: An overview. *Human Immunology*, 82(11). https://doi.org/10.1016/j.humimm.2021.02.012

igvteam. (2025). *igvteam/igv-docs*. GitHub. https://github.com/igvteam/igv-docs

Johnson, R., Ravenhall, M., Pickard, D., Dougan, G., Byrne, A., & Frankel, G. (2017). Comparison of *Salmonella enterica* Serovars Typhi and Typhimurium Reveals Typhoidal Serovar-Specific Responses to Bile. *Infection and Immunity*, 86(3). https://doi.org/10.1128/iai.00490-17

Knodler, L. A., & Elfenbein, J. R. (2019). *Salmonella enterica*. *Trends in Microbiology*, 27(11), 964–965. https://doi.org/10.1016/j.tim.2019.05.002

Kolmogorov, M. (2024). *Flye: De novo assembler for single-molecule sequencing reads*. GitHub. https://github.com/mikolmogorov/Flye

Li, H. (2022). *lh3/minimap2*. GitHub. https://github.com/lh3/minimap2

Oxford Nanopore Technologies. (2020). R10.3: the newest nanopore for high accuracy nanopore sequencing – now available in store. https://nanoporetech.com/news/news-r103-newest-nanopore-high-accuracy-nanopore-sequencing-now-available-store

samtools. (2020). *samtools/samtools*. GitHub. https://github.com/samtools/samtools

Schiffer, A. M., Rahman, A., Sutton, W., Putnam, M. L., & Weisberg, A. J. (2025). A comparison of short- and long-read whole-genome sequencing for microbial pathogen epidemiology. *mSystems*, 10(12), e01426-25. https://doi.org/10.1128/msystems.01426-25
