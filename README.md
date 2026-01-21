# Brainstorming

## Introduction (≈ 1-1.5 pages)

Explain what genome assembly is
- Explain why it matters (especially for bacteria & pathogens)

Describe the challenges:
- Compare Short-read vs long-read sequencing
- Pros & cons of long reads (Nanopore)

Introduce Salmonella enterica:
- Why it’s important (foodborne illness, outbreaks, surveillance)

Clearly state your goal:
- Assemble a genome from Nanopore reads
- Compare it to a reference genome
- Identify and visualize differences

-- Show understanding, discuss tradeoffs, Cite at least 5 real papers

## Proposed Methods (≈ 0.5 page)

Starting data:
- Oxford Nanopore FASTQ (R10 chemistry, Q20+, N50 5–15 kb)

Quality control:
- How i will inspect reads

Assembly:
- Tool (e.g., Flye)
- Why it’s appropriate
- Version
- Key parameters (e.g., --nano-hq, threads)

Reference genome:
- Download from NCBI
Alignment & comparison:
- Tool (e.g., minimap2)
- Variant calling approach

Visualization:
- How differences will be viewed
