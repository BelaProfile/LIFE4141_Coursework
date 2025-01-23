**Cochlearia excelsa Genomic Signatures of Adaptation: Understanding Metal Tolerance and Population Divergence**

We investigated the adaptation of three populations of Cochlearia excelsa (LAB, NEN and ODN) in different habitats. More specifically, it examines how stressors such as heavy metal contamination affect genetic diversity and divergence. PCA and Fst analyses revealed strong genetic differentiation among the populations. Two genes were found to be of particular interest—IRT1 and HMA2. Both are important in metal transport and adaptation. These results shed light on evolutionary dynamics, demonstrating the impact of environmental stressors in facilitating local adaptation in natural populations.

**Materials and Methods**
Study design
Examining the genetic diversity, population structure, and selective pressures seen in Cochlearia excelsa populations from several ecological contexts was the aim of this study. ODN, which was collected from a mine-adjacent area; NEN, which was collected from a metal-contaminated river; and LAB, a non-mine community, were the three populations that were analysed. By combining genomic approaches including Principal Component Analysis (PCA), Fst-based selection scans, and candidate gene identification, the goal was to understand genetic difference and adaptive evolution.

Data collection
  1.	Genomic Data: Three million base pairs of biallelic SNPs are included in the variant data in the VCF file LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz.
  2.	Genome Annotation: Gene locations and characteristics are included in the Cochlearia excelsa annotation file C_excelsa_V5_braker2_wRseq.gff3.
  3.	Functional Homologs: 1-2-1_hits_all_gene_descriptions.tsv contains information about Arabidopsis thaliana homologs that connect genes to established functions.

Data processing
  1.	Genomic Filtering: Following VCFtools' processing of the raw VCF file, SNPs with a minor allele frequency (MAF) ≥ 0.05 were kept.
    •	Avoid including locations if more than 50% of the data is missing.
    •	Get rid of people who have poor genotype call rates.
    •	A filtered dataset including 35,949 high-quality SNPs and 22 individuals was the end result.
  2.	Formatting: The SNP data that had been filtered was converted into formats that bedtools for genome intersection and R-based research could use.

Analytical Procedure
Population Structure Analysis
  1.	Principal Component Analysis (PCA):
    •	SNP genotypes were transformed into a genlight object for PCA computation using the adegenet package in R.
    •	Scatter plots were used to visualise clustering by population.

  2.	Hierarchical Clustering:
    •	Dendrograms and heatmaps showing the genetic links between populations were created using genetic distance matrices obtained from PCA.

Fst-Based Selection Scan
  1.	Pairwise Fst Calculation:
    •	For LAB vs. NENT and LAB vs. ODN, pairwise genome-wide Fst values were computed using VCFtools with 1 kb sliding windows.
    •	In outlier regions, or the top 1% of Fst values, candidate selective sweep loci were discovered.

  2.	Visualization:
    •	Manhattan plots displaying genome-wide Fst values with annotations in the upper sections were used to identify significant outliers.

Candidate Gene Identification
  1.	Gene Overlap:
    •	Using bedtools intersect, genes in C_excelsa_V5_braker2_wRseq.gff3 were overlapped with Fst outlier areas; the resulting intersect files contained genomic coordinates and corresponding gene IDs.
  2.	Gene Annotation:
    •	Gene IDs from the intersect files were annotated using homolog data from 1-2-1_hits_all_gene_descriptions.tsv in order to infer gene functions.
    •	R-based processes and Python scripts ensured accurate gene mapping and functional annotation.

Functional Insights
  1.	Metal-Related Gene Identification:
    •	Annotated candidate genes were checked for metal-related activities using descriptions from the homolog database.
    •	For extra validation, cross-referencing with the A. thaliana TAIR database was required.

Visualization
All of the data and visualisations were produced using custom R scripts, showing:
  1.	PCA scatter plots and heatmaps with hierarchical grouping.
  2.	Manhattan plots of the genome's Fst values.

Tools and Software
  •	R Packages: vcfR, ade4, adegenet, adegraphics, heatmaply, ggplot2.
  •	VCFtools: Genomic filtering and Fst calculation.
  •	bedtools: Gene overlap analysis.
  •	Python: Functional annotation and custom data processing scripts.

