##COMPLETED ASSESSMENT CODE###




# Load required libraries
library(vcfR)
library(ade4)
library(adegenet)
library(adegraphics)
library(heatmaply)
library(ggplot2)

#### Importing SNP Data ####

# Read the VCF file
vcf <- read.vcfR("LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf") 

# Extract genotypes from the VCF and transform to numeric values
df <- extract.gt(vcf)
df[df == "0|0" | df == "0/0"] <- 0
df[df == "0|1" | df == "1|0"] <- 1
df[df == "1|1" | df == "1/1"] <- 2
df <- data.frame(apply(df, 2, function(y) as.numeric(as.character(y))))

# Get sample IDs
sample_ID <- addID(vcf) 

# Convert the VCF data to a genlight object
vcf_gen <- vcfR2genlight(vcf)

# Set population groups for coloring PC axes
pop(vcf_gen) <- substr(indNames(vcf_gen), 1, 3) 

# Check populations
print(pop(vcf_gen))

#### Object Information ####

# Check ploidy for each sample (should be diploid)
print(ploidy(vcf_gen))

# Check sample names
print(indNames(vcf_gen))

# Verify positions match the VCF file
print(position(vcf_gen))

#### Principal Component Analysis (PCA) ####

# Perform PCA with 6 axes
pca_object <- glPca(vcf_gen, nf = 6)

# Scatter plot of eigenvalues (variance explained by each PC axis)
scatter(pca_object, ratio = 0.2, posi = 'bottomright') 

# Plot variable loadings for PCA
loadingplot(pca_object)

# Calculate and print the variance explained by PC1 and PC2
pc1_var <- pca_object$eig[1] / sum(pca_object$eig)
pc2_var <- pca_object$eig[2] / sum(pca_object$eig)
print(paste("Proportion of variance explained by PC1:", round(pc1_var * 100, 1), "%"))
print(paste("Proportion of variance explained by PC2:", round(pc2_var * 100, 1), "%"))

# Plot PCA colored by groups
plot1 <- s.class(
  pca_object$scores, pop(vcf_gen), 
  xax = 1, yax = 2, col = colors()[c(131, 133, 139)], 
  ellipseSize = 0, starSize = 0, ppoints.cex = 1, paxes.draw = TRUE, 
  pgrid.draw = FALSE, ppoints.pch = c(25, 22, 19),
  xlab = paste("PC1 (", round(pc1_var * 100, 1), "%)", sep = ""), 
  ylab = paste("PC2 (", round(pc2_var * 100, 1), "%)", sep = ""),
  plegend.drawKey = TRUE, plegend.size = 0.8
)
print(plot1)

# Plot PCA without coloring groups
plot2 <- s.label(
  pca_object$scores, xax = 1, yax = 2, ppoints.col = c("red"), 
  plabels = list(box = list(draw = FALSE), optim = TRUE), 
  paxes.draw = TRUE, pgrid.draw = FALSE, 
  plabels.cex = 1, plot = FALSE,
  xlab = paste("PC1 (", round(pc1_var * 100, 1), "%)", sep = ""), 
  ylab = paste("PC2 (", round(pc2_var * 100, 1), "%)", sep = "")
)
print(plot2)

# Combine and display both PCA plots in one layout
combined_plots <- ADEgS(c(plot1, plot2), layout = c(1, 2))
print(combined_plots)

#### Heatmap of PCA Distance Matrix ####

# Calculate a distance matrix using PCA scores
dist_matrix <- dist(pca_object$scores)

# Convert the distance matrix to a square matrix for heatmap
dist_matrix_square <- as.matrix(dist_matrix)

# Add sample names to rows and columns
rownames(dist_matrix_square) <- indNames(vcf_gen)
colnames(dist_matrix_square) <- indNames(vcf_gen)

# Plot heatmap using heatmaply
heatmaply(
  dist_matrix_square,
  dendrogram = "both", # Add hierarchical clustering to rows and columns
  xlab = "Samples",
  ylab = "Samples",
  main = "Heatmap of PCA Distance Matrix",
  colors = viridis::viridis(100), # Viridis color scheme
  margins = c(40, 40) # Adjust margins for better readability
)




#!/bin/bash

# Input VCF file
vcf_file="LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz"

# Population prefixes
populations=("LAB" "NEN" "ODN")

# Check if the VCF file is gzipped and decompress it if necessary
if [[ "$vcf_file" == *.gz ]]; then
echo "Decompressing VCF file..."
gunzip -k "$vcf_file"  # Keep the original gzipped file
vcf_file="${vcf_file%.gz}"
fi

# Extract the header line containing sample IDs
header=$(grep "^#CHROM" "$vcf_file")

# Loop through each population prefix
for pop in "${populations[@]}"; do
echo "Processing population: $pop"
# Extract sample IDs matching the population prefix and save to a file
echo "$header" | tr '\t' '\n' | grep "^$pop" > "${pop}_population.txt"
echo "File created: ${pop}_population.txt"
done

# FST Analysis for LAB vs NENT

# Output directory for LAB vs NENT FST results
output_dir_nent="LAB_NENT_output"

# Population files for LAB vs NENT
lab_population="LAB_population.txt"
nent_population="NEN_population.txt"

# Create output directory for LAB vs NENT if it does not exist
mkdir -p "$output_dir_nent"

# Perform FST analysis for LAB vs NENT
vcftools --vcf "$vcf_file" \
--max-missing 0.8 \
--maf 0.05 \
--weir-fst-pop "$lab_population" \
--weir-fst-pop "$nent_population" \
--fst-window-size 5000 \
--fst-window-step 5000 \
--out "$output_dir_nent/LABvNENT_WIN5K"

# Confirmation message for LAB vs NENT
echo "FST analysis completed for LAB vs NENT. Results saved in $output_dir_nent"

# FST Analysis for LAB vs ODN

# Output directory for LAB vs ODN FST results
output_dir_odn="LAB_ODN_output"

# Population files for LAB vs ODN
odn_population="ODN_population.txt"

# Create output directory for LAB vs ODN if it does not exist
mkdir -p "$output_dir_odn"

# Perform FST analysis for LAB vs ODN
vcftools --vcf "$vcf_file" \
--max-missing 0.8 \
--maf 0.05 \
--weir-fst-pop "$lab_population" \
--weir-fst-pop "$odn_population" \
--fst-window-size 5000 \
--fst-window-step 5000 \
--out "$output_dir_odn/LABvODN_WIN5K"

# Confirmation message for LAB vs ODN
echo "FST analysis completed for LAB vs ODN. Results saved in $output_dir_odn"



### FST ###

# FST Analysis for LAB-NENT

# Load required library
library(ggplot2)

#### Read and Process FST Results ####

# Read Fst results for LAB-NENT
fst_results_nent <- read.table("LAB_NENT_output/LABvNENT_WIN5K.windowed.weir.fst", header = TRUE)

# Inspect the data
cat("First few rows of LAB-NENT Fst results:\n")
print(head(fst_results_nent))

# Remove rows with NA values in WEIGHTED_FST
fst_results_nent <- fst_results_nent[!is.na(fst_results_nent$WEIGHTED_FST), ]

# Calculate the 99th percentile as the Fst cutoff
fst_cutoff_nent <- quantile(fst_results_nent$WEIGHTED_FST, 0.99)
cat("99th percentile Fst cutoff for LAB-NENT:", fst_cutoff_nent, "\n")

# Filter candidate regions based on the Fst cutoff
outlier_regions_nent <- fst_results_nent[fst_results_nent$WEIGHTED_FST >= fst_cutoff_nent, ]

# Display the candidate regions
cat("Candidate regions for LAB-NENT (top 1% Fst):\n")
print(outlier_regions_nent)

#### Visualization: Manhattan Plot ####

# Plot the Manhattan plot highlighting outlier regions
ggplot(fst_results_nent, aes(x = BIN_START, y = WEIGHTED_FST, color = as.factor(CHROM))) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_point(data = outlier_regions_nent, aes(x = BIN_START, y = WEIGHTED_FST), color = "red", size = 2) +
  geom_hline(yintercept = fst_cutoff_nent, color = "blue", linetype = "dashed") +
  labs(
    title = "Fst Outlier Regions (LAB-NENT)",
    x = "Genomic Position",
    y = "Fst"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# FST Analysis for LAB-ODN

#### Read and Process FST Results ####

# Read Fst results for LAB-ODN
fst_results_odn <- read.table("LAB_ODN_output/LABvODN_WIN5K.windowed.weir.fst", header = TRUE)

# Inspect the data
cat("First few rows of LAB-ODN Fst results:\n")
print(head(fst_results_odn))

# Remove rows with NA values in WEIGHTED_FST
fst_results_odn <- fst_results_odn[!is.na(fst_results_odn$WEIGHTED_FST), ]

# Calculate the 99th percentile as the Fst cutoff
fst_cutoff_odn <- quantile(fst_results_odn$WEIGHTED_FST, 0.99)
cat("99th percentile Fst cutoff for LAB-ODN:", fst_cutoff_odn, "\n")

# Filter candidate regions based on the Fst cutoff
outlier_regions_odn <- fst_results_odn[fst_results_odn$WEIGHTED_FST >= fst_cutoff_odn, ]

# Display the candidate regions
cat("Candidate regions for LAB-ODN (top 1% Fst):\n")
print(outlier_regions_odn)

#### Visualization: Manhattan Plot ####

# Plot the Manhattan plot highlighting outlier regions
ggplot(fst_results_odn, aes(x = BIN_START, y = WEIGHTED_FST, color = as.factor(CHROM))) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_point(data = outlier_regions_odn, aes(x = BIN_START, y = WEIGHTED_FST), color = "red", size = 2) +
  geom_hline(yintercept = fst_cutoff_odn, color = "blue", linetype = "dashed") +
  labs(
    title = "Fst Outlier Regions (LAB-ODN)",
    x = "Genomic Position",
    y = "Fst"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

write.table(outlier_regions_nent[, c("CHROM", "BIN_START", "BIN_END")],
            file = "LAB_NENT_outliers.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(outlier_regions_odn[, c("CHROM", "BIN_START", "BIN_END")],
            file = "LAB_ODN_outliers.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





###Intersect FST Outliers with GFF###
#!/bin/bash

# Activate the bedtools environment
conda activate bedtools

# Input files
gff_file="C_excelsa_V5_braker2_wRseq.gff3"
bed_file_nent="LAB_NENT_outliers.bed"
bed_file_odn="LAB_ODN_outliers.bed"

# Output files
intersect_nent="LAB_NENT_output/LABvNENT_WIN5K_gffoverlapped.bed"
intersect_odn="LAB_ODN_output/LABvODN_WIN5K_gffoverlapped.bed"

# LAB-NENT: Intersect FST outliers with GFF
bedtools intersect -a "$gff_file" -b "$bed_file_nent" -wa -wb > "$intersect_nent"
echo "LAB-NENT gene names overlap saved in $intersect_nent."

# LAB-ODN: Intersect FST outliers with GFF
bedtools intersect -a "$gff_file" -b "$bed_file_odn" -wa -wb > "$intersect_odn"
echo "LAB-ODN gene names overlap saved in $intersect_odn."

# Deactivate conda environment
conda deactivate



###matching gene names with homologs###
##labvsnent##
#!/bin/bash

# Input files
bed_file="LAB_NENT_output/LABvNENT_WIN5K_gffoverlapped.bed"  # Replace with LAB_ODN file for ODN analysis
homolog_table="1-2-1_hits_all_gene_descriptions.tsv"

# Output file
output_file="LAB_NENT_output/labvsnent_homologs.txt"  # Replace with LAB_ODN path for ODN analysis

# Extract gene IDs from the BED file
awk -F"[=;]" '/ID=/ {print $2}' "$bed_file" | sort | uniq > extracted_gene_ids.txt

# Match extracted gene IDs with the homolog table
grep -Ff extracted_gene_ids.txt "$homolog_table" > "$output_file"

# Check if the output file contains data
if [ -s "$output_file" ]; then
echo "Homologs identified and saved in $output_file."
else
  echo "No homologs found. Please check input files."
fi

# Clean up temporary file
rm extracted_gene_ids.txt


##labvsodn##
#!/bin/bash

# Input files
bed_file="LAB_ODN_output/LABvODN_WIN5K_gffoverlapped.bed"
homolog_table="1-2-1_hits_all_gene_descriptions.tsv"

# Output file
output_file="LAB_ODN_output/labvsodn_homologs.txt"

# Extract gene IDs from the BED file
awk -F"[=;]" '/ID=/ {print $2}' "$bed_file" | sort | uniq > extracted_gene_ids_odn.txt

# Match extracted gene IDs with the homolog table
grep -Ff extracted_gene_ids_odn.txt "$homolog_table" > "$output_file"

# Check if the output file contains data
if [ -s "$output_file" ]; then
echo "Homologs identified and saved in $output_file."
else
  echo "No homologs found. Please check input files."
fi

# Clean up temporary file
rm extracted_gene_ids_odn.txt


##IDENTIFYING HOMOLOGS###
#!/bin/bash

# Input files
bed_file_nent="LAB_NENT_output/LABvNENT_WIN5K_gffoverlapped.bed"
bed_file_odn="LAB_ODN_output/LABvODN_WIN5K_gffoverlapped.bed"
homolog_table="1-2-1_hits_all_gene_descriptions.tsv"

# Output files
output_homologs_nent="LAB_NENT_output/labvsnent_homologs.txt"
output_homologs_odn="LAB_ODN_output/labvsodn_homologs.txt"

# Extract gene IDs for LAB-NENT and match with homolog table
awk -F"[=;]" '/ID=/ {print $2}' "$bed_file_nent" | sort | uniq > extracted_gene_ids_nent.txt
grep -Ff extracted_gene_ids_nent.txt "$homolog_table" > "$output_homologs_nent"
echo "Homologs for LAB-NENT saved in $output_homologs_nent."

# Extract gene IDs for LAB-ODN and match with homolog table
awk -F"[=;]" '/ID=/ {print $2}' "$bed_file_odn" | sort | uniq > extracted_gene_ids_odn.txt
grep -Ff extracted_gene_ids_odn.txt "$homolog_table" > "$output_homologs_odn"
echo "Homologs for LAB-ODN saved in $output_homologs_odn."

# Clean up temporary files
rm extracted_gene_ids_nent.txt extracted_gene_ids_odn.txt



##IDENTIFY METAL RELATED CANDIDATES AND SHARED GENES###
#!/bin/bash

# Input homolog files
homologs_nent="LAB_NENT_output/labvsnent_homologs.txt"
homologs_odn="LAB_ODN_output/labvsodn_homologs.txt"

# Output files
metal_related_nent="LAB_NENT_output/labvsnent_metal_related.txt"
metal_related_odn="LAB_ODN_output/labvsodn_metal_related.txt"
shared_candidates="shared_candidates_between_nent_odn.txt"

# Keywords for metal-related functions
metal_keywords=("iron" "zinc" "copper" "manganese" "metal")

# LAB-NENT: Filter for metal-related homologs
grep -Ei "$(printf "|%s" "${metal_keywords[@]}")" "$homologs_nent" > "$metal_related_nent"
if [ -s "$metal_related_nent" ]; then
echo "Metal-related homologs for LAB-NENT saved in $metal_related_nent."
else
  echo "No metal-related homologs found for LAB-NENT."
fi

# LAB-ODN: Filter for metal-related homologs
grep -Ei "$(printf "|%s" "${metal_keywords[@]}")" "$homologs_odn" > "$metal_related_odn"
if [ -s "$metal_related_odn" ]; then
echo "Metal-related homologs for LAB-ODN saved in $metal_related_odn."
else
  echo "No metal-related homologs found for LAB-ODN."
fi

# Identify shared candidates between LAB-NENT and LAB-ODN
comm -12 <(sort "$homologs_nent") <(sort "$homologs_odn") > "$shared_candidates"
if [ -s "$shared_candidates" ]; then
echo "Shared candidates between LAB-NENT and LAB-ODN saved in $shared_candidates."
else
  echo "No shared candidates found between LAB-NENT and LAB-ODN."
fi









