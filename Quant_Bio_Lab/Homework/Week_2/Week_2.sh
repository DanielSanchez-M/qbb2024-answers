#!/usr/bin/bash

# Question 1
# Download the files, and genes/coding exons/cCREs from UCSC Genome Table Browser

# Q1.4 Use bedtools to sort and merge the files and saved them as a "sorted file"
bedtools sort -i genes_chr1.bed > genes_sorted_chr1.bed
bedtools sort -i coding_exons.bed > exons_sorted.bed
bedtools sort -i cCREs.bed > cCREs_sorted.bed

#Then I merged the sorted files to out put the finalizaed files of sort_merge
bedtools merge -i genes_sorted_chr1.bed > genes_sort_merge_chr1.bed
bedtools merge -i exons_sorted.bed > exons_sort_merge_chr1.bed
bedtools merge -i cCREs_sorted.bed > cCREs_sort_merge_chr1.bed

# Q1.5 use bedtools to create intron feature file using the gene_sort_merge and the exons_sort_merge
# genes - exons = introns
bedtools subtract -a genes_sort_merge_chr1.bed -b exons_sort_merge_chr1.bed > introns_chr1.bed

#1.6 Use bedtools and find intervals not covered by other features
# This will take the genome file downloaded from the zip, and subtract the exons, introns and cCREs leaving only other components
bedtools subtract -a genome_chr1.bed -b exons_sort_merge_chr1.bed -b introns_chr1.bed -b cCREs_sort_merge_chr1.bed > other_chr1.bed


# Q2
# Pseudo Code
# "Create the results file with a header"
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

# "Loop through each possible MAF value"

# Place the snp files together under maf_snp_files
# Also placing the sorted and merged exon file, cCREs file, intron, and other files together under feature_files
maf_snp_files=("chr1_snps_0.1.bed" "chr1_snps_0.2.bed" "chr1_snps_0.3.bed" "chr1_snps_0.4.bed" "chr1_snps_0.5.bed")
feature_files=("exons_sort_merge_chr1.bed" "introns_chr1.bed" "cCREs_sort_merge_chr1.bed" "other_chr1.bed")
genome_file=("genome_chr1.bed")     # Determined it was best to place the genome_chr1 as a sepaerate file 

# Loop through the MAF files
for maf_file in "${maf_snp_files[@]}"; do
    # Find the total SNP coverage for the whole chromosome
    bedtools coverage -a "$genome_file" -b "$maf_file" > temp.file.txt

    # Calculate total SNPs in the chromosome from column 5
    chrm_snps_total=$(awk '{s+=$5}END{print s}' temp.file.txt)
    
    # obtain the chromosome length by suming the values in column 6
    chromosome_length=$(awk '{s+=$6}END{print s}' temp.file.txt)   

    # Calculate background SNP density by dividing the chromosome snps total by chromosome length
    background_snp_density=$(echo "$chrm_snps_total / $chromosome_length" | bc -l)

    # Loop through the feature files
    for feature_file in "${feature_files[@]}"; do
        # Find the SNP coverage from the the feature file in reference to maf_file
        feature_coverage=$(bedtools coverage -a "$feature_file" -b "$maf_file")

        # Sum SNPs for the feature
        feature_snps=$(awk '{s+=$5}END{print s}' <<< "$feature_coverage")

        # Sum the total bases of the feature
        total_bases=$(awk '{s+=$6}END{print s}' <<< "$feature_coverage")

        # Calculate SNP density and enrichment by dividing the feature_snps / total bases and then dividng those values by the background snp density
        feature_density=$(bc -l <<< "$feature_snps / $total_bases")
        enrichment=$(bc -l <<< "$feature_density / $background_snp_density")

        # Get the feature name (without the file extension)
        feature_name=$(basename "$feature_file" .bed)

        # Extract the MAF value from the maf_file name
        maf_value=$(basename "$maf_file" .bed | cut -d'_' -f3)

        # Save the result to the snp_counts.txt file
        echo -e "$maf_value\t$feature_name\t$enrichment" >> snp_counts.txt
    done
done

# "Pseudo Code
# Create the results file with a header
# Loop through each possible MAF value
#   Use the MAF value to get the file name for the SNP MAF file
#   Find the SNP coverage of the whole chromosome
#   Sum SNPs from coverage
#   Sum total bases from coverage
#   Calculate the background
#   Loop through each feature name
#     Use the feature value to get the file name for the feature file
#     Find the SNP coverage of the current feature
#     Sum SNPs from coverage
#     Sum total bases from coverage
#     Calculate enrichment
#     Save result to results file