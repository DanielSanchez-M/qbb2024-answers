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
bedtools subtract -a genome_chr1.bed -b exons_sort_merge_chr1.bed introns_chr1.bed cCREs_sort_merge_chr1.bed > other_chr1.bed
