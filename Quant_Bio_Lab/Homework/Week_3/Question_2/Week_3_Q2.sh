#!/usr/bin/env bash

# Step 2.1: Download and index the sacCer3 genome
## Pseudo Code: 
wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz   # Download the genome file
gunzip sacCer3.fa.gz                                                            # Unzip the zipped file
bwa index sacCer3.fa                                                            # Index the genome file sacCer3.fa 

## Question 2.1: How many chromosomes are in the yeast genome?
## Code: 
grep ">" sacCer3.fa | wc -l
### grep ">" sacCer3.fa             # Search the lines in the sacCer3.fa that start with ">" (chromosome files)
### wc -l                           # Count the number of lines

### Answer: There are 17 chromosomes in the yeast genome


# Step 2.2 & Step 2.4: Align your reads to the reference and Format and index your alignments
## Check code used in Step2.4
## Code: 
for my_alignments in A01_*.fastq
### A for-loop that goes over current directory files with the name A01_*.fastq (* is wildcard to cover the different file numbers)
do
    my_alignments=`basename ${my_alignments} .fastq`
### Extract the basename of the file and removes the .fastq extension and then restored as my_alignments
    bwa mem -R "@RG\tID:${my_alignments}\tSM:${my_alignments}" sacCer3.fa ${my_alignments}.fastq > ${my_alignments}.sam
### bwa men used to align the sequences from the sacCer3.fa reference genome file with the my_alignment.fastq files and saved the alignment in SAM files
    samtools sort -@ 4 -O bam -o ${my_alignments}.bam ${my_alignments}.sam
### Sort the SAM files and converts them to BAM files
    samtools index ${my_alignments}.bam
### Index the my_alignments.bam files to out put the .bai files
done
### Terminate the For-loop

# Step 2.3: Sanity check your alignments
## Question 2.2: How many total read alignments are recorded in the SAM file?
## Code: 
grep -v "^@" A01_09.sam | wc -l
### grep -v "^@" A01_09.sam         # Search through the A01_09.sam file without the lines that start with "@" to have only alignments
### wc -l                           # Count the number of lines

### Answer: 669,548 (Consistent with Exercise 1)

## Question 2.3: How many of the alignments are to loci on chromosome III?
tail -n +21 A01_09.sam | cut -d $'\t' -f 3 | grep "chrIII" | wc -l
### tail -n +21 A01_09.sam          # Skips the first 20 lines that contain header material from A01_09.sam
### cut -d $'\t' -f 3               # Extract the 3rd column from the reference sequence (Chromosome Name)
### grep "chrIII"                   # grep searches through the tabs searching for "chrIII"
### wc -l                           # Count the number of lines

### Answer: 17,815

# Step 2.5: Visualize your alignments
## Code: igv in Terminal, then open the .bai file, change genome to sacCer3 and zoom into visualize alignments

## Question 2.4: Does the depth of coverage appear to match that which you estimated in Step 1.3? Why or why not?
### Answer: Overall, yes the depth of coverage apeears to match well with what I estimated in Step 1.3 (4). There are areas that have coverages of more than 4, but also some areas that were not covered. So it can be safe to assume that the average is around 4.

## Question 2.5: Set your window to chrI:113113-113343 (paste that string in the search bar and click enter). 
## How many SNPs do you observe in this window? Are there any SNPs about which you are uncertain? Explain your answer.
### Answer: There are 3 SNPs observed in this window. 2 of the SNPs at chrI:113,132 and chrI:113,207 have about 4 read coverages, while the other SNP at chrI:113,326 only has 1 coverage leaving me to speculate with a degree of uncertainty on this SNP.

## Question 2.6: Set your window to chrIV:825548-825931. 
## What is the position of the SNP in this window? Does this SNP fall within a gene?
### Answer: The SNP is positioned at 825,834 which falls between scc2 and sas4