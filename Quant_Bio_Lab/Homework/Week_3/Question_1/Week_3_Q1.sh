#!/usr/bin/env bash

# Question 1.1: How long are the sequencing reads?
## Code: 
head -n 2 A01_09.fastq | tail -n 1 | wc -c
### head -n 2 A01_09.fastq              # Grabs the 2nd line from the fastQ file (sequence)
### tail -n 1                           # extracts the read sequence
### wc - c                              # counts the number of characters = length of sequence

## Answer:  The sequence is 77bp long


# Question 1.2. How many reads are present within the file?
## Code:
wc -l A01_09.fastq | awk '{print $1/4}'
### wc -l A01_09.fastq                  # counts the number of lines in the file
### awk '{print $1/4}'                  # divides the total line count by 4 (Each read in a FASTQ file has 4 lines of information)

## Answer:  There are about 669,548 reads in the A01_09.fastq file


# Question 1.3: Given your answers to 1 and 2, as well as knowledge of the length of the S. cerevisiae reference genome, what is the expected average depth of coverage?
## Code:
genome_size = 12100000                  # S. cerevisiae reference genome length is about 12.1Mb
echo "(669548 * 77) / 12100000" | bc    # Outputs the string as a methematical operation due to "bc"

## Answer:  Expected average depth of coverage = 4


# Question 1.4: While you do not need to repeat for all samples, looking at the size of the files can give us information about whether we have similar amounts of data from other samples. Use the du command to check the file sizes of the rest of the samples. Which sample has the largest file size (and what is that file size, in megabytes)? Which sample has the smallest file size (and what is that file size, in megabytes)?
## Code:
du -h *.fastq | sort -hr
### du -h *.fastq                       # du = "Disk Usage" to find the size of all ___.fastq files present in the directory
### sort -hr                            # Sorts the lists by descending order

## Answer:  Largest File = A01_62.fastq (149M)      Smallest File = A01_27.fastq (110M)


# Question 1.5: Run the program FastQC on your samples (with default settings). Open the HTML report for sample A01_09.
## Code:
fastqc A01_09.fastq                     # generates a quality control report of the file

## Answer:  
### What is the median base quality along the read?
#### The median base quality along the reads is about 35-36

### How does this translate to the probability that a given base is an error?
#### From the formula given in class of P = 10 ^ (-Quality score / 10), then the probability that base is incorrect can be determined to be within 0.025% and 0.0316%, a very low rate for error.

### Do you observe much variation in quality with respect to the position in the read?
#### Yes, from the quality scores there seems to be more variation on either end of the reads, while the centers have what appears to be higher quality (smaller boxplot spread).