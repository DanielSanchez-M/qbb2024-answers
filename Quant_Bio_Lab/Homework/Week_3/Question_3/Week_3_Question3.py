#!/usr/bin/env python3

#import python packages
import sys
import numpy

# Use pseudo code as skeleton for developing code for Steps 3.1 Parse the VCF file, 3.2 Allele frequency spectrum, and 3.2 Read depth distribution:
# for line in open(<vcf_file_name>):
#     if line.startswith('#'):
#         continue
#     fields = line.rstrip('\n').split('\t')

#     # grab what you need from `fields`

# Step 3.1 and Step 3.2: Parse the VCF file and Allele frequency spectrum
allele_frequency_file = open('AF.txt', mode = 'w')                        # Create a new "currently empty" file named allele_frequency_file with mode set to 'write' which will be filled in from the for loop to nest the allele frequencies from the biallelic.vcf file
allele_frequency_file.write("Allele Frequencies\n")                       # Create the header "Allele Frequncies" with a new line character in the allele_frequency_file made from the previous line of code

# For Loop:
for line in open(sys.argv[1], "r"):                                       # For loop that will open the file listed in the terminal as the "2nd" argument (biallelic.vcf) and set it as a readable file
    if line.startswith('#'):                                              # cycle through the first few lines of the file that starts with the "#" to loop through the data sections
        continue
    fields = line.rstrip('\n').split('\t')                                # remove the new line character at the end, and splits the list of fields into tabs ('\t)
    info_field = fields[7].split(";")                                     # Extract the INFO column (fields[7]) and splits it at the ";" that contain each column
    allele_frequency = info_field[3].lstrip("AF=")                        # info_fields will output all Allele Frequencies as "AF= ___" in the 4th column ([3]), thereby needing the strip from the left side to remove the "AF=" leaving only the number value
    allele_frequency_file.write(allele_frequency + '\n')                  # Overwrite the allele_frequency_file made from before the for loop, and inputs the allele_frequency with a new line character to separate the inputs in the text file
allele_frequency_file.close()                                             # Closee the allele_frequency_file

## Question 3.1: 
### Interpret this figure in two or three sentences in your own words.
### Does it look as expected? 
### Why or why not? 
### Bonus: what is the name of this distribution?

#### This figure is slightly skewed to the right with the average allele frequency just under 0.50. This is expected as most variants tend to have low frequencies due to occurances of rare mutations.


# Step 3.2: Read depth distribution
depth_file = open('DP.txt', mode = 'w')                                   # Create a new "currently empty" file named depth_file with mode set to 'write' which will be filled in from the for loop to nest the read depths from the biallelic.vcf file
depth_file.write("Read Depth of Each Variant\n")                          # Create the header "Read Depth of Each Variant" with a new line character in the depth_file made from the previous line of code

for line in open(sys.argv[1], "r"):                                       # For loop that will open the file listed in the terminal as the "2nd" argument (biallelic.vcf) and set it as a readable file
    if line.startswith('#'):                                              # cycle through the first few lines of the file that starts with the "#" to loop through the data sections
        continue
    fields = line.rstrip('\n').split('\t')                                # remove the new line character at the end, and splits the list of fields into tabs ('\t)
    for format_field in fields[9:]:                                       # For loop to go through the format fields starting from column 9
        sample_data = format_field.split(':')                             # Split the format fields by the ":" to get individual data
        read_depth = sample_data[2]                                       # Extract the read depth that is in the 3rd field
        depth_file.write(read_depth + '\n')                               # Overwrite the allele_frequency_file made from before the for loop, and inputs the allele_frequency with a new line character to separate the inputs in the text file
depth_file.close()                                                        # Close the depth_file


## Question 3.2: 
### Interpret this figure in two or three sentences in your own words. 
### Does it look as expected? 
### Why or why not? 
### Bonus: what is the name of this distribution?

#### This figure shows the number of variants at each allele frequency and matches the coverage obtained from questions 1 and 2 where a majority of the variants are covered within the 3-5 allele frequency range with a relative average of 4.0x. This figure follows a skewed poisson distribution.