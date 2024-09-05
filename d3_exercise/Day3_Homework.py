#!/usr/bin/env python3

# open file
# skip 2 lines
# split column header by tabs and skip first 2 entries
# create way to hold gene names
# create way to hold expression values
# for each line
#     split the line
#     save field 0 into gene IDs
#     save field 1 into gene names
#     save fields 2+ into expression values

import sys

import numpy

fs = open(sys.argv[1], mode = "r")                      #Open specific document at position 1 in 'read' mode

fs.readline()                                           #Skips line 1
fs.readline()                                           #Skips line 2
#Now on line 3 w/ usefull information
line = fs.readline()
fields = line.strip("\n").split("\t")                   #Strips the new line character inputted at the end | Splits the line by tabs
tissues = fields[2:]                                    #Skips the first two column entries

gene_names = []                                         #Create list for gene names
gene_IDs = []                                           #Create list for gene IDs
expression = []                                         #Create list for expression

for line in fs:
    fields = line.strip("\n").split("\t")               #Strips the new line character inputted at the end | Splits the line by tabs
    gene_IDs.append(fields[0])                          #Append the information for gene IDs
    gene_names.append(fields[1])                        #Append the information for gene names
    expression.append(fields[2:])                       #Append the information for expression values
#Now each gene_names | gene_IDs | and expression values are placed into separate lists

fs.close()

#Now we need to combine the lists to design a matrix using numpy
tissues = numpy.array(tissues)                          #Create an array for tissues (They are strings)
gene_IDs = numpy.array(gene_IDs)                        #Create an array for gene_IDs (They are strings)
gene_names = numpy.array(gene_names)                    #Create an array for gene_names (They are strings)
expression = numpy.array(expression, dtype = float)     #Create an array for expression (They are default strings need to convert to floats)

# Q4 - Use numpy arrays to calculate the mean expression values for the first 10 genes
expression_row_mean = numpy.mean(expression[:10], axis = 1)     #Will calculate the mean from each expression row up to row 10 | row is considered axis 1
#print(expression_row_mean)        # Used to check the expression means

# Q5 - Calculate the mean and median of expression values from entire data set
expression_mean = numpy.mean(expression)
expression_median = numpy.median(expression)
#print(expression_mean)                          # Mean = 16.557814350910945
#print(expression_median)                        # Median = 0.0271075

# One can infer that the difference between these statistical values means that the data distribution is positively skewed.

# Q6 - Normalize range of expression values, apply a log-transformation to the data and check mean and median values
expression = expression + 1                     # Adding a pseudocount to avoid the "can't divide by 0" error messages
expression_normalized = numpy.log2(expression)
expression_normalized_mean = numpy.mean(expression_normalized)
expression_normalized_median = numpy.median(expression_normalized)

#print(expression_normalized_mean)               # Mean = 1.1150342022364093
#print(expression_normalized_median)             # Median = 0.03858718613570538

# The mean and median values in the normalized expression set are much closer in value and more comparable compared to the unnormalized values

# Q7 - Find the expression gap for each gene between their highest and next highest expression level to identify highly tissue specific genes
# create a copy of the expression array
# Sort this copy across tissues by expression values | specify an axis to sort along
# For each gene, find the difference between the highest expression value and the second highest expression value to identify genes that are highly tissue-specific.
copy_expression_norm = numpy.copy(expression_normalized)
sorted_expression_norm = numpy.sort(copy_expression_norm, axis = 1)
diff_array = sorted_expression_norm[:,-1] - sorted_expression_norm[:,-2]
#print(diff_array)

#Q8 - Print the number of genes whose difference between the highest and second highest tissue expression is greater than 10

print(numpy.sum(diff_array >= 10))      # 33 genes