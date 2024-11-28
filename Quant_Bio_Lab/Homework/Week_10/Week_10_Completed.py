#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly.express as px 
import plotly

na = np.newaxis

# 1: Loading the Image Data

## Load directory with images and define genes, fields, and channel
Week_10_Folder = "/Users/cmdb/qbb2024-answers/Quant_Bio_Lab/Homework/Week_10/"
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "PCNA", "nascentRNA"]

## Open an empty dictionary 
images = {} 

for gene in genes:
    for field in fields:
        image = np.zeros((520, 616, 3), np.uint16)
        for i, channel in enumerate(channels):
            image_name = Week_10_Folder + gene + '_' + field + '_' + channel + '.tif' 
            image[:, :, i] = imageio.v3.imread(image_name)
        
        image = image.astype(np.float32)
        for i in range(3):
            image[:, :, i] -= np.amin(image[:, :, i])
            image[:, :, i] /= np.amax(image[:, :, i])
    
        image = (np.minimum(255, np.floor(image * 256))).astype(np.uint8) #resize
        images[f"{gene}_{field}"] = image
        
        # 2.1: Creating binary mask from DAPI channel
        dapi_mean = np.mean(image[:, :, 0])
        mask = image[:, :, 0] >= dapi_mean

print(images.keys())

# Exercise 2.2: Find labels for each image based on DAPI mask
## From Live Code

def find_labels(mask):
    # Set initial label
    l = 0
    # Create array to hold labels
    labels = np.zeros(mask.shape, np.int32)
    # Create list to keep track of label associations
    equivalence = [0]
    # Check upper-left corner
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    # For each non-zero column in row 0, check back pixel label
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If back pixel has a label, use same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    # For each non-zero row
    for x in range(1, mask.shape[0]):
        # Check left-most column, up  and up-right pixels
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If up pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If up-right pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        # For each non-zero column except last in nonzero rows, check up, up-right, up-right, up-left, left pixels
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If up pixel has label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If not up but up-right pixel has label, need to update equivalence table
                    if mask[x - 1, y - 1]:
                        # If up-left pixel has label, relabel up-right equivalence, up-left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If left pixel has label, relabel up-right equivalence, left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither up-left or left pixels are labeled, use up-right equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If not up, or up-right pixels have labels but up-left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If not up, up-right, or up-left pixels have labels but left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, add new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        # Check last pixel in row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # if up pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # if not up but up-left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # if not up or up-left but left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = np.array(equivalence)
    # Go backwards through all labels
    for i in range(1, len(equivalence))[::-1]:
        # Convert labels to the lowest value in the set associated with a single object
        labels[np.where(labels == i)] = equivalence[i]
    # Get set of unique labels
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[np.where(labels == j)] = i
    return labels

labels = find_labels(mask)

# 2.3: Filter out labeled outliers based on size
## From Live Code 
sizes = np.bincount(labels.ravel())

for i in range(1, np.amax(labels)+1):
    where = np.where(labels == i)
    if sizes[i] < 100:
        labels[where] = 0


def filter_by_size(labels, minsize, maxsize):
    # Find label sizes
    sizes = np.bincount(labels.ravel())
    # Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        # If the number of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find all pixels for label
            where = np.where(labels == i)
            labels[where] = 0
    # Get set of unique labels
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[np.where(labels == j)] = i
    return labels

# 3: Find mean signal for each nucleus from PCNA and nascentRNA

results = []
## Loop through each gene and field to process the images
for gene in genes:
    for field in fields:
        image = images[f"{gene}_{field}"]
        dapi_mask = mask
        labels = find_labels(dapi_mask)
        labels = filter_by_size(labels, minsize=100, maxsize=10000)
        
        # Create empty lists to to store the mean signal of nascentRNA, PCNA, and the ratios between them
        mean_nascentRNA = []
        mean_pcna = []
        ratios = []

        # Loop through each label (representing a nucleus) to calculate the mean signal
        for i in range(1, np.amax(labels) + 1):
            where = np.where(labels == i)
            mean_nascentRNA.append(np.mean(image[where[0], where[1], 1]))
            mean_pcna.append(np.mean(image[where[0], where[1], 2]))
        # Calculate the ratio of nascentRNA to PCNA, avoiding division by zero
        for nRNA, pcna in zip(mean_nascentRNA, mean_pcna):
            if pcna != 0:
                ratios.append(np.log2(nRNA / pcna))
            else:
                ratios.append(0)
        # Append the results for the current gene to the results list
        for mnRNA, mPCNA, ratio in zip(mean_nascentRNA, mean_pcna, ratios):
            results.append([gene, mnRNA, mPCNA, ratio])

# Save the results as a text file
output_filename = "week_10_signals.txt"
with open(output_filename, 'w') as file:
    file.write("Gene\tMean_NascentRNA\tMean_PCNA\tRatio\n")
    # Loop through the results and write each line to the file
    for result in results:
        gene = result[0]
        mean_nascentRNA = result[1]
        mean_pcna = result[2]
        ratio = result[3]
        file.write(f"{gene}\t{mean_nascentRNA:.4f}\t{mean_pcna:.4f}\t{ratio:.4f}\n")
