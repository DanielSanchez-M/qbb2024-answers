#!/usr/bin/env python

# Import libraries
import numpy
import scipy
import matplotlib.pyplot as plt #import plotting package in python without having to bounce to R
import imageio
import plotly.express as px
import plotly

# Exercise 1: Loading the image data

## Observe a test image to obtain dimensions for image_Array
sample_image = imageio.v3.imread("APEX1_field0_DAPI.tif")
print(sample_image.shape)

## Make an empty list to contain the channel images
images = []

## Creating as gene list for loop
genenames = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
print(genenames)
## Creating a fields list for loop
fields = ["field0","field1"]
print(fields)
## Creating a channels list for loop
channels = ["DAPI", "PCNA", "nascentRNA"]
print(channels)

## Creating the for loop
for gene in genenames:
    for field in fields:
        image_array = numpy.zeros((sample_image.shape[0], sample_image.shape[1], 3), numpy.uint16)
        for i, channel in enumerate(channels):
            image_array[:, :, i] = imageio.v3.imread(f"{gene}_{field}_{channel}.tif")
        images.append(image_array)

plt.imshow(image_array)
plt.show()

# Exercise 2: Identifying individual cells

mask = []
for i in range(len(images[:])):
    mask.append(images[i][:, :, 0] >= numpy.mean(images[i][:, :, 0]))

label_array = []

def find_labels(mask):
    # Initialize label counter
    l = 0
    # Create array to hold the label values
    labels = numpy.zeros(mask.shape, numpy.int32)
    # Create list to track label equivalence
    equivalence = [0]
    
    # Check if the upper-left corner is part of an object (non-zero value)
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    
    # Loop through the first row to label connected components
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If the previous pixel has a label, use the same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # If no label, assign a new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    
    # Loop through each row (except the first row)
    for x in range(1, mask.shape[0]):
        # Check the left-most pixel in the row and the adjacent pixels (up and up-right)
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If the pixel above has a label, use the same label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If the pixel diagonally above-right has a label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # If no adjacent labeled pixel, assign a new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        
        # Loop through the rest of the row (except the last column)
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If the pixel above has a label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If the pixel diagonally above-right has a label, update equivalence
                    if mask[x - 1, y - 1]:
                        # If the pixel diagonally above-left has a label, relabel all connected components
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If the left pixel has a label, relabel accordingly
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither the left nor diagonal pixels have labels, use the diagonal equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If the pixel diagonally above-left has a label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If the left pixel has a label, use that label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, assign a new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        
        # Check the last pixel in the row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # If the pixel above has a label, use that label
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # If the pixel diagonally above-left has a label, use that label
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # If the left pixel has a label, use that label
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, assign a new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    
    # Convert equivalence list to a numpy array for efficient indexing
    equivalence = numpy.array(equivalence)
    
    # Go backwards through all equivalence labels to normalize labels
    for i in range(1, len(equivalence))[::-1]:
        # Relabel all pixels with a label from the equivalence list
        labels[numpy.where(labels == i)] = equivalence[i]
    
    # Get unique labels
    ulabels = numpy.unique(labels)
    
    # Normalize label values so they span from 1 to the number of labels
    for i, j in enumerate(ulabels):
        labels[numpy.where(labels == j)] = i
    
    return labels

# Function to filter labels by size
def filter_by_size(labels, minsize, maxsize):
    # Count the size of each label (connected component) in the label array
    sizes = np.bincount(labels.ravel())
    
    # Loop through each label (starting from 1 to skip the background label)
    for i in range(1, sizes.shape[0]):
        # Check if the size of the label is smaller than the minimum size or larger than the maximum size
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find the locations where the label appears in the label array
            where = np.where(labels == i)
            # Set the pixels of this label to 0 (removes the label from the output)
            labels[where] = 0
    
    # Get the unique remaining labels after filtering
    ulabels = np.unique(labels)
    
    # Normalize label values so they span from 1 to the number of labels
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    
    # Return the filtered label array
    return labels
plt.imshow(lab_array[0])
plt.show()


# Exercise 3: Score the PCNA and nascent RNA signal in each nucleus and plot them
## Step 3.1 Find the mean signal for each nucleus from the PCNA and nascent RNA channels



### Question: What do each of these values mean, based on the descriptions of what is being labeled?
### Answer:

## Step 3.2 Plot each set of data in a separate violin plot



### Question: What do each of these values mean, based on the descriptions of what is being labeled?
### Answer: 