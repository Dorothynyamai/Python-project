import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import numexpr as ne
from joblib import Parallel, delayed

# Set the number of threads (cores) you want to use
num_threads = 8  # You can adjust this value based on the number of available CPU cores
ne.set_num_threads(num_threads)

# Load the gene expression data
# Replace 'path/to/gene_expression_data.csv' with the actual file path.
# gene_expression_data = pd.read_csv('path/to/gene_expression_data.csv', index_col=0)


# Transpose the data to have samples as rows instead of columns
gene_expression_data_transposed = gene_expression_data.T

# Calculate the distance matrix using parallel computation with numexpr
distance_matrix = hierarchy.distance.pdist(gene_expression_data_transposed.values, metric='euclidean')

# Perform hierarchical clustering
linkage_matrix = hierarchy.linkage(distance_matrix, method='ward')

# Create a large figure to accommodate the sample names
fig, ax = plt.subplots(figsize=(10, 8))  # You can adjust the figure size as needed


# Create the dendrogram
dendrogram = hierarchy.dendrogram(linkage_matrix, labels=gene_expression_data_transposed.index, orientation='left')

# Set the font size of the labels
plt.tick_params(axis='y', labelsize=10)  # Adjust the font size as needed

# Remove the grid around the figure
plt.grid(False)


plt.xlabel('Distance')
plt.ylabel('Sample')
plt.title('Sample Clustering Dendrogram')

# Adjust the layout to prevent label overlapping
plt.tight_layout()

# Save the plot as a PNG file
# Replace 'path/to/sample_dendrogram.png' with the actual file path.
plt.savefig('path/to/sample_dendrogram.png')
