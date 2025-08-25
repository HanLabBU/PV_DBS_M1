#%%
import h5py
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from scipy.stats import trim_mean
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import zscore
import scipy.stats as stats

from matplotlib import lines

import pickle
from cliffs_delta import cliffs_delta

import pydeseq2 as pds

from pydeseq2 import dds
from pydeseq2 import ds
import csv

import time
import math
import os

import io
import contextlib
import re

#%% loading market gene data -- takes awhile to filter (few min)
f_sep = os.sep

# this file data is actually in the folder: "/home/pierfier/Projects/Pierre Fabris/PV DBS neocortex/Seq_Data/"
data_dir = "/home/pierfier/Projects/Pierre Fabris/PV DBS neocortex/Seq_Data/"
hdf5_file = r"expression_matrix.hdf5"
metadata_file = r"metadata.csv"

#savepath = r"/ad/eng/research/eng_research_handata/Pierre Fabris/GenePVPlots/"
savepath = "/home/pierfier/Dropbox/RKC-HanLab/Pierre PV DBS Project Dropbox/Materials/Plots/Gene_Exp/"

plt.rcParams['pdf.fonttype'] = 42
# Removed genes because they are zero "Kcna3", "Chrna1", , "Htr3a", "Best1"
marker_genes = [
    "Pvalb",
    "Kcnc1", "Kcnc2", "Kcnc3",
    "Kcnab1", "Kcnab2",
    "Kcna1", "Kcna2",
    "Kcna6",
    "Hcn1", "Hcn2", 
    "Scn1a", "Scn2a", "Scn9a",
    "Cacna1c", "Cacna1d", "Cacna1a", "Cacna1b",
    "Grin1", "Gria1", "Grik1", "Atp1a1",
    "Kcnma1"
]

# Percentage amount to remove for trimmed
pp_cut = 0.1
cm = 1/2.54
log_scale = np.array([2, 9])

#%% filter data/counts for marker genes iteratively 
print("Filtering data/counts for marker genes...")
marker_gene_data = []  

with h5py.File(data_dir + hdf5_file, "r") as f:
    gene_names = f["data/gene"][:].astype(str)
    
    for gene in marker_genes:
        if gene in gene_names:
            gene_index = np.where(gene_names == gene)[0][0]  
            print(f"Extracting data for gene: {gene} (index {gene_index})")
            gene_data = f["data/counts"][gene_index, :]  
            marker_gene_data.append(gene_data)  
        else:
            print(f"Gene {gene} not found in data/gene.")

marker_gene_data = np.vstack(marker_gene_data)
print(f"Marker gene data extracted. Shape: {marker_gene_data.shape}")

#%% Perform single cell normalization, by dividing total number of gene expression
with h5py.File(data_dir + hdf5_file, 'r') as hdf:
    print(hdf['data'].keys())
    # Get the gene names dataset object:
    dset = hdf['data']['gene']
    print (dset.shape, dset.dtype)
        
    gene_names = [x.decode('utf-8') for x in hdf['data']['gene']]
    print(len(gene_names))
    
    cell_names = hdf['data']['samples'][:].astype('U')    
    print(hdf['data']['samples'].shape)
    
    gene_intersect = np.intersect1d(marker_genes,gene_names,return_indices = True)
    gene_inds = np.array(sorted(gene_intersect[2]))
    matrix_df = [hdf['data']['counts'][x,:] for x in gene_inds]
    
    
    data_size = (hdf['data']['counts'].shape)
    print(f"Coutns shape: {data_size}")
    
    chunk_size = 10000
    array_len = data_size[1]
    cell_counts_total = np.empty(array_len)
    n_iters = math.ceil(data_size[1]/chunk_size)
    print(f"Num iterations {n_iters}")
    
    # start_time = time.time()
    # iter_15 = np.array(hdf['data']['counts'][:,30000:40000]).sum(axis=0)
    # print("--- %s seconds ---" % (time.time() - start_time))
    print('Starting chunk extraction')
    for cnt in range(n_iters):
        print('get slice #',str(cnt))

        start_time = time.time()
        if cnt == n_iters:
            cell_counts_total[cnt*chunk_size:array_len] = np.array(hdf['data']['counts'][:,cnt*chunk_size:array_len]).sum(axis=0) 
            break
        else:
            cell_counts_total[cnt*chunk_size:(cnt+1)*chunk_size] = np.array(hdf['data']['counts'][:,cnt*chunk_size:(cnt+1)*chunk_size]).sum(axis=0)
        
        print("--- %s seconds ---" % (time.time() - start_time))
        
    hdf.close() 


# rest of code is fast

#%% Calculate the trimmed means for the shared clusters between visual and motor cortex 
#V1_regs = ['VISal', 'VISl', 'VISli', 'VISam', 'VISpm', 'VISp', 'VISpl', 'VISpor']

# Vocab definitions
# clusters: distinct groups of cells that generally belong to a cell type, like PV cells (generally by cell type)
# samples: cell identification

# filter metadata for motor and visual cortex, and shared cluster labels
print("Filtering metadata for motor and visual cortex...")
metadata_df = pd.read_csv(data_dir + metadata_file)
metadata_filtered = metadata_df[
    ((metadata_df["region_label"] == "MOp") | (metadata_df["region_label"] == "VISp") ) &
    metadata_df["cluster_label"].str.contains("Pvalb", na=False)
][["sample_name", "region_label", "cluster_label"]]

# Remove stray brain region
metadata_filtered = metadata_filtered[metadata_filtered['region_label'] != 'SSs-GU-VISC-AIp']

print(f"Filtered metadata shape: {metadata_filtered.shape}")

# match metadata_filtered with data/sample
print("Mapping metadata_filtered to data/sample...")
with h5py.File(data_dir + hdf5_file, "r") as f:
    samples = f["data/samples"][:].astype(str)
sample_to_index = {sample: idx for idx, sample in enumerate(samples)}


# find indices of samples in data/samples that match metadata_filtered
metadata_filtered["sample_index"] = metadata_filtered["sample_name"].map(sample_to_index)
metadata_filtered = metadata_filtered.dropna(subset=["sample_index"])
metadata_filtered["sample_index"] = metadata_filtered["sample_index"].astype(int)
print(f"Number of valid samples after matching: {len(metadata_filtered)}")

# separate metadata into visual and motor cortex
visual_metadata = metadata_filtered[metadata_filtered["region_label"].str.contains("VIS", na=False)]
motor_metadata = metadata_filtered[metadata_filtered["region_label"].str.contains("MOp", na=False)]
print(f"Visual cortex samples: {len(visual_metadata)}")
print(f"Motor cortex samples: {len(motor_metadata)}")

shared_clusters = set(visual_metadata["cluster_label"]).intersection(set(motor_metadata["cluster_label"]))
visual_metadata = visual_metadata[visual_metadata["cluster_label"].isin(shared_clusters)]
motor_metadata = motor_metadata[motor_metadata["cluster_label"].isin(shared_clusters)]
print(f"Shared clusters: {shared_clusters}")

# use sample indices to filter marker gene data for each region
print("Filtering marker gene data for visual and motor cortex...")

# Sample here refers to cell, so here we are indexing all the cells within the corresponding region
# The metadata for each region is determined as the shared PV clusters between motor and visual cortex
visual_sample_indices = visual_metadata["sample_index"].values
motor_sample_indices = motor_metadata["sample_index"].values
visual_marker_counts = marker_gene_data[:, visual_sample_indices]   
motor_marker_counts = marker_gene_data[:, motor_sample_indices]

visual_total_counts = cell_counts_total[visual_sample_indices].reshape(1, -1)
motor_total_counts = cell_counts_total[motor_sample_indices].reshape(1, -1)

#%% Manually calculating the trimmed means values and storig them into a dataframe

# Trim the expression values for each group
# Uses the global proportion_cutoff value
def trim_group(group):
    # Remove the zeros from each group
    group = group.sort_values(by='Exp')
    cutoff_num = int(pp_cut * group.shape[0])

    # If the cutoff is 0 just return the whole group
    if cutoff_num == 0:
        return group

    # Splice the group from the front and back
    trimmed_data = group.iloc[cutoff_num:-cutoff_num]

    #DEBUG
    if (group['gene'].unique() == 'Scn9a') & (group['clust'].unique() == '117_Pvalb'):
        print(trimmed_data)

    # If the cutoff removes elements from 
    if len(trimmed_data) == 0:
        return group
    
    # Return the trimmed array
    return trimmed_data


motor_dict = {}
visual_dict = {}

# Intialize empty dataframes
motor_dict['gene_exp_vals'] = pd.DataFrame()
motor_dict['trimmed_exp_vals'] = pd.DataFrame()
motor_dict['counts_total'] = pd.DataFrame()
motor_dict['trimmed_means'] = pd.DataFrame()

visual_dict['gene_exp_vals'] = pd.DataFrame()
visual_dict['trimmed_exp_vals'] = pd.DataFrame()
visual_dict['counts_total'] = pd.DataFrame()
visual_dict['trimmed_means'] = pd.DataFrame()

# Loop through each shared cluser and save the expression data and calculate the trimmed means
for cluster in shared_clusters:
    # Grab the marker gene expression values for this cluster
    motor_cluster_data_idx = motor_metadata[motor_metadata['cluster_label'] == cluster]['sample_index'].values
    motor_cluster_gene_data = marker_gene_data[:, motor_cluster_data_idx]

    # Grab expression data for all genes of interest (gois) and map them to the cell index values
    temp_df = pd.DataFrame(motor_cluster_gene_data, index=marker_genes, columns=motor_cluster_data_idx)
    temp_df = temp_df.melt(ignore_index=False, var_name='cell_idx', value_name='Exp')
    temp_df = temp_df.rename_axis('gene').reset_index()
    temp_df['clust'] = cluster
    temp_df['clust'] = temp_df['clust'].astype('category')

    # Put everything into a dictionary
    motor_dict['gene_exp_vals'] = pd.concat([motor_dict['gene_exp_vals'], temp_df])
    
    # Calculate the values to remove for trimmed means
    temp_df = temp_df.groupby('gene')[temp_df.columns].apply(trim_group).reset_index(drop=True)

    #DEBUG
    #if cluster == '108_Pvalb':
    #    print(temp_df[temp_df['gene'] == 'Scn9a'])
    #    raise Exception('Creating the trimmed means values')

    # Create NaNs for expressions that are at 0
    temp_df[temp_df == 0] = np.nan
    
    #-- Normalize each expression value by the total for each cell
    
    # Get all of the total counts for the corresponding cells from the trimmed mean
    motor_clust_counts_total = cell_counts_total[motor_cluster_data_idx].reshape(1, -1)
    total_tmp_df = pd.DataFrame(motor_clust_counts_total, columns=motor_cluster_data_idx)
    total_tmp_df = total_tmp_df.melt(var_name = 'cell_idx', value_name='counts_total')

    #total_tmp_df['clust'] = cluster
    #total_tmp_df['clust'] = total_tmp_df['clust'].astype('category')

    motor_dict['counts_total'] = pd.concat([motor_dict['counts_total'], total_tmp_df])

    # Match the expression values with the cell counts total
    norm_df = temp_df.merge(total_tmp_df, on='cell_idx')

    # Normalize each expression value
    norm_df['norm'] = 100000*(norm_df['Exp']/norm_df['counts_total']) + 1

    # Add trimmed data with normalization to dictionary
    motor_dict['trimmed_exp_vals'] = pd.concat([motor_dict['trimmed_exp_vals'], norm_df])

    # Compute the trimmed means for all genes
    trim_df = norm_df.groupby('gene')['norm'].apply(np.nanmean).reset_index()
    trim_df['clust'] = cluster
    trim_df['clust'] = trim_df['clust'].astype('category')

    motor_dict['trimmed_means'] = pd.concat([motor_dict['trimmed_means'], trim_df])

    #-- Same thing but for visual PVs
    # Grab the marker gene expression values for this cluster
    visual_cluster_data_idx = visual_metadata[visual_metadata['cluster_label'] == cluster]['sample_index'].values
    visual_cluster_gene_data = marker_gene_data[:, visual_cluster_data_idx]

    # Grab expression data for all genes of interest (gois) and map them to the cell index values
    temp_df = pd.DataFrame(visual_cluster_gene_data, index=marker_genes, columns=visual_cluster_data_idx)
    temp_df = temp_df.melt(ignore_index=False, var_name='cell_idx', value_name='Exp')
    temp_df = temp_df.rename_axis('gene').reset_index()
    temp_df['clust'] = cluster
    temp_df['clust'] = temp_df['clust'].astype('category')

    # Put everything into a dictionary
    visual_dict['gene_exp_vals'] = pd.concat([visual_dict['gene_exp_vals'], temp_df])
    
    # Calculate the values to remove for trimmed means
    temp_df = temp_df.sort_values(by=['gene', 'Exp'])
    temp_df = temp_df.groupby('gene')[temp_df.columns].apply(trim_group).reset_index(drop=True)

    # Create NaNs for expressions that are at 0
    temp_df[temp_df == 0] = np.nan
    
    #-- Normalize each expression value by the total for each cell
    
    # Get all of the total counts for the corresponding cells from the trimmed mean
    visual_clust_counts_total = cell_counts_total[visual_cluster_data_idx].reshape(1, -1)
    total_tmp_df = pd.DataFrame(visual_clust_counts_total, columns=visual_cluster_data_idx)
    total_tmp_df = total_tmp_df.melt(var_name = 'cell_idx', value_name='counts_total')

    #total_tmp_df['clust'] = cluster
    #total_tmp_df['clust'] = total_tmp_df['clust'].astype('category')

    visual_dict['counts_total'] = pd.concat([visual_dict['counts_total'], total_tmp_df])

    # Match the expression values with the cell counts total
    norm_df = temp_df.merge(total_tmp_df, on='cell_idx')

    # Normalize each expression value
    norm_df['norm'] = 100000*(norm_df['Exp']/norm_df['counts_total']) + 1

    # Add trimmed data with normalization to dictionary
    visual_dict['trimmed_exp_vals'] = pd.concat([visual_dict['trimmed_exp_vals'], norm_df])

    # Compute the trimmed means for all genes
    trim_df = norm_df.groupby('gene')['norm'].apply(np.nanmean).reset_index()
    trim_df['clust'] = cluster
    trim_df['clust'] = trim_df['clust'].astype('category')

    visual_dict['trimmed_means'] = pd.concat([visual_dict['trimmed_means'], trim_df])


# %% Perform min-max normalization for each gene

# Perform the (x - min)/ (max - min) for each gene across both regions
m1_trimmed_means_df = motor_dict['trimmed_means'].copy()
v1_trimmed_means_df = visual_dict['trimmed_means'].copy()

# Reorganize the order of genes

max_min_df = pd.DataFrame()
for gene in order_genes:
    # Read in the expression value for both motor and visual brain region
    m1_gene_means = m1_trimmed_means_df[m1_trimmed_means_df['gene'] == gene]['norm']
    v1_gene_means = v1_trimmed_means_df[v1_trimmed_means_df['gene'] == gene]['norm']
    
    #DEBUG
    if gene == 'Scn9a':
        print(m1_gene_means)

    m1_gene_means.index = m1_trimmed_means_df[m1_trimmed_means_df['gene'] == gene]['clust']
    v1_gene_means.index = v1_trimmed_means_df[v1_trimmed_means_df['gene'] == gene]['clust']

    m1_gene_means.name = 'm1'
    v1_gene_means.name = 'v1'

    # Combine region values for single gene to perform total gene max-min normalization
    #TODO could just switch the order and have V1 then M1
    comb_df = pd.concat([m1_gene_means, v1_gene_means], axis=1)      

    min_val = np.nanmin(comb_df.values)
    max_val = np.nanmax(comb_df.values)

    comb_df = (comb_df - min_val) / (max_val - min_val)

    #DEBUG
    #if comb_df.isna().any().any():
    #    raise Exception()

    # Add the gene as the higher level for the Multi-Index
    reg_col = comb_df.columns
    mult_col = pd.MultiIndex.from_tuples([(gene, col) for col in reg_col])
    comb_df.columns = mult_col

    # Concatenate the values horizontally for both gene data
    max_min_df = pd.concat([max_min_df, comb_df], axis=1)
    
plt.figure(figsize=(8, 18))
sns.heatmap(max_min_df.T, cmap=sns.diverging_palette(250, 10, as_cmap=True), annot=False, vmin=0, vmax=1, cbar_kws={'label': 'Trimmed Mean'})
plt.hlines(y=np.arange(2, 50, 2), xmin=0, xmax=26, color='k')
plt.title("Individual Gene Range")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + 'Max_Min_Norm_M1_V1_per_gene.png', format='png')
plt.savefig(savepath + 'Max_Min_Norm_M1_V1_per_gene.pdf', format='pdf')
plt.show()

# Average across clusters
#avg_clust_gene = max_min_df.mean(axis=0)
#
#plt.figure(figsize=(18, 6))
#sns.heatmap(pd.DataFrame(avg_clust_gene).T, cmap=sns.diverging_palette(250, 10, as_cmap=True), annot=False, vmin=0, vmax=1, cbar_kws={'label': 'Trimmed Mean'})
#plt.vlines(x=np.arange(2, 50, 2), ymin=0, ymax=26, color='k')
#plt.title("Motor then Visual")
#plt.xlabel("Genes")
#plt.ylabel("Cluster Labels")
#plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
#plt.savefig(savepath + 'Max_Min_Norm_M1_V1_per_gene.png', format='png')
plt.show()

# %% Test the significance of each cluster between M1 and V1 for each gene
m1_trimmed_vals_df = motor_dict['trimmed_exp_vals'].copy()
v1_trimmed_vals_df = visual_dict['trimmed_exp_vals'].copy()

# Pre-emptively remove NaNs
m1_trimmed_vals_df = m1_trimmed_vals_df[~np.isnan(m1_trimmed_vals_df['norm']) ]
v1_trimmed_vals_df = v1_trimmed_vals_df[~np.isnan(v1_trimmed_vals_df['norm']) ]

m1_trimmed_vals_df['reg'] = 'm1'
m1_trimmed_vals_df['reg'] = m1_trimmed_vals_df['reg'].astype('category')
v1_trimmed_vals_df['reg'] = 'v1'
v1_trimmed_vals_df['reg'] = v1_trimmed_vals_df['reg'].astype('category')

# Concatenate the two brain region dataframes
consol_reg_df = pd.concat([m1_trimmed_vals_df, v1_trimmed_vals_df])

# Function to calculate the p-value between 
def calc_pval(gene_df):
    # Get V1 values
    v1_vals = gene_df[gene_df['reg'] == 'v1']['norm'].values
    m1_vals = gene_df[gene_df['reg'] == 'm1']['norm'].values
    
    # Remove NaNs
    #v1_vals = v1_vals[~np.isnan(v1_vals)]
    #m1_vals = m1_vals[~np.isnan(m1_vals)]

    # Replace empty values with 0
    if v1_vals.size == 0:
        v1_vals = np.array([0])
        print('Empty V1')    
    if m1_vals.size == 0:
        m1_vals = np.array([0])
        print('Empty M1')

    # Calculate the Mann Whitney between the brain regions
    u_stats, p_value = stats.mannwhitneyu(m1_vals, v1_vals, alternative='two-sided')
    
    dir = 0
    # Compare if V1 or M1 is higher
    if np.median(v1_vals) > np.median(m1_vals):
        dir = 1
    elif np.median(v1_vals) < np.median(m1_vals):
        dir = -1
    #print(p_value)

    return p_value, dir

stats_df = pd.DataFrame()
for gene in order_genes:
    gene_reg_df = consol_reg_df[consol_reg_df['gene'] == gene]

    # TODO delete this
    sam = gene_reg_df[gene_reg_df['clust'] == '108_Pvalb']
    

    # Group by cluster to determine how many are significant between the two regions
    clust_df = gene_reg_df.groupby('clust').apply(lambda x: calc_pval(x))
    clust_df = clust_df.to_frame()

    clust_df.columns = [gene]
    #print(clust_df)

    stats_df = pd.concat([stats_df, clust_df], axis=1)

# Output csv file
stats_df = stats_df.T
p_val_df = stats_df.map(lambda x: x[0]) 

p_val_df.to_csv(savepath + 'V1_M1_gene_p_values.csv')

# Determine for each cluster if V1 or M1 was higher
sig_diff_df = stats_df.map(lambda x: x[1] if x[0] < 0.05 else 0)
sig_diff_df = sig_diff_df.T

# Count how many clusters for each gene that are V1 or M1 significantly higher
m1_sum = sig_diff_df.apply(lambda col: (col == -1).sum())
v1_sum = sig_diff_df.apply(lambda col: (col == 1).sum())

m1_sum.name = 'm1'
v1_sum.name = 'v1'

# Significant clusters count
clust_count = pd.concat([m1_sum, v1_sum], axis=1)

clust_count.to_csv(savepath + 'V1_M1_cluster_count.csv')

# %% Expression difference for each PV cluster

# Flag to determine how to show heatmap
#opts = 'fold_change'
opts = 'scale_diff'

m1_trimmed_means_df = motor_dict['trimmed_means'].copy()
v1_trimmed_means_df = visual_dict['trimmed_means'].copy()

reg_diff_df = pd.DataFrame()

for gene in marker_genes:
    # Read in the expression value for both motor and visual brain region
    m1_gene_means = m1_trimmed_means_df[m1_trimmed_means_df['gene'] == gene]['norm']
    v1_gene_means = v1_trimmed_means_df[v1_trimmed_means_df['gene'] == gene]['norm']
    
    m1_gene_means.index = m1_trimmed_means_df[m1_trimmed_means_df['gene'] == gene]['clust']
    v1_gene_means.index = v1_trimmed_means_df[v1_trimmed_means_df['gene'] == gene]['clust']

    m1_gene_means.name = 'm1'
    v1_gene_means.name = 'v1'

    # Difference of V1 to M1
    if opts == 'scale_diff':
        diff_s = v1_gene_means.sub(m1_gene_means, fill_value=0).div(v1_gene_means.add(m1_gene_means, fill_value=0), fill_value=1)

        diff_df = diff_s.to_frame()
        diff_df.columns = [gene]

        reg_diff_df = pd.concat([reg_diff_df, diff_df], axis=1)
    elif opts == 'fold_change':
        # Calculating log2 (V1/M1) fold change
        fc_df = np.log2(v1_gene_means.div(m1_gene_means, fill_value=1))
        fc_df = fc_df.to_frame()
        fc_df.columns = [gene]
        reg_diff_df = pd.concat([reg_diff_df, fc_df], axis=1)

#plt.rcParams.update({'font.size': 4})
#plt.figure(figsize=(5*cm, 8.8*cm))
plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(8, 9))

# Only plot the significant genes
reg_diff_df = reg_diff_df.T
reg_diff_df = reg_diff_df.loc[sorted_genes]

max_chan = np.max([np.abs(np.min(reg_diff_df)), np.abs(np.max(reg_diff_df))])

sns.heatmap(reg_diff_df, cmap=sns.diverging_palette(250, 10, as_cmap=True), annot=False,\
    vmin=-max_chan, vmax=max_chan, cbar_kws={'label': opts})
plt.title("Gene Cluster difference")
plt.xlabel("Cluster Labels")
plt.ylabel("Genes")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + opts + '_Norm_M1_V1_sig_diff_per_cluster.png', format='png')
plt.savefig(savepath + opts + '_Norm_M1_V1_sig_diff_per_cluster.pdf', format='pdf')
plt.show()

# %%
# Gene and cluster expression values for both regions (grouped by cluster labels)

m1_trimmed_means_df = motor_dict['trimmed_means'].copy()
v1_trimmed_means_df = visual_dict['trimmed_means'].copy()

trim_exp_df = pd.DataFrame()

for cluster in m1_trimmed_means_df['clust'].unique():
    # Grab all of the norm gene data from the cluster
    m1_gene_counts = m1_trimmed_means_df[m1_trimmed_means_df['clust'] == cluster]['norm']
    m1_gene_counts.index = m1_trimmed_means_df[m1_trimmed_means_df['clust'] == cluster]['gene']
    v1_gene_counts = v1_trimmed_means_df[v1_trimmed_means_df['clust'] == cluster]['norm']
    v1_gene_counts.index = v1_trimmed_means_df[v1_trimmed_means_df['clust'] == cluster]['gene']

    # Sort the genes by array order
    m1_gene_counts = m1_gene_counts.loc[sorted_bydiff_genes]
    v1_gene_counts = v1_gene_counts.loc[sorted_bydiff_genes]

    m1_gene_counts.name = 'm1'
    v1_gene_counts.name = 'v1'

    # Combine the two series
    comb_df = pd.concat([v1_gene_counts, m1_gene_counts], axis=1)
    reg_col = comb_df.columns
    mult_col = pd.MultiIndex.from_tuples([(cluster, col) for col in reg_col])
    comb_df.columns = mult_col

    # Concatenate the values horizontally for both gene data
    trim_exp_df = pd.concat([trim_exp_df, comb_df], axis=1)

# Re-organize the genes order by the sorted difference between V1 and M1
trim_exp_df = trim_exp_df.loc[sorted_bydiff_genes]

# Plot the ordered heatmap
plt.rcParams.update({'font.size': 6})
plt.figure(figsize=(15*cm, 15*cm))
#plt.rcParams.update({'font.size': 12})
#plt.figure(figsize=(7, 12))

# Transform and log the expression values
log_exp = np.log2(trim_exp_df)

ax = sns.heatmap(log_exp, cmap=sns.color_palette('viridis', as_cmap=True), annot=False, vmin=np.min(log_scale),\
    vmax=np.max(log_scale), cbar_kws={'label': 'Log2 Trimmed Mean'})

# Store the colorbar axis values
cbar_clust_exp = ax.collections[0].colorbar

#plt.hlines(y=np.arange(2, 50, 2), xmin=0, xmax=26, color='k')

# Lines that separate between each cluster
plt.vlines(x=np.arange(2, 50, 2), ymin=0, ymax=23, linewidth=0.5, color='k')
max_y = log_exp.shape[0]
min_y = 0

# Plot dashed vertical lines
for v_x in np.arange(1, 50, 2):
    line = lines.Line2D([v_x, v_x], [min_y, max_y], color='k', linewidth=0.25,\
                        linestyle='--', dashes=(10, 5))
    plt.gca().add_line(line)

plt.title("Gene Expression")
plt.xlabel("PV Subclusters")
plt.ylabel("Genes")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + 'Norm_M1_V1_per_cluster_by_cluster.png', format='png')
plt.savefig(savepath + 'Norm_M1_V1_per_cluster_by_cluster.pdf', format='pdf')
plt.show()

# %%
# Gene with brain region median (regions are put side-by-side)

# Flag to determine how to show heatmap
#opts = 'fold_change'
opts = 'scale_diff'

# Compute the averages of each gene for both brain regions
med_motor_exp = motor_dict['trimmed_exp_vals'].groupby('gene')['norm'].apply(np.nanmedian)
med_motor_exp.name = 'm1'

med_visual_exp = visual_dict['trimmed_exp_vals'].groupby('gene')['norm'].apply(np.nanmedian)
med_visual_exp.name = 'v1'

# Re-order genes
#med_motor_exp = med_motor_exp.loc[order_genes]
#med_visual_exp = med_visual_exp.loc[order_genes]

# Consolidate into one dataframe
med_reg_df = pd.concat([med_visual_exp, med_motor_exp], axis=1, keys=['V1', 'M1'])

# Plotting the difference in averages
if opts == 'scale_diff':
    med_reg_diff_df = (med_visual_exp - med_motor_exp) / (med_visual_exp + med_motor_exp)
elif opts == 'fold_change':
    med_reg_diff_df = np.log2(med_visual_exp.div(med_motor_exp, fill_value=1))

med_reg_diff_df = med_reg_diff_df.to_frame()

#-- Sort the genes by V1 value 
# Sort by significant genes and then non-significant genes
sorted_byval_genes = []

m1_high_genes = med_reg_diff_df[med_reg_diff_df[0] < 0].index
v1_high_genes = med_reg_diff_df[med_reg_diff_df[0] > 0].index

# M1 high significant genes
m1_sig_genes = list(set(m1_high_genes) & set(des_sig_genes))
m1_non_genes = list(set(m1_high_genes) - set(des_sig_genes))
# V1 high significant genes
v1_sig_genes = list(set(v1_high_genes) & set(des_sig_genes))
v1_non_genes = list(set(v1_high_genes) - set(des_sig_genes))

# Sort based on the V1 values
m1_sig_genes = med_reg_df.loc[m1_sig_genes].sort_values(by='V1')
v1_sig_genes = med_reg_df.loc[v1_sig_genes].sort_values(by='V1')

# Sort the non-significant genes individually
#m1_non_genes = med_reg_df.loc[m1_non_genes].sort_values(by='V1')
#v1_non_genes = med_reg_df.loc[v1_non_genes].sort_values(by='V1')
## Create the final long index list
#sorted_genes = pd.concat([m1_sig_genes, m1_non_genes, v1_non_genes, v1_sig_genes]).index

# Sort the non-significant genes togehter
non_genes = med_reg_df.loc[list(set(v1_non_genes) | set(m1_non_genes))].sort_values(by='V1')
sorted_byval_genes = pd.concat([m1_sig_genes, non_genes, v1_sig_genes]).index

# -- Done sorting by V1 value

# -- Sort by difference value
sorted_bydiff_genes = []

# M1 high significant genes
m1_sig_genes = list(set(m1_high_genes) & set(des_sig_genes))
m1_non_genes = list(set(m1_high_genes) - set(des_sig_genes))
# V1 high significant genes
v1_sig_genes = list(set(v1_high_genes) & set(des_sig_genes))
v1_non_genes = list(set(v1_high_genes) - set(des_sig_genes))

# Sort based on the difference values
m1_sig_genes = med_reg_diff_df.loc[m1_sig_genes].sort_values(by=0)
v1_sig_genes = med_reg_diff_df.loc[v1_sig_genes].sort_values(by=0)

non_genes = med_reg_diff_df.loc[list(set(v1_non_genes) | set(m1_non_genes))].sort_values(by=0)
sorted_bydiff_genes = pd.concat([m1_sig_genes, non_genes, v1_sig_genes]).index

# -- Done sorting by diff value


# Re-organize by the sort of V1
#med_reg_df = med_reg_df.loc[sorted_byval_genes]
#med_reg_diff_df = med_reg_diff_df.loc[sorted_byval_genes]

# Re-organize by the sort of V1
med_reg_df = med_reg_df.loc[sorted_bydiff_genes]
med_reg_diff_df = med_reg_diff_df.loc[sorted_bydiff_genes]

# Plotting both brain region averages
plt.rcParams.update({'font.size': 6})
fig = plt.figure(figsize=(6*cm, 15*cm))

#max_chan = np.max([np.abs(np.min(med_reg_df)), np.abs(np.max(med_reg_df))])
log_exp = np.log2(med_reg_df)

ax = sns.heatmap(log_exp, cmap=sns.color_palette('viridis', as_cmap=True),\
    vmin=cbar_clust_exp.vmin, vmax=cbar_clust_exp.vmax, cbar_kws={'label': 'log2(expression)'})

# Set the colorbar the same as the per cluster
#cbar_reg = ax.collections[0].colorbar

# Save figure
plt.savefig(savepath + 'Median_Region.png', format='png')
plt.savefig(savepath + 'Median_Region.pdf', format='pdf')

plt.rcParams.update({'font.size': 6})
fig = plt.figure(figsize=(6*cm, 15*cm))

log_exp = med_reg_diff_df
max_chan = np.max([np.abs(np.min(log_exp)), np.abs(np.max(log_exp))])

sns.heatmap(log_exp, cmap=sns.diverging_palette(250, 10, as_cmap=True),\
    vmin=-max_chan, vmax=max_chan, cbar_kws={'label': opts})

# Save figure
plt.savefig(savepath + opts + '_Median_Region.png', format='png')
plt.savefig(savepath + opts + '_Median_Region.pdf', format='pdf')

plt.show()

# %%
# Calculate differential expression using DESeq2 package

# First consolidate all cells into a dataframe where rows are genes columns are cell names
motor_df = motor_dict['trimmed_exp_vals']
visual_df = visual_dict['trimmed_exp_vals']

motor_pv_counts = motor_df.pivot(index='gene', columns='cell_idx', values='Exp')
visual_pv_counts = visual_df.pivot(index='gene', columns='cell_idx', values='Exp')

# Concatenate to into one dataframe
gene_cell_df = pd.concat([motor_pv_counts, visual_pv_counts], axis=1).T
gene_cell_df[np.isnan(gene_cell_df)] = 0

# Construct col_data to correspond with cells
motor_meta = pd.DataFrame(['m1'] * len(motor_pv_counts.columns), index=motor_pv_counts.columns, columns=['Condition'])
visual_meta = pd.DataFrame(['v1'] * len(visual_pv_counts.columns), index=visual_pv_counts.columns, columns=['Condition'])

col_meta = pd.concat([motor_meta, visual_meta])

des = pds.dds.DeseqDataSet(adata=None, counts=gene_cell_df, metadata=col_meta, design='~ Condition')
des.deseq2()

#%% Save deseq object to pickle file
# Save the deseq output (It takes a really long time to run the analysis)
with open(data_dir + 'des_results.pkl', "wb") as f:
    pickle.dump(des.to_picklable_anndata(), f)

# Write des to pickle file

## This did not seem to work.....
## Save the deseq output (It takes a really long time to run the analysis)
#with open('des_results.pkl', "wb") as f:
#    pickle.dump(des, f)

# %% Read in pickle file
## This also misses some data needed for further analysis
with open('des_results.pkl', "rb") as f:
    red_data = pickle.load(f)
    #des = red_data

# %% Output the results of the differential expression dispersion and size factors
des_seq_stats = pds.ds.DeseqStats(des, contrast = ['Condition', 'm1', 'v1'])

des_seq_stats.summary()

# %% Save 'des_seq_stats'

# %% Output the dataframe to a csv file
#TODO continue execution here, need to find where sorted_genes is created
#des_seq_stats.results_df = des_seq_stats.results_df.loc[sorted_bydiff_genes]
des_seq_stats.results_df.to_csv(savepath + 'Deseq' +f_sep+ 'deseq_results.csv')

# This is the old way of outputting the results
# Display the deseq results into a file
#buffer = io.StringIO()
#with contextlib.redirect_stdout(buffer):
#    
#
#summary_str = buffer.getvalue()
#summary_str = re.sub(r'\s+', ' ', summary_str)
##summary_str = summary_str.replace(' ', ',')
#lines = summary_str.split('\n')
#
## Write to file
#with open('deseq_results.csv', 'w') as f:
#    writer = csv.writer(f)
#    for line in lines:
#        writer.writerow([line])


# Plot the histograms from the deseq counts and from their p-values
des_res_df = des_seq_stats.results_df
des_sig_genes = des_res_df[des_res_df['pvalue'] < 0.05].index

# Get the cell IDX's for each brain region
m1_idx = des.obs['Condition'] == 'm1'
v1_idx = des.obs['Condition'] == 'v1'

# Total cells for each brain region
m1_cell_nums = des.X[m1_idx, :].shape[0]
v1_cell_nums = des.X[v1_idx, :].shape[0]

# Plot all of the histograms for the significant genes
for gene in des_sig_genes:
    # Grab the m1 data for this gene
    gene_id = des.var_names.get_loc(gene)
    
    m1_counts = des.layers["normed_counts"][m1_idx, gene_id]
    v1_counts = des.layers["normed_counts"][v1_idx, gene_id]

    m1_counts = m1_counts[m1_counts != 0]
    v1_counts = v1_counts[v1_counts != 0]

    # Calculate the medians for each region
    m1_med_counts = np.nanmedian(m1_counts)
    v1_med_counts = np.nanmedian(v1_counts)

    # Determine which region has the higher median
    high_reg = ''
    if m1_med_counts > v1_med_counts:
        high_reg = 'm1'
    else:
        high_reg = 'v1'

    #m1_counts = m1_counts.astype(float)
    #v1_counts = v1_counts.astype(float)
#
    ## Take out the 0's in each array
    #m1_counts[m1_counts == 0] = np.nan
    #v1_counts[v1_counts == 0] = np.nan

    # Percentage of neurons expressing gene
    gene_cell_perc_v1 = 100*(v1_counts.shape[0]/v1_cell_nums)
    gene_cell_perc_m1 = 100*(m1_counts.shape[0]/m1_cell_nums)

    # Calculate the cliff's delta
    delta = cliffs_delta(v1_counts, m1_counts)

    # Store the p-value
    print(gene)
    if gene == 'Hcn2':
        print(des_res_df.loc[gene, 'pvalue'])
        raise Exception('Stop here')
    p_value = des_res_df.loc[gene, 'pvalue']

    # Calculate the histogram parts with the difference
    vis_hist, vis_bins = np.histogram(v1_counts, bins=100)
    mot_hist, mot_bins = np.histogram(m1_counts, bins=vis_bins) # Use the visual bins to keep it consistent
    vis_hist = vis_hist/v1_counts.shape[0]
    mot_hist = mot_hist/m1_counts.shape[0]

    # Calculate inter-quartile range
    q1 = np.percentile(np.concatenate((m1_counts, v1_counts)), 25)
    q3 = np.percentile(np.concatenate((m1_counts, v1_counts)), 75)

    iqr = q3 - q1
    lower = q1 - (1.5 * iqr)
    upper = q3 + (1.5 * iqr)

    # Use 0 as the lowest low component
    if lower < 0:
        lower = 0

    plt.rcParams.update({'font.size': 4})
    plt.figure(figsize=(3*cm, 2.2*cm))

    #plt.rcParams.update({'font.size': 6})
    #plt.figure(figsize=(10*cm, 8*cm))
    
    plt.bar(mot_bins[:-1], mot_hist, width=np.diff(mot_bins),\
            color='blue', alpha=0.5, label='M1 (' +str(round(gene_cell_perc_m1, 1))+ '%)')
    plt.bar(vis_bins[:-1], vis_hist, width=np.diff(vis_bins),\
            color='red', alpha=0.5, label='V1 (' +str(round(gene_cell_perc_v1, 1))+ '%)')

    # Set the cutoff points
    #plt.vlines([lower, upper], ymin=0, ymax=0.2, linestyle='-', color='k')
    plt.xlim(lower, upper)

    #plt.hist(mot_cpm_norm[0], bins=30, alpha=0.5, label='M1', color='blue')
    #plt.hist(vis_cpm_norm[0], bins=30, alpha=0.5, label='V1', color='red')
    
    plt.legend()
    plt.xlabel("Expression Value")
    plt.ylabel("Cell Percentage")
    plt.title(gene + '\ncliff\'s deltra: '\
                 + str(delta) + '\n p='+ str(round(p_value, 4)) +\
                '\nHigher Reg: ' + str(high_reg) )
        
    plt.savefig(savepath + 'Deseq' + f_sep + 'Gene' + f_sep + gene + '_exp_cell_perc.png', format='png')
    plt.savefig(savepath + 'Deseq' + f_sep + 'Gene' + f_sep + gene + '_exp_cell_perc.pdf', format='pdf')
    
plt.show()

# %% (Old way)
# Gene and cluster expression values for both regions (grouped by genes)
  
m1_trimmed_means_df = motor_dict['trimmed_means'].copy()
v1_trimmed_means_df = visual_dict['trimmed_means'].copy()

trim_exp_df = pd.DataFrame()
for gene in order_genes:
    # Read in the expression value for both motor and visual brain region
    m1_gene_means = m1_trimmed_means_df[m1_trimmed_means_df['gene'] == gene]['norm']
    v1_gene_means = v1_trimmed_means_df[v1_trimmed_means_df['gene'] == gene]['norm']

    # Maintain cluster indices
    m1_gene_means.index = m1_trimmed_means_df[m1_trimmed_means_df['gene'] == gene]['clust']
    v1_gene_means.index = v1_trimmed_means_df[v1_trimmed_means_df['gene'] == gene]['clust']

    m1_gene_means.name = 'm1'
    v1_gene_means.name = 'v1'

    comb_df = pd.concat([m1_gene_means, v1_gene_means], axis=1)

    # Add the gene as the higher level for the Multi-Index
    reg_col = comb_df.columns
    mult_col = pd.MultiIndex.from_tuples([(gene, col) for col in reg_col])
    comb_df.columns = mult_col

    # Concatenate the values horizontally for both gene data
    trim_exp_df = pd.concat([trim_exp_df, comb_df], axis=1)

plt.rcParams.update({'font.size': 6})
plt.figure(figsize=(8*cm, 15*cm))
#plt.rcParams.update({'font.size': 12})
#plt.figure(figsize=(7, 12))

# Transform and log the expression values
log_exp = np.log2(trim_exp_df.T)

sns.heatmap(log_exp, cmap=sns.color_palette('viridis', as_cmap=True), annot=False, vmin=np.min(log_exp),\
    vmax=np.max(log_exp), cbar_kws={'label': 'Log2 Trimmed Mean'})
plt.hlines(y=np.arange(2, 50, 2), xmin=0, xmax=26, color='k')
plt.hlines(y=np.arange(1, 50, 2), xmin=0, xmax=26, linewidth=0.5, color='k')
plt.title("Cluster Gene Expression")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + 'Norm_M1_V1_per_cluster_by_gene.png', format='png')
plt.savefig(savepath + 'Norm_M1_V1_per_cluster_by_gene.pdf', format='pdf')
plt.show()

# %% (old way)
# Gene with brain region average, both values and difference between each regions

# Flag to determine how to show heatmap
opts = 'fold_change'
#opts = 'scale_diff'

med_motor_exp = motor_dict['trimmed_exp_vals'].groupby('gene')['norm'].apply(np.nanmean)
avg_motor_exp.name = 'm1'

avg_visual_exp = visual_dict['trimmed_exp_vals'].groupby('gene')['norm'].apply(np.nanmean)
avg_visual_exp.name = 'v1'

# Sort all genes in a specified order

# Grab all of the Sodium, "Scn*" genes
scn_genes = [gene for gene in marker_genes if gene.startswith("Scn")]
scn_genes.append("Nav1")

# Potassium genes, "Kcn*"
kcn_genes = [gene for gene in marker_genes if gene.startswith("Kcn")]
# Calcium genes, "Cacn*" genes
cacn_genes = [gene for gene in marker_genes if gene.startswith("Cacn")]

# Order genes by average expression levels
scn_df = avg_visual_exp.loc[scn_genes].copy()
scn_df = scn_df.sort_values()
scn_genes = list(scn_df.index)

kcn_df = avg_visual_exp.loc[kcn_genes].copy()
kcn_df = kcn_df.sort_values()
kcn_genes = list(kcn_df.index)

cacn_df = avg_visual_exp.loc[cacn_genes].copy()
cacn_df = cacn_df.sort_values()
cacn_genes = list(cacn_df.index)

# Reorganize gene names
order_genes = ['Pvalb']
order_genes = order_genes + scn_genes + kcn_genes + cacn_genes
other_genes = [gene for gene in marker_genes if gene not in order_genes]

order_genes = order_genes + other_genes

avg_motor_exp = avg_motor_exp.loc[order_genes]
avg_visual_exp = avg_visual_exp.loc[order_genes]

# Check the expression of the calcium gene
print('Gene diff, v1 then m1')
print(avg_visual_exp.loc['Atp1a1'])
print(avg_motor_exp.loc['Atp1a1'])

# Combine V1 and M1 together with multi-index
# This way does not seem great
#cur_index = avg_motor_exp.index
#new_index = pd.MultiIndex.from_tuples([(idx, 'M1') for idx in cur_index])
#avg_motor_exp.index = new_index
#
#cur_index = avg_visual_exp.index
#new_index = pd.MultiIndex.from_tuples([(idx, 'V1') for idx in cur_index])
#avg_visual_exp.index = new_index

avg_reg_df = pd.concat([avg_visual_exp, avg_motor_exp], axis=0, keys=['V1', 'M1'])
avg_reg_df = avg_reg_df.swaplevel(0, 1)
avg_reg_df = avg_reg_df.loc[order_genes]
avg_reg_df = avg_reg_df.to_frame()

# Change the index level
#avg_reg_df = avg_reg_df.sort_index(kind='merge')
#avg_reg_df = avg_reg_df.swaplevel(0, 1).sort_index().to_frame()

if opts == 'scale_diff':
    avg_reg_diff_df = (avg_visual_exp - avg_motor_exp) / (avg_visual_exp + avg_motor_exp)
elif opts == 'fold_change':
    avg_reg_diff_df = np.log2(avg_visual_exp.div(avg_motor_exp, fill_value=1))

#print(v1_higher_genes)

plt.rcParams.update({'font.size': 6})
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6*cm, 15*cm))
#plt.rcParams.update({'font.size': 8})
#fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20, 8))

avg_exp = np.log2(avg_reg_df)
sns.heatmap(avg_exp, ax=axs[0], cmap=sns.color_palette('viridis', as_cmap=True), vmin=np.min(log_scale),\
    vmax=np.max(log_scale), cbar_kws={'label': ' Log2 Exp Values'})
axs[0].hlines(y=np.arange(2, 50, 2), xmin=0, xmax=26, color='k')
axs[0].hlines(y=np.arange(1, 50, 2), xmin=0, xmax=26, linewidth=0.5, color='k')
axs[0].set_title("Pooled Gene Expression")
axs[0].set_xlabel("Region")
axs[0].set_ylabel("Gene Labels")

avg_reg_diff_df = avg_reg_diff_df.to_frame()
avg_reg_diff_df = avg_reg_diff_df.loc[sig_genes]

max_chan = np.max([np.abs(np.min(avg_reg_diff_df)), np.abs(np.max(avg_reg_diff_df))])

sns.heatmap(avg_reg_diff_df, ax=axs[1], cmap=sns.diverging_palette(250, 10, as_cmap=True),\
    vmin=-max_chan, vmax=max_chan, cbar_kws={'label': opts})
axs[1].set_title("V1 and M1 " + opts)

axs[1].set_ylabel("Gene Labels")
plt.savefig(savepath + opts + '_Reg_Diff.png', format='png')
plt.savefig(savepath + opts+ '_Reg_Diff.pdf', format='pdf')
plt.show()

# %% Select significantly different genes between motor and visual cortex
m1_trim_vals = motor_dict['trimmed_exp_vals']
v1_trim_vals = visual_dict['trimmed_exp_vals']

sig_genes = []
# Loop through each gene and get the p-values between motor and visual cortex
for gene in marker_genes:
    m1_gene_vals = m1_trim_vals[m1_trim_vals['gene'] == gene]['norm'].values
    v1_gene_vals = v1_trim_vals[v1_trim_vals['gene'] == gene]['norm'].values

    # Remove cell expressions with NaNs
    m1_gene_vals = m1_gene_vals[~np.isnan(m1_gene_vals)]
    v1_gene_vals = v1_gene_vals[~np.isnan(v1_gene_vals)]

    u_stat, p_value = stats.mannwhitneyu(v1_gene_vals, m1_gene_vals, alternative='two-sided')
    print(p_value)

    if p_value < 0.05:
        sig_genes.append(gene)

#%%
# Plot individual histograms of the genes that are higher in the V1 than in M1
m1_trim_vals = motor_dict['trimmed_exp_vals']
v1_trim_vals = visual_dict['trimmed_exp_vals']

for gene in sig_genes: #v1_higher_genes:
    m1_gene_vals = m1_trim_vals[m1_trim_vals['gene'] == gene]['norm'].values
    v1_gene_vals = v1_trim_vals[v1_trim_vals['gene'] == gene]['norm'].values

    # Determine total number of cells
    m1_cell_total_nums = m1_gene_vals.shape[0]
    v1_cell_total_nums = v1_gene_vals.shape[0]

    # Remove cell expressions with NaNs
    m1_gene_vals = m1_gene_vals[~np.isnan(m1_gene_vals)]
    v1_gene_vals = v1_gene_vals[~np.isnan(v1_gene_vals)]

    u_stat, p_value = stats.mannwhitneyu(v1_gene_vals, m1_gene_vals, alternative='two-sided')

    # Calculate the Cliff's Delta
    delta = cliffs_delta(v1_gene_vals, m1_gene_vals)

    #DEBUG
    if gene == 'Atp1a1':
        print('Atp1a1: averages')
        print('V1 avg: ' + str(np.nanmean(v1_gene_vals)))
        print('M1 avg: ' + str(np.nanmean(m1_gene_vals)))

    # Calculate the histogram parts with the difference
    vis_hist, vis_bins = np.histogram(v1_gene_vals, bins=50)
    mot_hist, mot_bins = np.histogram(m1_gene_vals, bins=vis_bins) # Use the visual bins to keep it consistent
    vis_hist = vis_hist/v1_gene_vals.shape[0]
    mot_hist = mot_hist/m1_gene_vals.shape[0]
    
    # Check the histogram values
    print('Sum of histograms')
    print(np.sum(vis_hist))
    print(np.sum(mot_hist))

    gene_cell_perc_v1 = 100*(v1_gene_vals.shape[0]/v1_cell_total_nums)
    gene_cell_perc_m1 = 100*(m1_gene_vals.shape[0]/m1_cell_total_nums)
    
    plt.rcParams.update({'font.size': 4})

    plt.figure(figsize=(3*cm, 2.2*cm))
    plt.bar(mot_bins[:-1], mot_hist, width=np.diff(mot_bins), color='blue', alpha=0.5, label='M1')
    plt.bar(vis_bins[:-1], vis_hist, width=np.diff(vis_bins), color='red', alpha=0.5, label='V1')
        
    #plt.hist(mot_cpm_norm[0], bins=30, alpha=0.5, label='M1', color='blue')
    #plt.hist(vis_cpm_norm[0], bins=30, alpha=0.5, label='V1', color='red')
    
    plt.legend()
    plt.xlabel("Expression Value")
    plt.ylabel("Cell Percentage")
    plt.title(gene + ' V1: ' + str(round(gene_cell_perc_v1, 1)) +\
              ' M1: ' + str(round(gene_cell_perc_m1, 1)) + " \ncliff's deltra: " +\
                 str(delta) + " p="+ str(p_value))
        
    plt.savefig(savepath + 'Gene' + f_sep + gene + '_exp_cell_perc.png', format='png')
    plt.savefig(savepath + 'Gene' + f_sep + gene + '_exp_cell_perc.pdf', format='pdf')
    
plt.show()  


#%%
# compute trimmed means
print("Computing trimmed means for shared clusters...")

def compute_trimmed_means(args):
    cluster, sample_indices, marker_counts = args
    cluster_indices = np.isin(sample_indices, cluster["sample_index"].values)
    cluster_data = marker_counts[:, cluster_indices]
    
    cell_nums = cluster_data.shape[1]
    trim_num = int(pp_cut*cell_nums)
    
    # Grab elements that were used for the trim_mean()
    sorted_mat = np.sort(cluster_data, axis=1)

    if trim_num == 0:
        used_elements = sorted_mat
    else:
        used_elements = sorted_mat[:, trim_num:- trim_num]
    
    #DEBUG
    if cluster["cluster_label"].iloc[0] == '108_Pvalb':
        print('Num Cells ' + str(cell_nums))
        print(used_elements)


    #print(use_elements[0, :])
    #print(f"Elements kept {use_elements.size}, Total number of cells: {cell_nums}")
    result = [cluster["cluster_label"].iloc[0], trim_mean(cluster_data, proportiontocut=pp_cut, axis=1), int(cell_nums - 2*trim_num), used_elements]
    
    return result


# multiprocessing
visual_args = [
    (visual_metadata[visual_metadata["cluster_label"] == cluster], visual_sample_indices, visual_marker_counts)
    for cluster in shared_clusters
]
motor_args = [
    (motor_metadata[motor_metadata["cluster_label"] == cluster], motor_sample_indices, motor_marker_counts)
    for cluster in shared_clusters
]

visual_total_args = [
    (visual_metadata[visual_metadata["cluster_label"] == cluster], visual_sample_indices, visual_total_counts)
    for cluster in shared_clusters
]

motor_total_args = [
    (motor_metadata[motor_metadata["cluster_label"] == cluster], motor_sample_indices, motor_total_counts)
    for cluster in shared_clusters
]

with Pool(processes=10) as pool:
    motor_trimmed_means_info = pool.map(compute_trimmed_means, motor_args)
    visual_trimmed_means_info = pool.map(compute_trimmed_means, visual_args)
    motor_total_counts_info = pool.map(compute_trimmed_means, motor_total_args)
    visual_total_counts_info = pool.map(compute_trimmed_means, visual_total_args)
    
    mot_clust_cell_nums = [tup[2] for tup in motor_trimmed_means_info]
    vis_clust_cell_nums = [tup[2] for tup in visual_trimmed_means_info]
    
    vis_trimmed_means = [tup[0:2] for tup in visual_trimmed_means_info]
    mot_trimmed_means = [tup[0:2] for tup in motor_trimmed_means_info]
    vis_clust_total_counts = [tup[0:2] for tup in visual_total_counts_info]
    mot_clust_total_counts = [tup[0:2] for tup in motor_total_counts_info]
    
    # Storing the single cell expression levels for each cluster from the trimmed mean
    mot_cell_exp = [(tup[0], tup[3]) for tup in motor_trimmed_means_info]
    vis_cell_exp = [(tup[0], tup[3]) for tup in visual_trimmed_means_info]
    mot_cell_exp = dict(mot_cell_exp)
    vis_cell_exp = dict(vis_cell_exp)
    
    # Calculates the trimmed means for each cluster
    visual_trimmed_means = dict(vis_trimmed_means)
    motor_trimmed_means = dict(mot_trimmed_means)
    
    # These total counts are for each PV cluster as a trimmed_mean
    vis_total_counts = dict(vis_clust_total_counts)
    mot_total_counts = dict(mot_clust_total_counts)
    
    
motor_trimmed_means_df = pd.DataFrame.from_dict(motor_trimmed_means, orient="index", columns=marker_genes)
visual_trimmed_means_df = pd.DataFrame.from_dict(visual_trimmed_means, orient="index", columns=marker_genes)


#%% plotting with just trimmed means
output_pdf = "/projectnb/hanlab/genesPierre/trimmed_means_heatmaps.pdf"

with PdfPages(output_pdf) as pdf:
    plt.figure(figsize=(10, 6))
    sns.heatmap(visual_trimmed_means_df, cmap="viridis", annot=False, vmin=0, vmax=70, cbar_kws={'label': 'Trimmed Mean'})
    plt.title("Trimmed Means, Visual Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    #pdf.savefig()  
    plt.show()  

    plt.figure(figsize=(10, 6))
    sns.heatmap(motor_trimmed_means_df, cmap="viridis", annot=False, vmin=0, vmax=70, cbar_kws={'label': 'Trimmed Mean'})
    plt.title("Trimmed Means, Motor Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    #pdf.savefig()  
    plt.show()  

#%% CPM normalization with trimmed means heatmaps
# Note: that the vis_total_counts is value of the trimmed_mean of total gene expression of cells in that cluster

# Construct the total_gene expression normalized
motor_trimmed_means_norm = {clust:\
                            100000*(motor_trimmed_means[clust]/mot_total_counts[clust]) + 1\
                                for clust in motor_trimmed_means}
visual_trimmed_means_norm = {clust:\
                             100000*(visual_trimmed_means[clust]/vis_total_counts[clust]) + 1\
                             for clust in visual_trimmed_means}

motor_trimmed_means_norm_df = pd.DataFrame.from_dict(motor_trimmed_means_norm, orient="index", columns=marker_genes)
visual_trimmed_means_norm_df = pd.DataFrame.from_dict(visual_trimmed_means_norm, orient="index", columns=marker_genes)


plt.figure(figsize=(8, 6))
sns.heatmap(np.log2(motor_trimmed_means_norm_df), cmap="viridis", annot=False, vmin=0, vmax=10, cbar_kws={'label': 'Trimmed Mean'})
plt.title("Norm Trimmed Means, Motor Cortex")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45)
#pdf.savefig()

plt.figure(figsize=(8, 6))
sns.heatmap(np.log2(visual_trimmed_means_norm_df), cmap="viridis", annot=False, vmin=0, vmax=10, cbar_kws={'label': 'Trimmed Mean'})
plt.title("Norm Trimmed Means, Visual Cortex")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45)
#pdf.savefig()

plt.show()

# %% Calculate the weighted averages for each gene and plot the weighted average based on individual clusters
# Calculating the weighted averages based on each cluster and trimmed means

# TODO need to re-change this for the new dictionary structures
weighted_avg_m1 = pd.Series(
    np.average(motor_trimmed_means_norm_df, axis=0, weights=mot_clust_cell_nums),
    index=motor_trimmed_means_df.columns,
    name='M1'
)

weighted_avg_v1 = pd.Series(
    np.average(visual_trimmed_means_norm_df, axis=0, weights=vis_clust_cell_nums),
    index=visual_trimmed_means_df.columns,
    name='V1'
)

weighted_df = pd.concat([weighted_avg_m1, weighted_avg_v1], axis=1)
weighted_diff_df = (weighted_avg_v1 - weighted_avg_m1)/(weighted_avg_v1 + weighted_avg_m1)

# Print the genes that have higher expression in V1
print('Genes greater in V1 than M1')
print(weighted_diff_df[weighted_diff_df > 0].index.tolist())

heatmap_df = pd.concat([weighted_df, weighted_diff_df], axis=1)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
sns.heatmap(weighted_df, ax=axs[0], cmap=sns.diverging_palette(250, 10, as_cmap=True), vmin=1, vmax=240, cbar_kws={'label': 'Exp Values'})
axs[0].set_title("Weighted Expression")
axs[0].set_xlabel("Genes")
axs[0].set_ylabel("Gene Labels")
#axs[0].tick_params(axis='x', rotation=45)
#axs[0].set_xticklabels(axs[0].get_xticklabels(), ha='right')

weighted_diff_df = weighted_diff_df.to_frame()
sns.heatmap(weighted_diff_df, ax=axs[1], cmap=sns.diverging_palette(250, 10, as_cmap=True), vmin=-0.06, vmax=0.06, cbar_kws={'label': 'Exp Values'})
axs[1].set_title("Weighted Differences")
axs[1].set_xlabel("Genes")
axs[1].set_ylabel("Gene Labels")

plt.savefig(savepath + 'CPM_Norm_Weighted.png', format='png')
plt.savefig(savepath + 'CPM_Norm_Weighted.pdf', format='pdf')
plt.show()

# %% Pools all of the expression levels together and performs CPM normalization
# Calculate the histogram for a specific gene betweeen M1 and V1
# Need to use per cell because trimmed means does not work

# Conslidate all of the expression values from each cluster
# This is wrong because we need to normalize by CELL not cluster
# mot_cell_exp_df = [ mot_cell_exp[clust] / mot_total_counts[clust] for clust in mot_cell_exp]
# vis_cell_exp_df = [ vis_cell_exp[clust] / vis_total_counts[clust] for clust in vis_cell_exp]

# mot_cell_exp_df = pd.DataFrame(np.hstack(mot_cell_exp_df).T, columns=marker_genes)
# vis_cell_exp_df = pd.DataFrame(np.hstack(vis_cell_exp_df).T, columns=marker_genes)

# mot_cell_exp_df.replace(0, np.nan, inplace=True)
# vis_cell_exp_df.replace(0, np.nan, inplace=True)


# Total Cell Number for each region
vis_cell_num = visual_total_counts.shape[1]
mot_cell_num = motor_total_counts.shape[1]

# Perform CPM normalization
vis_cpm_norm = 100000*(visual_marker_counts/visual_total_counts) + 1
mot_cpm_norm = 100000*(motor_marker_counts/motor_total_counts) + 1

# Replace 0's with NaN
vis_cpm_norm[vis_cpm_norm == 1] = np.nan
mot_cpm_norm[mot_cpm_norm == 1] = np.nan

# Count the number of cells that have an expression value (greter than 0) per gene
vis_cellw_exp_num = np.count_nonzero(~np.isnan(vis_cpm_norm), axis=1)
mot_cellw_exp_num = np.count_nonzero(~np.isnan(mot_cpm_norm), axis=1)

#vis_trim_num = int(pp_cut*vis_cell_num)
#mot_trim_num = int(pp_cut*mot_cell_num)
#vis_cpm_norm = np.sort(vis_cpm_norm, axis=1)[:, vis_trim_num:]
#mot_cpm_norm = np.sort(mot_cpm_norm, axis=1)[:, mot_trim_num:]


# Store all of the normalized CPM into a pickle file
norm_cpm_dict = {}

#goi = ['Scn2a', 'Kcna1', 'Kcna2', 'Kcnc2']
#goi = ['Kcnc3', 'Kcna1', 'Hcn2', 'Grin1', 'Scn9a', 'Cacna1a', 'Cacna1b']
goi = ['Kcna1', 'Hcn2', 'Scn9a', 'Cacna1a', 'Cacna1b', 'Grin1', 'Atp1a1']
for gene in goi:
   
    g_index = marker_genes.index(gene)
    #TODO manually change bin width
    # Normalize the expression frequencies by the total number of cells in each region
    vis_cpm_data = vis_cpm_norm[g_index, ~np.isnan(vis_cpm_norm[g_index, :])]
    mot_cpm_data = mot_cpm_norm[g_index, ~np.isnan(mot_cpm_norm[g_index, :])]
    
    # Add values to the dictionary
    norm_cpm_dict[gene] = {}
    norm_cpm_dict[gene]['v1'] = vis_cpm_data
    norm_cpm_dict[gene]['m1'] = mot_cpm_data
        
    # Calculate the P-value for the CPM difference
    u_stat, p_value = stats.mannwhitneyu(vis_cpm_data, mot_cpm_data, alternative='two-sided')
    
    # Calculate the Cliff's Delta
    delta = cliffs_delta(vis_cpm_data, mot_cpm_data)

    vis_hist, vis_bins = np.histogram(vis_cpm_data, bins=50)
    mot_hist, mot_bins = np.histogram(mot_cpm_data, bins=vis_bins)
    
    vis_hist = vis_hist/vis_cellw_exp_num[g_index]
    mot_hist = mot_hist/mot_cellw_exp_num[g_index]
    
    gene_cell_perc_v1 = 100*(vis_cellw_exp_num[g_index]/vis_cell_num)
    gene_cell_perc_m1 = 100*(mot_cellw_exp_num[g_index]/mot_cell_num)
    
    plt.figure(figsize=(10, 6))
    plt.bar(mot_bins[:-1], mot_hist, width=np.diff(mot_bins), color='blue', alpha=0.5, label='M1')
    plt.bar(vis_bins[:-1], vis_hist, width=np.diff(vis_bins), color='red', alpha=0.5, label='V1')
        
    #plt.hist(mot_cpm_norm[0], bins=30, alpha=0.5, label='M1', color='blue')
    #plt.hist(vis_cpm_norm[0], bins=30, alpha=0.5, label='V1', color='red')
    
    plt.legend()
    plt.xlabel("Expression Value")
    plt.ylabel("Cell Percentage")
    plt.title(gene + ' V1: ' + str(round(gene_cell_perc_v1, 1)) +\
              ' M1: ' + str(round(gene_cell_perc_m1, 1)) + " cliff's deltra: " +\
                 str(delta) + " p="+ str(p_value))
        
    plt.savefig(savepath + 'Gene' + f_sep + gene + '_exp_cell_perc.png', format='png')
    plt.savefig(savepath + 'Gene' + f_sep + gene + '_exp_cell_perc.pdf', format='pdf')
    
plt.show()

# Save the values from the dictionary
#with open(savepath + "cell_counts_goi.pkl", 'wb') as f:
#    pickle.dump(norm_cpm_dict, f)

#%% Heatmap with total cell averages and difference between two (NOT using trimmed mean values from each cluster)
# This section requires the CPM normalization lines 364-365
avg_gene_exp_m1 = pd.Series(np.nanmean(mot_cpm_norm, axis=1) )
avg_gene_exp_m1.index = marker_genes
avg_gene_exp_v1 = pd.Series(np.nanmean(vis_cpm_norm, axis=1) )
avg_gene_exp_v1.index = marker_genes

# Expressions for both M1 and V1
avg_exp_df = pd.concat([avg_gene_exp_v1, avg_gene_exp_m1], axis=1)
avg_diff_exp_df = (avg_gene_exp_v1 - avg_gene_exp_m1)/(avg_gene_exp_v1 + avg_gene_exp_m1)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
sns.heatmap(avg_exp_df, ax=axs[0], cmap=sns.diverging_palette(250, 10, as_cmap=True), vmin=1, vmax=240, cbar_kws={'label': 'Exp Values'})
axs[0].set_title("Cell Average Expression")
axs[0].set_xlabel("Genes")
axs[0].set_ylabel("Gene Labels")
#axs[0].tick_params(axis='x', rotation=45)
#axs[0].set_xticklabels(axs[0].get_xticklabels(), ha='right')

avg_diff_exp_df = avg_diff_exp_df.to_frame()
sns.heatmap(avg_diff_exp_df, ax=axs[1], cmap=sns.diverging_palette(250, 10, as_cmap=True), vmin=-0.06, vmax=0.06, cbar_kws={'label': 'Exp Values'})
axs[1].set_title("Average Expression Difference")
axs[1].set_xlabel("Genes")
axs[1].set_ylabel("Gene Labels")

plt.savefig(savepath + 'CPM_Norm_Cell_Avg.png', format='png')
plt.savefig(savepath + 'CPM_Norm_Cell_Avg.pdf', format='pdf')
plt.show()

#%% Make violin plots from the top normalization dataframes
goi = ['Scn2a', 'Kcna1', 'Kcna2', 'Kcnc2']
for gene in goi:
   
    g_index = marker_genes.index(gene)
    #TODO manually change bin width
    # Normalize the expression frequencies by the total number of cells in each region
    vis_cpm_data = vis_cpm_norm[g_index, ~np.isnan(vis_cpm_norm[g_index, :])].reshape(-1, 1)
    mot_cpm_data = mot_cpm_norm[g_index, ~np.isnan(mot_cpm_norm[g_index, :])].reshape(-1, 1)
        
    # Pad with NaNs
    max_len = max(len(vis_cpm_data), len(vis_cpm_data))
    vis_cpm_data_pd = np.pad(vis_cpm_data, ((0, max_len - len(vis_cpm_data) ), (0, 0)), constant_values=np.nan)
    mot_cpm_data_pd = np.pad(mot_cpm_data, ((0, max_len - len(mot_cpm_data) ), (0, 0)), constant_values=np.nan)
    
    # Calculate the average expression level
    avg_v1 = np.median(vis_cpm_data)
    avg_m1 = np.median(mot_cpm_data)
    
    data = [vis_cpm_data.flatten(), mot_cpm_data.flatten()]
    
    plt.figure(figsize=(10, 6))
    sns.violinplot(data, inner='box', alpha=0.5)
    plt.xticks([0, 1], ["V1", "M1"])
    plt.hlines(y=[avg_v1, avg_m1], xmin=0, xmax=2, color='k')
    
    
    plt.legend()
    plt.xlabel("")
    plt.ylabel("Expression Vals")
    plt.title(gene + ' V1: ' + str(round(gene_cell_perc_v1, 1)) +\
              ' M1: ' + str(round(gene_cell_perc_m1, 1)))
    
plt.show()

#%% Heatmap with genes on Y-axis and brain region on X-axis
m1_df = motor_trimmed_means_norm_df.T
v1_df = visual_trimmed_means_norm_df.T

m1_df = m1_df.loc[order_genes]
v1_df = v1_df.loc[order_genes]

cur_col = m1_df.columns
new_col = pd.MultiIndex.from_tuples([('M1', col) for col in cur_col])
m1_df.columns = new_col

cur_col = v1_df.columns
new_col = pd.MultiIndex.from_tuples([('V1', col) for col in cur_col])
v1_df.columns = new_col

heatmap_df = pd.concat([m1_df, v1_df], axis=1)

plt.figure(figsize=(14, 5))
sns.heatmap(np.log2(heatmap_df), cmap="viridis", annot=False, vmin=0, vmax=10, cbar_kws={'label': ' Log2(Trimmed Mean)'})
plt.vlines(x=[m1_df.shape[1], ], ymin=0, ymax=26, color='k')
plt.title("M1 and V1 Gene Expression")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45, ha='right')
plt.savefig(savepath + 'CPM_Norm_M1_V1_Group.png', format='png')
plt.savefig(savepath + 'CPM_Norm_M1_V1_Group.pdf', format='pdf')
plt.show()

#%% Merge Motor and Visual Cortex together
m1_df = motor_trimmed_means_norm_df
v1_df = visual_trimmed_means_norm_df

reg_merg_df = pd.DataFrame()

#existing_columns = df.columns

# Create the new column level (for example, 'group1')
#new_columns = pd.MultiIndex.from_tuples([('group1', col) for col in existing_columns])

for gene in m1_df.columns:
    
    m1_series = m1_df[gene]
    m1_series.name = 'm1'
    
    v1_series = v1_df[gene]
    v1_series.name = 'v1'
    
    comb_df = pd.concat([m1_series, v1_series], axis=1)
    exist_col = comb_df.columns
    new_columns = pd.MultiIndex.from_tuples([(gene, col) for col in exist_col])
    comb_df.columns = new_columns 
    
    reg_merg_df = pd.concat([reg_merg_df, comb_df], axis=1)

plt.rcParams['pdf.fonttype'] = 42

plt.figure(figsize=(18, 6))
sns.heatmap(np.log2(reg_merg_df), cmap=sns.diverging_palette(250, 10, as_cmap=True), annot=False, vmin=0, vmax=10, cbar_kws={'label': 'Trimmed Mean'})
plt.vlines(x=np.arange(2, 50, 2), ymin=0, ymax=26, color='k')
plt.title("Motor then Visual")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + 'Gene' f 'CPM_Norm_M1_V1_per_gene.png', format='png')
plt.savefig(savepath + 'Gene' f 'CPM_Norm_M1_V1_per_gene.pdf', format='pdf')
plt.show()


# %% Normalization with pvalb, updated way

m1_df = motor_trimmed_means_norm_df
v1_df = visual_trimmed_means_norm_df

m1_df = m1_df.div(m1_df['Pvalb'], axis=0)
v1_df = v1_df.div(v1_df['Pvalb'], axis=0)

#m1_df = m1_df/m1_pval
#v1_df = v1_df/v1_pval

#m1_pval = m1_df['Pvalb']
#m1_pval.index = m1_df.index

#v1_pval = v1_df['Pvalb']
#v1_pval.index = v1_df.index

reg_merg_df = pd.DataFrame()

#existing_columns = df.columns

# Create the new column level (for example, 'group1')
#new_columns = pd.MultiIndex.from_tuples([('group1', col) for col in existing_columns])

for gene in m1_df.columns:
    m1_series = m1_df[gene]
    m1_series.name = 'm1'
    
    v1_series = v1_df[gene]
    v1_series.name = 'v1'
    
    comb_df = pd.concat([m1_series, v1_series], axis=1)
    exist_col = comb_df.columns
    new_columns = pd.MultiIndex.from_tuples([(gene, col) for col in exist_col])
    comb_df.columns = new_columns 
    
    reg_merg_df = pd.concat([reg_merg_df, comb_df], axis=1)

plt.figure(figsize=(18, 6))
sns.heatmap(np.log2(reg_merg_df), cmap=sns.diverging_palette(250, 10, as_cmap=True), annot=False, vmin=0, vmax=1, cbar_kws={'label': 'Trimmed Mean'})
plt.vlines(x=np.arange(2, 50, 2), ymin=0, ymax=26, color='k')
plt.title("Motor then Visual")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + 'PValb_Norm_M1_V1_per_gene.png', format='png')
plt.show()

# %% Normalization with SCN1A (Sodium Channel), updated way

m1_df = motor_trimmed_means_norm_df
v1_df = visual_trimmed_means_norm_df

m1_df = m1_df.div(m1_df['Scn1a'], axis=0)
v1_df = v1_df.div(v1_df['Scn1a'], axis=0)

#m1_df = m1_df/m1_pval
#v1_df = v1_df/v1_pval

#m1_pval = m1_df['Pvalb']
#m1_pval.index = m1_df.index

#v1_pval = v1_df['Pvalb']
#v1_pval.index = v1_df.index

reg_merg_df = pd.DataFrame()

#existing_columns = df.columns

# Create the new column level (for example, 'group1')
#new_columns = pd.MultiIndex.from_tuples([('group1', col) for col in existing_columns])

for gene in m1_df.columns:
    m1_series = m1_df[gene]
    m1_series.name = 'm1'
    
    v1_series = v1_df[gene]
    v1_series.name = 'v1'
    
    comb_df = pd.concat([m1_series, v1_series], axis=1)
    exist_col = comb_df.columns
    new_columns = pd.MultiIndex.from_tuples([(gene, col) for col in exist_col])
    comb_df.columns = new_columns 
    
    reg_merg_df = pd.concat([reg_merg_df, comb_df], axis=1)

plt.figure(figsize=(18, 6))
sns.heatmap(np.log2(reg_merg_df), cmap=sns.diverging_palette(250, 10, as_cmap=True), annot=False, vmin=0, vmax=1, cbar_kws={'label': 'Trimmed Mean'})
plt.vlines(x=np.arange(2, 50, 2), ymin=0, ymax=26, color='k')
plt.title("Motor then Visual")
plt.xlabel("Genes")
plt.ylabel("Cluster Labels")
plt.xticks(rotation=45, ha='right')

#plt.tick_params(axis='x', which='both', labelbottom=True, pad=10)
plt.savefig(savepath + 'Scn1a_Norm_M1_V1_per_gene.png', format='png')
plt.show()



# ------- ----- other ways and attempts ------------------------

# %% Read in metadata and keep the motor and visual cortices separate clusters (-- This does not seem to matter --)

# filter metadata for motor and visual cortex
print("Filtering metadata for motor and visual cortex...")
metadata_df = pd.read_csv(metadata_file)
metadata_filtered = metadata_df[
    metadata_df["region_label"].str.contains("VIS|MOp", na=False) &
    metadata_df["cluster_label"].str.contains("Pvalb", na=False)
][["sample_name", "region_label", "cluster_label"]]
print(f"Filtered metadata shape: {metadata_filtered.shape}")

# match metadata_filtered with data/sample
print("Mapping metadata_filtered to data/sample...")
with h5py.File(hdf5_file, "r") as f:
    samples = f["data/samples"][:].astype(str)
sample_to_index = {sample: idx for idx, sample in enumerate(samples)}

# find indices of samples in data/samples that match metadata_filtered
metadata_filtered["sample_index"] = metadata_filtered["sample_name"].map(sample_to_index)
metadata_filtered = metadata_filtered.dropna(subset=["sample_index"])
metadata_filtered["sample_index"] = metadata_filtered["sample_index"].astype(int)
print(f"Number of valid samples after matching: {len(metadata_filtered)}")


# separate metadata into visual and motor cortex
visual_metadata = metadata_filtered[metadata_filtered["region_label"].str.contains("VIS", na=False)]
motor_metadata = metadata_filtered[metadata_filtered["region_label"].str.contains("MOp", na=False)]
print(f"Visual cortex samples: {len(visual_metadata)}")
print(f"Motor cortex samples: {len(motor_metadata)}")

# use sample indices to filter marker gene data for each region
print("Filtering marker gene data for visual and motor cortex...")
visual_sample_indices = visual_metadata["sample_index"].values
motor_sample_indices = motor_metadata["sample_index"].values
visual_marker_counts = marker_gene_data[:, visual_sample_indices]   
motor_marker_counts = marker_gene_data[:, motor_sample_indices]

visual_total_counts = cell_counts_total[visual_sample_indices].reshape(1, -1)
motor_total_counts = cell_counts_total[motor_sample_indices].reshape(1, -1)

# compute trimmed means 
print("Computing trimmed means for cortical-specific clusters...")
def compute_trimmed_means(args):
    cluster, sample_indices, marker_counts = args
    cluster_indices = np.isin(sample_indices, cluster["sample_index"].values)
    cluster_data = marker_counts[:, cluster_indices]
    return cluster["cluster_label"].iloc[0], trim_mean(cluster_data, proportiontocut=0.1, axis=1)


# multiprocessing
visual_args = [
    (visual_metadata[visual_metadata["cluster_label"] == cluster], visual_sample_indices, visual_marker_counts)
    for cluster in visual_metadata["cluster_label"].unique()
]
motor_args = [
    (motor_metadata[motor_metadata["cluster_label"] == cluster], motor_sample_indices, motor_marker_counts)
    for cluster in motor_metadata["cluster_label"].unique()
]

visual_total_args = [
    (visual_metadata[visual_metadata["cluster_label"] == cluster], visual_sample_indices, visual_total_counts)
    for cluster in visual_metadata["cluster_label"].unique()
]

motor_total_args = [
    (motor_metadata[motor_metadata["cluster_label"] == cluster], motor_sample_indices, motor_total_counts)
    for cluster in motor_metadata["cluster_label"].unique()
]

with Pool(processes=10) as pool:
    visual_trimmed_means = dict(pool.map(compute_trimmed_means, visual_args))
    motor_trimmed_means = dict(pool.map(compute_trimmed_means, motor_args))
    vis_total_trimmed_means = dict(pool.map(compute_trimmed_means, visual_total_args))
    mot_total_trimmed_means = dict(pool.map(compute_trimmed_means, motor_total_args))


motor_trimmed_means_df = pd.DataFrame.from_dict(motor_trimmed_means, orient="index", columns=marker_genes)
visual_trimmed_means_df = pd.DataFrame.from_dict(visual_trimmed_means, orient="index", columns=marker_genes)

#%% normalization with pvalb


print("Normalizing trimmed means across genes (Z-score normalization)...")

visual_trimmed_means_df_normalized = visual_trimmed_means_df.apply(zscore, axis=1, nan_policy='omit')
motor_trimmed_means_df_normalized = motor_trimmed_means_df.apply(zscore, axis=1, nan_policy='omit')

print("Normalizing by Pvalb...")
pvalb_index = visual_trimmed_means_df.columns.get_loc("Pvalb")
visual_trimmed_means_df_normalized_by_pvalb = visual_trimmed_means_df_normalized.div(visual_trimmed_means_df_normalized["Pvalb"], axis=0)
motor_trimmed_means_df_normalized_by_pvalb = motor_trimmed_means_df_normalized.div(motor_trimmed_means_df_normalized["Pvalb"], axis=0)

#output_pdf = "/projectnb/hanlab/genesPierre/normalized_trimmed_means_by_pvalb_heatmaps.pdf"
output_pdf = "/ad/eng/research/eng_research_handata/Pierre Fabris/GenePVPlots"

with PdfPages(output_pdf + 'norm_pval_trim_heatmap') as pdf:
    plt.figure(figsize=(10, 6))
    sns.heatmap(visual_trimmed_means_df_normalized_by_pvalb, cmap="viridis", annot=False, vmin=-15, vmax=15, cbar_kws={'label': 'Normalized Value'})
    plt.title("Normalized Trimmed Means (by Pvalb), Visual Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    #pdf.savefig()  
    plt.show()  

    plt.figure(figsize=(10, 6))
    sns.heatmap(motor_trimmed_means_df_normalized_by_pvalb, cmap="viridis", annot=False, vmin=-15, vmax=15, cbar_kws={'label': 'Normalized Value'})
    plt.title("Normalized Trimmed Means (by Pvalb), Motor Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    #pdf.savefig()  
    plt.show()  

print(f"Heatmaps saved to {output_pdf}!")

#%% normalization with Scn1a


print("Normalizing trimmed means across genes (Z-score normalization)...")

visual_trimmed_means_df_normalized = visual_trimmed_means_df.apply(zscore, axis=1, nan_policy='omit')
motor_trimmed_means_df_normalized = motor_trimmed_means_df.apply(zscore, axis=1, nan_policy='omit')

print("Normalizing by Scn1a...")
#pvalb_index = visual_trimmed_means_df.columns.get_loc("Pvalb")
visual_trimmed_means_df_normalized_by_scn1a = visual_trimmed_means_df_normalized.div(visual_trimmed_means_df_normalized["Scn1a"], axis=0)
motor_trimmed_means_df_normalized_by_scn1a = motor_trimmed_means_df_normalized.div(motor_trimmed_means_df_normalized["Scn1a"], axis=0)

#output_pdf = "/projectnb/hanlab/genesPierre/normalized_trimmed_means_by_pvalb_heatmaps.pdf"
output_pdf = "/ad/eng/research/eng_research_handata/Pierre Fabris/GenePVPlots"

with PdfPages(output_pdf + 'norm_pval_trim_heatmap') as pdf:
    plt.figure(figsize=(10, 6))
    sns.heatmap(visual_trimmed_means_df_normalized_by_pvalb, cmap="viridis", annot=False, vmin=-15, vmax=15, cbar_kws={'label': 'Normalized Value'})
    plt.title("Normalized Trimmed Means (by Scn1a), Visual Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    #pdf.savefig()  
    plt.show()  

    plt.figure(figsize=(10, 6))
    sns.heatmap(motor_trimmed_means_df_normalized_by_pvalb, cmap="viridis", annot=False, vmin=-15, vmax=15, cbar_kws={'label': 'Normalized Value'})
    plt.title("Normalized Trimmed Means (by Scn1a), Motor Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    #pdf.savefig()  
    plt.show()  

print(f"Heatmaps saved to {output_pdf}!")


#%% normalization without pvalb

print("Normalizing trimmed means across genes (Z-score normalization)...")
visual_trimmed_means_df_normalized = visual_trimmed_means_df.apply(zscore, axis=1, nan_policy='omit')
motor_trimmed_means_df_normalized = motor_trimmed_means_df.apply(zscore, axis=1, nan_policy='omit')

output_pdf = "/projectnb/hanlab/genesPierre/normalized_trimmed_means_heatmaps.pdf"

with PdfPages(output_pdf) as pdf:
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(visual_trimmed_means_df_normalized, cmap="viridis", annot=False, vmin=0, vmax=4, 
                      cbar_kws={'label': 'Normalized Value'}, rasterized=False)  # Ensure vectorization
    plt.title("Normalized Trimmed Means, Visual Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    pdf.savefig()
    plt.close() 

    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(motor_trimmed_means_df_normalized, cmap="viridis", annot=False, vmin=0, vmax=4, 
                      cbar_kws={'label': 'Normalized Value'}, rasterized=False)  # Ensure vectorization
    plt.title("Normalized Trimmed Means, Motor Cortex")
    plt.xlabel("Genes")
    plt.ylabel("Cluster Labels")
    plt.xticks(rotation=45)
    pdf.savefig() 
    plt.close() 
print(f"Heatmaps saved to {output_pdf}!")
