#%%
import pickle
import os
from cliffs_delta import cliffs_delta

#%%
f_sep = os.sep

savepath = "/home/pierfier/handata_server/Pierre Fabris/GenePVPlots"

# Read in pickle file with cells data
with open(savepath +f_sep+ 'cell_counts_goi.pkl', 'rb') as file:
    data = pickle.load(file)

# %%
# Loop through and print out cliff's delta values
with open(savepath+f_sep+ 'cliffs_delta.txt', 'w') as file:


    for gene in data:
        delta = cliffs_delta(data[gene]['v1'], data[gene]['m1'])

        file.write(str("Gene: " + gene+ " delta:" + str(delta) +"\n"))
