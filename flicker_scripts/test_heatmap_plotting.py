import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.transforms import Bbox
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap
import sys
from glob import glob
from pymatreader import read_mat
from scipy import signal
from scipy import stats
from scipy.ndimage import uniform_filter1d

import itertools

# Create a simple data frame
# Specify the number of repetitions
num_repetitions = 5

# Create a list of repeated numbers from 1 to 5
data = [[i] * num_repetitions for i in range(1, 6)]

# Create a DataFrame from the list
df = pd.DataFrame(data)
timeline = df.columns.values
row_idx = np.arange(df.shape[0]) + 1

plt.figure()
pcol_plot = plt.pcolormesh(timeline, row_idx, df.values)
plt.colorbar(pcol_plot)
plt.show()

print(df)