# %%
import matlab.engine
import numpy as np

eng = matlab.engine.start_matlab()

# %%
data1 = np.random.rand(10) * 1
data2 = np.random.rand(10) * 3
data = np.concatenate((data1, data2))

labels = ['d1' for _ in range(10)]
labels.extend(['d2' for _ in range(10)])

print(labels)
# %%

matlab_exp = "opts.ShowMedian = true;\
            opts.ShowMean = false;\
            opts.MedianColor = [1 1 1];\
            opts.MarkerSize = 5;\
            opts.MedianMarkerSize = 5;\
            opts.BoxWidth = 0.1;\
            opts.BoxColor = [0 0 0];\
            opts.ViolinAlpha = {[0.3], [0.3]};"

# Specifying options for a clean violin plot
opts = {'ShowMedian':True, 'ShowMean':False, 'MedianColor':np.array([1, 0, 0]), 'MarkerSize':5,\
        'MedianMarkerSize':10, 'BoxWidth':0.1, 'BoxColor':np.array([0, 0, 0]), 'ViolinAlpha':[np.array([0.3]), np.array([0.3])]}

# Passing variables into workspace
eng.eval(matlab_exp, nargout=0)
eng.workspace['data'] = data
eng.workspace['labels'] = labels

eng.figure()
#eng.eval("data", nargout=0)
eng.eval("v = violinplot(data, labels, opts);", nargout=0)
eng.eval("v(1).ViolinColor = {[0, 0, 1]};", nargout=0)
eng.eval("v(2).ViolinColor = {[1, 0, 1]};", nargout=0)

# %%
