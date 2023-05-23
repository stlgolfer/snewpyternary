import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# so we have an array that gets turned into a histogram
arr = [3,3,3,1,5]

hist, bin_edges = np.histogram(np.array(arr), bins=5)

# for now, just make the binning an interval of the other binning for a sanity check
hist_fine, bin_edges_fine = np.histogram(np.array(arr), bins=10)

# ok we ball, we just gonna get a nice kernel function and multiply them shits together
# nah, this shit way too hard

plt.figure()
plt.plot(bin_edges[:-1], hist, label='Regular')
plt.plot(bin_edges_fine[:-1], hist_fine, label='Fine')
# print(hist.copy().resize(hist_fine.copy().shape))
# plt.plot(bin_edges[:-1], np.resize(hist_fine, hist.shape), label='Hist Fine Reshaped')


# ok so start with bigger bin and go piece by piece
mult_bin_edges = bin_edges_fine # also for now assume they have the same range, but we would take the smaller one
def find_num_in_range(edges, graph, num) -> float:
    index = -1
    for i, e in enumerate(edges):
        if num >= e:
            index = i
    return graph[index]
mult = np.zeros_like(hist_fine)
for bin_index, small_bins in enumerate(bin_edges_fine):
    if bin_index < len(bin_edges_fine) - 1:
        val = hist_fine[bin_index]
        mult[bin_index] = find_num_in_range(bin_edges, hist, val) * val
plt.plot(mult_bin_edges[:-1], mult, label="Multiply")
plt.legend()
plt.show()

#TODO: now we just need the range to be the same, which will just be a pre-processing thing
print("yuh")