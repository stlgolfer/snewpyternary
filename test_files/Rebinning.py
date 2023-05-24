import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# so we have an array that gets turned into a histogram
arr = [3,3,3,1,5]
wider_arr = [-1,0,3,3,7]
hist, bin_edges = np.histogram(np.array(wider_arr), bins=5)

# for now, just make the binning an interval of the other binning for a sanity check
hist_fine, bin_edges_fine = np.histogram(np.array(arr), bins=10)

plt.figure()
plt.bar(bin_edges[:-1], hist, width=2, label='Regular', align='center')
plt.bar(bin_edges_fine[:-1], hist_fine, label='Fine', align='center')
plt.legend()
plt.show()
# print(hist.copy().resize(hist_fine.copy().shape))
# plt.plot(bin_edges[:-1], np.resize(hist_fine, hist.shape), label='Hist Fine Reshaped')


# ok so start with bigger bin and go piece by piece
# TODO: now we just need the range to be the same, which will just be a pre-processing thing
# TODO: also need the case when the large bin is smaller than the fine bin
# we want to map to the fine one, so cut off values
mult_bin_edges = bin_edges_fine # also for now assume they have the same range, but we would take the smaller one
def find_num_in_range(edges, graph, num) -> float:
    for i, e in enumerate(edges):
        if num <= e:
            return graph[i-1]
    return 0

plt.figure()
mult = np.zeros_like(hist_fine)
for bin_index, small_bin in enumerate(bin_edges_fine):
    if bin_index < len(bin_edges_fine) - 1 and small_bin >= min(bin_edges) and small_bin <= max(bin_edges):
        small_val = hist_fine[bin_index]
        big_val = find_num_in_range(bin_edges, hist, small_bin)
        mult[bin_index] = big_val * small_val
        print("dank")
plt.bar(mult_bin_edges[:-1], mult, label="Multiply")
plt.show()
print("yuh")
