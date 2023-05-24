import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# so we have an array that gets turned into a histogram
# arr = [3,3,3,1,5]
# wider_arr = [-1,0,3,3,7]
# hist, bin_edges = np.histogram(np.array(wider_arr), bins=5)

# for now, just make the binning an interval of the other binning for a sanity check
# hist_fine, bin_edges_fine = np.histogram(np.array(arr), bins=10)

def order_histograms(hist1, bins1, hist2, bins2):
    width1 = (bins1[-1]-bins1[0])/len(bins1)
    width2 = (bins2[-1] - bins2[0]) / len(bins2)
    return hist1, bins1, hist2, bins2 if width1 > width2 else hist2, bins2, hist1, bins1

def histogram_mult(hist1, bins1, hist2, bins2):
    '''
    Histogram multiplication interpolation. Must first determine which has the 'finer' bins and plug in accordingly.
    Note that bounds might not be preserved
    Parameters
    ----------
    hist Bigger-bin histogram values
    bin_edges Bin edge locations of bigger histogram
    hist_fine Smaller-bin histogram values
    bin_edges_fine Smaller-bin edges

    Returns
    -------
    mult, mult_edges
    '''
    # plt.figure()
    # plt.bar(bin_edges[:-1], hist, width=2, label='Regular', align='center')
    # plt.bar(bin_edges_fine[:-1], hist_fine, label='Fine', align='center')
    # plt.legend()
    # plt.show()

    # ok so start with bigger bin and go piece by piece
    # now we just need the range to be the same, which will just be a pre-processing thing
    # also need the case when the large bin is smaller than the fine bin
    # we want to map to the fine one, so cut off values

    hist, bin_edges, hist_fine, bin_edges_fine = order_histograms(hist1, bins1, hist2, bins2)

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
    return mult, bin_edges_fine
    # plt.bar(mult_bin_edges[:-1], mult, label="Multiply")
    # plt.show()
