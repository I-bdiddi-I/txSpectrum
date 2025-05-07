# data-mc comparasions
# comparisons: log energy, x and y hybrid cores, core proximity (how close hybrid core is to sd core)
# histogram both
# normalize mc to area of data hist (separate from do_datamc)
# plot (root plot idealy) (root should have x, y errorbars, but if can't do that, use stderror on mean)
# break into energy bins (low, mid, high) and compare
# (if have time) do KS test

import logging
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp
import flux

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

#def norm_datamc(data, mc):
    
def create_histogram(data, bin_width, weights=None):
    """Creates histogram bins based on a fixed bin width."""
    min_val, max_val = min(data), max(data)
    bins = np.arange(min_val, max_val + bin_width, bin_width)  # Create bins with specified width
    hist, bin_edges = np.histogram(data, bins=bins, weights=weights)
    return hist, bin_edges


def bin_data(data, bin_edges):
    bin_indicies = np.digitize(data, bin_edges) - 1
    return [data[bin_indicies == i] for i in range(len(bin_edges)-1)]
    
    
def comp_plot(data, mc, name, size, weights):
    """Normalizes MC and creates a comparison histogram between Data and MC"""
    bin_width = size  # Adjust this as needed
    data_hist, d_bin_edges = create_histogram(data, bin_width)
    mc_hist_raw, bin_edges = create_histogram(mc, bin_width, weights=weights)  # Use different bin edges for MC (have different max energies)
    binned_data = bin_data(data, d_bin_edges)
    #binned_mc = bin_data(mc, bin_edges)
    #weights = flux.reweight_true_energy(mc.Energy_true)
    #print(binned_data)

    # Normalize MC to data area
    mc_hist_weights = mc_hist_raw * (sum(data_hist) / sum(mc_hist_raw))

    # Error calculations
    data_std_y = []
    
    for i in binned_data:
        data_std_y.append(np.std(i))
    
    data_err_y = np.sqrt(data_hist) # Poisson errors for Y
    data_err_x = bin_width / 2  # X errors are half bin width

    # Bin centers for plotting
    d_bin_centers = (d_bin_edges[:-1] + d_bin_edges[1:]) / 2
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Plot
    plt.figure(figsize=(8, 6))
    plt.hist(bin_centers, bins=bin_edges, weights=mc_hist_weights, alpha=0.5, label="Normalized and Weighted MC Histogram")
    plt.errorbar(d_bin_centers, data_hist, xerr=data_err_x, yerr=data_err_y, fmt='o', label="Data with X & Y Error Bars")
    plt.xlabel(f"{name}")
    plt.ylabel("Counts")
    plt.legend()
    plt.title(f"Data vs MC Comparison (Bin Width: {bin_width})")
    plt.show()
    
    #ks_stat, ks_pval = ks_2samp(data, mc)
    #print(f"KS Test {name} Statistic: {ks_stat}, P-Value: {ks_pval}")
 
    
def do_datamc(data, mc):
    #logger.info(f"min data LE = {min(data.LogEnergy)}, max data LE = {max(data.LogEnergy)}")
    #logger.info(f"min mc LE = {min(mc.LogEnergy)}, max mc LE = {max(data.LogEnergy)}")
    weights = mc.Weights
    comp_plot(data.LogEnergy, mc.LogEnergy, "Log Energy", 0.1, weights)
    comp_plot(data.HybridCoreX, mc.HybridCoreX, "Hybrid Core (X)", 1000, weights)
    comp_plot(data.HybridCoreY, mc.HybridCoreY, "Hybrid Core (Y)", 1000, weights)
    comp_plot(data.CoreProximity, mc.CoreProximity, "Core Proximity", 0.1, weights)
    

