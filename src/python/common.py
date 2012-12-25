import sys

import numpy as np

import matplotlib

# required in order to work without X-server running
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.mlab import csv2rec, rec_groupby

def read_table(filename, converter_dict=None, **kwargs):
    return csv2rec(filename, delimiter='\t', converterd=converter_dict, **kwargs)

def plot_intensity_distribution(data, out_fn, xlabel="Intensity"):
    approx_density = (data['count'] / 0.01) / np.sum(data['count'])

    intensities = data['intensity'] / 100.0

    plt.plot(intensities, approx_density, 'b-')
    plt.fill_between(intensities, approx_density, color="grey", alpha=0.8)
    plt.xticks(np.arange(0.5, np.max(intensities)))
    plt.xlabel(xlabel)
    plt.ylabel("Probability density")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.savefig(out_fn)
