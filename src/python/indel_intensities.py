from common import *

def plot(in_fn, out_fn, xlabel="Intensity"):
    data = read_table(in_fn)
    plot_intensity_distribution(data, out_fn, xlabel)
