from common import *

def plot(in_fn, out_fn):
    data = rec_groupby(read_table(in_fn),
                       ('intensity', ),
                       [('count', sum, 'count')])

    plot_intensity_distribution(data, out_fn, "Signal intensity")
