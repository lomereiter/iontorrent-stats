from common import *

def plot(in_fn, out_fn):
    data = read_table(in_fn)
    data = rec_groupby(data, ('intensity',), 
                             [('count', sum, 'count'), 
                              ('ins', sum, 'ins'), 
                              ('del', sum, 'del')])

    data = data[data['count'] > 10]

    i_freq = data['ins'] * 1.0 / data['count']
    d_freq = data['del'] * 1.0 / data['count']

    intensities = data['intensity'] / 100.0

    plt.plot(intensities, i_freq, 'b.-', label="Overcall frequency")
    plt.plot(intensities, d_freq, 'g.-', label="Undercall frequency")

    plt.xlabel("Signal intensity")
    plt.ylabel("Error frequencies")

    plt.xticks(np.arange(0.5, np.max(intensities)))

    plt.legend(loc='upper left')

    plt.savefig(out_fn)
