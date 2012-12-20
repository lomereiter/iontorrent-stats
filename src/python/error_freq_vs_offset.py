from common import *

def plot(in_fn, out_fn):
    data = read_table(in_fn, dict(zip(range(6), [np.float] * 6)))

    m_freq = data['mismatches'] / data['total']
    i_freq = data['insertions'] / data['total']
    d_freq = data['deletions'] / data['total']

    plt.plot(data['offset'], m_freq, 'b-', label="mismatches")
    plt.plot(data['offset'], i_freq, 'r-', label="insertions")
    plt.plot(data['offset'], d_freq, 'g-', label="deletions")

    plt.xlabel("Offset, bp")
    plt.ylabel("Error frequencies")

    plt.legend()

    plt.savefig(out_fn)
