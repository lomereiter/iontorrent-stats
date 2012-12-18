# Plots frequency rates of mismatches, insertions, and deletions
# as functions of offset from the read start
import sys

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.mlab import csv2rec

def print_usage():
    print("Usage: ", sys.argv[0], " <input.dat> <output.png>")

def plot(in_fn, out_fn):
    data = csv2rec(in_fn, 
                   delimiter='\t',
                   converterd=dict(zip(range(6), [np.float] * 6)))

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

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print_usage()
        sys.exit(0)

    input = sys.argv[1]
    output = sys.argv[2]
    plot(input, output)
