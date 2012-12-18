#!/usr/bin/env python
from common import *

@make_console_app
def plot(in_fn, out_fn):
    data = csv2rec(in_fn, delimiter='\t')
    data['intensity']

    counts = {}
    for record in data:
        if record[2] not in counts:
            counts[record[2]] = 0
        counts[record[2]] += record[3]

    intensities = np.sort(counts.keys())

    plt.xticks(np.arange(0.5, np.max(intensities / 100.0)))
    plt.hist(intensities / 100.0, 
             bins=200, 
             weights=[counts[int] for int in intensities],
             normed=True)

    plt.xlabel('Signal intensity')
    plt.ylabel('Probability density')

    plt.savefig(out_fn)
