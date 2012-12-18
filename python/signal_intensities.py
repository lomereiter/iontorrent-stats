#!/usr/bin/env python
from common import *

@make_console_app
def plot(in_fn, out_fn):
    data = rec_groupby(csv2rec(in_fn, delimiter='\t'),
                       ('intensity', ),
                       [('count', sum, 'count')])

    intensities = data['intensity'] / 100.0

    plt.xticks(np.arange(0.5, np.max(intensities)))
    plt.hist(intensities,
             bins=200, 
             weights=data['count'],
             normed=True)

    plt.xlabel('Signal intensity')
    plt.ylabel('Probability density')

    plt.savefig(out_fn)
