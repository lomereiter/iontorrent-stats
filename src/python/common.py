import sys

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.mlab import csv2rec, rec_groupby

def read_table(filename, converter_dict=None, **kwargs):
    return csv2rec(filename, delimiter='\t', converterd=converter_dict, **kwargs)
