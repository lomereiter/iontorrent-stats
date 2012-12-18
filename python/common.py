import sys

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.mlab import csv2rec

def make_console_app(fn):
    def print_usage():
        print("Usage: %s <input.dat> <output.png>" % sys.argv[0])

    if len(sys.argv) < 3:
        print_usage()
        sys.exit(0)

    input = sys.argv[1]
    output = sys.argv[2]
    fn(input, output)
