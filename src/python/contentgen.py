#!/usr/bin/env python
import sys
import os

from common import plt

import error_freq_vs_offset
import error_freq_vs_intensity
import signal_intensities
import indel_intensities

class ContentGenerator:
    def __init__(self, report_dir, dest_dir):
        self.report_dir = report_dir
        self.dest_dir = dest_dir

        if not os.path.exists(self.dest_dir):
            os.makedirs(self.dest_dir)

    def _make_plot(self, module, input_fn, output_fn, *args, **kwargs):
        plt.figure()
        module.plot(os.path.join(self.report_dir, input_fn),
                    os.path.join(self.dest_dir, output_fn),
                    *args, **kwargs)

    def run(self):
        self._make_plot(error_freq_vs_offset, "offsets.dat", "err_freq_vs_oft.png")

        self._make_plot(error_freq_vs_intensity, "flows.dat", "err_freq_vs_int.png")

        self._make_plot(signal_intensities, "flows.dat", "intensities.png")

        self._make_plot(indel_intensities, "overcall.intensities.dat", "overcalls.png", 
                        xlabel="Overcall intensity")

        self._make_plot(indel_intensities, "undercall.intensities.dat", "undercalls.png",
                        xlabel="Undercall intensity")

if __name__ == '__main__':
    gen = ContentGenerator(sys.argv[1], sys.argv[2])
    gen.run()
