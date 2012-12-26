#!/usr/bin/env python
import sys
import os

from django.conf import settings
from django.template.loader import render_to_string

cur_dir = os.path.dirname(__file__)

settings.configure(TEMPLATE_DIRS=(
                       os.path.join(cur_dir, 'templates'), 
                   ),
                   TEMPLATE_DEBUG=True)

class HtmlGenerator:
    def __init__(self, outfile):
        self.outfile = outfile
        self.dest_dir = os.path.dirname(outfile)

        if self.dest_dir != '' and not os.path.exists(self.dest_dir):
            os.makedirs(self.dest_dir)

    def run(self):
        with open(self.outfile, 'w+') as f:
            f.write(render_to_string('iontorrent-stats_block.djt', { 
                                            'overcall_plot': 
                                                os.path.join('images', 'overcalls.png')
                                          , 'undercall_plot':
                                                os.path.join('images', 'undercalls.png')
                                          , 'intensity_plot':
                                                os.path.join('images', 'intensities.png')
                                          , 'error_frequency_vs_offset_plot':
                                                os.path.join('images', 'err_freq_vs_oft.png')
                                          , 'error_frequency_vs_intensity_plot':
                                                os.path.join('images', 'err_freq_vs_int.png')
                                                  }))

if __name__ == '__main__':
    gen = HtmlGenerator(sys.argv[1])
    gen.run()
