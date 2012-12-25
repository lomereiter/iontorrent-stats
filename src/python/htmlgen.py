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
    def __init__(self, content_dir, outfile):
        self.content_dir = content_dir
        self.outfile = outfile
        self.dest_dir = os.path.dirname(outfile)

        if not os.path.exists(self.dest_dir):
            os.makedirs(self.dest_dir)

    def run(self):
        with open(self.outfile, 'w+') as f:
            f.write(render_to_string('iontorrent-stats_block.djt', { 
                                            'overcall_plot': 
                                                os.path.join(self.content_dir, 'overcalls.png')
                                          , 'undercall_plot':
                                                os.path.join(self.content_dir, 'undercalls.png')
                                          , 'intensity_plot':
                                                os.path.join(self.content_dir, 'intensities.png')
                                          , 'error_frequency_plot':
                                                os.path.join(self.content_dir, 'error_frequencies.png')
                                                  }))

if __name__ == '__main__':
    gen = HtmlGenerator(sys.argv[1], sys.argv[2])
    gen.run()
