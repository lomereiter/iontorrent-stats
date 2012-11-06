#!/usr/bin/env ruby
`Rscript -e 'source("knitpost.R"); KnitPost("#{ARGV[0]}", "/iontorrent-stats/");'`
