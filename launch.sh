#!/usr/bin/env sh

# for testing on local machine
#
#export TSP_FILEPATH_BAM='~/bam/B7-295.100k.bam'
#export RESULTS_DIR='/tmp'
#export TSP_FLOWORDER='TACGTACGTCTGAGCATCGATCGATGTACAGC'
#export TSP_NUM_FLOWS=1024
#export TSP_LIBRARY_KEY='TCAG'
#export TSP_RUN_NAME='9IKNG'

./build/iontorrent-stats -d $RESULTS_DIR/$TSP_RUN_NAME/reports $TSP_FILEPATH_BAM 
./src/python/contentgen.py $RESULTS_DIR/$TSP_RUN_NAME/reports $RESULTS_DIR/$TSP_RUN_NAME/images
