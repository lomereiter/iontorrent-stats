#!/usr/bin/env sh

VERSION="0.1"

# for testing on local machine
#
#export PLUGINNAME='iontorrent-stats'
#export TSP_FILEPATH_BAM='~/bam/B7-295.100k.bam'
#export RESULTS_DIR='/tmp'
#export TSP_FLOWORDER='TACGTACGTCTGAGCATCGATCGATGTACAGC'
#export TSP_NUM_FLOWS=1024
#export TSP_LIBRARY_KEY='TCAG'


OUTFILE=${RESULTS_DIR}/${PLUGINNAME}_block.html

${DIRNAME}/build/iontorrent-stats -e -C -d $RESULTS_DIR/reports $TSP_FILEPATH_BAM 
${DIRNAME}/src/python/contentgen.py $RESULTS_DIR/reports $RESULTS_DIR/images
${DIRNAME}/src/python/htmlgen.py $OUTFILE

# copy JS/CSS
cp -R ${DIRNAME}/js/ $RESULTS_DIR
cp -R ${DIRNAME}/css/ $RESULTS_DIR
