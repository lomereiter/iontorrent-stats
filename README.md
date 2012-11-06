Exploring dataset [B7-295.bam](https://s3.amazonaws.com/iontorrent/datasets/B7-295/B7-295.bam.gz)
produced by IonTorrent 318 chip on August 26, 2012.

---------------------------------------------

### Description of project files

* U00096.2.fas
    Reference sequence for the dataset.

The following files require [sambamba](https://github.com/lomereiter/sambamba) library and latest stable [DMD](http://dlang.org/download.html) compiler to be installed:

* flowindices.d
    Gets distribution of flow indices contributing to actual base calls.

* counts.d
    Fetches information about numbers of mismatches, insertions, and deletions, per offset from the beginning of a read.

* deletions.d
    Collects more specific information about deletions, such as deleted sequence, 
    surrounding nucleotides and corresponding flow signal intensities, offset from the beginning, 
    and a few other, less important, variables. Produces ~200MB text file.
