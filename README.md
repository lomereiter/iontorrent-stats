IonTorrent-Stats
================

Tool for collecting various statistics about BAM files produced by
recent versions of IonTorrent Suite.

## Installation

[DMD](http://dlang.org/download) compiler must be installed in order to build the application.

```sh
git clone --recursive https://github.com/lomereiter/iontorrent-stats
cd iontorrent-stats
make
```

## Usage

```sh
./build/iontorrent-stats <input.bam> -d <output directory>
```
(Mapped reads in the BAM file must contain ZF, FZ, and MD tags.)

Output directory is automatically created if it doesn't exist.

In case of successful execution, several tab-delimited files will be
written to there. All of them start with comments describing the variables,
followed by a header line. They can be easily loaded into R via standard
`read.table` function.

Directory `src/python` contains Python scripts for visualizing these
files by means of Matplotlib library. 

## License

IonTorrent-Stats is distributed under GPLv2+ license.

## TODO

* Collect mismatch statistics
* Add more plot generation scripts
* Add Python script for generating HTML pages (via Django templates)
* Add `launch.sh` script for integration with Torrent Browser Plugin Store
