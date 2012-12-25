IonTorrent-Stats
================

Tool for collecting various statistics about BAM files produced by
recent versions of IonTorrent Suite.

## Installation

[LDC](https://github.com/ldc-developers/ldc/downloads) compiler must be installed in order to build the application.

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

## Examples of generated plots

![Error frequency along the read](http://dl.dropbox.com/u/7916095/iontorrent-stats/plots/error_frequencies.png)

![Call signal intensity distribution](http://dl.dropbox.com/u/7916095/iontorrent-stats/plots/intensities.png)

![Undercall signal intensity distribution](http://dl.dropbox.com/u/7916095/iontorrent-stats/plots/undercalls.png)

![Overcall signal intensity distribution](http://dl.dropbox.com/u/7916095/iontorrent-stats/plots/overcalls.png)

## License

IonTorrent-Stats is distributed under GPLv2+ license.

## TODO

* Collect mismatch statistics
* Add more plot generation scripts
* _Test everything on Torrent Server Virtual Machine_
