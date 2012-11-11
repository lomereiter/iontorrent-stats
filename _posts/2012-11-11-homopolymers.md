---
layout: post
title: Homopolymer distribution
---

## Reference sequence (E. coli MG1655)

For the next section to be meaningful, let's first assess reference sequence
homopolymer distribution.


{% highlight r %}
library(ggplot2)
{% endhighlight %}



{% highlight text %}
## Loading required package: methods
{% endhighlight %}



{% highlight r %}
library(data.table)

refhoms <- data.table(read.table("_data/referencehomopolymers.dat", header=T))

ggplot(data=refhoms, 
       aes(x=factor(len), 
           y=log10(count), 
           fill=factor(nuc, c('a', 't', 'g', 'c')))) +
geom_bar() +
theme_bw() +
xlab("Homopolymer length") +
ylab(expression("lg" ~~ bgroup("(", "number of homopolymers", ")"))) +
labs(fill="Nucleotide")
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-11-homopolymers/referencehoms.png) 


Average homopolymer lengths for different nucleotides explain why indels related
to A/T bases occur about twice more often than those related to G/C:


{% highlight r %}
print(refhoms[, list(mean=weighted.mean(len, count)), by=nuc])
{% endhighlight %}



{% highlight text %}
##    nuc  mean
## 1:   a 1.420
## 2:   c 1.299
## 3:   g 1.298
## 4:   t 1.424
{% endhighlight %}


## BAM data

The file `homopolymers.dat` contains information about flow base calls,
i.e. how many times each homopolymer was called from raw flow intensity data.
That is, for each read, called bases are matched with corresponding 
flow intensities, and for each possible combination of nucleotide, called length, and 
intensity value (which is stored as an integer equal to 
`round(100.0 * normalized signal intensity)` in BAM file)
the number of such flow base calls is counted.


{% highlight r %}
library(lattice)
library(data.table)

homs <- data.table(read.table("_data/homopolymers.dat", header=T))
homs[, base := factor(base, c('T', 'C', 'A', 'G'))]
{% endhighlight %}



{% highlight text %}
##        base length intensity  count
##     1:    A      1        51 411172
##     2:    A      1        52 240347
##     3:    A      1        53 129832
##     4:    A      1        54 134885
##     5:    A      1        55 146169
##    ---                             
## 11332:    T     11      1117      1
## 11333:    T     11      1119      1
## 11334:    T     11      1124      1
## 11335:    T     11      1132      1
## 11336:    T     11      1149      2
{% endhighlight %}



{% highlight r %}

intensity.counts <- homs[, list(n=sum(count)), by=list(intensity, base)]

print(head(intensity.counts))
{% endhighlight %}



{% highlight text %}
##    intensity base      n
## 1:        51    A 411172
## 2:        52    A 240347
## 3:        53    A 129832
## 4:        54    A 134885
## 5:        55    A 146169
## 6:        56    A 162345
{% endhighlight %}



{% highlight r %}

plotIntensityDistribution <- function(df) {
    xyplot(n ~ I(intensity / 100) | base, 
           df,
           xlab="Signal intensity",
           ylab="Number of flow base calls",
           type='h',
           layout=c(2,2))
}
{% endhighlight %}


### Overall picture

{% highlight r %}
plotIntensityDistribution(intensity.counts)
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-11-homopolymers/overall.png) 


### Short homopolymers (1-2 bases)

{% highlight r %}
plotIntensityDistribution(intensity.counts[intensity < 250])
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-11-homopolymers/short.png) 


### Long homopolymers (>2 bases)

{% highlight r %}
plotIntensityDistribution(intensity.counts[intensity >= 250])
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-11-homopolymers/long.png) 


## Noticeably abnormal things about BAM data:

* None of flow base calls has intensity equal to 0.50, 1.50, 2.50, etc.
  In fact, that holds not only for those flow intensities from which 
  read bases were actually called, but for all of them.

* Peaks/pits at intensities 1.00, 2.00, etc.:

{% highlight r %}
plotIntensityDistribution(homs[length == 1 & intensity < 150, 
                               list(n=count),
                               by=list(intensity, base)])
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-11-homopolymers/length1.png) 

