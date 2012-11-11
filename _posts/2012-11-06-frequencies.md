---
layout: post
title: Frequencies of error types along the length of a read
---


{% highlight r %}
library(lattice)
library(data.table)

counts <- data.table(read.table("_data/counts.dat", header=T))
{% endhighlight %}



{% highlight r %}
xyplot(I(mismatches/reads) + I(insertions/reads) + I(deletions/reads) ~ offset, 
       counts[mismatches > 0 | insertions > 0 | deletions > 0], 
       key=simpleKey(c("Mismatches", "Insertions", "Deletions")), 
       pch=19, 
       xlab="Offset, bp",
       ylab="Error frequencies")
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-06-frequencies/frequencies.png) 



{% highlight r %}
xyplot(I((insertions + deletions) / mismatches) ~ offset,
       counts[mismatches > 0],
       pch=19,
       xlab="Offset, bp",
       ylab="Indel/mismatch ratio")
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-06-frequencies/ratio.png) 

