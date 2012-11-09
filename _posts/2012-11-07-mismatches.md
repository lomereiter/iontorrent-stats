---
layout: post
title: Mismatches
---


{% highlight r %}

library(lattice)
library(data.table)
load("~/Dropbox/Public/iontorrent-stats/data/mismatches.RData")

print(head(mismatches))
{% endhighlight %}



{% highlight text %}
##    pos strand ref.nuc nuc intensity offset
## 1: 278      +       A   G      1.08    282
## 2: 279      +       G   A      0.84    283
## 3:  55      +       G   A      1.80     87
## 4: 169      +       T   A      5.96    242
## 5: 171      +       A   T      2.20    244
## 6:  53      +       G   A      1.54     69
{% endhighlight %}



{% highlight r %}

homs <- data.table(read.table("_data/homopolymers.dat", header = T))

intensity.counts <- homs[, list(total = sum(length * count)), by = intensity]
setkey(intensity.counts, "intensity")
{% endhighlight %}


### Mismatch frequency as a function of intensity


{% highlight r %}
mismatch.counts <- mismatches[, list(mismatches = nrow(.SD)), by = intensity]
mismatch.counts$intensity <- mismatch.counts$intensity * 100

xyplot(I(mismatches/total) ~ I(intensity/100), intensity.counts[mismatch.counts][intensity < 
    550], pch = 19, col = "#00DDFF", scales = list(x = list(at = seq(0.5, 5.5, 
    0.5))), xlab = "Signal intensity", ylab = "Mismatch frequency")
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-07-mismatches/mismatchfreq.png) 



{% highlight r %}
xyplot(I(log(mismatches/total)) ~ I(intensity/100), intensity.counts[mismatch.counts], 
    pch = 19, col = "#00DDFF", scales = list(x = list(at = seq(0.5, 11.5, 0.5))), 
    xlab = "Signal intensity", ylab = expression("ln" ~ ~bgroup("(", "mismatch frequency", 
        ")")))
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-07-mismatches/logmismatchfreq.png) 

