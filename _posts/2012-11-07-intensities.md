---
layout: post
title: Intensity distribution
---


{% highlight r %}
library(data.table)
library(weights)

homs <- data.table(read.table("_data/homopolymers.dat", header = T))
print(head(homs))
{% endhighlight %}



{% highlight text %}
##    base length intensity  count
## 1:    A      1        51 411172
## 2:    A      1        52 240347
## 3:    A      1        53 129832
## 4:    A      1        54 134885
## 5:    A      1        55 146169
## 6:    A      1        56 162345
{% endhighlight %}



{% highlight r %}
intensity.counts <- homs[, list(n = sum(length * count)), by = intensity]

wtd.hist(intensity.counts$intensity/100, weight = intensity.counts$n, breaks = seq(0, 
    12, 0.05), xlab = "Signal intensity", prob = T, main = "Intensity distribution", 
    xaxp = c(0.5, 11.5, 22))
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-07-intensities/intensities.png) 

