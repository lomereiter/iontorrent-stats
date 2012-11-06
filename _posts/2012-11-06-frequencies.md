---
layout: layout
title: Frequencies of error types along the length of a read
---


{% highlight r %}

library(lattice)

counts <- read.table("_data/counts.dat", header = T)

plotErrorTypeFrequency <- function() {
    counts <- with(counts, counts[mismatches > 0 | insertions > 0 | deletions > 
        0, ])
    
    xyplot(I(mismatches/reads) + I(insertions/reads) + I(deletions/reads) ~ 
        offset, counts, key = simpleKey(c("Mismatches", "Insertions", "Deletions")), 
        pch = 19, xlab = "Offset, bp", ylab = "Error frequencies")
}

plotErrorTypeRatio <- function() {
    counts <- with(counts, counts[mismatches > 0, ])
    
    plot(I((insertions + deletions)/mismatches) ~ offset, counts, pch = 19, 
        xlab = "Offset, bp", ylab = "Indel/mismatch ratio")
}
{% endhighlight %}



{% highlight r %}
plotErrorTypeFrequency()
{% endhighlight %}

![center](figures/2012-11-06-frequencies_rmd/frequencies.png) 



{% highlight r %}
plotErrorTypeRatio()
{% endhighlight %}

![center](figures/2012-11-06-frequencies_rmd/ratio.png) 

