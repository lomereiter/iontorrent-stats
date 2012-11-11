---
layout: post
title: Intensity distribution
---


{% highlight r %}
library(data.table)
library(weights)

homs <- data.table(read.table("_data/homopolymers.dat", header=T))
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


### Signal intensity distribution


{% highlight r %}
intensity.counts <- homs[, list(n=sum(length * count)), by=intensity]

wtd.hist(intensity.counts$intensity / 100,
         weight=intensity.counts$n,
         breaks=seq(0, 12, 0.05),
         xlab="Signal intensity", 
         prob=T,
         main="Intensity distribution",
         xaxp=c(0.5, 11.5, 22))
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-07-intensities/intensities.png) 


### Homopolymer length distribution


{% highlight r %}
len.dist <- homs[, list(n=sum(count)), by=length]
print(len.dist)
{% endhighlight %}



{% highlight text %}
##     length          n
##  1:      1 1138183694
##  2:      2  304151304
##  3:      3   72436829
##  4:      4   18786233
##  5:      5    5916739
##  6:      6    1613363
##  7:      7     317912
##  8:      8      75862
##  9:      9      14470
## 10:     10       1427
## 11:     11        106
{% endhighlight %}



{% highlight r %}

len.dist.lm <- lm(log(n) ~ length, len.dist)
print(summary(len.dist.lm))
{% endhighlight %}



{% highlight text %}
## 
## Call:
## lm(formula = log(n) ~ length, data = len.dist)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.3155 -0.2866 -0.0185  0.5222  0.6342 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   22.932      0.414    55.4  1.0e-12 ***
## length        -1.541      0.061   -25.3  1.1e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 0.64 on 9 degrees of freedom
## Multiple R-squared: 0.986,	Adjusted R-squared: 0.985 
## F-statistic:  638 on 1 and 9 DF,  p-value: 1.15e-09
{% endhighlight %}



{% highlight r %}

plot(log(n) ~ length, len.dist)
abline(len.dist.lm)
{% endhighlight %}

![center](/iontorrent-stats/figures/2012-11-07-intensities/lengths.png) 

