setwd('~/Dropbox/Dissertation_2014/DAS_Paird/DASPairedSimulation')
library(boot)

one <- read.table('test.txt', header = F, sep = '\t')
cairo_pdf('histHellDist.pdf', width = 6, height = 4, onefile = T, family = 'Lekton')
hist(one[one[,2]>0,2])
dev.off()

hista<-hist(one[one[,2]>0,2],plot = F)
length(hista$breaks[1:(length(hista$breaks)-1)])
length(hista$density)

x<-(hista$breaks[1:(length(hista$breaks)-1)]+hista$breaks[2:(length(hista$breaks))])/2
y<-hista$density

ab<-lm(y~log(x))
b<-ab$coefficients[2]
a<-ab$coefficients[1]

plot(x,(b*log(x)+a)/(a-b))
plot(x,y)
points(x,(b*log(x)+a)/(a-b),type='l')

inv.cdf<-function(p,a,b){
  exp((lambert_Wm1(-p*(a-b)*exp(b/a-1)/a)*a+a-b)/a)
}

#find top 20% and define it as DAS.
