#setwd("C:/Users/Yishi/Dropbox/NICE_marker/VICC/threshold");
library(locfdr)
library(stats)
library(splines)
d1=read.csv('all.csv', header=F)
d1=as.numeric(d1$V1)
zn=mylocfdr(d1, bre = 120, df = 30, pct = 0, pct0 = 1/4, nulltype = 2, type =0, plot = 1,   sw = 0, all=1)

#Only keep the mu_0 and sigma_0 given by locfdr on the overall data and calculate p0 as follows
mu_hat=zn$fp0[3,1]
sigma_hat=zn$fp0[3,2]
N0=length(d1[which(d1 <= mu_hat-sigma_hat)])
N=length(d1)
p0=N0/N/(pnorm(-1,0,1))

d1.in=read.csv('in.csv',header=F)
d1.in=as.numeric(d1.in$V1)
#zn.in=mylocfdr(d1.in, bre = 120, df = 30, pct = 0, pct0 = 1/4, nulltype = 2, type =0, plot = 1, sw = 0, all=0)
N0.in=length(d1.in[which(d1.in <= mu_hat-sigma_hat)])
N.in=length(d1.in)
p0.in=N0.in/N.in/(pnorm(-1,0,1))

d1.out=read.csv('out.csv',header=F)
d1.out=as.numeric(d1.out$V1)
#zn.out=mylocfdr(d1.out, bre = 120, df = 30, pct = 0, pct0 = 1/4, nulltype = 2, type =0, plot = 1,  sw = 0,all=0)
N0.out=length(d1.out[which(d1.out <= mu_hat-sigma_hat)])
N.out=length(d1.out)
p0.out=N0.out/N.out/(pnorm(-1,0,1))

czn=c(zn$fdr)
write.table(c(czn,p0,p0.in,p0.out),file='fdr_p0.csv',row.names = F,col.names =F)