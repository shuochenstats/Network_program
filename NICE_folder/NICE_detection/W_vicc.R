
 library(locfdr)
 library(pracma)
 library(glasso)
 
cor.f=read.csv("cor_fvicc.csv",header=F)
str(cor.f)
 w <- locfdr(cor.f$V1)
 str(w)
 cor.fW=1-w$fdr
  windows()
 plot(cor.f$V1 ,cor.fW )


 write.csv(file="w.csv",cor.fW,row.names = F,
            col.names = F)
 
 #Glasso
 install.packages("glasso")
 
 library(glasso)
 
 install.packages("pracma")
 
 library(pracma) 
 set.seed(100)
x<-matrix(rnorm(50*20),ncol=20)
s<- var(x)
a<-glasso(s, rho=.01)

 cor.raw=(exp(2*cor.f$V1)-1)/ (exp(2*cor.f$V1)+1)
 W=squareform(cor.raw)
 for(i in 1:dim(W)[1])
 {
   W[i,i]=1
 }
 
 str(W)
 
 output=glasso(W,0.1)
 
 str(output)
 
 write.csv(file="inverse_wVICC.csv",output$wi,row.names = F,
           col.names = F)
 