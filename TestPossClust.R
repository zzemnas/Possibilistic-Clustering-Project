
###################APPLICATION OF POSSIBILISTIC CLUSTERING TO THE X12 AND X16 DATA
#The methods applied here are proposed in the paper "A Possibilistic Approach to Clustering#
#written by Krishnapuram and Keller in 1993

library(ppclust)
library(fclust)
#?x12
data(x12) 

#proto is basically the initialization of the matrices U (membership degree)
#As suggested by the authors of the paper, we will use the Fuzzy K-means algorithm
#in order to have good starting points for the possibilistic clustering algorithm, 
#whose outcome heavily depends on them.
#Here we will use the FKM function, which is part of the fclust package.

#proto<-FKmeans(x12,k=2,m=2,RS=60,conv=1e-5,maxit=1e-6)
proto<-FKM(x12,k=2,m=2,RS=60,conv=1e-5,maxit=1e-6)

#possclust is the core function of the possibilistic clustering algorithm. (see EtaFunz for more details)
#It works in synergy with the FKM function, whose output is the proto variable.
myp<-possclust(X=x12, k=2, m=2, conv=1e-6, maxit=1e+9,RS=60,proto=proto)

#plot of x12 points.
plot(x12)

#points of the centroid, highlighted in red.
points(myp$H[1,1], myp$H[1,2],col="red",pch=16)
points(myp$H[2,1], myp$H[2,2],col="red",pch=16)

#application to x16 dataset
y<-as.matrix(x16[,1:2])
proto2<-FKM(y, k=2, m=2, conv=1e-6, maxit=1e+7,RS=60)
myp2<-possclust(y, k=2, m=2, conv=1e-6, maxit=1e+6, proto=proto2,RS=60)

plot(y)
points(myp2$H[1,1], myp2$H[1,2],col="red",pch=16)
points(myp2$H[2,1], myp2$H[2,2],col="red",pch=16)


