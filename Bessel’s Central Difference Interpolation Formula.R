# Bessel’s Central Difference Interpolation Formula
Bessel<-function(x, y, x0){
  h<-x[2]-x[1]
  n<-length(which(x<x0))-1
  u<-(x0-x[n+1])/h
  U<-rep(1,n+1)
  v<-u-0.5
  for(i in 1:n){
    U[i+1]<-U[i]*(v^2-((2*i-1)^2)/4)
  }
  S1<-0
  for(i in 1:n){
    del1<-diff.default(y,lag=1,differences=2*i-1)
    S1<-S1+(v*U[i]/((factorial(2*i-1)))*del1[n+2-i])
  }
  S2<-0
  for(i in 1:(n-1)){
    del2<-diff.default(y,lag=1,differences=2*i)
    S2<-S2+(U[i+1]/(factorial(2*i)))*(del2[n+1-i]+del2[n+2-i])/2
  }
  return(((y[n+1]+y[n+2])/2)+S1+S2)
}
# An example
x<-seq(.51,.57,.01)
y<-c(.5292437,.5378987,.5464641,.5549392,.5633233,.5716157, .5798158)
Bessel(x,y,0.55) 	# More appropriate than Stirling’s
# R output
[1]	 0.5633233
Bessel(x,y,0.5437)	# Stirling’s is more appropriate	
# R output
[1] 	0.5580519613
