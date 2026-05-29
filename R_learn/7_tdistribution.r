#t distribution in R

#1.pt(t,n-1)-->cumilative probability from (-infinity to t) for n-1

n=12
t=2.5
pvalue=pt(t,n-1)
print(paste("pvalue=",pvalue))

#2.qt(pvalue,n-1)--> returns the t value upto which the cumulative probability is pvalue 

n=12
pvalue=0.9
t=qt(pvalue,n-1)
print(paste("tvalue=",t))

#3. dt(t,n-1)-->probility desnity at t

t=-1.8
n=12
pdense=dt(t,n-1)
print(paste("pdense=",pdense))

#4. rt(m,n-1)--> returns m random deviates from a t distribution with n-1 degree of freedom 
 
 n=12
 m=20
 t=rt(m,n-1)
 print(t)
 
# 5. plot t distribution curve
 
 t=seq(-4,4,0.1)
 n=12
 pdense=dt(t,n-1)
 plot(t,pdense,col="red",lwd=2,type="l")
 
 #plot histogram of samples 
 n=12
 t=rt(10000,n-1)
 hist(t,breaks=30,col="green")
 
 
 
 
