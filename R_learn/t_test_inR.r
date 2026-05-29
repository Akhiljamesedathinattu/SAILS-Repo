#one sample t-test
x=c(183, 152, 178, 157, 194, 163, 144, 114, 178, 152, 118, 158, 172, 138)
res=t.test(x,alternative="two.sided",mu=165,conf.level=0.95)
print(res)



res1=t.test(x,alternative="less",mu=165,conf.level=0.95)
res2=t.test(x,alternative="greater",mu=165,conf.level=0.95)


#2 sample independent t_test with equal variance 
x=c(......)
y=c(......)
res=t.test(x,y,alternative="less",var.equal=TRUE,paired=FALSE,conf.level=0.95)


#2 sample welshes sample t-test
x=c()
y=c()
res=t.test(x,y,alternative="less",var.equal=FALSE,paired=FALSE,conf.level=0.95)

res1=t.test(x,y,alternative="greater",var.equal=FALSE,paired=FALSE,conf.level=0.95)
res2=t.test(x,y,alternative="two.sided",var.equal=FALSE,paired=FALSE,conf.level=0.95)

#4 paired t test

x=c(22,63.3,96,9.2,3.1,50,33,69,64,18.8,0,34)


 y=c(63.5,91.5,59,37.8,10.1,19.6,41,87.8,86,55,88,40)


res=t.test(x,y,paired=TRUE,alternative="less",conf.level=0.95)




#example 1 for one_sample t test

> x = c(12.8, 13.1, 13.9, 14.0, 12.7, 13.4, 13.6, 13.2, 14.1, 13.0, 12.9, 13.8, 13.5, 13.3, 12.6)
> res=t.test(x,alternative="two.sided",mu=13.5,conf.level=0.99)
> res

	One Sample t-test

data:  x
t = -2.2822, df = 15, p-value = 0.01874
alternative hypothesis: true mean is less than 60
95 percent confidence interval:
     -Inf 58.82618
sample estimates:
mean of x 
  54.9375 

