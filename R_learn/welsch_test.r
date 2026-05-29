x=c(131,115,124,131,122,117,88,114,150,169)
y=c(60,150,130,180,163,130,121,119,130,148)
x_bar=mean(x)
y_bar=mean(y)
sd_x=sd(x)
#sd_x2=sd(x)*sd(x)
sd_y=sd(y)
#sd_y2=sd(y)*sd(y)
n=length(x)
m=length(y)


#welsch t test 

w = (x_bar - y_bar) / sqrt((sd_x^2 / n) + (sd_y^2 / m))

r = (((sd_x^2/n) + (sd_y^2/m))^2) /
    ((((sd_x^2/n)^2)/(n-1)) + (((sd_y^2/m)^2)/(m-1)))
    
alpha=0.10
a=aplha/2
tc=qt(1-a,r)
p_value=pt(w,r)


#very imp
 bp = c(183, 152, 178, 157, 194, 163, 144, 114, 178, 152, 118, 158, 172, 138)
> res=t.test(bp,alternative="two.sided",mu=165,conf.level=0.95)
> res

	One Sample t-test

data:  bp
t = -1.2432, df = 13, p-value = 0.2358
alternative hypothesis: true mean is not equal to 165
95 percent confidence interval:
 143.6845 170.7441
sample estimates:
mean of x 
 157.2143 
 
 
 
 res=t.test(bp,alternative="less",mu=165,conf.level=0.95)
> res

	One Sample res=t.test(x,y,paired=TRUE,alternative="less",conf.level=0.95)t-test

data:  bp
t = -1.2432, df = 13, p-value = 0.1179
alternative hypothesis: true mean is less than 165
95 percent confidence interval:
     -Inf 168.3052
sample estimates:
mean of x 
 157.2143 












