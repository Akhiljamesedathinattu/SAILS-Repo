#the following are the systolic blood preausre of 14 patient undergoing some drug thrapy 
populationcan we conclude the data of population is different from 165? let alpha be 0.05 




 x_data=c(14.5,12.9,14.0,16.1,12.0,17.5,14.1,12.9,17.9,12.0,16.4,24.2,12.2,14.4,17.0,10.0,18.5,20.8,16.2,14.9,19.6,22.3,17.8,12.1)
> l=length(x_data)
> l
[1] 14
> x_mean=mean(x_data)
> sd_x_data=sd(x_data)
> x_mean
[1] 157.2143
> sd_x_data
[1] 23.43298
> a1=x_mean-165
> a2=
> a1=x_mean-165
> a2=sd_x_data/l
> t_value=a1/a2
> t_value
[1] -4.651564
> a2=sd_x_data/sqrt(l)
> t_value=a1/a2
> t_value
[1] -1.243183
> p_value=qt(t_value,l-1)
Warning message:
In qt(t_value, l - 1) : NaNs produced
> p_value=qt(t_value,13)
Warning message:
In qt(t_value, 13) : NaNs produced
> p_value=pt(t_value,l-1)
> p_value
[1] 0.1178774
> tc=qt(1-(0.05/2),l-1)
> tc
[1] 2.160369
> ci=
> 
> ci=x_mean+(tc*(sd_x_data/sqrt(l))
+ 
> ci=(x_mean+(tc*(sd_x_data/sqrt(l)))
+ ci=(x_mean+(tc*(sd_x_data/sqrt(l))))
Error: unexpected symbol in:
"ci=(x_mean+(tc*(sd_x_data/sqrt(l)))
ci"
> ci=(x_mean+(tc*(sd_x_data/sqrt(l)))
+ 
> ci=x_mean+(tc*(sd_x_data/sqrt(l)))
> ci_neg=x_mean-(tc*(sd_x_data/sqrt(l)))
> ci
[1] 170.7441
> ci_neg
[1] 143.6845
> print(paste(ci,"+ or -",ci_neg)
+ print(paste(ci,"+ or -",ci_neg))
Error: unexpected symbol in:
"print(paste(ci,"+ or -",ci_neg)
print"
> 
> print(paste(ci,ci_neg))
[1] "170.744083633215 143.684487795356"
> print(ci,ci_neg)
Error in print.default(ci, ci_neg) : invalid printing digits 143
> print(paste(ci,ci_neg))
[1] "170.744083633215 143.684487795356"
> 
 # a condition with >165
 
 
 
 
 
 > dat=c(14.5,12.9,14.0,16.1,12.0,17.5,14.1,12.9,17.9,12.0,16.4,24.2,12.2,14.4,17.0,10.0,18.5,20.8,16.2,14.9,19.6,22.3,17.8,12.1)
> mean(dat)
[1] 15.84583
> sd(dat)
[1] 3.547317
> a1=mean(dat)-14
> a2=sd(dat)/sqrt(length(dat))
> t_value=a1/a2
> t_value
[1] 2.549166
> p_value=qt(t_value,length(dat)-1)
Warning message:
In qt(t_value, length(dat) - 1) : NaNs produced
> l=length(dat)
> l
[1] 24
> p_value=qt(t_value,23)
Warning message:
In qt(t_value, 23) : NaNs produced
> p_value=pt(t_value,23)
> p_value
[1] 0.991034
> tc=qt((1-0.05),23)
> tc
[1] 1.713872
> 1-p_value
[1] 0.008966007


