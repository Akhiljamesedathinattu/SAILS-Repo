x=c(131,115,124,131,122,117,88,114,150,169)
y=c(60,150,130,180,163,130,121,119,130,148)
x_bar=mean(x)
y_bar=mean(y)
sd_x=sd(x)
sd_x2=sd(x)*sd(x)
sd_y=sd(y)
sd_y2=sd(y)*sd(y)
n=length(x)
m=length(y)

# for t test 


a1=x_bar-y_bar

a111=9*sd_x2
a112=9*sd_y2
a113=10+10-2
a114=(a111+a112)/a113
a115=sqrt(a114)

a12=1/m+1/n

a14=sqrt(a12)
a15=a115*a14
t_value=a1/a15

tc=qt(1-0.05,n+m-2)
p_value=pt(t_value,n+m-2)

t_value = (x_bar - y_bar) / 
          sqrt( (((n-1)*sd_x^2 + (m-1)*sd_y^2) / (n+m-2)) * ((1/n) + (1/m)) )





# VERY IMP





#taking a data set
x=c(183,152,178,157,194,163,144,114,178,152,118,158,172,138)

#or x=c(14.5,12.9,14.0,16.1,12.0,17.5,14.1,12.9,17.9,12.0,16.4,24.2,12.2,14.4,17.0,10.0,18.5,16.2,14.9,19.6,22.3,17.8,12.1,20.8)


#finding mean
xbar=mean(x)
#[1] 157.2143

#standard deviation
s=sd(x)
#[1] 23.43298

#t test
t=(xbar-165)/(s/sqrt(length(x)))
print(t)
#[1] -1.243183

#length
n=14

#pval if the t value is negtive then do not use 1- use 1- when t value is postive
pval=1-pt(t,n-1)
print (pval)

pval=pt(t,n-1)
print (pval)
#[1] 0.1178774 for this data set t is negtive so we are useing the second formula


#giving alpha value
aplha=0.05

#if statment to reject and accept
if (pval<aplha/2) {
  print("accept null")
} else {
  print("reject null")
}


#finding t critical when ur taking both side alpha/2 or  if ur taking for 1 side take alpha only
tc=qt(1-(aplha/2),n-1)
print (tc)
#[1] 2.160369

#finding confidence intervel

ci=tc*s/sqrt(length(x))
print(ci)
#[1] 13.5298

posconi=xbar+ci
negconi=xbar-ci

print(posconi)
#[1] 170.7441
print(negconi)
#[1] 143.6845

#y is the valuetaken as a generalized mean
y=165

# now checking if y is in the range of posconi and negconi so that we can accpet or reject them
if (y > negconi & y < posconi) {
  print("accept null")
} else {
  print("reject null")
}
#[1] "accept null"

#res
res=t.test(x,alternative="two.sided",mu=165,conf.level=0.95)
print(res)
