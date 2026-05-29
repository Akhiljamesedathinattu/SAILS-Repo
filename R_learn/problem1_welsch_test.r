x=c(10.2,9.5,10.1,9.8,10.9,11.4,10.8,9.7,10.4)
y=c(11.0,11.2,10.1,11.4,11.7,11.2,10.8,11.6,10.9,10.9)
x_bar=mean(x)
y_bar=mean(y)
sd_x=sd(x)
#sd_x2=sd(x)*sd(x)
sd_y=sd(y)
#sd_y2=sd(y)*sd(y)
n=length(x)
m=length(y)
alpha=0.10


#welsch t test 

w = (x_bar - y_bar) / sqrt((sd_x^2 / n) + (sd_y^2 / m))

r = (((sd_x^2/n) + (sd_y^2/m))^2) /
    ((((sd_x^2/n)^2)/(n-1)) + (((sd_y^2/m)^2)/(m-1)))
    
alpha=0.10
a=alpha/2
tc=qt(1-a,r)
p_value=pt(w,r)
ci=(x_bar - y_bar)+tc*sqrt((sd_x^2 / n) + (sd_y^2 / m))
 cii=(x_bar - y_bar)-tc*sqrt((sd_x^2 / n) + (sd_y^2 / m))


print(paste(ci, "\u00B1", cii))


dependent

x=c(22,63.3,96,9.2,3.1,50,33,69,64,18.8,0,34)
y <- c(63.5,91.5,59,37.8,10.1,19.6,41,87.8,86,55,88,40)

d <- c()

for(i in 1:length(x)) {
  d[i] <- x[i] - y[i]
}

print(d)

