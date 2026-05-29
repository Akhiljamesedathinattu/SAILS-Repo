#Analysis of Hieght weight data 
dat=read.csv(file="SOCR_height_weight_data.csv",header=TRUE)
print(colnames(dat))
print(dim(dat))

#convert height to cm ,weight to kgs

Height=dat$height.Inches*2.5
weight=dat$weight.Pounds*0.45

#stastical summary of height and weight 
print(summary(Height))
print(summary(weight))
mu_height=mean(Height)
sigma_height=sd(Height)
mu_weight=mean(weight)
sigma_weight=sd(weight)

#compute 2 variable 
z_height=(Height-mu_height)/sigma_height
z_weight=(weight-mu_weight)/sigma_weight

#plot the distribution 
par(mfrow=c(2,2))
hist(Height,breaks=30,col="blue")
hist(z_height,breaks=30,col="blue")
hist(weight,breaks=30,col="red")
hist(z_weight,breaks=30,col="red")

