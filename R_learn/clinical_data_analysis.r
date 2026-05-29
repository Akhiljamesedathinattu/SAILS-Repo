Questions for SLE data analysis

1. Read the data into a data frame called “dat”
2. how many columns and how many patiests are there?
3. data summary of following columns : TLC, Platelet, Neutrophils, Lymphocyte, creatinine,ALT, Bilirubin, Proteins, ALP, Calcium, Albumin
4. Divide the screen into 4 and plot the histograms of : TLC, Platelets, Neutrophil, Lymphosite
5. Divide the screen into 4 and plot the histograms of : Bilirubin, Proteins, Calcium and Albumin
6. Divide the data into two subsets: subset 1 with Lupus_anticoagulant_binary=1
subset 2 with Lupus_anticoagulant_binary = 0
7. Between subset 1 and subset 2, compare the statistical summaries of Neutrophils.
Between subset 1 and subset 2, compare the box plots of Bilirubin.



dat=read.csv(file="SLE_clinical_data.csv",header=TRUE)
View(dat)
print(dim(dat))
summary(dat[, c("TLC","Platelet","Neutrophils","Lymphocyte","Creatinine","ALT","Bilirubin","Proteins","ALP","Calcium","Albumin")])

par(mfrow = c(2,2))
#TLC contain some bad reading remove them 
TLC=subset(dat$TLC,dat$TLC < 60)


hist(TLC,main = "Histogram of TLC",xlab = "TLC",xlim=c(0,50),col="blue")


min(dat$Platelet)
max(dat$Platelet)

platelet=subset(dat$Platelet,dat$Platelet <1000)

hist(platelet,main = "Histogram of Platelet",col="red",xlab = "Platelet",xlim=c(0,500))


hist(dat$Neutrophils,main = "Histogram of Neutrophils",col="lightpink",xlab = "Neutrophils")

hist(dat$Lymphocyte,main = "Histogram of Lymphocyte",col="lightgreen",xlab = "Lymphocyte")

par(mfrow = c(2,2))

Bilirubin=subset(dat$Bilirubin,dat$Bilirubin < 5)
hist(Bilirubin,main = "Histogram of Bilirubin",breaks=30,col="orange",xlab = "Bilirubin",xlim=c(0,2))
protein=subset(dat$Proteins,dat$Proteins< 11)
hist(protein,main = "Histogram of Proteins",col="yellow",xlab = "Proteins")
Calcium=subset(dat$Calcium,dat$Calcium< 11)
hist(Calcium,main = "Histogram of Calcium",col="black",xlab = "Calcium")

hist(dat$Albumin,main = "Histogram of Albumin",col="purple",xlab = "Albumin")
subset1=subset(dat, Lupus_anticoagulant_binary == 1)
subset2=subset(dat, Lupus_anticoagulant_binary == 0)
summary(subset1$Neutrophils)

summary(subset2$Neutrophils)


x11()
par(mfrow = c(1,2))

#boxplot(subset1$Bilirubin,
        #main = "Bilirubin_1",col="green",
        #xlab = "Bilirubin")

#boxplot(subset2$Bilirubin,
        #main = "Bilirubin_0",col="red",
        #xlab = "Bilirubin")
# Remove extreme outliers if needed
bil1 = subset(subset1$Bilirubin, subset1$Bilirubin < 5)
bil2 = subset(subset2$Bilirubin, subset2$Bilirubin < 5)

# Mean values


# Bar plot
barplot(bil1)
barplot(bil2)


















