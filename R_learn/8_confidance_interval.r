


#getting a confidence interval for the given data

#dataset of biopsy samples
x = c(1009, 1280, 1180, 1255, 1547, 2352, 1956, 1080, 1776, 1767, 1680, 2050, 1452, 2857, 3100, 1621)
#size of dataset
n = length(x)
alpha = 0.05
conf = qt((1-(alpha/2)), n-1)*sd(x)/sqrt(n)
CI_upper = mean(x) + conf
CI_lower = mean(x) - conf
CI = 1 - alpha
print(paste((CI)*100,"%","CI =(",round(CI_lower), ",", round(CI_upper),")"))
