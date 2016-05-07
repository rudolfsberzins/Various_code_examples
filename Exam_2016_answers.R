#### Exam 2016 Answers by Rudolfs Berzins (ckj464) ####

#### Problem 1 Question 1 ####

hist(data, main = "Histogram of dataset named data")
d <- density(data)
plot(d, main = "Density plot of dataset named data")

#### Problem 1 Question 2 ####

plot(ecdf(data), main = "")

#### Problem 1 Question 3 ####

quantile(data, probs = c(0.25, 0.5, 0.75))

#### Problem 1 Question 4 ####

plot(ecdf(data))
curve(quantile(data,x), from = 0, to = 1, main = "Empirical Quantile Function plot")

#### Problem 1 Question 5 ####
# pp-plot
plot (x=pnorm((data), mean = mean(data), sd = sd(data)), y=ecdf(data)(data), main = "PP-plot of data", xlab = "Theoretical N(mu,sigma) distribution", ylab = "Emperical distribution function")

#### Problem 1 Question 6 ####

# qq-plot
qqnorm(data)

#### Problem 1 Question 7 ####

loglike <- function(x) -sum(log(dt(data, df = x)))

#### Problem 1 Question 8 ####

a <- c(1:10)
plot(sapply(a,loglike), main = "Plot of a minus log likelihood function") # It's 2

#### Problem 1 Question 9 ####

plot (x = qt((1:length(data)-0.5)/length(data), df = 2), y = sort(data), main = 'QQ-plot of data', xlab = "Theoretical N(mu, sigma) quantiles", ylab = "Sample Quantiles")
plot (x = pt(data, df = 2), y=ecdf(data)(data), main = "PP-plot of data", xlab = "Theoretical N(mu,sigma) distribution", ylab = "Emperical distribution function")

#### Problem 1 Question 10 ####

#### Problem 2 Question 1 ####

a <- subset(diabetes, diabetes=='TRUE', select = diabetes)
nrow(a)

tru <- subset(diabetes, diabetes=='TRUE', select = waist.hip)
tru_mean <- mean(tru$waist.hip)
tru_mean
fal <- subset(diabetes, diabetes=='FALSE', select = waist.hip)
fal_mean <- mean(fal$waist.hip)
fal_mean

#### Problem 2 Question 2 ####

v <- var(tru$waist.hip)
sem <- sqrt(v/length(tru$waist.hip))
t.975 <- qt(0.975, df=length(tru$waist.hip)-1)
c(tru_mean-t.975*sem, tru_mean+t.975*sem)

#### Problem 2 Question 3 ####

v <- var(fal$waist.hip)
sem <- sqrt(v/length(fal$waist.hip))
t.975 <- qt(0.975, df=length(fal$waist.hip)-1)
c(fal_mean-t.975*sem, fal_mean+t.975*sem)

#### Problem 2 Question 4 ####

# Based on the calculated data, we can determin that the waist.hip 
# ratio is different between the subjects with diagnosis of diabetes
# and without diagnosis of diabetes.

#### Problem 2 Question 5 ####

R <- function(x) (nrow(subset(diabetes, diabetes == 'TRUE' & waist.hip>=x)))/(nrow(subset(diabetes, waist.hip>=x)))

R(1)
  

#### Problem 2 Question 6 ####

x <- c(0.8,0.85,0.9,0.95,1.0,1.05,1.1) 

plot(sapply(x,R),x, main = "Function R(x) against x")

# The value number 1.0 jumps out

#### Problem 2 Question 7 ####

fit <- glm(diabetes ~ waist.hip, family = binomial('logit'), data = diabetes) 
summary(fit)

#### Problem 2 Question 8 ####

fit1 <- glm(diabetes ~ age + gender + factor(age.group):waist.hip, family = binomial('logit'), data = diabetes)
summary(fit1)

#### Problem 2 Question 9 ####

#### Problem 3 Question 1 ####

diet_1 <- subset(cw, Time==10.00 & Diet == 1, select = weight)
mean_diet_1 <- mean(diet_1$weight)
mean_diet_1

#### Problem 3 Question 2 ####

diet_2 <- subset(cw, Time==10.00 & Diet == 2, select = weight)
mean_diet_2 <- mean(diet_2$weight)
mean_diet_2

#### Problem 3 Question 3 ####

boot.samples <- replicate(10000, sample(cw$weight[cw$Time==10.00 & cw$Diet==1], replace = T), simplify = F)
mean <- sapply(boot.samples, mean)
sem <- sd(mean)
z <- qnorm(0.95)
mean_diet_1 + c(-z*sem, z*sem)

#### Problem 3 Question 4 ####

boot.samples <- replicate(10000, sample(cw$weight[cw$Time==10.00 & cw$Diet==2], replace = T), simplify = F)
mean <- sapply(boot.samples, mean)
sem <- sd(mean)
z <- qnorm(0.95)
mean_diet_2 + c(-z*sem, z*sem)

#### Problem 3 Question 5 ####

mean_diet_1
mean_diet_2

# The chicks raised with diet 2 will have a larger mean weight at
# day 10 than the chicks raised with diet 1. ....

#### Problem 3 Question 6 ####

k <- function(i) min(subset(cw, weight>=100 & Chick == i, select = Time)) 

#### Problem 3 Question 7 ####

chicks_on_1 <- subset(cw, Diet == 1 & Time == 0, select=Chick)

Ek1 <- mean(sapply(chicks_on_1$Chick, k))

#### Problem 3 Question 8 ####

chicks_on_2 <- subset(cw, Diet == 2 & Time == 0, select=Chick)

Ek2 <- mean(sapply(chicks_on_2$Chick, k))

#### Problem 3 Question 9 ####

boot.samples <- replicate(10000, sample(cw$Chick[cw$Time==0 & cw$Diet==1], replace = T), simplify = F)
mean <- sapply(boot.samples, mean)
sem <- sd(mean)
z <- qnorm(0.95)
Ek1 + c(-z*sem, z*sem)

#### Problem 3 Question 10 ####

boot.samples <- replicate(10000, sample(cw$Chick[cw$Time==0 & cw$Diet==2], replace = T), simplify = F)
mean <- sapply(boot.samples, mean)
sem <- sd(mean)
z <- qnorm(0.95)
Ek2 + c(-z*sem, z*sem)

#### Problem 3 Question 11 ####

# It takes around 2 days more to get the chicks to reach weight of 100 
# on diet 1 than on diet 2

#### Problem 4 Question 2 ####

curve(pexp(x^2,rate=1),type='l', from = 0, to = 3, ylab = 'Fn(y)', xlab = 'y', main = "Exponential distribution function of FY(y)")


