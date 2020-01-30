# PRACTICE_EXAM_1

rm(list = ls()) ## clean workspace
load("exam1819.RData") ## change the path to your location
ls()


### 1
# 1.1

d201 <- radiation$dose == 201
d220 <- radiation$dose == 220
d243 <- radiation$dose == 243
d260 <- radiation$dose == 260
treat_stm <- radiation$treatment == 1
treat_ctrl <- radiation$treatment == 0
dead <- radiation$dead == 1
alive <- radiation$dead == 0
doses_v <- c(201, 220, 243, 260) 

sum(treat_stm)
probability <- function(dose, treatment, dead){
 return(sum(dose & treatment & dead)/sum(dose & treatment)) 
}

c201_stm <- probability(d201, treat_stm, dead)
c220_stm <- probability(d220, treat_stm, dead)
c243_stm <- probability(d243, treat_stm, dead)
c260_stm <- probability(d260, treat_stm, dead)
stm_prob <- c(c201_stm, c220_stm, c243_stm, c260_stm)

c201_ctrl <- probability(d201, treat_ctrl, dead)
c220_ctrl <- probability(d220, treat_ctrl, dead)
c243_ctrl <- probability(d243, treat_ctrl, dead)
c260_ctrl <- probability(d260, treat_ctrl, dead)
ctrl_prob <- c(c201_ctrl, c220_ctrl, c243_ctrl, c260_ctrl)

plot(doses_v, ctrl_prob, type = "l", col = 'red', ylim = c(0, 1), ylab = 'prob')
lines(doses_v, stm_prob, col = 'blue')


# 1.2

chi_sq <- function(dosi){
  table <- table(radiation[radiation$dose == dosi, c(2,3)])
  chisq.test(table)
}
chi_sq(201)
chi_sq(220)
chi_sq(243)
chi_sq(260)

table201 <- table(radiation[radiation$dose == 201, c(2, 3)])
chisq.test(table201)
# the p-value isn't smaller than the significance level of 0.05 so we don't reject 
# the null hypothesis, the 2 variables are independent

table220 <- table(radiation[radiation$dose == 220, c(2, 3)])
chisq.test(table220)
# the p-value is smaller than the significance level of 0.05 so we reject 
# the null hypothesis, the 2 variables are dependent

table243 <- table(radiation[radiation$dose == 243, c(2, 3)])
chisq.test(table243)
# the p-value is smaller than the significance level of 0.05 so we reject 
# the null hypothesis, the 2 variables are dependent

table260 <- table(radiation[radiation$dose == 260, c(2, 3)])
chisq.test(table260)
# the p-value isn't smaller than the significance level of 0.05 so we don't reject 
# the null hypothesis, the 2 variables are independent

# 1.3
# I compute the bootstrap since I'm not able to obtain a closed form expression

delta_bt <- function(dosi){ 
  radiation_bt <- radiation[sample(seq_len(nrow(radiation)), nrow(radiation), replace = T),]
      # to apply bootstrap to a dataframe
  dead_bt <- radiation_bt$dead == 1
  treat_stm_bt <- radiation_bt$treatment == 1
  treat_ctrl_bt <- radiation_bt$treatment == 0
  dose_bt <- radiation_bt$dose == dosi
 
  probabilities_stm_bt <- sum(dose_bt & treat_stm_bt & dead_bt)/sum(dose_bt & treat_stm_bt)
  probabilities_ctrl_bt <- sum(dose_bt & treat_ctrl_bt & dead_bt)/sum(dose_bt & treat_ctrl_bt)
  
  delta <- probabilities_stm_bt - probabilities_ctrl_bt
  return(delta)
}

delta_201 <- c201_stm - c201_ctrl
delta_220 <- c220_stm - c220_ctrl
delta_243 <- c243_stm - c243_ctrl
delta_260 <- c260_stm - c260_ctrl

delta_bt_201 <- replicate(1000, delta_bt(201))
delta_bt_220 <- replicate(1000, delta_bt(220))
delta_bt_243 <- replicate(1000, delta_bt(243))
delta_bt_260 <- replicate(1000, delta_bt(260))
delta_bt_v <- c(delta_bt_201, delta_bt_220, delta_bt_243, delta_bt_260)

# one-sided Wald test
alpha <- 0.05
z <- qnorm(1 - alpha)
w_201 <- delta_201/sd(delta_bt_201)
w_201 < -z

w_220 <- delta_220/sd(delta_bt_220)
w_220 < -z

w_243 <- delta_243/sd(delta_bt_243)
w_243 < -z

w_260 <- delta_260/sd(delta_bt_260)
w_260 < -z

# p-values
p_201 <- 1 - pnorm(-w_201)
p_201
p_220 <- 1 - pnorm(-w_220)
p_220
p_243 <- 1 - pnorm(-w_243)
p_243
p_260 <- 1 - pnorm(-w_260)
p_260

# 1.4

treated <- radiation[radiation$treatment == 1,]     
              # dataframe with just the streptomycin treated mice

logit_model_1 <- glm(dead ~ dose, data = treated, family = 'binomial')

n <- c(201,220,243,260)
predictions <- sapply(n,function(i){
  predict(logit_model_1, type = 'response', data.frame(dose = i))
})

plot(doses_v, predictions, type = "l", col = 'red', 
                      ylab = 'prob', ylim = c(0.05,0.8))
lines(doses_v, stm_prob, col = 'blue')
legend("topleft", legend = c('predicted', 'estimated'), col = c('red', 'blue'), 
                                  lty = 1, cex = 0.75)

logit_model_2 <- glm(dead ~ dose + I(dose^2), data = treated, family = 'binomial')
logit_model_3 <- glm(dead ~ dose + I(dose^2) + I(dose^3), data = treated, 
                     family = 'binomial')

est_c_mod_1 <- logit_model_1$coefficients
est_c_mod_2 <- logit_model_2$coefficients
est_c_mod_3 <- logit_model_3$coefficients
coefficients(logit_model_1)
coefficients(logit_model_2)
coefficients(logit_model_3)

AIC(logit_model_1, logit_model_2, logit_model_3) 
                # logit_model_2 performs better (lower AIC)
anova(logit_model_1, logit_model_2, test = 'LRT')
  # The LRT tests whether the reduction in the residual sum of squares (RSS) is 
    # statistically significant or not. The null hypothesis is that the reduction 
    # is not significant, p-value > ?? so we don't reject the null hypothesis 
    # and we should use the simpler model (logit_model_1).

### 2
# 2.1
model_fitall <- lm(fat ~ ., data = bodyfat)
summary(model_fitall)
  # the features that seem to be the most relevant to predict the 
    # percentage of body fat are the abdomen circumference, the forearm and 
      # the wrist (low p-value)
# 
# We can't reject that the coefficient for knee (knee circumference) is equal to 0
# (null hypothesis) because its p-value is extremely high.

# 2.2
model_intercept <- lm(fat ~ 1, data=bodyfat)
fw_model <- step(model_intercept, scope = formula(model_fitall), 
                 direction = "forward", trace = 0, k = log(nrow(bodyfat)))
summary(fw_model)
bw_model <- step(model_fitall, direction = "backward", trace = 0, 
                 k = log(nrow(bodyfat)))
summary(bw_model)

# the 2 models select the same covariates, so we can select just one of the 2.

# 2.3
bodyfat$BMI <- bodyfat$weight/(bodyfat$height)^2
BMI_model <- lm(fat ~ BMI + age, data = bodyfat)
BMI_model_coefficients <- summary(BMI_model)$coefficients
BMI_model_coefficients

# 2.4
# bootstrap for the se
M <- 1000
n <- nrow(bodyfat)
coefficients_bt <- replicate(M, {
  BMI_model_bt <- lm(formula(BMI_model), bodyfat[sample(1:n, replace = T),])
  return(coefficients(BMI_model_bt))
  })

intercept <- coefficients_bt[1,]
BMI <- coefficients_bt[2,]
age <- coefficients_bt[3,]

alpha <- 0.05

quantile(intercept, probs = c(alpha/2, 1 - alpha/2))
quantile(BMI, probs = c(alpha/2, 1 - alpha/2))
quantile(age, probs = c(alpha/2, 1 - alpha/2))

confint(BMI_model)

# 2.5
BIC(fw_model, BMI_model)
# Using BIC, the forward model (or the backward) is selected

# LOOCV

cv <- function(model){
  error = 0                   #T positive + T negative
  for(i in 1:nrow(bodyfat)){
    # 1- fit model with all except row i (with subset of data)   -> train the model
    less_i <- bodyfat[-i,]     
    model_i <- lm(formula(model), data = less_i)
    # 2- predict value for row i                                     -> test the model
    y_pred <- predict(model_i, newdata = bodyfat[i,])
    # 3- compare predicted with true
    error <- error + (bodyfat[i,]$fat - y_pred)^2
  }     # the list of TRUE and FALSE that I get from (coris[i,]$chd == chd_pred), if
# TRUE means that using the trained data to predict x value correspond to
# real value, when FALSE the trained data didn't work properly in predicting
# the correct value
# 4- take average
  return(error)
}

cv_fw <- cv(fw_model)

cv_bmi <- cv(BMI_model)

### 3
# 3.1
# PDF
dgumbel <- function(x, mu = 0, b = 1, log = F){
  temp <- (1/b) * exp((mu-x)/b - exp((mu - x)/b))
  if (log){
    return(log(temp))
  }
  else{
    return(temp)
  }
}

# CDF
pgumbel <- function(q, mu = 0, b = 1, log.p = FALSE){ 
  temp <- exp(- exp((mu - q)/ b)) 
  if (log.p){
    return(log(1 - temp))       # pnorm return the p(x) <= xo
  }else{                          # so it is equal to 1 - prob = 1 - temp
    return(1 - temp)
  }
}

# Quantile F
qgumbel <- function(p, mu = 0, b = 1){
  temp <- mu - b * log(-log(p))
}

# Sampling F
rgumbel <- function(n = 1, mu = 0, b = 1){
  qgumbel(runif(n), mu = mu, b = b)   # random uniform (from 0 to 1)
}

# Test
test <- rgumbel(1000)
hist(test, probability = T, breaks = 'FD')
curve(dgumbel(x), add = T, col = 'green')

# Check that dgumbel is positive
check_positive <- function(pdensityf){
  for (mu in c(-100:100)){
    for (b in c(1:10)){
      return(all(pdensityf(-1000:1000, mu = mu, b = b) >= 0))
    }  
  }
}
check_positive(dgumbel)      # OK, it is positive

# Check that dgumbel integrate to 1
check_integr.to1 <- function(densità){
  for (mu in c(-10:10)){
    for (b in c(1:10)){
      return(abs(integrate(densità, -Inf, Inf, mu = mu, b= b)$value - 1) < 1e-7)
    }  
  }
}
check_integr.to1(dgumbel)       # OK, it integrates to 1

# 3.2     # da scrivere le formuleeee!!!!!!!!!!!

y <- -digamma(1)
b_est <- (sd(wind)*sqrt(6))/pi
mu_est <- mean(wind) - (b_est*y)

# plot
hist(wind, probability = T)
curve(dgumbel(x, mu = mu_est, b = b_est), add = T, col = 'orange')

# qq plot
qq_perc <- qgumbel(ppoints(wind), mu = mu_est, b = b_est)
plot(sort(wind), qq_perc, xlab = "empirical", ylab = "theoretical", 
     main = "Q-Q plot (Gumbel)")
abline(0, 1, col = "red")
    # the Gumbel distribution fits the data

# 3.3
gumbel_mle <- function(pars, data){
  return(-sum(dgumbel(data, mu = pars[1], b = pars[2], log = TRUE))) 
  }
gumbel_parameters <- optim(fn = gumbel_mle, par = c(mu_est, b_est), 
                           data = wind)$par
gumbel_parameters
# the parameters are slightly different from those calculated with method of moments
# plot
hist(wind, probability = T)
curve(dgumbel(x, mu = mu_est, b = b_est), add = T, col = 'red')
curve(dgumbel(x, mu = gumbel_parameters[1], b = gumbel_parameters[2]), 
      add = T, col = 'blue', lty = 2)
legend('topright',legend = c('mom_density', 'mle_density'), lty = c(1,2), 
       col = c('orange', 'blue'))

# the 2 densities are slightly different, as expected from the previous calculations

# 3.4
mean_w <- mean(wind)
sd_w <- sd(wind)
hist(wind, probability = T)
curve(dnorm(x, mean = mean_w, sd = sd_w), col = 'purple', add = T)

qqnorm(wind)
qqline(wind, col='red')

# Gaussian mle
gaussian_value <- -sum(dnorm(wind, mean_w, sd_w, log = TRUE))

# Calculate AIC and BIC
gumbel_value <- gumbel_mle(pars = c(mu_est, b_est), data = wind)

k = 2
AIC_gumbel <- 2 * gumbel_value + 2 * k 
AIC_gaussian <- 2 * gaussian_value + 2 * k
BIC_gumbel <- 2 * gumbel_value + k * log(length(wind))
BIC_gaussian <- 2 * gaussian_value + k * log(length(wind))

data.frame(row.names = c("gumbel", "gauss"),
           mll = c(gumbel_value, gaussian_value),
           aic = c(AIC_gumbel, AIC_gaussian),
           bic = c(BIC_gumbel, BIC_gaussian))

# AIC and BIC select the gumbel model, as expected.

# 3.5
M <- 1000
gumbel_mle_bt <- replicate(M, {
  wind_bt <- sample(wind, size = length(wind), replace = T)
  optim(fn = gumbel_mle, par = c(mu_est, b_est), data = wind_bt)$par
})

mu_se <- sd(gumbel_mle_bt[1,])
b_se <- sd(gumbel_mle_bt[2,])

# confidence interval
alpha <- 0.05

z <- qnorm(1 - alpha/2 ) ## alpha/2 upper quantile 
a_mu <- mu_est - mu_se * z 
b_mu <- mu_est + mu_se * z 
c(a_mu, b_mu)

a_b <- b_est - b_se * z
b_b <- b_est + b_se * z
c(a_b, b_b)

# percentile confidence interval
mu_cf <- quantile(gumbel_mle_bt[1,], probs = c(alpha/2, 1- alpha/2))
mu_cf
b_cf <- quantile(gumbel_mle_bt[2,], probs = c(alpha/2, 1 - alpha/2))
b_cf
