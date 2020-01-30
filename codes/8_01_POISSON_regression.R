crab <- read.table('crab.txt', header = T)
head(crab)
dim(crab)
model_Wt <- glm(Sa ~ Wt, data = crab, family = 'poisson')
summary(model_Wt)

# we have a very strong evidence to reject the nul hypothesis that coefficient for the weight is = 0
  # the weight has a great influence on the number of satellites
      # H0 : beta1 = 0 (Estimate Wt)


# Now we need to categorize the colors as factors so that they are not evaluated for
    # how close they are (1-2 closer than 1-3), but as independent unique things.
crab1 <- crab
crab1$C <- as.factor(crab$C)
levels(crab1$C)
levels(crab1$C) <- c('blue', 'red', 'yellow', 'green')
model_c <- glm(Sa ~ C, data = crab1, family = 'poisson')
summary(model_c)
summary(crab1$C)
# blue is the base, we should put as base the most common color (summary(crab1$C)) -> red 
# yellow has positive influence with respect of blue
  # green same

# Spine condition is also convenient to put it as factor
crab2 <- crab2
crab2$S <- as.factor(crab$S)
levels(crab2$S)
levels(crab2$C) <- c('')

# not part of the course but if we want to take into account 2 variables together
    # to see how they influence the data, we do:
model <- glm(Sa ~ C * S, family = 'poisson', data = crab1)      # the * is not a product
                                                           # for the product we would need I(C * S)