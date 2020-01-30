load("golub_data.RData")

# what gene (column) allow us to classify different tumors?
fit <- glm(cl ~ ., family = 'binomial', data = golub_data)
summary(fit)
  # after 38 parameters the other columns are useless because they have few observations
# when we have a model like this, horizontal dataframe with few observations but many variables.

fit2 <- glm(cl ~ AFFX-HUMISGF3A/M97935_MA_at, family = 'binomial' , data = golub_data)
# with this code doesn't work because the name of columns have characters that R doesn't like
# Hot to CHANGE COLUMNS NAMES:
colnames(golub_data) <- c(paste0("G", 1:3049), "cl")

golub_data$cl <- as.factor(golub_data$cl)# to change the values of cl column in names
levels(golub_data$cl) <- c("ALL", "AML")
golub_data$cl

# if we don't want to change the names of the columns we can refer to them by index
fit2 <- glm(golub_data$cl ~ golub_data[,1], family = 'binomial')
summary(fit2)

# the fit2 model is the correct one: we need to do feature selection
      # selction of which predictors are important

# let's suppose we want to use just one column as predictor for one kind of leukemia
# I use the AIC to see which model is best to predict leukemia gene

score() <- c()
models() <- c()
for (variable in colnames(golub_data)[-3050]){
  tmp_model <- glm(paste("cl ~ ", variable), data = golub_data, family = 'binomial')
  scores <- c(scores, AIC(tmp_model))
  models <- c(models, variable)
}

df <- data.frame(scores = scores, models = models)
idx <- sort(df$scores, decreasing = F, index.return = T)   # to find the minimum
sort(df[,1], decreasing = F)
df[idx$ix,][1:5,]


n <- ncol(golub_data) - 1
res <- lapply(1:n, function(i){
  model <- glm(golub_data$cl ~ golub_data[,i], family = 'binomial')
  return(model)
})

#sapply(res, function(x) x$aic)
#ii <- which.min(sapply(res, funtion(x) x$aic))
plot(sapply(res, function(x) AIC(x)))

best_gene <- which.min(sapply(res, function(x) x$aic))

# it could take some minutes
# backward selection is not possible to do it because the full model doesn't make sense
# we could do forwars selection

model2 <- glm(cl ~ G394 + G829, family = 'binomial', data = golub_data)
summary(model2)

plot(golub_data$G394, golub_data$G829)

# lot of colinearity btw the 2 variables (visible from plot) so we use one of the 2 variables and another one

model2 <- glm(cl ~ G829 + G670, family = 'binomial', data = golub_data)
plot(golub_data$G670, golub_data$G829, col = golub_data$cl +1)
legend(legend = c(ALL, AML))

### P(cl = 1| G829, G670)