### Cross Validation Practice

coris <- read.table("coris.dat", skip = 4, sep = ",",
                    col.names = c("row.names", "sbp", "tobacco",
                                  "ldl", "adiposity",
                                  "famhist", "typea", "obesity",
                                  "alcohol",
                                  "age", "chd"))[,-1]     
                                    # I take all rows and all columns except the first
class(coris$chd)
fit_all <- glm(chd ~ ., data = coris, family ='binomial')

# cross validation: 1- fit model with subset of data
#                   2- predict on new data
#                   3- compare predicted with true
#                   4- take average

# when you do predict depending on the tipe argument you choose

predict(fit_all, newdata = coris[1,]) > 0
predict(fit_all, newdata = coris[1,], type = 'response') > 0.5
fit <- glm()

# if I add as.numeric it is gonna give 1 if TRUE and 0 if FALSE
as.numeric(predict(fit_all, newdata = coris[1,]) > 0)
as.numeric(predict(fit_all, newdata = coris[1,], type = 'response') > 0.5)

curve(exp(x)/(1+ exp(x)), from = -10, to = 10)

### LOOCV - Leave One Out Cross Validation
correct = 0                   #T positive + T negative
for(i in 1:nrow(coris)){
  # 1- fit model with all except row i (with subset of data)   -> train the model
  coris_i <- coris[-i,]     
  model_i <- glm(chd ~ ., data = coris_i, family = 'binomial')
  # 2- predict value for row i                                     -> test the model
  chd_pred <- as.numeric(predict(model_i, newdata = coris[i,]) > 0)
  chd_pred <- as.numeric(predict(model_i, type = "response", newdata = coris[i,]) > 0.5)
  # 3- compare predicted with true
  correct <- correct + (coris[i,]$chd == chd_pred)
}     # the list of TRUE and FALSE that I get from (coris[i,]$chd == chd_pred), if
      # TRUE means that using the trained data to predict patient x value correspond to
      # real value, when FALSE the trained data didn't work properly in predicting
      # the correct value
# 4- take average
correct/nrow(coris)
# 


### K-FOLD Cross Validation

k <- 10
folds <- list()
m <- nrow(coris) %/% k
for (i in 1:k){
  folds[[i]] <-  ((i-1)*m + 1) : (i * m)
}

folds[[1]]
folds[[2]]
folds[[3]]
folds[[4]]

coris_group <- coris[-folds[[i]],]

tot <- c()                   #T positive + T negative
for(i in 1:k){
  # 1- fit model with all but i group (with subset of data)
  model_group <- glm(chd ~ ., data = coris_group, family = 'binomial')
  # 2- predict on group i
  chd_pred_group <- as.numeric(predict(model_group, newdata = coris[folds[[i]],]) > 0)
  # chd_pred_group <- as.numeric(predict(model_group, type = "response", newdata = coris[i,]) > 0.5)
  # 3- compare predicted with true
  correct_group <- chd_pred_group == coris[folds[[i]], 'chd']
  tot[i] <- sum(correct_group)
}     

mean(tot/m)

# if I want to add the shuffle i write before the function 
# new_coris <- coris[sample(1:nrow(coris)),]