library(glmnet)
getwd()
x <- read.table(file = "expression3-Kim.txt", header = TRUE, sep = "",row.names = 1)
x= as.matrix(x)
y <- read.table(file = "Outcome-Kim.txt", header = TRUE, sep = "",row.names = 1)
y <- y$Response
fit <- glmnet(x, y, family="binomial", nlambda=100, alpha=1)
plot(fit, xvar="lambda", label=TRUE)
cv.fit <- cv.glmnet(x, y, family = "binomial",nfolds = 10)
plot(cv.fit)
min = cv.fit$lambda.min 
se = cv.fit$lambda.1se
coef(cv.fit$glmnet.fit,s=min,exact = F)

