set.seed(1234)


ipf_fit <- cvr.adaptive.ipflasso(dat, condition,family = 'binomial', type.measure="auc",
                     standardize=TRUE, alpha = 1, type.step1="sep",
                     blocks=groups,nfolds=10,ncv=10)

pl_fit <- prioritylasso(X = dat, Y = condition, family = "binomial", type.measure = "auc", 
                         blocks = groups, nfolds = 10,
                         standardize = FALSE, lambda.type = "lambda.min", block1.penalization = T)

pred_test <- ipflasso.predict(ipf_fit, dat_test)

pred_test_pl <- predict(pl_fit, newdata = dat_test, type = "response")

preds <- data.frame(ipf = pred_test$probabilitiestest, pl = pred_test_pl)
preds$avg <- (preds$ipf + preds$pl)/2

test_auc <- roc(factor(condition_3), preds$avg)

auc(test_auc)
plot.roc(test_auc, grid=0.1)
