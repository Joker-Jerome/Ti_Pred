library(xgboost)

load("../data/cad_1000_ahasig_en.RData")

df <- as.data.frame(feature_mtx)
n_sample <- dim(feature_mtx)[1]
set.seed(7)
train_idx <- sample(1:n_sample, 1000, replace = F)
test_idx <- setdiff(1:n_sample, train_idx)

train_labels <- df$feature_mtx[train_idx]
test_labels <- df$feature_mtx[test_idx]

new_tr <- model.matrix(~.+0, data = df[train_idx, -1], with = F)
new_ts <- model.matrix(~.+0, data = df[test_idx, -1], with = F)

dtrain <- xgb.DMatrix(data = new_tr, label = train_labels)
dtest <- xgb.DMatrix(data = new_ts, label = test_labels)

params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)

xgb1 <- xgb.train (params = params, data = dtrain, nrounds = 79, watchlist = list(val=dtest,train=dtrain), print_every_n = 10, early_stop_round = 10, maximize = F , eval_metric = "error")