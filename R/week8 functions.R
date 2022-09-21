# tune rf model, input x, y are data matrix
# ntree will be tuned!
rf_model_selection_parallel <- function(x, y, TEST_SIZE=0.3, ntrees=c(100),
                               mtrys=seq(3, ncol(x), by=1),
                               nodesizes=seq(1, 5)){
  
  # split the examples into training and test
  test_indices <- sample(1:nrow(x), size=as.integer(TEST_SIZE*nrow(x)), replace=FALSE)
  x_train <- x[-test_indices,]
  y_train <- y[-test_indices]
  x_test <- x[test_indices,]
  y_test <- y[test_indices]
  
  # grid over which we will perform the hyperparameter search:
  hparam_grid <- as.data.frame(expand.grid(ntree = ntrees,
                                           mtrys=mtrys, 
                                           nodesize=nodesizes))
  
  # to store the OOB estimates of the MSE
  oob_mses <- rep(0.0, nrow(hparam_grid))
  
  
  # tic(sprintf("Iteration begins: Iteration %i/%i, (ntree, mtry, nodesize): (%i, %.3f, %.3f)",
  #             1, nrow(hparam_grid),
  #             hparam_grid[1, 1], hparam_grid[1, 2], hparam_grid[1, 3]))
  
  support_parallel = function(rf_hparam){
    # train candidate model
    this_ntree <- hparam_grid[rf_hparam, 1]
    this_mtry <- hparam_grid[rf_hparam, 2]
    this_nodesize <- hparam_grid[rf_hparam, 3]
    
    # print(sprintf("Iteration %i/%i, (ntree, mtry, nodesize): (%i, %.3f, %.3f)",
    #               rf_hparam,nrow(hparam_grid), this_ntree, this_mtry, this_nodesize))
    
    rf <- randomForest(x_train, y_train, 
                       ntree = this_ntree,
                       mtry=this_mtry, 
                       nodesize=this_nodesize)
    
    if(class(y)=="factor"){
      # calculate OOB MSE for classification models
      return(rf$err.rate[this_ntree,1])
    } else{
      # calculate OOB MSE for classification models
      return(mse(y_train, predict(rf)))
    }
  }
  
  # perform the gridsearch
  oob_mses = foreach (hparam_idx = 1:nrow(hparam_grid), .combine=rbind) %dopar%{
    support_parallel(hparam_idx)
  }
  # toc()
  
  # select the best model (that which has the minimum OOB MSE)
  best_hparam_set <- hparam_grid[which.min(oob_mses),]
  print(best_hparam_set)
  
  # train a model on the whole training set with the selected hyperparameters
  rf_final <- randomForest(x_train, y_train,
                           ntree = best_hparam_set$ntree,
                           mtry=best_hparam_set$mtrys,
                           nodesize = best_hparam_set$nodesize,
                           importance=TRUE)
  
  # the test performance of the final model
  yhat_test <- predict(rf_final, newdata=x_test)
  
  # default hyperparmaeter model
  rf_default <- randomForest(x_train, y_train)
  yhat_test_default <- predict(rf_default, newdata=x_test)
  
  # MSEs
  test_mse <- mse(as.numeric(y_test), as.numeric(yhat_test))
  test_mse_default <- mse(as.numeric(y_test), as.numeric(yhat_test_default))
  
  cat(sprintf("Test MSE with default hyperparameters: %.6f, Test MSE with OOB-tuned hyperparameters: %.6f\n", 
              test_mse_default, test_mse))
  
  return(list(final_model = rf_final, 
              mse_grid = cbind(hparam_grid, oob_mses),
              test_result = c(test_mse_default, test_mse)))
}

get_param_parallel <- function(x, y, SUBSET_SIZE = 1000,
                               ntrees = c(200,300,400,500),
                               mtrys=seq(3, ncol(x), by=1),
                               nodesizes=c(1,3,5)){
  subset_idx = sample(nrow(x), SUBSET_SIZE)
  x <- x[subset_idx, ]
  y <- y[subset_idx]
  
  # a quite wide grid search

  rf_select_result = rf_model_selection_parallel(data.matrix(x), y,
                                                 TEST_SIZE = 0.3, 
                                                 ntrees = ntrees,
                                                 mtrys=mtrys,
                                                 nodesizes=nodesizes)
  
  rf_select_result %>% return()
}

mse <- function(y_true, y_obs) { mean((y_true-y_obs)**2.0) }
