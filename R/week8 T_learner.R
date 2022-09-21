# T-learner ---------------------------------------------------
# to get the CATE we simply run the function T_learner

# DOCUMENT HERE
#' Remove duplicated strings
#'
#' `T_learner()` computes CATE for given treatment and control data
#'
#' @param x0 data.frame: control data for predictive variables
#' @param y0 vector: control response variable
#' @param x1 data.frame: treatment data for predictive variables
#' @param y1 vector: treatment response variable
#' @param k_fold int: number of cross-validation
#' @param hparam0;hparam1 list of hyper-parameter if pre-trained value available
#' @returns CATE: a data.frame with each folumn for one fold data.
#' @examples
#' T_learner(x0, y0, x1, y1, hparam0=param_0, hparam1=param_1)
#' @export


T_learner = function(x0, y0, x1, y1, # input data X, y for control and treatment
                     k_fold=5,
                     hparam0=NULL, hparam1=NULL){
  # find best hyper-parameter if it's not pretrained
  if(is.null(hparam0)){
    rf_select_result_0 = get_param_parallel(x0, y0, SUBSET_SIZE = floor(nrow(x0)/10),
                                            ntrees = c(100, 200, 300,400,500),
                                            mtrys=seq(20, ncol(x0), by=1),
                                            nodesizes=seq(1,5))
    param_0 = rf_select_result_0$mse_grid %>% arrange() %>% .[1,] %>% return()
  } else{
    param_0 = hparam0
  }
  if(is.null(hparam1)){
    rf_select_result_1 = get_param_parallel(x1, y1, SUBSET_SIZE = floor(nrow(x1)/10),
                                            ntrees = c(100, 200, 300, 400, 500),
                                            mtrys=seq(20, ncol(x1), by=1),
                                            nodesizes=seq(1,5))
    param_1 = rf_select_result_1$mse_grid %>% arrange() %>% .[1,] %>% return()
  } else{
    param_1 = hparam1
  }

  # k-fold CV
  CATE = foreach(k=1:k_fold, .combine = rbind) %dopar% {{
    # split x0, y0
    test_idx0 <- sample(1:nrow(x0), size=as.integer(nrow(x0)/k_fold), replace=FALSE)
    x0_train <- x0[-test_idx0,]
    y0_train <- y0[-test_idx0]
    x0_test <- x0[test_idx0,]
    y0_test <- y0[test_idx0]
    # fit M0
    tic()
    model.0 <- randomForest(data.matrix(x0_train), y0_train,
                            ntree=param_0$ntree,
                            mtry=param_0$mtrys,
                            nodesize = param_0$nodesize)
    toc() 
    
    # split x1, y1
    test_idx1 <- sample(1:nrow(x1), size=as.integer(nrow(x1)/k_fold), replace=FALSE)
    x1_train <- x1[-test_idx1,]
    y1_train <- y1[-test_idx1]
    x1_test <- x1[test_idx1,]
    y1_test <- y1[test_idx1]
    # fit M1
    tic()
    model.1 <- randomForest(data.matrix(x1_train), y1_train,
                            ntree=param_1$ntree,
                            mtry=param_1$mtrys,
                            nodesize = param_1$nodesize)
    toc() 
    
    # compute CATE with train data
    # x_test <- rbind(x0_test, x1_test)
    
    compute_CATE <- function(x_data){
      y0_star = predict(model.0, newdata=x_data)
      y1_star = predict(model.1, newdata=x_data)
      return(mean(y1_star - y0_star))
    }
    # we have four inputs: x0_train, x1_train, x0_test, x1_test
    CATE0_train = compute_CATE(x0_train)
    CATE0_test = compute_CATE(x0_test)
    CATE1_train = compute_CATE(x1_train)
    CATE1_test = compute_CATE(x1_test)
  
    
    
    list(CATE0_train=CATE0_train,
         CATE0_test=CATE0_test,
         CATE1_train=CATE1_train,
         CATE1_test=CATE1_test)
  }
  }
  return(CATE)
}
