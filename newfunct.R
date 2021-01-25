k_fullsample_1l <- function(data, 
                            n_steps, 
                            n_feat = 1, 
                            # model_compiled,
                            nodes = 50,
                            size_batch = 1,
                            epochs = 2000,
                            ES = F,
                            keepBest = F,
                            view_loss = F,
                            internal_validation = FALSE,
                            data_val = NULL){
  
  # major IF
  if (is.null(data_val)){
    # do something when user does not provide data val
    
    #' *when n_sample is prime ditch oldest obs!*
    
    if (sad){ 
    }
    
  }else{
    
    # user provides val data, need to check
    # consistent batch/sample/validation sizes
    # raise error consistently with other options
    internal_validation <- FALSE
    
    if (size_batch == 'auto'){
      # do our thingy
    }else if (size_batch == 1){
      # easy green light, go ahead
    }else if (is.numeric(size_batch) && (size_batch%%1==0)){
      # check consistency 
      
    }
    
  }
  
  # batcher
  
  
}


y <- 3
t2 <- T
if (is.null(y) && t2){
  print('went tru')
}else{
  print('failed')
}


if (is.null(data_val)){
  if (internal_validation){
    
    # n_sample to account val
    n_train <- n_sample - batch_prime
    
    # this ensures that mod(val_sample)==mod(n_sample)
    val_sample <- tail(data_lagged, batch_prime)  
    train_sample <- data_lagged[-((n_train + 1):n_sample),]
    
    
    
    # recast in arrays:
    # training target - 2d
    y_data_arr <- array(data = train_sample[, 1],
                        dim = c(n_train, n_feat))
    # training features - 3d
    x_data_arr <- array(data = train_sample[, -1],
                        dim = c(n_train, n_steps, n_feat))
    
    # validation target
    y_val_arr <- array(data = val_sample[, 1],
                       dim = c(batch_prime, n_feat))
    x_val_arr <- array(data = val_sample[, -1],
                       dim = c(batch_prime, n_steps, n_feat))
    
    # fine for univariate case
  }
}
