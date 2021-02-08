k_fullsample_1l <- function(data, 
                            data_val = NULL,
                            n_feat = 1, 
                            n_steps, 
                            nodes = 50,
                            epochs = 2000,
                            size_batch = 1,
                            internal_validation = FALSE,
                            ES = F,
                            keepBest = F,
                            view_loss = F,
                            ){
  
  n_feat <- ncol(data)
  n_sample <- nrow(data) - n_steps
  
  usr_batch <- (is.numeric(size_batch) &&
                size_batch%%1 == 0 &&
                size_batch != 1)
  
  
  if (usr_batch &&
      !(n_sample%%size_batch == 0)){
    
    # user provided size_batch, check consistency
    stop('\nProvided incompatible {size_batch}, must evenly divide sample size less lags {n_steps}.\n')
  }
  
    
  # batcher gates
  if (numbers::isPrime(n_sample)){
    warning('\n{data} and lags {n_steps} produce prime sample: ditching oldest observation.\nVary {n_steps} to avoid such behaviour.\n')
    
    # if n_sample is prime ditch oldest obs
    data <- data[-1,]
    n_sample <- nrow(data)
  }
  
  
  # compute unique prime factors
  prime_fs <- numbers::primeFactors(n_sample) %>% 
    unique()
  batch_prime <- tail(prime_fs,1)
    
  if (batch_prime>.25*n_sample){
    
    # take second last prime factor
    batch_prime <- tail(prime_fs, 2) %>% .[1]
  }
      
    
  data_lagged <- embed(x = as.matrix(data), dimension = (n_steps+1))
  
  
  # major IF
  if (is.null(data_val)){
    if (internal_validation==F){
      
      # error
      stop('\nProvide either compatible validation data {data_val} or activate {internal_validation = T}\n')
      
    }else if (internal_validation){
      # internally determine validation data,
      # by construction sizes are consistent and evenly divide
      # train and validation.
      
      #' *needs data_lagged first!*
      
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
    }

  }else if (!is.null(data_val)){
    if (nrow(data_val)<=n_steps) stop('\n{data_val} is too short wtr lags {n_steps}, adjust either of them.\n')
    
    # user provides val data, need to check
    # consistent batch/sample/validation sizes
    # raise error consistently with other options
    # if batch size is not auto, 1 or else problematic, check if it mods val and train
    
    data_val_n_sample <- nrow(data_val) - n_steps
    
    # check consistency with val&train
    if (usr_batch){
      
      # test
      data_val_gate <- data_val_n_sample%%size_batch == 0
      if (!data_val_gate) stop('\nProvided incompatible {size_batch} or {data_val}, must evenly divide sample sizes less lags {n_steps}.\n')
      
      data_train_gate <- n_sample%%size_batch == 0
      if (!data_train_gate) stop('\nProvided incompatible {size_batch} or {data}, must evenly divide sample sizes less lags {n_steps}.\n')
    }else if (size_batch == 'auto'){
      
      # test: even division with val data
      val_test <- data_val_n_sample%%batch_prime==0
      
      # test: size_batch in both 
      primes_cand <- intersect(numbers::primeFactors(n_sample),
                               numbers::primeFactors(data_val_n_sample))
      if (length(primes_cand)==0) stop('\nNo common factor in {data} and {data_val}: vary {n_steps} or activate {internal_validation = T}\n')
      
      if (!(batch_prime %in% primes_cand)) stop('\nNo common factor in {data} and {data_val}: vary {n_steps} or activate {internal_validation = T}\n')
      
    }
    
  }
  

  
}


# # main structure for gates
# 
# declare batcher fn
# 
# if is null data_va 
#   if internal_validation is F
#     error!
#   
#   else
#     gate on size_batch
#     batcher
#     
# elseif data_val not null (user provides data)
#   internal_validation <- F
#   
#   if batch_size is integer/is %in% primes(data)
#     great, cool
#   else if batch == auto
#     full run of batcher
#   else if batch == 1
#     do nothing, carry on