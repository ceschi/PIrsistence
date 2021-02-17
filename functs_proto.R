lstm_1l <- function(data, 
                    data_val = NULL,
                    n_feat = 1, 
                    n_steps, 
                    nodes = 50,
                    epochs = 2000,
                    size_batch = 1,
                    internal_validation = TRUE,
                    ES = F,
                    keepBest = F,
                    view_loss = F){
  
  # Function to fit a model on the whole sample of data;
  # it takes care of lagging & reshaping the data according to parameters
  # passed and to declare a model with one layer of LSTM and a final
  # dense layer. It requires 'keras' and pipes. It preserves the 
  # time dimension dropping the NAs generated after lagging.
  # Optionally it accommodates the early stopping callback
  # for which the epochs are then the patience limit.
  # Output contains history, model, and optionally time index.
  
  # Keras usually run in parallel by default, call this function 
  # sequentially to avoid nested parallelism.
  suppressWarnings(require(magrittr))
  suppressWarnings(require(keras))
  suppressWarnings(require(tictoc))
  suppressWarnings(require(numbers))
  
  # data come in as a simple TS,
  # it comes then with n of observations (n_sample)
  # and must be lagged according to n_steps.

  # it will drop as many obs as many lags we include
  n_sample <- nrow(data) - n_steps
  n_feat <- ncol(data)
  
  # preserve the time index of data
  # for later use
  if (is.xts(data)) time_index <- time(data)[(n_steps+1):length(data)]
  
  # shorthand condition
  usr_batch <- (is.numeric(size_batch) &&
                  size_batch%%1 == 0 &&
                  size_batch != 1)
  
  # when the sample is prime raise err
  if (numbers::isPrime(n_sample)){
    warning('\n{data} and lags {n_steps} produce prime effective sample: ditching oldest observation.\nVary {n_steps} to avoid such behaviour.\n')
    
    # if n_sample is prime ditch oldest obs,
    # preserve og data
    data_og <- data
    data <- data[-1,, drop = F]
    n_sample <- nrow(data) - n_steps
  }
  
  if (usr_batch &&
      !(n_sample%%size_batch == 0)){
    
    # user provided size_batch, check consistency
    stop('\nProvided incompatible {size_batch}, must evenly divide sample size less lags {n_steps}.\n')
  }
  
  
  # compute unique prime factors
  prime_fs <- numbers::primeFactors(n_sample) %>% 
    unique()
  batch_prime <- tail(prime_fs,1)
  
  # check batch too big
  if (batch_prime>.25*n_sample){
    
    # alternative:
    # combn(primeFactors(204), 2, FUN = prod) %>% unique() %>% which(x = (./204<.25))
    # take second last prime factor
    batch_prime <- tail(prime_fs, 2) %>% head(1)
  }
   
  # lag training data
  data_lagged <- embed(x = as.matrix(data), dimension = (n_steps+1))
  
  # cross check conditions
  if (is.null(data_val)){
    if (internal_validation==F){
      # error
      stop('\nProvide either compatible validation data {data_val} or activate {internal_validation = T}\n')
      
      }else if (internal_validation){
        # internally determine validation data,
        # by construction sizes are consistent and evenly divide
        # train and validation.
        
        # n_sample to account val
        n_train <- n_sample - batch_prime
        n_val <- batch_prime
        
        # this ensures that mod(val_sample)==mod(n_sample)
        val_sample <- tail(data_lagged, batch_prime)  
        train_sample <- data_lagged[-((n_train + 1):n_sample),] 
      }
    
    }else if (!is.null(data_val)){ # usr provides data
      if (nrow(data_val)<=n_steps) stop('\n{data_val} is too short wrt lags {n_steps}, adjust either of them.\n')
      
      # user provides val data, need to check
      # consistent batch/sample/validation sizes
      # raise error consistently with other options
      # if batch size is not auto, 1 or else problematic, check if it mods val and train
      
      n_val <- nrow(data_val) - n_steps
      n_train <- n_sample
      
      # check consistency with val & train
      if (usr_batch){ # usr provided batch size, check compat, lag data
        
        # test on validation
        if (!(n_val%%size_batch == 0)) stop('\nProvided incompatible {size_batch} or {data_val}, must evenly divide sample sizes less lags {n_steps}.\n')
        
        # test on train
        if (!(n_sample%%size_batch == 0)) stop('\nProvided incompatible {size_batch} or {data}, must evenly divide sample sizes less lags {n_steps}.\n')
        
        # batch size ok, assign 
        batch_prime <- size_batch
        
        # lag validation
        val_sample <- embed(as.matrix(data_val), n_steps + 1)
        
      }else if (size_batch == 'auto'){
        # size_batch in both 
        primes_cand <- intersect(numbers::primeFactors(n_sample),
                                 numbers::primeFactors(n_val))
        
        if (length(primes_cand)==0 ||
            !(batch_prime %in% primes_cand) ||
            !(n_val%%batch_prime==0)) stop('\nNo common factor in {data} and {data_val}: vary {n_steps} or activate {internal_validation = T}\n')
        
        size_batch <- batch_prime <- tail(primes_cand, 1)
      }
    
    }
  
  ##### data reshaping, inherits from previous gates
  
  ## train samples
  # target: 2D array (nobs; covars)
  y_data_arr <- array(train_sample[,1],
                      dim = c(n_train, n_feat))
  
  # features: 3D array (nobs; lags; covars)
  x_data_arr <- array(train_sample[,-1],
                      dim = c(n_train, n_steps, n_feat))
  
  ## validation samples
  # target
  y_val_arr <- array(val_sample[,1],
                     dim = c(n_val, n_feat))
  # features
  x_val_arr <- array(val_sample[,-1],
                     dim = c(n_val, n_steps, n_feat))
  
  vld <- list(x_val_arr, y_val_arr)
  
  #### Model set up
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  model_compiled <- keras_model_sequential() %>%
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = F,
               stateful = T,
               kernel_regularizer = regularizer_l2(l = 0.01),
               batch_size = size_batch
               ) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  model_online <- keras_model_sequential() %>% 
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = F,
               stateful = T,
               batch_size = 1
    ) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  ## callbacks list
  clbks <- list()
  
  if (ES){
    clbks <- c(clbks,
               callback_early_stopping(monitor = 'val_loss',
                                       mode = 'auto',
                                       patience = floor(epochs*.2),
                                       min_delta = 1e-5, 
                                       restore_best_weights = keepBest)
               )
  }
  
  clbks <- c(clbks,
             callback_reduce_lr_on_plateau(
               monitor = "val_loss",
               factor = 2.1,
               patience = 10,
               verbose = 1,
               mode = 'min',
               min_delta = 1e-03,
               cooldown = 0,
               min_lr = 0)
             )
  
  
  
  tictoc::tic('\n\nModel estimation\n')
  history <- fit(object = model_compiled, 
                 y = y_data_arr, 
                 x = x_data_arr, 
                 epochs = epochs, 
                 verbose = 2,
                 shuffle = F,
                 validation_data = vld,
                 batch_size = size_batch,
                 view_metrics = view_loss,
                 callbacks = clbk)
  tictoc::toc()
  
  out <- list()
  out[['model_fitted']] <- model_compiled
  out[['history']] <- history
  out[['model_batch']] <- batch_prime
  out[['model_weights']] <- keras::get_weights(model_compiled)
  out[['model_online']] <- keras::set_weights(object = model_online, 
                                              weights = out[['model_weights']])
  if (is.xts(data)){
    out[['time_index']] <- time_index
  }
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  return(out)
}