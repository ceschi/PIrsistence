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
  
  # when the sample is prime raise err,
  # surgery on data
  if (numbers::isPrime(n_sample)){
    warning('\n{data} and lags {n_steps} produce prime sample: ditching oldest observation.\nVary {n_steps} to avoid such behaviour.\n')
    
    # if n_sample is prime ditch oldest obs
    data <- data[-1,]
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
    
    # take second last prime factor
    batch_prime <- tail(prime_fs, 2) %>% head(1)
  }
  
  # lag training data
  data_lagged <- embed(x = as.matrix(data), dimension = (n_steps+1))
  
  
  #' *GATES (OF HELL) FOR VALIDATION*
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
    
  }else if (!is.null(data_val)){ # usr provides data: lagged or not? ASSUMING NOT LAGGED
    if (nrow(data_val)<=n_steps) stop('\n{data_val} is too short wtr lags {n_steps}, adjust either of them.\n')
    
    # user provides val data, need to check
    # consistent batch/sample/validation sizes
    # raise error consistently with other options
    # if batch size is not auto, 1 or else problematic, check if it mods val and train
    
    val_n_sample <- nrow(data_val) - n_steps
    
    # check consistency with val & train
    #' when usr provides *val & batch*
    if (usr_batch){ # usr provided batch size, check compat, lag data
      
      # test on validation
      data_val_gate <- val_n_sample%%size_batch == 0
      if (!data_val_gate) stop('\nProvided incompatible {size_batch} or {data_val}, must evenly divide sample sizes less lags {n_steps}.\n')
      
      # test on train
      data_train_gate <- n_sample%%size_batch == 0
      if (!data_train_gate) stop('\nProvided incompatible {size_batch} or {data}, must evenly divide sample sizes less lags {n_steps}.\n')
      
      # lag validation
      val_lagged <- embed(as.matrix(data_val), n_steps + 1)
      
      # validation target
      y_val_arr <- array(val_lagged[,1],
                          dim = c(val_n_sample, n_feat))
      x_val_arr <- array(val_lagged[,-1],
                         dim = c(val_n_sample, n_steps, n_feat))
      
      
    }else if (size_batch == 'auto'){
      
      # test: even division with val data
      val_test <- val_n_sample%%batch_prime==0
      
      # test: size_batch in both 
      primes_cand <- intersect(numbers::primeFactors(n_sample),
                               numbers::primeFactors(val_n_sample))
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


#TODO: HALFLIFE SHOCK DECAY ####################################################

#TODO: HALFLIFE FUNCTION WITH ESTIMATED MODEL

ysin <- tibble(y = 2*sin(2*base::pi*1:300))

yy <- data_prepper(ysin, train = .8, .2)


steps <- 15
feat <- 1
batch <- 5


yy_lag_t <- embed(yy$train$train_norm, steps+1)
yy_lag_v <- embed(yy$test$test_norm, steps+1)

t_sample <- nrow(yy_lag_t)
v_sample <- nrow(yy_lag_v)

yy_arr_ty <- array(yy_lag_t[,1], dim = c(t_sample, feat))
yy_arr_tx <- array(yy_lag_t[,-1], dim = c(t_sample, steps, feat))

yy_arr_vy <- array(yy_lag_v[,1], dim = c(t_sample, feat))
yy_arr_vx <- array(yy_lag_v[,-1], dim = c(t_sample, steps, feat))
  
mdl_tr <- keras_model_sequential() %>% 
  layer_lstm(units = 20,
             input_shape = c(steps, feat), 
             return_sequences = F,
             # batch_input_shape = c(5, 1),
             stateful = T,
             batch_size = batch
             ) %>% 
  layer_dense(units = 1) %>% 
  compile(optimizer= 'adam', loss = 'mse')


his <- fit(mdl_tr,
           x = yy_arr_tx,
           y = yy_arr_ty,
           shuffle = F,
           # validation_split = .5,
           validation_data = list(yy_arr_vx, yy_arr_vy),
           epochs = 1300,
           batch_size = batch
           )  


foca <- rep(0, v_sample-1) %>% c(.,1) %>% array(dim = c(v_sample,1))



mdl_tr$get_input_shape_at(0L) %>% sapply(c)


mdl_1 <- keras_model_sequential() %>% 
  layer_lstm(units = 20,
             input_shape = c(steps, feat), 
             return_sequences = F,
             # batch_input_shape = c(5, 1),
             stateful = T,
             batch_size = 1
  ) %>% 
  layer_dense(units = 1) %>% 
  compile(optimizer= 'adam', loss = 'mse')

# this one works inplace
set_weights(mdl_1, weights = get_weights(mdl_tr))

input_lagged <- embed(foca, steps)

in_arr <- array(input_lagged,
                dim = c(nrow(input_lagged), ncol(input_lagged),1))

sha <-  mdl_1$get_input_shape_at(0L) %>% sapply(c)

llast <- nrow(in_arr)

h <- 99
# not sure it immediately converges to above .5
# while (h>=.5){
#   h <- predict(mdl_1,
#                x = array(data = in_arr[llast, , ], dim = sha),
#                batch_size = 1)
# }

for (i in 1:100){
  
  input_lagged <- embed(foca, steps)
  
  in_arr <- array(input_lagged,
                  dim = c(nrow(input_lagged), ncol(input_lagged),1))
  
  llast <- nrow(in_arr)
  
  h <- predict(mdl_1,
               x = array(data = in_arr[llast, , ], dim = sha),
               batch_size = 1)
  
  foca <- c(foca, h)
  
  
}

ert <- bind_rows(tibble(yy_arr_vy) %>% add_column(lab = 'real'),
                 predict(mdl_tr, x = yy_arr_vx[nrow(yy_arr_vx), , ], batch_size = 5) %>%tibble() %>%  add_column(lab = 'pred')
                 )