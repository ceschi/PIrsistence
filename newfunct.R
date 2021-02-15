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
  compile(optimizer= optimizer_adam(lr = 10, decay = 1),
          loss = 'mse')


his <- fit(mdl_tr,
           x = yy_arr_tx,
           y = yy_arr_ty,
           shuffle = F,
           # validation_split = .5,
           validation_data = list(yy_arr_vx, yy_arr_vy),
           epochs = 30,
           batch_size = batch, callbacks = uot
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


