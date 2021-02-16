##### Callbacks


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
  compile(optimizer= optimizer_adam(),
          loss = 'mse')


clbk <- R6::R6Class(classname = 'bump_plateau',
                    # ideally, a custom call back that 
                    # after 'trigger_epoch' periods on high
                    # loss bumps the lr to a higher level
                    # for some periods, then starts to
                    # converge back to the starting 
                    # lr in an exponential decayng way
                    inherit = KerasCallback,
                    public = list(start_lr = NULL,
                                  high_lr = NULL,
                                  i_epoch = NULL,
                                  i_high_lr = NULL,
                                  flag_high_lr = NULL,
                                  i_loss = NULL,
                                  trigger_epoch = NULL,
                                  epoch = NULL,
                                  explore = NULL,
                                  current_lr = NULL,
                                  current_loss = NULL,
                                  current_val_loss = NULL,
                                  threshold_loss = NULL,
                                  threshold_val_loss = NULL,
                                  hist_i_loss = NULL,
                                  hist_i_high_lr = NULL,
                                  hist_flag = NULL,
                                  hist_lr = NULL,
                                  outprod = NULL,
                                  initialize = function(high_lr = 2,
                                                        explore = 20,
                                                        threshold_loss = 2,
                                                        threshold_val_loss = 6,
                                                        trigger_epoch = 100){
                                    
                                    # bounce parameters
                                    self$high_lr <- high_lr
                                    self$explore <- explore
                                    # thresholds
                                    self$threshold_loss <- threshold_loss
                                    self$threshold_val_loss <- threshold_loss*3
                                    # counters
                                    self$i_epoch <- 
                                    self$i_high_lr <- 0
                                    self$i_loss <- 0
                                    self$trigger_epoch <- trigger_epoch
                                    self$epoch <- c()
                                    # flag
                                    self$flag_high_lr <- FALSE
                                    # store starting lr
                                    self$start_lr <- k_get_value(self$model$optimizer$lr)
                                    self$current_lr <- k_get_value(self$model$optimizer$lr)
                                  },
                                  
                                  on_epoch_end = function(epoch, logs = list()){
                                    
                                    # count epoch above threshold
                                    if (logs[['loss']]>self$threshold_loss){
                                      self$i_loss <- self$i_loss + 1
                                    }
                                    
                                    # if too much epochs, trigger high lr
                                    if (self$i_loss>self$trigger_epoch){
                                      
                                      # bump up lr
                                      k_set_value(self$model$optimizer$lr,
                                                  self$high_lr)
                                      
                                      # reset loss counter
                                      self$i_loss <- 0
                                      self$i_high_lr <- 0
                                      
                                      # set flag
                                      self$flag_high_lr <- TRUE
                                    }
                                    
                                    if (self$flag_high_lr){
                                      # update counter
                                      self$i_high_lr <- self$i_high_lr + 1
                                    }
                                    
                                    # after exploration, restore lr
                                    if (self$i_high_lr>self$explore){
                                      k_set_value(self$model$optimizer$lr,
                                                  self$start_lr)
                                    }
                                    
                                    
                                    # store current lr, loss, val
                                    self$current_lr <- k_get_value(self$model$optimizer$lr)
                                    self$current_loss <- c(self$current_loss,
                                                           logs[['loss']])
                                    self$current_val_loss <- c(self$current_val_loss,
                                                               logs[['val_loss']])
                                    
                                    # keep track of epochs
                                    self$i_epoch <- self$i_epoch + 1
                                    self$epoch <- c(self$epoch,
                                                    self$i_epoch)
                                    # keep track of stuff
                                    self$hist_i_loss <- c(self$hist_i_loss,
                                                          self$i_loss)
                                    self$hist_i_high_lr <- c(self$hist_i_high_lr,
                                                             self$i_high_lr)
                                    self$hist_flag <- c(self$hist_flag,
                                                        self$flag_high_lr)
                                    self$hist_lr <- c(self$hist_lr,
                                                      self$current_lr)
                                    
                                  },
                                  on_train_end = function(logs=list()){
                                    
                                    # assemble df with stats and behaviour
                                    self$outprod <- data.frame(
                                      # epoch = self$epoch,
                                      loss = self$current_loss,
                                      val_loss = self$current_val_loss,
                                      lr = self$hist_lr,
                                      loss_above = self$hist_i_loss,
                                      high_lr = self$hist_i_high_lr,
                                      flag_high_lr = self$hist_flag)
                                  }
)
)

## testing 
u <- clbk$new()


his <- fit(mdl_tr,
           x = yy_arr_tx,
           y = yy_arr_ty,
           shuffle = F,
           # validation_split = .5,
           validation_data = list(yy_arr_vx, yy_arr_vy),
           epochs = 5,
           batch_size = batch, 
           callbacks = list(u)
           )  


u$outprod
