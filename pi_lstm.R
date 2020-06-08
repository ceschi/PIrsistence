##### LSTM part ###############################################################

# todo list
# - function to prep data                                                 DONE
# - function to fit model on whole sample                                 DONE
# - predictions                                                           DONE
# - function to slice data two ways                                       DONE
#   + rolling window                                                      DONE
#   + increasing width                                                    DONE
#   + chunks with 0 offset                                                DONE
# - reuse function to fit models                                          DONE?
# - use stored models to make predictions two ways
#   + indirectly, by iterating on previous forecasts
#   + directly, by specifying an appropriate model
# - compute persistence
# - plot all of the above

# restructure list for having a big node for LSTM stuff







inflation[['lstm_data']] <- future_pmap(.l = list(data = sapply(pi, list),
                                                  train = sapply(rep(1, n), list)
                                                  ),
                                        .f = data_prepper
                                        )


if (!keras::is_keras_available()){
  keras::install_keras()
}


tic('Full Loop: 1 layer LSTM')
sink(file = './log_lstm_full.txt', split = T, append = F)
for (i in 1:n){
  inflation[['lstm_fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation[['lstm_data']][[i]]$train$train_norm,
                                                            # either twice the BIC lags or 9 quarters to prevent
                                                            # too much sample shrinking
                                                            n_steps = min(inflation[['aropti']][[i]]*2,9),
                                                            n_feat = 1,
                                                            # baseline for one single layer
                                                            nodes = 75,
                                                            # online model with one batch, workaround needed
                                                            size_batch = 1,
                                                            # either the max epochs or patience
                                                            epochs = 40,
                                                            ES = T,
                                                            keepBest = T)
  
  # save the fitted model (with max batch size optionally)
  keras::save_model_hdf5(object = inflation[['lstm_fullsample_1l']][[i]]$model_fitted,
                         filepath = file.path(models_dir,
                                              paste0(inflation[['names']][[i]],
                                                     '_1l_fullsample.h5'))
  )
  
  # ### tester
  # inflation[['lstm_fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation[['lstm_data']][[i]]$train$train_norm,
  #                                                     # either twice the BIC lags or 9 quarters to prevent
  #                                                     # too much sample shrinking
  #                                                     n_steps = 1,
  #                                                     n_feat = 1,
  #                                                     # baseline for one single layer
  #                                                     nodes = 1,
  #                                                     # online model with one batch, workaround needed
  #                                                     size_batch = 1,
  #                                                     # either the max epochs or patience
  #                                                     epochs = 2,
  #                                                     ES = F,
  #                                                     keepBest = F)
  
  # todo improvements:
  #   - let batch size depend on highest prime factor in n_sample - done but critical
  #   - critical error when using validation set and prime factor batch size:
  #       it's possible to def outside fit() the validation data w/ size of 
  #       precisely one batch, but still there is some dimensions mismatch + it's
  #       really a lot of data that the model does not see and fit on..
}
toc()
sink(NULL)

sink(file = './log_lstm_2l.txt', split = T, append = F)
tic('Full loop: 2 layers LSTM')
for (i in 1:n){
  # fit model
  inflation[['lstm_fullsample_2l']][[i]] <- 
    k_fullsample_2l(data = inflation[['lstm_data']][[i]]$train$train_norm, 
                    n_steps = min(inflation[['aropti']][[i]]*2, 12), 
                    n_feat = 1, 
                    nodes = 75, 
                    size_batch = 1, 
                    epochs = 40, 
                    ES = T, 
                    keepBest = T)
  # save model somewhere on disk
  save_model_hdf5(object = inflation[['lstm_fullsample_2l']][[i]]$model_fitted, 
                  filepath = file.path(models_dir,
                                       paste0(inflation[['names']][[i]],
                                              '_2l_fullsample.h5')
                  )
  )
}
toc()
sink(NULL)

##### If models are fitted externally, load in those files

# for (i in 1:n){
#   inflation[['lstm_fullsample']][[i]]$model_fitted <-
#     keras::load_model_hdf5(filepath = file.path(paste0(models_dir,'_4k_n75/'),
#                                                 paste0(inflation[['names']][[i]],
#                                                        ' fullsample.h5')),
#                                                 compile = T)
# }



##### Online predictions #######################################################
# predictions for 2L models
for (i in 1:n){
  inflation[['lstm_online_pred_1l']][[i]] <- online_pred(model_fitted = inflation[['lstm_fullsample_1l']][[i]], 
                                                         model_type = 'model_online',
                                                         data_train = inflation[['lstm_data']][[i]],
                                                         horizon = 40)
  
  inflation[['lstm_online_pred_2l']][[i]] <- online_pred(model_fitted = inflation[['lstm_fullsample_2l']][[i]], 
                                                         model_type = 'model_online',
                                                         data_train = inflation[['lstm_data']][[i]],
                                                         horizon = 40)
  
  inflation[['plot_lstm_full_1l']][[i]] <- ggplot(data = inflation[['lstm_online_pred_1l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ' 1L')) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))
  
  inflation[['plot_lstm_full_2l']][[i]] <- ggplot(data = inflation[['lstm_online_pred_2l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ' 2L')) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))
  
  
  plot(inflation[['plot_lstm_full_1l']][[i]])
  plot(inflation[['plot_lstm_full_2l']][[i]])
  
}


##### Split data in chuncks for backtesting ####################################

# this section assumes the lags from the top analysis. It shall be tested whether
# it is better to have homogeneous lags across series, say 15 lags, so to train
# the LSTMs on even grounds.

# 10y splits w/o overlap
inflation[['lstm_chunk_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                       initial = sapply(inflation[['aropti']], FUN = function(x) x + 40),
                                                       assess = fm_apply(0, n),
                                                       skip = sapply(inflation[['aropti']], FUN = function(x) x + 40 - 1),
                                                       cumulative = fm_apply(F, n)
),
.f = rsample::rolling_origin)

# 10y rolling windows, moves on by one quarter
inflation[['lstm_wind_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                      initial = fm_apply(4*10, n),
                                                      assess = fm_apply(0, n),
                                                      cumulative = fm_apply(F, n),
                                                      skip = fm_apply(0, n),
                                                      lag = inflation[["aropti"]]), 
                                            .f = rsample::rolling_origin)

# incremental splits: they grow over time incorporating more obs
inflation[['lstm_increm_splits']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                           initial = fm_apply(4*10, n),
                                                           assess = fm_apply(0, n),
                                                           cumulative = fm_apply(F, n),
                                                           skip = fm_apply(0, n),
                                                           lag = inflation[["aropti"]]), 
                                                 .f = rsample::rolling_origin)


##### LSTM on 10y of data ######################################################

# this section takes the data chuncks above, rescales, train an LSTM, makes
# predictions, computes persistence on train+forecast data. In addition it puts 
# together the resulting dataframes for plotting. 

# Nested loops to take care of differing subsamples. Possibly a more efficient
# way exists.

chunks <- list()

for (i in 1:n){
  # preallocate for results
  chunks[[i]] <- list()
  
  # process data chunks all at once
  prepped_chunks <- inflation[['lstm_chunk_10y']][[i]]$splits %>% 
    lapply(FUN = rsample::analysis) %>% 
    lapply(FUN = data_prepper)
  
  # store number of chunks
  len_chunks <- length(prepped_chunks)
  
  for (s in 1:len_chunks){
    
    
    # pull out data
    prepped_data <- prepped_chunks[[s]]
    # fit model
    lstm_list <- k_fullsample_1l(data = prepped_data$train$train_norm, 
                                 n_steps = inflation[['aropti']][[n]], 
                                 nodes = 5, 
                                 epochs = 10, 
                                 ES = F, 
                                 keepBest = F,
                                 size_batch = 'auto')
    # make predictions: horizon small to avoid overestimates
    # see paper and make point clear for flatlining preds
    predictions <- online_pred(model_fitted = lstm_list, 
                               model_type = 'model_online', 
                               data_train = prepped_data, 
                               horizon = 5*4)
    # store id for this chunk in this series
    # to add label 
    id <- paste0(names(pi)[i], '_chunk_', s)
    predictions <- add_column(predictions, data_chunk = id)
    
    # store predictions
    chunks[[i]]$predictions[[s]] <- predictions
    
    # compute simple AR(1)
    
    # compute AR(1) with rolling window
    
    # compute AR(k) with given k, rolling window, sum of coefficients
  }
  
  
  # simple AR(1)
  inflation$lstm$chunks[[i]][['ar1']] <- 
    future_pmap(.l = list(data = chunks[[i]]$predictions,
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg)
  
  # AR(1) w rolling window
  inflation$lstm$chunks[[i]][['ar1_wind']] <- 
    future_pmap(.l = list(df = chunks[[i]]$predictions,
                          window = fm_apply(20, len_chunks),
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = rolloop)
  
  # simple AR(3) - SOC
  inflation$lstm$chunks[[i]][['ar3']] <- 
    future_pmap(.l = list(data = chunks[[i]]$predictions,
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg.sum)
  
  # AR(3) - rolling SOC
  inflation$lstm$chunks[[i]][['ar3_wind']] <- 
    future_pmap(.l = list(df = chunks[[i]]$predictions,
                          window = fm_apply(20, len_chunks),
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = rolloop.sum)
  
  inflation$lstm$chunks[[i]]$predictions <- bind_rows(chunks[[i]]$predictions)
  inflation$lstm$chunks[[i]]$plot_hair <- 
    inflation$lstm$chunks[[i]]$predictions %>% ggplot() + 
        geom_line(aes(x = date, y = value, colour = label, group = data_chunk))
}

rm(splits_temp, id)