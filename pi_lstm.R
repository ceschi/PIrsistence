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







inflation$lstm[['data']] <- future_pmap(.l = list(data = sapply(pi, list),
                                                  train = sapply(rep(1, n), list)
                                                  ),
                                        .f = data_prepper
                                        )

library(keras)
if (!keras::is_keras_available()){
  keras::install_keras()
}


tic('Full Loop: 1 layer LSTM')
sink(file = './log_lstm_full.txt', split = T, append = F)
for (i in 1:n){
  inflation$lstm[['fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation$lstm[['data']][[i]]$train$train_norm,
                                                            # either twice the BIC lags or 9 quarters to prevent
                                                            # too much sample shrinking
                                                            n_steps = 15,
                                                            n_feat = 1,
                                                            # baseline for one single layer
                                                            nodes = 100,
                                                            # online model with one batch, workaround needed
                                                            size_batch = 'auto',
                                                            # either the max epochs or patience
                                                            epochs = fit_epochs,
                                                            ES = F,
                                                            keepBest = T)
  
  # save the fitted model (with max batch size optionally)
  keras::save_model_hdf5(object = inflation$lstm[['fullsample_1l']][[i]]$model_fitted,
                         filepath = file.path(models_dir,
                                              paste0(inflation[['names']][[i]],
                                                     '_1l_fullsample.h5'))
  )
  
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
  inflation$lstm[['fullsample_2l']][[i]] <- 
    k_fullsample_2l(data = inflation$lstm[['data']][[i]]$train$train_norm, 
                    n_steps = 15, 
                    n_feat = 1, 
                    nodes = 100, 
                    size_batch = 'auto', 
                    epochs = fit_epochs, 
                    ES = F, 
                    keepBest = T)
  # save model somewhere on disk
  save_model_hdf5(object = inflation$lstm[['fullsample_2l']][[i]]$model_fitted, 
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
  
  # double plots
  inflation$lstm[['online_pred_1l']][[i]] <- online_pred(model_fitted = inflation$lstm[['fullsample_1l']][[i]], 
                                                         model_type = 'model_online',
                                                         data_train = inflation$lstm[['data']][[i]],
                                                         horizon = fore_horiz)
  
  inflation$lstm[['online_pred_2l']][[i]] <- online_pred(model_fitted = inflation$lstm[['fullsample_2l']][[i]], 
                                                         model_type = 'model_online',
                                                         data_train = inflation$lstm[['data']][[i]],
                                                         horizon = fore_horiz)
  
  # prepare canvases for plots
  inflation$lstm$plots[['full_1l']][[i]] <- ggplot(data = inflation$lstm[['online_pred_1l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ' 1L: online forecasts')) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))
  
  inflation$lstm$plots[['full_2l']][[i]] <- ggplot(data = inflation$lstm[['online_pred_2l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ' 2L: online forecasts')) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))
  
  
  # write out plots
  ggsave(filename = file.path(graphs_dir, 
                              paste0(inflation$names[[i]], '_1l_forecast.pdf')),
         plot = inflation$lstm$plots[['full_1l']][[i]],
         device = 'pdf',
         width = 8,
         height = 9*8/16,
         units = 'in')
  
  ggsave(filename = file.path(graphs_dir, 
                              paste0(inflation$names[[i]], '_1l_forecast.pdf')),
         plot = inflation$lstm$plots[['full_1l']][[i]],
         device = 'pdf',
         width = 8,
         height = 9*8/16,
         units = 'in')
  
  # display plots
  plot(inflation$lstm$plots[['full_1l']][[i]])
  plot(inflation$lstm$plots[['full_2l']][[i]])
  
}


##### Split data in chuncks for backtesting ####################################

# this section assumes the lags from the top analysis. It shall be tested whether
# it is better to have homogeneous lags across series, say 15 lags, so to train
# the LSTMs on even grounds.

# 10y splits w/o overlap
inflation$lstm[['chunk_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                       initial = sapply(inflation[['aropti']], FUN = function(x) x + 40),
                                                       assess = fm_apply(0, n),
                                                       skip = sapply(inflation[['aropti']], FUN = function(x) x + 40 - 1),
                                                       cumulative = fm_apply(F, n)),
                                             .f = rsample::rolling_origin)

# 10y rolling windows, moves on by one quarter
inflation$lstm[['wind_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                      initial = fm_apply(4*10, n),
                                                      assess = fm_apply(0, n),
                                                      cumulative = fm_apply(F, n),
                                                      skip = fm_apply(0, n),
                                                      lag = inflation[["aropti"]]), 
                                            .f = rsample::rolling_origin)

# incremental splits: they grow over time incorporating more obs
inflation$lstm[['increm_splits']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                           initial = fm_apply(4*10, n),
                                                           assess = fm_apply(0, n),
                                                           cumulative = fm_apply(T, n),
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

tic('Chunks nested loop')
for (i in 1:n){
  # preallocate for results
  chunks[[i]] <- list()
  
  # process data chunks all at once
  prepped_chunks <- #inflation$lstm[['increm_splits']][[i]]$splits %>% 
    inflation$lstm[['chunk_10y']][[i]]$splits %>%
    lapply(FUN = rsample::analysis) %>% 
    lapply(FUN = data_prepper)
  
  # store number of chunks
  len_chunks <- length(prepped_chunks)
  
  for (s in 1:len_chunks){
    
    # pull out data
    prepped_data <- prepped_chunks[[s]]
    # fit model
    lstm_list <- k_fullsample_1l(data = prepped_data$train$train_norm, 
                                 n_steps = 10,                                  #inflation[['aropti']][[n]], 
                                 nodes = 50, 
                                 epochs = fore_epochs, 
                                 ES = F,                                        # F: because there's so little data 
                                 keepBest = F,
                                 size_batch = 'auto')                           # 'auto' is also an alternative but needs testing
    # make predictions: horizon small to avoid overestimates
    # see paper and make point clear for flatlining preds
    predictions <- online_pred(model_fitted = lstm_list, 
                               model_type = 'model_online', 
                               data_train = prepped_data, 
                               horizon = fore_horiz)
    # store id for this chunk in this series
    # to add label 
    id <- paste0(names(pi)[i], '_chunk_', s)
    predictions <- add_column(predictions, data_chunk = id)
    
    # store predictions
    chunks[[i]]$predictions[[s]] <- predictions
    
    # checks
    plot(lstm_list$history)
    
    # dump model to avoid learning spillover
    rm(lstm_list)
  }
  
  # preallocate to avoid annoying behaviour
  inflation[['lstm']][['chunks']][[i]] <- list()
  
  # convert prediction tbl to xts faster
  chunks[[i]]$predictions_xts <- lapply(X = chunks[[i]]$predictions,
                                       FUN = tbl_xts)
  
  # store all good stuff in the main list
  # simple AR(1)
  inflation[['lstm']][['chunks']][[i]][['ar1']] <- 
    future_pmap(.l = list(data = chunks[[i]]$predictions_xts,
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg)
  
  # AR(1) w rolling window
  inflation$lstm$chunks[[i]][['ar1_wind']] <- 
    future_pmap(.l = list(df = chunks[[i]]$predictions_xts,
                          window = fm_apply(20, len_chunks),
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = rolloop)
  
  # simple AR(3) - SOC
  inflation$lstm$chunks[[i]][['ar3']] <- 
    future_pmap(.l = list(data = chunks[[i]]$predictions_xts,
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg.sum)
  
  # AR(3) - rolling SOC
  inflation$lstm$chunks[[i]][['ar3_wind']] <- 
    future_pmap(.l = list(df = chunks[[i]]$predictions_xts,
                          window = fm_apply(20, len_chunks),
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = rolloop.sum)
  
  # stitch all chunks back together with forecasts
  inflation$lstm$chunks[[i]]$predictions <- bind_rows(chunks[[i]]$predictions)
  
  # hairplot 
  inflation$lstm$chunks[[i]]$plot_hair <- 
    inflation$lstm$chunks[[i]]$predictions %>% ggplot() + 
        geom_line(aes(x = date, y = value, colour = label, group = data_chunk))+
        theme_minimal() + xlab(label = element_blank()) + 
        ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ': forecasts on ',len_chunks, ' chunks' )) + 
        theme(legend.position = 'bottom', 
              legend.title = element_blank())+
        guides(colour = guide_legend(nrow = 1))
  
  # filename title
  tt <- paste0(inflation$names[[i]], '_', len_chunks, '_chunks_forecasts.pdf')
  
  # save plots
  ggsave(filename = file.path(graphs_dir, tt),
         plot = inflation$lstm$increm_chunks[[i]]$plot_hair, 
         device = 'pdf', 
         width = 8, 
         height = 9*8/16, 
         units = 'in')
  
  # display
  plot(inflation$lstm$chunks[[i]]$plot_hair)
  
  # some housekeeping
  rm(tt, prepped_data, predictions, id)
  invisible(gc())
  
  cat('Done with model on ', inflation$names[[i]])
}
toc()

rm(chunks, len_chunks)


##### Now increasing subsamples ################################################

incre_win <- list()

tic('Chunks nested loop')
for (i in 1:n){
  # preallocate for results
  incre_win[[i]] <- list()
  
  # process data chunks all at once
  prepped_chunks <- #inflation$lstm[['increm_splits']][[i]]$splits %>% 
    inflation$lstm[['increm_splits']][[i]]$splits %>%
    lapply(FUN = rsample::analysis) %>% 
    lapply(FUN = data_prepper)
  
  # store number of chunks
  len_chunks <- length(prepped_chunks)
  
  for (s in 1:len_chunks){
    
    # pull out data
    prepped_data <- prepped_chunks[[s]]
    # fit model
    lstm_list <- k_fullsample_1l(data = prepped_data$train$train_norm, 
                                 n_steps = 10,                                  #inflation[['aropti']][[n]], 
                                 nodes = 50, 
                                 epochs = fore_epochs, 
                                 ES = F,                                        # F: because there's so little data 
                                 keepBest = F,
                                 size_batch = 'auto')                           # 'auto' is also an alternative but needs testing
    # make predictions: horizon small to avoid overestimates
    # see paper and make point clear for flatlining preds
    predictions <- online_pred(model_fitted = lstm_list, 
                               model_type = 'model_online', 
                               data_train = prepped_data, 
                               horizon = fore_horiz)
    # store id for this chunk in this series
    # to add label 
    id <- paste0(names(pi)[i], '_increm_win_', s)
    predictions <- add_column(predictions, data_chunk = id)
    
    # store predictions
    incre_win[[i]]$predictions[[s]] <- predictions
    
    # checks
    plot(lstm_list$history)
    
    # dump model to avoid learning spillover
    rm(lstm_list)
  }
  
  # preallocate to avoid annoying behaviour
  inflation[['lstm']][['increm_chunks']][[i]] <- list()
  
  # convert prediction tbl to xts faster
  incre_win[[i]]$predictions_xts <- lapply(X = incre_win[[i]]$predictions,
                                        FUN = tbl_xts)
  
  # store all good stuff in the main list
  # simple AR(1)
  inflation[['lstm']][['increm_chunks']][[i]][['ar1']] <- 
    future_pmap(.l = list(data = incre_win[[i]]$predictions_xts,
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg)
  
  # AR(1) w rolling window
  inflation$lstm$increm_chunks[[i]][['ar1_wind']] <- 
    future_pmap(.l = list(df = incre_win[[i]]$predictions_xts,
                          window = fm_apply(20, len_chunks),
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = rolloop)
  
  # simple AR(3) - SOC
  inflation$lstm$increm_chunks[[i]][['ar3']] <- 
    future_pmap(.l = list(data = incre_win[[i]]$predictions_xts,
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg.sum)
  
  # AR(3) - rolling SOC
  inflation$lstm$increm_chunks[[i]][['ar3_wind']] <- 
    future_pmap(.l = list(df = incre_win[[i]]$predictions_xts,
                          window = fm_apply(20, len_chunks),
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = rolloop.sum)
  
  # stitch all chunks back together with forecasts
  inflation$lstm$increm_chunks[[i]]$predictions <- bind_rows(incre_win[[i]]$predictions)
  
  # hairplot 
  inflation$lstm$increm_chunks[[i]]$plot_hair <- 
    inflation$lstm$increm_chunks[[i]]$predictions %>% ggplot() + 
    geom_line(aes(x = date, y = value, colour = label, group = data_chunk))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ': forecasts on ',len_chunks, ' chunks' )) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))
  
  # filename title
  tt <- paste0(inflation$names[[i]], '_', len_chunks, '_increm_chunks_forecasts.pdf')
  
  # save plots
  ggsave(filename = file.path(graphs_dir, tt),
         plot = inflation$lstm$increm_chunks[[i]]$plot_hair, 
         device = 'pdf', 
         width = 8, 
         height = 9*8/16, 
         units = 'in')
  
  # display
  plot(inflation$lstm$increm_chunks[[i]]$plot_hair)
  
  # some housekeeping
  rm(tt, prepped_data, predictions, id)
  invisible(gc())
  
  cat('Done with model on ', inflation$names[[i]])
}
toc()

rm(increm_chunks, len_chunks)