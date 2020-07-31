##### LSTM part ################################################################

# todo list
# - function to prep data                                                 DONE
# - function to fit model on whole sample                                 DONE
# - predictions                                                           DONE
# - function to slice data two ways                                       DONE
#   + rolling window                                                      DONE
#   + increasing width                                                    DONE
#   + chunks with 0 offset                                                DONE
# - reuse function to fit models                                          DONE
# - use stored models to make predictions two ways
#   + indirectly, by iterating on previous forecasts                      DONE
#   + directly, by specifying an appropriate model
# - compute persistence                                                   DONE
# - plot all of the above                                                 DONE




##### Data prep for ALL models #################################################

# full sample
inflation$lstm[['data']] <- future_pmap(.l = list(data = sapply(pi, list),
                                                  train = sapply(rep(1, n), list)
                                                  ),
                                        .f = data_prepper
                                        )

##### Split data in chuncks for backtesting ####################################

# this section assumes the lags from the top analysis. It shall be tested whether
# it is better to have homogeneous lags across series, say 15 lags, so to train
# the LSTMs on even grounds.
# We select 10y lenghts since when joined with forecasts lenght will be 14y, as 
# in the first section of the paper.

# 10y splits w/o overlap
inflation$lstm[['chunk_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                       initial = sapply(inflation[['aropti']], FUN = function(x) x + 4*10),
                                                       assess = fm_apply(0, n),
                                                       skip = sapply(inflation[['aropti']], FUN = function(x) x + 4*10 - 1),
                                                       cumulative = fm_apply(F, n)),
                                             .f = rsample::rolling_origin)

# 10y rolling windows, moves on by one quarter
inflation$lstm[['wind_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                      initial = sapply(inflation[['aropti']], FUN = function(x) x + 4*10),
                                                      assess = fm_apply(0, n),
                                                      cumulative = fm_apply(F, n),
                                                      skip = fm_apply(0, n)), 
                                            .f = rsample::rolling_origin)

# incremental splits: they grow over time incorporating more obs
inflation$lstm[['increm_splits']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                           initial = sapply(inflation[['aropti']], FUN = function(x) x + 4*10),
                                                           assess = fm_apply(0, n),
                                                           cumulative = fm_apply(T, n),
                                                           skip = fm_apply(0, n)), 
                                                 .f = rsample::rolling_origin)


####' *TIME SAVING BACKSTOP*'
n <- 5

# safety check for keras
library(keras)
if (!keras::is_keras_available()){
  keras::install_keras()
}

##### ONE layer LSTM on full sample ############################################
tic('Full Loop: 1 layer LSTM')
sink(file = './log_lstm_full.txt', split = T, append = F)
for (i in 1:n){
  inflation$lstm[['fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation$lstm[['data']][[i]]$train$train_norm,
                                                            # either twice the BIC lags or 9 quarters to prevent
                                                            # too much sample shrinking
                                                            n_steps = 15,
                                                            n_feat = 1,
                                                            # baseline for one single layer
                                                            nodes = 1000,
                                                            # online model with one batch, workaround needed
                                                            size_batch = 'auto',
                                                            # either the max epochs or patience
                                                            epochs = fit_epochs,
                                                            ES = F,
                                                            keepBest = T)
  
  # save the fitted model (with max batch size optionally)
  keras::save_model_hdf5(object = inflation$lstm[['fullsample_1l']][[i]]$model_fitted,
                         filepath = file.path(models_dir,
                                              paste0(inflation[['names']][[i]] %>% noms,
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


##### TWO layer LSTM on full sample ############################################
sink(file = './log_lstm_2l.txt', split = T, append = F)
tic('Full loop: 2 layers LSTM')
for (i in 1:n){
  # fit model
  inflation$lstm[['fullsample_2l']][[i]] <- 
    k_fullsample_2l(data = inflation$lstm[['data']][[i]]$train$train_norm, 
                    n_steps = 15, 
                    n_feat = 1, 
                    nodes = 750, 
                    size_batch = 'auto', 
                    epochs = fit_epochs, 
                    ES = F, 
                    keepBest = T)
  # save model somewhere on disk
  save_model_hdf5(object = inflation$lstm[['fullsample_2l']][[i]]$model_fitted, 
                  filepath = file.path(models_dir,
                                       paste0(inflation[['names']][[i]] %>% noms,
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
# predictions for both models
for (i in 1:n){
  
  # double predictions
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
    guides(colour = guide_legend(nrow = 1))+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black'))
  
  inflation$lstm$plots[['full_2l']][[i]] <- ggplot(data = inflation$lstm[['online_pred_2l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ' 2L: online forecasts')) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black'))
  
  
  # write out plots
  ggsave(filename = file.path(graphs_dir, 
                              paste0(inflation$names[[i]] %>% noms,
                                     '_1l_forecast.pdf')),
         plot = inflation$lstm$plots[['full_1l']][[i]],
         device = 'pdf',
         width = 8,
         height = 9*8/16,
         units = 'in')
  
  ggsave(filename = file.path(graphs_dir, 
                              paste0(inflation$names[[i]] %>% noms,
                                     '_2l_forecast.pdf')),
         plot = inflation$lstm$plots[['full_2l']][[i]],
         device = 'pdf',
         width = 8,
         height = 9*8/16,
         units = 'in')
  
  # display plots
  plot(inflation$lstm$plots[['full_1l']][[i]])
  plot(inflation$lstm$plots[['full_2l']][[i]])
  
}





##### LSTM on 10y of data ######################################################
tic('\n10 years chunks (with lags)')
source('pi_lstm_10y.R')
toc()


##### LSTM on rolling sample ###################################################
tic('\n10y rolling samples')
source('pi_lstm_wind_10y.R')
toc()


##### increasing samples #######################################################
tic('Increasing samples')
source('pi_lstm_increm.R')
toc()

# save all elements to disk: needs to be adapted
# saveRDS(mylist, file = "MYLIST.Rds")
# read.all <- readRDS("MYLIST.Rds")

##### Save results to disk #####################################################
saveRDS(object = inflation$lstm,
        file = './lstm_list')
