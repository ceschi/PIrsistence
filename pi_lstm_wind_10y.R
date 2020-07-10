##### LSTM on rolling samples of 10y ###########################################

rolling_wind <- list()


for (i in 1:n){
  # preallocate for results
  rolling_wind[[i]] <- list()
  
  # process data chunks all at once
  prepped_chunks <- #inflation$lstm[['increm_splits']][[i]]$splits %>% 
    inflation$lstm[['wind_10y']][[i]]$splits %>%
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
                                 nodes = 50,                                    # tune he nodes to max 75
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
    id <- paste0(names(pi)[i] %>% noms, '_rolling_wind_', s)
    predictions <- add_column(predictions, data_chunk = id)
    
    # store predictions
    rolling_wind[[i]]$predictions[[s]] <- predictions
    
    # store stuff for hairplot
    if (s==1){
      rolling_wind[[i]]$selected_predictions[[s]] <- predictions
    }else{
      # to avoid double counts due to increasing/moving sample
      rolling_wind[[i]]$selected_predictions[[s]] <- tail(x = predictions,
                                                       n = fore_horiz+1)
    }
    
    # checks
    # plot(lstm_list$history)
    
    # dump model to avoid learning spillover
    rm(lstm_list)
  }
  
  # preallocate to avoid annoying behaviour
  inflation[['lstm']][['rolling_wind']][[i]] <- list()
  
  # convert prediction tbl to xts faster
  rolling_wind[[i]]$predictions_xts <- lapply(X = incre_win[[i]]$predictions,
                                           FUN = tbl_xts)
  
  # store all good stuff in the main list
  # simple AR(1)
  inflation[['lstm']][['rolling_wind']][[i]][['ar1']] <- 
    future_pmap(.l = list(data = rolling_wind[[i]]$predictions_xts,
                          lags = fm_apply(1, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg)
  
  # # AR(1) w rolling window
  # inflation$lstm$rolling_wind[[i]][['ar1_wind']] <- 
  #   future_pmap(.l = list(df = rolling_wind[[i]]$predictions_xts,
  #                         window = fm_apply(20, len_chunks),
  #                         lags = fm_apply(1, len_chunks),
  #                         interc = fm_apply(intercep, len_chunks)),
  #               .f = rolloop)
  
  # simple AR(3) - SOC
  inflation$lstm$rolling_wind[[i]][['ar3']] <- 
    future_pmap(.l = list(data = rolling_wind[[i]]$predictions_xts,
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg.sum)
  
  # # AR(3) - rolling SOC
  # inflation$lstm$rolling_wind[[i]][['ar3_wind']] <- 
  #   future_pmap(.l = list(df = rolling_wind[[i]]$predictions_xts,
  #                         window = fm_apply(20, len_chunks),
  #                         lags = fm_apply(3, len_chunks),
  #                         interc = fm_apply(intercep, len_chunks)),
  #               .f = rolloop.sum)
  
  # stitch all chunks back together with forecasts
  inflation$lstm$rolling_wind[[i]]$predictions <- bind_rows(rolling_wind[[i]]$selected_predictions)
  
  # some patchwork for the plot
  inflation$lstm$rolling_wind[[i]]$predictions <- bind_rows(inflation$lstm$rolling_wind[[i]]$predictions %>% 
                                                               filter(label == 'train') %>% 
                                                               mutate(data_chunk = '-'),
                                                             inflation$lstm$rolling_wind[[i]]$predictions %>% 
                                                               filter(label == 'forecast')
                                                            )
  
  # hairplot 
  #' *needs fixing labels*
  inflation$lstm$rolling_wind[[i]]$plot_hair <- 
    inflation$lstm$rolling_wind[[i]]$predictions %>% ggplot() + 
    geom_line(aes(x = date, y = value, colour = label, group = interaction(label, data_chunk), alpha = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ': forecasts on 10y rolling window' )) + 
    theme(legend.position = 'bottom',
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    guides(colour = guide_legend(nrow = 1))+
    scale_alpha_discrete(range = c(.03, 1))+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black'))
  
  # filename title
  tt <- paste0(inflation$names[[i]] %>% noms,
               '_',
               len_chunks, 
               '_rolling_forecasts.pdf')
  
  # save plots
  ggsave(filename = file.path(graphs_dir, tt),
         plot = inflation$rolling_wind[[i]]$plot_hair, 
         device = 'pdf', 
         width = 8, 
         height = 9*8/16, 
         units = 'in')
  
  # display
  plot(inflation$lstm$rolling_window[[i]]$plot_hair)
  
  
  # additional analysis
  inflation$lstm$rolling_window[[i]][['freq_stats']] <- 
    chunk_rolling(regs_list = inflation[['lstm']][['rolling_window']][[i]][['ar1']],
                  regs_list_sum = inflation$lstm$rolling_window[[i]][['ar3']], 
                  ar_lags_sum = 3,
                  fore_horiz = fore_horiz)
  
  # save and print plots
  inflation$lstm$rolling_window[[i]][['plots_freq_stats']] <- 
    plot_rollregs_lines(chunk_regs_obj = inflation$lstm$rolling_window[[i]][['freq_stats']],
                        graphs_dir. = graphs_dir, 
                        name = inflation$names[[i]])
  
  # some housekeeping
  rm(tt, prepped_data, predictions, id, prepped_chunks)
  invisible(gc())
  
  cat('\n\n\n\nDone with model on ', inflation$names[[i]])
  
}


rm(incre_win, len_chunks)
