##### Now increasing subsamples ################################################

incre_win <- list()

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
    id <- paste0(names(pi)[i] %>% noms, '_increm_win_', s)
    predictions <- add_column(predictions, data_chunk = id)
    
    # store predictions
    incre_win[[i]]$predictions[[s]] <- predictions
    
    # store stuff for hairplot
    if (s==1){
      incre_win[[i]]$selected_predictions[[s]] <- predictions
    }else{
      # to avoid double counts due to increasing sample
      incre_win[[i]]$selected_predictions[[s]] <- tail(x = predictions,
                                                       n = fore_horiz+1)
    }
    
    # checks
    # plot(lstm_list$history)
    
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
  inflation$lstm$increm_chunks[[i]]$predictions <- bind_rows(incre_win[[i]]$selected_predictions)
  
  # some patchwork for the plot
  inflation$lstm$increm_chunks[[i]]$predictions <- bind_rows(inflation$lstm$increm_chunks[[i]]$predictions %>% 
                                                               filter(label == 'train') %>% 
                                                               mutate(data_chunk = '-'),
                                                             inflation$lstm$increm_chunks[[i]]$predictions %>% 
                                                               filter(label == 'forecast')
                                                            )
  
  # hairplot 
  inflation$lstm$increm_chunks[[i]]$plot_hair <- 
    inflation$lstm$increm_chunks[[i]]$predictions %>% ggplot() + 
    geom_line(aes(x = date, y = value, colour = label, group = interacton(label, data_chunk)))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ': forecasts on ',len_chunks, ' chunks' )) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank())+
    guides(colour = guide_legend(nrow = 1))
  
  # filename title
  tt <- paste0(inflation$names[[i]] %>% noms,
               '_',
               len_chunks, 
               '_increm_chunks_forecasts.pdf')
  
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
  
  cat('\n\n\n\nDone with model on ', inflation$names[[i]])
  
}


rm(increm_chunks, len_chunks)