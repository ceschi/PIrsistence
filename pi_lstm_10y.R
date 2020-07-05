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
                                 n_steps = 10,                                  #inflation[['aropti']][[n]],  # here, 10q lags but this could vary
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
    id <- paste0(names(pi)[i] %>% noms, '_chunk_', s)
    predictions <- add_column(predictions, data_chunk = id)
    
    # store predictions
    chunks[[i]]$predictions[[s]] <- predictions
    
    # checks
    # plot(lstm_list$history)
    
    # dump model to avoid learning spillover
    rm(lstm_list)
  }
  
  # preallocate to avoid annoying behaviour
  inflation[['lstm']][['chunks']][[i]] <- list()
  
  
  # convert prediction tbl to xts faster
  chunks[[i]]$predictions_xts <- lapply(X = chunks[[i]]$predictions,
                                        FUN = tbl_xts)
  
  # store all the predictions and data, raw
  inflation[['lstm']][['chunks']][[i]][['raw']] <- chunks[[i]]$predictions
  inflation[['lstm']][['chunks']][[i]][['raw_xts']] <- chunks[[i]]$predictions_xts
  
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
  
  inflation$lstm$chunks[[i]]$predictions_xts <- inflation$lstm$chunks[[i]]$predictions %>% 
    group_by(data_chunk) %>% nest() %>% map(tbl_xts)
  
  d_vline <- inflation$lstm$chunks[[i]]$predictions %>% 
                filter(label == 'train') %>% 
                group_by(data_chunk) %>% 
                mutate(ll = last(date)) %>% 
                ungroup() %>% 
                distinct(ll)
  
  # hairplot 
  inflation$lstm$chunks[[i]]$plot_hair <- 
    inflation$lstm$chunks[[i]]$predictions %>% ggplot() + 
    geom_line(aes(x = date, y = value, colour = label, group = interaction(label, data_chunk)))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]], ': forecasts on ',len_chunks, ' chunks' )) + 
    theme(legend.position = 'bottom', 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    guides(colour = guide_legend(nrow = 1))+
    geom_vline(xintercept = d_vline$ll, linetype = 'dashed', alpha = .5)+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black'))
  
  # filename title
  tt <- paste0(inflation$names[[i]] %>% noms,
               '_',
               len_chunks,
               '_chunks_forecasts.pdf')
  
  # save plots
  ggsave(filename = file.path(graphs_dir, tt),
         plot = inflation$lstm$chunks[[i]]$plot_hair, 
         device = 'pdf', 
         width = 8, 
         height = 9*8/16, 
         units = 'in')
  
  # display
  plot(inflation$lstm$chunks[[i]]$plot_hair)
  
  # some housekeeping
  rm(tt, prepped_data, predictions, id, d_vline)
  invisible(gc())
  
  cat('\n\n\n\nDone with model on ', inflation$names[[i]])
}


rm(chunks, len_chunks)
# plus change colours to red/black for all series in the ML part
# + 
# scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))


