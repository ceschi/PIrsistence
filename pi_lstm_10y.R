# this section takes the data chuncks, rescales, train an LSTM, makes
# predictions, computes persistence on train+forecast data. In addition it puts 
# together the resulting dataframes for plotting. 

# Nested loops to take care of differing subsamples. Possibly a more efficient
# way exists.

chunks <- list()
chunks_histories <- list()
chunks_net_weights <- list()

for (i in 1:n){
  # preallocate for results
  chunks[[i]] <- list()
  chunks_histories[[i]] <- list()
  chunks_net_weights[[i]] <- list()
  
  # process data chunks all at once
  prepped_chunks <-  
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
                                 n_steps = 15,
                                 nodes = 500,
                                 epochs = fore_epochs,
                                 ES = F,
                                 keepBest = T,
                                 size_batch = 'auto')
    # make predictions: horizon small to avoid overestimates
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
    
    # storing
    chunks_histories[[i]][[s]] <- lstm_list$history$metrics
    chunks_net_weights[[i]][[s]] <- lstm_list$model_weights
    # plot(lstm_list$history)
    
    # dump model to avoid learning spillover
    rm(lstm_list)
    
    cat('\n\nJust ended iteration ', s, ' of ', len_chunks, ' in the inner loop. 
         \nOuter loop is at iteration ', i, ' of the total ', n, ' for the 10y samples.\n')
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
  
  # simple AR(3) - SOC
  inflation$lstm$chunks[[i]][['ar3']] <- 
    future_pmap(.l = list(data = chunks[[i]]$predictions_xts,
                          lags = fm_apply(3, len_chunks),
                          interc = fm_apply(intercep, len_chunks)),
                .f = auto.reg.sum)
  
  # stitch all chunks back together with forecasts
  inflation$lstm$chunks[[i]]$predictions <- bind_rows(chunks[[i]]$predictions)
  
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
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]] %>% noms_tt(), ': forecasts on ',len_chunks, ' chunks' )) + 
    guides(colour = guide_legend(nrow = 1))+
    geom_vline(xintercept = d_vline$ll, linetype = 'dashed', alpha = .5)+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black')) +
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5), 
          legend.position = 'bottom', 
          legend.title = element_blank())
  
  # filename title
  tt <- paste0(inflation$names[[i]] %>% noms,
               '_',
               len_chunks,
               '_chunks_forecasts.pdf')
  
  # save plots
  ggsave(filename = file.path(chunks_dir, tt),
         plot = inflation$lstm$chunks[[i]]$plot_hair, 
         device = 'pdf', 
         width = 8, 
         height = 9*8/16, 
         units = 'in')
  
  write.csv(x = inflation$lstm$chunks[[i]]$predictions,
            file = file.path(data_dir, paste0(tt,'.csv')),
            row.names = F)
  
  # display
  plot(inflation$lstm$chunks[[i]]$plot_hair)
  
  # additional analysis
  inflation$lstm$chunks[[i]][['freq_stats']] <- chunk_regs(regs_list = inflation[['lstm']][['chunks']][[i]][['ar1']],
                                                           regs_list_sum = inflation$lstm$chunks[[i]][['ar3']], 
                                                           ar_lags_sum = 3,
                                                           fore_horiz = fore_horiz)
  # save and print plots
  inflation$lstm$chunks[[i]][['plots_freq_stats']] <- 
    plot_chunkregs_bar(chunk_regs_obj = inflation$lstm$chunks[[i]][['freq_stats']],
                        graphs_dir. = chunks_dir, 
                        name = inflation$names[[i]])
  
  # save and print tables in TeX
  chunk_stargazer(ar1 = inflation[['lstm']][['chunks']][[i]][['ar1']], 
                  chunk_out = inflation$lstm$chunks[[i]][['freq_stats']], 
                  name = inflation$names[[i]], 
                  pathout = tables)
  
  # some housekeeping
  rm(tt, prepped_data, predictions, id, d_vline)
  invisible(gc())
  
  cat('\n\n\n\nDone with model on ', inflation$names[[i]])
}

rm(chunks, len_chunks, prepped_chunks)



##### Save results to disk #####################################################
saveRDS(object = inflation$lstm$chunks,
        file = file.path(rds_dir,'/lstm_chunks_10y_list.rds'))
