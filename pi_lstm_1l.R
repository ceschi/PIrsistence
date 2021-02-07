##### LSTM - 1L ################################################################

n <- 5 # only pch series

l1_histories <- list()
l1_net_weights <- list()

##### ONE layer LSTM on full sample ############################################
tic('Full Loop: 1 layer LSTM')
for (i in 1:n){
  inflation$lstm[['fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation$lstm[['data']][[i]]$train$train_norm,
                                                            # either twice the BIC lags or 9 quarters to prevent
                                                            # too much sample shrinking
                                                            n_steps = 15,
                                                            n_feat = 1,
                                                            # baseline for one single layer
                                                            # nodes = 2,
                                                            nodes = 20,
                                                            # online model with one batch, workaround needed
                                                            size_batch = 'auto',
                                                            # either the max epochs or patience
                                                            epochs = fit_epochs,
                                                            ES = F,
                                                            keepBest = F)
  
  # save the fitted model (with max batch size optionally)
  keras::save_model_hdf5(object = inflation$lstm[['fullsample_1l']][[i]]$model_fitted,
                         filepath = file.path(models_dir,
                                              paste0(inflation[['names']][[i]] %>% noms,
                                                     '_1l_fullsample.h5'))
  )
  
  l1_histories[[i]] <- inflation$lstm$fullsample_1l[[i]]$history$metrics
  l1_net_weights[[i]] <- inflation$lstm$fullsample_1l[[i]]$model_weights
  
  plot(inflation$lstm$fullsample_1l[[i]]$history)
  cat('\nJust done with ', inflation$names[[i]] %>% noms_tt(),
      ' model on iteration ', i,' of one-layer LSTM.\n\n')
  
  # compute predictions
  inflation$lstm[['online_pred_1l']][[i]] <- online_pred(model_fitted = inflation$lstm[['fullsample_1l']][[i]], 
                                                         model_type = 'model_online',
                                                         data_train = inflation$lstm[['data']][[i]],
                                                         horizon = fore_horiz)
  
  # compute RW-ARk on extended series
  inflation$lstm$ark_1l[[i]] <- inflation$lstm[['online_pred_1l']][[i]] %>% 
    select(-label) %>% 
    tbl_xts() %>% 
    rolloop.sum(window = wind,
                lags = inflation[['aropti']][[i]], 
                interc = intercep) %>% 
    na.omit(.)
  
  # prepare canvases for plots
  inflation$lstm$plots[['full_1l']][[i]] <- ggplot(data = inflation$lstm[['online_pred_1l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]] %>% noms_tt(), ' 1L: online forecasts')) + 
    guides(colour = guide_legend(nrow = 1))+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black')) +
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom', 
          legend.title = element_blank())
  
  inflation$lstm$plots[['full_1l_rwark_coef']][[i]] <- 
            plot_autoregsum(df = inflation$lstm$ark_1l[[i]],
                            names = inflation$names[[i]], 
                            path = l1_dir, 
                            laags = inflation$aropti[[i]])
  
  inflation$lstm$plots[['full_1l_rwark_trend']][[i]] <- 
            plot_autoregsum(df = inflation$lstm$ark_1l[[i]],
                            names = inflation$names[[i]], 
                            path = l1_dir, 
                            laags = inflation$aropti[[i]], 
                            .slot = 2)
  
  # write out plots
  ggsave(filename = file.path(l1_dir, 
                              paste0(inflation$names[[i]] %>% noms,
                                     '_1l_forecast.pdf')),
         plot = inflation$lstm$plots[['full_1l']][[i]],
         device = 'pdf',
         width = 8,
         height = 9*8/16,
         units = 'in')
  
  # cleanup stuff
  invisible(gc())
}
toc()

saveRDS(object = inflation$lstm$fullsample_1l,
        file = file.path(rds_dir,'lstm_fullsample_1l_list.rds'))