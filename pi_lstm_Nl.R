# Little nodes, many layers

# set series
n <- 5

# how many layers?


#' *NEEDS ADAPTATION*

##### TWO layer LSTM on full sample ############################################
tic('Full loop: n layers LSTM')
for (i in 1:n){
  # fit model
  inflation$lstm[['fullsample_2l']][[i]] <- 
    k_fullsample_nl(data = inflation$lstm[['data']][[i]]$train$train_norm, 
                    n_steps = 15, 
                    n_feat = 1, 
                    # nodes = 7,
                    nodes = 5,
                    size_batch = 'auto', 
                    epochs = fit_epochs*2, 
                    ES = F, 
                    keepBest = T)
  # save model somewhere on disk
  save_model_hdf5(object = inflation$lstm[['fullsample_2l']][[i]]$model_fitted, 
                  filepath = file.path(models_dir,
                                       paste0(inflation[['names']][[i]] %>% noms,
                                              '_2l_fullsample.h5')
                  )
  )
  
  
  cat('\nJust done with ', inflation$names[[i]] %>% noms_tt(),
      ' model on iteration ', i,' of two-layer LSTM.\n\n')
  
  inflation$lstm[['online_pred_2l']][[i]] <- online_pred(model_fitted = inflation$lstm[['fullsample_2l']][[i]], 
                                                         model_type = 'model_online',
                                                         data_train = inflation$lstm[['data']][[i]],
                                                         horizon = fore_horiz)
  
  # compute RW-ARk on extended series
  inflation$lstm$ark_2l[[i]] <- inflation$lstm[['online_pred_2l']][[i]] %>% 
    select(-label) %>% 
    tbl_xts() %>% 
    rolloop.sum(window = wind,
                lags = inflation[['aropti']][[i]], 
                interc = intercep) %>% 
    na.omit(.)
  
  inflation$lstm$plots[['full_2l']][[i]] <- ggplot(data = inflation$lstm[['online_pred_2l']][[i]])+
    geom_line(aes(x = date, y = value, colour = label))+
    theme_minimal() + xlab(label = element_blank()) + 
    ylab(element_blank()) + ggtitle(paste0(inflation$names[[i]] %>% noms_tt(), ' 2L: online forecasts')) + 
    guides(colour = guide_legend(nrow = 1))+
    scale_colour_manual(labels = c('Forecast', 'Data'), values = c('red', 'black')) +
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          legend.position = 'bottom', 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  inflation$lstm$plots[['full_2l_rwark_coef']][[i]] <- 
    plot_autoregsum(df = inflation$lstm$ark_2l[[i]],
                    names = inflation$names[[i]], 
                    path = l2_dir, 
                    laags = inflation$aropti[[i]])
  
  inflation$lstm$plots[['full_2l_rwark_trend']][[i]] <- 
    plot_autoregsum(df = inflation$lstm$ark_2l[[i]],
                    names = inflation$names[[i]], 
                    path = l2_dir, 
                    laags = inflation$aropti[[i]], 
                    .slot = 2)
  
  # save plot
  ggsave(filename = file.path(l2_dir, 
                              paste0(inflation$names[[i]] %>% noms,
                                     '_2l_forecast.pdf')),
         plot = inflation$lstm$plots[['full_2l']][[i]],
         device = 'pdf',
         width = 8,
         height = 9*8/16,
         units = 'in')
  
  
  # cleanup stuff
  invisible(gc())
}
toc()

saveRDS(object = inflation$lstm$fullsample_2l,
        file = file.path(rds_dir,'lstm_fullsample_2l_list.rds'))