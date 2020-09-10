##### Freq. plots for differing window sizes ###################################

# setup window widths in years

windows <- c(5, 10, 20)
other_wind <- list()

# loop over widths to make graphs

for (w in 1:length(windows)){
  
  lbl <- paste0(windows[w], 'y_')

  other_wind[['rollark']] <- lapply(X = pi,
                                   FUN = rolloop,
                                   window = windows[w]*4,
                                   lags = 1,
                                   interc = intercep)

  other_wind[['aroptirollm']] <- future_pmap(.l = list(df = sapply(pi, list),
                                                      window = sapply(rep(windows[w]*4, n), list),
                                                      lags = inflation[['aropti']],
                                                      interc = sapply(rep(intercep, n), list)
                                                      ),
                                            .f = rolloop.sum)
  
  
  #### Plots et al ##########################################################
  
  # AR(1) rolling
  
  other_wind[['plot_rollm']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                                     names = paste0(lbl,inflation[['names']]),
                                                     path = sapply(rep(ar1_dir, n), list)),
                                           .f = plot_roller)
  
  
  # AR(k*) plots
  
  other_wind[['plot_aropti']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                                      names = paste0(lbl, inflation[['names']]),
                                                      laags = inflation[['aropti']],
                                                      path = sapply(rep(ark_dir, n), list)),
                                            .f = plot_autoregsum)
  
}

rm(other_wind, windows, lbl)