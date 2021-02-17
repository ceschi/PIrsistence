##### LSTM - Data Prep #########################################################

# here apply prep to all series
n <- length(names(pi))

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