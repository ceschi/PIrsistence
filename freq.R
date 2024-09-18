##### I - Frequentist analyses on inflation series #############################

# set number of series
n <- length(names(pi))

# runs ADF test with max llags, minimises BIC, computes pval
inflation[['unitroot']] <- lapply(X = na.omit(pi),
                                  FUN = urcabis::ur.df.ol,
                                  lags = llags,
                                  type = 'drift',
                                  selectlags = 'BIC')


# extract and store optimal lags
inflation[['aropti']] <- lapply(X = inflation[['unitroot']],
                                FUN = function(x) return(x@optilags))



#### II - AR(1) ################################################################

# one model on the whole sample
inflation[['ark']] <- lapply(X = pi,
                             FUN = auto.reg,
                             lags = 1,
                             interc = intercep)

# rolling window
inflation[['rollark']] <- lapply(X = pi,
                                 FUN = rolloop,
                                 window = wind,
                                 lags = 1,
                                 interc = intercep)



#### III - AR(k*) ##############################################################

# one fit with optilags
# on the whole sample
inflation[['aroptilm']] <- future_pmap(.l = list(data = sapply(pi, list),
                                                 lags = inflation[['aropti']],
                                                 interc = sapply(rep(intercep, n), list)
                                                 ),
                                       .f = auto.reg.sum)


# rolling window estimate
# width is preselected,
# optilag is computed in step I
inflation[['aroptirollm']] <- future_pmap(.l = list(df = sapply(pi, list),
                                                    window = sapply(rep(wind, n), list),
                                                    lags = inflation[['aropti']],
                                                    interc = sapply(rep(intercep, n), list)
                                                    ),
                                           .f = rolloop.sum)


# ridges plot material
inflation[['aroptiridges']] <- future_pmap(.l = list(tseries = sapply(pi, list),
                                                     window = sapply(rep(wind, n), list),
                                                     lags = inflation[['aropti']]),
                                           .f = persistence_ridges)