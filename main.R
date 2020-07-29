# new file for inflation analysis
# start anew with different structure
# each step is run in vectorised way
# on each series.


#' *MAIN SERIES ARE Q-O-Q ANNUALISED!*



#### 0 - Setup, data, flags ####

# exogenous lags
k <- 1

# intercept
intercep <- T


# rolling window width
wind <- 14*4


# flag on optimal lags:
# 1 for computing k*
# 0 for skipping
flag___optilag <- 1

# upperbound on optimal lags
llags <- 18
# if (flag___optilag == 1) llags <- 18

# flag on plotting
# 1 for plotting results
# 0 for skipping
flag___plot <- 0

# flags for LSTM scripts
# run at all that script?
flag___lstm <- T
# run with small n of epochs?
flag___epochs <- F


# directories, functions and data
source('functs.R')
# pick ahead to set how many quarters ahead 
# to consider for SPF forecasts:
# -1 for previous quarter estimates
# 0 for nowcast
# 1 for one quarter ahead -- default
# 2 for two quarters ahead
# 3 for three quarters ahead
# 4 for one year ahead
ahead <- 1
download.file(url = 'https://raw.githubusercontent.com/ceschi/us_macro_data/master/USdata_coll.R',
      destfile = 'temp.R',
      quiet = T)
source('temp.R')
unlink(x = 'temp.R')

plan(multiprocess)

# subselect data
pi <- merge(
  # qoq percentage change
  db_US$rev_cpi_pch,
  db_US$rev_cpi_fe_pch,
  db_US$rev_defl_pch,
  db_US$rev_pce_pch,
  db_US$rev_pce_fe_pch,
  
  # yoy percentage change
  db_US$rev_cpi_yoy,
  db_US$rev_cpi_fe_yoy,
  db_US$rev_defl_yoy,
  db_US$rev_pce_yoy,
  db_US$rev_pce_fe_yoy
)

# db_US not needed in its entirety
rm(db_US)

# reproducibility with Reis&Pivetta
# pi <- pi["/2002-12-31"]

# reshape data in long format
pi_long <- pi %>% 
            as_tibble %>%
            add_column(date = time(pi)) %>% 
            gather(1:(ncol(.)-1), key = 'index', value = 'rate', na.rm = T)

# write out to disk the series
write.zoo(x=pi, 
          file=file.path(data_dir, 'PI_data.csv'), 
          sep=';',
          row.names=F, 
          index.name='date')

write_csv(x = pi_long,
          path = file.path(data_dir, 'PI_long.csv'))

# scripts for further analysis in MATLAB
# source('matlab_exp.R')
# source('matlab_plot.R')

n=length(names(pi))

# this preallocated list will
# collect all results
inflation <- list(
  names=list(
    
    # forecasts/nowcasts
    # 'CPI nowcast',
    # 'PCE nowcast',
    # 'GDP deflator nowcast',
    # 'GDP deflator forecast',

    # continously compounded annual rate of change
    'Revised CPI pch',
    'Revised CPI, no FE pch',
    'Revised GDP deflator pch',
    'Revised PCE pch',
    'Revised PCE, no FE pch',
    
    # percentage change from a year ago
    'Revised CPI yoy',
    'Revised CPI, no FE yoy',
    'Revised GDP deflator yoy',
    'Revised PCE yoy',
    'Revised PCE, no FE yoy'
    ),
  
  # 1
  unitroot = list(),
  aropti = list(),
  
  # 2
  ark = list(),
  rollark = list(),
  
  # 3
  aroptilm = list(),
  aroptirollm = list(),
  aroptiridges = list(),
  
  # 3bis
  aropti_ms = list(),
  
  
  # 4+
  plot_rollm = list(),
  plot_aropti = list(),
  plot_ridges = list(),
  plot_aropti_ms = list()
)



#### I - ADF test and optimal lags #############################################

# runs ADF test with max llags, minimises BIC, computes pval
inflation[['unitroot']] <- lapply(X = na.omit(pi),
                                  FUN = urca::ur.df,
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







#### IV - Plots et al ##########################################################

# AR(1) rolling

inflation[['plot_rollm']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                            names = inflation[['names']] %>% noms,
                                            path = sapply(rep(graphs_dir, n), list)),
                                          .f = plot_roller
                                          )


# AR(k*) plots

inflation[['plot_aropti']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                             names = inflation[['names']],
                                             laags = inflation[['aropti']],
                                             path = sapply(rep(graphs_dir, n), list)),
                                           .f = plot_autoregsum
                                           )

# plotting ridges

inflation[["plot_ridges"]] <- future_pmap(.l = list(df = inflation[['aroptiridges']],
                                             nam = inflation[['names']] %>% noms,
                                             laags = inflation[['aropti']],
                                             path = sapply(rep(graphs_dir, n), list)),
                                    
                                           .f = plot_ridges
                                           )



##### LSTM #####################################################################
if (flag___lstm){
  if (flag___epochs){
    fit_epochs <- 4
    fore_epochs <- 2
    fore_horiz <- 1
  } else {
    fit_epochs <- 5000
    fore_epochs <- 2000
    fore_horiz <- 40
  }
  # 
  # fit_epochs <- 15
  # fore_epochs <- 5
  # fore_horiz <- 4
  # 
  tic('Machine learning fit and forecasts.\n')
  
  # workaround to limit CPUs usage @PSE
  library(tensorflow)
  tf$config$threading$set_intra_op_parallelism_threads(25L)
  tf$config$threading$set_inter_op_parallelism_threads(25L)
  tf$executing_eagerly()
  
  source('pi_lstm.R')
  toc()
}




##### Plots ####################################################################

tic('Plotting')
source('plotte.R')
toc()