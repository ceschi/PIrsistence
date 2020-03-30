# new file for inflation analysis
# start anew with different structure
# each step is run in vectorised way
# on each series.


#' *MAIN SERIES ARE Q-O-Q ANNUALISED!*



#### 0 - Setup, data, flags ####

# horizon for now/forecast
ahead <- 1

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

# flag on the MarkovS states,
# default 2
flag___ms <- 2


# directories, functions and data

source('directories.R')
source('functs.R')
source('USdatacoll.R')

plan(multiprocess)

# subselect data
pi <- merge(
  # nowcasts and forecasts
  # db_US$cpit,
  # db_US$coret,
  # db_US$deflt,
  # db_US$deflt1,

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

# write out to disk the series
write.zoo(x=pi, 
          file=file.path(data_dir, 'PI_data.csv'), 
          sep=';',
          row.names=F, 
          index.name='date')

# chunck for further analysis in MATLAB
source('matlab_exp.R')

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

# run ADF test with max llags, minimises BIC, computes pval
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


# inflation[['aroptirollm_var']]


# Markov Switching model on the k* lags
# on the whole sample


# NOT USEFUL NOW, NOT INFORMATIVE

# inflation[['aropti_ms']] <- future_pmap(.l = list(df = sapply(pi,FUN = function(x) list(as.data.frame(x))),
#                                            lags = inflation[['aropti']],
#                                            states = flag___ms
#                                            ),
#                                  .f = ms_aropti)




# ridges plot material

inflation[['aroptiridges']] <- future_pmap(.l = list(tseries = sapply(pi, list),
                                              window = sapply(rep(wind, n), list),
                                              lags = inflation[['aropti']]),
                                    .f = persistence_ridges)







#### IV - Plots et al ##########################################################

# AR(1) rolling

inflation[['plot_rollm']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                            names = inflation[['names']],
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
                                             nam = inflation[['names']],
                                             laags = inflation[['aropti']],
                                             path = sapply(rep(graphs_dir, n), list)),
                                    
                                   .f = plot_ridges
                                   )


# # plotting msm
# 
# inflation[['plot_aropti_ms']] <- future_pmap(.l = list(ms_model = inflation[['aropti_ms']],
#                                                 nam = inflation[['names']],
#                                                 laags = inflation[['aropti']],
#                                                 path = sapply(rep(graphs_dir, n), list)
#                                                 ),
#                                                 
#                                       .f = plot_msm
#                                         )


inflation[['plot_ts']] <- ggplot(pi["1945/2020"], aes(x = index(pi["1945/2020"]))) + 
      geom_line(aes(y = rev_cpi_pch, colour = 'CPI'), alpha = .75) + 
      geom_line(aes(y = rev_cpi_fe_pch, colour = 'CPI FE'), alpha = .75) + 
      geom_line(aes(y = rev_pce_pch, colour = 'PCE'), alpha = .75) + 
      geom_line(aes(y = rev_pce_fe_pch, colour = 'PCE FE'), alpha = .75)+
      geom_line(aes(y = rev_defl_pch, colour = 'Deflt.'), alpha = .75) +
      theme_minimal() + labs(colour = ' ') +
      ggtitle('Inflation series') + xlab(' ') + ylab(' ') +
      guides(colour=guide_legend(nrow = 1, byrow = T)) + 
      theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45))

ggsave(filename = file.path(graphs_dir, 'ts_plot.pdf'),
       plot = inflation[['plot_ts']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)