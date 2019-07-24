# new file for inflation analysis
# start anew with different structure
# each step is run in vectorised way
# on each series.

#### 0 - Setup, data, flags ####

# horizon for now/forecast
ahead <- 1

# exogenous lags
k <- 1

# coefficient selector
r <- 1

# rolling window width
wind <- 14*4


# flag on optimal lags:
# 1 for computing k*
# 0 for skipping
flag___optilag <- 1

# upperbound on optimal lags
llags <- 8
if (flag___optilag == 1) llags <- 18

# flag on plotting
# 1 for plotting results
# 0 for skipping
flag___plot <- 0


# directories, functions and data

source('directories.R')
source('functs.R')
source('USdatacoll.R')


# subselect data
pi <- merge(# db_US$cpit,
  # db_US$coret,
  # db_US$deflt,
  # db_US$deflt1,
  db_US$rev_cpi,
  db_US$rev_cpi_fe,
  db_US$rev_defl,
  db_US$rev_pce,
  db_US$rev_pce_fe
)

n=length(names(pi))

# this preallocated list will
# collect all results
inflation <- list(
  names=list(# 'CPI nowcast',
    # 'PCE nowcast',
    # 'GDP deflator nowcast',
    # 'GDP deflator forecast',
    'Revised CPI',
    'Revised CPI, no FE',
    'Revised GDP deflator',
    'Revised PCE',
    'Revised PCE, no FE'),
  
  # 1
  unitroot = list(),
  aropti = list(),
  
  # 2
  ark = list(),
  rollark = list(),
  
  # 3
  aroptilm = list(),
  aroptirollm = list(),
  
  
  # 4+
  plot_aropti = list(),
  rollridges = list(),
  plot_rollm = list(),
  plot_ridges = list()
)



#### I - ADF test and optimal lags #############################################

# run ADF test with max llags, minimises BIC, computes pval
inflation[['unitroot']] <- lapply(X = na.omit(pi),
                                  FUN = urca::ur.df,
                                  lags = llags,
                                  selectlags = 'BIC')


# extract and store optimal lags
inflation[['aropti']] <- lapply(X = inflation[['unitroot']],
                                FUN = function(x) return(x@optilags))



#### II - AR(1) ################################################################

# one model on the whole sample
inflation[['ark']] <- lapply(X = pi,
                             FUN = auto.reg,
                             lags = 1,
                             interc = T)

# rolling window
inflation[['rollark']] <- lapply(X = pi,
                                 FUN = rolloop,
                                 window = wind,
                                 lags = 1,
                                 interc = T)



#### III - AR(k*) ##############################################################

# one fit with optilags
# on the whole sample
inflation[['aroptilm']] <- pmap(.l = list(data = sapply(pi, list),
                                          lags = inflation[['aropti']],
                                          interc = sapply(rep(T, n), list)
                                          ),
                                .f = auto.reg.sum)


# rolling window estimate
# width is preselected,
# optilag is computed in step I
inflation[["aroptirollm "]] <- pmap(.l = list(df = sapply(pi, list),
                                              window = sapply(rep(wind, n), list),
                                              lags = inflation[['aropti']],
                                              interc = sapply(rep(T, n), list)
                                              ),
                                    .f = rolloop.sum)








#### IV - Plots et al ##########################################################














