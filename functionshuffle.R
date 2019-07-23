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
  unitroot=list(),
  ark=list(),
  rollark=list(),
  aropti=list(),
  aroptilm=list(),
  plot_aropti=list(),
  rollm=list(),
  rollridges=list(),
  plot_rollm=list(),
  plot_ridges=list()
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

inflation[['ark']] <- lapply(X = pi,
                             FUN = auto.reg,
                             lags = 1,
                             interc = F)

inflation[['rollark']] <- lapply(X = pi,
                                 FUN = rolloop,
                                 window = wind,
                                 lags = 1,
                                 interc = F)




#### III - AR(k*) ####



