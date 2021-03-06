#' *Inflation persistence analyses*
#' 
#' This scripts runs analyses on inflation persistence
#' in the US economy from the eponymous paper.
#' After setting up some flags that prevents parts to run and
#' detail some operations, it first fetch all necessary data 
#' to run the analyses. Importantly, some observations are dropped 
#' from the sample, eg those prior to 1966 for core CPI, since they
#' are produced artificially by the data provider (BEA); optionally,
#' also COVID-19 observations can be ditched, since they display
#' largest, unprecedented swings in inflation.
#' 
#' After collecting the data and loading up custom functions, 
#' the script compiles an 'inflation' list object that will collect 
#' all the results and plots. Then, the script exports files for Matlab
#' bayesian analysis.
#' 
#' First, frequentist analyses are run, importantly finding the optimal lag
#' number for each series that will used subsequently.
#' 
#' Then, the Machine Learning operation is set up in a separate
#' tree of scripts. This part does really take much time to complete.
#' 
#' Lastly, some more plots are produced in the last call to external script.

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
source('data_pi.R')

# reproducibility with Reis&Pivetta
# pi <- pi["/2002-12-31"]

# drop COVID19 obs?
# pi <- pi["/2020-02-29"]




# scripts for further analysis in MATLAB
# source('matlab_exp.R')
# source('matlab_plot.R')

n <- length(names(pi))

# this preallocated list will
# collect all results
inflation <- list(
  names=list(

    # continously compounded annual rate of change
    'CPI pch',
    'CPI, no FE pch',
    'PCE pch',
    'PCE, no FE pch',
    'GDP deflator pch',
    
    # percentage change from a year ago
    'CPI yoy',
    'CPI, no FE yoy',
    'PCE yoy',
    'PCE, no FE yoy',
    'GDP deflator yoy'
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
  
  # 4+
  plot_rollm = list(),
  plot_aropti = list(),
  plot_ridges = list()
)

# parallel
plan(multisession)

##### I - ADF, optilags, AR1, ARk  #############################################
tic('Frequentist part')
source('freq.R')
toc()


##### III - LSTM ###############################################################
if (flag___lstm){
  if (flag___epochs){
    fit_epochs <- 4
    fore_epochs <- 5
    fore_horiz <- 6
  } else {
    fit_epochs <- 2500
    fore_epochs <- 750
    fore_horiz <- 40
  }
  
  tic('Machine learning fit and forecasts.\n')
  
  # make sure that keras finds proper conda env
  library(reticulate)
  reticulate::use_condaenv('r-reticulate', required = T)
  
  # workaround to limit CPUs usage @PSE to 5*4 threads
  # library(tensorflow)
  # tf$config$threading$set_intra_op_parallelism_threads(7L)
  # tf$executing_eagerly()
  
  source('pi_lstm.R')
  toc()
}




##### Plots ####################################################################
tic('Plotting')
source('plotte.R')
toc()


##### TODO #####################################################################

# + chop out matlab part in separate script
#   - possibly in Julia
# + look into trend changes via the intercept moves.