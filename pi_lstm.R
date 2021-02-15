##### LSTM part ################################################################

# This subscript runs ML LSTM analyses in three parts:
# - first, networks are trained on the whole sample lenghts in two flavours,
#   with one or two layers; then predictions from those trained networks are
#   produced and stored.
# 
# - secondly, LSTMs are trained and predictions made on non-overlapping subsamples
#   of the data; these predictions are then used to compute metrics of inertia changes
#   
# - third, a rolling window framework is used to track inflation persistence changes
#   over time, with the same approach as in bayesian analysis - that is train a model
#   and use it to produce several forecasts ahead, stitch data and forecasts and
#   compute some metrics 




##### Data prep for ALL models #################################################
tic('LSTM data prep')
source('pi_lstm_dataprep.R')
toc()

# safety check for keras
library(reticulate)
reticulate::use_condaenv('r-reticulate', required = T)
library(keras)
invisible(keras::is_keras_available())
if (!keras::is_keras_available()){
  keras::install_keras()
}

####' *TIME SAVING BACKSTOP*'
n <- 5 # only pch series

##### LSTM: 1L, full sample ####################################################
source('pi_lstm_1l.R')

##### LSTM: 2L, full sample ####################################################
source('pi_lstm_2l.R')

##### LSTM on 10y of data ######################################################
tic('\n10 years chunks (with lags)')
source('pi_lstm_10y.R')
toc()


##### LSTM on rolling sample ###################################################
tic('\n10y rolling samples')
source('pi_lstm_wind_10y.R')
toc()


##### increasing samples #######################################################
#' *to date it does not provide meaningful output on persistence!*
# tic('Increasing samples')
# source('pi_lstm_increm.R')
# toc()

# save all elements to disk: needs to be adapted
# saveRDS(mylist, file = "MYLIST.Rds")
# read.all <- readRDS("MYLIST.Rds")
#' *a bit problematic for keras output!*


##### Save results to disk #####################################################
saveRDS(object = inflation$lstm,
        file = file.path(rds_dir,'lstm_list.Rds'))

#### TODO ######################################################################
# - store trained model for flex reuse
# - implement HL measure w/ while condition
# - plot losses ggridges
# - review 'increm' script