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
source('functs.R')
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

# chunck for further analysis in MATLAB
# source('matlab_exp.R')

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




##### LSTM part ###############################################################

# todo list
  # - function to prep data                           DONE
  # - function to fit model on whole sample           DONE
  # - predictions
  # - function to slice data two ways                 
  #   + rolling window
  #   + increasing width
  # - reuse function to fit models
  # - use stored models to make predictions two ways
  #   + indirectly, by iterating on previous forecasts
  #   + directly, by specifying an appropriate model
  # - compute persistence
  # - plot all of the above



inflation[['lstm_data']] <- future_pmap(.l = list(data = sapply(pi, list),
                                                  train = sapply(rep(1, n), list)
                                                  ),
                                        .f = data_prepper
                                        )


inflation[['lstm_fullsample']] <- list()

if (!keras::is_keras_available()){
  keras::install_keras()
}

tic('Full Loop')
sink(file = './log_lstm_full.txt', split = T, append = F)
for (i in 1:n){
  inflation[['lstm_fullsample']][[i]] <- k_fullsample(data = inflation[['lstm_data']][[i]]$train$train_norm,
                                                      # either twice the BIC lags or 9 quarters to prevent 
                                                      # too much sample shrinking
                                                      n_steps = min(inflation[['aropti']][[i]]*2,9), 
                                                      n_feat = 1, 
                                                      # baseline for one single layer
                                                      nodes = 75, 
                                                      # online model with one batch, workaround needed
                                                      size_batch = 1, 
                                                      # either the max epochs or patience
                                                      epochs = 2000, 
                                                      ES = T)
  
  keras::save_model_hdf5(object = inflation[['lstm_fullsample']][[i]]$model_fitted,
                         filepath = file.path(models_dir,
                                              paste0(inflation[['names']][[i]], 
                                                     ' fullsample.h5'))
                        )
  
  # todo improvements:
  #   - this is online fit, batch = 1
  #   - let batch size depend on highest prime factor in n_sample
  #   - fit model and then copy weights into new model for online forecasts
}
toc()
sink(NULL)


##### If models are fitted externally, load in those files

for (i in 1:n){
  inflation[['lstm_fullsample']][[i]]$model_fitted <- 
    keras::load_model_hdf5(filepath = file.path(models_dir,
                                                paste0(inflation[['names']][[i]], 
                                                       ' fullsample.h5')),
                                                compile = T)
}


online_pred <- function(model_fitted, data_train, horizon = 4*10){
  
  # This function produces iterative, indirect predictions with 
  # a previously trained model. It copies weights and model structure
  # and resets batch to 1 so to make online predictions easily
  # and consistently. 'horizon' gives the nomber of indirect
  # predictions to produce. If data are TS then also dates are generated.
  
  require(keras)
  require(dplyr)
  
  # data_train is a list from data_prepper function!
  if (!is.list(data_train)) error('Provide list from "data_prepper" function')
  
  # preallocate array with results
  pred <- array(NA, dim = c(horizon, 1))
  
  if (is.xts(data_train$train[['train_norm']])){
    time_preds <- seq(from = end(input),
                      length.out = (horizon+1),
                      by = periodicity(input)$label)
    time_preds <- tail(time_preds, n = horizon)
    pred <- xts(pred, order.by = time_preds)
  }
  
  # retrieve input shape
  in_shape <- get_input_shape_at(object = model_fitted[['model_online']],
                                 node_index = 0)
  in_shape <- sapply(in_shape, FUN = c)
  
  # retrieve input data and fitted model
  input <- data_train$train[['train_norm']]
  model_online <- model_fitted[['model_online']]
  

  
  for (h in 1:horizon){
    input_lagged <- embed(input, in_shape[2])
    input_arr <- array(data = input_lagged,
                       dim = c(nrow(input_lagged), ncol(input_lagged),1))
    last_row <- nrow(input_arr)
    
    pred[h, ] <- predict(model_online,
                         x = array(data = input_arr[last_row, , ],
                                   dim = in_shape),
                         batch_size = 1)
    input <- rbind(input, pred[h,])
    }
  
  if (is.xts(data_train$train[['train_norm']])){
    
    time_preds <- seq(from = end(input),
                      length.out = (horizon+1),
                      by = periodicity(input)$label)
    time_preds <- tail(time_preds, n = horizon)
    
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train', 
                                   date = time(data_train$train[['train_norm']])) %>% 
                        rename(value = V1) %>% 
                        mutate(value = as.numeric(value)) %>% 
                        xts(order.by = .$date), 
                      pred %>% 
                        as_tibble %>% 
                        add_column(label = 'forecast', 
                                   date = time_preds) %>% 
                        rename(value = V1) %>% 
                        mutate(value = as.numeric(value)) %>% 
                        xts(order.by = .$date))
  } else {
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train'), 
                      pred %>% 
                        as_tibble %>% 
                        add_column(label = 'forecast'))
  }
  
  forecast$date <- NULL
  forecast$value <- as.numeric(forecast$value)
  forecast <- forecast*data_train$train[['sd']] + data_train$train[['mean']]
  
  
  
  
  return(forecast)
  
}