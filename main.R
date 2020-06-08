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
  # - predictions                                     DONE
  # - function to slice data two ways                 DONE
  #   + rolling window                                DONE
  #   + increasing width                              DONE?
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
  inflation[['lstm_fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation[['lstm_data']][[i]]$train$train_norm,
                                                      # either twice the BIC lags or 9 quarters to prevent
                                                      # too much sample shrinking
                                                      n_steps = min(inflation[['aropti']][[i]]*2,9),
                                                      n_feat = 1,
                                                      # baseline for one single layer
                                                      nodes = 75,
                                                      # online model with one batch, workaround needed
                                                      size_batch = 1,
                                                      # either the max epochs or patience
                                                      epochs = 40,
                                                      ES = T,
                                                      keepBest = T)
  
  # save the fitted model (with max batch size optionally)
  keras::save_model_hdf5(object = inflation[['lstm_fullsample']][[i]]$model_fitted,
                         filepath = file.path(models_dir,
                                              paste0(inflation[['names']][[i]],
                                                     '_fullsample.h5'))
                        )
  
  # ### tester
  # inflation[['lstm_fullsample_1l']][[i]] <- k_fullsample_1l(data = inflation[['lstm_data']][[i]]$train$train_norm,
  #                                                     # either twice the BIC lags or 9 quarters to prevent
  #                                                     # too much sample shrinking
  #                                                     n_steps = 1,
  #                                                     n_feat = 1,
  #                                                     # baseline for one single layer
  #                                                     nodes = 1,
  #                                                     # online model with one batch, workaround needed
  #                                                     size_batch = 1,
  #                                                     # either the max epochs or patience
  #                                                     epochs = 2,
  #                                                     ES = F,
  #                                                     keepBest = F)
  
  # todo improvements:
  #   - let batch size depend on highest prime factor in n_sample - done but critical
}
toc()
sink(NULL)

sink(file = './log_lstm_2l.txt', split = T, append = F)
tic('2 layers loop')
for (i in 1:n){
  # fit model
  inflation[['lstm_fullsample_2l']][[i]] <- 
    k_fullsample_2l(data = inflation[['lstm_data']][[i]]$train$train_norm, 
                    n_steps = min(inflation[['aropti']][[i]]*2, 12), 
                    n_feat = 1, 
                    nodes = 75, 
                    size_batch = 1, 
                    epochs = 40, 
                    ES = T, 
                    keepBest = T)
  # save model somewhere on disk
  save_model_hdf5(object = inflation[['lstm_fullsample_2l']][[i]]$model_fitted, 
                  filepath = file.path(models_dir,
                                       paste0(inflation[['names']][[i]],
                                              '_2l_fullsample.h5')
                                       )
                  )
}
toc()
sink(NULL)

##### If models are fitted externally, load in those files

# for (i in 1:n){
#   inflation[['lstm_fullsample']][[i]]$model_fitted <-
#     keras::load_model_hdf5(filepath = file.path(paste0(models_dir,'_4k_n75/'),
#                                                 paste0(inflation[['names']][[i]],
#                                                        ' fullsample.h5')),
#                                                 compile = T)
# }



##### Online predictions #######################################################

for (i in 1:n){
  inflation[['lstm_online_pred']][[i]] <- online_pred(model_fitted = inflation[['lstm_fullsample']][[i]], 
                                                      model_type = 'model_fitted',
                                                      data_train = inflation[['lstm_data']][[i]],
                                                      horizon = 20)
  
  inflation[['plot_lstm_full']][[i]] <- ggplot(data = inflation[['lstm_online_pred']][[i]])+
                                          geom_line(aes(x = date, y = value, colour = label))+
                                          theme_minimal() + xlab(label = element_blank()) + 
                                          ylab(element_blank()) + ggtitle(inflation$names[[i]]) + 
                                          theme(legend.position = 'bottom', 
                                                legend.title = element_blank())+
                                          guides(colour = guide_legend(nrow = 1))
  
  plot(inflation[['plot_lstm_full']][[i]])
  
}


##### Split data in 10y chuncks ################################################


inflation[['lstm_splits_10y']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                        initial = fm_apply(4*10, n),
                                                        assess = fm_apply(0, n),
                                                        cumulative = fm_apply(F, n),
                                                        skip = fm_apply(0, n),
                                                        lag = inflation[["aropti"]]), 
                                              .f = rsample::rolling_origin)

inflation[['lstm_increm_splits']] <- future_pmap(.l = list(data = sapply(pi, FUN = function(x) {list(na.omit(x))}),
                                                           initial = fm_apply(4*10, n),
                                                           assess = fm_apply(0, n),
                                                           cumulative = fm_apply(F, n),
                                                           skip = fm_apply(0, n),
                                                           lag = inflation[["aropti"]]), 
                                                 .f = rsample::rolling_origin)