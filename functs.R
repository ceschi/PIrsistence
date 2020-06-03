##### Specifically designed functions ####



##### I - create directories ###################################################
working_directory <- getwd()
temp_dir <- 'downloaded_files'
data_dir <- 'processed_data'
graphs_dir <- 'plots'
models_dir <- 'models'

temp_dir <- file.path(working_directory, temp_dir)
data_dir <- file.path(working_directory, data_dir)
graphs_dir <- file.path(working_directory, graphs_dir)
models_dir <- file.path(working_directory, models_dir)

options(warn=-1) # turns off warnings momentarily
dir.create(temp_dir)
dir.create(data_dir)
dir.create(graphs_dir)
dir.create(models_dir)
options(warn=0) # turns warnings back on




# A file to gather all home made functions with relative descriptions


instant_pkgs <- function(pkgs) { 
  ## Function loading or installing packages in
  ## current R instance.
  ## Developed by Jaime M. Montana Doncel - V1

  
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss, dependencies = T, Ncpus = 2, INSTALL_opts = '--no-lock')
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss, dependencies = T, Ncpus = 2, INSTALL_opts = '--no-lock')
  }
  
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) suppressPackageStartupMessages(library(need_to_attach[i], character.only = TRUE))
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded!\n")
  }
}


rollm <- function(df, formula){
  
  # function to extract and store coefficients 
  # and double SD in a named row tibble
  
  
  # estimates the linear model
  lmod <- summary(lm(data=df, formula=formula))
  
  # extracts point estimates and 2*SD (+- 95%),
  # put info in named row tibble dropping 
  # intercept info from first column
  
  cofs <- as_tibble(coefficients(lmod)[2:(lmod %>% coefficients() %>% 
                                            t() %>% ncol()),1] %>% t(),
                    .name_repair = 'minimal')
  SD2 <- as_tibble(2*coefficients(lmod)[2:(lmod %>% coefficients() %>% 
                                            t() %>% ncol()),2] %>% t(),
                   .name_repair = 'minimal')
  
  # adds suffix for bands
  names(SD2) <- paste0(names(SD2), '.SD2')
  
  # merges in one row with names
  estim <- cbind(cofs, SD2)
  
  # outputs
  return(estim)
}

rolloop <- function(df, window=8, lags=1, interc = T){
  
  # width of the rolling window
  window <- as.integer(window)
  
  # select lags 
  k <- as.integer(lags)
  
  # lags the time series, names it, cuts out NAs
  df <- df %>% lagger(laag=k, na.cut=T)
  # and creates related formula
  formulae <- formula.maker(df, 
                            df %>%  names(.) %>% first(),
                            intercept = interc)
  
  # computes point estimates and 2SD
  # stocks in a dataframe for convenience
  regs <-rollapply(as.data.frame(df),
                   width=window,
                   by.column = F,
                   FUN=function(x, formula) rollm(df=as.data.frame(x), formula=formulae))
  
  # converts and dates the regressions
  regs <- xts(regs, frequency=4, 
              order.by=index(df)[window:length(index(df))])
  
 return(regs)
}

make_stars <- function(x){
  # ancillary function for
  # printing stars alongside
  # with converted parameters
  

  # pre-allocate 
  signif <- NULL
  
  if (x < .001) {
    signif <- as.factor('***')
  }else if (x < .01 & x >= .001){
    signif <- as.factor('**')
  }else if (x < .05 & x >= .01){
    signif <- as.factor('*')
  }else if (x < .1 & x >= .05){
    signif <- as.factor('.')
  }else if (x>=.1){
    signif <- as.factor('')
  }
  
  return(signif)
}

lagger <- function(series, laag, na.cut=F){
  # Takes a time series and creates a matrix with given number
  # of lags, also generating appropriate names
  
  
  matrix <- as.data.frame(matrix(ncol=laag+1, nrow=nrow(series)))
  for (i in 1:laag+1){
    matrix[,i] <- stats::lag(series, k=(i-1))
  }
  names(matrix) <- c(names(series), paste(names(series), 1:laag, sep='.'))
  matrix[, 1] <- series
  matrix <- as.xts(matrix, order.by=index(series))
  
  # conditional to remove NAs from output
  if (na.cut){
    matrix <- na.omit(matrix)
  }
  
  # output
  return(matrix)
}

# lagger_bis benchmarks better
# roughly ten times faster

lagger_bis <- function(series, lag, na.cut=F){
  # Takes a time series and creates a matrix with given number
  # of lags, also generating appropriate names
  # 
  matrix <- embed(as.matrix(series), lag+1)
  matrix <- as.data.frame(matrix)
  names(matrix) <- c(names(series), paste(names(series), 1:lag, sep='.'))
  
  # conditional to remove NAs from output
  if (na.cut){
    matrix <- na.omit(matrix)
  }
  
  # output
  return(matrix)
}

formula.maker <- function(df, y, intercept = T){
  # provided with a df and a dependent variable name
  # this generates a formula for estimation in R, y is the 
  # dependent variable, all the others are considered
  # independent and explanatory ones
  
  if (intercept){
    fomu <- as.formula(paste(y, 
                           paste(names(df)[names(df)!=y], collapse='+'),
                           sep='~'))
  } else {
    fomu <- as.formula(paste(y,
                             paste(c(0,names(df)[names(df)!=y]), collapse='+'),
                             sep='~'))
                }
  
  
  attr(fomu, which='.Environment') <- .GlobalEnv
  return(fomu)
}


############ PIRSISTENCE FUNCTIONS #############################################

auto.reg <- function(data, lags = 1, interc = T){
  
  # function to estimate AR(lags)
  
  transf_data <- lagger(series = data,
                        laag = lags,
                        na.cut = F)
  
  model_formula <- formula.maker(df = transf_data,
                                 y = first(names(transf_data)),
                                 intercept = interc)
  linear_model <- lm(formula = model_formula,
                     data = transf_data)
  
  return(linear_model)
}

auto.reg.sum <- function(data, lags = 1, interc = T){
  
  if (!require(broom)) {install.packages('broom'); library(broom)}
  # not necessary as already in tidyverse
  # if (!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
  # if (!require(magrittr)) {install.packages('magrittr'); library(magrittr)}
  
  # function to estimate AR(lags) and sum over parameters

  transf_data <- lagger(series = data,
                        laag = lags,
                        na.cut = F)
  

  
  model_formula <- formula.maker(df = transf_data,
                                 y = first(names(transf_data)),
                                 intercept = interc)
  
  linear_model <- lm(formula = model_formula,
                     data = transf_data)
  
  output <- broom::tidy(linear_model)
  
  coef_sum <- output %>% filter(term != '(Intercept)') %>% dplyr::select(estimate) %>%  sum()
  
  return(coef_sum)
}




rolloop.sum <- function(df, window, lags = 1, interc = T){
  
  # remove troublesome NAs
  df_na <- na.omit(df)
  
  # computes point estimates
  # stocks in a dataframe for convenience
  regs <-rollapply(df_na,
                   # as.data.frame(df),
                   width=window,
                   by.column = F,
                   FUN = auto.reg.sum,
                   lags = lags,
                   interc = interc)
  
  # # converts and dates the regressions
  # regs <- xts(regs, frequency=4, 
  #             order.by=index(df_na)[window:length(index(df_na))])
}




# TRACKING PERSISTENCE OVER TIME #
persistence_ridges <- function(tseries, window = 24, lags = 8){
  # requires zoo, broom
  if (!require(zoo))    {install.packages('zoo');   library(zoo)}
  if (!require(broom))  {install.packages('broom'); library(broom)}
  
  # check out the nature of the input
  # throw an error if it's not time series class
  if (!(class(tseries)=='ts' || class(tseries)=='xts' || class(tseries) == 'zoo')) error('Wrong object, please provide a time series object (ts, zoo, xts).')
  if (window<=lags*2) warning('\nWrong window/lag sizes: \nto get meaningful estimates window width should be at least twice the lags.')
  
  # define function to be applied rolling over
  bloc_ar <- function(tseries, lags = 8, interc = F, last){
    
    # save out last observation of the series
    # will be the identifier later on
    # last <- time(tseries)[length(tseries)]
    
    # generate a matrix with lags+1 columns
    # to have original series + lagged cols
    # 
    # It outputs a flat matrix, its 
    # length is cut down by lags
    mat_lag <- embed(tseries, lags+1)
    
    
    # estimate linear model without intercept,
    # store the results
    estlm <- lm(data = as.data.frame(mat_lag),
                formula = formula.maker(as.data.frame(mat_lag), 'V1', intercept = interc))
    
    # flip in tidy format the lm output
    # and delete the "statistic" col
    est_tidy <- broom::tidy(estlm)
    est_tidy$statistic <- NULL
    
    # gather all in a dated dataframe:
    # lengths = 8
    # width   = 5
    # names = c('last.date', 'term', 'estimate', 'std.error', 'p.value')
    col <- data.frame(last.date = rep(last, lags),
                      est_tidy)
    
    # output the resulting df
    return(col)
  }
  
  # remove all NAs - experimental
  tseries <- na.omit(tseries)
  
  # this object out_fin will accommodate 
  # the results, iteration by iteration
  # add names and preallocate cells
  out_fin <- matrix(nrow = (length(tseries)-window+1)*lags, ncol = 5)
  out_fin <- as.data.frame(out_fin)
  names(out_fin) <- c('last.date', 'term', 'estimate', 'std.error', 'p.value')
  
  for (i in 1:(length(tseries)-window+1)){
    
    last_date <- time(tseries)[(i+window-1)]
    
    col_fin <- bloc_ar(tseries = tseries[i:(i+window-1)],
                       lags = lags,
                       interc = F,
                       last = last_date)
    
    out_fin[((i-1)*lags + 1):((i)*lags),] <- col_fin
    # out_df <- rbind(out_df,col)
  }
  
  out_fin$term <- as.numeric(
    gsub(pattern = '[V]',
         x = out_fin$term,
         replacement = '')
  ) - 1
  
  return(out_fin)
  
}


# Markov Switching model with optimal lags
ms_aropti <- function(df, lags, states){
  # adapt the dataset creating lags
  data <-  lagger_bis(series = df, 
                      lag = lags)
  
  # estimate a linear model
  l_model <- lm(data = data)
  
  # run the MS
  estimate <- msmFit(#data = data,
    object = l_model,
    k = states,
    sw = rep(T, 1+1+lags),)
  # output results
  return(estimate)
  
}


# plot rolling estimates for AR1
plot_roller <- function(df, names, path){
  po <- ggplot(data=df,
               aes(x=index(df),
                   y=df$Var.1))+
    # plot the above with line geom, in black
    geom_line(colour='black', size=1)+
    # adds upper confidence band in red
    geom_line(aes(y=(df$Var.1 + df$.SD2)),
              colour='red')+
    # adds lower confidence band in red
    geom_line(aes(y=(df$Var.1 - df$.SD2)),
              colour='red')+
    # adds unit root line
    geom_line(aes(y=1), colour='black', size=.8)+
    geom_line(aes(y=0), colour='black', size=.8)+
    # plot makeup
    geom_smooth(method='loess', colour='blue', formula = 'y~x')+
    scale_x_yearqtr(format='%Y Q%q')+theme_minimal()+
    scale_y_continuous()+xlab(' ') + ylab(paste0('AR(1) coeff. estimates')) + 
    ggtitle(paste0(names, ' - 1 exogenous lag'))
  
  
  
  # saves the plots in given path
  ggsave(paste0(names, ' - AR(1) coeff estimates.pdf'),
         plot = po,
         device='pdf',
         path = path,
         height=8/2, width=14.16/2, units='in')
  
  
  return(po)
}

# plots summed coefficients of optimal AR
plot_autoregsum <- function(df, names, path, laags){
  po <- ggplot(data=df,
               aes(x=index(df),
                   y=df[,1]))+
    # plot the above with line geom, in black
    geom_line(colour='black', size=1)+
    # adds unit root line
    geom_line(aes(y=1), colour='black', size=.8)+
    geom_line(aes(y=0), colour='black', size=.8)+
    # plot makeup
    geom_smooth(method='loess', colour='blue', formula = 'y~x')+
    scale_x_yearqtr(format='%Y Q%q')+theme_minimal()+
    scale_y_continuous()+xlab(' ') + ylab(paste0('AR(',laags,') coeff. estimates sum')) + 
    ggtitle(paste0(names, ' - ', laags, ' optimal lags: sum of coefficients')) 
  
  
  # save plot
  
  ggsave(paste0(names, ' - AR(',laags,') coeff estimates sum.pdf'),
         plot = po,
         device='pdf',
         path = path,
         height=8/2, width=14.16/2, units='in')
  
  return(po)
  
}

# plots ridges for AR
plot_ridges <- function(df, nam, laags, path){
  out <- ggplot(data = df)+
    geom_ridgeline_gradient(aes(x = term,
                                y = as.factor(last.date),
                                height = estimate,
                                group = as.factor(last.date),
                                fill = p.value),
                            min_height = -2) +
    scale_fill_viridis(name = 'P-values',option = "C", direction = 1) +
    ggtitle(paste0('Evolving persistence - ',
                   nam,
                   ' ',
                   laags,
                   ' end. lags')) +
    xlab('Lag order') + ylab(' ') + theme_minimal()
  
  
  ggsave(paste0(nam, ' - AR(',laags,') acf.pdf'),
         plot = out,
         device = 'pdf',
         path = path,
         # extra height needed for full display
         height = 100,
         width = 14.16,
         units = 'in',
         limitsize = FALSE
  )
  
  return(out)
}


# plot Markov Switching results
plot_msm <- function(ms_model, nam, laags, path){
  
  # setting device size
  # mar sets margings
  # cex.main scales title to 70%
  par(mar = c(1,1,2.85,1), cex.main = .70)
  
  # store actual plot
  # it's automatically printed
  plot_out <- plotProb(ms_model, which = 2)
  
  # fix title
  title(paste0(flag___ms, '-state MS regimes for ', nam, ' with ', laags, ' lags'), line = 2.3)
  
  # copy dev output to file (pdf)
  dev.copy(pdf, height=8/1.5, width=14.6/1.5,
           file.path(path,
                     paste0(nam, ' ', flag___ms, '-state MSM.pdf')
           )
  ) %>% invisible() # just to remove annoying output
  
  # shut down device, comment for keeping the plot
  invisible(dev.off())
  
  # output
  return(plot_out)
}


# LSTM functions ----------------------------------------------------------

data_prepper <- function(data, train = 1, test = NULL){
  # function to rescale data to normal values
  # that are required for the LSTM model fit.
  # It stores first and second moments for later
  # scaling back to observed data. All outputs are
  # stored in a list, usually TS properties (index) are 
  # preserved.
  
  
  # train and test shall be expressed in percentage terms
  if (train>1){stop('Train should be expressed in percentage not intergers.')}
  
  # preallocate list
  output <- list()
  
  # remove NAs
  data <- data[!is.na(data)]
  if (is.null(dim(data))){
    data <- array(data = data, dim = c(length(data), 1))
  }
  
  # train/rep
  len <- dim(data)[1]
  l_train <- floor(len*train)
  l_test <- ifelse(is.null(test), 0, (len-l_train))
  
  # split data
  data_train <- data[1:l_train, ]
  
  # rescale train
  output$train[['mean']] <- mean(data_train)
  output$train[['sd']] <- sd(data_train)
  output$train[['train_norm']] <- array(data = (data_train - output$train[['mean']])/output$train[['sd']],
                                        dim = c(l_train, 1))
  output$train[['n_sample']] <- l_train

  
  if (train != 1){
    # condition on test subsample, whether present
    data_test <- data[(l_train+1):len, ]
    # rescale test
    output$test[['mean']] <- mean(data_test)
    output$test[['sd']] <- sd(data_test)
    output$test[['test_norm']] <- array(data = (data_test - output$test[['mean']])/output$test[['sd']],
                                        dim = c(l_test, 1))
    output$train[['n_sample']] <- l_test
  }
  
  return(output)
}


k_fullsample <- function(data, 
                         n_steps, 
                         n_feat = 1, 
                         # model_compiled,
                         nodes = 50,
                         # size_batch = 1, 
                         epochs = 2000,
                         ES = F,
                         keepBest = F){
  
  # Function to fit a model on the whole sample of data;
  # it takes care of lagging & reshaping the data according to parameters
  # passed and to declare a model with one layer of LSTM and a final
  # dense layer. It requires 'keras' and pipes. It preserves the 
  # time dimension dropping the NAs generated after lagging.
  # Optionally it accommodates the early stopping callback
  # for which the epochs are then the patience limit.
  # Output contains history, model, and optionally time index.
  
  # Keras usually run in parallel by default, call this function 
  # sequentially to avoid nested parallelism.
  
  
  
  require(magrittr)
  require(keras)
  require(tictoc)
  require(numbers)
  
  # data come in as a simple TS,
  # it comes then with n of observations (n_sample)
  # and must be lagged according to n_steps.
  # The number of features is 1 by default, can be varied tho
  
  # we will drop as many obs as many lags we include
  n_sample <- nrow(data) - n_steps
  n_feat <- ncol(data)
  
  # preserve the time index of data
  # for later use
  if (is.xts(data)){
    time_index <- time(data)[(n_steps+1):length(data)]
    
  }
  
  # highest prime factor in the number of obs
  batch_prime <- numbers::primeFactors(n_sample)
  batch_prime <- batch_prime[length(batch_prime)]
  
  # KIM batch must mod train and test samples!
  size_batch <- batch_prime
  
  # embed automates lags and turns into lower
  # object matrix/array: first col is original series
  # second to end are lags
  data_lagged <- embed(x = as.matrix(data), dimension = (n_steps+1))
  
  
  # NB: y must be 2D array/matrix
  y_data_arr <- array(data = data_lagged[,1],
                      dim = c(n_sample, n_feat))
  
  # X must be 3D for a stateful LSTM
  x_data_arr <- array(data = data_lagged[,-1],
                      dim = c(n_sample, n_steps, n_feat))
  
  model_compiled <- keras_model_sequential()
  model_compiled %>%
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = F,
               stateful = T,
               batch_size = size_batch,
              ) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  model_online <- keras_model_sequential() %>% 
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = F,
               stateful = T,
               batch_size = 1,
              ) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  
  tictoc::tic('Model estimation')
  if (ES){
    # estimate with early stopping
    history <- fit(object = model_compiled, 
                   y = y_data_arr, 
                   x = x_data_arr,
                   verbose = 2,
                   shuffle = F,
                   callbacks = list(
                     callback_early_stopping(monitor = 'val_loss',
                                             mode = 'auto',
                                             patience = epochs,
                                             min_delta = 1e-5, 
                                             restore_best_weights = keepBest)
                   ),
                   epochs = epochs,
                   validation_split = .1,
                   batch_size = size_batch)
  } else {
    # estimate with given number of epochs
    history <- fit(object = model_compiled, 
                   y = y_data_arr, 
                   x = x_data_arr, 
                   epochs = epochs, 
                   verbose = 2,
                   shuffle = F,
                   batch_size = size_batch)
  }
  tictoc::toc()
  
  out <- list()
  out[['model_fitted']] <- model_compiled
  out[['history']] <- history
  out[['model_batch']] <- batch_prime
  out[['model_weights']] <- keras::get_weights(model_compiled)
  out[['model_online']] <- keras::set_weights(object = model_online, 
                                              weights = out[['model_weights']])
  if (is.xts(data)){
    out[['time_index']] <- time_index
  }
  
  return(out)
}

extra_layers <- function(nodes_list, options){
  
  # function to automate the layering of models 
  # in keras. It outputs a model ready to compile.
  
  require(keras)
  require(magrittr)
  
  model <- keras::keras_model_sequential()
  
  for (l in 1:length(nodes_list)){
    # the first layer needs extra info on input
    # and batch size.
    if (l==1){
      model %>% 
        layer_lstm(
          input_shape = options$input_shape,
          batch_size = options$size_batch,
          units = nodes_list[[l]],
          return_sequences = options$ret_sequences[[i]],
          stateful = options$stateful[[i]]
          )
      
    }else{
      model %>% 
        layer_lstm(
          units = nodes_list[[l]],
          return_sequences = options$ret_sequences[[i]],
          stateful = options$stateful[[i]]
          )
    }
  }
  
  model %>% layer_dense(units = 1)
  
  return(model)
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
  
  # retrieve input data and fitted model
  input <- data_train$train[['train_norm']]
  model_online <- model_fitted[['model_online']]
  
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
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train', 
                                   date = time(data_train$train[['train_norm']])) %>% 
                        rename(value = V1) %>% 
                        mutate(value = as.numeric(value)),
                      
                      pred %>% 
                        as_tibble %>% 
                        add_column(label = 'forecast', 
                                   date = time_preds) %>% 
                        rename(value = V1) %>% 
                        mutate(value = as.numeric(value)))
  } else {
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train') %>% 
                        rename(value = V1), 
                      pred %>% 
                        as_tibble %>% 
                        add_column(label = 'forecast') %>% 
                        rename(value = V1))
  }
  
  # reconversion to values
  forecast$value <- forecast$value*data_train$train[['sd']] + data_train$train[['mean']]
  
  return(forecast)
  
}




##### Packages Loader #####

pkgs <- c(
  'broom',
  'devtools',
  'furrr', 
  'future',
  'ggridges', 
  'glue',
  'MSwM',
  'stargazer',
  'strucchange',
  'tictoc',
  'viridis',
  'keras',
  'numbers'
  )
# fill pkgs with names of the packages to install

instant_pkgs(pkgs)

# devtools::install_github('sboysel/fredr')
devtools::install_github('ceschi/urcabis')
# devtools::install_version("readxl", version = "1.0.0")
# library(urcabis) # for when the package will be duly updated (pull request)



#### housekeeping ####
rm(pkgs)
