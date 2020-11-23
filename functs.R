##### Specifically designed functions ####



##### I - create directories ###################################################
working_directory <- getwd()
temp_dir <- 'downloaded_files'
data_dir <- 'processed_data'
graphs_dir <- 'plots'
models_dir <- 'models'
rds_dir <- 'rds_files'

# detailed folders for LSTM
lstm_dir <- 'plots_lstm'
l1_dir <- 'plots_1lfore'
l2_dir <- 'plots_2lfore'
chunks_dir <- 'plots_chunks'
tables <- 'tables'
increm_dir <- 'plots_increm'
rolling_dir <- 'plots_rollwind'

ar1_dir <- 'plots_ar1'
ark_dir <- 'plots_ark'
vars_dir <- 'plots_vars'
acf_dir <- 'plots_acf'

# bayes
lines_dir <- 'plots_lines'
distro_dir <- 'plots_distro'




temp_dir <- file.path(working_directory, temp_dir)
data_dir <- file.path(working_directory, data_dir)
graphs_dir <- file.path(working_directory, graphs_dir)
models_dir <- file.path(working_directory, models_dir)
l1_dir <- file.path(graphs_dir, lstm_dir, l1_dir)
l2_dir <- file.path(graphs_dir, lstm_dir, l2_dir)
chunks_dir <- file.path(graphs_dir, lstm_dir, chunks_dir)
rolling_dir <- file.path(graphs_dir, lstm_dir, rolling_dir)
tables <- file.path(graphs_dir, lstm_dir, tables)
increm_dir <- file.path(graphs_dir, lstm_dir, increm_dir)
ar1_dir <- file.path(graphs_dir, ar1_dir)
ark_dir <- file.path(graphs_dir, ark_dir)
vars_dir <- file.path(graphs_dir, vars_dir)
acf_dir <- file.path(graphs_dir, acf_dir)
rds_dir <- file.path(working_directory, rds_dir)
lines_dir <- file.path(graphs_dir, lines_dir)
distro_dir <- file.path(graphs_dir, distro_dir)

options(warn=-1) # turns off warnings momentarily
dir.create(temp_dir)
dir.create(data_dir)
dir.create(graphs_dir)
dir.create(models_dir)
dir.create(ar1_dir)
dir.create(ark_dir)
dir.create(acf_dir)
dir.create(vars_dir)
dir.create(rds_dir)

dir.create(file.path(graphs_dir, lstm_dir))
dir.create(l1_dir)
dir.create(l2_dir)
dir.create(chunks_dir)
dir.create(rolling_dir)
dir.create(tables)
dir.create(increm_dir)
dir.create(lines_dir)
dir.create(distro_dir)
options(warn=0) # turns warnings back on

##### II - general purpose functions ###########################################

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
  #' *Add here storage of intercept with related SDs*
  
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

rollmbis <- function(.df, .formula){
  invisible(require(dplyr))
  invisible(require(broom))
  
  lmod <- broom::tidy(lm(data = .df,
                         formula = .formula))
  
  # check whether .formula has intercept (T=1)
  # attr(terms(ff), 'intercept')
  
  out <- lmod %>% 
    dplyr::select(term, estimate, std.error) %>%
    dplyr::mutate(term = tolower(gsub(pattern = '\\(',
                                      replacement =  '', 
                                      x = gsub('\\)', '', term)
                                      ))) %>% 
    tidyr::pivot_wider(id_cols = term,
                       names_from = term,
                       values_from = c('estimate', 'std.error'),
                       names_sort = T)
  
  return(out)
  
}

rolloop <- function(df, window=8, lags=1, interc = T){
  
  #' *This should be already able to accomodate intercept*
  
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

fm_apply <- function(foo, n){
  # shorthand to use in future_pmap
  # calls for atom foo
  return(sapply(rep(foo, n), list))
}


############ II - frequentist part #############################################

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
  
  invisible(require(broom))
  invisible(require(dplyr))
  invisible(require(magrittr))
  
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
  
  #' *hand this part for storing intecept*
  coef_sum <- output %>% 
    filter(term != '(Intercept)') %>% 
    dplyr::select(estimate) %>%  
    sum()
  
  if (interc){
    coef_sum_se <- linear_model %>% 
                    vcov() %>%
                    .[-1,-1] %>% 
                    sum() %>%
                    sqrt()
  }else{
    coef_sum_se <- linear_model %>% 
      vcov() %>%
      sum() %>% 
      sqrt()
  }
  
  coef_sum <- cbind(coef_sum, coef_sum_se)
  
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

persistence_ridges <- function(tseries, window = 24, lags = 8){
  # requires zoo, broom
  if (!invisible(require(zoo)))    {install.packages('zoo');   invisible(library(zoo))}
  if (!invisible(require(broom)))  {install.packages('broom'); invisible(library(broom))}
  
  # check out the nature of the input
  # throw an error if it's not time series class
  if (!(class(tseries)=='ts' || class(tseries)=='xts' || class(tseries) == 'zoo')) warning('Wrong object, please provide a time series object (ts, zoo, xts).')
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
               aes(x=index(df) %>% as.yearqtr(),
                   y=df$Var.1))+
    # plot the above with line geom, in black
    geom_line(colour='black', size=1)+
    # add confidence interval
    geom_ribbon(aes(ymax = (df$Var.1 + df$.SD2),
                    ymin = (df$Var.1 - df$.SD2)),
                colour = 'grey',
                size = .25,
                alpha = .1) +
    # adds unit root line
    geom_line(aes(y=1), colour='black', size=.5)+
    geom_line(aes(y=0), colour='black', size=.5)+
    # plot makeup
    geom_smooth(method='loess', colour='blue', formula = 'y~x', se = F)+
    scale_x_yearqtr(format='%Y Q%q')+theme_minimal()+
    scale_y_continuous()+xlab(' ') + ylab(paste0('AR(1) coeff. estimates')) + 
    ggtitle(paste0(names %>% noms_tt(), ' - 1 exogenous lag'))+
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5))
  
  
  
  # saves the plots in given path
  ggsave(paste0(names %>% noms(), '_AR(1)_coeff.pdf'),
         plot = po,
         device='pdf',
         path = path,
         units='in',
         width = 8,
         height = 8*9/16)
  
  
  return(po)
}

# plots summed coefficients of optimal AR
plot_autoregsum <- function(df, names, path, laags){
  po <- ggplot(data=df,
               aes(x=index(df) %>% as.yearqtr(),
                   y=df[,1]))+
    # plot the above with line geom, in black
    geom_line(colour='black', size=1)+
    # adds unit root line
    geom_line(aes(y=1), colour='black', size=.5)+
    geom_line(aes(y=0), colour='black', size=.5)+
    # plot makeup
    geom_smooth(method='loess', colour='blue', formula = 'y~x', se = F)+
    # scale_x_yearqtr(format='%Y Q%q')+ 
    # scale_x_date(date_labels = '%Y', breaks = '10 years') + 
    scale_x_yearqtr(format = '%Y Q%q')+
    theme_minimal()+
    scale_y_continuous()+xlab(' ') + ylab(paste0('AR(',laags,') coeff. estimates sum')) + 
    ggtitle(paste0(names %>% noms_tt(), ' - ', laags, ' optimal lags: sum of coefficients')) +
    # add ribbon style standard errors
    geom_ribbon(aes(ymin = (df[,1] - df[,2]),
                    ymax = (df[,1] + df[,2])),
                size = .25, colour = 'grey', alpha = .1)+
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5))
  
  
  # save plot
  
  ggsave(paste0(names %>% noms(), '_AR(',laags,')_coeffsum.pdf'),
         plot = po,
         device='pdf',
         path = path,
         units='in',
         width = 8,
         height = 8*9/16)
  
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
                   nam %>% noms_tt(),
                   ' ',
                   laags,
                   ' end. lags')) +
    xlab('Lag order') + ylab(' ') + theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(nam %>% noms(), '_AR(',laags,')acf.pdf'),
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

# nicer names from human readable strings
noms <- function(x){
  # if (!is.character(x)) stop('\nNot a string')
  invisible(require(magrittr))
  x %>% 
    tolower() %>% 
    gsub(x = ., 'revised ', '') %>% 
    gsub(x = ., 'rev_', '') %>% 
    gsub(x = ., ',', '') %>% 
    gsub(x = ., ' ', '_') %>% 
    trimws() %>% 
    return()
}

noms_tt <- function(x){
  invisible(require(magrittr))
  
  x %>% 
    gsub(x = ., '^[0-9]{1,2}y_', '') %>% 
    gsub(x = ., 'Revised ', '') %>% 
    gsub(x = ., ', no FE ', ' core ') %>% 
    gsub(x = ., 'pch', '') %>% 
    trimws() %>% 
    return()
}

##### III - treat LSTM predictions output ######################################

chunk_regs <- function(regs_list, regs_list_sum, ar_lags_sum, fore_horiz){
  # create a tibble from chunks' regressions
  # and labels for start/end
  
  invisible(require(magrittr))
  
  tidyout <- function(onereg, fore_horiz){
    # extract results from lm and 
    # add label. Does the heavy lifting.
    
    # start date
    start <- onereg %>% 
      model.frame() %>% 
      rownames() %>% 
      head(1) %>% 
      as.Date()
    # adjust for one lag
    # true only for this AR1 case
    start <- start - base::months(3, F)
    
    end <- onereg %>% 
      model.frame() %>% 
      rownames() %>% 
      tail(fore_horiz + 1) %>% 
      head(1) %>% 
      as.Date()
    
    labl <- paste0(lubridate::year(start), quarters(start), ' - ', lubridate::year(end), quarters(end))
    
    out <- broom::tidy(onereg) %>% tibble::add_column(chunk = labl)
    
    return(out)
  }
  
  
  tidycoefs <- lapply(X = regs_list, FUN = tidyout, fore_horiz)
  
  out_ar1 <- tidycoefs %>% dplyr::bind_rows()
  
  
  # second part for regsum
  tidyout_sum <- function(ar1, regs_list_sum, ar_lags_sum, fore_horiz){
    # extract results from lm and 
    # add label. Does the heavy lifting.
    
    # start date
    start <- ar1 %>% 
      model.frame() %>% 
      rownames() %>% 
      head(1) %>% 
      as.Date()
    # adjust for lags wrt ar1
    start <- start - base::months(3, F)*(ar_lags_sum-1)
    
    end <- ar1 %>% 
      model.frame() %>% 
      rownames() %>% 
      tail(fore_horiz + 1) %>% 
      head(1) %>% 
      as.Date()
    
    labl <- paste0(lubridate::year(start), quarters(start), ' - ', lubridate::year(end), quarters(end))
    
    out <- tibble::tibble(ar_sum = regs_list_sum[,1],
                          ar_sum_se = regs_list_sum[,2],
                          chunk = labl,
                          k_lags = ar_lags_sum)
    
    return(out)
  }
  
  # out_ark_sum <- furrr::future_pmap(.l = list(ar1 = regs_list, 
  #                                             ar_lags_sum = fm_apply(ar_lags_sum, length(regs_list_sum)), 
  #                                             fore_horiz = fm_apply(fore_horiz, length(regs_list_sum)), 
  #                                             regs_list_sum = regs_list_sum), 
  #                                   .f = tidyout_sum) %>% 
  #   dplyr::bind_rows()
  
  out_ark_sum <- purrr::pmap(.l = list(ar1 = regs_list, 
                                       ar_lags_sum = fm_apply(ar_lags_sum, length(regs_list_sum)), 
                                       fore_horiz = fm_apply(fore_horiz, length(regs_list_sum)), 
                                       regs_list_sum = regs_list_sum), 
                             .f = tidyout_sum) %>% 
    dplyr::bind_rows()
  
  
  out <- list(ar1 = out_ar1,
              ark_sum = out_ark_sum)
  
  return(out)
  
}

plot_chunkregs_bar <- function(chunk_regs_obj, graphs_dir. = graphs_dir, name){
  
  # bar plot/save for the AR1 models on chunks of data
  # can be used also for rolling windows and increasing windows
  # but DOES require lm object - hence need adaptation for autosum fct
  
  # part 1: bar plot for AR1 on each chunk
  
  # number of chunks
  len <- chunk_regs_obj$ar1 %>% 
    select(chunk) %>% 
    unique() %>% 
    nrow()
  
  # setup title
  tt <- paste0(name %>% noms_tt(),
               ': ',
               len,
               ' chunks')
  
  jj <- name %>% 
    noms() %>% 
    paste0('_ar1_chunks.pdf')
  
  # make plot, filtering out intercept
  plt <- chunk_regs_obj$ar1 %>% 
    filter(term != '(Intercept)') %>% 
    ggplot(aes(x = chunk, 
               y = estimate, 
               group = chunk))+
    geom_col(alpha = .5) +
    geom_errorbar(aes(ymin = (estimate - std.error), 
                      ymax = (estimate + std.error)), 
                  width = .2)+
    ggtitle(tt) + theme_minimal() + ylab('AR(1) coefficient') + xlab('Time periods') + 
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0,
                                     size = rel(.65)))
  
  
  
  ggsave(plot = plt, 
         filename = file.path(graphs_dir., 
                              jj),
         device = 'pdf', 
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  
  # Part 2: barplot for each chunk, sum of coefficients
  
  labely <- chunk_regs_obj$ark_sum %>% 
    select(k_lags) %>% 
    unique() %>% 
    paste0('sum of AR(', ., ') coefficients')
  
  tt <- paste0(name %>% noms_tt(),
               ': ', len, ' chunks')
  
  jj <- name %>% 
    noms() %>% 
    paste0('_ark_sum_chunks.pdf')
  
  
  plt_sum <- chunk_regs_obj$ark_sum %>% 
    ggplot(aes(x = chunk, y = ar_sum, group = chunk)) + 
    geom_col(alpha = .5) + theme_minimal() + ylab(labely) + xlab('Time periods') + 
    ggtitle(tt) + 
    geom_errorbar(aes(ymin = (ar_sum - ar_sum_se), ymax = (ar_sum + ar_sum_se)), width = .2)+
    theme(axis.text = element_text(size = rel(1.5)), 
          legend.text = element_text(size = rel(1.5)), 
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0,
                                     size = rel(.65)))
  
  ggsave(plot = plt_sum, 
         filename = file.path(graphs_dir., 
                              jj),
         device = 'pdf', 
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  plt_list <- list(ar1 = plt, 
                   ark_sum = plt_sum)
  
  
  return(plt_list)
}

chunk_stargazer <- function(ar1, chunk_out, name, pathout = graphs_dir){
  
  # function to output regressions results with stargazer from 
  # tidy object and lm one
  
  # store periods in character form
  mod_labels <- chunk_out$ar1 %>% 
    dplyr::select(chunk) %>% 
    unique() %>% 
    unlist() 
  
  # create file name 
  filename <- name %>% noms() %>% paste0('_chunkreg_tab.tex')
  name_tt <- name %>% noms_tt()
  
  # turn into appropriate path
  destination <- file.path(pathout, filename)
  
  # regex stuff
  strin <- paste0('\\multicolumn\\{', length(mod_labels), '\\}\\{c\\}\\{', name, '\\}')
  repla <- paste0('\\multicolumn\\{', length(mod_labels), '\\}\\{c\\}\\{', name_tt, '\\}')
  
  # ll <- length(mod_labels)
  # notes <- paste0('(', 1:ll, '): ', mod_labels) %>%
  #   paste(collapse = '; ')
  
  mod_labels <- gsub(pattern = ' - ', 
                     replacement = ' ',  
                     x = mod_labels)
  
  mod_num <- paste0('& (', 1:length(mod_labels), ') ')%>% 
    paste0(collapse = '') %>% trimws(which = 'right')
  patt_mod_num <- paste0('\\[-1.8ex] ', mod_num, '\\') 
  
  mod_labels <- paste0('\\multirow{2}{1cm}{', mod_labels, '}')
  repl_mod_labels <- paste0('\\[-1.8ex] ', rep('& ', length(mod_labels)) %>% paste0(collapse = ''), '\\') 
  
  
  # produce table and suppress console output
  sink('oo')
  # stargazer formatting
  tabtex <- stargazer::stargazer(ar1, 
                                 type = 'latex', 
                                 covariate.labels = c('\\nth{1} lag', NA), 
                                 dep.var.labels = name,
                                 no.space = T,
                                 ci = F,
                                 font.size = 'small',
                                 initial.zero = F,
                                 column.labels = mod_labels,
                                 # notes = notes,
                                 header = F, 
                                 model.numbers = T,
                                 column.sep.width = '1pt',
                                 omit.stat = c('ser', 'res.dev')) %>% 
    # regex replace to fix cell width
    gsub(pattern = strin,
         replacement = repla,
         x = .) %>% 
    gsub(pattern = patt_mod_num,
         replacement = repl_mod_labels, 
         x = ., 
         fixed = T) %>% 
    write(x = .,
          file = destination) %>% 
    capture.output()
  
  sink(NULL)
  unlink(x = 'oo')

}

chunk_rolling <- function(regs_list, regs_list_sum, ar_lags_sum, fore_horiz){
  # function to extract and manipulate regressions made on rolling windows with 
  # lstm predictions
  invisible(require(magrittr))
  
  ### For the ar1 rolling window results
  tidyout <- function(onereg, fore_horiz){
    # ancillary function to extract and label coefficients and dates
    
    # extract only the end date, as the window moves on
    end <- onereg %>% 
      model.frame() %>% 
      rownames() %>% 
      tail(fore_horiz + 1) %>% 
      head(1) %>% 
      as.Date()
    
    # extract coefficients, add enddate, turn in time format
    out <- broom::tidy(onereg) %>% 
      tibble::add_column(enddate = end) %>% 
      tibbletime::as_tbl_time(index = enddate)
    
    return(out)
  }
  
  # apply extractor to all regressions
  tidycoefs <- lapply(X = regs_list, FUN = tidyout, fore_horiz)
  
  # rowbind all different results
  out_ar1 <- tidycoefs %>% dplyr::bind_rows()
  
  ### For the ark rolling window results
  tidyout_sum <- function(ar1, regs_list_sum, ar_lags_sum, fore_horiz){
    # extract results from lm and 
    # add label. Does the heavy lifting.
    
    end <- ar1 %>% 
      model.frame() %>% 
      rownames() %>% 
      tail(fore_horiz + 1) %>% 
      head(1) %>% 
      as.Date()
    
    out <- tibble::tibble(ar_sum = regs_list_sum[,1],
                          ar_sum_se = regs_list_sum[,2],
                          enddate = end,
                          k_lags = ar_lags_sum) %>% 
      tibbletime::as_tbl_time(index = enddate)
    
    return(out)
  }
  
  
  # out_ark_sum <- furrr::future_pmap(.l = list(ar1 = regs_list, 
  #                                             ar_lags_sum = ar_lags_sum, 
  #                                             fore_horiz = fm_apply(fore_horiz, length(regs_list_sum)), 
  #                                             regs_list_sum = regs_list_sum), 
  #                                   .f = tidyout_sum) %>% 
  #   dplyr::bind_rows()
  
  out_ark_sum <- purrr::pmap(.l = list(ar1 = regs_list, 
                                       ar_lags_sum = ar_lags_sum, 
                                       fore_horiz = fm_apply(fore_horiz, length(regs_list_sum)), 
                                       regs_list_sum = regs_list_sum), 
                             .f = tidyout_sum) %>% 
    dplyr::bind_rows()
  
  out <- list(ar1_roll = out_ar1, 
              ark_sum_roll = out_ark_sum)
  
  return(out)  
}

plot_rollregs_lines <- function(chunk_regs_obj, graphs_dir. = graphs_dir, name){
  
  # line plot/save for the AR1 models on rolling window of data
  
  # part 1: bar plot for AR1 on each chunk
  
  # number of windows
  len <- chunk_regs_obj$ar1_roll %>% 
    select(enddate) %>% 
    unique() %>% 
    nrow()
  
  # setup title
  tt <- paste0(name %>% noms_tt(),
               ': rolling window')
  
  jj <- name %>% 
    noms() %>% 
    paste0('_ar1_rollwind.pdf')
  
  # make plot, filtering out intercept
  plt <- chunk_regs_obj$ar1_roll %>% 
    filter(term != '(Intercept)') %>% 
    ggplot(aes(x = enddate, y = estimate))+
    geom_line(size = .75) +
    geom_ribbon(aes(ymin = (estimate - std.error), ymax = (estimate + std.error)), alpha = .25)+
    geom_hline(yintercept = 0:1, size = .25, linetype = 2, colour = 'black') +
    ggtitle(tt) + theme_minimal() + ylab('AR(1) coefficient') + xlab('Sample end date') + 
    theme(plot.title = element_text(hjust = 0.5,
                                    size = rel(1.5)),
          axis.title = element_text(size = rel(1.5)),
          axis.text = element_text(size = rel(1.5))) + 
    geom_smooth(method = 'loess', se = F, colour = 'blue', formula = 'y~x', size = .5)
    
  
  
  
  ggsave(plot = plt, 
         filename = file.path(graphs_dir., 
                              jj),
         device = 'pdf', 
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  
  # Part 2: barplot for each chunk, sum of coefficients
  
  labely <- chunk_regs_obj$ark_sum_roll %>% 
    select(k_lags) %>% 
    unique() %>% 
    paste0('sum of AR(', ., ') coefficients')
  
  tt <- paste0(name %>% noms_tt(),
               ': rolling window')
  
  jj <- name %>% 
    noms() %>% 
    paste0('_ar3_rollwind.pdf')
  
  
  plt_sum <- chunk_regs_obj$ark_sum_roll %>% 
    ggplot(aes(x = enddate, y = ar_sum)) + 
    geom_line(size = .75) + theme_minimal() + ylab(labely) + xlab('Sample end date') + 
    geom_hline(yintercept = 0:1, size = .25, linetype = 2, colour = 'black') +
    geom_ribbon(aes(ymin = (ar_sum - ar_sum_se), ymax = (ar_sum + ar_sum_se)), alpha = .5) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = rel(1.5)), 
          axis.text = element_text(size = rel(1.5))) + ggtitle(tt) + 
    geom_smooth(method = 'loess', se = F, colour = 'blue', formula = 'y~x', size = .5)
  
  ggsave(plot = plt_sum, 
         filename = file.path(graphs_dir., 
                              jj),
         device = 'pdf', 
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  plt_list <- list(ar1 = plt, 
                   ark_sum = plt_sum)
  
  
  return(plt_list)
}

chunk_increm <- function(regs_list, regs_list_sum, ar_lags_sum, fore_horiz){
  # function to extract and manipulate regressions made on rolling windows with 
  # lstm predictions
  
  invisible(require(magrittr))
  
  ### For the ar1 rolling window results
  tidyout <- function(onereg, fore_horiz){
    # ancillary function to extract and label coefficients and dates
    
    # extract only the end date, as the window moves on
    end <- onereg %>% 
      model.frame() %>% 
      rownames() %>% 
      tail(fore_horiz + 1) %>% 
      head(1) %>% 
      as.Date()
    
    # extract coefficients, add enddate, turn in time format
    out <- broom::tidy(onereg) %>% 
      tibble::add_column(enddate = end) %>% 
      tibbletime::as_tbl_time(index = enddate)
    
    return(out)
  }
  
  # apply extractor to all regressions
  tidycoefs <- lapply(X = regs_list, FUN = tidyout, fore_horiz)
  
  # rowbind all different results
  out_ar1 <- tidycoefs %>% dplyr::bind_rows()
  
  ### For the ark rolling window results
  tidyout_sum <- function(ar1, regs_list_sum, ar_lags_sum, fore_horiz){
    # extract results from lm and 
    # add label. Does the heavy lifting.
    
    end <- ar1 %>% 
      model.frame() %>% 
      rownames() %>% 
      tail(fore_horiz + 1) %>% 
      head(1) %>% 
      as.Date()
    
    out <- tibble::tibble(ar_sum = regs_list_sum,
                          enddate = end,
                          k_lags = ar_lags_sum) %>% 
      tibbletime::as_tbl_time(index = enddate)
    
    return(out)
  }
  
  
  out_ark_sum <- furrr::future_pmap(.l = list(ar1 = regs_list, 
                                              ar_lags_sum = ar_lags_sum, 
                                              fore_horiz = fm_apply(fore_horiz, length(regs_list_sum)), 
                                              regs_list_sum = regs_list_sum), 
                                    .f = tidyout_sum) %>% 
    dplyr::bind_rows()
  
  out <- list(ar1_roll = out_ar1, 
              ark_sum_roll = out_ark_sum)
  
  return(out)  
}

plot_increm_lines <- function(chunk_regs_obj, graphs_dir. = graphs_dir, name){
  
  # line plot/save for the AR1 models on rolling window of data
  
  # part 1: bar plot for AR1 on each chunk
  
  # number of windows
  len <- chunk_regs_obj$ar1_roll %>% 
    select(enddate) %>% 
    unique() %>% 
    nrow()
  
  # setup title
  tt <- paste0(name %>% noms_tt(),
               ': increasing sample with forecasts')
  
  jj <- name %>% 
    noms() %>% 
    paste0('_ar1_increm.pdf')
  
  # make plot, filtering out intercept
  plt <- chunk_regs_obj$ar1_roll %>% 
    filter(term != '(Intercept)') %>% 
    ggplot(aes(x = enddate, y = estimate))+
    geom_line(size = .5) +
    geom_ribbon(aes(ymin = (estimate - std.error), ymax = (estimate + std.error)), alpha = .5)+
    geom_hline(yintercept = 0:1, size = .25, linetype = 2, colour = 'black') +
    ggtitle(tt) + theme_minimal() + ylab('AR(1) coefficient') + xlab('Sample end date') + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_smooth(method = 'loess', se = F, colour = 'blue', formula = 'y~x', size = .5)
  
  
  
  ggsave(plot = plt, 
         filename = file.path(graphs_dir., 
                              jj),
         device = 'pdf', 
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  
  # Part 2: barplot for each chunk, sum of coefficients
  
  labely <- chunk_regs_obj$ark_sum_roll %>% 
    select(k_lags) %>% 
    unique() %>% 
    paste0('sum of AR(', ., ') coefficients')
  
  tt <- paste0(name %>% noms_tt(),
               ': increasing sample with forecasts - ',
               labely)
  
  jj <- name %>% 
    noms() %>% 
    paste0('_ar3_increm.pdf')
  
  
  plt_sum <- chunk_regs_obj$ark_sum_roll %>% 
    ggplot(aes(x = enddate, y = ar_sum)) + 
    geom_line(size = .5) + theme_minimal() + ylab(labely) + xlab('Sample end date') + 
    geom_hline(yintercept = 0:1, size = .25, linetype = 2, colour = 'black') +
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(tt) + 
    geom_smooth(method = 'loess', se = F, colour = 'blue', formula = 'y~x', size = .5)
  
  ggsave(plot = plt_sum, 
         filename = file.path(graphs_dir., 
                              jj),
         device = 'pdf', 
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  plt_list <- list(ar1 = plt, 
                   ark_sum = plt_sum)
  
  
  return(plt_list)
}

chunk_increm_window <- function(ar1, ark, lags, name, graphs_dir. = graphs_dir){
  
  
  # AR1 part: take data and do a plot
  
  # labels and data
  
  ids <- 1:length(ar1) %>% 
    paste0('chunk_', .) %>% 
    fm_apply()
  
  ar1_tbl <- furrr::future_map(.x = ar1, 
                               .f = tbl2xts::xts_tbl)
  
  data_tbltime <- furrr::future_pmap_dfr(.l = list(.data = ar1_tbl,
                                                   chunk_id = ids),
                                         .f = tibble::add_column)
  
  tt <- name %>% noms_tt() %>%  paste0(., ' - AR(1) rolling window on increasing sample')
  jj <- name %>% noms %>% paste0(., '_ar1_rolling_within_increm.pdf')
  
  ar1_plt <- data_tbltime %>% 
    ggplot(aes(x = date, y = Var.1, group = chunk_id))+
    geom_ribbon(aes(ymin = Var.1 - .SD2, ymax = Var.1 + .SD2, group = chunk_id),
                colour = 'grey', alpha = .01, size = .1)+
    geom_line(colour = 'red', size = .5, alpha = .1)+
    geom_hline(yintercept = 0:1, size = .25, linetype = 2, colour = 'black') +
    theme_minimal() + ggtitle(tt) + 
    ylab('AR(1) coefficient') + xlab(element_blank())
  theme(plot.title = element_text(hjust = .5))
  
  ggsave(plot = ar1_plt,
         filename = file.path(graphs_dir., jj),
         device = 'pdf',
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  # AR_k part, just sum over basically
  
  ark_tbltime <- furrr::future_map(.x = ark,
                                   .f = tbl2xts::xts_tbl) %>% 
    furrr::future_map(.f = na.omit) 
  
  ark_tbltime <- furrr::future_pmap_dfr(.l = list(.data = ark_tbltime,
                                                  chunk_id = ids,
                                                  lag_k = fm_apply(lags, length(ark_tbltime))),
                                        .f = tibble::add_column) %>% 
    dplyr::rename(ark_sum=coredata.xts.)
  
  
  
  tt <- name %>% noms_tt() %>%  paste0(., ' - AR(',lags, ') sum rolling window on increasing sample')
  jj <- name %>% noms %>% paste0(., '_ark_sum_rolling_within_increm.pdf')
  
  ark_plt <- ark_tbltime %>% 
    ggplot(aes(x = date, y = ark_sum, group = chunk_id))+
    geom_line(size = .5, alpha = .25)+ theme_minimal() + 
    geom_hline(yintercept = 0:1, size = .25, linetype = 2, colour = 'black') +
    ylab('Sum of coefficients') + xlab(element_blank())+
    theme(plot.title = element_text(hjust = .5))+
    ggtitle(tt)
  
  ggsave(plot = ark_plt,
         filename = file.path(graphs_dir., jj),
         device = 'pdf',
         units = 'in', 
         width = 8, 
         height = 9*8/16)
  
  out <- list(ar1_plot = ar1_plt,
              ark_sum_plot = ark_plt)
  
  return(out)
  
}


##### IV - LSTM functions ######################################################

data_prepper <- function(data, train = 1, test = NULL){
  # function to rescale data to normal values
  # that are required for the LSTM model fit.
  # It stores first and second moments for later
  # scaling back to observed data. All outputs are
  # stored in a list, usually TS properties (index) are 
  # preserved.
  
  
  # train and test shall be expressed in percentage terms
  if (train>1){stop('\nTrain should be expressed in percentage not intergers.')}
  
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
  output$train <- list() # preallocate for compatibility
  output$train[['mean']] <- mean(data_train)
  output$train[['sd']] <- sd(data_train)
  output$train[['n_sample']] <- l_train
  # storage
  store_train <- ((data_train - output$train[['mean']])/output$train[['sd']])
  if (is.xts(store_train)){
    output$train[['train_norm']] <- store_train
  }else{
    output$train[['train_norm']] <- array(data = store_train,
                                          dim = c(l_train, 1))
  }
  
  
  if (train != 1){
    # condition on test subsample, whether present
    data_test <- data[(l_train+1):len, ]
    # rescale test
    output$test <- list() # preallocate for compatibility
    output$test[['mean']] <- mean(data_test)
    output$test[['sd']] <- sd(data_test)
    output$test[['n_sample']] <- l_test
    # storage
    store_test <- ((data_test - output$test[['mean']])/output$test[['sd']])
    if (is.xts(store_test)){
      output$test[['test_norm']] <- store_test
    }else{
      output$test[['test_norm']] <- array(data = store_test,
                                          dim = c(l_test, 1))
    }
    
  }
  
  return(output)
}

k_fullsample_1l <- function(data, 
                            n_steps, 
                            n_feat = 1, 
                            # model_compiled,
                            nodes = 50,
                            size_batch = 1,
                            epochs = 2000,
                            ES = F,
                            keepBest = F,
                            view_loss = F){
  
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
  
  
  
  invisible(require(magrittr))
  invisible(require(keras))
  invisible(require(tictoc))
  invisible(require(numbers))
  
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
  # not working as of now, needs impro
  if (size_batch == 'auto' && ES == F){
    size_batch <- batch_prime
  }
  
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
  
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  model_compiled <- keras_model_sequential()
  model_compiled %>%
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = F,
               stateful = T,
               kernel_regularizer = regularizer_l2(l = 0.01),
               batch_size = size_batch
    ) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  model_online <- keras_model_sequential() %>% 
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = F,
               stateful = T,
               batch_size = 1
    ) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  
  tictoc::tic('\n\nModel estimation')
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
                                             patience = floor(epochs*.2),
                                             min_delta = 1e-5, 
                                             restore_best_weights = keepBest)
                   ),
                   epochs = epochs,
                   validation_split = .1,
                   batch_size = size_batch,
                   view_metrics = view_loss)
  } else {
    # estimate with given number of epochs
    history <- fit(object = model_compiled, 
                   y = y_data_arr, 
                   x = x_data_arr, 
                   epochs = epochs, 
                   verbose = 2,
                   shuffle = F,
                   batch_size = size_batch,
                   view_metrics = view_loss)
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
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  return(out)
}

k_fullsample_2l <- function(data, 
                            n_steps, 
                            n_feat = 1, 
                            # model_compiled,
                            nodes = 50,
                            size_batch = 1,
                            epochs = 2000,
                            ES = F,
                            keepBest = F,
                            view_loss = F){
  
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
  
  #' *this version uses only TWO layers*
  
  
  
  invisible(require(magrittr))
  invisible(require(keras))
  invisible(require(tictoc))
  invisible(require(numbers))
  
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
  # not working as of now, needs impro
  if (size_batch == 'auto' && ES == F){
    size_batch <- batch_prime
  }
  
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
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  model_compiled <- keras_model_sequential()
  model_compiled %>%
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = T,
               stateful = T,
               batch_size = size_batch,
               kernel_regularizer = regularizer_l2(l = 0.01)
    ) %>% 
    layer_lstm(units = nodes,
               return_sequences = F,
               stateful = T,
               kernel_regularizer = regularizer_l2(l = 0.01)) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  model_online <- keras_model_sequential() %>% 
    layer_lstm(units = nodes,
               input_shape = c(n_steps, n_feat),
               return_sequences = T,
               stateful = T,
               batch_size = 1
    ) %>% 
    layer_lstm(units = nodes,
               return_sequences = F,
               stateful = T) %>% 
    layer_dense(units = 1) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  
  tictoc::tic('\n\nModel estimation')
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
                                             patience = floor(epochs*.2),
                                             min_delta = 1e-5, 
                                             restore_best_weights = keepBest)
                   ),
                   epochs = epochs,
                   validation_split = .1,
                   batch_size = size_batch,
                   view_metrics = view_loss)
  } else {
    # estimate with given number of epochs
    history <- fit(object = model_compiled, 
                   y = y_data_arr, 
                   x = x_data_arr, 
                   epochs = epochs, 
                   verbose = 2,
                   shuffle = F,
                   batch_size = size_batch,
                   view_metrics = view_loss)
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
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  return(out)
}

k_fullsample_nl <- function(data, 
                            n_steps, 
                            n_feat = 1, 
                            # model_compiled,
                            nodes = 50,
                            size_batch = 1,
                            epochs = 2000,
                            ES = F,
                            keepBest = F,
                            nodes_list,
                            options,
                            view_loss = F){
  
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
  
  #' *this version adapts for multiple layers*
  
  # costum, scope specific function to automate layering
  extra_layers <- function(nodes_list, options){
    
    # function to automate the layering of models 
    # in keras. It outputs a model ready to compile.
    
    # bunch of checks
    if (!is.list(nodes_list)) stop('\nNumber of nodes per layer must be declared in a list.')
    if (!is.list(options)) stop('\nOptions must be provided in a list.')
    if (is.null(options$input_shape)) stop('\nOptions must contain input shape ONLY FOR THE FIRST LAYER.')
    if (is.null(options$ret_sequences)) stop('\nOptions must declare whether or not to return sequences from one layer to the following as a list.\nThe last must NOT return sequences unless direct forecstas.')
    if (is.null(options$stateful)) stop('\nOptions must declare whether each layer is stateful (and in case set batch size accordingly).')
    
    
    invisible(require(keras))
    invisible(require(magrittr))
    
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
            return_sequences = options$ret_sequences[[l]],
            stateful = options$stateful[[l]]
          )
        
      }else{
        model %>% 
          layer_lstm(
            units = nodes_list[[l]],
            return_sequences = options$ret_sequences[[l]],
            stateful = options$stateful[[l]],
            kernel_regularizer = regularizer_l2(l = 0.01)
          )
      }
    }
    
    model %>% layer_dense(units = 1)
    
    return(model)
  }
  
  
  # packages required
  invisible(require(magrittr))
  invisible(require(keras))
  invisible(require(tictoc))
  invisible(require(numbers))
  
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
  # not working as of now, needs impro
  if (size_batch == 'auto' && ES == F){
    size_batch <- batch_prime
  }
  
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
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  # set batch_size in options list
  options$size_batch <- size_batch
  options$input_shape <- c(n_steps, n_feat)
  
  # setup and compile main model
  model_compiled <- extra_layers(nodes_list,
                                 options) %>% 
    compile(optimizer = 'adam', 
            loss = 'mse')
  
  # setup and compile online model
  options_online <- options
  options_online$size_batch <- 1
  model_online <- extra_layers(nodes_list,
                               options_online) %>% 
    compile(optimizer = 'adam',
            loss = 'mse')
  
  # train model
  tictoc::tic('\n\nModel estimation')
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
                                             patience = floor(epochs*.2),
                                             min_delta = 1e-5, 
                                             restore_best_weights = keepBest)
                   ),
                   epochs = epochs,
                   validation_split = .1,
                   batch_size = size_batch,
                   view_metrics = view_loss)
  } else {
    # estimate with given number of epochs
    history <- fit(object = model_compiled, 
                   y = y_data_arr, 
                   x = x_data_arr, 
                   epochs = epochs, 
                   verbose = 2,
                   shuffle = F,
                   batch_size = size_batch,
                   view_metrics = view_loss)
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

online_pred <- function(model_fitted, 
                        data_train, 
                        model_type = 'model_online', 
                        horizon = 4*10){
  
  # This function produces iterative, indirect predictions with 
  # a previously trained model. It copies weights and model structure
  # and resets batch to 1 so to make online predictions easily
  # and consistently. 'horizon' gives the nomber of indirect
  # predictions to produce. If data are TS then also dates are generated.
  
  # to use the 'online' model (same weights, batch set to 1) use the default option
  # 'model_online'; to use other versions of the model, specify it in the model_type
  # option, eg 'model_fitted'.
  
  invisible(require(keras))
  invisible(require(dplyr))
  
  
  # data_train is a list from data_prepper function!
  if (!is.list(data_train)) warning('\nProvide list from "data_prepper" function')
  
  # preallocate array with results
  pred <- array(NA, dim = c(horizon, 1))
  
  # retrieve input data and fitted model
  input <- data_train$train[['train_norm']]
  model_online <- model_fitted[[model_type]]
  
  if (is.xts(data_train$train[['train_norm']])){
    time_preds <- seq(from = end(input),
                      length.out = (horizon+1),
                      by = periodicity(input)$label)
    time_preds <- tail(time_preds, n = horizon)
    pred <- xts(pred, order.by = time_preds)
  }
  
  # retrieve input shape
  in_shape <- get_input_shape_at(object = model_fitted[[model_type]],
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
  
  va <- names(data_train$train[['train_norm']])
  if (is.xts(data_train$train[['train_norm']])){
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train', 
                                   date = time(data_train$train[['train_norm']])) %>% 
                        rename(value = all_of(va)) %>%
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
                        rename(value = all_of(va)), 
                      pred %>% 
                        as_tibble %>% 
                        add_column(label = 'forecast') %>% 
                        rename(value = V1))
  }
  
  # reconversion to values
  forecast$value <- forecast$value*data_train$train[['sd']] + data_train$train[['mean']]
  # names(forecast)[1] <- va
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  return(forecast)
}


multi_online <- function(model_fitted, 
                         data_train, 
                         model_type = 'model_online', 
                         horizon = 4*10,
                         reps = 100){
  
  # This function produces iterative, indirect predictions with 
  # a previously trained model. It copies weights and model structure
  # and resets batch to 1 so to make online predictions easily
  # and consistently. 'horizon' gives the nomber of indirect
  # predictions to produce. If data are TS then also dates are generated.
  
  # to use the 'online' model (same weights, batch set to 1) use the default option
  # 'model_online'; to use other versions of the model, specify it in the model_type
  # option, eg 'model_fitted'.
  
  invisible(require(keras))
  invisible(require(dplyr))
  
  
  # data_train is a list from data_prepper function!
  if (!is.list(data_train)) warning('\nProvide list from "data_prepper" function')
  
  # preallocate array with results
  pred <- array(NA, dim = c(horizon, reps))
  
  # retrieve input data and fitted model
  input <- data_train$train[['train_norm']]
  model_online <- model_fitted[[model_type]]
  
  if (is.xts(data_train$train[['train_norm']])){
    time_preds <- seq(from = end(input),
                      length.out = (horizon+1),
                      by = periodicity(input)$label)
    time_preds <- tail(time_preds, n = horizon)
    pred <- xts(pred, order.by = time_preds)
  }
  
  # retrieve input shape
  in_shape <- get_input_shape_at(object = model_fitted[[model_type]],
                                 node_index = 0)
  in_shape <- sapply(in_shape, FUN = c)
  
  pred_out <- NULL
  
  for (i in 1:reps){
    input <- data_train$train[['train_norm']]
    for (h in 1:horizon){
      input_lagged <- embed(input, in_shape[2])
      input_arr <- array(data = input_lagged,
                         dim = c(nrow(input_lagged), ncol(input_lagged),1))
      last_row <- nrow(input_arr)
      
      pred[h, i] <- predict(model_online,
                            x = array(data = input_arr[last_row, , ],
                                      dim = in_shape),
                            batch_size = 1)
      input <- rbind(input, pred[h,i])
    }
    pred_out <- bind_rows(pred_out, 
                          pred[, i] %>% 
                            as_tibble() %>% 
                            add_column(label = 'forecast',
                                       date = time_preds) %>% 
                            rename(value = V1) %>% 
                            mutate(value = as.numeric(value))
    )
  }
  
  va <- names(data_train$train[['train_norm']])
  if (is.xts(data_train$train[['train_norm']])){
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train', 
                                   date = time(data_train$train[['train_norm']])) %>% 
                        rename(value = all_of(va)) %>%
                        mutate(value = as.numeric(value)),
                      pred_out)
  } else {
    forecast <- rbind(data_train$train[['train_norm']] %>% 
                        as_tibble %>% 
                        add_column(label = 'train') %>% 
                        rename(value = all_of(va)), 
                      pred_out)
  }
  
  # reconversion to values
  forecast$value <- forecast$value*data_train$train[['sd']] + data_train$train[['mean']]
  # names(forecast)[1] <- va
  
  # wipe out mem from previous runs
  keras::k_clear_session()
  
  return(forecast)
  
}




##### Packages Loader #####

pkgs <- c(
  'ggplot2',
  'magrittr',
  'dplyr',
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
  'numbers',
  'rsample',
  'tbl2xts',
  'tibbletime',
  'xts',
  'tibble',
  'cowplot'
)
# fill pkgs with names of the packages to install

instant_pkgs(pkgs)

devtools::install_github('sboysel/fredr')
devtools::install_github('ceschi/urcabis')
# devtools::install_version("readxl", version = "1.0.0")
# library(urcabis) # for when the package will be duly updated (pull request)



#### housekeeping ####
rm(pkgs,lstm_dir)
