##### Specifically designed functions ####

# A file to gather all home made functions with relative descriptions


instant_pkgs <- function(pkgs) { 
  ## Function loading or installing packages in
  ## current R instance.
  ## Developed by Jaime M. Montana Doncel - V1

  
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
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



subfilter <- function(df){
  # function to convert a df with multiple observations per unit
  # of time in a df with one observation per unit of time,
  # namely the last one among those previously present
  
  
  indice <- as.character(unique(df$date))
  len <- length(indice)
  outp <- matrix(NA, ncol=ncol(df), nrow=len)
  outp <- data.frame(outp)
  names(outp) <- names(df)
  for (i in 1:len){
    supp <- indice[i]
    ram <- subset(df, date==supp)
    outp[i,] <- ram[nrow(ram),]
    outp[i,1] <- indice[i]
  }
  return(outp)
}


subfilter.mean <- function(df){
  # function to convert a df with multiple observations per unit
  # of time in a df with one observation per unit of time,
  # namely the mean of those previously present
  
  
  indice <- as.character(unique(df$date))
  len <- length(indice)
  outp <- matrix(NA, ncol=ncol(df), nrow=len)
  outp <- data.frame(outp)
  names(outp) <- names(df)
  for (i in 1:len){
    supp <- indice[i]
    ram <- subset(df, date==supp)
    outp[i,] <- c(0, as.numeric(apply(ram[,-1], 2, mean)))
  }
  outp[,1] <- indice
  return(outp)
}


trendev<-function(mat){
  # for multiple observation in particular shape, this function
  # estimates a quadratic trend on the available series and consider
  # the deviation from the trend in the last observation. This deviation
  # is put into another time series. The purpose of this function is to
  # extract real time output gap from Philadelphia dataset.
  
  
  matdat<-mat[,2:ncol(mat)]
  temp<-1:nrow(mat)
  temp2<-temp^2
  regr<-function(x){
    dta<-data.frame(x, temp, temp2)
    names(dta)<-c('x', 'temp', 'temp2')
    model<-lm(x~temp+temp2, data=dta)
    GAPS<-(model$residuals/(x-model$residuals))
    gaps<-as.matrix(na.omit(GAPS))
    gap<-gaps[nrow(gaps)]
    return(gap)
  }
  outcome<-apply(matdat, 2, regr)
  outcome<-as.matrix(outcome)
  return(outcome*100)
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


spf_funct <-  function(filnam, typs, ahead=1) {
  # this function imports the files, reformats,
  # renames, saves in raw format and produces
  # aggregate statistics in XTS format
  
  # read in xlsx files and reshape w\ spread
  # this block selects one quarter ahead forecasts
  # but adjusting 'ahead' parameter below one can
  # extract other values
  
  # ad-hoc function inconsistent w/ external use
  # typs is one of CPI, CORECPI, PCE, COREPCE
  
  
  # 'ahead' allows to select the horizon of 
  # forecasts one wishes to extract:
  # -1 for previous quarter estimates
  # 0 for nowcast
  # 1 for one quarter ahead -- default
  # 2 for two quarters ahead
  # 3 for three quarters ahead
  # 4 for one year ahead
  
  typ=tolower(typs)
  
  colu=c(rep('numeric',3),  # picks year, quarter, ID
         rep('skip', 2+ahead),	 # skips industry
         'numeric',				 # moving target picking 'ahead' horizon
         rep('skip', 7-ahead)	 # skips the rest
  )
  
  df=read_excel(file.path(temp_dir,filnam), 
                na='#N/A', col_types=colu) %>%
    spread(ID, paste0(typs,ahead+2)) %>% 
    ts(start=c(1968, 4), frequency=4) %>%
    as.xts()
  
  pst=paste0(typ,'_')
  if (ahead==-1){
    pst=paste0(pst,'b1')
  } 	else {
    pst=paste0(pst,'h') %>% paste0(ahead)
  }
  
  names(df)=c('year', 'quarter', paste(pst, (1:(ncol(df)-2)), sep='_'))
  
  df$year <- df$quarter <- NULL
  
  # saving in txt csv format the raw data
  write.zoo(df, file.path(data_dir, paste(paste0('SPF_IND_',pst),'txt', sep='.')), sep=';', row.names=F, index.name='time')
  
  
  iqr <- apply(df, 1, IQR, na.rm=TRUE) %>% ts(start=c(1968, 4), frequency=4) %>% as.xts()
  stand<-apply(df, 1, var, na.rm=T) %>% sqrt()%>% ts(start=c(1968, 4), frequency=4) %>% as.xts()
  mean<-apply(df, 1, mean, na.rm=T)%>% ts(start=c(1968, 4), frequency=4) %>% as.xts()
  mean[is.nan(mean)] <- NA
  
  lab <- paste0('spf_', pst)
  
  df_stat=merge(iqr, stand, mean)
  names(df_stat)=paste(lab, c('iqr', 'sd', 'mean'), sep='_')
  
  
  return(df_stat)
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
  
  coef_sum <- output %>% filter(term != '(Intercept)') %>% select(estimate) %>%  sum()
  
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
    geom_smooth(method='loess', colour='blue')+scale_x_yearqtr(format='%Y Q%q')+theme_minimal()+
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
    geom_smooth(method='loess', colour='blue')+scale_x_yearqtr(format='%Y Q%q')+theme_minimal()+
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



##### Packages Loader #####

pkgs <- c('vars', 'glue', 'lazyeval',
          'quantreg', 'tidyverse', 'devtools',
          'tseries', 'dynlm', 'stargazer',
          'dyn', 'strucchange', 'xts', 'httr',
          'MASS', 'car', 'rvest', 'viridis',
          'mFilter', 'fredr','ggridges', 'MSwM',
          'readr', 'quantmod','broom', 'fredr',
          'devtools', 'lubridate', 'ggridges',
          'readxl', 'tbl2xts', 'tictoc', 'gmm',
          'future', 'furrr', 'broom')
# fill pkgs with names of the packages to install

instant_pkgs(pkgs)

# devtools::install_github('sboysel/fredr')
devtools::install_github('ceschi/urcabis')
# devtools::install_version("readxl", version = "1.0.0")
# library(urcabis) # for when the package will be duly updated (pull request)



#### housekeeping ####
rm(pkgs)
