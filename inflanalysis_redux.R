##### Code for Inflation analysis #####
##### On revised time series #####



##### Preamble #####
# Sets up stuff and collects
# data with old files.
source('directories.R')
source('functs.R')


ahead = 1
source('USdatacoll.R')


##### Parameters and al ####



# exogenous lag
k <- 5

# selector for coefficient
# AR(r) will be plotted
# MUST be <k
r <- 1

# select window width for
# rolling estimates, pick <80
# to get interesting results
wind <- 56 # as Pivetta and Reis, 14y. Alter: 5y


flag___optilag <- 1
flag___plot <- 0

##### Data subselection ######


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

# looping upper bound
n=length(names(pi))

# results collector
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
  aropti=list(),
  aroptilm=list(),
  plot_aropti=list(),
  rollm=list(),
  rollridges=list(),
  plot_rollm=list(),
  plot_ridges=list()
)

##### Unit root tests #####
if (flag___optilag==1){
  llags <- 18
} else if (flag___optilag==0){
  llags <- 8
}

for (i in 1:n){
  # saves output for every series
  sink(file=paste0(file.path(graphs_dir, inflation$names[[i]]), ' AR inflation results.txt'),
       append=F,
       split=T)
  
  ##### Unit root tests #####
  inflation[['unitroot']][[i]] <- ur.df(na.omit(pi[,i]),
                                         # DF test, max lag to consider
                                        lags=llags,
                                         # lags selection criterion, min BIC
                                         selectlags='BIC')
  # storing optimal lags -- forked urca pkg
  
  if (flag___optilag==1) inflation[['aropti']][[i]]<- inflation[['unitroot']][[i]]@optilags
  
  if (inflation[['unitroot']][[i]]@teststat>
      inflation[['unitroot']][[i]]@cval[3]) {
              cat(paste0('For ', names(pi)[i], ' it is not possible to reject \nthe null hypothesis of unit root.\n\n'))
  }else if ((inflation[['unitroot']][[i]]@teststat<
            inflation[['unitroot']][[i]]@cval[3]) &
            (inflation[['unitroot']][[i]]@teststat>
             inflation[['unitroot']][[i]]@cval[2])) {
              cat(paste0('For ', names(pi)[i], ' it is possible to reject \nthe null hypothesis of unit root at 90%.\n\n'))
  }else if((inflation[['unitroot']][[i]]@teststat<
            inflation[['unitroot']][[i]]@cval[2]) &
           (inflation[['unitroot']][[i]]@teststat>
            inflation[['unitroot']][[i]]@cval[1])) {
              cat(paste0('For ', names(pi)[i], ' it is possible to reject \nthe null hypothesis of unit root at 95%.\n\n'))
  }else if(inflation[['unitroot']][[i]]@teststat<
           inflation[['unitroot']][[i]]@cval[1]) {
              cat(paste0('For ', names(pi)[i], ' it is possible to reject \nthe null hypothesis of unit root at 99%.\n\n'))
  }
  
  
  
  ##### AR regression with fixed lags k=5 #####
  # k = 5 as standard practice but this might vary
  inflation[['ark']][[i]] <- lm(data= (pi[,i] %>% lagger(lag=k)),
                                formula=formula.maker(df=pi[,i] %>% lagger(lag=k),
                                                      y= pi[,i] %>% lagger(lag=k) %>% 
                                                        names(.) %>% first(),
                                                      intercept = F)
                                ) #%>% summary() %>% coef()
  cat('\n')
  print(paste0(inflation$names[[i]],',  ', k, ' exogenously defined lags'))
  print(summary(inflation[['ark']][[i]]))
  
  ##### AR regression with optimal lags #####
  if (flag___optilag==1){
    inflation[['aroptilm']][[i]] <- rolloop(df = pi[,i], window = wind, lags = inflation[['aropti']][[i]])

    inflation[['plot_aropti']][[i]] <- ggplot(data=inflation[['aroptilm']][[i]],
                                              aes(x=index(inflation[['aroptilm']][[i]]),
                                                  y=inflation[['aroptilm']][[i]][,r]))+
      # plot the above with line geom, in black
      geom_line(colour='black', size=1)+
      # adds upper confidence band in red
      geom_line(aes(y=inflation[['aroptilm']][[i]][,r]+inflation[['aroptilm']][[i]][,inflation[['aropti']][[i]]+r]),
                colour='red')+
      # adds lower confidence band in red
      geom_line(aes(y=inflation[['aroptilm']][[i]][,r]-inflation[['aroptilm']][[i]][,inflation[['aropti']][[i]]+r]),
                colour='red')+
      # adds unit root line
      geom_line(aes(y=1), colour='black', size=.8)+
      # plot makeup
      geom_smooth(method='loess', colour='blue')+scale_x_yearqtr(format='%Y Q%q', n=20)+theme_bw()+
      scale_y_continuous()+xlab(' ') + ylab(paste0('AR(',r,') coeff. estimates')) +
      ggtitle(paste0(inflation$names[[i]],' - ', inflation[['aropti']][[i]], ' optimal lags'))

    if  (flag___plot==0) plot(inflation[['plot_aropti']][[i]])

    # saves graphs in proper directory with names
    ggsave(paste0('AR(',r,') coeff. estimates ', inflation[['names']][[i]],
                  ' - ', inflation[['aropti']][[i]], ' optimal lags.pdf'),
           inflation[['plot_aropti']][[i]],
           device='pdf',
           graphs_dir,
           height=8, width=14.16, units='in')
  }
  ##### AR(k) rolling window regressions #####
  ## consider using map() or apply() function families
  ## to vectorise code and make it faster. This requires
  ## wrapping loop's actions in an ad-hoc function
  ## about window length: 58 is used by Pivetta&Reis as well as Fuhrer, trade off
  ## is about precision and extension of the resulting estimates
  inflation[['rollm']][[i]] <- rolloop(df = pi[,i], window = wind, lags = k)
  
  ##### AR(k) while keeping all coefs ####
  inflation[['rollridges']][[i]] <- persistence_ridges(tseries = pi[,i],
                                                       window = wind,
                                                       lags = inflation[['aropti']][[i]])
  
  
  ##### Plots registration #####
  # attaches data, select first column of estimates and time
  inflation[['plot_rollm']][[i]] <- ggplot(data=inflation[['rollm']][[i]],
                                           aes(x=index(inflation[['rollm']][[i]]),
                                               y=inflation[['rollm']][[i]][,r]))+
    # plot the above with line geom, in black
    geom_line(colour='black', size=1)+
    # adds upper confidence band in red
    geom_line(aes(y=inflation[['rollm']][[i]][,r]+inflation[['rollm']][[i]][,k+r]),
              colour='red')+
    # adds lower confidence band in red
    geom_line(aes(y=inflation[['rollm']][[i]][,r]-inflation[['rollm']][[i]][,k+r]),
              colour='red')+
    # adds unit root line
    geom_line(aes(y=1), colour='black', size=.8)+
    # plot makeup
    geom_smooth(method='loess', colour='blue')+scale_x_yearqtr(format='%Y Q%q', n=20)+theme_bw()+
    scale_y_continuous()+xlab(' ') + ylab(paste0('AR(',r,') coeff. estimates')) + 
    ggtitle(paste0(inflation$names[[i]], ' - ', k, ' exogenous lags'))
  
  inflation[['plot_ridges']][[i]] <- ggplot(data = inflation[['rollridges']][[i]])+
                                         geom_ridgeline_gradient(aes(x = term,
                                                            y = as.factor(last.date),
                                                            height = estimate,
                                                            group = as.factor(last.date),
                                                            fill = p.value),
                                                        min_height = -2) +
                                          scale_fill_viridis(name = 'P-values',option = "C", direction = -1) +
                                          ggtitle(paste0('Evolving persistence - ',
                                                         inflation$names[[i]],
                                                         ' ',
                                                         inflation[['aropti']][[i]],
                                                         ' end. lags')) +
                                          xlab('Lag order') + ylab(' ')
  ###########!!!!!! make up is missing!
  
  if  (flag___plot==0) plot(inflation[['plot_rollm']][[i]])
  
  # saves graphs in proper directory with names
  ggsave(paste0('AR(',r,') coeff. estimates ', inflation[['names']][[i]], '.pdf'),
         inflation[['plot_rollm']][[i]],
         device='pdf',
         graphs_dir,
         height=8, width=14.16, units='in')
  
  ggsave(paste0(i, ' - AR(',inflation[['aropti']][[i]],') acf ', inflation[['names']][[i]], '.pdf'),
         inflation[['plot_ridges']][[i]],
         device = 'pdf',
         graphs_dir,
         # extra height needed for full display
         height = 100,
         width = 14.16,
         units = 'in',
         limitsize = FALSE
         )
  
  cat('\n\n\n')
  # stopping printing
  sink()
}


lapply(inflation[['ark']], FUN = function(x) summary(x) %>% coefficients() %>% sum())


##### Housekeeping ####
rm(pi, n, i, llags)

if (flag___singular == 1) rm(r, k, wind)







