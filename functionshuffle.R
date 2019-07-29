# new file for inflation analysis
# start anew with different structure
# each step is run in vectorised way
# on each series.



#### 0 - Setup, data, flags ####

# horizon for now/forecast
ahead <- 1

# exogenous lags
k <- 1


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
  aroptiridges = list(),
  
  
  # 4+
  plot_rollm = list(),
  plot_aropti = list(),
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


# ridges plot material

inflation[['aroptiridges']] <- pmap(.l = list(tseries = sapply(pi, list),
                                              window = sapply(rep(wind, n), list),
                                              lags = inflation[['aropti']]),
                                    .f = persistence_ridges)







#### IV - Plots et al ##########################################################

# AR(1) rolling

inflation[['plot_rollm']] <- pmap(.l = list(df = inflation[['rollark']],
                                            names = inflation[['names']],
                                            path = sapply(rep(graphs_dir, n), list)),
                                  .f = function(df, names, path){
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
                                      # plot makeup
                                      geom_smooth(method='loess', colour='blue')+scale_x_yearqtr(format='%Y Q%q', n=20)+theme_bw()+
                                      scale_y_continuous()+xlab(' ') + ylab(paste0('AR(1) coeff. estimates')) + 
                                      ggtitle(paste0(names, ' - 1 exogenous lag'))
                                    
                                    
                                    
                                    # saves the plots in given path
                                    ggsave(paste0(names, ' - AR(1) coeff. estimates.pdf'),
                                           plot = po,
                                           device='pdf',
                                           path = path,
                                           height=8, width=14.16, units='in')
                                    
                                    
                                    return(po)
                                    
                                  }
                                  )


# AR(k*) plots

inflation[['plot_aropti']] <- pmap(.l = list(df = inflation[['rollark']],
                                             names = inflation[['names']],
                                             laags = inflation[['aropti']],
                                             path = sapply(rep(graphs_dir, n), list)),
                                   .f = function(df, names, path, laags){
                                     po <- ggplot(data=df,
                                                  aes(x=index(df),
                                                      y=df[,1]))+
                                       # plot the above with line geom, in black
                                       geom_line(colour='black', size=1)+
                                       # adds unit root line
                                       geom_line(aes(y=1), colour='black', size=.8)+
                                       # plot makeup
                                       geom_smooth(method='loess', colour='blue')+scale_x_yearqtr(format='%Y Q%q', n=20)+theme_bw()+
                                       scale_y_continuous()+xlab(' ') + ylab(paste0('AR(',laags,') coeff. estimates sum')) + 
                                       ggtitle(paste0(names, ' - ', laags, ' optimal lags: sum of coefficients'))
                                     
                                     
                                     # save plot
                                     
                                     ggsave(paste0(names, ' - AR(',laags,') coeff. estimates sum.pdf'),
                                            plot = po,
                                            device='pdf',
                                            path = path,
                                            height=8, width=14.16, units='in')
                                     
                                     return(po)
                                     
                                   }
                                   )

# plotting ridges

inflation[["plot_ridges"]] <- pmap(.l = list(df = inflation[['aroptiridges']],
                                             nam = inflation[['names']],
                                             laags = inflation[['aropti']],
                                             path = sapply(rep(graphs_dir, n), list)),
                                    
                                   .f = function(df, nam, laags, path){
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
                                        xlab('Lag order') + ylab(' ')
                                      
                                      
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
                                   )



