# new file for inflation analysis
# start anew with different structure
# each step is run in vectorised way
# on each series.



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
llags <- 8
if (flag___optilag == 1) llags <- 18

# flag on plotting
# 1 for plotting results
# 0 for skipping
flag___plot <- 0

# flag on the MarkovS states,
# default 2
flag___ms <- 2


# directories, functions and data

source('directories.R')
source('functs.R')
source('USdatacoll.R')


# subselect data
pi <- merge(
  # db_US$cpit,
  # db_US$coret,
  # db_US$deflt,
  # db_US$deflt1,
  db_US$rev_cpi_pch,
  db_US$rev_cpi_fe_pch,
  db_US$rev_defl_pch,
  db_US$rev_pce_pch,
  db_US$rev_pce_fe_pch,
  db_US$rev_cpi_yoy,
  db_US$rev_cpi_fe_yoy,
  db_US$rev_defl_yoy,
  db_US$rev_pce_yoy,
  db_US$rev_pce_fe_yoy
)

n=length(names(pi))

# this preallocated list will
# collect all results
inflation <- list(
  names=list(# 'CPI nowcast',
    # 'PCE nowcast',
    # 'GDP deflator nowcast',
    # 'GDP deflator forecast',
    'Revised CPI pch',
    'Revised CPI, no FE pch',
    'Revised GDP deflator pch',
    'Revised PCE pch',
    'Revised PCE, no FE pch',
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
inflation[['aroptilm']] <- pmap(.l = list(data = sapply(pi, list),
                                          lags = inflation[['aropti']],
                                          interc = sapply(rep(intercep, n), list)
                                          ),
                                .f = auto.reg.sum)


# rolling window estimate
# width is preselected,
# optilag is computed in step I
inflation[['aroptirollm']] <- pmap(.l = list(df = sapply(pi, list),
                                              window = sapply(rep(wind, n), list),
                                              lags = inflation[['aropti']],
                                              interc = sapply(rep(intercep, n), list)
                                              ),
                                    .f = rolloop.sum)


# Markov Switching model on the k* lags
# on the whole sample

inflation[['aropti_ms']] <- pmap(.l = list(df = sapply(pi,FUN = function(x) list(as.data.frame(x))),
                                           lags = inflation[['aropti']],
                                           states = flag___ms
                                           ),
                                 .f = function(df, lags, states){
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
                                   
                                 })




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
                                    ggsave(paste0(names, ' - AR(1) coeff estimates.pdf'),
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
                                     
                                     ggsave(paste0(names, ' - AR(',laags,') coeff estimates sum.pdf'),
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


# plotting msm

inflation[['plot_aropti_ms']] <- pmap(.l = list(ms_model = inflation[['aropti_ms']],
                                                nam = inflation[['names']],
                                                laags = inflation[['aropti']],
                                                path = sapply(rep(graphs_dir, n), list)
                                                ),
                                                
                                      .f = function(ms_model, nam, laags, path){
                                        
                                        # setting device size
                                        # mar sets margings
                                        # cex.main scales title to 70%
                                        par(mar = c(1,1,2.85,1), cex.main = .70)
                                        
                                        # store actual plot
                                        plot_out <- plotProb(ms_model, which = 2)
                                        
                                        
                                        
                                        # print to device
                                        print(plot_out)
                                        
                                        # fix title
                                        title(paste0(flag___ms, '-state MS regimes for ', nam, ' with ', laags, ' lags'), line = 2.3)
                                        # copy dev output to file (pdf)
                                        dev.copy(pdf, height=8/1.5, width=14.6/1.5,
                                                 file.path(path,
                                                           paste0(nam, ' ', flag___ms, '-state MSM.pdf')))
                                        
                                        # shut down device, comment for keeping the plot
                                        invisible(dev.off())
                                        
                                        # output
                                        return(plot_out)
                                      }
                                        )


