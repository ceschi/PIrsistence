

for (i in 1:n){
  # saves output for every series
  sink(file=paste0(file.path(graphs_dir, inflation$names[[i]]), ' AR inflation results.txt'),
       append=F,
       split=T)
  
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
  
  cat('\n')
  print(paste0(inflation$names[[i]],',  ', k, ' exogenously defined lags'))
  print(summary(inflation[['ark']][[i]]))
  
  cat('\n\n\n')
  # stopping printing
  sink()
}


##### Housekeeping ####
rm(pi, n, i, llags)

if (flag___singular == 1) rm(r, k, wind)







