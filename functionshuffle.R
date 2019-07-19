# new file for inflation analysis
# all is just a funtion and then 
# one can apply it to the series DB

PI_test <- function(serie, coll){
	# serie: an individual TS with the inflation qtrly obs
	# coll: a strictly structured list that collects all
	#		results and stores them. It is pre-declared in
	#		the inflanalysis.R file

	##### Part I: printed output is saved to disk
	sink(file=paste0(file.path(graphs_dir, coll$names[[i]]), ' AR inflation results.txt'),
       append=F,
       split=T)

	##### Part II: unit root tests

	# 1 - run ADF test, compute last significant lag
	# llags is the upper limit of lags
	coll[['unitroot']][[i]] <- ur.df(na.omit(pi[,i]),
                                         # DF test, max lag to consider
                                        lags=llags,
                                         # lags selection criterion, min BIC
                                         selectlags='BIC')

	# tells if we will also use optimal lags in the
	# rest of the estimations: if yes then we store 
	# computed optilags
	if (flag___optilag==1) coll[['aropti']][[i]] <- coll[['unitroot']][[i]]@optilags

	# 2 - Unit root test results:
	# print results with p-values on the unit-root test
	# modified urca pkg
	if (coll[['unitroot']][[i]]@teststat>
      coll[['unitroot']][[i]]@cval[3]) {
              cat(paste0('For ', names(pi)[i], ' it is not possible to reject \nthe null hypothesis of unit root.\n\n'))
  	}else if ((coll[['unitroot']][[i]]@teststat<
            coll[['unitroot']][[i]]@cval[3]) &
            (coll[['unitroot']][[i]]@teststat>
             coll[['unitroot']][[i]]@cval[2])) {
              cat(paste0('For ', names(pi)[i], ' it is possible to reject \nthe null hypothesis of unit root at 90%.\n\n'))
  	}else if((coll[['unitroot']][[i]]@teststat<
            coll[['unitroot']][[i]]@cval[2]) &
           (coll[['unitroot']][[i]]@teststat>
            coll[['unitroot']][[i]]@cval[1])) {
              cat(paste0('For ', names(pi)[i], ' it is possible to reject \nthe null hypothesis of unit root at 95%.\n\n'))
  	}else if(coll[['unitroot']][[i]]@teststat<
           coll[['unitroot']][[i]]@cval[1]) {
              cat(paste0('For ', names(pi)[i], ' it is possible to reject \nthe null hypothesis of unit root at 99%.\n\n'))
  	}


  	# 3 - AR(5) regression: sum over the estimated coefficients
  	coll[['ark']]


}




