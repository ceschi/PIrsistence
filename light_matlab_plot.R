plot_light_draws <- function(main_path, var, name, lines_dir, ttil){
  invisible(require(R.matlab))
  invisible(require(ggplot2))
  invisible(require(ggridges))
  invisible(require(dplyr))
  invisible(require(lubridate))
  invisible(require(reshape2))
  
  
  full_path <- file.path(main_path, var)
  
  bsum <- R.matlab::readMat(file.path(full_path, 'BSUM.mat')) %>% 
          as.data.frame() 
  
  temp <- pi[, grep(x = names(pi),
                    pattern = var,
                    value = F)
             ]
  
  temp <- temp[!is.na(temp),]
  
  col_k <- bsum %>% ncol()
  row_k <- temp %>% nrow()
  k <- row_k - col_k
  
  names(bsum) <- temp[-(1:k), ] %>% index()

  
  bsum_collapsed <- bsum %>% 
    reshape2::melt() %>% 
    mutate(variable = lubridate::ymd(variable)) %>% 
    dplyr::group_by(variable) %>% 
    dplyr::mutate(median = median(value), 
           mean = mean(value), 
           q5 = quantile(value, .05), 
           q95 = quantile(value, .95)
           ) %>% 
    dplyr::select(-value) %>% 
    dplyr::distinct() 

  
  tt <- paste0(name, ttil)
  
  p2 <- bsum_collapsed %>% 
    ggplot() +
    geom_line(aes(x = variable,
                  y = mean),
              size = 1, 
              alpha = .1,
              colour = 'black') + 
    geom_ribbon(aes(ymin = q5,
                    ymax = q95,
                    y = median,
                    x = variable),
              colour = 'red', 
              alpha = .1,
              size = .5) +
    geom_hline(yintercept = 0:1, colour = 'black', size = .75)+
    theme_ridges() + ggtitle(tt) + labs(colour = 'Qtls') + 
    theme(axis.text.x = element_text(angle = 0),
          legend.position = 'none',
      	  plot.title = element_text(hjust = 0.5), 
      	  axis.title = element_blank()) +
    ylab(' ') + xlab(' ')
  
  
  # saves the plots in given path
  ggsave(paste0(var, '_lines.pdf'),
         plot = p2,
         device='pdf',
         path = lines_dir,
         width = 8,
         height = 8*9/16, 
         units='in')
  
  return(p2)
  invisible(gc())
}

names <- c(
  'Revised CPI pch',
  'Revised CPI, no FE pch',
  'Revised GDP deflator pch',
  'Revised PCE pch',
  'Revised PCE, no FE pch'
)


# multipath ---------------------------------------------------------------

# get folders
par_dir <- list.dirs(path = '/home/e.franceschi/Bureau/brute_force_shortdata',
                     recursive = F)

# store numbers
draws <- strsplit(x = par_dir, split = '/draws') %>% 
  sapply(X = ., FUN = `[[`, 2)



dir_lines <- sapply(X = file.path(lines_dir, draws),
                    FUN = dir.create) %>% 
  names()



so_many_plots <- list()
interator <- 1

for (dd in 1:length(draws)){
  for (v in 1:length(var)){
    
    
    templot <- plot_light_draws(main_path = par_dir[dd], 
                                var = var[v],
                                name = names[v],
                                lines_dir = dir_lines[dd],
                                ttil = draws[dd]
                                )
    
    so_many_plots[[iterator]] <- templot
    iterator <- 1 + iterator
  }
}

