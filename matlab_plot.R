plot_draws <- function(main_path, var, name){
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
  
  
  bsum <- bsum %>% 
    reshape2::melt() %>% 
    dplyr::mutate(variable = lubridate::ymd(variable)) %>% 
    dplyr::group_by(variable) %>% 
    dplyr::mutate(Median = median(value)) %>% 
    dplyr::ungroup()
  
  bsum_collapsed <- bsum %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(median = median(value), 
           mean = mean(value), 
           q5 = quantile(value, .05), 
           q95 = quantile(value, .95)
           ) %>% 
    dplyr::select(-value) %>% 
    dplyr::distinct() 
  
  # test chunck
  # bsum <- bsum[sample(1:nrow(bsum), 10000),]
  
  p <- bsum %>% 
    ggplot(aes(x = value, y = variable, group = variable, fill = Median)) +
    stat_density_ridges(quantile_lines = TRUE, 
                        quantiles = c(0.05, .5, 0.95),
                        scale = 5, 
                        rel_min_height = .05) + 
    coord_flip() + ylab(' ') + 
    xlab('Per period draws') + ggtitle(name) +
    theme(axis.text.x = element_text(angle = 0)) +
    geom_vline(xintercept = 1, colour = 'black', size = .75) +
    geom_vline(xintercept = 0, colour = 'black', size = .75) +
    scale_fill_viridis(option = 'C') + theme_ridges() + 
    theme(legend.position = 'bottom', 
    	legend.justification = 'center',
    	plot.title = element_text(hjust = 0.5)) +
    guides(colour=guide_legend(nrow = 1, byrow = T))
  
  
  p2 <- bsum_collapsed %>% 
    ggplot() +
    geom_line(aes(x = variable,
                  y = median),
              size = 1.5, 
              alpha = 1,
              colour = 'black') + 
    geom_ribbon(aes(ymin = q5,
                    ymax = q95,
                    y = median,
                    x = variable),
              colour = 'red', 
              alpha = .2,
              size = .5) +
    geom_hline(yintercept = 0, colour = 'black', size = .75)+
    geom_hline(yintercept = 1, colour = 'black', size = .75)+
    theme_ridges() + ggtitle(name) + labs(colour = 'Qtls') + 
    theme(axis.text.x = element_text(angle = 0),
          legend.position = 'none',
    	  plot.title = element_text(hjust = 0.5)) +
    ylab(' ') + xlab(' ')
  
  
  # saves the plots in given path
  ggsave(paste0(var, '_distro.pdf'),
         plot = p,
         device='pdf',
         path = graphs_dir,
         width = 8,
         height = 8*9/16, 
         units='in')
  
  # saves the plots in given path
  ggsave(paste0(var, '_lines.pdf'),
         plot = p2,
         device='pdf',
         path = graphs_dir,
         width = 8,
         height = 8*9/16, 
         units='in')
  
  
  plots <- list(p, p2)
  
  return(plots)
  invisible(gc())
}

var <- names(pi) %>% grep(pattern = '_pch', x = ., value = T)
main_path <- 'C:/Users/emanu/Desktop'
main_path <- 'D:/emanu/OneDrive/R/pirsistence/matlab/'
names <- c(
  'Revised CPI pch',
  'Revised CPI, no FE pch',
  'Revised GDP deflator pch',
  'Revised PCE pch',
  'Revised PCE, no FE pch'
)

mega_plots <- list()
names <- names %>% gsub('Revised ', '', .) %>% gsub(', no FE pch', ' core', .) %>% gsub(' pch', ' headline', .)


pb <- txtProgressBar(min = 0, 
                     max = length(var),
                     style = 3)

for (i in 1:length(var)){
  
  mega_plots[[i]] <- plot_draws(main_path, var[i], names[i])
  setTxtProgressBar(pb, i)
}

rm(pb)
