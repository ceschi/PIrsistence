library(R.matlab)
library(ggridges)

bsum <- R.matlab::readMat(file.choose()) %>% as.data.frame()

temp <- pi$rev_defl_pch
temp <- temp[!is.na(temp),]
names(bsum) <- temp[-(1:40), ] %>% index()

rm(temp)

bsum <- bsum %>% reshape2::melt() %>% 
  group_by(variable) %>% mutate(med = median(value),
                                mea = mean(value)) %>% 
  ungroup()

p <- bsum %>% 
  as_tibble() %>% 
  ggplot(aes(x = value, y = variable, group = variable)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.05, .5, 0.95),
                      scale = .95, rel_min_height = .05) + 
  coord_flip() + theme_ridges() + scale_y_yearqtr(n = 15) +
  geom_line(aes(x = med, y = value), size = 2) + 
  geom_line(aes(x = mea, y = value), size = 2)

format(object.size(bsum), 'Gb')


plot_draws <- function(main_path, var, name){
  require(R.matlab)
  require(ggplot2)
  require(ggridges)
  require(dplyr)
  require(lubridate)
  require(reshape2)
  
  
  full_path <- file.path(main_path, var)
  
  bsum <- R.matlab::readMat(file.path(full_path, 'BSUM.mat')) %>% 
          as.data.frame()
  
  temp <- pi[, grep(x = names(pi),
                    pattern = var,
                    value = F)
             ]
  
  temp <- temp[!is.na(temp),]
  
  names(bsum) <- temp[-(1:40), ] %>% index()
  
  
  bsum <- bsum %>% 
    reshape2::melt() %>% 
    mutate(variable = lubridate::yq(variable))
  
  bsum_collapsed <- bsum %>% group_by(variable) %>% 
    mutate(median = median(value), 
           mean = mean(value), 
           q5 = quantile(value, .05), 
           q95 = quantile(value, .95)
           ) %>% 
    select(-value) %>% 
    distinct() 

  
  p <- bsum %>% 
    ggplot(aes(x = value, y = variable, group = variable)) +
    stat_density_ridges(quantile_lines = TRUE, 
                        quantiles = c(0.05, .5, 0.95),
                        scale = 5, 
                        rel_min_height = .05) + 
    coord_flip() + theme_ridges() + 
    ylab(' ') + xlab('Per period draws') + ggtitle(name)+
    theme(axis.text.x = element_text(angle = 45)) +
    geom_vline(xintercept = 1, colour = 'black', size = .75) +
    geom_vline(xintercept = 0, colour = 'black', size = .75)
  
  p2 <- bsum_collapsed %>% 
    ggplot() +
    geom_line(aes(x = variable,
                  y = median),
              size = 1.5, 
              alpha = .5,
              colour = 'black') + 
    geom_line(aes(x = variable,
                  y = q5),
              size = .75,
              colour = 'red') +
    geom_line(aes(x = variable,
                  y = q95),
              size = .75,
              colour = 'red') +
    geom_hline(yintercept = 0, colour = 'black', size = .75)+
    geom_hline(yintercept = 1, colour = 'black', size = .75)+
    theme_minimal() + ggtitle(name) + labs(colour = 'Qtls') + 
    theme(axis.text.x = element_text(angle = 45),
          legend.position = 'none') +
    ylab(' ') + xlab(' ')
  
  
  # saves the plots in given path
  ggsave(paste0(var, '_distro.pdf'),
         plot = p,
         device='pdf',
         path = graphs_dir,
         height=8/.8, width=14.16/.8, units='in')
  
  # saves the plots in given path
  ggsave(paste0(var, '_lines.pdf'),
         plot = p2,
         device='pdf',
         path = graphs_dir,
         height=8/.8, width=14.16/.8, units='in')
  
  
  plots <- list(p, p2)
  
  return(plots)
  invisible(gc())
}

var <- names(pi) %>% grep(pattern = '_pch', x = ., value = T)
main_path <- 'C:/Users/emanu/Desktop'
names <- c(
  'Revised CPI pch',
  'Revised CPI, no FE pch',
  'Revised GDP deflator pch',
  'Revised PCE pch',
  'Revised PCE, no FE pch'
)

mega_plots <- list()
pb <- txtProgressBar(min = 0, 
                     max = length(var),
                     style = 3)
for (i in 1:length(var)){
  
  mega_plots[[i]] <- plot_draws(main_path, var[i], names[i])
  setTxtProgressBar(pb, i)
}