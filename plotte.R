###### Plotting script ########################################################

# general plot with raw data
inflation$plots[['plot_ts']] <- ggplot(pi["1945/2020"], aes(x = index(pi["1945/2020"]))) + 
  geom_line(aes(y = rev_cpi_pch, colour = 'CPI'), alpha = .75) + 
  geom_line(aes(y = rev_cpi_fe_pch, colour = 'CPI FE'), alpha = .75) + 
  geom_line(aes(y = rev_pce_pch, colour = 'PCE'), alpha = .75) + 
  geom_line(aes(y = rev_pce_fe_pch, colour = 'PCE FE'), alpha = .75)+
  geom_line(aes(y = rev_defl_pch, colour = 'Deflt.'), alpha = .75) +
  theme_minimal() + theme(legend.title = element_blank()) +
  ggtitle('Inflation series') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 45),
    plot.title = element_text(hjust = 0.5, size = rel(3.5)))

ggsave(filename = fvars_dir, 'ts_plot.pdf'),
       plot = inflation$plots[['plot_ts']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

##### CPI & PCE headline vs core ##############################################
inflation$plots[['cpis']] <- pi_long %>% 
  filter(index == 'rev_cpi_pch' | index == 'rev_cpi_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index))+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('CPI: core vs headline') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 45),
    plot.title = element_text(hjust = 0.5, size = rel(3.5)))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))


ggsave(filename = file.path(vars_dir, 'cpis_plot.pdf'),
       plot = inflation$plots[['cpis']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['pces']] <- pi_long %>% 
  filter(index == 'rev_pce_pch' | index == 'rev_pce_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index))+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('PCE: core vs headline') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 45),
    plot.title = element_text(hjust = 0.5, size = rel(3.5)))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))

ggsave(filename = file.path(vars_dir, 'pces_plot.pdf'),
       plot = inflation$plots[['pces']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['cpis_zoom']] <- pi_long %>% filter(date>as.Date('1989-12-31')) %>% 
  filter(index == 'rev_cpi_pch' | index == 'rev_cpi_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index))+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('CPI: core vs headline since 1990') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 45),
    plot.title = element_text(hjust = 0.5, size = rel(3.5)))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))


ggsave(filename = file.path(vars_dir, 'cpis_zoomed_plot.pdf'),
       plot = inflation$plots[['cpis_zoom']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['pces_zoom']] <- pi_long %>% filter(date>as.Date('1989-12-31')) %>% 
  filter(index == 'rev_pce_pch' | index == 'rev_pce_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index))+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('PCE: core vs headline since 1990') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 45),
    plot.title = element_text(hjust = 0.5, size = rel(3.5)))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))

ggsave(filename = file.path(vars_dir, 'pces_zoomed_plot.pdf'),
       plot = inflation$plots[['pces_zoom']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['headcores']] <- cowplot::plot_grid( 
  # full cpi
  inflation$plots[['cpis']],
  # full pce
  inflation$plots[['pces']],
  # zoomed cpi
  inflation$plots[['cpis_zoom']],
  # zoomed pce
  inflation$plots[['pces_zoom']],
  # arrangement
  ncol = 2, nrow = 2
)

ggsave(filename = file.path(vars_dir, 'headcores_plot.pdf'),
       plot = inflation$plots[['headcores']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)


oil_p <- fredr_series_observations(series_id = 'DCOILWTICO',
                                    frequency = 'q',
                                    aggregation_method = 'eop') %>% 
  select(-series_id) %>% 
  tibbletime::as_tbl_time(date) %>% 
  rename(oil=value) %>% mutate(oil_p=c(NA, 400*diff(log(oil))))


inflation$plots[['oil']] <- oil_p %>% ggplot()+
  geom_line(aes(x = date, y = oil)) + theme_minimal()+
  ggtitle('Crude oil spot price: WTI level')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = rel(3.5))) + 
  xlab(' ') + ylab(' ')

inflation$plots[['oil_p']] <- oil_p %>% ggplot()+
  geom_line(aes(x = date, y = oil_p)) + theme_minimal()+
  ggtitle('Crude oil spot price: WTI price change')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = rel(3.5))) + 
  xlab(' ') + ylab(' ')

inflation$plots[['oil_grid']] <- cowplot::plot_grid(inflation$plots[['oil']],
                                                    inflation$plots[['oil_p']],
                                                    nrow = 2)

ggsave(filename = file.path(vars_dir, 'wti_pi_plot.pdf'),
       plot = inflation$plots[['oil_grid']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

#### all commodities

comms <- fredr_series_observations(series_id = 'PALLFNFINDEXM',
                                   frequency = 'q',
                                   aggregation_method = 'eop') %>% 
  select(-series_id) %>% 
  tibbletime::as_tbl_time(date) %>% 
  mutate(value_p=c(NA, 400*diff(log(value))))

inflation$plots[['comm']] <- comms %>% ggplot()+
  geom_line(aes(x = date, y = value)) + theme_minimal()+
  ggtitle('Commodities GPI: level')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = rel(3.5))) + 
  xlab(' ') + ylab(' ')

inflation$plots[['comm_p']] <- comms %>% ggplot()+
  geom_line(aes(x = date, y = value_p)) + theme_minimal()+
  ggtitle('Commodities GPI: price change')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = rel(3.5))) + 
  xlab(' ') + ylab(' ')

inflation$plots[['comm_grid']] <- cowplot::plot_grid(inflation$plots[['comm']],
                                                    inflation$plots[['comm_p']],
                                                    nrow = 2)

ggsave(filename = file.path(vars_dir, 'comm_pi_plot.pdf'),
       plot = inflation$plots[['comm_grid']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)



##### histograms ###############################################################
inflation$plots$histo <- pi_long %>%
  filter(str_detect(index, 'pch')) %>% 
  ggplot(data = ., aes(x = rate, group = index)) + geom_density(fill = 'grey') + 
  facet_wrap(~index) + theme_minimal() + ggtitle('Kernel densities') + 
  theme(plot.title = element_text(hjust = .5))


##### other windows width for robustness #######################################
source('other_win.R')

# housekeeping
rm(oil_p, comms)


