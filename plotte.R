###### Plotting script ########################################################

# set n back to all variables
n <- length(names(pi))

# general plot with raw data
inflation$plots[['plot_ts']] <- ggplot(pi["1945/2020"], aes(x = index(pi["1945/2020"]))) + 
  geom_line(aes(y = rev_cpi_pch, colour = 'CPI'), alpha = .75) + 
  geom_line(aes(y = rev_cpi_fe_pch, colour = 'CPI Core'), alpha = .75) + 
  geom_line(aes(y = rev_pce_pch, colour = 'PCE'), alpha = .75) + 
  geom_line(aes(y = rev_pce_fe_pch, colour = 'PCE Core'), alpha = .75)+
  geom_line(aes(y = rev_defl_pch, colour = 'GDP Deflator'), alpha = .75) +
  theme_minimal() + theme(legend.title = element_blank()) +
  xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_viridis_d(option = 'A', end = .8)
  

ggsave(filename = file.path(vars_dir, 'ts_plot.pdf'),
       plot = inflation$plots[['plot_ts']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

##### CPI & PCE headline vs core ##############################################
inflation$plots[['cpis']] <- pi_long %>% 
  filter(index == 'rev_cpi_pch' | index == 'rev_cpi_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index), alpha = .5)+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('CPI: core vs headline') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'none', 
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))


ggsave(filename = file.path(vars_dir, 'cpis_plot.pdf'),
       plot = inflation$plots[['cpis']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['pces']] <- pi_long %>% 
  filter(index == 'rev_pce_pch' | index == 'rev_pce_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index), alpha = .5)+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('PCE: core vs headline') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'none', 
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))

ggsave(filename = file.path(vars_dir, 'pces_plot.pdf'),
       plot = inflation$plots[['pces']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['cpis_zoom']] <- pi_long %>% filter(date>as.Date('1989-12-31')) %>% 
  filter(index == 'rev_cpi_pch' | index == 'rev_cpi_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index), alpha = .5)+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('CPI: core vs headline since 1990') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5))+ 
  scale_colour_manual(labels = c("core", "headline"), values = c("darkblue", "red"))


ggsave(filename = file.path(vars_dir, 'cpis_zoomed_plot.pdf'),
       plot = inflation$plots[['cpis_zoom']], 
       device = 'pdf',
       width = 8,
       units = 'in', 
       height = 9*8/16)

inflation$plots[['pces_zoom']] <- pi_long %>% filter(date>as.Date('1989-12-31')) %>% 
  filter(index == 'rev_pce_pch' | index == 'rev_pce_fe_pch') %>% 
  ggplot() + geom_line(aes(x = date, y = rate, colour = index), alpha = .5)+
  theme_minimal() + theme(legend.title = element_blank()) + 
  ggtitle('PCE: core vs headline since 1990') + xlab(' ') + ylab(' ') +
  guides(colour=guide_legend(nrow = 1, byrow = T)) + 
  theme(legend.position = 'bottom', 
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5))+ 
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
  rename(oil=value) %>% 
  mutate(oil_p=c(NA, 400*diff(log(oil))))

# imf food indx 
#' *WIP*
download.file('https://api.db.nomics.world/v22/series/IMF/COMMP/M.W0.PFANDB_Index.xlsx',
              destfile = file.path(data_dir, 'imffood.xlsx'),
              method = 'curl')

food_p <- readxl::read_excel(path = file.path(data_dir, 'imffood.xlsx'),
                             col_names = c("date", "food"),
                             skip = 1) %>%
  mutate(date = lubridate::ym(date),
         date_q = lubridate::quarter(date, with_year = T)) %>%
  group_by(date_q) %>% 
  mutate(food_q = last(food)) %>% 
  select(date_q, food_q) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(food_q_g = 4*log(food_q)/lag(log(food_q)))

inflation$plots[['food']] <- food_p %>% 
  ggplot(aes(x = date_q, y = food_q)) +
  geom_line() + theme_minimal() + 
  ggtitle('Food Commodities Index - IMF') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = .5))+ 
  xlab(' ') + ylab(' ')


inflation$plots[['food_p']] <- food_p %>% 
  ggplot(aes(x = date_q, y = food_q_g)) +
  geom_line() + theme_minimal() + 
  ggtitle('Food Commodities Index Quarterly Change - IMF') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = .5))+ 
  xlab(' ') + ylab(' ')

inflation$plots[['food_grid']] <- 
  cowplot::plot_grid(inflation$plots[["food"]],
                     inflation$plots[["food_p"]],
                     nrow = 2)


inflation$plots[['oil']] <- oil_p %>% ggplot()+
  geom_line(aes(x = date, y = oil)) + theme_minimal()+
  ggtitle('Crude oil spot price: WTI level')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5)) + 
  xlab(' ') + ylab(' ')

inflation$plots[['oil_p']] <- oil_p %>% ggplot()+
  geom_line(aes(x = date, y = oil_p)) + theme_minimal()+
  ggtitle('Crude oil spot price: WTI price change')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5)) + 
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
    plot.title = element_text(hjust = 0.5)) + 
  xlab(' ') + ylab(' ')

inflation$plots[['comm_p']] <- comms %>% ggplot()+
  geom_line(aes(x = date, y = value_p)) + theme_minimal()+
  ggtitle('Commodities GPI: price change')+
  theme(legend.position = 'none',
    plot.title = element_text(hjust = 0.5)) + 
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
  filter(stringr::str_detect(index, 'pch')) %>% 
  group_by(index) %>% 
  mutate(mean = mean(rate),
         median = median(rate)) %>% 
  ungroup() %>% 
  ggplot(data = ., aes(x = rate, group = index)) + geom_density(fill = 'grey') + 
  geom_vline(xintercept = 0) + 
  facet_wrap(~index) + theme_minimal() + ggtitle('Kernel densities') + 
  theme(plot.title = element_text(hjust = .5))

##### Plots for frequentist part ###############################################

# AR(1) rolling
inflation$plots[['plot_rollm']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                                   names = inflation[['names']],
                                                   path = sapply(rep(ar1_dir, n), list)),
                                         .f = plot_roller
                                         )

# AR(1) trend rolling
inflation$plots[['plot_rollm_trend']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                                   names = inflation[['names']],
                                                   path = sapply(rep(ar1_dir, n), list),
                                                   .slot = fm_apply(2, n)),
                                         .f = plot_roller
                                         )


# AR(k*) plots
inflation$plots[['plot_aropti']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                                    names = inflation[['names']],
                                                    laags = inflation[['aropti']],
                                                    path = sapply(rep(ark_dir, n), list)),
                                          .f = plot_autoregsum
                                          )

# AR(k*) plots trend
inflation$plots[['plot_aropti_trend']] <- future_pmap(.l = list(df = inflation[['rollark']],
                                                          names = inflation[['names']],
                                                          laags = inflation[['aropti']],
                                                          path = sapply(rep(ark_dir, n), list),
                                                          .slot = fm_apply(2, n)),
                                                .f = plot_autoregsum
                                                )

# plotting ridges
inflation$plots[["plot_ridges"]] <- future_pmap(.l = list(df = inflation[['aroptiridges']],
                                                    nam = inflation[['names']],
                                                    laags = inflation[['aropti']],
                                                    path = sapply(rep(acf_dir, n), list)),
                                          .f = plot_ridges
                                          )


##### other windows width for robustness #######################################
source('other_win.R')

# housekeeping
rm(oil_p, comms)


