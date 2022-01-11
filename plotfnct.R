

plot_layers <- function(.name, .l1, .l2, .destpath){

  all <- bind_rows(.l1 %>%
                     mutate(label = case_when(label == 'train' ~ 'Data',
                                              label == 'forecast' ~ 'One Layer')),
                   .l2 %>%
                     mutate(label = case_when(label == 'train' ~ 'Data',
                                              label == 'forecast' ~ 'Two Layers'))
                   ) %>%
    distinct()

  tt <- paste0(noms_tt(.name), ' LSTMs Forecasts')
  ff <- paste0(noms(.name), '_multiforecasts.pdf')
  vls <- c('One Layer' = 'firebrick1',
           'Two Layers' = 'navyblue',
           'Data' = 'black')

  plt <- all %>%
    ggplot(aes(x = date, y = value, group = label))+
    geom_line(aes(colour = label)) + theme_bw() +
    xlab(element_blank()) + ylab(element_blank())+
    ggtitle(tt) + guides(colour = guide_legend(nrow = 1)) +
    scale_colour_manual(values = vls) +
    theme(axis.text = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1.5)),
          title = element_text(size = rel(1.5)),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'bottom',
          legend.title = element_blank(),
          panel.background = element_blank())

  ggsave(filename = file.path(.destpath, ff),
         plot = plt,
         device = 'pdf',
         width = 8,
         height = 9*8/16)

}


pmap(.l = list(.l1 = inflation$lstm[['online_pred_1l']],
               .l2 = inflation$lstm[['online_pred_2l']],
               .name = inflation$names[1:n],
               .destpath = file.path(graphs_dir, 'plots_lstm')),
     .f = plot_layers)