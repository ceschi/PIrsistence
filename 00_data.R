### Data collection


# key for the FRED API
fredr_set_key('5d4b5f1e6667727ee4ea90affbad1e6a')


# 1 - Quarterly series EOP --------------------------------------------------------

pi_q_eop <-
  bind_rows(
    rev_pci = fredr_series_observations(series_id='CPIAUCSL',
                                        frequency='q',
                                        aggregation_method='eop') %>%
      add_column(id = "cpi"),

    rev_pci_core = fredr_series_observations(series_id='CPILFESL',
                                        frequency='q',
                                        aggregation_method='eop') %>%
      add_column(id = "cpi_core"),

    rev_pce = fredr_series_observations(series_id='PCEPI',
                                        frequency='q',
                                        aggregation_method='eop') %>%
      add_column(id = "pce"),

    rev_pce_core = fredr_series_observations(series_id='PCEPILFE',
                                        frequency='q',
                                        aggregation_method='eop') %>%
      add_column(id = "pce_core"),

    rev_defl = fredr_series_observations(series_id='GDPDEF',
                                        frequency='q',
                                        aggregation_method='eop') %>%
      add_column(id = "defl")
  ) %>%
  add_column(freq = "q", agg = "eop") %>%
  select(-starts_with('real'), -series_id)


### Series are monthly by default

pi_m <- 
  bind_rows(
    rev_pci = fredr_series_observations(series_id='CPIAUCSL') %>% 
      add_column(id = "cpi"),
    
    rev_pci_core = fredr_series_observations(series_id='CPILFESL') %>% 
      add_column(id = "cpi_core"),
    
    rev_pce = fredr_series_observations(series_id='PCEPI') %>% 
      add_column(id = "pce"),
    
    rev_pce_core = fredr_series_observations(series_id='PCEPILFE') %>% 
      add_column(id = "pce_core"),
    
    # GDP defl is only quarterly
    # rev_defl = fredr_series_observations(series_id='GDPDEF') %>% 
    #   add_column(id = "defl")
  ) %>% 
  add_column(freq = "m", agg = NA) %>% 
  select(-starts_with('real'), -series_id)



pi_m %>% 
  group_by(id, q_date = paste0(year(date), "-", quarter(date)) %>% yq) %>% 
  summarise(
    lin_eop = last(value),
    lin_avg = mean(value)
  )
