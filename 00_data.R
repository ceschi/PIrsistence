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

pi <- 
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
  add_column(agg = NA) %>% 
  select(-starts_with('real'), -series_id) %>% 
  mutate(q_date = paste0(year(date), "-", quarter(date)) %>% yq,
         y_date = year(date))


pi_mlybase <- pi %>% 
  group_by(id) %>% 
  mutate(
    mom = 1200*(log(value) - log(lag(value))),
    qoq = 1200*(log(value) - log(lag(value, n = 3))),
    yoy = 100*(log(value) - log(lag(value, n = 12)))
    )


pi_qlybase <- pi %>% 
  group_by(id, q_date) %>% 
  summarise(eop = last(value),
            avg = mean(value)) %>% 
  mutate(
    qoq_eop = 400*(log(eop) - log(lag(eop))),
    qoq_avg = 400*(log(avg) - log(lag(avg))),
    yoy_eop = 100*(log(eop) - log(lag(eop, 4))),
    yoy_avg = 100*(log(avg) - log(lag(avg, 4))),
    ) %>% 
  select(-eop, -avg) %>% 
  pivot_longer(cols = -(1:2), names_to = "vars", values_to = "vals") %>% 
  separate(col = vars, 
           into = c("freq", "agg"), 
           sep = '_', 
           remove = T) %>% 
  rename(date = q_date)

pi_ylybase <- pi %>% 
  group_by(id, y_date) %>% 
  summarise(
    eop = last(value),
    avg = mean(value)
  ) %>% 
  mutate(
    yoy_eop = 100*(log(eop) - log(lag(eop))),
    yoy_avg = 100*(log(avg) - log(lag(avg))),
  ) %>% 
  select(-eop, -avg) %>% 
  pivot_longer(cols = -(1:2), names_to = "vars", values_to = "vals") %>% 
  separate(col = vars, 
           into = c("freq", "agg"), 
           sep = '_', 
           remove = T) %>% 
  mutate(date = as_date(paste0(y_date, "-01-01")),
         y_date = NULL)

  

