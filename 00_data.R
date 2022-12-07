### Data collection


# key for the FRED API
fredr_set_key('5d4b5f1e6667727ee4ea90affbad1e6a')


# 1 - Fetch and treat series ---------------------------------------------------

##### Get and transform GDP deflator data
# get the data and parse year - it's only quarterly
defl <- fredr_series_observations(series_id='GDPDEF') %>% 
  add_column(id = "defl") %>% 
  select(date, value, id) %>% 
  mutate(y_date = year(date))

# from qtly basis compute yoy and qoq
defl_qlybase <- defl %>% 
  select(-y_date) %>% 
  mutate(qoq = 400*(log(value) - log(lag(value))),
         yoy = 100*(log(value) - log(lag(value, n = 4))),
         base = 'quarter',
         value = NULL,
         agg = "eop",
         ) %>% 
  pivot_longer(c(qoq, yoy), 
               names_to = "freq", 
               values_to = "vals")

# aggregate up to yearly two ways: eop and avg
# for each compute the yoy
defl_ylybase <- defl %>% 
  group_by(id, y_date) %>% 
  summarise(eop = last(value),
            avg = mean(value)) %>% 
  mutate(yoy_eop = 100*(log(eop) - log(lag(eop))),
         yoy_avg = 100*(log(avg) - log(lag(avg))),
         base = "year",
         date = as_date(paste0(y_date, "-01-01")),
         y_date = NULL,
         eop = NULL,
         avg = NULL,
         ) %>% 
  pivot_longer(cols = starts_with("yoy_"), 
               names_to = 'freq', 
               values_to = "vals") %>% 
  separate(freq, 
           into = c("freq", "agg"), 
           sep = "_")
  
##### Get and transform CPI, PCE data
### Series are monthly by default

# first off all data at monthly freq
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
  add_column(agg = "eop") %>% 
  select(-starts_with('real'), -series_id) %>% 
  # add quarter dates and years
  mutate(q_date = paste0(year(date), "-", quarter(date)) %>% yq,
         y_date = year(date))

#' from a monthly frequency compute
#' mom, qoq, yoy,
#' for both avg and eop aggregations
pi_mlybase <- pi %>% 
  select(-q_date, -y_date) %>% 
  group_by(id) %>% 
  mutate(
    mom = 1200*(log(value) - log(lag(value))),
    qoq = 400*(log(value) - log(lag(value, n = 3))),
    yoy = 100*(log(value) - log(lag(value, n = 12)))
    ) %>% 
  pivot_longer(c(mom, qoq, yoy), 
               names_to = "freq", 
               values_to = "vals") %>% 
  mutate(base = "month", 
         value = NULL)

#' aggregate up the frequency to
#' qtrly via eop and avg
#' then compute qoq, yoy
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
  rename(date = q_date) %>% 
  add_column(base = "quarter")

#' last, up to yearly
#' and compute yoy
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
         base = "year",
         y_date = NULL)

##### Join all together

allflation <- 
  bind_rows(
    defl_qlybase,
    defl_ylybase,
    pi_mlybase,
    pi_qlybase,
    pi_ylybase
    )

## keep only those that make sense
inflation <- allflation %>% 
  filter(base == "quarter" &
           agg == "eop" &
           date >= "1960-01-01")


inflation_monthly <- allflation %>% 
  filter(base == 'month' &
           id %in% c("pce", "pce_core"))

# N - housekeeping --------------------------------------------------------

rm(defl, defl_qlybase, defl_ylybase,
   pi, pi_mlybase, pi_qlybase, pi_ylybase)


# N2 - Add other series at industry level ---------------------------------


