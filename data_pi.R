#### data fetching #############################################################

# required pkgs
invisible(require(magrittr))
invisible(require(xts))
invisible(require(tbl2xts))
invisible(require(devtools))
if (!('fredr' %in% installed.packages()[,1])){
  devtools::install_github('sboysel/fredr', force = TRUE)
}
invisible(require(fredr))



# pick ahead to set how many quarters ahead
# to consider for SPF forecasts:
# -1 for previous quarter estimates
# 0 for nowcast
# 1 for one quarter ahead -- default
# 2 for two quarters ahead
# 3 for three quarters ahead
# 4 for one year ahead
if (!base::exists(x = "ahead", envir = .GlobalEnv)){stop('Provide how many quarters ahead in an "ahead" variable.')}

# key for the FRED API
fredr_set_key('5d4b5f1e6667727ee4ea90affbad1e6a')


##### Actual data fetching
# Section to get historical, revised inflation TS

rev_hist_pch <- merge(

  # Consumer Price Index for All Urban Consumers: All Items
  rev_pci = fredr_series_observations(series_id='CPIAUCSL',
                                      frequency='q',
                                      aggregation_method='eop',
                                      units='cca') %>% tbl_xts(),

  # Consumer Price Index for All Urban Consumers: All Items Less Food and Energy
  rev_pci_fe  = fredr_series_observations(series_id='CPILFESL',
                                          frequency='q',
                                          aggregation_method='eop',
                                          units='cca') %>% tbl_xts(),

  # Gross Domestic Product: Implicit Price Deflator
  rev_defl = fredr_series_observations(series_id='GDPDEF',
                                       frequency='q',
                                       aggregation_method='eop',
                                       units='cca') %>% tbl_xts(),

  # Personal Consumption Expenditures including Food and Energy
  rev_pce  = fredr_series_observations(series_id='PCEPI',
                                       frequency='q',
                                       aggregation_method='eop',
                                       units='cca') %>% tbl_xts(),

  # Personal Consumption Expenditures Excluding Food and Energy
  rev_pce_fe  = fredr_series_observations(series_id='PCEPILFE',
                                          frequency='q',
                                          aggregation_method='eop',
                                          units='cca') %>% tbl_xts()
)

# renames variables
names(rev_hist_pch) <-  c('rev_cpi_pch', 'rev_cpi_fe_pch', 'rev_defl_pch',
                          'rev_pce_pch', 'rev_pce_fe_pch')

##### Year on Year instead of from previous quarter
rev_hist_yoy <- merge(

  # Consumer Price Index for All Urban Consumers: All Items
  rev_pci = fredr_series_observations(series_id='CPIAUCSL',
                                      frequency='q',
                                      aggregation_method='eop',
                                      units='pc1') %>% tbl_xts(),

  # Consumer Price Index for All Urban Consumers: All Items Less Food and Energy
  rev_pci_fe  = fredr_series_observations(series_id='CPILFESL',
                                          frequency='q',
                                          aggregation_method='eop',
                                          units='pc1') %>% tbl_xts(),

  # Gross Domestic Product: Implicit Price Deflator
  rev_defl = fredr_series_observations(series_id='GDPDEF',
                                       frequency='q',
                                       aggregation_method='eop',
                                       units='pc1') %>% tbl_xts(),

  # Personal Consumption Expenditures including Food and Energy
  rev_pce  = fredr_series_observations(series_id='PCEPI',
                                       frequency='q',
                                       aggregation_method='eop',
                                       units='pc1') %>% tbl_xts(),

  # Personal Consumption Expenditures Excluding Food and Energy
  rev_pce_fe  = fredr_series_observations(series_id='PCEPILFE',
                                          frequency='q',
                                          aggregation_method='eop',
                                          units='pc1') %>% tbl_xts()
)

# renames variables
names(rev_hist_yoy) <-  c('rev_cpi_yoy', 'rev_cpi_fe_yoy', 'rev_defl_yoy',
                          'rev_pce_yoy', 'rev_pce_fe_yoy')


#### merge series to single object ############################################
# some obs are dropped because they are not reliable
pi <- merge(
  # QoQ
  rev_hist_pch$rev_cpi_pch,
  rev_hist_pch$rev_cpi_fe_pch['1966-01-01/'],
  rev_hist_pch$rev_pce_pch,
  rev_hist_pch$rev_pce_fe_pch,
  rev_hist_pch$rev_defl_pch,

  #YoY
  rev_hist_yoy$rev_cpi_yoy,
  rev_hist_yoy$rev_cpi_fe_yoy,
  rev_hist_yoy$rev_pce_yoy,
  rev_hist_yoy$rev_pce_fe_yoy,
  rev_hist_yoy$rev_defl_yoy
)

# reshape data in long format
pi_long <- pi %>%
  as_tibble %>%
  add_column(date = time(pi)) %>%
  pivot_longer(
    cols = -date,
    names_to = "index",
    values_to = "rate"
  )

# write out to disk the series
write.zoo(x=pi,
          file=file.path(data_dir, 'PI_data.csv'),
          sep=';',
          row.names=F,
          index.name='date')

write.csv(x = pi_long,
          file = file.path(data_dir, 'PI_long.csv'))

##### housekeeping
rm(rev_hist_pch, rev_hist_yoy)
