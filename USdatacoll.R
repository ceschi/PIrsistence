# This code collects and scraps core data for the analysis
# and puts it in xts format


#### Scraping US data ####


#### FEDERAL INTEREST RATE ####

fredr_set_key('5d4b5f1e6667727ee4ea90affbad1e6a')
# key for the FRED API

ffr <- fredr_series_observations(series_id='FEDFUNDS', frequency='m') %>% tbl_xts()
ffr <- as.xts(aggregate(ffr, as.yearqtr(as.yearmon(time(ffr))), last))
# aggregates up to quarters picking quarter's last month value

ffrb <- lag(ffr)

ffrate <- merge(ffr, ffrb)

#### INFLATION FORECASTS & REVISED ####

# downloads the big xlsx Greenbook file in
# a specifically created folder
download.file('https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/greenbook-data/documentation/gbweb_row_format.xls?la=en',
              file.path(temp_dir,'Greenbook_allvar_row.xlsx'), mode='wb',
              extra='--no-check-certificate',
              quiet = T)

# reads the single interesting sheets and
# imports them in df format

classi <- c('text', rep('numeric', 14), 'text')

cpi_greenbook <- read_excel(file.path(temp_dir,'Greenbook_allvar_row.xlsx'), 
                            sheet='gPCPI', col_types=classi, na='#N/D')
core_greenbook <- read_excel(file.path(temp_dir,'Greenbook_allvar_row.xlsx'),
                             sheet='gPCPIX', col_types=classi, na='#N/D')
deflator_greenbook <- read_excel(file.path(temp_dir,'Greenbook_allvar_row.xlsx'),
                                 sheet='gPGDP', col_types=classi, na='#N/D')

# # replace NAs
# cpi_greenbook[cpi_greenbook=='NaN'] <- NA
# core_greenbook[core_greenbook=='NaN'] <- NA
# deflator_greenbook[deflator_greenbook=='NaN'] <- NA

# drop useless columns
cpi_greenbook <- cpi_greenbook[,-c(2:5, 16, 15)]
core_greenbook <- core_greenbook[,-c(2:5, 16, 15)]
deflator_greenbook <- deflator_greenbook[,-c(2:5, 16, 15)]

# name columns
names(cpi_greenbook) <- c('date', 'cpit', paste(rep('cpit', 7), 1:8, sep=''))
names(core_greenbook) <- c('date', 'coret', paste(rep('coret', 7), 1:8, sep=''))
names(deflator_greenbook) <- c('date', 'deflt', paste(rep('deflt', 7), 1:8, sep=''))


# drop useless observations via subfilter fncts
cpi <- subfilter(cpi_greenbook)
cpi.mean <- subfilter.mean(cpi_greenbook)

core <- subfilter(core_greenbook)
core.mean <- subfilter.mean(core_greenbook)

defl <- subfilter(deflator_greenbook)
defl.mean <- subfilter.mean(deflator_greenbook)


# time series conversion

cpi <- as.xts(ts(cpi, start=c(1967, 1), frequency = 4))
cpi$date <- NULL

core <- as.xts(ts(core, start=c(1967, 1), frequency = 4))
core$date <- NULL

defl <- as.xts(ts(defl, start=c(1967, 1), frequency = 4))
defl$date <- NULL

rates <- merge(ffrate, cpi, core, defl)

# same, but for mean series

cpi.mean <- as.xts(ts(cpi.mean, start=c(1967, 1), frequency = 4))
cpi.mean$date <- NULL

core.mean <- as.xts(ts(core.mean, start=c(1967, 1), frequency = 4))
core.mean$date <- NULL

defl.mean <- as.xts(ts(defl.mean, start=c(1967, 1), frequency = 4))
defl.mean$date <- NULL

rates.mean <- merge(ffrate, cpi.mean, core.mean, defl.mean)

# Section to get historical, revised inflation TS

rev_hist_pch <- merge(

          # Consumer Price Index for All Urban Consumers: All Items 
          rev_pci = 4*fredr_series_observations(series_id='CPIAUCSL', 
                                          frequency='q', 
                                        # aggregation_method='eop',
                                        units='pch') %>% tbl_xts(), 
          
          # Consumer Price Index for All Urban Consumers: All Items Less Food and Energy
          rev_pci_fe  = 4*fredr_series_observations(series_id='CPILFESL', 
                                            frequency='q', 
                                            # aggregation_method='eop', 
                                            units='pch') %>% tbl_xts(),
          
          # Gross Domestic Product: Implicit Price Deflator
          rev_defl = 4*fredr_series_observations(series_id='GDPDEF', 
                                          frequency='q', 
                                         # aggregation_method='eop', 
                                         units='pch') %>% tbl_xts(),
          
          # Personal Consumption Expenditures including Food and Energy
          rev_pce  = 4*fredr_series_observations(series_id='PCEPI', 
                                             frequency='q', 
                                         # aggregation_method='eop', 
                                         units='pch') %>% tbl_xts(),
          
          # Personal Consumption Expenditures Excluding Food and Energy
          rev_pce_fe  = 4*fredr_series_observations(series_id='PCEPILFE', 
                                            frequency='q', 
                                            # aggregation_method='eop', 
                                            units='pch') %>% tbl_xts()
) 

# renames variables
names(rev_hist_pch) <-  c('rev_cpi_pch', 'rev_cpi_fe_pch', 'rev_defl_pch',
                      'rev_pce_pch', 'rev_pce_fe_pch')

##### Year on Year instead of from previous quarter
rev_hist_yoy <- merge(
  
  # Consumer Price Index for All Urban Consumers: All Items 
  rev_pci = fredr_series_observations(series_id='CPIAUCSL', 
                                        frequency='q', 
                                        # aggregation_method='eop',
                                        units='pc1') %>% tbl_xts(), 
  
  # Consumer Price Index for All Urban Consumers: All Items Less Food and Energy
  rev_pci_fe  = fredr_series_observations(series_id='CPILFESL', 
                                            frequency='q', 
                                            # aggregation_method='eop', 
                                            units='pc1') %>% tbl_xts(),
  
  # Gross Domestic Product: Implicit Price Deflator
  rev_defl = fredr_series_observations(series_id='GDPDEF', 
                                         frequency='q', 
                                         # aggregation_method='eop', 
                                         units='pc1') %>% tbl_xts(),
  
  # Personal Consumption Expenditures including Food and Energy
  rev_pce  = fredr_series_observations(series_id='PCEPI', 
                                         frequency='q', 
                                         # aggregation_method='eop', 
                                         units='pc1') %>% tbl_xts(),
  
  # Personal Consumption Expenditures Excluding Food and Energy
  rev_pce_fe  = fredr_series_observations(series_id='PCEPILFE', 
                                            frequency='q', 
                                            # aggregation_method='eop', 
                                            units='pc1') %>% tbl_xts()
)

# renames variables
names(rev_hist_yoy) <-  c('rev_cpi_yoy', 'rev_cpi_fe_yoy', 'rev_defl_yoy',
                          'rev_pce_yoy', 'rev_pce_fe_yoy')



## UNEMPLOYMENT METRICS ####

claims <- fredr_series_observations(series_id='ICSA', frequency='q', aggregation_method='sum') %>% tbl_xts()
# initial claims, number

natural_unemp_short <- fredr_series_observations(series_id='NROUST', frequency='q') %>% tbl_xts()
# natural employment on the short run

natural_unemp_long <- fredr_series_observations(series_id='NROU', frequency='q') %>% tbl_xts()
# longer term natural unemployment rate

current_unemp <- fredr_series_observations(series_id='UNRATE', frequency='q') %>% tbl_xts()
# current unemployment rate

tot_emp <- fredr_series_observations(series_id='PAYEMS', frequency='q') %>% tbl_xts() %>% `*`(.,1000)
# total employed, thousands

## Unemployment manipulation

short_long_diff <- natural_unemp_short - natural_unemp_long
layoffs <- 100*claims/tot_emp
employment_fluct <- current_unemp - natural_unemp_long

## merging

unemployment <- merge(layoffs, employment_fluct, current_unemp, short_long_diff)
names(unemployment) <- c('layoffs', 'employment_fluct', 'unempl_rate', 'unemp_sh_lng')


#### OUTPUT GAPS ####
# expost gap

capacity <- fredr_series_observations(series_id='GDPPOT', frequency='q') %>% tbl_xts()
# real installed capacity, 2009 chained dollars

actual <- fredr_series_observations(series_id='GDPC1', frequency='q') %>% tbl_xts()
# actual gdp

gap_expost <- (actual-capacity)*100/capacity 

# real time gap

download.file('https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/real-time-data/data-files/files/xlsx/routputqvqd.xlsx?la=en',
              file.path(temp_dir,'PhilFed_realtime_realgdp.xlsx'), mode='wb',
              extra='--no-check-certificate',
              quiet = T)

gdp_waves <- read_excel(file.path(temp_dir,'PhilFed_realtime_realgdp.xlsx'), 
                        sheet='ROUTPUT', na='#N/A')
cols <- ncol(gdp_waves)

options(warn=-1) # line below produces more than 50 warnings as it produces NAs, which I want




y_real_gap <- as.xts(ts(trendev(gdp_waves), start=c(1965, 4), frequency = 4))

gap_output <- merge(y_real_gap, gap_expost)
names(gap_output) <- c('realtime_gap', 'expost_gap')
                      # philly and st louis gaps, respectively
options(warn=0) # reactivates warnings

cfnai <- fredr_series_observations(series_id = 'CFNAI', frequency = 'q') %>% tbl_xts()

##### Consumption #####
# work in progress
# real_cons_exp <-  fredr_series_observations(series_id = 'PCECC96', frequency = 'q') %>% 

###### List of additional series ########
# mnemonics	desc
# PCECC96	Real personal consumption expenditures
# PCEDG	personal consumption expenditures, durable goods
# PCND	personal consumption expenditures, non durable goods
# PCESV	personal consumption expenditures, services
# DPCCRC1M027SBEA	personal consumption expenditures ex food and energy
# GCEC1	Real Government Consumption Expenditures and Gross Investment
# FDEFX	Federal Government: National Defense Consumption Expenditures and Gross Investment
# GPDI	Gross Private Domestic Investment
# GPDIC1	Real Gross Private Domestic Investment
# PNFI	Private Nonresidential Fixed Investment
# PRFI	Private Residential Fixed Investment
# DRSFRMACBS	Delinquency Rate on Single-Family Residential Mortgages, Booked in Domestic Offices, All Commercial Banks
# EMRATIO	Civilian Employment-Population Ratio
# PAYEMS	All Employees: Total Nonfarm Payrolls
# CES0500000003	Average Hourly Earnings of All Employees: Total Private
# CE16OV	Civilian Employment Level
# LREM25TTUSQ156S	Employment Rate: Aged 25-54: All Persons for the United States
# CES9091000001	All Employees: Government: Federal
# LNS12300060	Employment Population Ratio: 25 - 54 years
# U6RATE	Total unemployed, plus all marginally attached workers plus total employed part time for economic reasons
# UEMPMEAN	Average (Mean) Duration of Unemployment
# AWHNONAG	Average Weekly Hours of Production and Nonsupervisory Employees: Total private
# AWHMAN	Average Weekly Hours of Production and Nonsupervisory Employees: Manufacturing
# CES0600000007	Average Weekly Hours of Production and Nonsupervisory Employees: Goods-Producing 
# CES4300000007	Average Weekly Hours of Production and Nonsupervisory Employees: Transportation and Warehousing
# CEU4200000007	Average Weekly Hours of Production and Nonsupervisory Employees: Retail Trade
# CES4200000007	Average Weekly Hours of Production and Nonsupervisory Employees: Retail Trade  seasonally adj
# CEU3100000007	Average Weekly Hours of Production and Nonsupervisory Employees: Durable Goods




##### SPREADS ####

## BAA 10Y bonds        !!! - DISCONTINUED BY FRED - !!!
spread_baa <- fredr_series_observations(series_id='BAA10Y', frequency='q') %>% tbl_xts()

## BAA 20+ year bonds rate
baa <- fredr_series_observations(series_id = 'BAA', frequency = 'q') %>% tbl_xts()

## AAA 20+ year bonds rate
aaa <- fredr_series_observations(series_id = 'AAA', frequency = 'q') %>% tbl_xts()

## 3 months Tbill rate
tbill_rate_3m <- fredr_series_observations(series_id='TB3MS',frequency='q') %>% tbl_xts()

## 1 year
tbill_rate_1y <- fredr_series_observations(series_id='DGS1', frequency='q') %>% tbl_xts()

## 10 years 
tbill_rate_10y <- fredr_series_observations(series_id='DGS10',frequency='q') %>% tbl_xts()




## spread btw 3m tbill and FFR
tbill3_ffr <- fredr_series_observations(series_id='TB3SMFFM', 
                                        # aggregation_method='eop',
                                        frequency = 'q') %>% tbl_xts()



options("getSymbols.warning4.0"=FALSE)
options("getSymbols.yahoo.warning"=FALSE)# disables disclaimer about version update
# downloads daily prices time series through new Yahoo! API
# to be fixed sooner than later
sp_ret <- getSymbols(src='yahoo', Symbols='^GSPC',
                     from='1950-01-03',
                     to=format(Sys.Date()-1, '%Y-%m-%d'),
                     auto.assign = F)

# aggregating up to quarterly data
# through monthly upscaling
sp_ret <- to.monthly(sp_ret)

# adapts the order of the observations                        <- code for old API
# sp_ret <- sp_ret[order(-1:-nrow(sp_ret)),]
# sp_ret <- data.frame(sp_ret$sp_ret.Close)
# sp_ret <- as.xts(ts(sp_ret[1:(nrow(sp_ret)-1),], start=c(1950, 01), frequency=12))

sp_ret <- diff(log(sp_ret$sp_ret.Close))*100
sp_ret <- as.xts(aggregate(sp_ret, as.yearqtr(as.yearmon(time(sp_ret))), mean))

# one_year <- fredr_series_observations(series_id='DGS1', frequency='q') %>% tbl_xts()
# spread_sp <- (sp_ret - one_year)
spread_sp_3m <- sp_ret - tbill_rate_3m

# replicate spread from moody's
spread_baa_long <- baa - tbill_rate_10y

# quality spread
spread_baa_aaa <- baa - aaa

# corp vs treas spread
spread_aaa <- aaa - tbill_rate_10y



spreads <- merge(spread_baa, spread_sp_3m, 
                 tbill3_ffr, tbill_rate_3m,
                 tbill_rate_1y, tbill_rate_10y,
                 spread_baa_long, spread_baa_aaa,
                 spread_aaa, aaa, baa)
names(spreads) <- c('spread_baa', 'spread_sp_3m',
                    'tbill3_ffr', 'tbill_rate_3m',
                    'tbill_rate_1y', 'tbill_rate_10y',
                    'spread_baa_long', 'spread_baa_aaa',
                    'spread_aaa', 'aaa', 'baa')

options("getSymbols.warning4.0"=T) # activates disclaimer v0.4

##### Fiscal vars #####

#### DEFICIT as % OF GDP

surplus_season <- fredr_series_observations(series_id='M318501Q027NBEA', frequency='q') %>% tbl_xts()

inizio <- as.Date(min(time(surplus_season)), format='%Y-%m-%d')
fine <- as.Date(max(time(surplus_season)), format='%Y-%m-%d')

gdp <- fredr_series_observations(series_id='GDP', frequency='q',
                           observation_start= inizio,
                           observation_end= fine) %>% tbl_xts()

# to deseasonalize one needs to get back to ts format

# surplus_ratio <- 100*surplus/gdp
surplus.ts <- ts((100*surplus_season/gdp), frequency=4,
                       start=c(year(inizio), quarter(inizio)),
                       end=  c(year(fine), quarter(fine)))

surplus_gdp <- decompose(surplus.ts)$x - decompose(surplus.ts)$seasonal
surplus_gdp <- as.xts(surplus_gdp)
names(surplus_gdp) <- 'surplus_gdp'
 
#### DEBT 
# debt to gdp series
debt_gdp <- fredr_series_observations(series_id='GFDEGDQ188S', frequency='q') %>% tbl_xts()

# debt level, millions of $
debt_lev <- fredr_series_observations(series_id='GFDEBTN', frequency='q') %>% tbl_xts()

# debt growth rate to proxy deficit
debt_g <- diff(log(debt_lev))*100


# debt held by FED, billions of $
debt_fed <- fredr_series_observations(series_id = 'FDHBFRBN', frequency='q') %>% tbl_xts()


# percentage of debt held by FED
debt_fed_share <- 100*(debt_fed*1000/debt_lev)


fiscal <- merge(surplus_gdp, debt_g, debt_gdp, debt_fed, debt_fed_share)
names(fiscal) <- c('surplus_gdp', 'debt_growth', 'debt_gdp', 'debt_fed', 'debt_fed_share')

#### MONEY AGGREGATES #####

base <- fredr_series_observations(series_id='BOGMBASE', frequency='q') %>% tbl_xts() %>% `/`(.,1000)
m1 <- fredr_series_observations(series_id='M1SL', frequency='q') %>% tbl_xts()
m2 <- fredr_series_observations(series_id='M2SL', frequency='q') %>% tbl_xts()
m3 <- fredr_series_observations(series_id = 'MABMM301USM189S', frequency = 'q') %>% tbl_xts()

money <- merge(base, m1, m2, m3)
names(money) <- c('base', 'm1', 'm2', 'm3')

# monetary aggregates growth rates
money_g <- diff(log(money))*100
names(money_g) <- c('base_g', 'm1_g', 'm2_g', 'm3_g')

money <-  merge(money, money_g)


#### SPF DATA ####

# automate download of the xlsx file, import, run statistics and merge

# download CPI inflation rate raw file for individuals in the SPF
download.file('https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/files/individual_cpi.xlsx?la=en',
              file.path(temp_dir,'spf_ind_cpi_rate.xlsx'), mode='wb',
              extra='--no-check-certificate',
              quiet = T)

# download CORE CPI inflation rate raw file for individuals in the SPF
download.file('https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/files/individual_corecpi.xlsx?la=en',
              file.path(temp_dir,'spf_ind_corecpi_rate.xlsx'), mode='wb',
              extra='--no-check-certificate',
              quiet = T)

# download PCE inflation rate raw file for individuals in the SPF
download.file('https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/files/individual_pce.xlsx?la=en',
              file.path(temp_dir,'spf_ind_pce_rate.xlsx'), mode='wb',
              extra='--no-check-certificate',
              quiet = T)

# download CORE PCE inflation rate file for individuals in the SPF
download.file('https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/files/individual_corepce.xlsx?la=en',
              file.path(temp_dir,'spf_ind_corepce_rate.xlsx'), mode='wb',
              extra='--no-check-certificate',
              quiet = T)



spf_cpi <- spf_funct('spf_ind_cpi_rate.xlsx', 'CPI',
                     ahead=ahead)
spf_corecpi <- spf_funct('spf_ind_corecpi_rate.xlsx', 'CORECPI',
                         ahead=ahead)
spf_pce <- spf_funct('spf_ind_pce_rate.xlsx','PCE',
                     ahead=ahead)
spf_corepce <- spf_funct('spf_ind_corepce_rate.xlsx', 'COREPCE',
                         ahead=ahead)

spf <- merge(spf_cpi,spf_corecpi,spf_pce, spf_corepce)
# take back here to have tidy spf raw
# raw_spf_cpi <- sapply()



##### Wu-Xia Shadow FFR #####

# This shadow rate ends in 2015 when the ZLB period was terminated.
# using this shadow rate allows the inclusion of 
# non-Taylor rule interventions of the CB in the 
# quantity of money in the policy rate path

download.file(#url = 'http://faculty.chicagobooth.edu/jing.wu/research/data/policyrate.xls', old link
              url = 'https://drive.google.com/uc?export=download&id=1Z0UdvdPhXaIQNVXTZzYtvxDXwxA6YD2e',
              destfile = file.path(temp_dir, 'wuxia_dwnl.xls'),
              mode = 'wb',
              quiet = T)

shffr <- read_excel(path = file.path(temp_dir,'wuxia_dwnl.xls'),
                    col_names = c('date', 'shffr'),
                    sheet = 'Sheet1')

# the above lines download the xls file and dataframe it

#### conversion to xts ####

shffr <- shffr$shffr
shffr <- as.xts(ts(shffr, frequency = 12, start=c(1960, 1)))
shffr <- aggregate(shffr, as.yearqtr(as.yearmon(time(shffr))), mean) %>% as.xts()


# stitching
shffr_ext <- ffr
shffr_ext["2007/2015"] <- shffr["2007/2015"]

shffr <- merge(shffr_ext, lag(shffr_ext, 1))
names(shffr) <- c('shffr', 'shffrb')


# ##### Un. of Michigan Surveys of Consumers #######
# 
# # Glossary:
# # perch is PERcentage CHange; 
# # SOC(M) is Survey Of Consumers (Michigan);
# # g is for Growth rate
# # ind is for INDex
# # diff is for DIFFerence, levels when not otherwise stated
# create list with requests to send
# to the UoM website
req_fields <- list(#'User-Agent' = 'Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/77.0.3865.120 Mobile Safari/537.36',
  table = '32',
  year = '1960',
  qorm = "Q",
  order = 'asc',
  format = "View Below")

# consumers' expected inflation, 1 year ahead
socm_1y_inflation <- POST(url = 'https://data.sca.isr.umich.edu/data-archive/mine.php',
                          body = req_fields,
                          encode = 'form', verbose()
) %>% 
  content('parsed') %>%
  rvest::html_table() %>% 
  .[[2]] %>% 
  select(X13, X15) %>% 
  .[-1,] %>% 
  ts(start = c(1960, 1), frequency = 4) %>%
  xts(order.by = time(.) %>% as.Date()) %>% 
  rename(X13 = socm_1y_inflation_mean,
         X14 = socm_1y_inflation_sd)


# consumers' expected inflation, 5 years ahead
req_fields[['table']] <- '33'

socm_5y_inflation <- POST(url = 'https://data.sca.isr.umich.edu/data-archive/mine.php',
                          body = req_fields,
                          encode = 'form', verbose()
) %>% 
  content('parsed') %>%
  rvest::html_table() %>% 
  .[[2]] %>% 
  select(X13, X15) %>% 
  .[-1,] %>% 
  ts(start = c(1960, 1), frequency = 4) %>%
  xts(order.by = time(.) %>% as.Date()) %>% 
  rename(X13 = socm_5y_inflation_mean,
         X14 = socm_5y_inflation_var)


# consumer confidence indexes + transformations
req_fields[['table']] <- '5'

socm_indexes <- POST(url = 'https://data.sca.isr.umich.edu/data-archive/mine.php',
                     body = req_fields,
                     encode = 'form', verbose()
) %>% 
  content('parsed') %>%
  rvest::html_table() %>% 
  .[[2]] %>% 
  select(X8, X9) %>% 
  .[-1,] %>% 
  ts(start = c(1960, 1), frequency = 4) %>%
  xts(order.by = time(.) %>% as.Date()) %>% 
  rename(X13 = socm_actual_ind,
         X14 = socm_expected_ind) %>% 
  
  # create deviations from expected index
  # percentage deviations from previous level
  # of the indexes, actual and expected
  mutate(soc_diff_ind = socm_actual_ind - socm_expected_ind,
         socm_perch_actual_ind = 100*diff(log(socm_actual_ind)),
         socm_perch_expected_ind = 100*diff(log(socm_expected_ind)))


# merge all datafiles together
SOC_Michigan <- merge(socm_1y_inflation,
                      socm_5y_inflation,
                      socm_indexes)


##### EPU indexes #####
# nb: the US EPU contains data from SPF
# that we already have: using EPU + SPF
# together is not advised, use rather 
# categorical series with the exclusion of
# economic forecasters disagreement

# download excel files
download.file(url = 'http://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.xlsx',
              destfile = file.path(temp_dir, 'EPU_aggregate.xlsx'), 
              mode = 'wb',
              quiet = T)

download.file(url = 'http://www.policyuncertainty.com/media/Categorical_EPU_Data.xlsx',
              destfile = file.path(temp_dir, 'EPU_categories.xlsx'), 
              mode = 'wb',
              quiet = T)

# import data, drop last row
epu_aggregate <- read_excel(path = file.path(temp_dir, 'EPU_aggregate.xlsx'),
                            col_names = T, 
                            sheet = 1,
                            col_types = 'numeric'
                            )
epu_aggregate <- epu_aggregate[-nrow(epu_aggregate),]

epu_aggregate_comp <- read_excel(path = file.path(temp_dir, 'EPU_aggregate.xlsx'),
                                 col_names = T, 
                                 sheet = 2, 
                                 col_types = 'numeric'
                                 )
epu_aggregate_comp <- epu_aggregate_comp[-nrow(epu_aggregate_comp),]

epu_cat <- read_excel(path = file.path(temp_dir, 'EPU_categories.xlsx'),
                      col_names = T,
                      sheet = 1,
                      col_types = c('date', rep('numeric', 12))
                      )

epu_cat <- epu_cat[-c((nrow(epu_cat)-1),nrow(epu_cat)),]

# convert data to xts
epu_aggregate_ts <- ts(data = epu_aggregate, 
                       start=c(epu_aggregate$Year[1], 
                               epu_aggregate$Month[1]),
                       frequency = 12) %>% 
                    aggregate(FUN = mean, nfrequency = 4) %>%
                    xts(x = ., order.by = as.Date(time(.)))



epu_aggregate_comp_ts <- ts(data = epu_aggregate_comp, 
                            start=c(epu_aggregate_comp$Year[1], 
                                    epu_aggregate_comp$Month[1]),
                            frequency = 12) %>% 
                          aggregate(FUN = mean, nfrequency = 4) %>%
                          xts(x = ., order.by = as.Date(time(.)))

epu_cat_ts <- ts(data = epu_cat,
                 start = c(1985, 01),
                 frequency = 12) %>% 
              aggregate(FUN = mean, nfrequency = 4) %>% 
              xts(x = ., order.by = as.Date(time(.)))
# epu_cat_ts <- xts(x = epu_cat, frequency = 12, order.by = epu_cat$Date)

epu_aggregate_ts$Month <- epu_aggregate_ts$Year <- NULL
epu_aggregate_comp_ts$Month <- epu_aggregate_comp_ts$Year <- epu_aggregate_comp_ts$News_Based_Policy_Uncert_Index <- NULL
epu_cat_ts$Date <- NULL

names(epu_aggregate_ts) <- c('epu_baseline_index', 'epu_news_index')
names(epu_aggregate_comp_ts) <- c('epu_public_purchases', 'epu_cpi_spf_disagr', 'epu_tax_expir')
names(epu_cat_ts) <- c('epu_econpol', 'epu_monpol', 'epu_fiscal', 'epu_taxes', 
                       'epu_gvt_spending', 'epu_healthcare', 'epu_natsec',
                       'epu_entitle', 'epu_regul', 'epu_finregul', 'epu_tradepol',
                       'epu_debt_curr')

epu <- merge(epu_aggregate_ts,
             epu_aggregate_comp_ts,
             epu_cat_ts, fill = NA)

#### Quarterly dummies to account for seasonality (JBC)






#### Merge to dataset ####

db_US <- merge(rates, 
               rev_hist_pch,
               rev_hist_yoy,
               unemployment,
               gap_output,
               cfnai,
               spreads,
               money,
               fiscal,
               spf,
               shffr,
               SOC_Michigan,
               epu)

dums <- xts(x = data.frame(q1 = rep(c(1,0,0,0), floor(length(index(db_US)))/4),
                           q2 = rep(c(0,1,0,0), floor(length(index(db_US)))/4),
                           q3 = rep(c(0,0,1,0), floor(length(index(db_US)))/4)),
            order.by = index(db_US)
            )

db_US <- merge(db_US, dums)



write.zoo(x=db_US, 
          file=file.path(data_dir, 'US_data.txt'), 
          sep=';', 
          row.names=F, 
          index.name='time')


# this snippet normalises all vars to standard normal
# JBC suggested this but it really messes up results
# and changes the estimated parameters in an important way.

# db_US <- apply(db_US, FUN = function(x){y=(x-mean(na.omit(x)))/sd(na.omit(x)); return(y)}, MARGIN = 2) %>% as.data.frame()



#### Other countries ####

# UK
# Canada
# Norway
# Denmark
# Japan
# Israel
# Mexico
# Finland
# S Korea
# Sweden





## housekeeping
rm(ffr, classi, core_greenbook, cpi_greenbook, deflator_greenbook,
cpi, core, defl, cpi.mean, core.mean, defl.mean,shffr_ext,
claims, natural_unemp_long, natural_unemp_short,
current_unemp, tot_emp, layoffs, employment_fluct,
cols, gdp_waves, rates, ffrate, 
unemployment, gap_output, cfnai,

spreads, sp_ret, spread_baa, spread_sp_3m,
tbill_rate_3m, tbill_rate_10y, tbill_rate_1y,ffrb,
spread_baa_long, spread_baa_aaa, spread_aaa, aaa, baa,

actual, capacity, y_real_gap, gap_expost, rates.mean,
base, m1, m2, m3, money, money_g, gdp,
inizio, fine, surplus.ts, debt_fed,
debt_fed_share, debt_g, debt_gdp, debt_lev, fiscal,
surplus_gdp, surplus_season, spf, spf_corecpi,
spf_corepce, spf_cpi, spf_pce, rev_hist_yoy, rev_hist_pch,
tbill3_ffr, shffr,

socm_inflation, socm_indexes, socm_indexes, req_fields,
socm_1y_inflation, socm_5y_inflation, SOC_Michigan,

short_long_diff, epu_aggregate, epu_aggregate_comp,
epu_aggregate_comp_ts, epu_aggregate_ts,
epu_cat, epu_cat_ts, epu,

dums)

