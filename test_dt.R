# script to create synthetic DT to test ACC stuff

library(data.table)
library(tidyverse)
library(Hmisc)

# total observations
n <- 400000000

# years  <- 2000:2015
# month <- 1:12
# siren <- runif(50, 100, 1000)  %>% floor()
# id_bank <- sample(letters, 10, T) %>% paste0(., sample(letters, 10, T)) %>% paste0(., runif(10, 1, 10) %>% floor())
# credit <- rcauchy(n, 10000, 1) %>% abs() %>% `*`(., (1:n)*10/1000) # last part adds some growth in credit
# risk <- 1:6


dt <- list(year = 2000:2020 %>%  sample(replace = T, ., n),
           month = 1:12%>%  sample(replace = T, ., n),
           siren = runif(50, 100, 1000)  %>% floor() %>%  sample(replace = T, ., n),
           id_bank = sample(letters, 1000, T) %>% paste0(., sample(letters, 1000, T)) %>% paste0(., runif(1000, 1, 10) %>% floor()) %>%  sample(replace = T, ., n),
           credit =  rcauchy(n, 10000, 1) %>% abs() %>% `*`(., (1:n)*10/1000), # last part adds some growth in credit,
           risk = 1:6 %>% sample(replace = T, ., n)
) %>% setDT

# rm(credit, id_bank, month, risk, siren, years)


##### do something to pick up later ############################################

dt[risk == 4][year>2010, credit := 1:.N + credit]

# split dt by year and write to disk

for (i in 2000:2020){
  temp <- dt[year == i]
  fwrite(x = temp,
         file = file.path("D:/dt", paste0(i, "_dt_yearly.rds")),
         sep = ";",
         append = F)
  rm(temp)
}

rm(dt, n, i)

gc()

##### Import data back with loop  ##############################################
folder_data <- file.path("D:/dt/")
files_list <- list.files(folder_data, full.names = T)

tic('for loop')
dt <- fread(file = files_list[1], nrows = 0, sep = ';', header = T)

for (i in 1:length(files_list)){
  temp <- fread(file = files_list[i], sep = ';', header = T)
  dt <- rbind(dt, temp)
  rm(temp)
}
toc()

tic('pmap')
# plan(multiprocess)
list_DT <- future_pmap(.l = list(file = as.list(files_list),
                                  sep = rep(list(';'), length(files_list)),
                                  header = rep(list(T), length(files_list))),
                    .f = fread)
toc()



##### now some data wrangling ##################################################
dt[, keyby = c("year", "month", "id_bank"), tot_cr := sum(credit)
   ][
     risk == 4, by = c("year", "month", "id_bank"), tot4_cr := sum(credit)
     ][
       is.na(tot4_cr), tot4_cr := 0
       ][
         , keyby = c("year", "month", "id_bank"), share4 := tot4_cr/tot_cr
         ]
dt[, keyby = c("year", "month"), p50_thres := quantile(share4, .5)]
dt[, dummy := 0]
dt[share4>p50_thres, dummy := 1][, keyby = c("dummy", "year", "month"), a_tot := log(sum(as.numeric(tot_cr)))]
dt[, tot4_cr := NULL]


# id_dummy <- dt[dummy == 1, by = c("year", "month", "dummy"), .(id_bank)]

# ddd <- dt[id_bank %in% id_dummy$id_bank] # too few banks


##### Growth rates #############################################################

dt[, keyby = c('year', 'month', 'id_bank', 'risk'), cr_byrisk := sum(credit)]
dt[order(year, month, id_bank, risk)]

dt[, keyby = c('id_bank', 'risk'), g_cr_risk := c(NA, 100*diff(log(cr_byrisk)))]

for (i in unique(dt$banks)){
  for (r in unique(dt$risk)){
    temp <- dt[id_bank == i][risk==r][order(year, month, risk)][, .(year, month, id_bank, risk, cr_byrisk)][, g_cr_risk := c(NA, 100*diff(log(cr_byrisk)))][, .(year, month, id_bank, risk, cr_byrisk, g_cr_risk)]
  }
}



dt[, keyby = c('year', 'month', 'id_bank', 'risk'), .(mean(credit))]





dt[, keyby = c('year', 'month', 'id_bank', 'risk'), .(cane = sum(credit))] -> bra
bra[, keyby = c('id_bank', 'risk'), .(year, month, ert = c(NA, 100*diff(log(cane))))]