##### Export to MATLAB #####
# quick dirty use of R&P07 codes


write.table(x = pi[, grep(pattern = '_pch', names(pi))], 
            file = 'D:/emanu/OneDrive/Matlab/infl_pers_07/brute_force/pi_data.csv', 
            row.names = F, 
            quote = F, 
            na = "NaN", 
            col.names = F)


write.table(x = pi %>% names() %>% grep(pattern = '_pch', value = T), 
            file = 'D:/emanu/OneDrive/Matlab/infl_pers_07/brute_force/pi_names.csv', 
            row.names = F, 
            quote = F, 
            na = "NaN", 
            col.names = F)

nn <- list()
for (i in 1:length(names(pi))){
  temp <- pi[, i]
  temp <- temp[!is.na(temp)]
  ind <- index(temp) %>% dplyr::first()
  nn[[i]] <- ind %>% year()
}

write.table(x = matrix(nn),
            file = 'D:/emanu/OneDrive/Matlab/infl_pers_07/brute_force/pi_dates.csv',
            row.names = F,
            col.names = F)

rm(temp, ind, i, nn)
