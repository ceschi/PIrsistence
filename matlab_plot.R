library(R.matlab)
library(ggridges)

bsum <- R.matlab::readMat(file.choose()) %>% as.data.frame()

temp <- pi$rev_defl_pch
temp <- temp[!is.na(temp),]
names(bsum) <- temp[-(1:40), ] %>% index()

rm(temp)

bsum <- bsum %>% reshape2::melt() %>% 
  group_by(variable) %>% mutate(med = median(value),
                                mea = mean(value)) %>% 
  ungroup()

p <- bsum %>% 
  as.tibble() %>% 
  ggplot(aes(x = value, y = variable, group = variable)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.05, .5, 0.95),
                      scale = .95, rel_min_height = .05) + 
  coord_flip() + theme_ridges() + scale_y_yearqtr(n = 15) +
  geom_line(aes(x = med, y = value), size = 2) + 
  geom_line(aes(x = mea, y = value), size = 2)