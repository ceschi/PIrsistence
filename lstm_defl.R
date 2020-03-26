# Core Tidyverse
library(tidyverse)
library(glue)
library(forcats)

# Time Series
library(timetk)
library(tidyquant)
library(tibbletime)

# Visualization
library(cowplot)

# Preprocessing
library(recipes)

# Sampling / Accuracy
library(rsample)
library(yardstick) 

# Modeling
library(keras)

# adaptation from lsmt to defl

sun_spots <- pi %>% 
  tk_tbl() %>%
  dplyr::select(rev_defl_pch) %>% 
  mutate(index = pi %>% index() %>% as_date(),
         value = rev_defl_pch,
         rev_defl_pch = NULL) %>% 
  as_tbl_time(index = index) %>% 
  dplyr::filter(!is.na(value))

p1 <- sun_spots %>%
  ggplot(aes(index, value)) +
  geom_point(color = palette_light()[[1]], alpha = 0.5) +
  theme_tq() +
  labs(
    title = "From 1749 to 2013 (Full Data Set)"
  )

p2 <- sun_spots %>%
  filter_time("start" ~ "1980") %>%
  ggplot(aes(index, value)) +
  geom_line(color = palette_light()[[1]], alpha = 0.5) +
  geom_point(color = palette_light()[[1]]) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE) +
  theme_tq() +
  labs(
    title = "1749 to 1800 (Zoomed In To Show Cycle)",
    caption = "datasets::sunspot.month"
  )

p_title <- ggdraw() + 
  draw_label("Sunspots", size = 18, fontface = "bold", colour = palette_light()[[1]])

plot_grid(p_title, p1, p2, ncol = 1, rel_heights = c(0.1, 1, 1))

tidy_acf <- function(data, value, lags = 0:20) {
  acf_values <- data %>%
    pull(value) %>%
    acf(lag.max = tail(lags, 1), plot = FALSE) %>%
    .$acf %>%
    .[,,1]
  
  ret <- tibble(acf = acf_values) %>%
    rowid_to_column(var = "lag") %>%
    mutate(lag = lag - 1) %>%
    filter(lag %in% lags)
  
  return(ret)
}

sun_spots %>%
  tidy_acf(value, lags = 0:20)

sun_spots %>%
  tidy_acf(value, lags = 1:20) %>%
  ggplot(aes(lag, acf)) +
  geom_segment(aes(xend = lag, yend = 0), size = 4, color = palette_light()[[1]]) +
  geom_vline(xintercept = 5, size = 1, color = palette_light()[[2]]) +
  annotate("text", label = "5 qtr Mark", x = 5, y = 0.8, 
           color = palette_light()[[2]], size = 3, hjust = 0) +
  theme_tq() +
  labs(title = "ACF: Sunspots")

optimal_lag_setting <- sun_spots %>%
  tidy_acf(value, lags = 1:12) %>%
  filter(acf == max(acf)) %>%
  pull(lag)

