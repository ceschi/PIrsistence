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

defl <- pi %>% 
  tk_tbl() %>%
  dplyr::select(rev_defl_pch) %>% 
  mutate(index = pi %>% index() %>% as_date(),
         value = rev_defl_pch,
         rev_defl_pch = NULL) %>% 
  as_tbl_time(index = index) %>% 
  dplyr::filter(!is.na(value))

p1 <- defl %>%
  ggplot(aes(index, value)) +
  geom_point(color = palette_light()[[1]], alpha = 0.5) +
  theme_tq() +
  labs(
    title = "From 1749 to 2013 (Full Data Set)"
  )

p2 <- defl %>%
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

defl %>%
  tidy_acf(value, lags = 1:20)

defl %>%
  tidy_acf(value, lags = 1:20) %>%
  ggplot(aes(lag, acf)) +
  geom_segment(aes(xend = lag, yend = 0), size = 4, color = palette_light()[[1]]) +
  geom_vline(xintercept = 5, size = 1, color = palette_light()[[2]]) +
  annotate("text", label = "5 qtr Mark", x = 5, y = 0.8, 
           color = palette_light()[[2]], size = 3, hjust = 0) +
  theme_tq() +
  labs(title = "ACF: Sunspots")

optimal_lag_setting <- defl %>%
  tidy_acf(value, lags = 1:12) %>%
  filter(acf == max(acf)) %>%
  pull(lag)


# set up the quantities for the blocks

# train on 10y
periods_train <- 4 * 15
# test on 2y
periods_test  <- 4 * 5
# skip_span to 5y
skip_span     <- 4 * 10

rolling_origin_resamples <- rolling_origin(
  defl,
  initial    = periods_train,
  assess     = periods_test,
  cumulative = FALSE,
  skip       = skip_span
)

# Plotting function for a single split
plot_split <- function(split, expand_y_axis = TRUE, alpha = 1, size = 1, base_size = 14) {
  
  # Manipulate data
  train_tbl <- training(split) %>%
    add_column(key = "training") 
  
  test_tbl  <- testing(split) %>%
    add_column(key = "testing") 
  
  data_manipulated <- bind_rows(train_tbl, test_tbl) %>%
    as_tbl_time(index = index) %>%
    mutate(key = fct_relevel(key, "training", "testing"))
  
  # Collect attributes
  train_time_summary <- train_tbl %>%
    tk_index() %>%
    tk_get_timeseries_summary()
  
  test_time_summary <- test_tbl %>%
    tk_index() %>%
    tk_get_timeseries_summary()
  
  # Visualize
  g <- data_manipulated %>%
    ggplot(aes(x = index, y = value, color = key)) +
    geom_line(size = size, alpha = alpha) +
    theme_tq(base_size = base_size) +
    scale_color_tq() +
    labs(
      title    = glue("Split: {split$id}"),
      subtitle = glue("{train_time_summary$start} to {test_time_summary$end}"),
      y = "", x = ""
    ) +
    theme(legend.position = "none") 
  
  if (expand_y_axis) {
    
    sun_spots_time_summary <- defl %>% 
      tk_index() %>% 
      tk_get_timeseries_summary()
    
    g <- g +
      scale_x_date(limits = c(sun_spots_time_summary$start, 
                              sun_spots_time_summary$end))
  }
  
  return(g)
}

# Plotting function that scales to all splits 
plot_sampling_plan <- function(sampling_tbl, expand_y_axis = TRUE, 
                               ncol = 3, alpha = 1, size = 1, base_size = 14, 
                               title = "Sampling Plan") {
  
  # Map plot_split() to sampling_tbl
  sampling_tbl_with_plots <- sampling_tbl %>%
    mutate(gg_plots = map(splits, plot_split, 
                          expand_y_axis = expand_y_axis,
                          alpha = alpha, base_size = base_size))
  
  # Make plots with cowplot
  plot_list <- sampling_tbl_with_plots$gg_plots 
  
  p_temp <- plot_list[[1]] + theme(legend.position = "bottom")
  legend <- get_legend(p_temp)
  
  p_body  <- plot_grid(plotlist = plot_list, ncol = ncol)
  
  p_title <- ggdraw() + 
    draw_label(title, size = 18, fontface = "bold", colour = palette_light()[[1]])
  
  g <- plot_grid(p_title, p_body, legend, ncol = 1, rel_heights = c(0.05, 1, 0.05))
  
  return(g)
  
}

rolling_origin_resamples$splits[[1]] %>%
  plot_split(expand_y_axis = TRUE) +
  theme(legend.position = "bottom")

rolling_origin_resamples %>%
  plot_sampling_plan(expand_y_axis = T, ncol = 3, alpha = 1, size = 1, base_size = 10, 
                     title = "Backtesting Strategy: Rolling Origin Sampling Plan")



# first LSTM on a single slice of data, the most recent
split    <- rolling_origin_resamples$splits[[6]]
split_id <- rolling_origin_resamples$id[[6]]

plot_split(split, expand_y_axis = FALSE, size = 0.5) +
  theme(legend.position = "bottom") +
  ggtitle(glue("Split: {split_id}"))

# get data ready:
# create tibble with key identifying train vs test
df_trn <- training(split)
df_tst <- testing(split)

df <- bind_rows(
  df_trn %>% add_column(key = "training"),
  df_tst %>% add_column(key = "testing")
) %>% 
  as_tbl_time(index = index)