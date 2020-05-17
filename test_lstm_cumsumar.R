##### LSTM for fools ###########################################################
library(recipes)
library(tbl2xts)
library(timetk)
library(tidyverse)
library(rlang)
library(lubridate)
library(xts)
library(tibbletime)
library(rsample)
library(keras)
library(tictoc)



##### Some more complicated data ##############################################


fake <- arima.sim(model = list(ar = c(.8, .1, 0.02)), n = 10000)  %>% cumsum 

# make it time series
fake <-fake %>% 
  as_tibble() %>% 
  mutate(date = as_date(index(.))) %>% 
  # rename(value = x) %>% 
  as_tbl_time(date)

##### splitting ################################################################

t_train <- 7500
t_test <- 2500
skip_span <- 10000

# rolling origin creates instructions to 
# later split data, especially good 
# for backtesting, could be done in
# a rough way too
rolling_samples <- rolling_origin(fake,
                                  initial = t_train,
                                  assess = t_test, 
                                  cumulative = F, 
                                  skip = skip_span)

# pull out just the first split
# in this case it's all samples, too
split1 <- rolling_samples$splits[[1]]

# combine with labels in
# one df

df <- bind_rows(
  training(split1) %>% add_column(key = factor('train',
                                               c('train', 'test'))
  ),
  testing(split1) %>% add_column(key = factor('test',
                                              c('train', 'test'))
  )
) %>% 
  as_tbl_time(index = date)

# pull out the two groups 
df_train <- df %>% filter(key == 'train')
df_test <- df %>% filter(key == 'test')

# recipe to rescale, center data
rec_train <- recipe(value ~ ., df_train) %>% 
  step_center(value) %>% 
  step_scale(value) %>% 
  prep()

# this is done *SEPARATELY* for
# train and test parts
rec_test <- recipe(value ~ ., df_test) %>% 
  step_center(value) %>% 
  step_scale(value) %>% 
  prep()


# compose together the rescaled
# observations
df_proc <- bind_rows(
  bake(rec_train, df_train),
  bake(rec_test, df_test)
) %>% 
  as_tbl_time(index = date)

# store means & SDs
train_mean <- rec_train$steps[[1]]$means
train_sd <- rec_train$steps[[2]]$sds

test_mean <- rec_test$steps[[1]]$means
test_sd <- rec_test$steps[[2]]$sds


######' *NOW THE FUCKING LSTM* ############################################

# how much past use?
lag_set <- 50

# how many observations feed?
batch <- 50

# no idea here, really
train_length <- 6000

# still, not clear
tsteps <- 1

# for how long train the model?
epochs <- 500

lag_train <- df_proc %>% 
  mutate(value_lag = lag(value, lag_set)) %>% 
  filter(!is.na(value_lag)) %>% 
  filter(key == 'train') %>%
  tail(train_length)

x_train_vec <- lag_train$value_lag
x_train_arr <- array(data = x_train_vec,
                     dim = c(length(x_train_vec), tsteps, 1))

y_train_vec <- lag_train$value
y_train_arr <- array(data = y_train_vec, 
                     dim = c(length(y_train_vec), tsteps))


# Testing Set
lag_test <- df_proc %>%
  mutate(
    value_lag = lag(value, n = lag_set)
  ) %>%
  filter(!is.na(value_lag)) %>%
  filter(key == "test")

x_test_vec <- lag_test$value_lag
x_test_arr <- array(data = x_test_vec, dim = c(length(x_test_vec), tsteps, 1))

y_test_vec <- lag_test$value
y_test_arr <- array(data = y_test_vec, dim = c(length(y_test_vec), tsteps))


##### model stuff ##############################################################
if (exists('model')) rm(model) 


model <- keras_model_sequential()
model %>%
  layer_lstm(units = 50, 
             # activation = 'relu',
             input_shape = c(tsteps, 1), 
             batch_size = batch, 
             return_sequences = T, 
             stateful = T) %>%
  # layer_lstm(units = 50,
  #            return_sequences = F,
  #            stateful = T) %>%
  # layer_lstm(units =100,
  #            return_sequences = T,
  #            stateful = T) %>%
  # layer_gru(units = 20,
  #           return_sequences = T,
  #           stateful = T,
  #           batch_size = batch) %>%
  layer_lstm(units = 55,
             return_sequences = F,
             stateful = T) %>%
  layer_dense(units = 1)

model %>% compile(loss = 'mse',
                  optimizer = 'adam',
                  metrics = 'mae')


tic('model fit')
hist <- model %>% fit(x = x_train_arr,
                      y = y_train_arr,
                      batch_size = batch,
                      epochs = epochs,
                      verbose = 1,
                      shuffle = F)
toc()

##### predictions

pred <- model %>% 
  predict(x_test_arr, batch_size = batch) 

t_train <- tibble(value = df_proc %>% 
                    filter(key == 'train') %>% 
                    select(value) %>% 
                    mutate(value = value*train_sd + train_mean), 
                  date = df_proc %>% 
                    filter(key == 'train') %>% 
                    select(date)) %>% 
  add_column(key = 'actual_train')

t_test <- tibble(value = df_proc %>% 
                   filter(key == 'test') %>% 
                   select(value) %>% 
                   mutate(value = value*test_sd + test_mean), 
                 date = df_proc %>% 
                   filter(key == 'test') %>% 
                   select(date)) %>% 
  add_column(key = 'actual_test')

t_pred <- tibble(value = pred*test_sd + test_mean) %>% 
  add_column(key = 'predicted', 
             date = df_proc %>% 
               filter(key == 'test') %>% 
               select(date))

output <- rbind(t_test, t_train, t_pred)
output <- tibble(value = output$value$value, date = output$date$date, key = output$key)
output <- arrange(output, date, key)

output %>% ggplot(aes(x = date, colour = as.factor(key)))+geom_line(aes(y=value), size = 1, alpha = .5) + ggtitle('model2')
