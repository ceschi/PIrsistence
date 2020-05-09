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


##### create fake linear data ##################################################

# just a sequence
fake <- seq(1, 10000, 1)

# make it time series
fake <-fake %>% 
    as_tibble() %>% 
    mutate(date = as_date(index(.))) %>% 
    as_tbl_time(date)


# plot just because
fake %>% ggplot() + geom_line(aes(x = date, y = value))


##### splitting ################################################################

t_train <- 750
t_test <- 250
skip_span <- 10000

rolling_samples <- rolling_origin(fake,
                                  initial = t_train,
                                  assess = t_test, 
                                  cumulative = F, 
                                  skip = skip_span)


# pull out just the first split
split1 <- rolling_samples$splits[[1]]

# combine with labels in
# one df

df <- bind_rows(
  training(split1) %>% add_column(key = 'train'),
  testing(split1) %>% add_column(key = 'test')
) %>% 
  as_tbl_time(index = date)


# rescale and center
rec <- recipe(value ~ ., df) %>% 
  step_center(value) %>% 
  step_scale(value) %>% 
  prep()

df_proc <- bake(rec, df)

rec_mean <- rec$steps[[1]]$means
rec_sd <- rec$steps[[2]]$sds


######' *NOW THE FUCKING LSTM* ############################################


lag_set <- 10
batch <- 50
train_length <- 500
tsteps <- 1

epochs <- 300

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
set.seed(24051988)
model <- keras_model_sequential()
model %>%
  layer_lstm(units = 25, 
             # activation = 'relu',
             input_shape = c(tsteps, 1), 
             batch_size = batch, 
             return_sequences = T, 
             stateful = T) %>% 
  # layer_lstm(units = 35,
  #            return_sequences = T,
  #            stateful = T) %>%
  # layer_lstm(units =17,
  #            return_sequences = T,
  #            stateful = T) %>%
  # layer_gru(units = 20,
  #           return_sequences = T,
  #           stateful = T,
  #           batch_size = batch) %>%
  layer_lstm(units = 15,
             return_sequences = F,
             stateful = T) %>%
  layer_dense(units = 1)

model %>% compile(loss = 'mae', 
                  # metrics = 'mse', 
                  optimizer = 'adam')

for (i in 1:epochs+200){
model %>% fit(x = x_train_arr,
              y = y_train_arr,
              batch_size = batch,
              epochs = 1,
              verbose = 1,
              shuffle = F)
model %>% reset_states()

cat('\nIteration ', i, ' completed out of ', epochs,'.\n')
}

# these two should be equivalent,
# I guess

model_one <- keras::clone_model(model)
model_one %>% compile(loss = 'mae', 
                      # metrics = 'mse', 
                      optimizer = 'adam')
model_one %>% fit(x = x_train_arr,
              y = y_train_arr,
              batch_size = batch,
              epochs = epochs+200,
              verbose = 1,
              shuffle = F)

##### predictions

pred <- model %>% 
  predict(x_test_arr, batch_size = batch)

t_train <- tibble(value = df_proc %>% 
                    filter(key == 'train') %>% 
                    select(value) %>% 
                    mutate(value = value*rec_sd + rec_mean), 
                  date = df_proc %>% 
                    filter(key == 'train') %>% 
                    select(date)) %>% 
  add_column(key = 'actual_train')

t_test <- tibble(value = df_proc %>% 
                   filter(key == 'test') %>% 
                   select(value) %>% 
                   mutate(value = value*rec_sd + rec_mean), 
                 date = df_proc %>% 
                   filter(key == 'test') %>% 
                   select(date)) %>% 
  add_column(key = 'actual_test')

t_pred <- tibble(value = pred[,1]*rec_sd + rec_mean) %>% 
  add_column(key = 'predicted', 
             date = df_proc %>% 
               filter(key == 'test') %>% 
               select(date))

what <- rbind(t_test, t_train, t_pred)
what <- tibble(value = what$value$value, date = what$date$date, key = what$key)
what <- arrange(what, date, key)

what %>% ggplot(aes(x = date, colour = as.factor(key)))+geom_point(aes(y=value))
# what %>% ggplot(aes(x = date, colour = as.factor(key)))+geom_line(aes(y=value))

## sum-up: cool, fit out of the loop as it appears to give better results here
#          then experiment with epochs. Next step is to normalise and store 
#          seprately test and train, then figure out how to have longer horizon
#          of prediction with self feeding. then move to real data.