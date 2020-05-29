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
rm(list = ls())
gc(full = T, verbose = T)

##### Part I: linear data ######################################################

# just a sequence
# fake <- seq(1, 10000, 1) %>% log() %>% `+`(., rnorm(10000, 0, .10))
# fake <- cos(seq(1, 10000, 1)*.005)
# fake <- rep(1, 10000) + 
#   rnorm(n = 10000, 0, .01) + (1:10000)/10000
  # arima.sim(n = 10000, model = list(ar = c(.995)), sd = .001)

len <- 100000
seq <- seq(1, len, 1)
fake <- (cos(seq*.005*base::pi) +
  rnorm(n = len, 0, .01) -
  seq/len + 10 +
  sin(seq*.008*base::pi)) + log(seq)

fake <- seq
# fake <- seq + rnorm(n = len, 0, 10)

fake %>% ts.plot
# make it time series
fake <-fake %>% 
  as_tibble() %>% 
  mutate(date = as_date(index(.))) %>% 
  # rename(value = x) %>% 
  as_tbl_time(date)

##### splitting ################################################################

l_train <- floor(len*.75)
l_test <- floor(len*.25)
skip_span <- len + 2

# rolling origin creates instructions to 
# later split data, especially good 
# for backtesting, could be done in
# a rough way too
rolling_samples <- rolling_origin(fake,
                                  initial = l_train,
                                  assess = l_test, 
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
# for how long train the model?
epochs <- 300


# how much past use?
lag_set <- 1

# how many observations feed?
batch <- 1000

# no idea here, really
train_length <- 70000

# still, not clear
tsteps <- 20
# tsteps here seem to have a bad effect on overall 
# fit and prediction: lower values entail 
# faster fit and better predictions in general


lag_train <- df_proc %>% 
  filter(key == 'train') %>%
  mutate(value_lag = lag(value, lag_set)) %>% 
  filter(!is.na(value_lag)) %>% 
  tail(train_length)

x_train_vec <- lag_train$value_lag
x_train_arr <- array(data = x_train_vec,
                     dim = c(length(x_train_vec), tsteps, 1))

y_train_vec <- lag_train$value
y_train_arr <- array(data = y_train_vec, 
                     dim = c(length(y_train_vec), tsteps))


# Testing Set
lag_test <- df_proc %>%
  filter(key == 'test') %>% 
  mutate(
    value_lag = lag(value, n = lag_set)
  ) %>%
  filter(!is.na(value_lag))

x_test_vec <- lag_test$value_lag
x_test_arr <- array(data = x_test_vec, dim = c(length(x_test_vec), tsteps, 1))

y_test_vec <- lag_test$value
y_test_arr <- array(data = y_test_vec, dim = c(length(y_test_vec), tsteps))


##### model stuff ##############################################################
if (exists('model')) rm(model)


model <- keras_model_sequential()
model %>%
  layer_lstm(units = 25,
             # activation = 'tanh',
             # activation = 'relu',
             input_shape = c(tsteps, 1), 
             batch_size = batch,
             return_sequences = F, 
             stateful = T) %>% 
  # layer_lstm(units = 5,
  #            return_sequences = T,
  #            stateful = T) %>%
  # layer_lstm(units = 150,
  #            return_sequences = T,
  #            stateful = T) %>%
  # layer_lstm(units = 75,
  #            return_sequences = F,
  #            stateful = T) %>%
  layer_dense(units = 1)

model %>% compile(
                  loss = 'mae',
                  # loss = 'mse',
                  # metrics = 'mae',
                  optimizer = 'adam')
summary(model)

tic('model fit')
history <- model %>% fit(x = x_train_arr,
                         y = y_train_arr,
                         batch_size = batch,
                         epochs = epochs,
                         callbacks = list(
                           callback_early_stopping(monitor = 'val_loss',
                                                   mode = 'min',
                                                   patience = 50,
                                                   min_delta = .00001)
                         ),
                         validation_split = .1,
                         verbose = 1,
                         shuffle = F)
toc()

##### predictions

pred <- model %>% 
  predict(x_test_arr, batch_size = batch) %>% .[,1]

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
               filter(key == 'test', date != last(date)) %>% 
               select(date))

output <- rbind(t_test, t_train, t_pred)
output <- tibble(value = output$value$value, date = output$date$date, key = output$key)
output <- arrange(output, date, key)

tail(output, 7500) %>% ggplot(aes(x = date, colour = as.factor(key)))+
  geom_line(aes(y=value), size = 1.2, alpha = .5) +
  ggtitle('model1') + theme_minimal()#+geom_smooth(method = 'loess', aes(x = date, y = value, colour = as.factor(key)))


##### online estimation and forecast ###########################################

# horizon of forecast
h_forecast <- l_test + 1

pred_h <- array(data = NA, dim = c(h_forecast, 1))
pred_h[1,] <- x_train_arr[nrow(x_train_arr), ncol(x_train_arr), 1]


for (i in 2:h_forecast){
  pred_h[i,] <- predict(model,
                        x = array(data = pred_h[(i-1), 1],
                                  dim = c(1, tsteps, 1)),
                        batch_size = 1)
    # prediction plateaus after about 32
    # iterations, might well be model
    # depentent factor
}


  
pred_h_trans <- pred_h*test_sd+test_mean

pred_h_trans %>% ts.plot