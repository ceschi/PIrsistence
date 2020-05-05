###### Import data -------------------------------------------------------------
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

lags <- function(var, n=10){
  ## use as follows:
  ##    mutate( !!!lags(x, 3), !!!lags(y,3) )
  
  var <- enquo(var)
  
  indices <- seq_len(n)
  map( indices, ~quo(lag(!!var, !!.x)) ) %>% 
    set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
  
}

# adaptation from lsmt to defl
pi <- readr::read_delim(file = './processed_data/PI_data.csv',
                        delim = ';') %>% 
  dplyr::mutate(date = as.Date(as.yearqtr(date))) %>% 
  tbl2xts::tbl_xts()

defl <- pi %>%
  tk_tbl() %>%
  dplyr::select(rev_defl_pch) %>%
  mutate(
    index = pi %>% index() %>% as_date(),
    value = rev_defl_pch,
    rev_defl_pch = NULL
  ) %>%
  as_tbl_time(index = index) %>%
  dplyr::filter(!is.na(value))


##### do some stuff ------------------------------------------------------------

defl_splits <- rolling_origin(defl,
                               initial = 200,
                               assess = 50,
                               skip = 1000,
                               lag = 0,
                               cumulative = FALSE)

# recompose df
split <- defl_splits$splits[[1]]
split_id <- defl_splits$id[[1]]

df <- bind_rows(
  training(split) %>% add_column(key = 'train'),
  testing(split) %>% add_column(key = 'test')
) %>% as_tbl_time(index = index)


# preprocess via centering
recipe_df <- recipe(value ~ ., df) %>% 
  step_center(value) %>% 
  step_scale(value) %>% 
  prep()

df_processed <- bake(recipe_df, df)

# store previous values
df_mean <- recipe_df$steps[[1]]$means[['value']]
df_ssd <- recipe_df$steps[[2]]$sds[['value']]


##### *BUILDING UP THE LSTM* ###################################################
##### setup the data
# how many lags ?
lags <- 1
# how many obs to feed to training in one forward/backward pass before weight update?
batch_size <- 10

# 
train <- 200
testi <- 50

# number of lags included in the training/testing set
tsteps <- 1

# epochs are iterations
epochs <- 30

# setup the arrays
# X array: the predictors
#       it needs to be a 3d array, with [samples, timesteps, features]
#       dimensions. Again, the first dimension is the lenght
#       of values (?), the second the lags, the third is the
#       dimensionality of the predictors, ie 1 if univariate,
#       n if multivariate
# 
# Y array: the target/dependent variable
#       it needs to be a 2d array, with dimensions [samples, timesteps]
#       in this case it includes the y variables itself as lag?
#
# NB: are NAs values accepted in these arrays?


# training arrays
x_training_array <- array(
  data = df_processed %>% 
    filter(key == 'train') %>% 
    select(value) %>% 
    mutate(value1 = lag(value, n = 1),
           value2 = lag(value, n = 2),
           value3 = lag(value, n = 3),
           value4 = lag(value, n = 4),
           value5 = lag(value, n = 5)
           ) %>%  na.omit(.) %>%
    # drop time-t obs as it will 
    # be in the y array
    select(-value)%>% 
    pull(),
  dim = c(train,
          lags,
          tsteps)
)

y_training_array <- array(
  data = df_processed %>% 
    filter(key == 'train') %>% 
    select(value)%>% 
    pull(),
  dim = c(train,
          tsteps)
)

# testing arrays
x_test_array <- array(
  data = df_processed %>% 
    filter(key == 'test') %>% 
    select(value) %>% 
    mutate(value1 = lag(value, n = 1),
           value2 = lag(value, n = 2),
           value3 = lag(value, n = 3),
           value4 = lag(value, n = 4),
           value5 = lag(value, n = 5)
          ) %>% na.omit(.) %>% 
    # drop time-t obs as it will 
    # be in the y array
    select(-value) %>% 
    pull(),
  dim = c(testi,
          lags,
          tsteps)
)

y_test_array <- array(
  data = df_processed %>% 
    filter(key == 'test') %>% 
    select(value)%>% 
    pull(),
  dim = c(testi,
          tsteps)
)


##### SETUP tHE MODEL ##########################################################

# declare model
lstm <- keras::keras_model_sequential()

# stack layers - just one to start with
lstm %>% 
  keras::layer_lstm(units = 20,
                    input_shape = c(tsteps, lags),
                    batch_size = batch_size,
                    return_sequences = F,
                    stateful = T) %>% 
  # keras::layer_lstm(units = 20,
  #                   return_sequences = F,
  #                   stateful = T) %>%
  #                   
  # the last layer shall be 'dense'
  # with one unit, it outputs only
  # value that shall be the forecast
  # for next value
  keras::layer_dense(units = 1)


# compile the model
# also consider loss = 'mse' for mean square error
# instead of 'mae'
lstm %>% keras::compile(loss = 'mse', optimizer = 'adam')

summary(lstm)


##### now fit the model (finger crossed) ---------------------------------------


lstm %>% keras::fit(x = x_training_array, 
             y = y_training_array, 
             batch_size = batch_size, 
             epochs = 3000,#epochs, 
             verbose = T,
             shuffle = F)


##### Forecasts ----------------------------------------------------------------

pred <- lstm %>% predict(x_test_array,
                         batch_size = batch_size)

df_pred <- data.frame(test = pred[,1], true = x_test_array) %>% 
  apply(MARGIN = 2, FUN = function(x) x*df_ssd + df_mean) %>% 
  as_tibble()

  

df_pred %>% ggplot(aes(x = index(.))) + 
  geom_point(aes(y = test, colour = 'test')) + 
  geom_point(aes(y = true, colour = 'true'))


# sort of there, now lets try 
# with more units, one layer
############################################################################################


# declare model
lstm <- keras::keras_model_sequential()

# stack layers - just one to start with
lstm %>% 
  keras::layer_lstm(units = 100,
                    input_shape = c(tsteps, lags),
                    batch_size = batch_size,
                    return_sequences = F,
                    stateful = T) %>% 

  # the last layer shall be 'dense'
  # with one unit, it outputs only
  # value that shall be the forecast
  # for next value
  keras::layer_dense(units = 1)


# compile the model
# also consider loss = 'mse' for mean square error
# instead of 'mae'
lstm %>% keras::compile(loss = 'mae', optimizer = 'adam', metrics = 'accuracy')

summary(lstm)


##### now fit the model (finger crossed) ---------------------------------------

tic('lstm looong')
history <- lstm %>% keras::fit(x = x_training_array, 
                    y = y_training_array, 
                    batch_size = batch_size, 
                    epochs = 3000,#epochs, 
                    verbose = T,
                    shuffle = F)
toc()
plot(history)

##### Forecasts ----------------------------------------------------------------

pred <- lstm %>% predict(x_test_array,
                         batch_size = batch_size)

df_pred <- data.frame(test = pred[,1], true = x_test_array) %>% 
  apply(MARGIN = 2, FUN = function(x) x*df_ssd + df_mean) %>% 
  as_tibble()



df_pred %>% ggplot(aes(x = index(.))) + 
  geom_point(aes(y = test, colour = 'test'), size = 2) + 
  geom_line(aes(y = test, colour = 'test')) +
  geom_point(aes(y = true, colour = 'true'), size = 2) + 
  geom_line(aes(y = true, colour = 'true')) +
  theme_minimal()
