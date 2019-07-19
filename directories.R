##### Directories creation #####


working_directory <- getwd()
temp_dir <- 'Downloaded files'
data_dir <- 'Processed data'
graphs_dir <- 'Plots'

temp_dir <- file.path(working_directory, temp_dir)
data_dir <- file.path(working_directory, data_dir)
graphs_dir <- file.path(working_directory, graphs_dir)

options(warn=-1) # turns off warnings momentarily
dir.create(temp_dir)
dir.create(data_dir)
dir.create(graphs_dir)
options(warn=0) # turns warnings back on