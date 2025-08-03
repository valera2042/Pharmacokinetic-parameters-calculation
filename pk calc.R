# Step by step guide for the calculation of pharmacokinetic parameters

# install packages if required
requiredPackages <- c("here, writexl, pROC, readxl, tidyverse,
                      lubridate, hash, zoo, glue, gghighlight,
                      psych, ggplot2", 'openxlsx','DiagnosisMed',
                      'Metrics', 'caTools', 'gt')
for (package in requiredPackages) { #Installs packages if not yet installed
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(package)
}

# upload packges
library(here) # to replace absolute path on relative
library(writexl)
library(pROC) 
library(readxl)
library(tidyverse) # get for tibble
library(tidyverse) 
library(lubridate)
library(hash)
library(zoo) # for AUC calculation
library(glue)
library(gghighlight)
library("psych") # for geom mean calculation
library(gt)

# temporarily turn off warnings
options(warn=0)

# Import and data preprocessing -------------------------------------------

#set working directory
setwd('/Users/valer/Desktop/R_project/')
path_ <- getwd()
path_ <- "C:/Users/valer/Desktop/R_project/data.xlsx"

# loading data containing test and reference datasheets 
loading_data <- function(path_, sheet_) {
  
  # reading file
  data <- read_excel(path = path_, sheet = sheet_)
  
  # convert all types to numeric to enable further calculations 
  data <- data |> 
    mutate_all(as.numeric)
  sum(is.na(data))
  return (data)
  
}

# dataset for the test medicinal product
data_test <- loading_data(path_, "test_df")

# Calculation of linear regression parameters


# setting datastructues to store calculated data
# list of all r2 correlation coefficients
list_r2_test_product <- vector(mode = 'list', length = nrow(data_test))

# correlation coefficient
list_of_coeff_vs_r2_test_product <-
  vector(mode = 'list', length = ncol(data_test) - 1)

subject_list_data <- vector(mode = 'list', length = ncol(data_test))

# get the list of subjects to iterate over them later
subjects_list <- colnames(data_test)[2:ncol(data_test)]

# list of pk parameters
list_of_pharmacokinetic_parameters_test_product <-
  vector(mode = 'list', length = ncol(data_test) - 1)

# Looking for the best linear regression model 

for (subject in subjects_list) {
  
  data_work <- data_test[c('Time', subject)]
  
  window_start <- nrow(data_work[subject])
  
  # correction for the start
  Cmax <- max(data_work[subject], na.rm = TRUE)
  
  end <- which(data_work[[subject]] == Cmax) - 1 # need to substract one
  # to include the Cmax value itselfin the assessment
  
  # from start to iter start
  while (is.na(data_work[window_start,subject]) == TRUE){
    window_start <- window_start - 1
  }
  
  window_end <- window_start - 2 # at least 3 points need to be used
  # for the calculation of the linear model
  
  # to grow until it is equal to end (Cmax)
  r2max <- -1e5
  lm_coefficient_range <- NULL
  lm_r2_range <- NULL
  
  # iterate over each row of small data frame 
  for (row in (1: nrow(data_work))) {
    
    # immediately end iteration because iteration is done 
    # only to Cmax inclusive
    if (window_end == end) {
      break
    }
    
    # slice data frame to get y and x vectors for 
    # linear regression: x - time, y - concentration
    # linear regression fit starts with min of 3 points 
    # from the end of data frame
    # then iterate towards the beginning of data frame 
    # getting increasing x and y by 1
    
    # start with 3 points
    # time range
    x <- data_work[["Time"]][window_end:window_start] 
    
    # concentration range
    y <- data_work[[subject]][window_end:window_start] 
    # need to log transform y coordinates (concentration)
    y_log <- log(data_work[[subject]][window_end:window_start])
    
    # fit a model, get coefficient and r2 values
    model <- lm(formula = y_log ~ x, 
                data = data_work[window_end:window_start,])
    coefficient <- model$coefficients["x"]
    r2 <- summary(model)$r.squared
    
    # append coefficient to the coefficient range list
    lm_coefficient_range <- c(lm_coefficient_range, coefficient)
    # append r2 to the linear model r2 range list
    lm_r2_range <- c(lm_r2_range, r2)
    
    # use char as keys for list
    r2 <- as.character(r2)
    list_r2_test_product[[r2]] <- coefficient
    
    # increase the number of points for the current window of x and y
    # iteration from the end of data frame
    window_end <- window_end - 1
    
    # find the largest r2
    if (r2 > r2max) {
      r2max <- max(r2max, r2) # correct this code later
      
      # select x and y values that showed maximal r2
      xvalues_bestfit <- x
      yvalues_bestfit <- y
    }
    elimination_constant_optimized <- list_r2_test_product[[r2max]]
    
    output <- list(
      elimination_constant_optimized = elimination_constant_optimized,
      r2max = r2max,
      lm_coefficient_range = lm_coefficient_range,
      lm_r2_range = lm_r2_range,
      xvalues_bestfit = xvalues_bestfit,
      yvalues_bestfit = yvalues_bestfit,
      initial_df = data_work)
    
    subject_list_data[[subject]] <- output
  }
}


# this function computes 
pharmacokinetic_parameters_calculation <- function(data_temp) {
  
  # computing Cmax - max plasma concentration
  Cmax <- max(data_temp[subject])
  
  # computing the lowest concentration possible, it usually at the longest time
  # at the last row of data frame
  Clast <- data_temp[nrow(data_temp),][[subject]]
  
  # first timepoint (x coordinate) which showed the best fit that resulted in the
  # max eliminaiton constant
  Time_interval_start <- subject_list_data[[subject]]$xvalues_bestfit[1]
  
  # last timepoint (x coordinate) which showed the best fit that resulted in the
  # max eliminaiton constant
  Time_interval_end <- data_temp[nrow(data_temp),]["Time"]
  
  # Tmax - time at which Cmax was found
  Tmax <- which(data_temp[subject] == Cmax)
  
  # elimination constant was negative because of the negative trent of the 
  # elimination model. For furhter calculation it is needed to get absolute value
  elimination_constant_optimized <- abs(subject_list_data[[subject]]$elimination_constant_optimized)
  
  # elimination half-time, time at which 50 5 of drug is eliminated
  half_time <- round(log(2)/ elimination_constant_optimized, 2)
  
  # area under the curve (AUC) calculation
  # get coordinates for trapezoid
  x <- data_temp[["Time"]]
  y <- data_temp[[subject]]
  # order by index and get cumulative sume
  id <- order(x)
  AUC0_last <- sum(diff(x[id])*rollmean(y[id],2))
  
  # area under the curve at infinite time
  AUC0_inf <- AUC0_last + Clast/elimination_constant_optimized
  
  # residual area that is not covered (difference between AUC0_inf and AUC0_last)
  residual_area <- (AUC0_inf-AUC0_last)/AUC0_inf * 100
  
  # store calculated values for each subject
  list_of_pharmacokinetic_parameters_test_product[[subject]] <- 
    data.frame(Time_interval_start,
               Time_interval_end,
               elimination_constant_optimized,
               half_time, Clast, AUC0_last, 
               AUC0_inf, residual_area, Cmax, Tmax)
  
  return (list_of_pharmacokinetic_parameters_test_product)
}

# iterate over subjects
for (subject in subjects_list) {
  
  data_temp <- data_test |> 
    select(1, subject) |> 
    na.omit()
  list_of_pharmacokinetic_parameters_test_product <- pharmacokinetic_parameters_calculation(data_temp)
}

# to create dataframe
pharmac_param_ci_calculation <- 
  bind_rows(list_of_pharmacokinetic_parameters_test_product, .id = "Subject")

pk_df <- as.data.frame(pharmac_param_ci_calculation)
write_xlsx(pk_df, 'C:\\Users\\valer\\Desktop\\R\\data1.xlsx')



