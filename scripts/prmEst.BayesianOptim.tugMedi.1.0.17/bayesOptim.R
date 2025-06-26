#!/bin/Rscript --vanilla --slave

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: thisR filename 1([std], 2[try], 3[rndm]) \'c(3,4)\'(col num for log10) T(nuggetEstim) > hoge.txt 2>&1\n  nuggetEstim F (adjust the value) when Variance < Nugget in km() with T")
}
input_file <- args[1]
MODE       <- args[2] # 1 | 2 | 3 (random search)
LOG10      <- eval( parse( text=args[3] ) )
NUGGET.EST <- as.logical( args[4] ) # T or F, neglected in random search


dd <- read.table( input_file, head=T )
################################################
#rep     ABCdist m0      dN      M2      M3
#1828    0.46838 1e-08   1e-05   75      89
#1772    0.47664 1e-08   1e-04   84      96
#678     0.47687 1e-08   1e-05   81      98
#866     0.49527 1e-08   1e-04   70      76
#1696    0.50953 1e-08   1e-05   62      105
################################################
head(dd)


# dd <- read.table( "clipboard", head=T ) 


# MODE  <- 1 
# LOG10 <- c(3,4) 
# NUGGET.EST <- T 


data <- dd


# log10
forLog10 <- LOG10
data[ forLog10 ] <- lapply(data[ forLog10 ], log10)
head(data)


data <- data[ , -c(1) ]
head(data)
# 1: evaluation score; smaller, better
# 2, ...: variables


# remove constants
data <- data[, sapply(data, function(x) length( unique(x) ) > 1)]
head(data)


# remove duplicates
data <- data[!duplicated(data), ]



# ============================================================================


library(dplyr)


## response column
unique_data <- data %>%
  select(response = colnames(data)[1], everything())
print(head(unique_data))


## scale
scaled_unique_data <- scale(unique_data)
print(head(scaled_unique_data ))
means <- attr(scaled_unique_data, "scaled:center")  # Later used
sds   <- attr(scaled_unique_data, "scaled:scale")   # Later used


applied_data <- as.data.frame(scaled_unique_data)
summary(applied_data)



# Gaussian Process Regression =============================================


library(DiceKriging)


sttTime <- Sys.time()
ranges <- apply( as.matrix(applied_data[, -1, drop = FALSE]) , 2, function(x) diff(range(x)) )

if        ( MODE == 1 | MODE == 3 ) {
  my.covtype <- "matern5_2" 
  my.lower   <- ranges / 10 
  my.upper   <- ranges * 10 
} else if ( MODE == 2 ) {
  my.covtype <- "matern3_2"   # "matern5_2" or "matern3_2" 
  my.lower   <- rep( 1/10, length(ranges) )  
  my.upper   <- rep(   10, length(ranges) )  
}
gp_model <- km(
  formula  = ~1,
  design   = applied_data[, -1, drop = FALSE],  
  response = applied_data$response, 
  
  optim.method = "BFGS", 
  nugget = 1e-2, 
  nugget.estim = NUGGET.EST, 
  
  covtype = my.covtype,
  lower   = my.lower, 
  upper   = my.upper  
)

print(gp_model)
Sys.time() - sttTime



# =========================================================================


library(rBayesianOptimization)


# Evaluation Function Modeled by Gaussian Process  -------------------------
evaluate_function <- function(...) {
  input_parameters <- as.data.frame(list(...))
  colnames(input_parameters) <- colnames(applied_data[, -1, drop = FALSE]) 

  prediction <- predict(
    gp_model,
    newdata = input_parameters,
    type = "UK"
  )
  
  list(Score = -1 * prediction$mean) # For minimum
}


predictions_train <- predict(
  gp_model,
  newdata = applied_data[, -1, drop = FALSE],
  type = "UK"
)
rmse_train <- sqrt(mean((applied_data$response - predictions_train$mean)^2))
cat("Training RMSE: ", rmse_train, "\n")


cat("Log-likelihood: ", gp_model@logLik, "\n")



# Bayesian Optimization ------------------------------------------------------


# Set bounds 
bounds <- setNames(
  lapply(  applied_data[, -1, drop = FALSE], function(x) range(x)), 
  colnames(applied_data[, -1, drop = FALSE]) 
)


sttTime <- Sys.time()
result <- list()
### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if        ( MODE == 1 ) {
  my.acq   <- "ucb" 
  my.kappa <- 2.576       # def, 2.576 | fast, 2.5-2.0 | slow, 1.5-0.8   | sloth, <0.5 
  my.eps   <- 0           # def, 0     | fast, 0       | slow, 0.01-0.05 | sloth, >0.1 
} else if ( MODE == 2 ) {
  my.acq   <- "ucb"       # "ucb" or "ei" 
  my.kappa <- 2.576       # for ucb
  my.eps   <- 1e-4        # for ei
}
if        ( MODE == 1 | MODE == 2 ) {
result <- BayesianOptimization(
  FUN = function(...) { evaluate_function(...) }, 
  bounds = bounds, 
  
  init_points = 20,    # 
  n_iter      = 75,    # 
  verbose     = F, 
  
  acq   = my.acq, 
  kappa = my.kappa, 
  eps   = my.eps
)
### --------------------------------------------------------------
} else if ( MODE == 3 ) {
#nTry <- 1e1               #
nTry <- 1e4               #

random_samples <- as.data.frame(
  lapply(bounds, function(b) runif(nTry, min = b[1], max = b[2]))
)
colnames(random_samples) <- colnames(applied_data[, -1, drop = FALSE])

scores <- apply(random_samples, 1, function(params) {
  params_list <- as.list(params)  
  evaluate_result <- do.call(evaluate_function, params_list)  
  evaluate_result$Score  
})

best_index        <- which.min( scores )
result$Best_Par   <- as.numeric( random_samples[best_index, ] )
result$Best_Value <-             scores[        best_index]
names(result$Best_Par) <- colnames(applied_data[, -1, drop = FALSE])
}
### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Sys.time() - sttTime


cat("Best prms: ", result$Best_Par, "\n")    
cat("Best scr: ",  result$Best_Value, "\n") 



# ===========================================================================


# Helper function to restore original scale
rescale <- function(values, means, sds) {
  values * sds + means
}
original_best_par   <- rescale(      result$Best_Par,   means[names(result$Best_Par)], sds[names(result$Best_Par)] )
original_best_value <- rescale( -1 * result$Best_Value, means["response"],             sds["response"] )


cat("Original-scaled best prms: ", original_best_par,   "\n")
cat("Original-scaled best scr:  ", original_best_value, "\n")



