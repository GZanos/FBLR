

set.seed(12311)

## Options Input ##
#-----------------------------------------------------------#
# fuzzy.var <- runif(8,1,10) # assumed variance for fuzzification of crisp values into GFN
#-----------------------------------------------------------#


## Abdalla & Buckley 2007, Data ##
#-----------------------------------------------------------#
X1 <- c(2.00, 0.00, 1.13, 2.00, 2.19, 0.25, 0.75, 4.25)
X2 <- c(0.00, 5.00, 1.50, 1.25, 3.75, 3.50, 5.25, 2.00)
X3 <- c(15.25, 14.13, 14.13, 13.63, 14.75, 13.75, 15.25, 13.50)


Y <- matrix(c(2.27, 5.83, 9.39, 0.33, 0.85, 1.37, 5.43, 13.93, 22.43, 
              1.56, 4.00, 6.44, 0.64, 1.65, 2.66, 0.62, 1.58, 2.54,
              3.19, 8.18, 13.17, 0.72, 1.85, 2.98), nrow = 8, ncol = 3, byrow = TRUE)
#GC: Abdalla use 39% off vertex to fuzzify. So as per 21/11 meeting, we just take the vertex as the crisp value
Y_crisp <- Y[,2] 
#-----------------------------------------------------------#
fuzzy.var <-(Y[,3]-Y[,2])*(Y[,2]-Y[,1]) # Bhatia–Davis inequality     (Y[,3]-Y[,1])^2/4 # HD: Take the variance for fuzification from the range using Popoviciu's inequality: V(X)<= (M-m)^2/4, M = sup X, m = inf X. https://en.wikipedia.org/wiki/Popoviciu%27s_inequality_on_variances

## Required GFN Operations ##
#-----------------------------------------------------------#
# Addition 
GFN.add <- function(A, B) {
  mean <- A[1]+B[1] 
  variance <- A[2]+B[2] 
  return(c(mean, variance))
}

# Subtraction 
GFN.sub <- function(A, B) {
  mean <- A[1]-B[1]
  variance <- A[2]+B[2] 
  return(c(mean, variance))
}

# Multiplication 
GFN.multi <- function(A, B) {
  mean <- A[1]*B[1]
  variance <- (B[2]*A[1]^2)+(A[2]*B[1]^2)+(B[2]*A[2])
  return(c(mean, variance))
}

# Division 
GFN.div <- function(A, B) {
  mean <- A[1]*((1/B[1])+(B[2]/B[1]^3))
  variance <- (A[1]^2*(1/B[1]^4)*B[2])+((1/B[1]^2)*A[2])-(A[2]*(1/B[1]^4)*B[2])
  return(c(mean, variance))
}
#-----------------------------------------------------------#


## Function to fuzzify X into GFN ##
#-----------------------------------------------------------#
# fuzzify <- function(X, variance) {
#   matrix(c(X, rep(variance, length(X))), ncol = 2, 
#          dimnames = list(NULL, c("Mean", "Variance")))
# }

fuzzify <- function(X, variance) {
  matrix(c(X, variance), ncol = 2, 
         dimnames = list(NULL, c("Mean", "Variance")))
}

GFN_X1 <- fuzzify(X1, fuzzy.var) #GC: just assume a small variance here for Xi
GFN_X2 <- fuzzify(X2, fuzzy.var)
GFN_X3 <- fuzzify(X3, fuzzy.var)
GFN_Y <- fuzzify(Y_crisp, fuzzy.var)
#-----------------------------------------------------------#



## Error functions ##
#-----------------------------------------------------------#
#GC: these should be the same as in Dr. Duygu's paper

# MSE
MSEe_GFN <- function(Y_true, Y_pred) {
  n <- nrow(Y_true)
  mse <- c(0, 0) 
  
  for (i in 1:n) {
    diff <- GFN.sub(Y_pred[i, ], Y_true[i, ])
    squared_diff <- GFN.multi(diff, diff)    
    mse <- GFN.add(mse, squared_diff)       
  }
  
  mse <- GFN.div(mse, c(n, 0.001)) 
  return(mse)
}


# MPE
MPEe_GFN <- function(Y_true, Y_pred) {
  n <- nrow(Y_true)
  mpe <- c(0, 0) 
  
  for (i in 1:n) {
    diff <- GFN.sub(Y_pred[i, ], Y_true[i, ])
    perc_diff <- GFN.div(diff, Y_true[i, ])
    mpe <- GFN.add(mpe, perc_diff)
  }

  mpe <- GFN.div(mpe, c(n, 0.001))
  return(mpe)
}


# MAPE
MAPEe_GFN <- function(Y_true, Y_pred) {
  n <- nrow(Y_true)
  mape <- c(0, 0) 
  
  for (i in 1:n) {
    abs_diff <- GFN.sub(Y_pred[i, ], Y_true[i, ])
    abs_diff <- c(abs(abs_diff[1]), abs(abs_diff[2])) 
    abs_perc_diff <- GFN.div(abs_diff, abs(Y_true[i, ]))
    mape <- GFN.add(mape, abs_perc_diff)
  }

  mape <- GFN.div(mape, c(n, 0))
  mape <- GFN.multi(mape, c(100, 0.001)) 
  return(abs(mape))
}


# SMAPE
SMAPEe_GFN <- function(Y_true, Y_pred) {
  n <- nrow(Y_true)
  smape <- c(0, 0) 
  
  for (i in 1:n) {
    abs_diff <- GFN.sub(Y_pred[i, ], Y_true[i, ])
    abs_diff <- c(abs(abs_diff[1]), abs(abs_diff[2])) 
    mean_sum <- GFN.add(Y_pred[i, ], Y_true[i, ])   
    mean_avg <- GFN.div(mean_sum, c(2, 0))        
    sym_perc_diff <- GFN.div(abs_diff, mean_avg)    
    smape <- GFN.add(smape, sym_perc_diff)
  }
  
  smape <- GFN.div(smape, c(n, 0.001)) 
  return(abs(smape))
}


# MAE
MAEe_GFN <- function(Y_true, Y_pred) {
  n <- nrow(Y_true)
  mae <- c(0, 0) 
  
  for (i in 1:n) {
    abs_diff <- GFN.sub(Y_pred[i, ], Y_true[i, ]) 
    abs_diff <- c(abs(abs_diff[1]), abs(abs_diff[2])) 
    mae <- GFN.add(mae, abs_diff) 
  }
  
  mae <- GFN.div(mae, c(n, 0.001)) 
  return(mae)
}
#-----------------------------------------------------------#


## Added HD: load the TFN coefficients ##
#----------------------------------------------------------------------------

# Coef values from Icen & Demirhan paper
#A0 <- c(0.061, 0.316, 0.341)
#A1 <- c(-0.271, -0.268, -0.129)
#A2 <- c(-0.822, -0.727, -0.721)
#A3 <- c(0.259, 0.294, 0.336)

# Coef values from Icen & Gunay paper
A0 <- c(-1.3474, -1.1216, -0.8193)
A1 <- c(-0.6321, -0.6308, -0.6295)
A2 <- c(-1.5218, -1.5198, -1.5149)
A3 <- c(0.6687, 0.6714, 0.6726)

# Coef values from Abdalla & Buckley paper
#A0 <- c(-0.710, -0.539, -0.524)
#A1 <- c(-0.610, -0.473, -0.472)
#A2 <- c(-1.090, -1.089, -1.088)
#A3 <- c(0.459, 0.487, 0.680)

# HD: Directly use the optimised TFNs through Popoviciu's inequality to get the corresponding GFNs. No other optimisation required.

Beta_0_sigma2 <- (A0[3]-A0[2])*(A0[2]-A0[1]) # Bhatia–Davis inequality  #  ((A0[3]-A0[1])^2)/4
Beta_1_sigma2 <- (A1[3]-A1[2])*(A1[2]-A1[1]) #((A1[3]-A1[1])^2)/4
Beta_2_sigma2 <- (A2[3]-A2[2])*(A2[2]-A2[1]) #((A2[3]-A2[1])^2)/4
Beta_3_sigma2 <- (A3[3]-A3[2])*(A3[2]-A3[1]) #((A3[3]-A3[1])^2)/4


Beta_0_mode <- A0[2]
Beta_1_mode <- A1[2]
Beta_2_mode <- A2[2]
Beta_3_mode <- A3[2]

Beta <- as.matrix(data.frame(Mean = c(Beta_0_mode,Beta_1_mode,Beta_2_mode,Beta_3_mode), Variance = c(Beta_0_sigma2,Beta_1_sigma2,Beta_2_sigma2,Beta_3_sigma2)))

Beta_0 <- Beta[1, ] # Intercept
Beta_1 <- Beta[2, ] # Coefficient for X1
Beta_2 <- Beta[3, ] # Coefficient for X2
Beta_3 <- Beta[4, ] # Coefficient for X3

n <- nrow(GFN_X1)
Pred_Y <- matrix(0, nrow = n, ncol = 2)
# HD: Calculate predictions:
for (i in 1:n) {
  # Linear combination
  term1 <- GFN.multi(Beta_1, GFN_X1[i, ])
  #HD: Check if the resulting GFNs are symmetric. Otherwise, we cannot do GFN addition!
  delta_Beta_1 <- abs(Beta_1[1]/sqrt(Beta_1[2]))
  delta_GFN_X1 <- abs(GFN_X1[i,1]/ sqrt(GFN_X1[i,2]))
  if ((delta_Beta_1 < 3) & (delta_GFN_X1 < 3)){
    cat("Term1",i,"\n")
  }
  
  term2 <- GFN.multi(Beta_2, GFN_X2[i, ])
  #HD: Check if the resulting GFNs are symmetric. Otherwise, we cannot do GFN addition!
  delta_Beta_2 <- abs(Beta_2[1]/sqrt(Beta_2[2]))
  delta_GFN_X2 <- abs(GFN_X2[i,1]/ sqrt(GFN_X2[i,2]))
  if ((delta_Beta_2 < 3) & (delta_GFN_X2 < 3)){
    cat("Term2",i,"\n")
  }
  
  term3 <- GFN.multi(Beta_3, GFN_X3[i, ])
  #HD: Check if the resulting GFNs are symmetric. Otherwise, we cannot do GFN addition!
  delta_Beta_3 <- abs(Beta_3[1]/sqrt(Beta_3[2]))
  delta_GFN_X3 <- abs(GFN_X3[i,1]/ sqrt(GFN_X3[i,2]))
  if ((delta_Beta_3 < 3) & (delta_GFN_X3 < 3)){
    cat("Term3",i,"\n")
  }
  
  sum_terms <- GFN.add(Beta_0, GFN.add(term1, GFN.add(term2, term3)))
  Pred_Y[i, ] <- sum_terms
}

#HD: Check if the resulting prediction GFNs are symmetric. Otherwise, we cannot take the mean/mode directly!

smallDelta <- which(Pred_Y[,1]/sqrt(Pred_Y[,2]) < 0.4) # 0.4 is arbitrary depending on the variances




predicted1 <- Pred_Y[,1] # HD: just use mean
predicted2 <- Pred_Y[,1] + .1*Pred_Y[,2] # HD since all deltas are small, just shift the mean based on a fraction of variance. We don't know which direction to shift.
predicted3 <- predicted1
predicted3[smallDelta] <- Pred_Y[smallDelta,1] + .1*Pred_Y[smallDelta,2] # HD: Use mean for those with relatively large deltas and shift the mean based on a fraction of variance for those with small deltas.

mean(abs(Y_crisp - predicted1))
mean(abs(Y_crisp - predicted2))
mean(abs(Y_crisp - predicted3))

sqrt(mean((Y_crisp - predicted1)^2))
sqrt(mean((Y_crisp - predicted2)^2))
sqrt(mean((Y_crisp - predicted3)^2))

mean(abs(Y_crisp - predicted1)/Y_crisp)
mean(abs(Y_crisp - predicted2)/Y_crisp)
mean(abs(Y_crisp - predicted3)/Y_crisp)

predictedTFN <- c(6.64, 1.62, 4.85, 4.40, 2.18, 2.92, 1.57, 2.35)
mean(abs(Y_crisp - predictedTFN))
sqrt(mean((Y_crisp - predictedTFN)^2))
mean(abs(Y_crisp - predictedTFN)/Y_crisp)

#-----------------------------------------------------------#
#-----------------------------------------------------------#


#GC5: what if we do something like this to optimize that 0.1 (called 'm' here)? It runs super quick and gives better results.
#-----------------------------------------------------------#

# Define the range of m to test
m_range <- seq(0, 1, by = 0.01)


optimal_coeff <- NULL
min_error <- Inf

# Optimize m
for (m in m_range) {
  predicted1 <- Pred_Y[, 1] # Just use mean
  predicted2 <- Pred_Y[, 1] + m * Pred_Y[, 2] # Shift mean based on fraction of variance
  predicted3 <- predicted1
  predicted3[smallDelta] <- Pred_Y[smallDelta, 1] + m * Pred_Y[smallDelta, 2]
  
  # Calculate the errors
  error1 <- mean(abs(Y_crisp - predicted1))
  error2 <- mean(abs(Y_crisp - predicted2))
  error3 <- mean(abs(Y_crisp - predicted3))
  
  # Sum errors to find optimal m
  total_error <- error1 + error2 + error3
  
  if (total_error < min_error) {
    min_error <- total_error
    optimal_m<- m
  }
}

# Use the optimal m
cat("Optimal m:", optimal_m, "\n")
predicted1 <- Pred_Y[, 1] # Just use mean
predicted2 <- Pred_Y[, 1] + optimal_m * Pred_Y[, 2] # Shift mean based on fraction of variance
predicted3 <- predicted1
predicted3[smallDelta] <- Pred_Y[smallDelta, 1] + optimal_m * Pred_Y[smallDelta, 2]

# Print predictions
cat("Predicted1:\n", predicted1, "\n")
cat("Predicted2:\n", predicted2, "\n")
cat("Predicted3:\n", predicted3, "\n")
#-----------------------------------------------------------#

