# FBLR
Fuzzy Bayesian Logistic Regression using GFN Classification.

We work with crisp input data (Xij predictors and Yi binary response variable) and implement a RBLR JAGS model to extract posterior distributions. 
The information from these posteriors are then used to create GFN coefficients and compute GFN probabilities using fuzzy arithmetic.
The GFN probabilities are compared against a GFN threshold by use of z-score tests and binary output predictions are formed.
