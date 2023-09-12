################################################################################
require(cmdstanr)
require(rstan)
require(psych)
require(loo)
require(purrr)
require(ggplot2)
################################################################################
# Import the datasets 
  
  d_wide <- read.csv('./data/sb_wide.csv')
  d      <- read.csv('./data/sb_long.csv')

################################################################################
# Train/Test split

  set.seed(692023)

  # Randomly select 10% of data for test
  
  loc <- sample(1:nrow(d),nrow(d)*.10)
  
  d_test  <- d[loc,]
  d_train <- d[-loc,]

  # Although it is less likely, make sure that training data has observations
  # from all test takers

    sum(sort(unique(d_train$id))==sort(unique(d$id)))
    
      # should yield the number of unique test takers in the dataset

################################################################################
# data object for the Stan model 
# check the file "beta.stan" for the model syntax

  data_rt <- list(I              = length(unique(d_train$item)),
                  P              = length(unique(d_train$id)),
                  n_obs          = nrow(d_train),
                  id             = d_train$id,
                  item           = d_train$item,
                  W              = d_wide$w,
                  S              = d_wide$s,
                  r              = cor(d_wide[,'s'],d_wide[,'w']),
                  Y              = d_train$resp,
                  n_obs2         = nrow(d_test),
                  Y2             = d_test$resp,
                  item2          = d_test$item,
                  id2            = d_test$id)

################################################################################

# Variables to define the number of items and number of observations
    
  nit <- length(unique(d$item))
  N   <- length(unique(d$id))

# Starting Values (not necessary but may help under extreme situations)

  start <- list(ln_alpha=rnorm(nit,0,.2),
                ln_delta=rnorm(nit,0,.2),
                beta    =rnorm(nit,0,.5),
                theta   =rnorm(N),
                gamma0  =-abs(rnorm(nit,1,.2)),
                gamma1  = abs(rnorm(nit,1,.2)),
                a = 0.25,
                b = 0.25)

# Compile the model

  mod  <- cmdstan_model('./script/sb_irt/sb.stan')

# Fit the model

  fit <- mod$sample(data            = data_rt,
                    seed            = 1234,
                    chains          = 4,
                    parallel_chains = 4,
                    iter_warmup     = 250,
                    iter_sampling   = 750,
                    refresh         = 100,
                    init            = list(start,start,start,start))

# Save the output object

  fit$cmdstan_summary()
  
  stanfit <- rstan::read_stan_csv(fit$output_files())
  
  rm(list=ls()[!ls()%in%c("stanfit","d_wide","d","d_train","d_test")])
  
  save.image("./output/sb.RData")
  
  