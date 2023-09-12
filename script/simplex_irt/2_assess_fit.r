require(cmdstanr)
require(rstan)
require(psych)
require(loo)
require(purrr)
require(ggplot2)
################################################################################
# True Item Parameters

ipar <- data.frame(gamma0 = rep(c(-2,-.9,-2.6,-.8,-3,-2.5,-2.9,-2.6,-1.8,2),5),
                   gamma1 = rep(c(0.5,1.9,1.6,1,-1,-0.5,.9,.6,.8,3),5),
                   beta   = rep(c(2.5,2.2,2.6,1.5,1.2,2.5,2.2,2.6,1.5,1.2),5),
                   alpha  = rep(c(0.9,0.9,0.7,0.6,1.2,0.9,0.9,0.7,0.6,0.7),5),
                   phi    = rep(c(6.2,5.5,6.2,4.8,4.3,5.3,6.9,6.2,5.5,6.2),5))

################################################################################
# Import the datasets 

load("./output/simplex.RData")

################################################################################

# Extract the beta parameters

beta <- as.data.frame(summary(stanfit, 
                              pars = c("beta"), 
                              probs = c(0.025, 0.975))$summary)
beta$true_beta <- ipar$beta
round(beta,3)

ggplot(beta, aes(x = true_beta, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(0,3)+
  ylim(0,3)

# Extract the alpha parameters

alpha <- as.data.frame(summary(stanfit, 
                               pars = c("alpha"), 
                               probs = c(0.025, 0.975))$summary)

alpha$true_alpha <- ipar$alpha

round(alpha,3)

ggplot(alpha, aes(x = true_alpha, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(0,2)+
  ylim(0,2)

# Extract the phi parameters

phi <- as.data.frame(summary(stanfit, 
                             pars = c("phi"), 
                             probs = c(0.025, 0.975))$summary)

phi$true_phi <- ipar$phi

round(phi,3)

ggplot(phi, aes(x = true_phi, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(2,10)+
  ylim(2,10)

# Extract the gamma0 parameters

gamma0 <- as.data.frame(summary(stanfit, 
                                pars = c("gamma0"), 
                                probs = c(0.025, 0.975))$summary)

gamma0$true_gamma0 <- ipar$gamma0
round(gamma0,3)

ggplot(gamma0, aes(x = true_gamma0, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(-4,3)+
  ylim(-4,3)

# Extract the gamma1 parameters

gamma1 <- as.data.frame(summary(stanfit, 
                                pars = c("gamma1"), 
                                probs = c(0.025, 0.975))$summary)

gamma1$true_gamma1 <- ipar$gamma1

round(gamma1,3)

ggplot(gamma1, aes(x = true_gamma1, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(-1.5,3.5)+
  ylim(-1.5,3.5)

# Extract the theta parameters

theta <- as.data.frame(summary(stanfit, 
                               pars = c("theta"), 
                               probs = c(0.025, 0.975))$summary)

theta$true_theta <- d_wide$true_theta
theta

cor(theta$mean,theta$true_theta)

ggplot(theta, aes(x = true_theta, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(-5,4)+
  ylim(-5,4)

describe(theta$mean)
describe(theta$true_theta)

# Extract the standard deviation of theta estimates

summary(stanfit, pars = c("theta_sd"), probs = c(0.025, 0.975))$summary

# Extract the regression coefficient estimates 

# True values were 0.3 and 0.4

coefs2 <- summary(stanfit, pars = c("a","b"), probs = c(0.025, 0.975))$summary

round(coefs2,3)

# Extract pointwise log-likelihood from the Training data

loglik_1 <- extract_log_lik(stanfit, 
                            parameter_name = "log_lik", 
                            merge_chains = TRUE)

loo_1 <- loo(loglik_1,cores = 4)
print(loo_1)

ll <- rowSums(loglik_1)
hist(ll)

log(mean(exp(ll)))

max(ll)
max(ll) +  log(mean(exp(ll - max(ll))))


# Extract pointwise log-likelihood from the Test data

loglik_2 <- extract_log_lik(stanfit, 
                            parameter_name = "log_lik2", 
                            merge_chains = TRUE)

loo_2 <- loo(loglik_2, cores = 4)
print(loo_2)

ll2 <- rowSums(loglik_2)
hist(ll2)

log(mean(exp(ll2)))

max(ll2)
max(ll2) +  log(mean(exp(ll2 - max(ll2))))


################################################################################
################################################################################
# Extract the posterior distributions

posteriors <- extract(stanfit)

# Posterior Distribution of Model Generated Responses for Training Data
#  Yhat_train <- posteriors$Yhat_train

# Posterior Distribution of Model Generated Responses for Training Data
#  Yhat_test <- posteriors$Yhat_test
################################################################################
################################################################################

# Stan doesn't have a random number generator for the Simplex function
# so, we can't generate model-based responses under the generated quantities block

# We have to do this outside of Stan using the posterior draw of model parameters
# something we can extract from model fit.

# The availability of simplexreg pacage allow to generate data based on posterior
# draw of model parameters.

require(simplexreg)

# Modify some internal functions in the simplexreg package
# to avoid computational issues

rIG <- function(epsilon, Tau) {
  z <- rchisq(1, 1)
  ss <- sqrt(4 * epsilon * z / Tau + (epsilon * z)^2)
  z1 <- epsilon + (epsilon^2) * Tau * z / 2 - (epsilon * Tau / 2) * ss
  u1 <- runif(1, 0, 1)
  xxx <- z1
  if(u1 > (epsilon / (epsilon + z1)) ) {
    xxx <- (epsilon^2) / z1
  }
  return(as.numeric(xxx))      
}

rMIG <- function(epsilon, Tau, mu) {
  x1 <- rIG(epsilon, Tau)
  x2 <- rchisq(1, 1)     
  x3 <- x2 * Tau * (epsilon^2)
  u2 <- runif(1, 0, 1)
  xx <- x1
  if(u2 < mu) { xx <- x1 + x3 }
  return(as.numeric(xx))
}


rsimplex_ <- function (n, mu, sig) {
  Sigma <- sig^2
  epsilon <- mu/(1 - mu)
  Tau <- Sigma * ((1 - mu)^2)
  yy <- rep(0, n)
  for (i in 1:n) {
    x <- rMIG(epsilon[i], Tau[i], mu[i])
    yy[i] <- x/(1 + x)
  }
  return(as.vector(yy))
}


# This was run with 4 chains and 750 sampling iterations on each chain 
# (after 250 warmup iterations)
# So, there are 3000 draws for each parameter 

# Posterior samples for item parameters (3000 x 10)
# Each column represents the item, each row represents a draw

beta_posterior   <- posteriors$beta
alpha_posterior  <- posteriors$alpha
phi_posterior    <- posteriors$phi
gamma0_posterior <- posteriors$gamma0
gamma1_posterior <- posteriors$gamma1

# Posterior samples for person parameters (3000 x 250)
# Each column represents the item, each row represents a draw

theta_posterior <- posteriors$theta

# Vectorized item and person indices

p   <- d_test$id        
i   <- d_test$item              

# Matrix to save simulated values

Yhat_test <- matrix(nrow=3000,ncol=nrow(d_test))

# Progress bar

pb <- txtProgressBar(min = 0, max = 3000, style = 3)

# Iterate over 3000 draws

for(iter in 1:3000){
  
  # Probability of 0s and 1s (vectorized)
  # given the person and item parameters at the corresponding draw
  
  q0 <- 1/(1 + exp(-(gamma0_posterior[iter,i]-alpha_posterior[iter,i]*theta_posterior[iter,p])))
  q1 <- 1/(1 + exp(-(gamma1_posterior[iter,i]-alpha_posterior[iter,i]*theta_posterior[iter,p])))
  
  # simplex distribution parameters in case the response is neither 0 nor 1
  # (vectorized)
  
  mu_  <- 1/(1+exp(-(beta_posterior[iter,i]+alpha_posterior[iter,i]*theta_posterior[iter,p])))
  sig_ <- sqrt(phi_posterior[iter,i])
  
  # Simulate the continuous response (vectrorized) 
  
  y_   <- rsimplex_(n = length(mu_),mu = mu_, sig = sig_)
  
  # Sample 0,1, or 2
  # If it is 2 replace it with the continous response
  
  # Some preparation for the inputs required by the map function
  
  prob    <- data.frame(p0 = q0, 
                        p1 = 1-q1, 
                        p2 = q1-q0,
                        y = y_)
  
  prob$id <- 1:nrow(prob)
  by_id   <- split(prob,prob$id)
  
  # map is iterating over all observations (much faster then for loop)
  
  C <- map(.x = by_id,
           .f = function(x){
             tmp <- sample(0:2, 1, prob = c(x$p0,x$p1,x$p2))
             tmp <- ifelse(tmp==2,x$y,tmp)
           })
  
  # Save the simulated values at the correspond draw
  Yhat_test[iter,] <- unlist(C)
  
  setTxtProgressBar(pb, iter)
}
################################################################################

p   <- d_train$id        
i   <- d_train$item              

# Matrix to save simulated values

Yhat_train <- matrix(nrow=3000,ncol=nrow(d_train))

# Progress bar

pb <- txtProgressBar(min = 0, max = 3000, style = 3)

# Iterate over 3000 draws

for(iter in 1:3000){
  
  # Probability of 0s and 1s (vectorized)
  # given the person and item parameters at the corresponding draw
  
  q0 <- 1/(1 + exp(-(gamma0_posterior[iter,i]-alpha_posterior[iter,i]*theta_posterior[iter,p])))
  q1 <- 1/(1 + exp(-(gamma1_posterior[iter,i]-alpha_posterior[iter,i]*theta_posterior[iter,p])))
  
  # simplex distribution parameters in case the response is neither 0 nor 1
  # (vectorized)
  
  mu_  <- 1/(1+exp(-(beta_posterior[iter,i]+alpha_posterior[iter,i]*theta_posterior[iter,p])))
  sig_ <- sqrt(phi_posterior[iter,i])
  
  # Simulate the continuous response (vectrorized) 
  
  y_   <- rsimplex_(n = length(mu_),mu = mu_, sig = sig_)
  
  # Sample 0,1, or 2
  # If it is 2 replace it with the continous response
  
  # Some preparation for the inputs required by the map function
  
  prob    <- data.frame(p0 = q0, 
                        p1 = 1-q1, 
                        p2 = q1-q0,
                        y = y_)
  
  prob$id <- 1:nrow(prob)
  by_id   <- split(prob,prob$id)
  
  # map is iterating over all observations (much faster then for loop)
  
  C <- map(.x = by_id,
           .f = function(x){
             tmp <- sample(0:2, 1, prob = c(x$p0,x$p1,x$p2))
             tmp <- ifelse(tmp==2,x$y,tmp)
           })
  
  # Save the simulated values at the correspond draw
  Yhat_train[iter,] <- unlist(C)
  
  setTxtProgressBar(pb, iter)
}

################################################################################

# MODEL FIT

pb <- txtProgressBar(min = 0, max = 3000, style = 3)

# SSE BASELINE VS TEST SET

# Randomly sample a value from an observed distribution of an item
# in the training data to predict observations in the test set

SSE_test_baseline <- matrix(nrow=3000,ncol=1000)

for(i in 1:3000){
  for(j in 1:50){
    ind  <- d_test[j,]$item
    loc  <- which(d_test$item==j)
    sub  <- d_train[d_train$item==ind,]$resp
    SSE_test_baseline[i,loc] = sample(sub,length(loc),replace=TRUE)
  }
  setTxtProgressBar(pb, i)
}

SSE_baseline <- rep(0,3000)

for(i in 1:3000){
  tmp             <-  SSE_test_baseline[i,]
  SSE_baseline[i] <- sum((d_test$resp - tmp)^2)
  setTxtProgressBar(pb, i)
}

hist(SSE_baseline)
mean(SSE_baseline)

# Use the posterior predictive distribution to predict the observations in the
# test set

SSE_test <- rep(0,3000)

for(i in 1:3000){
  tmp          <- Yhat_test[i,]
  SSE_test[i]  <- sum((d_test$resp - tmp)^2)
}

hist(SSE_test)
mean(SSE_test)

# Proportional reduction in SSE

(mean(SSE_baseline) - mean(SSE_test))/mean(SSE_baseline)


# Recovery of Average Item Scores

obs_  <- c()   # The observed average item score in the test data
pred_ <- c()   # The posterior mean of average item scores in the test data

for(item in 1:50){
  
  loc1  <- which(d_test$item==item)
  obs   <- d_test$resp[loc1]
  obs_[item] <- mean(obs)
  
  pred <- c()
  
  for(j in 1:3000){
    pred[j] <- mean(Yhat_test[j,loc1])
  }
  
  pred_[item] <- mean(pred)
  
}


ggplot() +
  geom_point(aes(y = obs_, x = pred_)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("") +
  ylab("Observed Average Item Score") +
  xlab("Posterior Mean of Model Predicted Average Item Score")+
  xlim(0,1)+
  ylim(0,1)

# Recovery of Standard deviation of item scores

obs_  <- c()   # The observed standard deviation in the test data
pred_ <- c()   # The posterior mean of standard deviation in the test data

for(item in 1:50){
  
  loc1  <- which(d_test$item==item)
  obs   <- d_test$resp[loc1]
  obs_[item] <- sd(obs)
  
  pred <- c()
  
  for(j in 1:3000){
    pred[j] <- sd(Yhat_test[j,loc1])
  }
  
  pred_[item] <- mean(pred)
  
}


ggplot() +
  geom_point(aes(y = obs_, x = pred_)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("") +
  ylab("Observed Standard Deviation of Item Scores") +
  xlab("Posterior Mean of Model Predicted Standard Deviation of Item Scores")+
  xlim(0,1)+
  ylim(0,1)


# Total Score Distribution

d_train_wide <- reshape(data = d_train,
                        idvar = c('id','w','s','true_theta'),
                        timevar = 'item',
                        direction= 'wide')

obs_score          <- rowMeans(d_train_wide[,5:54],na.rm=TRUE)

score_posterior    <- c()

for(p in 1:1000){
  
  loc1  <- which(d_train$id==d_train_wide$id[p])
  
  tmp <- c()
  
  for(j in 1:3000){
    tmp[j] <- mean(Yhat_train[j,loc1])
  }
  
  score_posterior[p] <- mean(tmp)
  #print(p)
}

cor(obs_score,score_posterior)

ggplot() +
  geom_point(aes(y = obs_score, x = score_posterior)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("") +
  ylab("Observed Average Person Score") +
  xlab("Posterior Mean of Model Predicted Average Person Score")+
  xlim(0,1)+
  ylim(0,1)

score <- data.frame(type=c(rep('Observed',250),rep('Posterior',250)),
                    score = c(obs_score,score_posterior))

ggplot(score, aes(x=score, linetype=type)) + 
  geom_histogram(data=subset(score, type=="Observed"), 
                 aes(y=after_stat(density)),fill='white',color='black') + 
  geom_density(linewidth=0.1) + 
  theme_minimal() +
  labs(title="", 
       x="Total Score", y="Density") +  
  theme(legend.title = element_blank(),
        legend.key = element_blank())
