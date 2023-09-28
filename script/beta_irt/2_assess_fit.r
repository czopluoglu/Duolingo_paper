
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
                   omega  = rep(c(4.3,3.6,4.1,2.6,2.7,4.3,3.6,4.1,2.6,2.7),5))
################################################################################
# Import the datasets 

load("./output/beta.RData")

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
  xlim(0.5,3)+
  ylim(0.5,3)

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

# Extract the omega parameters

omega <- as.data.frame(summary(stanfit, 
                               pars = c("omega"), 
                               probs = c(0.025, 0.975))$summary)

omega$true_omega <- ipar$omega

round(omega,3)

ggplot(omega, aes(x = true_omega, y = mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Scatterplot of True Values vs. Mean Estimates") +
  xlab("True Values") +
  ylab("Mean Estimates")+
  xlim(1,5)+
  ylim(1,5)

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

Yhat_train <- posteriors$Yhat_train

# Posterior Distribution of Model Generated Responses for Training Data

Yhat_test <- posteriors$Yhat_test

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

for(item in 1:50){ # iterate over items
  
  loc1        <- which(d_test$item==item)
  obs         <- d_test$resp[loc1]
  obs_[item]  <- mean(obs)
  
  pred_[item] <- mean(rowMeans(Yhat_test[,loc1]))
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
  
  loc1        <- which(d_test$item==item)
  obs         <- d_test$resp[loc1]
  obs_[item]  <- sd(obs)
  
  Yhat_test[,loc1]
  
  # Var[E(Y|params)], variance of expected grade conditional on model parameters
  
  comp1 <- var(colMeans(Yhat_test[,loc1]))     
  
  # E[Var(Y|params)], expected value of variance of grades conditional on model parameters
  
  comp2 <- mean(apply(Yhat_test[,loc1],1,var)) 
  
  
  pred_[item] <- sqrt(comp1 + comp2)
  
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
  
  tmp <- rowMeans(Yhat_train[,loc1])
  
  score_posterior[p] <- sample(tmp,1)
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

score <- data.frame(type=c(rep('Observed',1000),rep('Posterior',1000)),
                    score = c(obs_score,score_posterior))

ggplot(score, aes(x=score, linetype=type)) + 
  geom_histogram(data=subset(score, type=="Observed"), 
                 aes(y=..density..),fill='white',color='black') + 
  geom_density(linewidth=0.1) + 
  theme_minimal() +
  labs(title="", 
       x="Total Score", y="Density") +  
  theme(legend.title = element_blank(),
        legend.key = element_blank())
