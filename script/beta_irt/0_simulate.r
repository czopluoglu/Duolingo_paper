require(MASS)
require(ggplot2)
################################################################################
# Data Generation

set.seed(9132022)

  N = 1000 # sample size

# True Item parameters for 100 hypothetical items

  ipar <- data.frame(gamma0 = rep(c(-2,-.9,-2.6,-.8,-3,-2.5,-2.9,-2.6,-1.8,2),5),
                     gamma1 = rep(c(0.5,1.9,1.6,1,-1,-0.5,.9,.6,.8,3),5),
                     beta   = rep(c(2.5,2.2,2.6,1.5,1.2,2.5,2.2,2.6,1.5,1.2),5),
                     alpha  = rep(c(0.9,0.9,0.7,0.6,1.2,0.9,0.9,0.7,0.6,0.7),5),
                     omega  = rep(c(4.3,3.6,4.1,2.6,2.7,4.3,3.6,4.1,2.6,2.7),5))

# Generate hypothetical Writing and speaking scores (standardized)

  r = 0.5 # correlation between writing and speaking scores
  
  Sigma <- matrix(c(1,r,r,1),2,2)
  
  ws <- mvrnorm(N,c(0,0),Sigma)
  
  w <- (ws[,1]-mean(ws[,1]))/sd(ws[,1])
  s <- (ws[,2]-mean(ws[,2]))/sd(ws[,2])

# Hypotheticalregression weights for writing and speaking scores

  b1 <- 0.3
  b2 <- 0.4

  # Standard deviation of errors, so the variance of true theta is equal to 1
  
  theta_sd <- sqrt(1-(b1^2 + b2^2 + 2*b1*b2*cor(w,s)))

# Generate the true person parameters, 
# based on the latent regression, 
  
  theta <- c()
  
  for(i in 1:N){
    theta[i] <- rnorm(1,b1*w[i]+b2*s[i],theta_sd)
  }
  
  mean(theta)  # should be close to 0
  sd(theta)    # should be close to 1

# Generate Item Response Data

d <- as.data.frame(matrix(nrow=N,nco=nrow(ipar)))

for(p in 1:N){
  for(i in 1:nrow(ipar)){
    
    q0 <- plogis(ipar$gamma0[i]-ipar$alpha[i]*theta[p])
    q1 <- plogis(ipar$gamma1[i]-ipar$alpha[i]*theta[p])
    
    p0 <- q0
    p1 <- 1-q1
    p2 <- q1-q0
    
    C=sample(0:2, 1, prob = c(p0,p1,p2))
    
    if(C==2){
      
      sh1  <- exp((ipar$alpha[i]*theta[p]+ipar$beta[i]+ipar$omega[i])/2)     # Eq4, p5
      sh2  <- exp((-(ipar$alpha[i]*theta[p]+ipar$beta[i])+ipar$omega[i])/2)  # Eq5, p5
      
      d[p,i] <- rbeta(1,sh1,sh2)
      
    } else {
      d[p,i] = C
    }
  }
}


# Remove 80% of data so each person responds only 10 items out of 50

for(i in 1:nrow(d)){
  tmp <- sample(1:50,10)
  d[i,-tmp] = NA
}

# Rehape it to Long format

d$id   <- 1:nrow(d)
d$w    <- w
d$s    <- s
d$true_theta <- theta

d_long <- reshape(data = d,
                  idvar = 'id',
                  varying = colnames(d)[1:nrow(ipar)],
                  v.names = 'resp',
                  timevar = 'item',
                  direction= 'long')

d_long <- na.omit(d_long)

i = 4

ggplot(d_long[d_long$item==i,], aes(x=resp)) + 
  geom_histogram() + 
  geom_density(linewidth=0.1) + 
  theme_minimal()


write.csv(d_long,'./data/beta_long.csv',row.names=FALSE)
write.csv(d,'./data/beta_wide.csv',row.names=FALSE)



