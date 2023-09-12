data{
  int  I;                 //  number of items
  int  P;                 //  number of individuals
  int  n_obs;             //  number of observed responses
  int  item[n_obs];       //  item id
  int  id[n_obs];         //  person id
  real W[P];              //  writing scores, standardized input
  real S[P];              //  speaking scores, standardized input
  real r;                 // correlation between writing and speaking
  real Y[n_obs];          //  vector outcome
  int  n_obs2;             //  number of observed responses in test set
  real Y2[n_obs2];         //  vector outcome for test set
  int  item2[n_obs2];       //  item id for test set
  int  id2[n_obs2];         //  person id for test set
}

parameters {
  real beta[I];         // vector of beta parameters for  I items
  real ln_alpha[I];     // vector of alpha parameters for I items     
  real ln_delta[I];     // vector of delta parameters for I items     
  real gamma0[I];       // vector of gamma0 parameters for I items     
  real gamma1[I];       // vector of gamma1 parameters for I items
  real a;               // regression slope for writing
  real b;               // regression slope for speaking
  real theta[P];        // vector of theta parameters for P individuals
}

transformed parameters{
  
  real<lower=0> delta[I];
  real<lower=0> alpha[I];
  
  for(i in 1:I){
       alpha[i] = exp(ln_alpha[i]);
       delta[i] = exp(ln_delta[i]);
  } 
  
  real<lower=0> theta_sd = sqrt(1 - (a*a+b*b+2*a*b*r));
}
   
model{
    
  beta      ~ normal(0,10);
  ln_alpha  ~ normal(0,10);
  ln_delta  ~ normal(0,10);
  gamma0    ~ normal(0,10);
  
  a         ~ normal(0,1);
  b         ~ normal(0,1);

  for(i in 1:I){
    gamma1[i] ~ normal(0,10)T[gamma0[i],];
  }
  
  for(i in 1:P){
    theta[i] ~ normal(a*W[i] + b*S[i], theta_sd);      
  }


  for(i in 1:n_obs) {
    
    real   k0 = inv_logit(gamma0[item[i]] - alpha[item[i]]*theta[id[i]]); 
    real   k1 = inv_logit(gamma1[item[i]] - alpha[item[i]]*theta[id[i]]);
  
    if (Y[i] == 0) {
      target += log(k0);
    } else if (Y[i] == 1) {
      target += log1m(k1);                       
    } else {
      real   y_star = logit(Y[i]);               
      real   p  = alpha[item[i]]*theta[id[i]] + beta[item[i]];
      real   q  = sqrt(delta[item[i]]);
      
      target += log(k1 - k0) + normal_lpdf(y_star | p, q);   
    }	

  }
}


generated quantities {
  
  vector[n_obs] log_lik;
  vector[n_obs2] log_lik2;
  
  vector[n_obs]  Yhat_train;
  vector[n_obs2] Yhat_test;
  
  for(i in 1:n_obs) {
    real   k0 = inv_logit(gamma0[item[i]] - alpha[item[i]]*theta[id[i]]); 
    real   k1 = inv_logit(gamma1[item[i]] - alpha[item[i]]*theta[id[i]]);
    
    real   y_star = logit(Y[i]);               
    real   p  = alpha[item[i]]*theta[id[i]] + beta[item[i]];
    real   q  = sqrt(delta[item[i]]);
    
    real probs_array[3] = {k0, 1-k1, k1-k0};
    vector[3] probs = to_vector(probs_array);
      
    if (Y[i] == 0) {
      log_lik[i]= log(k0);
    } else if (Y[i] == 1) {
      log_lik[i]= log1m(k1);                       
    } else {
      log_lik[i] = log(k1 - k0) + normal_lpdf(y_star | p, q);   
    }
    
    int tmp = categorical_rng(probs);
    if (tmp == 3) {
      Yhat_train[i] = inv_logit(normal_rng(p,q));
    } else {
      Yhat_train[i] = tmp - 1;
    }
  }
  
  
  for(i in 1:n_obs2) {
    real   k0 = inv_logit(gamma0[item2[i]] - alpha[item2[i]]*theta[id2[i]]); 
    real   k1 = inv_logit(gamma1[item2[i]] - alpha[item2[i]]*theta[id2[i]]);

    real   y_star2 = logit(Y2[i]);               
    real   p  = alpha[item2[i]]*theta[id2[i]] + beta[item2[i]];
    real   q  = sqrt(delta[item2[i]]);
    
    real probs_array[3] = {k0, 1-k1, k1-k0};
    vector[3] probs = to_vector(probs_array);
    
    if (Y2[i] == 0) {
      log_lik2[i]= log(k0);
    } else if (Y2[i] == 1) {
      log_lik2[i]= log1m(k1);                       
    } else {
      log_lik2[i] = log(k1 - k0) + normal_lpdf(y_star2 | p, q);   
    }
    
    int tmp = categorical_rng(probs);
    if (tmp == 3) {
      Yhat_test[i] = inv_logit(normal_rng(p,q));
    } else {
      Yhat_test[i] = tmp - 1;
    }
  }
}
