data {
    int<lower=0> N; //number of patients
    int<lower=N> stacked_N; //number of rows in the stacked data of all outcomes (e.g. N*2 for two outcome)
    
    int<lower=0> P; //number of parameters
    matrix[stacked_N, P] X; //design matrix pulled from the stacked dataset
    
    int IDs[stacked_N]; //IDs in the stacked dataset saying which pt. the i-th row of the stacked data corresponds
    
    vector<lower=0>[stacked_N] times; //observation times for the i-th row of the stacked dataset
    vector[stacked_N] status; //vector of status flags indicating if the i-th row of the stacked data corresponds to censoring (0) or event (1)
    int outcome[stacked_N]; //denotes which outcome the i-th row of stacked dataset corresponds to
}

parameters {
    vector[P] betas_A;        
    vector[P] betas_B;        
    
    real intercept_A;
    real intercept_B;
    
    real<lower=0> shape_A;
    real<lower=0> shape_B;

    vector[N] omega;
    real<lower=0> frail_param;
}


model {
    betas_A ~ normal(0,2);          
    betas_B ~ normal(0,2);          
    
    intercept_A ~ normal(-3,2);   
    intercept_B ~ normal(-3,2);     
        
    shape_A ~ normal(1,2);   
    shape_B ~ normal(1,2);   

    omega ~ normal(0,frail_param);
    frail_param ~ cauchy(0, 2.5);
    
    for (i in 1:stacked_N) {
      if (status[i] == 0){ //if censored
        if (outcome[i] == 1) {
          target += weibull_lccdf(times[i] | shape_A, exp(omega[IDs[i]] + intercept_A + (X[i,]*betas_A)));
        }else if(outcome[i] == 2) {
          target += weibull_lccdf(times[i] | shape_B, exp(omega[IDs[i]] + intercept_B + (X[i,]*betas_B)));  
        }
      }
      else if (status[i] == 1){ //if event/uncensored
        if (outcome[i] == 1) {
          target += weibull_lpdf(times[i] | shape_A, exp(omega[IDs[i]] + intercept_A + (X[i,]*betas_A)));
        }else if(outcome[i] == 2) {
          target += weibull_lpdf(times[i] | shape_B, exp(omega[IDs[i]] + intercept_B + (X[i,]*betas_B))); 
        }
      }
    }
}
