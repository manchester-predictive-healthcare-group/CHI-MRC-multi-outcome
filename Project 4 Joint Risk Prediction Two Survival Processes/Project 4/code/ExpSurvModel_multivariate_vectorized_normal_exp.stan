data {
    int<lower=0> N; //number of patients
    int<lower=1> N_cens_A; //number of patients
    int<lower=1> N_uncens_A; //number of patients
    int<lower=1> N_cens_B; //number of patients
    int<lower=1> N_uncens_B; //number of patients
    
    int<lower=0> P; //number of parameters
    matrix[N_cens_A, P] X_cens_A; //design matrix pulled from the stacked dataset
    matrix[N_uncens_A, P] X_uncens_A; //design matrix pulled from the stacked dataset
    matrix[N_cens_B, P] X_cens_B; //design matrix pulled from the stacked dataset
    matrix[N_uncens_B, P] X_uncens_B; //design matrix pulled from the stacked dataset
    
    int IDs_cens_A[N_cens_A]; //IDs in the stacked dataset saying which pt. the i-th row of the stacked data corresponds
    int IDs_uncens_A[N_uncens_A]; //IDs in the stacked dataset saying which pt. the i-th row of the stacked data corresponds
    int IDs_cens_B[N_cens_B]; //IDs in the stacked dataset saying which pt. the i-th row of the stacked data corresponds
    int IDs_uncens_B[N_uncens_B]; //IDs in the stacked dataset saying which pt. the i-th row of the stacked data corresponds
    
    vector<lower=0>[N_cens_A] times_cens_A; //observation times for the i-th row of the stacked dataset
    vector<lower=0>[N_uncens_A] times_uncens_A; //observation times for the i-th row of the stacked dataset
    vector<lower=0>[N_cens_B] times_cens_B; //observation times for the i-th row of the stacked dataset
    vector<lower=0>[N_uncens_B] times_uncens_B; //observation times for the i-th row of the stacked dataset
}

parameters {
    vector[P] betas_A;        
    vector[P] betas_B;        
    
    real intercept_A;
    real intercept_B;
    
    vector[N] omega;
    real<lower=0> frail_param;
}


model {
    betas_A ~ normal(0,1);          
    betas_B ~ normal(0,1);          
    
    intercept_A ~ normal(-7.5,2);   
    intercept_B ~ normal(-7.5,2);   
    
    omega ~ normal(0,frail_param);
    frail_param ~ gamma(1, 1);
    
    target += exponential_lpdf(times_uncens_A | exp(omega[IDs_uncens_A] + intercept_A + (X_uncens_A*betas_A)));
    target += exponential_lccdf(times_cens_A | exp(omega[IDs_cens_A] + intercept_A + (X_cens_A*betas_A)));
    target += exponential_lpdf(times_uncens_B | exp(omega[IDs_uncens_B] + intercept_B + (X_uncens_B*betas_B)));
    target += exponential_lccdf(times_cens_B | exp(omega[IDs_cens_B] + intercept_B + (X_cens_B*betas_B)));

}
