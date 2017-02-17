data{
    int N;
    int N_pops;
    int<lower=0,upper=1> y[N]; // outcome
    int pop_id[N];
    real age[N];
    vector[5] age_meta[N]; // adult,min,max,mean,sd
    int N_missing_ages;
    real age_miss_min[N_missing_ages];
    real age_miss_max[N_missing_ages];
    int miss_idx[N];
}
transformed data{
    vector[2] zeros;
    for ( i in 1:2 ) zeros[i] = 0;
}
parameters{
    real a; // overall intercept
    real b_age; // average age slope
    vector[2] v_pop[N_pops]; // intercept and slope on age for each population
    vector<lower=0>[2] sigma_pop; // scale parameters for population effects
    corr_matrix[2] Rho_pop; // correlation matrix for population effects
    vector<lower=0,upper=1>[N_missing_ages] age_impute_raw;
}
transformed parameters{
    vector[N_missing_ages] age_impute;
    for ( i in 1:N_missing_ages ) {
        age_impute[i] = age_impute_raw[i]*(age_miss_max[i]-age_miss_min[i]) + age_miss_min[i];
    }
}
model{
    vector[N] logit_p;
    matrix[2,2] SIGMA;
    vector[N] age_merge;

    // priors
    a ~ normal(0,10);
    b_age ~ normal(0,1);
    sigma_pop ~ exponential(1);
    Rho_pop ~ lkj_corr(4);
    SIGMA = quad_form_diag(Rho_pop,sigma_pop);
    v_pop ~ multi_normal( zeros , SIGMA );

    // age imputation
    for ( i in 1:N ) {
        if ( age[i]<0 ) {
            // age missing
            // uniform bounds should be done from parameter declaration block
            // so figure out which extra information we have to impute with
            // prior

            // insert right parameter
            age_merge[i] = age_impute[ miss_idx[i] ];
        } else {
            // age observed
            // so just insert into age_merge
            age_merge[i] = age[i];
        }
    }

    // likelihood
    for ( i in 1:N ) {
        logit_p[i] = a + v_pop[pop_id[i],1] + 
                    (b_age + v_pop[pop_id[i],2])*age_merge[i];
    }
    y ~ bernoulli_logit( logit_p );
}
