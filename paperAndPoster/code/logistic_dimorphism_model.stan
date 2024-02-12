data {
    // large sex data
    int<lower = 0> N_large; // number of males
    vector[N_large] age_large;
    vector[N_large] mass_large;

    // small sex data
    int<lower = 0> N_small;
    vector[N_small] age_small;
    vector[N_small] mass_small;

    // priors
    real prior_L_mean;
    real prior_L_sd;
    real prior_k_mean;
    real prior_k_sd;
    real prior_q_mean;
    real prior_q_sd;
    real prior_sigma_mean;
    real prior_sigma_sd;
}

parameters {
    real L_large;
    real k_large;
    real q_large;

    real L_small;
    real k_small;
    real q_small;

    real mass_sd;
}

model {
    mass_sd ~ normal(prior_sigma_mean, prior_sigma_sd);

    // large sex model
    L_large ~ normal(prior_L_mean, prior_L_sd);
    k_large ~ normal(prior_k_mean, prior_k_sd);
    q_large ~ normal(prior_q_mean, prior_q_sd);
    mass_large ~ normal(L_large * inv((1 + exp(q_large + k_large * age_large))), mass_sd);
    
    // small sex model
    L_small ~ normal(prior_L_mean, prior_L_sd);
    k_small ~ normal(prior_k_mean, prior_k_sd);
    q_small ~ normal(prior_q_mean, prior_q_sd);
    mass_small ~ normal(L_small * inv((1 + exp(q_small + k_small * age_small))), mass_sd);
}

generated quantities {
    real difference;
    difference = L_large - L_small;

    real relative_difference;
    relative_difference = L_large / L_small;
}
