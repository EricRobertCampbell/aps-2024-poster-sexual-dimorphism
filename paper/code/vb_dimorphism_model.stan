functions {
    vector vonBertalanffy(vector age, real L, real A, real K) {
        int N = size(age);
        vector[N] means = L * (1 - A * exp(-K * age));
        return means;
    }
    vector orZero(vector v) {
        int N = size(v);
        vector[N] newV;
        for (i in 1:N) {
            if (v[i] <= 0) {
                newV[i] = 0;
            } else {
                newV[i] = v[i];
            }
        }
        return newV;
    }
}
data {
    // male data
    int<lower = 0> N_male; // number of samples
    vector<lower=0>[N_male] age_male; // vector of the ages; index correlates with length (so length[i] happened at age[i])
    vector<lower=0>[N_male] length_male; // a vector of the lengths

    // female data
    int<lower = 0> N_female;
    vector<lower=0>[N_female] age_female;
    vector<lower=0>[N_female] length_female;

    // population data
    // prior data (mean and standard deviation)
    vector[2] prior_L;
    vector[2] prior_A;
    vector[2] prior_K;
    vector[2] prior_sigma;
}

parameters {
    real<lower=0> L_male;
    real<lower=0> A_male;
    real<lower=0> K_male;
    real<lower=0> sigma_male;

    real<lower=0> L_female;
    real<lower=0> A_female;
    real<lower=0> K_female;
    real<lower=0> sigma_female;
}

model {
    vector[N_male] mean_length_male;
    vector[N_male] sd_length_male;
    vector[N_female] mean_length_female;
    vector[N_female] sd_length_female;

    // male part of the model
    L_male ~ normal(prior_L[1], prior_L[2]);
    K_male ~ normal(prior_K[1], prior_K[2]);
    A_male ~ normal(prior_A[1], prior_A[2]);
    sigma_male ~ normal(prior_sigma[1], prior_sigma[2]);
    mean_length_male = vonBertalanffy(age_male, L_male, A_male, K_male);
    length_male ~ normal(mean_length_male, sigma_male);

    // female part of the model
    L_female ~ normal(prior_L[1], prior_L[2]);
    K_female ~ normal(prior_K[1], prior_K[2]);
    A_female ~ normal(prior_A[1], prior_A[2]);
    sigma_female ~ normal(prior_sigma[1], prior_sigma[2]);
    mean_length_female = vonBertalanffy(age_female, L_female, A_female, K_female);
    length_female ~ normal(mean_length_female, sigma_female);
}

generated quantities {
    real diff;
    diff = L_male - L_female;
}
