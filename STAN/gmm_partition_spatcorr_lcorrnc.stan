/*
 * Model to partition total residuals
 * estimate a coreg model at the same time
 */

#include functions.stan

data {
  int N;
  int NEQ;
  int NSTAT;
  int M;
  
  vector[N] Y; // total residuals

  array[N] int<lower=1,upper=NEQ> eq;
  array[N] int<lower=1,upper=NSTAT> stat;

  array[N] vector[2] X;
  array[NEQ] vector[2] X_e;

  array[NEQ, M] int ind_eq;
  array[NEQ] int len_eq;
}

parameters {
  real<lower=0> sigma_rec;
  real<lower=0> sigma_eq;
  real<lower=0> sigma_stat;

  vector[NEQ] eqterm;
  vector[NSTAT] statterm;

  real ic; // intrcept for each period
  real mu_log_ell;
  real<lower=0> sigma_log_ell;
  real<lower=0> range_ell; // length scale for ell

  vector[NEQ] z_log_ell;

  real<lower=0,upper=1> omega; // variance weights for ph_SS
}

transformed parameters {
  real theta = sqrt(omega .* square(sigma_rec));
  real theta2 = sqrt((1 - omega) .* square(sigma_rec));

  vector[NEQ] ell = exp(mu_log_ell + calc_f_M1(X_e, NEQ, sigma_log_ell, range_ell, z_log_ell));
}

model {
  // priors
  ic ~ normal(0,0.2);
  sigma_rec ~ normal(0,0.5);
  sigma_eq ~ normal(0,0.5);
  sigma_stat ~ normal(0,0.5);

  omega ~ beta(2,2);

  mu_log_ell ~ normal(0, 10);
  sigma_log_ell ~ normal(0, 10);
  range_ell ~ inv_gamma(2.1, 34);
  z_log_ell ~ std_normal();

  eqterm ~ normal(0, sigma_eq);
  statterm ~ normal(0, sigma_stat);

  vector[N] mu = ic + eqterm[eq] + statterm[stat];

  for(i in 1:NEQ) {
    matrix[len_eq[i], len_eq[i]] L;
    L = calc_L2_M1(X[ind_eq[i, 1:len_eq[i]]], len_eq[i], theta, ell[i], theta2);

    Y[ind_eq[i, 1:len_eq[i]]] ~ multi_normal_cholesky(mu[ind_eq[i, 1:len_eq[i]]], L);
  }
}
