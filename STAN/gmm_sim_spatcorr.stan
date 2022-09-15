/*
 * Model to simulate total residuals
 */

#include functions.stan

data {
  int N;
  int NEQ;
  int NSTAT;
  int M;
  
  array[N] int<lower=1,upper=NEQ> eq;
  array[N] int<lower=1,upper=NSTAT> stat;

  array[N] vector[2] X;

  array[NEQ, M] int ind_eq;
  array[NEQ] int len_eq;

  real<lower=0> sigma_rec;
  real<lower=0> sigma_eq;
  real<lower=0> sigma_stat;
  real<lower=0> ell; // length scales
  real<lower=0,upper=1> omega; // variance weights for ph_SS
}

transformed data {
  real theta = sqrt(omega .* square(sigma_rec));
  real theta2 = sqrt((1 - omega) .* square(sigma_rec));
}

parameters {

}

model {
}

generated quantities {
  array[NEQ] real eqterm = normal_rng(rep_vector(0, NEQ), sigma_eq);
  array[NSTAT] real statterm = normal_rng(rep_vector(0, NSTAT), sigma_stat);;
  vector[N] Y_sim;

  {
    vector[N] mu = to_vector(eqterm)[eq] + to_vector(statterm)[stat];

    for(i in 1:NEQ) {
      matrix[len_eq[i], len_eq[i]] L;
      L = calc_L2_M1(X[ind_eq[i, 1:len_eq[i]]], len_eq[i], theta, ell, theta2);

      Y_sim[ind_eq[i, 1:len_eq[i]]] = multi_normal_cholesky_rng(mu[ind_eq[i, 1:len_eq[i]]], L);
  }

  }
}
