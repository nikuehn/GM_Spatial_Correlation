functions {

/*
 *
 * */
  matrix calc_L_M1(array[] vector X, int N, real theta, real rho) {
    matrix[N, N] K;
    matrix[N, N] L;
    real delta = 1e-9;

    for(i in 1:N) {
      for(j in (i + 1):N) {
        K[i,j] = square(theta) * exp(- distance(X[i], X[j]) / rho);
        K[j,i] = K[i,j];
      }
      K[i,i] = square(theta) + delta;
    }
    L = cholesky_decompose(K);
    return L;
  }

/*
 *
 * */
  matrix calc_L2_M1(array[] vector X, int N, real theta, real rho, real theta2) {
    matrix[N, N] K;
    matrix[N, N] L;
    real delta = 1e-9;

    for(i in 1:N) {
      for(j in (i + 1):N) {
        K[i,j] = square(theta) * exp(- distance(X[i], X[j]) / rho);
        K[j,i] = K[i,j];
      }
      K[i,i] = square(theta) + square(theta2) + delta;
    }
    L = cholesky_decompose(K);
    return L;
  }

/*
 *
 * */
  matrix calc_Sigma2_M1(array[] vector X, int N, real theta, real rho, real theta2) {
    matrix[N, N] K;
    real delta = 1e-9;

    for(i in 1:N) {
      for(j in (i + 1):N) {
        K[i,j] = square(theta) * exp(- distance(X[i], X[j]) / rho);
        K[j,i] = K[i,j];
      }
      K[i,i] = square(theta) + square(theta2) + delta;
    }
    return K;
  }

/*
 *
 * */
  matrix calc_L2_M1_tau(array[] vector X, int N, real theta, real rho, real theta2, real tau) {
    matrix[N, N] K;
    matrix[N, N] L;
    real delta = 1e-9;

    for(i in 1:N) {
      for(j in (i + 1):N) {
        K[i,j] = square(theta) * exp(- distance(X[i], X[j]) / rho) + square(tau);
        K[j,i] = K[i,j];
      }
      K[i,i] = square(theta) + square(theta2) + square(tau) + delta;
    }
    L = cholesky_decompose(K);
    return L;
  }


/*
 *
 * */
  vector calc_f_M1(array[] vector X, int N, real theta, real rho, vector z) {
    matrix[N, N] K;
    matrix[N, N] L;
    vector[N] f;
    real delta = 1e-9;

    for(i in 1:N) {
      for(j in (i + 1):N) {
        K[i,j] = square(theta) * exp(- distance(X[i], X[j]) / rho);
        K[j,i] = K[i,j];
      }
      K[i,i] = square(theta) + delta;
    }
    L = cholesky_decompose(K);
    f = L * z;
    return f;
  }


}
