functions {
  // Compute the cumulative log likelihood via the Kalman filter
  real kalman_log_likelihood_SOE(
    matrix F,         // State transition matrix
    matrix Q,         // Process noise covariance
    matrix A_t,  
    matrix A, // Regression coefficients (if applicable)
    matrix H_t, 
    matrix H, // Observation matrix
    matrix R,         // Observation noise covariance
    vector kappa,     // Time-varying scaling (length T)
    matrix y,         // Observations (T x 4)
    matrix x,         // Covariates (T x 11)
    vector xi0,       // Initial state vector (n x 1)
    matrix P0,         // Initial state covariance (nxn)
    int T,
    int n
  ) {
    real ll_cum = 0;
    vector[n] xi_tt;
    xi_tt = xi0;

    matrix[n, n] P_tt;
    P_tt = P0;
    matrix[n, n] t_F;
    t_F = F';    // Transpose of F
    // For speed, we precompute a constant:
    real log_2_pi;
    log_2_pi = 1.8378770664093454835606594728112353;
    
   for (t in 1:T) {
      // Predict state
      vector[n] xi_ttm1 = F * xi_tt;
      matrix[n, n] P_ttm1 = F * P_tt * t_F + Q;
      
      // Compute prediction error:
      // Note: In Stan, a row of y is a row vector; we convert it to a column vector via transpose.
    
      vector[11] x_t = x[t,]';
      vector[4] y_t = y[t,]';
      vector[4] pred_err = y_t - A_t * x_t - H_t * xi_ttm1;
      
      // Compute the innovation covariance, kappa does not affect qgap:
      real kappa_t_s = square(kappa[t]);
      matrix[4,4] kappa_mat = diag_matrix(to_vector([kappa_t_s, kappa_t_s, 1,1]));
      matrix[4, 4] HPHR = H_t * P_ttm1 * H + kappa_mat * R;
     // matrix[4, 4] HPHR = H_t * P_ttm1 * H + square(kappa[t]) * R;
      
      // Solve for the weighted error: HPHR \ pred_err
      vector[4] inv_HPHR_hat = HPHR \ pred_err;
      
      // Update cumulative log likelihood:
      ll_cum += - (4 / 2.0) * log_2_pi - 0.5 * log_determinant(HPHR)
                - 0.5 * dot_product(pred_err, inv_HPHR_hat);
      
      // Update the state estimate:
      xi_tt = xi_ttm1 + P_ttm1 * H * inv_HPHR_hat;
      
      // Update the state covariance:
      P_tt = P_ttm1 - P_ttm1 * H * (HPHR \ (H_t * P_ttm1));
      
    }
    return ll_cum;
  }
}


data {
  int<lower=1> T;         // Number of time points
  int<lower=1> n;         // Dimension of observation vector y_t
  int use_kappa_int; //use covid vars or no
  matrix[T, 4] y;         // Observations y, pi, q, r
  matrix[T, 11] x;         // lagged variables + CSI (exogeneous)
  array[3] int varying_kappa_indices_2020;
  array[4] int varying_kappa_indices_2021;
  array[4] int varying_kappa_indices_2022;
  // Initial state and covariance:
  vector[n] xi0;
  matrix[n, n] P0;
}

parameters {
  real a_y1;
  real a_y2;
 // real <upper=0> a_r;
real<upper=0> a_r;
 real a_q;
  real b_pi;
  real b_y;
  real b_q;
  real<lower=0.001, upper=50> sigma_ygap_s;
  real<lower=0.001, upper=50> sigma_pi_s;
  real<lower=0.001, upper=50> sigma_ystar_s;
  real phi;
  real c;
  real<lower=0.0001, upper=1> rho_u;
 // real<lower=0> m;
 real m;
  real delta_1;
  real delta_2;
  
  real<lower=0.001, upper=50> sigma_qgap_s;
  real<lower=0.001, upper=50> sigma_qstar_s;
  real<lower=0.001, upper=50> sigma_u_s;
  real<lower=0.001, upper=50> sigma_z_s;
  real<lower=0.0001, upper=50> sigma_g_s;
  
  real<lower=1, upper=50> kappa2020;
  real<lower=1, upper=50> kappa2021;
  real<lower=1, upper=50> kappa2022;
//  real<lower=1> kappa2023;
  
}

transformed parameters {
  matrix[11,4] A;  // Create a 4x11 matrix
  matrix[4, 11] A_t;
  // Initialize A to zeros:
  A_t = rep_matrix(0.0, 4, 11);
  
  // Row 1
  A_t[1, 1] = a_y1;
  A_t[1, 2] = a_y2;
  A_t[1, 3] = a_r / 2;
  A_t[1, 4] = a_r / 2;
  A_t[1, 5] = a_q/2;
  A_t[1, 6] = a_q/2;
  A_t[1, 9] = phi;
  A_t[1, 10] = -a_y1 * phi;
  A_t[1, 11] = -a_y2 * phi;
  
  // Row 2
  A_t[2, 1] = b_y;
  A_t[2, 5] = b_q;
  A_t[2, 6] = -b_q;
  A_t[2, 7] = b_pi;
  A_t[2, 8] = 1 - b_pi;
  A_t[2, 10] = -b_y * phi;
  
  // Row 3
  A_t[3, 5] = delta_1;
  A_t[3, 6] = delta_2;
  
  // Row 4
  A_t[4, 5] = m;
  
  A= transpose(A_t);

  matrix[15, 4] H;
  matrix[4, 15] H_t; 
  H_t = rep_matrix(0.0, 4, 15);

  // Row 1
  H_t[1, 1] = 1;
  H_t[1, 2] = -a_y1;
  H_t[1, 3] = -a_y2;
  H_t[1, 5] = -c * a_r * 2;
  H_t[1, 6] = -c * a_r * 2;
  H_t[1, 8] = -a_r / 2;
  H_t[1, 9] = -a_r / 2;
  H_t[1, 11] = -a_q/2; 
  H_t[1, 12] = -a_q/2; 

  // Row 2
  H_t[2, 2] = -b_y;

  // Row 3
  H_t[3, 10] = 1;
  H_t[3, 11] = -delta_1;
  H_t[3, 12] = -delta_2;

  // Row 4
  H_t[4, 5] = c * 4;
  H_t[4, 8] = 1;
  H_t[4, 11] = -m;
  H_t[4, 14] = 1;

  // Compute transposed H
H = transpose(H_t);
  
   // Define and initialize R (4x4 diagonal matrix)
  matrix[4, 4] R;
R = diag_matrix(to_vector([
  sigma_ygap_s,
  sigma_pi_s,
  sigma_qgap_s,0]));
  // Define and initialize Q (15x15)
  matrix[15, 15] Q;
  Q = rep_matrix(0.0, 15, 15);

  Q[1, 1] = sigma_ystar_s;
  Q[4, 4] = sigma_g_s;
  Q[7, 7] = sigma_z_s;
  Q[10, 10] = sigma_qstar_s;
  Q[13, 13] = sigma_u_s;

  // Define and initialize F (15x15)
  matrix[15, 15] F;
  F = rep_matrix(0.0, 15, 15);

  // Assign individual elements for F
  F[1, 1] = 1;
  F[1, 4] = 1;
  F[2, 1] = 1;
  F[3, 2] = 1;
  F[4, 4] = 1;
  F[5, 4] = 1;
  F[6, 5] = 1;
  F[7, 7] = 1;
  F[8, 7] = 1;
  F[9, 8] = 1;
  F[10, 10] = 1;
  F[11, 10] = 1;
  F[12, 11] = 1;
  F[13, 13] = rho_u;
  F[14, 13] = 1;
  F[15, 14] = 1;
  
    vector[T] kappa;  
  
  // Initialize kappa with 1 (for all time points where kappa is fixed)
  kappa = rep_vector(1.0, T);  
  
if (use_kappa_int == 1 )  {
  for (i in 1:size(varying_kappa_indices_2020))
  kappa[varying_kappa_indices_2020[i]] = kappa2020;
  for (i in 1:size(varying_kappa_indices_2021))
  kappa[varying_kappa_indices_2021[i]] = kappa2021;
  for (i in 1:size(varying_kappa_indices_2022))
  kappa[varying_kappa_indices_2022[i]] = kappa2022;
 // for (i in 1:size(varying_kappa_indices_2023))
 // kappa[varying_kappa_indices_2023[i]] = kappa2023;
}
  
}


model {

  a_y1 ~ normal(1.5, 0.25);
  a_y2 ~ normal(-0.7, 0.25);
  a_r  ~ normal(-0.1, 0.1);
  a_q  ~ normal(0, 0.15);
  b_pi ~ normal(0.5, 0.25);
  
  if (use_kappa_int == 0 )  {
     phi ~ normal(0, 0.001);
    kappa2020 ~ normal(1, 0.001);
    kappa2021 ~ normal(1, 0.001);
    kappa2022 ~ normal(1, 0.001);
  } else { 
     phi ~ normal(-0.1, 0.1);}
    
  delta_1 ~ normal(1.5, 0.25);
  delta_2 ~ normal(-0.7, 0.25);
  b_y  ~ normal(0.5, 0.25);
m    ~ normal(0, 0.5);
  c    ~ normal(1, 0.1);
  b_q  ~ normal(-0.25, 0.1);
  rho_u ~ normal(0.5, 0.25);
  
  //kappas priors implicit (uniform)
  
  //sigma_ystar_s ~ gamma(0.0125, 0.05);
  // sigma_qgap_s ~ gamma(1.8, 0.6);
  // sigma_qstar_s ~ gamma(1.8, 0.6);
  // sigma_pi_s ~ gamma(1.8, 0.6);
//  sigma_g_s  ~ gamma(0.0125, 0.05);
 // sigma_ygap_s ~ gamma(4, 4);
 sigma_ystar_s ~ gamma(3, 20);
//   //sigma_ystar_s ~ gamma(1, 4);
sigma_qgap_s ~ gamma(9, 3);
sigma_qstar_s ~ gamma(9,12);
//   sigma_pi_s ~ gamma(3, 1);
// sigma_u_s ~ gamma(1.5, 0.5);
// //sigma_u_s ~ gamma(9, 3);
// sigma_z_s ~ normal(0, 0.1);
//   sigma_g_s  ~ normal(0, 0.001);
sigma_g_s  ~ gamma(1, 4);
 sigma_z_s ~ gamma(3, 20);
  
  // orignal TB values
   sigma_ygap_s ~ gamma(4, 8);
  sigma_ystar_s ~ gamma(4, 16);
 sigma_qgap_s ~ gamma(4.2, 1.4);
 //sigma_qstar_s ~ gamma(4.2, 1.4);
  sigma_pi_s ~ gamma(4.2, 1.4);
sigma_u_s ~ gamma(4, 8);
 //sigma_z_s ~ gamma(4, 8);
 // sigma_g_s  ~ gamma(4, 16);
  
  
  // Add the Kalman filter log likelihood to the target
  target += kalman_log_likelihood_SOE(F, Q, A_t, A, H_t, H, R, kappa, y, x, xi0, P0, T, n);
}

