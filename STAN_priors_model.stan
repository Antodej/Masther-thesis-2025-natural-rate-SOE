data {
  int use_kappa_int; //use covid vars or no
}

parameters {
  real a_y1;
  real a_y2;
  real a_r;
  //real<upper=0> a_q;
 real a_q;
  real b_pi;
  real b_y;
  real b_q;
  real<lower=0.0001> sigma_ygap_s;
  real<lower=0.0001> sigma_pi_s;
  real<lower=0.0001> sigma_ystar_s;
  real phi;
  real c;
  real<lower=0, upper=1> rho_u;
 // real<lower=0> m;
  real m;
  real delta_1;
  real delta_2;
  
  real<lower=0.0001> sigma_qgap_s;
  real<lower=0.0001> sigma_qstar_s;
  real<lower=0.0001> sigma_u_s;
  real<lower=0.0001> sigma_z_s;
  real<lower=0.0001> sigma_g_s;
  
  real<lower=1, upper=50> kappa2020;
  real<lower=1, upper=50> kappa2021;
  real<lower=1, upper=50> kappa2022;
//  real<lower=1, upper=50> kappa2023;
  
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
  
 sigma_g_s  ~ gamma(1, 4);
 sigma_z_s ~ gamma(3, 20);
   sigma_ygap_s ~ gamma(4, 8);
   
    sigma_ystar_s ~ gamma(3, 20);
sigma_qstar_s ~ gamma(9, 12);
  sigma_qgap_s ~ gamma(9, 3);
  sigma_pi_s ~ gamma(4.2, 1.4);
sigma_u_s ~ gamma(4, 8);
}

