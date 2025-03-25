
#============
# This code uses many functions from the code of Holston, Laubach and Williams (2023)
# It has been augmented to estimate r* with Bayesian method through the cmdstanr package
# This new estimation follows the model of Berger and Kempa (2014)
# It works alongside the files sourced here
# Antoine DEJEAN May 2025
#==========

rm(list=ls())

# =================
# DEFINE DIRECTORIES
# =================

# This directory should contain
#   - an 'inputData' folder with the data
#   - an 'output' folder to store estimation results
working.dir <- 'C:/R/master thesis/HLW_edited_code/code_thesis_DEJEAN_1303'

# Location of model code files
code.dir    <- 'C:/R/master thesis/HLW_edited_code/code_thesis_DEJEAN_1303'

# =================
# LOAD R PACKAGES
# =================

library("tis") # Time series package
library("mFilter") # HP filter
library("openxlsx") # Input from and write to Excel
library("cmdstanr") #load and run STAN models
library("posterior")
library("bayesplot") 
library("ggplot2") #graphs
library("dplyr") #data manip
library("tidyr")
library("purrr")
# ==================
# LOAD CODE PACKAGES
# ==================

setwd(code.dir)
source("kalman.log.likelihood.buncic.R") #from Buncic (202X), edited
source("kalman.states.R") #from HLW (2023)
source("kalman.states.wrapper.SOE.R") 
source("log.likelihood.wrapper.SOE.R") 
source("unpack.parameters.SOE.R") 
source("utilities.R") #from HLW (2023)
source("prior.posterior.plot.R") #to create a prior posterior plot
source("kalman.states.posterior.wrapper.SOE.R") #to output possible states
source("percentiles.SOE.R") #to output means and 90% CI interval for posterior parameters

# Set working directory back to output location
setwd(working.dir)


# Set the start and end dates of the estimation sample (format is c(year,quarter))
sample.start <- c(1972,1)
sample.end   <- c(2024,4)

# The estimation process uses data beginning 4 quarters prior to the sample start
data.start    <- shiftQuarter(sample.start,-4)

# Set start index for g.pot series; used in state vector initialization
g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(data.start,'quarterly')

# =================
# COVID-ADJUSTED MODEL SETTINGS
# =================

# Set to TRUE if using a sample covering COVID years
# Must specify kappa.inputs if TRUE
use.kappa <- TRUE #TRUE

kappa.inputs <- data.frame('name'=c('kappa2020Q2-Q4','kappa2021','kappa2022'),
                           'year'=c(2020,2021,2022),
                           'T.start'=c(NA,NA,NA),
                           'T.end'=c(NA,NA,NA),
                           'init'=c(5,1,1),
                           'lower.bound'=c(1,1,1),
                           'upper.bound'=c(Inf,Inf,Inf),
                           'theta.index'=c(NA,NA,NA),
                           't.stat.null'=c(1,1,1))

# NOTE: Sets Q1-Q4 of years provided
if (use.kappa) {
  # Number of kappas introduced
  n.kappa <- dim(kappa.inputs)[1]
  for (k in 1:n.kappa) {
    # Indexing to start of y_t vector
    covid.variance.start.yq <- c(kappa.inputs$year[k],1) - sample.start
    
    kappa.inputs$T.start[k] <- max(covid.variance.start.yq[1]*4 + covid.variance.start.yq[2] +1,0)
    
    covid.variance.end.yq <- c(kappa.inputs$year[k],4) - sample.start
    
    kappa.inputs$T.end[k] <- max(covid.variance.end.yq[1]*4 + covid.variance.end.yq[2] +1,0)
    
    rm(covid.variance.start.yq, covid.variance.end.yq)
    
    # Manual adjustment to start Kappa_2020 in second quarter
    # Comment out under alternative specifications
    if (kappa.inputs$year[k]==2020) {
      kappa.inputs$T.start[k] <- kappa.inputs$T.start[k] + 1
    }
  }
}

# =================
# INPUT DATA
# =================

# Read input data from folder
soe.data <- read.xlsx("inputData/DEJEAN_thesis_data_ER.xlsx", sheet="BE1971Q12024Q4",
                      na.strings = ".", colNames=TRUE, rowNames=FALSE, detectDates = TRUE)

log.output             <- soe.data$gdp.log
inflation              <- soe.data$inflation
inflation.expectations <- soe.data$inflation.expectations
nominal.interest.rate  <- soe.data$interest
real.interest.rate     <- nominal.interest.rate - inflation.expectations
covid.indicator        <- soe.data$covid.ind
exchange.rate         <- soe.data$exchange.rate


#----------------------------------------------------------------------------#
# Obtain initial parameter values
#----------------------------------------------------------------------------#

# Data must start 4 quarters before the estimation period
t.end  <- length(log.output) - 4

# Estimate log-linear trend in GDP via OLS
x.og <- cbind(rep(1,t.end+4), 1:(t.end+4))
y.og <- log.output
output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100

g.pot <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend # data.start : sample.end  
g.pot.diff <- diff(g.pot) # (data.start+1) : sample.end
  
#A: inititialization of the state vectors with estimated output and growth values
xi.00 <- c(100*g.pot[(g.pot.start.index+2):(g.pot.start.index)],100*g.pot.diff[(g.pot.start.index-1+2):((g.pot.start.index-1))],0,0,0,exchange.rate[5],exchange.rate[4],exchange.rate[3],0,0,0)

#A: many new params! When possible, starting value from Pedersen (2015)

# Initial parameter values
#initial.parameters <- c(b.is["a_1"], b.is["a_2"], b.is["a_r"], b.is["a_q"], b.ph[1], b.ph[3],b.ph[4], s.is, s.ph, 0.7, b.is["phi"], 1,
                       # 0.5,0.03,1.5,-0.7,3,3,0.5, 0.5, 0.1)
# Set parameter labels
param.num <- c("a_y1"=1,"a_y2"=2,"a_r"=3, "a_q" = 4, "b_pi"=5,"b_y"=6,
               "b_q"=7, "sigma_ygap_s"=8,"sigma_pi_s"=9,"sigma_ystar_s"=10,"phi"=11,"c"=12,
               "rho_u"=13, "m"=14, "delta_1"=15, "delta_2"=16, "sigma_qgap_s"=17,
               "sigma_qstar_s"=18, "sigma_u_s"=19, "sigma_z_s"=20, "sigma_g_s"=21)
n.params <- length(param.num)


#----------------------------------------------------------------------------#
# Build data matrices
#----------------------------------------------------------------------------#
y.data <- cbind(100 * log.output[5:(t.end+4)],
                inflation[5:(t.end+4)], exchange.rate[5:(t.end+4)],real.interest.rate[5:(t.end+4)])

x.data <- cbind(100 * log.output[4:(t.end+3)],
                100 * log.output[3:(t.end+2)],
                real.interest.rate[4:(t.end+3)],
                real.interest.rate[3:(t.end+2)],
                exchange.rate[4:(t.end+3)],
                exchange.rate[3:(t.end+2)],
                inflation[4:(t.end+3)],
                (inflation[3:(t.end+2)]+inflation[2:(t.end+1)]+inflation[1:t.end])/3,
                covid.indicator[5:(t.end+4)],
                covid.indicator[4:(t.end+3)],
                covid.indicator[3:(t.end+2)])

# Define Kappa
if (use.kappa) {
  
  n.kappa <- dim(kappa.inputs)[1]
  
  for (k in 1:n.kappa) {
    theta.ind <- n.params + k
    kappa.inputs$theta.index[k] <- theta.ind # Store index position within theta in kappa.inputs
    param.num[kappa.inputs$name[k]] <- theta.ind # Store index position within theta in param.num
    # Print statement for kappa value
    if (kappa.inputs$lower.bound[k]==kappa.inputs$upper.bound[k]) {
      print(paste0("Fixing ",kappa.inputs$name[k]," at ",as.character(kappa.inputs$lower.bound[k])))
    } else {
      print(paste0("Initializing ",kappa.inputs$name[k]," at ",as.character(kappa.inputs$init[k])))
    }
    
  }
}  else {
  param.num <- c("a_y1"=1,"a_y2"=2,"a_r"=3, "a_q" = 4, "b_pi"=5,"b_y"=6,
                 "b_q"=7, "sigma_ygap_s"=8,"sigma_pi_s"=9,"sigma_ystar_s"=10,"phi"=11,"c"=12,
                 "rho_u"=13, "m"=14, "delta_1"=15, "delta_2"=16, "sigma_qgap_s"=17,
                 "sigma_qstar_s"=18, "sigma_u_s"=19, "sigma_z_s"=20, "sigma_g_s"=21, "NOT USED"=22, "NOT USED"=23, "NOT USED"=24)#, "NOT USED"=25)
}

#initial covariance matrix
P.00=diag(0.2,length(xi.00))

kappa_fixed <- rep(1, length(log.output)-4)  

if (!is.na(kappa.inputs$T.start[1])) {
  varying_kappa_indices_2020 <- kappa.inputs$T.start[1]:kappa.inputs$T.end[1]   # For Q2, Q3, Q4 in 2020
  varying_kappa_indices_2021 <- kappa.inputs$T.start[2]:kappa.inputs$T.end[2]     # For all of 2021
  varying_kappa_indices_2022 <- kappa.inputs$T.start[3]:kappa.inputs$T.end[3]   # For all of 2022
 # varying_kappa_indices_2023 <- kappa.inputs$T.start[4]:kappa.inputs$T.end[4]  #look at Kappa inputs
  
  for (i in c(varying_kappa_indices_2020, varying_kappa_indices_2021, varying_kappa_indices_2022)) {
    kappa_fixed[i] <- NA  
  }
  
} else {
  varying_kappa_indices_2020 <- c(0,0,0)  
  varying_kappa_indices_2021 <- c(0,0,0,0)
  varying_kappa_indices_2022 <- c(0,0,0,0)
}

if (use.kappa==TRUE) { #communicate to STAN whether we estimate covid parameters
  use_kappa_int=1
} else {use_kappa_int=0
}

#STAN data prep
# Data list for Stan
data_list <- list(
  T = length(log.output)-4,
  n = length(xi.00),
  y = y.data,
  x = x.data,
  kappa = rep(1, length(log.output)-4),
  xi0 = xi.00,
  P0 = P.00,
  varying_kappa_indices_2020 = varying_kappa_indices_2020,  
  varying_kappa_indices_2021 = varying_kappa_indices_2021,  
  varying_kappa_indices_2022 = varying_kappa_indices_2022,  
  use_kappa_int=use_kappa_int
)

#starting values: mean of prior distribution, 1 for kappas
init_fun <- function() {
  list(
    a_y1=1.5,
    a_y2=-0.7,
    a_r=-0.1,
    a_q=-0.001,
    b_pi=0.5,
    b_y=0.5,
    b_q=-0.25,
    m=0.001,
    phi=-0.001,
    c=1,
    delta_1=1.5,
    delta_2=-0.7,
    rho_u=0.5,
    kappa2020=1.00001,
    kappa2021=1.00001,
    kappa2022=1.00001,
    sigma_ygap_s = 0.5,
    sigma_ystar_s = 0.25,
    sigma_qgap_s = 3,
    sigma_qstar_s = 1,
    sigma_pi_s = 3,
    sigma_u_s = 0.5,
    sigma_z_s = 0.5,
    sigma_g_s = 0.25 )
}

#STAN execution
# Compile the model
compiled_model = cmdstan_model("STAN_thesis_model.stan")
fit <- compiled_model$sample(
  data = data_list, 
  init = init_fun, 
  iter_warmup = 10000, 
  iter_sampling = 10000,  
  chains = 2,  # 
  parallel_chains = 3,
  refresh = 10,
  max_treedepth=12
)

#get states with the posterior draws
fit_results <- as_draws_df(fit$draws())
fit_results= fit_results[,2:25]
state.names = c("ystar_t","ystar_t1","ystar_t2","g_t","g_t1","g_t2","z_t","z_t1","z_t2","qstar_t","qstar_t1","qstar_t2","u_t","u_t1","u_t2")

states.bayesian.estimation <- kalman.states.posterior.wrapper.SOE(posterior.draws=fit_results, y.data=y.data, x.data=x.data,
                                  xi.00=xi.00, P.00=P.00, use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num,
                                  state.names=state.names, start=sample.start, end=sample.end, country.code="BE")
summary.smoothed=states.bayesian.estimation$tidy.smoothed.states
summary.filtered=states.bayesian.estimation$tidy.filtered.states

#90% interval (5th-95th percentile)
table_CI_mean <- percentiles.SOE(fit_results)

#PRIOR VS POSTERIOR
priors_model = cmdstan_model("STAN_priors_model.stan")
fit_priors <- priors_model$sample(seed=121, data = list(use_kappa_int=use_kappa_int), 
                                  iter_warmup = 1000, iter_sampling = 10000,  chains = 4, parallel_chains = 4, refresh = 10000)

PLOT1 = prior.posterior.plot(fit_results, fit_priors)

#facet plot of state estimates
PLOT2 <- ggplot(summary.smoothed, aes(x = time)) +
  geom_line(aes(y = mean), linewidth = 1.2) +  # Bold mean line
  geom_line(aes(y = lower), linetype = "dotted") +  # Lower CI as dotted line
  geom_line(aes(y = upper), linetype = "dotted") +  # Upper CI as dotted line
  facet_wrap(~ state, scales = "free_y") +
  labs(title = "Smoothed State Estimates with Credible Bands",
       x = "Time", y = "State Value")
PLOT3 <- ggplot(summary.filtered, aes(x = time)) +
  geom_line(aes(y = mean), linewidth = 1.2) +  # Bold mean line
  geom_line(aes(y = lower), linetype = "dotted") +  # Lower CI as dotted line
  geom_line(aes(y = upper), linetype = "dotted") +  # Upper CI as dotted line
  facet_wrap(~ state, scales = "free_y") +
  labs(title = "Filtered State Estimates with Credible Bands",
       x = "Time", y = "State Value")

#outpout plot
output_plot_df = summary.smoothed %>% filter(state == "ystar_t")
output_plot_df2 = data.frame(time = output_plot_df$time, y_t = y.data[,1])

PLOT4 <- ggplot(output_plot_df, aes(x = time)) +
  geom_line(aes(y = mean), linewidth = 1.2) +  # Bold mean line
  geom_line(aes(y = lower), linetype = "dotted") +  # Lower CI as dotted line
  geom_line(aes(y = upper), linetype = "dotted") +
  geom_line(data = output_plot_df2, aes(y = y_t), color = "red") +  # Add y.data values
  labs(title = "Time Series Plot", x = "Time", y = "Values") +
  theme_minimal()



#plot for thesis
