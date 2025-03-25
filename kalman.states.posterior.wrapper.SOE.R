#------------------------------------------------------------------------------#
# File:        kalman.states.posterior.wrapper.SOE.R
#
# Description: This is a function that 
#               loops over selected posteriors to obtain state estimates 
# Antoine DEJEAN (2025)
#------------------------------------------------------------------------------#
kalman.states.posterior.wrapper.SOE <- function(posterior.draws, y.data, x.data,
                                                 xi.00, P.00, use.kappa, kappa.inputs, param.num, state.names, start, end, country.code){
  n.draws=nrow(fit_results)
  state.smoothed.estimates <- vector("list", n.draws)
  state.filtered.estimates <- vector("list", n.draws)
 #create a timed y vector for later addition to state variables
   dates = seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
  y.data.time = as.data.frame(y.data)
  colnames(y.data.time) <- c("y_t", "pi_t", "q_t", "r_t")
  y.data.time <- y.data.time %>%
    mutate(time = dates) %>%
    select(time, everything())
  
  for (i in 1:n.draws) {
    draw <- unlist(posterior.draws[i, ], use.names = TRUE)
    out <- unpack.parameters.SOE(parameters=draw, y.data=y.data, x.data=x.data,
                                        xi.00=xi.00, P.00=P.00,
                                        use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)
      for (n in names(out)) {
        eval(parse(text=paste0(n, "<-out$", n)))
      }
    t.end <- dim(y.data)[1]
    # Run Kalman filter and smoother 
    states <- kalman.states(xi.tm1tm1=xi.00, P.tm1tm1=P.00, F=F, Q=Q, A=A, H=H, R=R, kappa=kappa.vec, y=y.data, x=x.data)
    state.smoothed.estimates[[i]] <- states$smoothed$xi.tT
    state.filtered.estimates[[i]] <- states$filtered$xi.tt
  }
  
  # Convert each draw's matrix into long format and add a draw identifier.
  tidy.smoothed.list <- lapply(seq_along(state.smoothed.estimates), function(i) {
    df <- as.data.frame(state.smoothed.estimates[[i]])
    colnames(df) <- state.names
    df$time <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
    df_long <- pivot_longer(df, cols = all_of(state.names),
                            names_to = "state", values_to = "value")
    df_long$draw <- i
    df_long
  })

  tidy.smoothed.states <- bind_rows(tidy.smoothed.list) %>% select(time, state, draw, value)
  
  summary.smoothed.states <- tidy.smoothed.states %>%
    group_by(time, state) %>%
    summarize(mean = mean(value),
              lower = quantile(value, 0.05),
              upper = quantile(value, 0.95),
              .groups = "drop")  
  
  summary.smoothed.wide <- summary.smoothed.states %>% #add y.data for computation
    pivot_wider(
      names_from = state,
      values_from = c(mean, lower, upper),
      names_sep = "_"
    ) %>%
    left_join(y.data.time, by = "time")
  
  summary.smoothed.wide <- summary.smoothed.wide %>% #compute composite states (r*) and indicators (output gap...)
    mutate(
      mean_rstar_t = 4 * mean_g_t1 + mean_z_t1,
      lower_rstar_t  = 4 * lower_g_t1  + lower_z_t1,
      upper_rstar_t  = 4 * upper_g_t1  + upper_z_t1,
      mean_gann_t = 4 * mean_g_t,
      lower_gann_t  = 4 * lower_g_t,
      upper_gann_t  = 4 * upper_g_t,
      mean_rgap_t = r_t - mean_rstar_t,
      lower_rgap_t  = r_t - lower_rstar_t,
      upper_rgap_t  = r_t - upper_rstar_t,
      mean_qgap_t = q_t - mean_qstar_t,
      lower_qgap_t  = q_t - lower_qstar_t,
      upper_qgap_t  = q_t - upper_qstar_t,
      mean_ygap_t = ((y_t - mean_ystar_t)),
      lower_ygap_t  = ((y_t - lower_ystar_t)),
      upper_ygap_t  = ((y_t - upper_ystar_t))
    )  %>%
    select(time, matches("_t$"))  %>%
    select(-r_t, -q_t, -y_t, -pi_t, -lower_g_t, -mean_g_t, - upper_g_t)
  
 summary.smoothed.states <- summary.smoothed.wide %>%
    pivot_longer(
      cols = -time,  
      names_to = "combined",  
      values_to = "value"
    ) %>%
   separate(
     col = combined, 
     into = c("statistic", "state"), 
     sep = "_", 
     extra = "merge"  # Keeps everything after the first underscore in "statistic"
   ) %>%
   pivot_wider(names_from = statistic, values_from = value) %>%
   mutate(Country = country.code)
  
  #filtered
  
  tidy.filtered.list <- lapply(seq_along(state.filtered.estimates), function(i) {
    df <- as.data.frame(state.filtered.estimates[[i]])
    colnames(df) <- state.names
    df$time <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
    df_long <- pivot_longer(df, cols = all_of(state.names),
                            names_to = "state", values_to = "value")
    df_long$draw <- i
    df_long
  })
  
  tidy.filtered.states <- bind_rows(tidy.filtered.list) %>% select(time, state, draw, value)
  
  summary.filtered.states <- tidy.filtered.states %>%
    group_by(time, state) %>%
    summarize(mean = mean(value),
              lower = quantile(value, 0.05),
              upper = quantile(value, 0.95),
              .groups = "drop") 
  summary.filtered.wide <- summary.filtered.states %>% #add y.data for computation
    pivot_wider(
      names_from = state,
      values_from = c(mean, lower, upper),
      names_sep = "_"
    ) %>%
    left_join(y.data.time, by = "time")
  summary.filtered.wide <- summary.filtered.wide %>% #compute composite states (r*) and indicators (output gap...)
    mutate(
      mean_rstar_t = 4 * mean_g_t1 + mean_z_t1,
      lower_rstar_t  = 4 * lower_g_t1  + lower_z_t1,
      upper_rstar_t  = 4 * upper_g_t1  + upper_z_t1,
      mean_gann_t = 4 * mean_g_t,
      lower_gann_t  = 4 * lower_g_t,
      upper_gann_t  = 4 * upper_g_t,
      mean_rgap_t = r_t - mean_rstar_t,
      lower_rgap_t  = r_t - lower_rstar_t,
      upper_rgap_t  = r_t - upper_rstar_t,
      mean_qgap_t = q_t - mean_qstar_t,
      lower_qgap_t  = q_t - lower_qstar_t,
      upper_qgap_t  = q_t - upper_qstar_t,
      mean_ygap_t = ((y_t - mean_ystar_t)),
      lower_ygap_t  = ((y_t - lower_ystar_t)),
      upper_ygap_t  = ((y_t - upper_ystar_t))
    )  %>%
    select(time, matches("_t$")) %>%
    select(-r_t, -q_t, -y_t, -pi_t,  -lower_g_t, -mean_g_t, - upper_g_t)
  
  summary.filtered.states <- summary.filtered.wide %>%
    pivot_longer(
      cols = -time,  
      names_to = "combined",  
      values_to = "value"
    ) %>%
    separate(
      col = combined, 
      into = c("statistic", "state"), 
      sep = "_", 
      extra = "merge"  # Keeps everything after the first underscore in "statistic"
    ) %>%
    pivot_wider(names_from = statistic, values_from = value) %>%
    mutate(Country = country.code)
  
  
  
  list.states = list(tidy.filtered.states=  summary.filtered.states,   tidy.smoothed.states=  summary.smoothed.states)
  
  return(list.states)
}