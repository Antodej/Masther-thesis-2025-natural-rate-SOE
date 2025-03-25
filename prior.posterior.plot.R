prior.posterior.plot <- function(fit_posteriors, fit_priors){
  
  # Create a combined long-format data frame
  prior_df <- as_draws_df(fit_priors)
  posterior_df <- as_draws_df(fit_posteriors)
  
  prior_long <- prior_df %>%
    select(all_of(names(fit_posteriors))) %>%
    mutate(sample = "Prior") %>%
    pivot_longer(cols = all_of(names(fit_posteriors)), names_to = "Parameter", values_to = "Value")
  
  posterior_long <- posterior_df %>%
    select(all_of(names(fit_posteriors))) %>%
    mutate(sample = "Posterior") %>%
    pivot_longer(cols = all_of(names(fit_posteriors)), names_to = "Parameter", values_to = "Value")
  
  plot1_df <- bind_rows(prior_long, posterior_long)
  
  plot1 = ggplot(plot1_df, aes(x = Value, color = sample, linetype = sample)) +
    geom_density(linewidth = 5) +
    facet_wrap(~ Parameter, scales = "free", ncol=6) +
    scale_color_manual(values = c("Prior" = "blue", "Posterior" = "red")) +
    scale_linetype_manual(values = c("Prior" = "dotted", "Posterior" = "solid")) +
    theme_minimal(base_size=40) +
    labs(x = "Parameter Value", y = "Density") +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
  return(plot1)
}