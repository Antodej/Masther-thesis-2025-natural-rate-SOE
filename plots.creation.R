rm(list=ls())
library("tis") # Time series package
library("mFilter") # HP filter
library("nloptr") # Optimization
library("openxlsx") # Input from and write to Excel
library("cmdstanr") #load and run STAN models
library("posterior")
library("bayesplot") 
library("ggplot2") #graphs
library("dplyr") #data manip
library("tidyr")
library("purrr")
library("ggh4x")

#import summary.smoothed for all countries
setwd('C:/R/master thesis/HLW_edited_code/code_thesis_DEJEAN_1303/output')

load("summary.smoothed.BE.RData")
summary.smoothed.BE <- summary.smoothed
summary.smoothed <- summary.smoothed %>%
  filter(state == "rstar_t")
time.vector=summary.smoothed[,1]
rm(summary.smoothed)
load("summary.smoothed.NL.RData")
summary.smoothed.NL <- summary.smoothed
rm(summary.smoothed)
load("summary.smoothed.DK.RData")
summary.smoothed.DK <- summary.smoothed
rm(summary.smoothed)
summary.smoothed.all <- bind_rows(summary.smoothed.BE, 
                                  summary.smoothed.NL, 
                                  summary.smoothed.DK)
rm(summary.smoothed.BE, summary.smoothed.DK, summary.smoothed.NL )

setwd('C:/R/master thesis/HLW_edited_code/code_thesis_DEJEAN_1303')

use_kappa_int = 1

#R* VS R, 3 countries facet plot
plot_df_IR <- summary.smoothed.all %>%
  filter(state == "rstar_t")

file_path <- "inputData/DEJEAN_thesis_data_ER.xlsx"
sheets <- list(BE = "BE1971Q12024Q4",
               NL = "NL1971Q12024Q4",
               DK = "DK1971Q12024Q4")

get_expost_IR <- function(country, sheet, file, time_vector) {
  raw_data <- read.xlsx(file, sheet = sheet, 
                        na.strings = ".", 
                        colNames = TRUE, 
                        rowNames = FALSE, 
                        detectDates = TRUE)
  raw_expost <- raw_data$interest - raw_data$inflation
  expost <- raw_expost[5:length(raw_expost)]
  tibble(
    time = time_vector,
    expost_real_IR = expost,
    Country = country
  )
}

expost_IR_all <- imap_dfr(sheets, ~get_expost_IR(.y, .x, file_path, time.vector))%>%
  mutate(time = time[[1]]) 

plot_df_IR <- left_join(plot_df_IR, expost_IR_all, by = c("time", "Country"))
plot_df_IR <- plot_df_IR %>% 
  mutate(rgaplow = expost_real_IR - upper, rgap = expost_real_IR - mean, rgapup = expost_real_IR - lower)

plot_IR_compared <- ggplot(plot_df_IR, aes(x = time, colour=Country)) +
  geom_line(aes(y = mean), linewidth = 1) +  # Bold mean line
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) + ylim(-4, 8) +
  theme_linedraw(base_size = 18) + ylab("r* (%)") + 
  theme(axis.title.x = element_blank())

plot_IR_vs_rstar <- ggplot(plot_df_IR, aes(x = time)) +
  geom_line(aes(y = expost_real_IR), linewidth = 1) +
   geom_line(aes(y = mean), linewidth = 2.5, col="red") +  # Bold mean line
  geom_line(aes(y = lower), linetype = "dashed", linewidth = 1, col="red") +  # Lower CI as dotted line
  geom_line(aes(y = upper), linetype = "dashed", linewidth = 1, col="red") +  # Upper CI as dotted line
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) + ylim(-10,15) +
labs(y = "% per annum") + theme_linedraw(base_size = 40) +
  facet_wrap(~ Country, ncol=1, scales = "free") +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_IR_gap <- ggplot(plot_df_IR, aes(x = time)) +
  geom_line(aes(y = rgap), linewidth = 1.5, col="red") +  # Bold mean line
  geom_line(aes(y = rgaplow), linetype = "dashed", linewidth = 0.5, col="red") + 
  geom_line(aes(y = rgapup), linetype = "dashed", linewidth = 0.5, col="red") +  
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) +
  labs(y = "% per annum") + theme_linedraw(base_size = 40) + ylim(-10,10) +
  facet_wrap(~ Country, ncol=1, scales = "free") +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#g 
plot_df_g <- summary.smoothed.all %>%
  filter(state == "gann_t")

plot_g_compared <- ggplot(plot_df_g, aes(x = time, colour=Country)) +
  geom_line(aes(y = mean), linewidth = 1, linetype="dashed") + 
  geom_line(data=plot_df_IR, aes(y = mean), linewidth = 1.5) + 
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) + ylab("g (dashed) and r* (full) (%)") + 
  scale_color_manual(values = c("BE" = "darkgreen", "DK"="red","NL" = "black")) +
  scale_y_continuous(breaks=seq(-5,10,2.5)) +
  theme_linedraw(base_size = 40) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#z 
plot_df_z <- summary.smoothed.all %>%
  filter(state == "z_t")

plot_z_compared <- ggplot(plot_df_z, aes(x = time, colour=Country)) +
  geom_line(aes(y = mean), linewidth = 1) + 
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) + ylim(-6, 5) +
  theme_linedraw(base_size = 18) + ylab("z (%)") + 
  theme(axis.title.x = element_blank())

#ER
plot_df_ER <- summary.smoothed.all %>%
  filter(state == "qstar_t")

file_path <- "inputData/DEJEAN_thesis_data_ER.xlsx"
sheets <- list(BE = "BE1971Q12024Q4",
               NL = "NL1971Q12024Q4",
               DK = "DK1971Q12024Q4")

get_expost_ER <- function(country, sheet, file, time_vector) {
  raw_data <- read.xlsx(file, sheet = sheet, 
                        na.strings = ".", 
                        colNames = TRUE, 
                        rowNames = FALSE, 
                        detectDates = TRUE)
  raw_er <- raw_data$exchange.rate 
  er <- raw_er[5:length(raw_er)]
  tibble(
    time = time_vector,
    ER = er,
    Country = country
  )
}

ER_all <- imap_dfr(sheets, ~get_expost_ER(.y, .x, file_path, time.vector))%>%
  mutate(time = time[[1]]) 

plot_df_ER <- left_join(plot_df_ER, ER_all, by = c("time", "Country"))

plot_df_ER <- plot_df_ER %>% 
  mutate(qgaplow = ER - upper, qgap = ER - mean, qgapup = ER - lower)


plot_ER_vs_qstar <- ggplot(plot_df_ER, aes(x = time)) +
  geom_line(aes(y = mean), linewidth = 2.5, col="red") +  # Bold mean line
  geom_line(aes(y = lower), linetype = "dashed", linewidth = 1, col="red") +  # Lower CI as dotted line
  geom_line(aes(y = upper), linetype = "dashed", linewidth = 1, col="red") +  # Upper CI as dotted line
  geom_line(aes(y = ER), linewidth = 1) +
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y")  +
  labs(y = "Index") + theme_linedraw(base_size = 40) +
  facet_wrap(~ Country, ncol=1, scales = "free") +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#output
plot_df_y = summary.smoothed.all %>% filter(state == "ystar_t")

get_y <- function(country, sheet, file, time_vector) {
  raw_data <- read.xlsx(file, sheet = sheet, 
                        na.strings = ".", 
                        colNames = TRUE, 
                        rowNames = FALSE, 
                        detectDates = TRUE)
  raw_y <- raw_data$gdp.log*100 
  y <- raw_y[5:length(raw_y)]
  tibble(
    time = time_vector,
    y = y,
    Country = country
  )
}

Y_all <- imap_dfr(sheets, ~get_y(.y, .x, file_path, time.vector))%>%
  mutate(time = time[[1]]) 

plot_df_y <- left_join(plot_df_y, Y_all, by = c("time", "Country"))

plot_y_vs_ystar <- ggplot(plot_df_y, aes(x = time)) +
  geom_line(aes(y = y), linewidth = 1) +
  geom_line(aes(y = mean), linewidth = 1.5, col="red") + 
  geom_line(aes(y = lower), linetype = "dashed", linewidth = 0.5, col="red") + 
  geom_line(aes(y = upper), linetype = "dashed", linewidth = 0.5, col="red") +  
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  labs(y = "log(output) x 100") + theme_linedraw(base_size = 40) +
  facet_wrap(~ Country, scales="free", ncol=1) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")

#ygap
#g 
plot_df_ygap <- summary.smoothed.all %>%
  filter(state == "ygap_t")

plot_ygap_compared <- ggplot(plot_df_ygap, aes(x = time)) +
  geom_line(aes(y = mean), linewidth = 1.5, col="red") + 
  geom_line(aes(y = lower), linetype = "dashed", linewidth = 0.5, col="red") + 
  geom_line(aes(y = upper), linetype = "dashed", linewidth = 0.5, col="red") +  
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) + 
  labs(y = "Output gap in %") + theme_linedraw(base_size = 40) +
  facet_wrap(~ Country, ncol=1, scales = "free_x") +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#decomposition of rstar
c = c("BE"=1.00, "NL"=1.01, "DK"=0.96) #my estimates
plot_df_dec <- summary.smoothed.all %>%
  filter(state %in% c("z_t", "gann_t")) %>%
  select(-lower, -upper)

plot_df_dec <- plot_df_dec %>% 
  mutate(mean = if_else(state == "gann_t", mean * c[as.character(Country)], mean))

plot_dec <- ggplot(plot_df_dec, aes(x = time, y=mean, fill = state)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_line(data=filter(plot_df_IR, state=="rstar_t") ,aes(y = mean), linewidth = 2.5, col="red") + 
  facet_wrap(~ Country, scales = "free_y") +  
  labs(y = "% per annum",
       fill = "Drivers of r*") +
  theme_linedraw(base_size = 40) +
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2025-01-01") , by = "10 years"), date_labels = "%y") +
  facet_wrap(~ Country, ncol=3, scales = "free_x") +
  scale_fill_manual(values = c("gann_t" = "steelblue", "z_t" = "black"),
                    labels = c("gann_t" = "Trend growth (g)", 
                               "z_t" = "Other factors (z)")) +
  scale_y_continuous(breaks=seq(-5,10,2.5)) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#rgap vs qgap
plot_IRER_gaps <- ggplot(plot_df_IR, aes(x = time)) +
  geom_line(aes(y = rgap), linewidth = 1.5, col="red") +  
  geom_line(data= plot_df_ER, aes(y = qgap), linewidth = 1.5, col="black") +  
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "10 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) +
  labs(y = "real interest rate and exchange rate gaps") + theme_linedraw(base_size = 40) + 
  facet_wrap(~ Country, ncol=3, scales = "free") +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#rgap vs qgap
plot_IRER_gaps <- ggplot(plot_df_IR, aes(x = time)) +
  geom_line(aes(y = rgap), linewidth = 1.5, col="red") +  
  geom_line(data= plot_df_ER, aes(y = qgap), linewidth = 1.5, col="black") +  
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "10 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) +
  labs(y = "real interest rate and exchange rate gaps") + theme_linedraw(base_size = 40) + 
  facet_wrap(~ Country, ncol=3, scales = "free") +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_df_gaps3 <- left_join(plot_df_IR, plot_df_ER, by = c("time", "Country")) %>% 
  select(time, Country, rgap, qgap) %>%
  pivot_longer(cols = c("rgap", "qgap"), 
               names_to = "Gaps", 
               values_to = "value")
plot_df_gaps3$time <- as.Date(plot_df_gaps3$time)

plot_IRER_gaps2 <- ggplot(plot_df_gaps3, aes(x = time, y = value, fill = Gaps)) +
  geom_col(position = position_dodge2(width = 2)) +   # side-by-side bars
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "10 years"), date_labels = "%y") +
  geom_hline(yintercept=0, linetype = "dotted", linewidth = 1) +
  labs(y = "real interest rate (%) and exchange rate (index) gaps") + theme_linedraw(base_size = 30) + 
  facet_wrap(~ Country, ncol=3) + ylim(-10,15)+
  scale_fill_manual(values = c("qgap" = "darkgreen", "rgap" = "darkblue"),
                    labels = c("qgap" = "Exchange rate gap", 
                               "rgap" = "Interest rate gap")) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#debt
raw_debt <- read.xlsx(file_path, sheet = "IMFDEBT", 
                   na.strings = ".", 
                   colNames = TRUE, 
                   rowNames = FALSE, 
                   detectDates = TRUE)
plot_df_debt <- pivot_longer(raw_debt, cols = c(BE, DK, NL),
                             names_to = "Country", values_to = "Debt")

plot_debt <- ggplot(plot_df_debt, aes(x=year, y=Debt, colour=Country)) +
  geom_line(linewidth=2) +  theme_linedraw(base_size = 20) + 
  labs(y = "Debt (% of GDP)") +
  scale_y_continuous(breaks=seq(0,150,25)) +
  scale_color_manual(values = c("BE" = "darkgreen", "DK"="red","NL" = "black")) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#plot data
get_raw_data <- function(country, sheet, file, time_vector) {
  raw_data <- read.xlsx(file, sheet = sheet, 
                        na.strings = ".", 
                        colNames = TRUE, 
                        rowNames = FALSE, 
                        detectDates = TRUE)
  nominal_IR <- raw_data$interest[5:length(raw_data$interest)] 
  inflation <- raw_data$inflation[5:length(raw_data$interest)] 
  gdp <- raw_data$gdp.log[5:length(raw_data$interest)] 
  er <- raw_data$exchange.rate[5:length(raw_data$interest)] 
  tibble(
    time = time_vector,
    nominal_IR = nominal_IR,
    inflation= inflation,
    gdp = gdp,
    er =er,
    Country = country
  )
}

plot_data_raw <- imap_dfr(sheets, ~get_raw_data(.y, .x, file_path, time.vector))%>%
  mutate(time = time[[1]]) 

plot_data_IR <- plot_data_raw %>% 
  select(-inflation, -gdp, -er)
plot_data_inf <- plot_data_raw %>% 
  select(-nominal_IR, -gdp, -er)
plot_data_ER <- plot_data_raw %>% 
  select(-inflation, -gdp, -nominal_IR)
plot_data_output <- plot_data_raw %>% 
  select(-inflation, -nominal_IR, -er)

plot_raw_y <- ggplot(plot_data_output, aes(x = time, colour=Country)) +
  geom_line(aes(y = gdp), linewidth = 1.5) +
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  labs(y = "log(output) x 100") + theme_linedraw(base_size = 40) +
  scale_color_manual(values = c("BE" = "darkgreen", "DK"="red","NL" = "black")) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")

plot_raw_IR <- ggplot(plot_data_IR, aes(x = time, colour=Country)) +
  geom_line(aes(y = nominal_IR), linewidth = 1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  labs(y = "Short-term nominal interest rates (%)") + theme_linedraw(base_size = 40) +
  scale_color_manual(values = c("BE" = "darkgreen", "DK"="red","NL" = "black")) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")

plot_raw_ER <- ggplot(plot_data_ER, aes(x = time, colour=Country)) +
  geom_line(aes(y = er), linewidth = 1.5) +
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  labs(y = "Real Effective Exchange Rate (Index)") + theme_linedraw(base_size = 40) +
  scale_color_manual(values = c("BE" = "darkgreen", "DK"="red","NL" = "black")) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")

plot_raw_inf <- ggplot(plot_data_inf, aes(x = time, colour=Country)) +
  geom_line(aes(y = inflation), linewidth = 1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_date(breaks = seq(as.Date("1975-01-01"),as.Date("2020-01-01") , by = "5 years"), date_labels = "%y") +
  labs(y = "Inflation (QoQ, %)") + theme_linedraw(base_size = 40) +
  scale_color_manual(values = c("BE" = "darkgreen", "DK"="red","NL" = "black")) +
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")
