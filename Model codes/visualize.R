library(tidyverse)

burn_in <- 1000
mcmc_file <- "output_1_1_7_7_1_5_C15_S56_S76.dat"
sim_file <- "output_1_1_7_7_1_5_C15_S56_S76_sims.txt"
trees_file <- "output_1_1_7_7_1_5_C15_S56_S76_trees.csv"

params <- c("phi", "psi", "alphaSZ", "beta_ext", "beta0", "beta1", "beta2",
            "beta3")
cnames <- c("ll", params, str_c("ar_", params))
mcmc_output <- read_delim(mcmc_file, delim = " ", col_names = cnames) %>%
  mutate(iter = 1:n())

estimates <- mcmc_output[ , c("iter", params)] %>%
  filter(iter > burn_in) %>%
  gather(-iter, key = "variable", value = "value") %>%
  group_by(variable) %>%
  summarize(mean = mean(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))

# For paper
inverse_estimates <- mcmc_output[ , c("iter", "phi", "psi", "alphaSZ")] %>%
  filter(iter > burn_in) %>%
  mutate(inverse_phi = 1 / phi,
         inverse_psi = 1 / psi,
         inverse_alphaSZ = 1 / alphaSZ) %>%
  gather(-iter, key = "variable", value = "value") %>%
  group_by(variable) %>%
  summarize(mean = mean(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))

################################################################################
# Epi curve
################################################################################
# Load real data
cnames <- c("id", "department", "insee", "type", "x", "y", "date_detection",
            "date_culling", "symptomatic", "preventively_culled", "delay_empty",
            "op", "is_pag", "in_zrp")
data <- read_delim("data.txt", delim = " ", col_names = cnames)

cnames <- c("sim", "id", "department", "insee", "date_detection",
            "date_culling", "symptomatic", "preventively_culled")
sims <- read_delim(sim_file, delim = " ", col_names = cnames)
tot_sims <- max(sims$sim)

get_real_epi <- function(data, dep = -1) {
  df <- data
  if (dep > 0) {
    df <- df %>%
      filter(department == dep)
  }
  df %>%
    filter(date_detection > -1) %>%
    group_by(date_detection) %>%
    summarize(n_infections = n())
}

real_epi <- get_real_epi(data)
real_epi$date <- as.Date("2016-11-28") + real_epi$date_detection
real_epi_40 <- get_real_epi(data, 40)
real_epi_40$date <- as.Date("2016-11-28") + real_epi_40$date_detection

################################################################################
# Figure 1
################################################################################
breaks <- seq(as.Date("2016-12-01"), as.Date("2017-04-01"), by = "1 month")
labels <- c("12/16", "01/17", "02/17", "03/17", "04/17")

font_size <- 16
ggplot(real_epi) +
  geom_vline(xintercept = as.Date("2017-01-04"), linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = as.Date("2017-02-02"), linetype = "dashed", size = 0.3) +
  geom_bar(aes(x = date, y = n_infections), stat = "identity") +
  scale_x_date("", breaks = breaks, labels = labels) +
  scale_y_continuous("Number of IPs", breaks = seq(0, 25, by = 5),
                     limits = c(0, 25)) +
  theme_bw() +
  theme(text = element_text(size = font_size))

################################################################################
# Figure 3
################################################################################
xs <- seq(min(real_epi$date), max(real_epi$date), by = "1 day")
s1 <- as.Date("2016-11-28") + 56
s2 <- as.Date("2016-11-28") + 76
ys <- rep(estimates$mean[which(estimates$variable == "beta3")], length(xs))
ys[xs < s1] <- estimates$mean[which(estimates$variable == "beta1")]
ys[xs >= s1 & xs < s2] <- estimates$mean[which(estimates$variable == "beta2")]
ys_lower <- rep(estimates$lower[which(estimates$variable == "beta3")],
                length(xs))
ys_lower[xs < s1] <- estimates$lower[which(estimates$variable == "beta1")]
ys_lower[xs >= s1 & xs < s2] <- estimates$lower[
  which(estimates$variable == "beta2")]
ys_upper <- rep(estimates$upper[which(estimates$variable == "beta3")],
                length(xs))
ys_upper[xs < s1] <- estimates$upper[which(estimates$variable == "beta1")]
ys_upper[xs >= s1 & xs < s2] <- estimates$upper[
  which(estimates$variable == "beta2")]
y0s <- rep(estimates$mean[which(estimates$variable == "beta0")], length(xs))
y0s_lower <- rep(estimates$lower[which(estimates$variable == "beta0")],
                 length(xs))
y0s_upper <- rep(estimates$upper[which(estimates$variable == "beta0")],
                 length(xs))
df <- tibble(
  date = xs,
  beta = ys,
  beta_lower = ys_lower,
  beta_upper = ys_upper,
  beta0 = y0s,
  beta0_lower = y0s_lower,
  beta0_upper = y0s_upper
) %>% left_join(real_epi_40, by = "date")

breaks <- seq(as.Date("2016-12-01"), as.Date("2017-04-01"), by = "1 month")
labels <- c("12/16", "01/17", "02/17", "03/17", "04/17")

font_size <- 16
ggplot(df) +
  geom_ribbon(aes(x = date, ymin = beta_lower, ymax = beta_upper),
              fill = "dodgerblue3", alpha = 0.3) +
  geom_line(aes(x = date, y = beta), color = "dodgerblue4") +
  geom_line(aes(x = date, y = beta0)) +
  geom_line(aes(x = date, y = beta0_lower), linetype = "dashed") +
  geom_line(aes(x = date, y = beta0_upper), linetype = "dashed") +
  scale_x_date("", breaks = breaks, labels = labels) +
  scale_y_continuous(expression(beta)) +
  theme_bw() +
  theme(text = element_text(size = font_size))

################################################################################
# Figure 5A
################################################################################
inflate <- function(vec, n_sims) {
  inflated <- rep(0, n_sims)
  inflated[1:length(vec)] <- vec
  return(inflated)
}

sim_epi <- sims %>%
  filter(date_detection %in% 0:max(real_epi$date_detection)) %>%
  group_by(date_detection, sim) %>%
  summarize(n_infections = n()) %>%
  summarize(mean_n_infections = mean(inflate(n_infections, tot_sims)),
            median_n_infections = median(inflate(n_infections, tot_sims)),
            lower_n_infections = quantile(inflate(n_infections, tot_sims),
                                          probs = 0.025),
            upper_n_infections = quantile(inflate(n_infections, tot_sims),
                                          probs = 0.975))
sim_epi$date <- as.Date("2016-11-28") + sim_epi$date_detection

breaks <- seq(as.Date("2016-12-01"), as.Date("2017-04-01"), by = "1 month")
labels <- c("12/16", "01/17", "02/17", "03/17", "04/17")

font_size <- 16
ggplot() +
  geom_ribbon(data = sim_epi, aes(x = date, ymin = lower_n_infections,
                                  ymax = upper_n_infections),
              fill = "#54278f", alpha = 0.3) +
  geom_line(data = sim_epi, aes(x = date, y = mean_n_infections),
            color = "#54278f") +
  geom_point(data = real_epi, aes(x = date, y = n_infections)) +
  scale_x_date("", breaks = breaks, labels = labels) +
  scale_y_continuous("Number of IPs", breaks = seq(0, 25, by = 5)) +
  theme_bw() +
  theme(text = element_text(size = font_size))

################################################################################
# Figure 5B
################################################################################
get_n_cases_by_department <- function(df, departments) {
  # Remove no infections
  tmp <- df[df$date_detection >= 0, ]
  res <- tibble(
    department = departments,
    n_cases = 0)
  for (curr in 1:length(departments)) {
    curr_dep <- departments[curr]
    res$n_cases[curr] <- sum(tmp$department == curr_dep)
  }
  return(res)
}

departments <- unique(data$department)

df_real <- get_n_cases_by_department(data, departments)
df_sim <- sims %>%
  filter(date_detection %in% 0:max(real_epi$date_detection)) %>%
  group_by(department, sim) %>%
  summarize(n_infections = sum(date_detection > -1)) %>%
  summarize(mean_n_infections = mean(inflate(n_infections)),
            median_n_infections = median(inflate(n_infections, tot_sims)),
            lower_n_infections = quantile(inflate(n_infections, tot_sims),
                                          probs = 0.025),
            upper_n_infections = quantile(inflate(n_infections, tot_sims),
                                          probs = 0.975))

df <- df_real %>%
  left_join(df_sim, by = "department") %>%
  arrange(desc(n_cases)) %>%
  mutate(department = factor(department, levels = department))

font_size <- 16
ggplot(df) +
  geom_pointrange(aes(x = department, y = mean_n_infections,
                      ymin = lower_n_infections, ymax = upper_n_infections,
                      color = department)) +
  geom_point(aes(x = department, y = n_cases), shape = 15, size = 2.5) +
  xlab("Department") +
  scale_y_continuous("Number of IPs") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = font_size))

################################################################################
# Transmission trees
################################################################################
t_trees <- read_csv(trees_file)

tTrees <- t_trees
# R and d over time (smoothed)
the_range <- range(tTrees$InfectionDate[!is.na(tTrees$InfectionDate)])
n_trees <- max(tTrees$Sim)
infection_dates <- seq(the_range[1], the_range[2], by = 1)

infected <- data[data$date_detection > -1, ]
infectious_at <- data.frame(
  ID = rep(infected$id, each = length(infection_dates)),
  Date = rep(infection_dates, nrow(infected)),
  Infectious = F)
latence <- 1
incubation <- 7
delayDetection <- 1
delayAsymptomatic <- 5
for (ID in infected$id) {
  ind_old <- which(infected$id == ID)
  if (length(ind_old) > 1) {
    print("Wooooo")
  }
  if (infected$symptomatic[ind_old] == 1) {
    delay <- incubation + delayDetection - latence
  } else {
    delay <- delayAsymptomatic - latence
  }
  start <- infected$date_detection[ind_old] - delay
  stop <- infected$date_culling[ind_old] - 1
  inds <- infectious_at$ID == ID & infectious_at$Date >= start &
    infectious_at$Date <= stop
  infectious_at$Infectious[inds] <- T
}

res_R_d <- data.frame(
  Sim = rep(1:n_trees, each = length(infection_dates)),
  Date = rep(infection_dates, n_trees),
  R = NA,
  Distance = NA
)

half_window <- 3
for (sim in unique(tTrees$Sim)) {
  print(sim)
  tmp <- tTrees[tTrees$Sim == sim & !is.na(tTrees$Infector), ]
  for (curr_date in infection_dates) {
    lower <- curr_date - half_window
    upper <- curr_date + half_window
    in_window <- tmp$InfectionDate >= lower & tmp$InfectionDate <= upper
    
    ind <- which(res_R_d$Sim == sim & res_R_d$Date == curr_date)
    if (length(ind) > 1) {
      print("Oh-oh!")
    }
    
    # Retrieve all infectious in time window
    infectiousIDs <- unique(infectious_at$ID[
      infectious_at$Date >= (curr_date - half_window) &
        infectious_at$Date <= (curr_date + half_window) &
        infectious_at$Infectious])
    
    # Compute number infected by each infector
    tmp2 <- data.frame(ID = infectiousIDs, NInfected = 0)
    for (ID in infectiousIDs) {
      tmp2$NInfected[tmp2$ID == ID] <- sum(tmp$Infector == ID)
    }
    
    # R = mean number of infections for each infectious individual
    res_R_d$R[ind] <- mean(tmp2$NInfected)
    res_R_d$Distance[ind] <- mean(tmp$Distance[in_window])
  }
}
rm(tmp, tmp2)

tmp <- res_R_d %>%
  group_by(Date) %>%
  summarize(MeanR = mean(R),
            LowerR = quantile(R, probs = 0.025), 
            UpperR = quantile(R, probs = 0.975),
            MeanD = mean(Distance),
            LowerD = quantile(Distance, probs = 0.025),
            UpperD = quantile(Distance, probs = 0.975))

to_plot <- tmp %>%
  filter(Date >= 0 & Date <= max(data$date_detection)) %>%
  mutate(date = Date + as.Date("2016-11-28"))

breaks <- seq(as.Date("2016-12-01"), as.Date("2017-04-01"), by = "1 month")
labels <- c("12/16", "01/17", "02/17", "03/17", "04/17")
font_size <- 16
font_size_label <- 10

################################################################################
# Figure 4
################################################################################
ggplot(to_plot) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  geom_ribbon(aes(x = date, ymin = LowerR, ymax = UpperR),
              fill = "indianred3", alpha = 0.2) +
  geom_line(aes(x = date, y = MeanR), color = "indianred4") +
  annotate("text", x = as.Date("2016-12-01"), y = 2.48, label = "A",
           size = font_size_label) +
  scale_x_date("", breaks = breaks, labels = labels) +
  scale_y_continuous("Reproduction Number", breaks = seq(0, 3.0, by = 0.5),
                     limits = c(0, 2.5)) +
  theme_bw() +
  theme(text = element_text(size = font_size))


ggplot(to_plot) +
  geom_hline(yintercept = 8.49, linetype = "dashed") +
  geom_ribbon(aes(x = date, ymin = LowerD, ymax = UpperD),
              fill = "dodgerblue3", alpha = 0.2) +
  geom_line(aes(x = date, y = MeanD), color = "dodgerblue4") +
  annotate("text", x = as.Date("2016-12-01"), y = 13.8, label = "B",
           size = font_size_label) +
  scale_x_date("", breaks = breaks, labels = labels) +
  scale_y_continuous("Infection Distance (km)",
                     breaks = seq(0, 15, by = 2), limits = c(0, 14)) +
  theme_bw() +
  theme(text = element_text(size = font_size))
