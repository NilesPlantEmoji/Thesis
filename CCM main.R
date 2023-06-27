###This is the code used for the analysis in the bachelor's thesis of Niels Kaptein###
###Please note that this code uses an older version of the "rEDM" package"
###The data used can be found on PANGEA. See thesis for the links to download them

#Packages
library(ggplot2)
library(remotes)
require(remotes)
install_version("rEDM", version = "0.7.3", repos = "http://cran.us.r-project.org")
library(rEDM)
library(latex2exp)
library(ggpubr)
library(ggpmisc)
library(zoo)

#Reading data. These files are cleaned up versions of the original datasets, removing all text so that R can read it.
methane_data <- read.table("ch4nat-noaa_CLEAN.txt", sep = "\t", header = TRUE)
CO2_data <- read.table("Vostok_CH4_CO2_age_CLEAN.tab", sep = "\t", header = TRUE)
d18O_data <- read.table("lr04_CLEAN.txt", sep = "\t", header = TRUE)

#Converting from years to ky
methane_data$age_ka <- methane_data$gas_ageBP/1000

#Synchronising all time series 
time_axis <- seq(1, 417, by = 1)
CH4_interp <- approx(methane_data$age_ka, methane_data$CH4, xout = time_axis, method = "linear", rule = 2)$y
CO2_interp <- approx(CO2_data$Age..ka.BP., CO2_data$CO2..ppmv., xout = time_axis, method = "linear", rule = 2)$y
d18O_interp <- approx(d18O_data$Time..ka., d18O_data$Benthic.d18O..per.mil., xout = time_axis, 
                      method = "linear", rule = 2)$y
combined_data <- data.frame(time_axis = time_axis, d18O_interp, CH4_interp, CO2_interp)

#Plotting time series
ggplot(combined_data, aes(x = time_axis)) +
  geom_line(aes(y = CO2_interp), color = "blue") +
  scale_x_reverse() +
  labs(x = "Time (ka)", y = TeX("$CO_{2}$ (ppmv)")) +
  theme_minimal()
ggsave("CO2_ts.png")

methane_ts <- ggplot(combined_data, aes(x = time_axis)) +
  geom_line(aes(y = CH4_interp), color = "green") +
  scale_x_reverse() +
  labs(x = "Time (ka)", y = TeX("$CH_{4}$ (ppbv)")) +
  theme_minimal()
d18O_ts <- ggplot(combined_data, aes(x = time_axis)) +
  geom_line(aes(y = d18O_interp), color = "darkorchid2") +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Time (ka)", y = TeX("$\\delta^{18}O$ (‰)")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
ggarrange(d18O_ts, methane_ts, ncol = 1, nrow = 2, labels = c("a", "b"))
ggsave("Time_Series.png")

#Simplex projection
ts1 <- d18O_interp
lib1 <- c(1, length(ts1))
pred1 <- c(1, length(ts1))
simplex_output1 <- simplex(ts1, lib1, pred1, silent = TRUE)

ts2 <- CH4_interp
lib2 <- c(1, length(ts2))
pred2 <- c(1, length(ts2))
simplex_output2 <- simplex(ts2, lib2, pred2, silent = TRUE)

ts3 <- CO2_interp
lib3 <- c(1, length(ts3))
pred3 <- c(1, length(ts3))
simplex_output3 <- simplex(ts3, lib3, pred3, silent = TRUE)

simplex_d18O <- ggplot(simplex_output1, aes(x = E, y = rho)) +
  geom_line() +
  scale_x_continuous(breaks = seq(min(simplex_output1$E), max(simplex_output1$E), by = 1),
                     limits = c(min(simplex_output1$E), max(simplex_output1$E))) +
  labs(x = "E", y = TeX("$\\rho$")) +
  theme_minimal()
simplex_CH4 <- ggplot(simplex_output2, aes(x = E, y = rho)) +
  geom_line() +
  scale_x_continuous(breaks = seq(min(simplex_output2$E), max(simplex_output2$E), by = 1),
                     limits = c(min(simplex_output2$E), max(simplex_output2$E))) +
  theme_minimal() +
  labs(x = "E", y = TeX("$\\rho$")) 
ggarrange(simplex_d18O, simplex_CH4, ncol = 2, nrow = 1, labels = c("a", "b"))
ggsave("Simplex.png")

ggplot(simplex_output3, aes(x = E, y = rho)) +
  geom_line() +
  scale_x_continuous(breaks = seq(min(simplex_output3$E), max(simplex_output3$E), by = 1),
                     limits = c(min(simplex_output3$E), max(simplex_output3$E))) +
  labs(x = "E", y = TeX("$\\rho$")) +
  theme_minimal()
ggsave("SimplexCO2.png")

#Main CCM analysis
max_lib_size <- seq(10, 100, by = 5)
d18O_xmap_CH4 <- ccm_means(ccm(combined_data, E = 8, random_libs = TRUE, lib_column = "d18O_interp", 
                               target_column = "CH4_interp", lib_sizes = max_lib_size, num_samples = 300, 
                               silent = TRUE))
CH4_xmap_d18O <- ccm_means(ccm(combined_data, E = 8, random_libs = TRUE, lib_column = "CH4_interp",
                               target_column = "d18O_interp", lib_sizes = max_lib_size, num_samples = 300,
                               silent = TRUE))

d18O_xmap_CO2 <- ccm_means(ccm(combined_data, E = 8, random_libs = TRUE, lib_column = "d18O_interp",
                               target_column = "CO2_interp", lib_sizes = max_lib_size, num_samples = 300,
                               silent = TRUE))
CO2_xmap_d18O <- ccm_means(ccm(combined_data, E = 8, random_libs = TRUE, lib_column = "CO2_interp",
                               target_column = "d18O_interp", lib_sizes = max_lib_size, num_samples = 300,
                               silent = TRUE))

CH4_xmap_CO2 <- ccm_means(ccm(combined_data, E = 5, random_libs = TRUE, lib_column = "CH4_interp",
                              target_column = "CO2_interp", lib_sizes = max_lib_size, num_samples = 300,
                              silent = TRUE))
CO2_xmap_CH4 <- ccm_means(ccm(combined_data, E = 5, random_libs = TRUE, lib_column = "CO2_interp",
                              target_column = "CH4_interp", lib_sizes = max_lib_size, num_samples = 300,
                              silent = TRUE))

#CCM on surrogate datasets
num_surr <- 100
surr_d18O <- make_surrogate_data(combined_data$d18O_interp, method = "ebisuzaki", 
                                 num_surr = num_surr)
surr_CH4 <- make_surrogate_data(combined_data$CH4_interp, method = "ebisuzaki", 
                                num_surr = num_surr)
surr_CO2 <- make_surrogate_data(combined_data$CO2_interp, method = "ebisuzaki",
                                num_surr = num_surr)

ccm_rho_surr <- data.frame(d18O = numeric(num_surr), CH4 = numeric(num_surr))

d18O_xmap_surrCH4 <- data.frame(max_lib_size)
column_names <- paste0("surr", 1:num_surr)
for (col_name in column_names) {
  d18O_xmap_surrCH4[[col_name]] <- NA
}
for (i in 1:num_surr) {
  result <- ccm_means(ccm(cbind(combined_data$d18O_interp, surr_CH4[, i]), E = 8,random_libs = TRUE, 
                          lib_column = 1, target_column = 2, lib_sizes = max_lib_size, num_samples = 300, 
                          silent = TRUE))
  d18O_xmap_surrCH4[, i+1] <- result$rho
}
d18O_xmap_surrCH4$CI_upper <- apply(d18O_xmap_surrCH4[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.975))
d18O_xmap_surrCH4$CI_lower <- apply(d18O_xmap_surrCH4[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.025))


CH4_xmap_surrd18O <- data.frame(lib_size = max_lib_size)
column_names <- paste0("surr", 1:num_surr)
for (col_name in column_names) {
  CH4_xmap_surrd18O[[col_name]] <- NA
}
for (i in 1:num_surr) {
  result <- ccm_means(ccm(cbind(combined_data$CH4_interp, surr_d18O[, i]), E = 8, 
                          random_libs = TRUE, lib_column = 1, target_column = 2, 
                          lib_sizes = max_lib_size, num_samples = 300, silent = TRUE))
  CH4_xmap_surrd18O[, i+1] <- result$rho
}
CH4_xmap_surrd18O$CI_upper <- apply(CH4_xmap_surrd18O[, 2:num_surr+1], 1, function(row) quantile(row, probs = 0.975))
CH4_xmap_surrd18O$CI_lower <- apply(CH4_xmap_surrd18O[, 2:num_surr+1], 1, function(row) quantile(row, probs = 0.025))

d18O_xmap_surrCO2 <- data.frame(max_lib_size)
column_names <- paste0("surr", 1:num_surr)
for (col_name in column_names) {
  d18O_xmap_surrCO2[[col_name]] <- NA
}
for (i in 1:num_surr) {
  result <- ccm_means(ccm(cbind(combined_data$d18O_interp, surr_CO2[, i]), E = 8,random_libs = TRUE, 
                          lib_column = 1, target_column = 2, lib_sizes = max_lib_size, num_samples = 300, 
                          silent = TRUE))
  d18O_xmap_surrCO2[, i+1] <- result$rho
}
d18O_xmap_surrCO2$CI_upper <- apply(d18O_xmap_surrCO2[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.975))
d18O_xmap_surrCO2$CI_lower <- apply(d18O_xmap_surrCO2[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.025))

CO2_xmap_surrd18O <- data.frame(max_lib_size)
column_names <- paste0("surr", 1:num_surr)
for (col_name in column_names) {
  CO2_xmap_surrd18O[[col_name]] <- NA
}
for (i in 1:num_surr) {
  result <- ccm_means(ccm(cbind(combined_data$CO2_interp, surr_d18O[, i]), E = 8,random_libs = TRUE, 
                          lib_column = 1, target_column = 2, lib_sizes = max_lib_size, num_samples = 300, 
                          silent = TRUE))
  CO2_xmap_surrd18O[, i+1] <- result$rho
}
CO2_xmap_surrd18O$CI_upper <- apply(CO2_xmap_surrd18O[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.975))
CO2_xmap_surrd18O$CI_lower <- apply(CO2_xmap_surrd18O[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.025))

CO2_xmap_surrCH4 <- data.frame(max_lib_size)
column_names <- paste0("surr", 1:num_surr)
for (col_name in column_names) {
  CO2_xmap_surrCH4[[col_name]] <- NA
}
for (i in 1:num_surr) {
  result <- ccm_means(ccm(cbind(combined_data$CO2_interp, surr_CH4[, i]), E = 5,random_libs = TRUE, 
                          lib_column = 1, target_column = 2, lib_sizes = max_lib_size, num_samples = 300, 
                          silent = TRUE))
  CO2_xmap_surrCH4[, i+1] <- result$rho
}
CO2_xmap_surrCH4$CI_upper <- apply(CO2_xmap_surrCH4[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.975))
CO2_xmap_surrCH4$CI_lower <- apply(CO2_xmap_surrCH4[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.025))

CH4_xmap_surrCO2 <- data.frame(max_lib_size)
column_names <- paste0("surr", 1:num_surr)
for (col_name in column_names) {
  CH4_xmap_surrCO2[[col_name]] <- NA
}
for (i in 1:num_surr) {
  result <- ccm_means(ccm(cbind(combined_data$CH4_interp, surr_CO2[, i]), E = 5,random_libs = TRUE, 
                          lib_column = 1, target_column = 2, lib_sizes = max_lib_size, num_samples = 300, 
                          silent = TRUE))
  CH4_xmap_surrCO2[, i+1] <- result$rho
}
CH4_xmap_surrCO2$CI_upper <- apply(CH4_xmap_surrCO2[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.975))
CH4_xmap_surrCO2$CI_lower <- apply(CH4_xmap_surrCO2[,2:num_surr+1], 1, function(row) quantile(row, probs = 0.025))

#Combining CCM results
xmap_combined <- rbind(d18O_xmap_CH4, CH4_xmap_d18O)
xmap_combined$direction <- ifelse(xmap_combined$lib_column == "d18O_interp", "d18O xmap to CH4", "CH4 xmap to d18O")
xmap_combined$CI_upper <- ifelse(xmap_combined$lib_column == "d18O_interp", d18O_xmap_surrCH4$CI_upper, CH4_xmap_surrd18O$CI_upper)
xmap_combined$CI_lower <- ifelse(xmap_combined$lib_column == "d18O_interp", d18O_xmap_surrCH4$CI_lower, CH4_xmap_surrd18O$CI_lower)

xmap_combined2 <- rbind(d18O_xmap_CO2, CO2_xmap_d18O)
xmap_combined2$direction <- ifelse(xmap_combined$lib_column == "d18O_interp", "d18O xmap to CO2", "CO2 xmap to d18O")
xmap_combined2$CI_upper <- ifelse(xmap_combined$lib_column == "d18O_interp", d18O_xmap_surrCO2$CI_upper, CO2_xmap_surrd18O$CI_upper)
xmap_combined2$CI_lower <- ifelse(xmap_combined$lib_column == "d18O_interp", d18O_xmap_surrCO2$CI_lower, CO2_xmap_surrd18O$CI_lower)

xmap_combined3 <- rbind(CH4_xmap_CO2, CO2_xmap_CH4)
xmap_combined3$direction <- ifelse(xmap_combined$lib_column == "CH4_interp", "CH4 xmap to CO2", "CO2 xmap to CH4")
xmap_combined3$CI_upper <- ifelse(xmap_combined$lib_column == "CH4_interp", CH4_xmap_surrCO2$CI_upper, CO2_xmap_surrCH4$CI_upper)
xmap_combined3$CI_lower <- ifelse(xmap_combined$lib_column == "CH4_interp", CH4_xmap_surrCO2$CI_lower, CO2_xmap_surrCH4$CI_lower)

#Plotting CCM results
ggplot() +
  geom_line(data = d18O_xmap_CH4, aes(x = max_lib_size, y = rho, color = "d18O xmap to CH4")) +
  geom_line(data = CH4_xmap_d18O, aes(x = max_lib_size, y = rho, color = "CH4 xmap to d18O")) +
  labs(x = "Library Size", y = TeX("$\\rho$"), color = "Direction") +
  scale_color_manual(values = c("red", "blue")) +
  guides(color = guide_legend()) +
  theme_minimal()
  

ggplot(xmap_combined, aes(x = lib_size, y = rho, color = direction, fill = direction)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = "Library size", y = TeX("$\\rho$")) +
  theme(legend.justification = "right") +
  theme_minimal()
ggsave("CCM.png")

CCM2 <- ggplot(xmap_combined2, aes(x = lib_size, y = rho, color = direction, fill = direction)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = "Library size", y = TeX("$\\rho$")) +
  theme(legend.justification = "right") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme_minimal()
CCM3 <- ggplot(xmap_combined3, aes(x = lib_size, y = rho, color = direction, fill = direction)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = "Library size", y = TeX("$\\rho$")) +
  theme(legend.justification = "right") +
  theme_minimal()
ggarrange(CCM2, CCM3, ncol = 1, nrow = 2, labels = c("a", "b"))
ggsave("CCM_supp.png")

#Post-hoc time lag analysis
vars <- names(combined_data)[2:3]
params <- expand.grid(lib_column = vars, target_column = vars, tp = -10:10)
params <- params[params$lib_column != params$target_column, ]
params$E <- 8
output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
  ccm_means(ccm(combined_data, E = params$E[i], lib_sizes = NROW(combined_data), 
                random_libs = FALSE, lib_column = params$lib_column[i], target_column = params$target_column[i], 
                tp = params$tp[i], silent = TRUE))
}))
output$direction <- paste(output$lib_column, "xmap to\n", output$target_column)

vars2 <- names(combined_data)[3:4]
params2 <- expand.grid(lib_column = vars2, target_column = vars2, tp = -10:10)
params2 <- params2[params2$lib_column != params2$target_column, ]
params2$E <- 5
output2 <- do.call(rbind, lapply(seq_len(NROW(params2)), function(i) {
  ccm_means(ccm(combined_data, E = params2$E[i], lib_sizes = NROW(combined_data), 
                random_libs = FALSE, lib_column = params2$lib_column[i], target_column = params2$target_column[i], 
                tp = params2$tp[i], silent = TRUE))
}))
output2$direction <- paste(output2$lib_column, "xmap to\n", output2$target_column)

vars3 <- names(combined_data)[c(2, 4)]
params3 <- expand.grid(lib_column = vars3, target_column = vars3, tp = -10:10)
params3 <- params3[params3$lib_column != params3$target_column, ]
params3$E <- 8
output3 <- do.call(rbind, lapply(seq_len(NROW(params3)), function(i) {
  ccm_means(ccm(combined_data, E = params3$E[i], lib_sizes = NROW(combined_data), 
                random_libs = FALSE, lib_column = params3$lib_column[i], target_column = params3$target_column[i], 
                tp = params3$tp[i], silent = TRUE))
}))
output3$direction <- paste(output3$lib_column, "xmap to\n", output3$target_column)

#Plotting post hoc results
ggplot(output, aes(x = tp, y = rho, color = direction)) + 
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Time lag", y = TeX("$\\rho$")) +
  theme_minimal()
ggsave("post-hoc.png")

phsupp1 <- ggplot(output2, aes(x = tp, y = rho, color = direction)) + 
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Time lag", y = TeX("$\\rho$")) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.2), legend.justification = c(1, 1))
phsupp2 <- ggplot(output3, aes(x = tp, y = rho, color = direction)) + 
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Time lag", y = TeX("$\\rho$")) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.2), legend.justification = c(1, 1))
ggarrange(phsupp2, phsupp1, ncol = 2, nrow = 1, labels = c("a", "b"))
ggsave("post-hoc_supp.png")

#Moving window correlation calculation
r_squared_calculation <- function(x, y) {
  r_squared <- cor(x, y)^2
  return(r_squared)
}
window_size <- 50
correlation_series <- rollapply(data = combined_data[, 2:3], width = window_size, 
                                FUN = function(x) r_squared_calculation(x[, 1], x[, 2]), 
                                by.column = FALSE, align = "left", fill = NA)
correlation_ts <- data.frame(time_axis, correlation_series)

#Moving window dynamic CCM
d18O_xmap_CH4_ts <- data.frame(time_axis, rho = NA)
for (i in 1:(NROW(time_axis) - window_size)) {
  result <- ccm_means(ccm(combined_data[i:(i + window_size), ], E = 8, random_libs = TRUE, lib_column = "d18O_interp",
                          target_column = "CH4_interp", lib_sizes = window_size, num_samples = 300, silent = TRUE))
  d18O_xmap_CH4_ts$rho[i] <- result$rho
}

CH4_xmap_d18O_ts <- data.frame(time_axis, rho = NA)
for (i in 1:(NROW(time_axis) - window_size)) {
  result <- ccm_means(ccm(combined_data[i:(i + window_size), ], E = 8, random_libs = TRUE, lib_column = "CH4_interp", 
                          target_column = "d18O_interp", lib_sizes = window_size, num_samples = 300, silent = TRUE))
  CH4_xmap_d18O_ts$rho[i] <- result$rho
}

#Plotting dynamic CCM (including alternative CH4 time series for the combination plot)
methane_ts2 <- ggplot(combined_data, aes(x = time_axis)) +
  geom_line(aes(y = CH4_interp), color = "green") +
  scale_x_reverse() +
  labs(x = "Time (ka)", y = TeX("$CH_{4}$ (ppbv)")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text = element_blank())
CCM_ts <- ggplot() +
  geom_line(data = d18O_xmap_CH4_ts, aes(x = time_axis, y = rho, color = "d18O xmap to CH4")) +
  geom_line(data = CH4_xmap_d18O_ts, aes(x = time_axis, y = rho, color = "CH4 xmap to d18O")) +
  scale_x_reverse() +
  labs(x = "Time (ka)", y = TeX("$\\rho$"), color = "Direction") +
  scale_color_manual(values = c("red", "blue")) +
  guides(color = guide_legend()) +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.1), legend.justification = c(0, 0))
corr_ts <- ggplot(correlation_ts, aes(x = time_axis, y = correlation_series)) +
  geom_line() + 
  scale_x_reverse() +
  labs(x = "Time (ka)", y = TeX("$R^{2}$")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
ggarrange(d18O_ts, methane_ts2, corr_ts, CCM_ts, ncol = 1, nrow = 4, labels = c("a", "b", "c", "d"))
ggsave("CCM_ts.png")

#Alternative plot formatting for poster
ggplot(output, aes(x = tp, y = rho, color = direction)) + 
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Time lag", y = "Forecast skill") +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.1), legend.justification = c(0, 0))
ggsave("post-hoc_POSTER.png")

ggplot(xmap_combined, aes(x = lib_size, y = rho, color = direction, fill = direction)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = "Library size", y = TeX("$\\rho$")) +
  theme(legend.justification = "right") +
  theme_minimal() +
  theme(legend.position = c(0.95, 0.05), legend.justification = c(1, 0))
ggsave("CCM_POSTER.png")