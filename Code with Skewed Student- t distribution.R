# SP500-45 volatility modeling using GARCH-family models with skewed Student-t innovations.
# The lagged VIX-up dummy is used as an external regressor in each variance equation.
# VIX_up_dummy_t equals 1 when VIX_{t-1} > VIX_{t-2}; otherwise, it equals 0.
# The script evaluates EGARCH, GARCH, and GJR-GARCH models using 1% and 5% VaR and QLIKE.
# All generated outputs are written to output_dir.

# 0. Environment setup
rm(list = ls())
cat("\014")   # Clear the RStudio console

# 1. Output directories
output_dir  <- ""
figures_dir <- file.path(output_dir, "Figures")
tables_dir  <- file.path(output_dir, "Tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

setwd(output_dir)

# 2. Required packages
library(quantmod)
library(fBasics)
library(jtools)
library(stats)
library(forecast)
library(rugarch)
library(xts)
library(PerformanceAnalytics)
library(openxlsx)

# 3. Data collection and return construction
ticker <- "^SP500-45"
vix_ticker <- "^VIX"

raw_data <- getSymbols(
  Symbols = ticker,
  src = "yahoo",
  from = "2010-01-01",
  to = "2024-12-31",
  periodicity = "daily",
  auto.assign = FALSE
)

vix_data <- getSymbols(
  Symbols = vix_ticker,
  src = "yahoo",
  from = "2010-01-01",
  to = "2024-12-31",
  periodicity = "daily",
  auto.assign = FALSE
)

prices <- Ad(raw_data)
colnames(prices) <- "Adjusted"

vix <- Ad(vix_data)
colnames(vix) <- "VIX"

# Compute daily log returns in percentage terms.
logprices <- log(prices)
returns <- diff(logprices) * 100
returns <- na.omit(returns)
colnames(returns) <- "Returns"

write.xlsx(
  data.frame(Date = index(prices), coredata(prices)),
  file = file.path(tables_dir, "prices_SP500_45.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(returns), coredata(returns)),
  file = file.path(tables_dir, "returns_SP500_45.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(vix), coredata(vix)),
  file = file.path(tables_dir, "VIX_levels.xlsx"),
  rowNames = FALSE
)

# 4. Lagged VIX-up dummy regressor
# VIX_up_dummy_t equals 1 when VIX_{t-1} > VIX_{t-2}; otherwise, it equals 0.
merged_data <- na.omit(merge(returns, vix))
colnames(merged_data) <- c("Returns", "VIX")

# Create the VIX lags used to define the dummy variable.
merged_data$VIX_lag1 <- lag(merged_data$VIX, 1)
merged_data$VIX_lag2 <- lag(merged_data$VIX, 2)

# Dummy variable indicating an increase in VIX on the previous trading day.
merged_data$VIX_up_dummy_lag1 <- ifelse(
  merged_data$VIX_lag1 > merged_data$VIX_lag2, 1, 0
)

merged_data <- na.omit(merged_data)

returns_final <- merged_data$Returns
colnames(returns_final) <- "Returns"

vix_dummy_lag1 <- merged_data$VIX_up_dummy_lag1
colnames(vix_dummy_lag1) <- "VIX_up_dummy_lag1"

# External regressor matrix used in the variance equation.
vix_reg_mat <- as.matrix(vix_dummy_lag1)

write.xlsx(
  data.frame(
    Date = index(merged_data),
    Returns = as.numeric(merged_data$Returns),
    VIX = as.numeric(merged_data$VIX),
    VIX_lag1 = as.numeric(merged_data$VIX_lag1),
    VIX_lag2 = as.numeric(merged_data$VIX_lag2),
    VIX_up_dummy_lag1 = as.numeric(merged_data$VIX_up_dummy_lag1)
  ),
  file = file.path(tables_dir, "merged_returns_vix_dummy.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "head_tail_returns.txt"))
cat("HEAD OF FINAL ALIGNED DATA\n")
print(head(merged_data))
cat("\nTAIL OF FINAL ALIGNED DATA\n")
print(tail(merged_data))
sink()

# Descriptive statistics for the final return series.
stats_returns <- basicStats(returns_final)

write.xlsx(
  as.data.frame(stats_returns),
  file = file.path(tables_dir, "basicStats_SP500_45.xlsx"),
  rowNames = TRUE
)

sink(file.path(output_dir, "basicStats_SP500_45.txt"))
cat("Basic Statistics for SP500-45 returns\n\n")
print(stats_returns)
sink()

# 5. In-sample and out-of-sample split
# In-sample period ends on 31 December 2021.
returns_is <- returns_final["/2021-12-31"]

T1 <- NROW(returns_is)
T2 <- NROW(returns_final["2022-01-01/"])

dates_oos <- index(returns_final[(T1 + 1):(T1 + T2)])

sink(file.path(output_dir, "sample_info.txt"))
cat("In-sample observations (T1): ", T1, "\n", sep = "")
cat("Out-of-sample observations (T2): ", T2, "\n\n", sep = "")
cat("Head of in-sample returns:\n")
print(head(returns_is))
cat("\nTail of in-sample returns:\n")
print(tail(returns_is))
sink()

# Out-of-sample evaluation starts on 1 January 2022.
realizedret <- returns_final[(T1 + 1):(T1 + T2)]
colnames(realizedret) <- "Realized_Return"

write.xlsx(
  data.frame(Date = index(realizedret), coredata(realizedret)),
  file = file.path(tables_dir, "realized_returns_oos.xlsx"),
  rowNames = FALSE
)

# 6. Model specifications
# The same lagged VIX-up dummy is included in each variance equation.

# Model 1: EGARCH(1,1) with skewed Student-t innovations.
spec_egarch11_dummy <- ugarchspec(
  variance.model = list(
    model = "eGARCH",
    garchOrder = c(1, 1),
    external.regressors = vix_reg_mat
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "sstd"
)

# Model 2: GARCH(1,1) with skewed Student-t innovations.
spec_garch11_dummy <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1),
    external.regressors = vix_reg_mat
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "sstd"
)

# Model 3: GJR-GARCH(1,1) with skewed Student-t innovations.
spec_gjr11t_dummy <- ugarchspec(
  variance.model = list(
    model = "gjrGARCH",
    garchOrder = c(1, 1),
    external.regressors = vix_reg_mat
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "sstd"
)

# 7. Model estimation
egarch11_dummy <- ugarchfit(
  spec = spec_egarch11_dummy,
  data = returns_final,
  out.sample = T2,
  solver = "hybrid"
)

garch11_dummy <- ugarchfit(
  spec = spec_garch11_dummy,
  data = returns_final,
  out.sample = T2,
  solver = "hybrid"
)

gjr11t_dummy <- ugarchfit(
  spec = spec_gjr11t_dummy,
  data = returns_final,
  out.sample = T2,
  solver = "hybrid"
)

# Save the complete estimation outputs.
sink(file.path(output_dir, "egarch11_dummy_output.txt"))
cat("EGARCH(1,1)-sstd with lagged VIX-up dummy in variance equation\n\n")
show(egarch11_dummy)
sink()

sink(file.path(output_dir, "garch11_dummy_output.txt"))
cat("GARCH(1,1)-sstd with lagged VIX-up dummy in variance equation\n\n")
show(garch11_dummy)
sink()

sink(file.path(output_dir, "gjr11t_dummy_output.txt"))
cat("GJR-GARCH(1,1)-sstd with lagged VIX-up dummy in variance equation\n\n")
show(gjr11t_dummy)
sink()

# 8. Diagnostic plots
# The ACF of squared standardized residuals is exported for each model.
png(file.path(figures_dir, "squared_residual_egarch11_dummy.png"),
    width = 1200, height = 800, res = 150)
plot(egarch11_dummy, which = 11)
dev.off()

png(file.path(figures_dir, "squared_residual_garch11_dummy.png"),
    width = 1200, height = 800, res = 150)
plot(garch11_dummy, which = 11)
dev.off()

png(file.path(figures_dir, "squared_residual_gjr11t_dummy.png"),
    width = 1200, height = 800, res = 150)
plot(gjr11t_dummy, which = 11)
dev.off()

# 9. In-sample model comparison
ic_egarch <- infocriteria(egarch11_dummy)
ic_garch  <- infocriteria(garch11_dummy)
ic_gjr    <- infocriteria(gjr11t_dummy)

coef_egarch <- coef(egarch11_dummy)
coef_garch  <- coef(garch11_dummy)
coef_gjr    <- coef(gjr11t_dummy)

model_comparison <- data.frame(
  Model = c("EGARCH(1,1)",
            "GARCH(1,1)",
            "GJR-GARCH(1,1)"),
  LogLikelihood = c(likelihood(egarch11_dummy),
                    likelihood(garch11_dummy),
                    likelihood(gjr11t_dummy)),
  AIC = c(ic_egarch[1], ic_garch[1], ic_gjr[1]),
  BIC = c(ic_egarch[2], ic_garch[2], ic_gjr[2]),
  Shibata = c(ic_egarch[3], ic_garch[3], ic_gjr[3]),
  HannanQuinn = c(ic_egarch[4], ic_garch[4], ic_gjr[4]),
  Dummy_Coefficient = c(
    ifelse("vxreg1" %in% names(coef_egarch), coef_egarch["vxreg1"], NA),
    ifelse("vxreg1" %in% names(coef_garch),  coef_garch["vxreg1"],  NA),
    ifelse("vxreg1" %in% names(coef_gjr),    coef_gjr["vxreg1"],    NA)
  )
)

write.xlsx(
  model_comparison,
  file = file.path(tables_dir, "model_comparison_insample.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "model_comparison_insample.txt"))
cat("IN-SAMPLE MODEL COMPARISON\n\n")
print(model_comparison)
sink()

# 10. Rolling one-step-ahead forecasts
# VaR is calculated at the 1% and 5% levels.
roll_egarch11_dummy <- ugarchroll(
  spec_egarch11_dummy,
  data = returns_final,
  n.start = T1,
  refit.every = 5,
  refit.window = "moving",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

roll_garch11_dummy <- ugarchroll(
  spec_garch11_dummy,
  data = returns_final,
  n.start = T1,
  refit.every = 5,
  refit.window = "moving",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

roll_gjr11t_dummy <- ugarchroll(
  spec_gjr11t_dummy,
  data = returns_final,
  n.start = T1,
  refit.every = 5,
  refit.window = "moving",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

# 11. Rolling volatility forecasts
vol_egarch11_dummy_roll <- xts(
  roll_egarch11_dummy@forecast$density$Sigma,
  order.by = dates_oos
)
colnames(vol_egarch11_dummy_roll) <- "Vol_EGARCH11_DUMMY"

vol_garch11_dummy_roll <- xts(
  roll_garch11_dummy@forecast$density$Sigma,
  order.by = dates_oos
)
colnames(vol_garch11_dummy_roll) <- "Vol_GARCH11_DUMMY"

vol_gjr11t_dummy_roll <- xts(
  roll_gjr11t_dummy@forecast$density$Sigma,
  order.by = dates_oos
)
colnames(vol_gjr11t_dummy_roll) <- "Vol_GJR11t_DUMMY"

write.xlsx(
  data.frame(Date = index(vol_egarch11_dummy_roll), coredata(vol_egarch11_dummy_roll)),
  file = file.path(tables_dir, "vol_egarch11_dummy_roll.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(vol_garch11_dummy_roll), coredata(vol_garch11_dummy_roll)),
  file = file.path(tables_dir, "vol_garch11_dummy_roll.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(vol_gjr11t_dummy_roll), coredata(vol_gjr11t_dummy_roll)),
  file = file.path(tables_dir, "vol_gjr11t_dummy_roll.xlsx"),
  rowNames = FALSE
)

png(file.path(figures_dir, "vol_egarch11_dummy_roll.png"),
    width = 1200, height = 800, res = 150)
plot(vol_egarch11_dummy_roll, main = "Rolling Volatility Forecast - EGARCH(1,1) ")
dev.off()

png(file.path(figures_dir, "vol_garch11_dummy_roll.png"),
    width = 1200, height = 800, res = 150)
plot(vol_garch11_dummy_roll, main = "Rolling Volatility Forecast - GARCH(1,1) ")
dev.off()

png(file.path(figures_dir, "vol_gjr11t_dummy_roll.png"),
    width = 1200, height = 800, res = 150)
plot(vol_gjr11t_dummy_roll, main = "Rolling Volatility Forecast - GJR-GARCH(1,1)")
dev.off()

vol_comparison <- na.omit(merge(vol_egarch11_dummy_roll, vol_garch11_dummy_roll, vol_gjr11t_dummy_roll))

write.xlsx(
  data.frame(Date = index(vol_comparison), coredata(vol_comparison)),
  file = file.path(tables_dir, "volatility_forecasts_comparison.xlsx"),
  rowNames = FALSE
)

png(file.path(figures_dir, "volatility_forecasts_comparison.png"),
    width = 1200, height = 800, res = 150)
plot(vol_comparison, main = "Volatility Forecast Comparison")
dev.off()


# 11A. QLIKE loss function comparison
# Ratio-based QLIKE is computed as h_t / hhat_t - log(h_t / hhat_t) - 1.
# Lower Mean QLIKE indicates better out-of-sample volatility forecasting.

# Realized variance proxy based on squared out-of-sample returns.
realized_var_oos <- xts(
  as.numeric(realizedret)^2,
  order.by = index(realizedret)
)

colnames(realized_var_oos) <- "Realized_Variance"

# Forecasted conditional variance is obtained by squaring the rolling volatility forecasts.
var_egarch11_dummy_roll <- vol_egarch11_dummy_roll^2
colnames(var_egarch11_dummy_roll) <- "Var_EGARCH11_sstd_DUMMY"

var_garch11_dummy_roll <- vol_garch11_dummy_roll^2
colnames(var_garch11_dummy_roll) <- "Var_GARCH11_sstd_DUMMY"

var_gjr11t_dummy_roll <- vol_gjr11t_dummy_roll^2
colnames(var_gjr11t_dummy_roll) <- "Var_GJR11_sstd_DUMMY"

# Add a small constant for numerical stability.
eps_qlike <- 1e-8

# Calculate the h_t / hhat_t ratio for each model.
ratio_egarch <- (realized_var_oos + eps_qlike) /
  (var_egarch11_dummy_roll + eps_qlike)

ratio_garch <- (realized_var_oos + eps_qlike) /
  (var_garch11_dummy_roll + eps_qlike)

ratio_gjr <- (realized_var_oos + eps_qlike) /
  (var_gjr11t_dummy_roll + eps_qlike)

# Compute the QLIKE series for each model.
qlike_egarch11_dummy <- ratio_egarch - log(ratio_egarch) - 1
colnames(qlike_egarch11_dummy) <- "QLIKE_EGARCH11_sstd_DUMMY"

qlike_garch11_dummy <- ratio_garch - log(ratio_garch) - 1
colnames(qlike_garch11_dummy) <- "QLIKE_GARCH11_sstd_DUMMY"

qlike_gjr11t_dummy <- ratio_gjr - log(ratio_gjr) - 1
colnames(qlike_gjr11t_dummy) <- "QLIKE_GJR11_sstd_DUMMY"

qlike_comparison_daily <- na.omit(merge(
  realized_var_oos,
  var_egarch11_dummy_roll,
  var_garch11_dummy_roll,
  var_gjr11t_dummy_roll,
  qlike_egarch11_dummy,
  qlike_garch11_dummy,
  qlike_gjr11t_dummy
))

# Build the QLIKE summary table.
qlike_summary <- data.frame(
  Model = c(
    "EGARCH(1,1)-sstd",
    "GARCH(1,1)-sstd",
    "GJR-GARCH(1,1)-sstd"
  ),
  
  Mean_QLIKE = c(
    mean(as.numeric(qlike_egarch11_dummy), na.rm = TRUE),
    mean(as.numeric(qlike_garch11_dummy), na.rm = TRUE),
    mean(as.numeric(qlike_gjr11t_dummy), na.rm = TRUE)
  ),
  
  Median_QLIKE = c(
    median(as.numeric(qlike_egarch11_dummy), na.rm = TRUE),
    median(as.numeric(qlike_garch11_dummy), na.rm = TRUE),
    median(as.numeric(qlike_gjr11t_dummy), na.rm = TRUE)
  ),
  
  SD_QLIKE = c(
    sd(as.numeric(qlike_egarch11_dummy), na.rm = TRUE),
    sd(as.numeric(qlike_garch11_dummy), na.rm = TRUE),
    sd(as.numeric(qlike_gjr11t_dummy), na.rm = TRUE)
  ),
  
  Min_QLIKE = c(
    min(as.numeric(qlike_egarch11_dummy), na.rm = TRUE),
    min(as.numeric(qlike_garch11_dummy), na.rm = TRUE),
    min(as.numeric(qlike_gjr11t_dummy), na.rm = TRUE)
  ),
  
  Max_QLIKE = c(
    max(as.numeric(qlike_egarch11_dummy), na.rm = TRUE),
    max(as.numeric(qlike_garch11_dummy), na.rm = TRUE),
    max(as.numeric(qlike_gjr11t_dummy), na.rm = TRUE)
  )
)

# Rank models by Mean QLIKE; lower values are preferred.
qlike_summary <- qlike_summary[order(qlike_summary$Mean_QLIKE), ]
qlike_summary$Ranking <- seq_len(nrow(qlike_summary))

qlike_summary <- qlike_summary[, c(
  "Ranking",
  "Model",
  "Mean_QLIKE",
  "Median_QLIKE",
  "SD_QLIKE",
  "Min_QLIKE",
  "Max_QLIKE"
)]

write.xlsx(
  list(
    "QLIKE_Summary" = qlike_summary,
    "QLIKE_Daily" = data.frame(
      Date = index(qlike_comparison_daily),
      coredata(qlike_comparison_daily)
    )
  ),
  file = file.path(tables_dir, "QLIKE_model_comparison_sstd_ratio.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "QLIKE_model_comparison_sstd_ratio.txt"))

cat("QLIKE MODEL COMPARISON - Skewed Student-t Models\n")
cat("Ratio-based QLIKE: h_t / hhat_t - log(h_t / hhat_t) - 1\n")
cat("Lower Mean QLIKE = better out-of-sample volatility forecasting\n\n")

print(qlike_summary)

sink()

cat("\n=====================================\n")
cat("RATIO-BASED QLIKE MODEL COMPARISON COMPLETED\n")
cat("Best model based on Mean QLIKE:\n")
cat(as.character(qlike_summary$Model[1]), "\n")
cat("Excel file saved at:\n")
cat(file.path(tables_dir, "QLIKE_model_comparison_sstd_ratio.xlsx"), "\n")
cat("=====================================\n")


# 12. VaR forecasts at the 1% level
var1_egarch11_dummy <- xts(
  roll_egarch11_dummy@forecast$VaR[, "alpha(1%)"],
  order.by = dates_oos
)
colnames(var1_egarch11_dummy) <- "VaR_1pct_EGARCH11_DUMMY"

var1_garch11_dummy <- xts(
  roll_garch11_dummy@forecast$VaR[, "alpha(1%)"],
  order.by = dates_oos
)
colnames(var1_garch11_dummy) <- "VaR_1pct_GARCH11_DUMMY"

var1_gjr11t_dummy <- xts(
  roll_gjr11t_dummy@forecast$VaR[, "alpha(1%)"],
  order.by = dates_oos
)
colnames(var1_gjr11t_dummy) <- "VaR_1pct_GJR11t_DUMMY"

write.xlsx(
  data.frame(Date = index(var1_egarch11_dummy), coredata(var1_egarch11_dummy)),
  file = file.path(tables_dir, "VaR_1pct_egarch11_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var1_garch11_dummy), coredata(var1_garch11_dummy)),
  file = file.path(tables_dir, "VaR_1pct_garch11_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var1_gjr11t_dummy), coredata(var1_gjr11t_dummy)),
  file = file.path(tables_dir, "VaR_1pct_gjr11t_dummy.xlsx"),
  rowNames = FALSE
)

var1_comparison <- na.omit(merge(var1_egarch11_dummy, var1_garch11_dummy, var1_gjr11t_dummy))

write.xlsx(
  data.frame(Date = index(var1_comparison), coredata(var1_comparison)),
  file = file.path(tables_dir, "VaR_1pct_forecasts_comparison.xlsx"),
  rowNames = FALSE
)

# Plot the 1% VaR forecasts against realized returns.
png(file.path(figures_dir, "VaRplot_1pct_egarch11_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.01, actual = realizedret, VaR = var1_egarch11_dummy)
title(main = "1% VaR Plot - EGARCH(1,1) ")
dev.off()

png(file.path(figures_dir, "VaRplot_1pct_garch11_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.01, actual = realizedret, VaR = var1_garch11_dummy)
title(main = "1% VaR Plot - GARCH(1,1) ")
dev.off()

png(file.path(figures_dir, "VaRplot_1pct_gjr11t_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.01, actual = realizedret, VaR = var1_gjr11t_dummy)
title(main = "1% VaR Plot - GJR-GARCH(1,1)")
dev.off()

# 13. VaR forecasts at the 5% level
var5_egarch11_dummy <- xts(
  roll_egarch11_dummy@forecast$VaR[, "alpha(5%)"],
  order.by = dates_oos
)
colnames(var5_egarch11_dummy) <- "VaR_5pct_EGARCH11_DUMMY"

var5_garch11_dummy <- xts(
  roll_garch11_dummy@forecast$VaR[, "alpha(5%)"],
  order.by = dates_oos
)
colnames(var5_garch11_dummy) <- "VaR_5pct_GARCH11_DUMMY"

var5_gjr11t_dummy <- xts(
  roll_gjr11t_dummy@forecast$VaR[, "alpha(5%)"],
  order.by = dates_oos
)
colnames(var5_gjr11t_dummy) <- "VaR_5pct_GJR11t_DUMMY"

write.xlsx(
  data.frame(Date = index(var5_egarch11_dummy), coredata(var5_egarch11_dummy)),
  file = file.path(tables_dir, "VaR_5pct_egarch11_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var5_garch11_dummy), coredata(var5_garch11_dummy)),
  file = file.path(tables_dir, "VaR_5pct_garch11_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var5_gjr11t_dummy), coredata(var5_gjr11t_dummy)),
  file = file.path(tables_dir, "VaR_5pct_gjr11t_dummy.xlsx"),
  rowNames = FALSE
)

var5_comparison <- na.omit(merge(var5_egarch11_dummy, var5_garch11_dummy, var5_gjr11t_dummy))

write.xlsx(
  data.frame(Date = index(var5_comparison), coredata(var5_comparison)),
  file = file.path(tables_dir, "VaR_5pct_forecasts_comparison.xlsx"),
  rowNames = FALSE
)

# Plot the 5% VaR forecasts against realized returns.
png(file.path(figures_dir, "VaRplot_5pct_egarch11_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.05, actual = realizedret, VaR = var5_egarch11_dummy)
title(main = "5% VaR Plot - EGARCH(1,1) ")
dev.off()

png(file.path(figures_dir, "VaRplot_5pct_garch11_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.05, actual = realizedret, VaR = var5_garch11_dummy)
title(main = "5% VaR Plot - GARCH(1,1) ")
dev.off()

png(file.path(figures_dir, "VaRplot_5pct_gjr11t_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.05, actual = realizedret, VaR = var5_gjr11t_dummy)
title(main = "5% VaR Plot - GJR-GARCH(1,1)")
dev.off()

# 14. VaR backtesting at the 1% level
test1_egarch11_dummy <- VaRTest(
  alpha = 0.01,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var1_egarch11_dummy)
)

test1_garch11_dummy <- VaRTest(
  alpha = 0.01,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var1_garch11_dummy)
)

test1_gjr11t_dummy <- VaRTest(
  alpha = 0.01,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var1_gjr11t_dummy)
)

sink(file.path(output_dir, "VaR_tests_1pct_all_models.txt"))
cat("1% VaR Backtest Results\n\n")

cat("===== EGARCH(1,1) =====\n")
print(test1_egarch11_dummy)

cat("\n===== GARCH(1,1) =====\n")
print(test1_garch11_dummy)

cat("\n===== GJR-GARCH(1,1) =====\n")
print(test1_gjr11t_dummy)
sink()

var1_backtest_comparison <- data.frame(
  Model = c("EGARCH(1,1)",
            "GARCH(1,1)",
            "GJR-GARCH(1,1)"),
  Expected_Exceed = c(test1_egarch11_dummy$expected.exceed,
                      test1_garch11_dummy$expected.exceed,
                      test1_gjr11t_dummy$expected.exceed),
  Actual_Exceed = c(test1_egarch11_dummy$actual.exceed,
                    test1_garch11_dummy$actual.exceed,
                    test1_gjr11t_dummy$actual.exceed),
  UC_LRstat = c(test1_egarch11_dummy$uc.LRstat,
                test1_garch11_dummy$uc.LRstat,
                test1_gjr11t_dummy$uc.LRstat),
  UC_pvalue = c(test1_egarch11_dummy$uc.LRp,
                test1_garch11_dummy$uc.LRp,
                test1_gjr11t_dummy$uc.LRp),
  UC_Decision = c(test1_egarch11_dummy$uc.Decision,
                  test1_garch11_dummy$uc.Decision,
                  test1_gjr11t_dummy$uc.Decision),
  CC_LRstat = c(test1_egarch11_dummy$cc.LRstat,
                test1_garch11_dummy$cc.LRstat,
                test1_gjr11t_dummy$cc.LRstat),
  CC_pvalue = c(test1_egarch11_dummy$cc.LRp,
                test1_garch11_dummy$cc.LRp,
                test1_gjr11t_dummy$cc.LRp),
  CC_Decision = c(test1_egarch11_dummy$cc.Decision,
                  test1_garch11_dummy$cc.Decision,
                  test1_gjr11t_dummy$cc.Decision)
)

write.xlsx(
  var1_backtest_comparison,
  file = file.path(tables_dir, "VaR_1pct_backtest_comparison.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "VaR_1pct_backtest_comparison.txt"))
cat("COMPARATIVE 1% VaR BACKTEST RESULTS\n\n")
print(var1_backtest_comparison)
sink()

# 15. VaR backtesting at the 5% level
test5_egarch11_dummy <- VaRTest(
  alpha = 0.05,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var5_egarch11_dummy)
)

test5_garch11_dummy <- VaRTest(
  alpha = 0.05,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var5_garch11_dummy)
)

test5_gjr11t_dummy <- VaRTest(
  alpha = 0.05,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var5_gjr11t_dummy)
)

sink(file.path(output_dir, "VaR_tests_5pct_all_models.txt"))
cat("5% VaR Backtest Results\n\n")

cat("===== EGARCH(1,1) =====\n")
print(test5_egarch11_dummy)

cat("\n===== GARCH(1,1) =====\n")
print(test5_garch11_dummy)

cat("\n===== GJR-GARCH(1,1) =====\n")
print(test5_gjr11t_dummy)
sink()

var5_backtest_comparison <- data.frame(
  Model = c("EGARCH(1,1)",
            "GARCH(1,1)",
            "GJR-GARCH(1,1)"),
  Expected_Exceed = c(test5_egarch11_dummy$expected.exceed,
                      test5_garch11_dummy$expected.exceed,
                      test5_gjr11t_dummy$expected.exceed),
  Actual_Exceed = c(test5_egarch11_dummy$actual.exceed,
                    test5_garch11_dummy$actual.exceed,
                    test5_gjr11t_dummy$actual.exceed),
  UC_LRstat = c(test5_egarch11_dummy$uc.LRstat,
                test5_garch11_dummy$uc.LRstat,
                test5_gjr11t_dummy$uc.LRstat),
  UC_pvalue = c(test5_egarch11_dummy$uc.LRp,
                test5_garch11_dummy$uc.LRp,
                test5_gjr11t_dummy$uc.LRp),
  UC_Decision = c(test5_egarch11_dummy$uc.Decision,
                  test5_garch11_dummy$uc.Decision,
                  test5_gjr11t_dummy$uc.Decision),
  CC_LRstat = c(test5_egarch11_dummy$cc.LRstat,
                test5_garch11_dummy$cc.LRstat,
                test5_gjr11t_dummy$cc.LRstat),
  CC_pvalue = c(test5_egarch11_dummy$cc.LRp,
                test5_garch11_dummy$cc.LRp,
                test5_gjr11t_dummy$cc.LRp),
  CC_Decision = c(test5_egarch11_dummy$cc.Decision,
                  test5_garch11_dummy$cc.Decision,
                  test5_gjr11t_dummy$cc.Decision)
)

write.xlsx(
  var5_backtest_comparison,
  file = file.path(tables_dir, "VaR_5pct_backtest_comparison.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "VaR_5pct_backtest_comparison.txt"))
cat("COMPARATIVE 5% VaR BACKTEST RESULTS\n\n")
print(var5_backtest_comparison)
sink()

# 16. Final summary of VaR backtests
backtest_summary_all <- data.frame(
  Model = rep(c("EGARCH(1,1)",
                "GARCH(1,1)",
                "GJR-GARCH(1,1)"), each = 2),
  VaR_Level = rep(c("1%", "5%"), times = 3),
  Expected_Exceed = c(
    test1_egarch11_dummy$expected.exceed, test5_egarch11_dummy$expected.exceed,
    test1_garch11_dummy$expected.exceed,  test5_garch11_dummy$expected.exceed,
    test1_gjr11t_dummy$expected.exceed,   test5_gjr11t_dummy$expected.exceed
  ),
  Actual_Exceed = c(
    test1_egarch11_dummy$actual.exceed, test5_egarch11_dummy$actual.exceed,
    test1_garch11_dummy$actual.exceed,  test5_garch11_dummy$actual.exceed,
    test1_gjr11t_dummy$actual.exceed,   test5_gjr11t_dummy$actual.exceed
  ),
  UC_LRstat = c(
    test1_egarch11_dummy$uc.LRstat, test5_egarch11_dummy$uc.LRstat,
    test1_garch11_dummy$uc.LRstat,  test5_garch11_dummy$uc.LRstat,
    test1_gjr11t_dummy$uc.LRstat,   test5_gjr11t_dummy$uc.LRstat
  ),
  UC_pvalue = c(
    test1_egarch11_dummy$uc.LRp, test5_egarch11_dummy$uc.LRp,
    test1_garch11_dummy$uc.LRp,  test5_garch11_dummy$uc.LRp,
    test1_gjr11t_dummy$uc.LRp,   test5_gjr11t_dummy$uc.LRp
  ),
  UC_Decision = c(
    test1_egarch11_dummy$uc.Decision, test5_egarch11_dummy$uc.Decision,
    test1_garch11_dummy$uc.Decision,  test5_garch11_dummy$uc.Decision,
    test1_gjr11t_dummy$uc.Decision,   test5_gjr11t_dummy$uc.Decision
  ),
  CC_LRstat = c(
    test1_egarch11_dummy$cc.LRstat, test5_egarch11_dummy$cc.LRstat,
    test1_garch11_dummy$cc.LRstat,  test5_garch11_dummy$cc.LRstat,
    test1_gjr11t_dummy$cc.LRstat,   test5_gjr11t_dummy$cc.LRstat
  ),
  CC_pvalue = c(
    test1_egarch11_dummy$cc.LRp, test5_egarch11_dummy$cc.LRp,
    test1_garch11_dummy$cc.LRp,  test5_garch11_dummy$cc.LRp,
    test1_gjr11t_dummy$cc.LRp,   test5_gjr11t_dummy$cc.LRp
  ),
  CC_Decision = c(
    test1_egarch11_dummy$cc.Decision, test5_egarch11_dummy$cc.Decision,
    test1_garch11_dummy$cc.Decision,  test5_garch11_dummy$cc.Decision,
    test1_gjr11t_dummy$cc.Decision,   test5_gjr11t_dummy$cc.Decision
  )
)

write.xlsx(
  backtest_summary_all,
  file = file.path(tables_dir, "VaR_backtests_summary_all_models.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "VaR_backtests_summary_all_models.txt"))
cat("FINAL SUMMARY OF VaR BACKTESTS FOR ALL MODELS\n\n")
print(backtest_summary_all)
sink()

# 17. Save workspace
save.image(file = file.path(output_dir, "workspace_SP500_45_with_VIX_dummy_1pct_5pct.RData"))

cat("Script execution completed successfully.\n")
cat("VaR was calculated and backtested at the 1% and 5% levels.\n")
cat("All results were saved to:\n")
cat(output_dir, "\n")
