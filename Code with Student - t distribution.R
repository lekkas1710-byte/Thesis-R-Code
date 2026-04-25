# SP500-45 volatility modelling with Student-t GARCH-family models
# The lagged VIX-up dummy is used as an external regressor in each variance equation.
# VIX_up_dummy_t equals 1 when VIX_{t-1} > VIX_{t-2}; otherwise, it equals 0.
# The script evaluates EGARCH, GARCH, and GJR-GARCH models using 1% and 5% VaR, QLIKE, and final model rankings.
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
vix_reg_mat <- matrix(as.numeric(vix_dummy_lag1), ncol = 1)
colnames(vix_reg_mat) <- "VIX_up_dummy_lag1"

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

# Out-of-sample evaluation starts on 1 January 2022.
realizedret <- returns_final["2022-01-01/"]
colnames(realizedret) <- "Realized_Return"

T2 <- NROW(realizedret)
dates_oos <- index(realizedret)

sink(file.path(output_dir, "sample_info.txt"))
cat("In-sample observations (T1): ", T1, "\n", sep = "")
cat("Out-of-sample observations (T2): ", T2, "\n\n", sep = "")
cat("Head of in-sample returns:\n")
print(head(returns_is))
cat("\nTail of in-sample returns:\n")
print(tail(returns_is))
cat("\nHead of out-of-sample returns:\n")
print(head(realizedret))
cat("\nTail of out-of-sample returns:\n")
print(tail(realizedret))
sink()

write.xlsx(
  data.frame(Date = index(realizedret), coredata(realizedret)),
  file = file.path(tables_dir, "realized_returns_oos.xlsx"),
  rowNames = FALSE
)

# 6. Model specifications
# The same lagged VIX-up dummy is included in each variance equation.

# Model 1: EGARCH(1,1) with Student-t innovations.
spec_egarch11_std_dummy <- ugarchspec(
  variance.model = list(
    model = "eGARCH",
    garchOrder = c(1, 1),
    external.regressors = vix_reg_mat
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

# Model 2: GARCH(1,1) with Student-t innovations.
spec_garch11_std_dummy <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1),
    external.regressors = vix_reg_mat
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

# Model 3: GJR-GARCH(1,1) with Student-t innovations.
spec_gjr11_std_dummy <- ugarchspec(
  variance.model = list(
    model = "gjrGARCH",
    garchOrder = c(1, 1),
    external.regressors = vix_reg_mat
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

# 7. Model estimation
egarch11_std_dummy <- ugarchfit(
  spec = spec_egarch11_std_dummy,
  data = returns_final,
  out.sample = T2,
  solver = "hybrid"
)

garch11_std_dummy <- ugarchfit(
  spec = spec_garch11_std_dummy,
  data = returns_final,
  out.sample = T2,
  solver = "hybrid"
)

gjr11_std_dummy <- ugarchfit(
  spec = spec_gjr11_std_dummy,
  data = returns_final,
  out.sample = T2,
  solver = "hybrid"
)

# Check whether each model converged successfully.
if (convergence(egarch11_std_dummy) != 0) warning("EGARCH(1,1)-std fit did not converge.")
if (convergence(garch11_std_dummy)  != 0) warning("GARCH(1,1)-std fit did not converge.")
if (convergence(gjr11_std_dummy)    != 0) warning("GJR-GARCH(1,1)-std fit did not converge.")

# Save the complete estimation outputs.
sink(file.path(output_dir, "egarch11_std_dummy_output.txt"))
cat("EGARCH(1,1)-std with lagged VIX-up dummy in variance equation\n\n")
show(egarch11_std_dummy)
sink()

sink(file.path(output_dir, "garch11_std_dummy_output.txt"))
cat("GARCH(1,1)-std with lagged VIX-up dummy in variance equation\n\n")
show(garch11_std_dummy)
sink()

sink(file.path(output_dir, "gjr11_std_dummy_output.txt"))
cat("GJR-GARCH(1,1)-std with lagged VIX-up dummy in variance equation\n\n")
show(gjr11_std_dummy)
sink()

# 7A. Estimated coefficient tables
# Classical and robust coefficient statistics are exported for each model.

extract_garch_table <- function(fit_obj, model_name, robust = FALSE) {

  # Select the classical or robust coefficient matrix.
  coef_mat <- if (robust) fit_obj@fit$robust.matcoef else fit_obj@fit$matcoef

  # Return a clear note when the requested matrix is unavailable.
  if (is.null(coef_mat)) {
    return(data.frame(
      Model = model_name,
      Note = if (robust) {
        "Robust coefficient matrix not available"
      } else {
        "Coefficient matrix not available"
      }
    ))
  }

  out <- as.data.frame(coef_mat, check.names = FALSE)
  out$Parameter <- rownames(out)
  rownames(out) <- NULL

  out <- out[, c("Parameter", setdiff(names(out), "Parameter"))]

  out$Model <- model_name
  out <- out[, c("Model", setdiff(names(out), "Model"))]

  # Add significance stars when a p-value column is available.
  pcol <- grep("Pr\\(", names(out), value = TRUE)
  if (length(pcol) > 0) {
    pcol <- pcol[1]
    out$Significance <- ifelse(out[[pcol]] < 0.01, "***",
                               ifelse(out[[pcol]] < 0.05, "**",
                                      ifelse(out[[pcol]] < 0.10, "*", "")))
  }

  return(out)
}

# Classical coefficient tables.
coef_egarch_classic <- extract_garch_table(
  egarch11_std_dummy,
  "EGARCH(1,1)",
  robust = FALSE
)

coef_garch_classic <- extract_garch_table(
  garch11_std_dummy,
  "GARCH(1,1)",
  robust = FALSE
)

coef_gjr_classic <- extract_garch_table(
  gjr11_std_dummy,
  "GJR-GARCH(1,1)",
  robust = FALSE
)

# Robust coefficient tables.
coef_egarch_robust <- extract_garch_table(
  egarch11_std_dummy,
  "EGARCH(1,1)",
  robust = TRUE
)

coef_garch_robust <- extract_garch_table(
  garch11_std_dummy,
  "GARCH(1,1)",
  robust = TRUE
)

coef_gjr_robust <- extract_garch_table(
  gjr11_std_dummy,
  "GJR-GARCH(1,1)",
  robust = TRUE
)

# Summary table for the external variance regressor vxreg1.
get_vxreg1_row <- function(df, type_label) {
  if ("Parameter" %in% names(df) && any(df$Parameter == "vxreg1")) {
    tmp <- df[df$Parameter == "vxreg1", , drop = FALSE]
    tmp$Table_Type <- type_label
    return(tmp)
  } else {
    return(data.frame(
      Model = unique(df$Model)[1],
      Parameter = "vxreg1",
      Table_Type = type_label,
      Note = "vxreg1 not found"
    ))
  }
}

dummy_summary <- rbind(
  get_vxreg1_row(coef_egarch_classic, "Classic"),
  get_vxreg1_row(coef_garch_classic,  "Classic"),
  get_vxreg1_row(coef_gjr_classic,    "Classic"),
  get_vxreg1_row(coef_egarch_robust,  "Robust"),
  get_vxreg1_row(coef_garch_robust,   "Robust"),
  get_vxreg1_row(coef_gjr_robust,     "Robust")
)

# Export all coefficient tables to one Excel workbook.
write.xlsx(
  list(
    "EGARCH_Classic" = coef_egarch_classic,
    "GARCH_Classic"  = coef_garch_classic,
    "GJR_Classic"    = coef_gjr_classic,
    "EGARCH_Robust"  = coef_egarch_robust,
    "GARCH_Robust"   = coef_garch_robust,
    "GJR_Robust"     = coef_gjr_robust,
    "Dummy_vxreg1"   = dummy_summary
  ),
  file = file.path(tables_dir, "estimated_parameters_detailed.xlsx"),
  rowNames = FALSE
)

# Save the coefficient tables as a text report.
sink(file.path(output_dir, "estimated_parameters_detailed.txt"))

cat("===== EGARCH CLASSIC =====\n")
print(coef_egarch_classic)

cat("\n\n===== GARCH CLASSIC =====\n")
print(coef_garch_classic)

cat("\n\n===== GJR CLASSIC =====\n")
print(coef_gjr_classic)

cat("\n\n===== EGARCH ROBUST =====\n")
print(coef_egarch_robust)

cat("\n\n===== GARCH ROBUST =====\n")
print(coef_garch_robust)

cat("\n\n===== GJR ROBUST =====\n")
print(coef_gjr_robust)

sink()

# 8. Residual diagnostics
# In rugarch, which = 11 plots the ACF of squared standardized residuals.
png(file.path(figures_dir, "squared_residual_egarch11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
plot(egarch11_std_dummy, which = 11)
dev.off()

png(file.path(figures_dir, "squared_residual_garch11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
plot(garch11_std_dummy, which = 11)
dev.off()

png(file.path(figures_dir, "squared_residual_gjr11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
plot(gjr11_std_dummy, which = 11)
dev.off()

# 9. In-sample model comparison
ic_egarch <- infocriteria(egarch11_std_dummy)
ic_garch  <- infocriteria(garch11_std_dummy)
ic_gjr    <- infocriteria(gjr11_std_dummy)

coef_egarch <- coef(egarch11_std_dummy)
coef_garch  <- coef(garch11_std_dummy)
coef_gjr    <- coef(gjr11_std_dummy)

model_comparison <- data.frame(
  Model = c("EGARCH(1,1)",
            "GARCH(1,1)",
            "GJR-GARCH(1,1)"),
  LogLikelihood = c(likelihood(egarch11_std_dummy),
                    likelihood(garch11_std_dummy),
                    likelihood(gjr11_std_dummy)),
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
  file = file.path(tables_dir, "model_comparison_insample_std.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "model_comparison_insample_std.txt"))
cat("IN-SAMPLE MODEL COMPARISON - Student-t\n\n")
print(model_comparison)
sink()

# 10. Rolling one-step-ahead forecasts
# VaR is calculated at the 1% and 5% levels.
roll_egarch11_std_dummy <- ugarchroll(
  spec_egarch11_std_dummy,
  data = returns_final,
  n.start = T1,
  refit.every = 10,
  refit.window = "moving",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

roll_garch11_std_dummy <- ugarchroll(
  spec_garch11_std_dummy,
  data = returns_final,
  n.start = T1,
  refit.every = 10,
  refit.window = "moving",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

roll_gjr11_std_dummy <- ugarchroll(
  spec_gjr11_std_dummy,
  data = returns_final,
  n.start = T1,
  refit.every = 10,
  refit.window = "moving",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

# Check rolling estimation convergence.
if (convergence(roll_egarch11_std_dummy) != 0) warning("EGARCH(1,1)-std roll has non-converged windows.")
if (convergence(roll_garch11_std_dummy)  != 0) warning("GARCH(1,1)-std roll has non-converged windows.")
if (convergence(roll_gjr11_std_dummy)    != 0) warning("GJR-GARCH(1,1)-std roll has non-converged windows.")

# 11. Rolling volatility forecasts
vol_egarch11_std_dummy_roll <- xts(
  roll_egarch11_std_dummy@forecast$density$Sigma,
  order.by = dates_oos
)
colnames(vol_egarch11_std_dummy_roll) <- "Vol_EGARCH11_STD_DUMMY"

vol_garch11_std_dummy_roll <- xts(
  roll_garch11_std_dummy@forecast$density$Sigma,
  order.by = dates_oos
)
colnames(vol_garch11_std_dummy_roll) <- "Vol_GARCH11_STD_DUMMY"

vol_gjr11_std_dummy_roll <- xts(
  roll_gjr11_std_dummy@forecast$density$Sigma,
  order.by = dates_oos
)
colnames(vol_gjr11_std_dummy_roll) <- "Vol_GJR11_STD_DUMMY"

write.xlsx(
  data.frame(Date = index(vol_egarch11_std_dummy_roll), coredata(vol_egarch11_std_dummy_roll)),
  file = file.path(tables_dir, "vol_egarch11_std_dummy_roll.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(vol_garch11_std_dummy_roll), coredata(vol_garch11_std_dummy_roll)),
  file = file.path(tables_dir, "vol_garch11_std_dummy_roll.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(vol_gjr11_std_dummy_roll), coredata(vol_gjr11_std_dummy_roll)),
  file = file.path(tables_dir, "vol_gjr11_std_dummy_roll.xlsx"),
  rowNames = FALSE
)

png(file.path(figures_dir, "vol_egarch11_std_dummy_roll.png"),
    width = 1200, height = 800, res = 150)
plot(vol_egarch11_std_dummy_roll, main = "Rolling Volatility Forecast - EGARCH(1,1)")
dev.off()

png(file.path(figures_dir, "vol_garch11_std_dummy_roll.png"),
    width = 1200, height = 800, res = 150)
plot(vol_garch11_std_dummy_roll, main = "Rolling Volatility Forecast - GARCH(1,1)")
dev.off()

png(file.path(figures_dir, "vol_gjr11_std_dummy_roll.png"),
    width = 1200, height = 800, res = 150)
plot(vol_gjr11_std_dummy_roll, main = "Rolling Volatility Forecast - GJR-GARCH(1,1)")
dev.off()

vol_comparison <- na.omit(merge(
  vol_egarch11_std_dummy_roll,
  vol_garch11_std_dummy_roll,
  vol_gjr11_std_dummy_roll
))

write.xlsx(
  data.frame(Date = index(vol_comparison), coredata(vol_comparison)),
  file = file.path(tables_dir, "volatility_forecasts_comparison_std.xlsx"),
  rowNames = FALSE
)

png(file.path(figures_dir, "volatility_forecasts_comparison_std.png"),
    width = 1200, height = 800, res = 150)
plot(vol_comparison, main = "Volatility Forecast Comparison - Student-t Models")
dev.off()

# 11A. QLIKE volatility forecast evaluation
# Ratio-based QLIKE:
# QLIKE_t = (h_t / hhat_t) - log(h_t / hhat_t) - 1
# Lower Mean QLIKE indicates better volatility forecasting performance.
# where:
# hhat_t is the forecasted conditional variance, sigma_t^2.
# h_t is the realized variance proxy, r_t^2.

# Realized variance proxy based on out-of-sample returns.
realized_var_oos <- xts(
  as.numeric(realizedret)^2,
  order.by = index(realizedret)
)

colnames(realized_var_oos) <- "Realized_Variance"

# Convert forecasted volatility, Sigma, to forecasted variance.
var_egarch11_std_dummy_roll <- vol_egarch11_std_dummy_roll^2
colnames(var_egarch11_std_dummy_roll) <- "Var_EGARCH11_STD_DUMMY"

var_garch11_std_dummy_roll <- vol_garch11_std_dummy_roll^2
colnames(var_garch11_std_dummy_roll) <- "Var_GARCH11_STD_DUMMY"

var_gjr11_std_dummy_roll <- vol_gjr11_std_dummy_roll^2
colnames(var_gjr11_std_dummy_roll) <- "Var_GJR11_STD_DUMMY"

# Small constant used for numerical stability.
# This prevents division by zero and log(0).
eps_qlike <- 1e-8

# Ratio h_t / hhat_t.
ratio_egarch11_std_dummy <- (realized_var_oos + eps_qlike) /
  (var_egarch11_std_dummy_roll + eps_qlike)

colnames(ratio_egarch11_std_dummy) <- "Ratio_EGARCH11_STD_DUMMY"


ratio_garch11_std_dummy <- (realized_var_oos + eps_qlike) /
  (var_garch11_std_dummy_roll + eps_qlike)

colnames(ratio_garch11_std_dummy) <- "Ratio_GARCH11_STD_DUMMY"


ratio_gjr11_std_dummy <- (realized_var_oos + eps_qlike) /
  (var_gjr11_std_dummy_roll + eps_qlike)

colnames(ratio_gjr11_std_dummy) <- "Ratio_GJR11_STD_DUMMY"


# QLIKE series for each model.
qlike_egarch11_std_dummy <- ratio_egarch11_std_dummy -
  log(ratio_egarch11_std_dummy) - 1

colnames(qlike_egarch11_std_dummy) <- "QLIKE_EGARCH11_STD_DUMMY"


qlike_garch11_std_dummy <- ratio_garch11_std_dummy -
  log(ratio_garch11_std_dummy) - 1

colnames(qlike_garch11_std_dummy) <- "QLIKE_GARCH11_STD_DUMMY"


qlike_gjr11_std_dummy <- ratio_gjr11_std_dummy -
  log(ratio_gjr11_std_dummy) - 1

colnames(qlike_gjr11_std_dummy) <- "QLIKE_GJR11_STD_DUMMY"


# Daily QLIKE comparison table.
qlike_comparison_daily <- na.omit(merge(
  realized_var_oos,
  var_egarch11_std_dummy_roll,
  var_garch11_std_dummy_roll,
  var_gjr11_std_dummy_roll,
  ratio_egarch11_std_dummy,
  ratio_garch11_std_dummy,
  ratio_gjr11_std_dummy,
  qlike_egarch11_std_dummy,
  qlike_garch11_std_dummy,
  qlike_gjr11_std_dummy
))

# Summary QLIKE comparison table.
qlike_summary <- data.frame(
  Model = c(
    "EGARCH(1,1)",
    "GARCH(1,1)",
    "GJR-GARCH(1,1)"
  ),

  Mean_QLIKE = c(
    mean(as.numeric(qlike_egarch11_std_dummy), na.rm = TRUE),
    mean(as.numeric(qlike_garch11_std_dummy), na.rm = TRUE),
    mean(as.numeric(qlike_gjr11_std_dummy), na.rm = TRUE)
  ),

  Median_QLIKE = c(
    median(as.numeric(qlike_egarch11_std_dummy), na.rm = TRUE),
    median(as.numeric(qlike_garch11_std_dummy), na.rm = TRUE),
    median(as.numeric(qlike_gjr11_std_dummy), na.rm = TRUE)
  ),

  SD_QLIKE = c(
    sd(as.numeric(qlike_egarch11_std_dummy), na.rm = TRUE),
    sd(as.numeric(qlike_garch11_std_dummy), na.rm = TRUE),
    sd(as.numeric(qlike_gjr11_std_dummy), na.rm = TRUE)
  ),

  Min_QLIKE = c(
    min(as.numeric(qlike_egarch11_std_dummy), na.rm = TRUE),
    min(as.numeric(qlike_garch11_std_dummy), na.rm = TRUE),
    min(as.numeric(qlike_gjr11_std_dummy), na.rm = TRUE)
  ),

  Max_QLIKE = c(
    max(as.numeric(qlike_egarch11_std_dummy), na.rm = TRUE),
    max(as.numeric(qlike_garch11_std_dummy), na.rm = TRUE),
    max(as.numeric(qlike_gjr11_std_dummy), na.rm = TRUE)
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

sink(file.path(output_dir, "QLIKE_model_comparison_std.txt"))

cat("QLIKE MODEL COMPARISON - Student-t Models\n")
cat("Ratio-based QLIKE: h_t / hhat_t - log(h_t / hhat_t) - 1\n")
cat("Lower Mean QLIKE = better out-of-sample volatility forecasting\n\n")

print(qlike_summary)

sink()

write.xlsx(
  list(
    "QLIKE_Summary" = qlike_summary,
    "QLIKE_Daily" = data.frame(
      Date = index(qlike_comparison_daily),
      coredata(qlike_comparison_daily)
    )
  ),
  file = file.path(tables_dir, "QLIKE_model_comparison_std.xlsx"),
  rowNames = FALSE
)

# Plot the QLIKE series.
png(
  file.path(figures_dir, "QLIKE_comparison_std.png"),
  width = 1200,
  height = 800,
  res = 150
)

plot(
  qlike_comparison_daily[, c(
    "QLIKE_EGARCH11_STD_DUMMY",
    "QLIKE_GARCH11_STD_DUMMY",
    "QLIKE_GJR11_STD_DUMMY"
  )],
  main = "QLIKE Comparison - Student-t Models"
)

dev.off()

cat("\n=====================================\n")
cat("QLIKE MODEL COMPARISON COMPLETED - Student-t\n")
cat("Best model based on Mean QLIKE:\n")
cat(as.character(qlike_summary$Model[1]), "\n")
cat("Excel file saved at:\n")
cat(file.path(tables_dir, "QLIKE_model_comparison_std.xlsx"), "\n")
cat("=====================================\n")

# 12. VaR forecasts at the 1% level
var1_egarch11_std_dummy <- xts(
  roll_egarch11_std_dummy@forecast$VaR[, "alpha(1%)"],
  order.by = dates_oos
)
colnames(var1_egarch11_std_dummy) <- "VaR_1pct_EGARCH11_STD_DUMMY"

var1_garch11_std_dummy <- xts(
  roll_garch11_std_dummy@forecast$VaR[, "alpha(1%)"],
  order.by = dates_oos
)
colnames(var1_garch11_std_dummy) <- "VaR_1pct_GARCH11_STD_DUMMY"

var1_gjr11_std_dummy <- xts(
  roll_gjr11_std_dummy@forecast$VaR[, "alpha(1%)"],
  order.by = dates_oos
)
colnames(var1_gjr11_std_dummy) <- "VaR_1pct_GJR11_STD_DUMMY"

write.xlsx(
  data.frame(Date = index(var1_egarch11_std_dummy), coredata(var1_egarch11_std_dummy)),
  file = file.path(tables_dir, "VaR_1pct_egarch11_std_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var1_garch11_std_dummy), coredata(var1_garch11_std_dummy)),
  file = file.path(tables_dir, "VaR_1pct_garch11_std_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var1_gjr11_std_dummy), coredata(var1_gjr11_std_dummy)),
  file = file.path(tables_dir, "VaR_1pct_gjr11_std_dummy.xlsx"),
  rowNames = FALSE
)

var1_comparison <- na.omit(merge(
  var1_egarch11_std_dummy,
  var1_garch11_std_dummy,
  var1_gjr11_std_dummy
))

write.xlsx(
  data.frame(Date = index(var1_comparison), coredata(var1_comparison)),
  file = file.path(tables_dir, "VaR_1pct_forecasts_comparison_std.xlsx"),
  rowNames = FALSE
)

# Plot 1% VaR exceedances.
png(file.path(figures_dir, "VaRplot_1pct_egarch11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.01, actual = realizedret, VaR = var1_egarch11_std_dummy)
title(main = "1% VaR Plot - EGARCH(1,1)")
dev.off()

png(file.path(figures_dir, "VaRplot_1pct_garch11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.01, actual = realizedret, VaR = var1_garch11_std_dummy)
title(main = "1% VaR Plot - GARCH(1,1)")
dev.off()

png(file.path(figures_dir, "VaRplot_1pct_gjr11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.01, actual = realizedret, VaR = var1_gjr11_std_dummy)
title(main = "1% VaR Plot - GJR-GARCH(1,1)")
dev.off()

# 13. VaR forecasts at the 5% level
var5_egarch11_std_dummy <- xts(
  roll_egarch11_std_dummy@forecast$VaR[, "alpha(5%)"],
  order.by = dates_oos
)
colnames(var5_egarch11_std_dummy) <- "VaR_5pct_EGARCH11_STD_DUMMY"

var5_garch11_std_dummy <- xts(
  roll_garch11_std_dummy@forecast$VaR[, "alpha(5%)"],
  order.by = dates_oos
)
colnames(var5_garch11_std_dummy) <- "VaR_5pct_GARCH11_STD_DUMMY"

var5_gjr11_std_dummy <- xts(
  roll_gjr11_std_dummy@forecast$VaR[, "alpha(5%)"],
  order.by = dates_oos
)
colnames(var5_gjr11_std_dummy) <- "VaR_5pct_GJR11_STD_DUMMY"

write.xlsx(
  data.frame(Date = index(var5_egarch11_std_dummy), coredata(var5_egarch11_std_dummy)),
  file = file.path(tables_dir, "VaR_5pct_egarch11_std_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var5_garch11_std_dummy), coredata(var5_garch11_std_dummy)),
  file = file.path(tables_dir, "VaR_5pct_garch11_std_dummy.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  data.frame(Date = index(var5_gjr11_std_dummy), coredata(var5_gjr11_std_dummy)),
  file = file.path(tables_dir, "VaR_5pct_gjr11_std_dummy.xlsx"),
  rowNames = FALSE
)

var5_comparison <- na.omit(merge(
  var5_egarch11_std_dummy,
  var5_garch11_std_dummy,
  var5_gjr11_std_dummy
))

write.xlsx(
  data.frame(Date = index(var5_comparison), coredata(var5_comparison)),
  file = file.path(tables_dir, "VaR_5pct_forecasts_comparison_std.xlsx"),
  rowNames = FALSE
)

# Plot 5% VaR exceedances.
png(file.path(figures_dir, "VaRplot_5pct_egarch11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.05, actual = realizedret, VaR = var5_egarch11_std_dummy)
title(main = "5% VaR Plot - EGARCH(1,1)")
dev.off()

png(file.path(figures_dir, "VaRplot_5pct_garch11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.05, actual = realizedret, VaR = var5_garch11_std_dummy)
title(main = "5% VaR Plot - GARCH(1,1)")
dev.off()

png(file.path(figures_dir, "VaRplot_5pct_gjr11_std_dummy.png"),
    width = 1200, height = 800, res = 150)
VaRplot(alpha = 0.05, actual = realizedret, VaR = var5_gjr11_std_dummy)
title(main = "5% VaR Plot - GJR-GARCH(1,1)")
dev.off()

# 14. VaR backtesting at the 1% level
test1_egarch11_std_dummy <- VaRTest(
  alpha = 0.01,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var1_egarch11_std_dummy)
)

test1_garch11_std_dummy <- VaRTest(
  alpha = 0.01,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var1_garch11_std_dummy)
)

test1_gjr11_std_dummy <- VaRTest(
  alpha = 0.01,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var1_gjr11_std_dummy)
)

sink(file.path(output_dir, "VaR_tests_1pct_all_models_std.txt"))
cat("1% VaR Backtest Results - Student-t Models\n\n")

cat("===== EGARCH(1,1) =====\n")
print(test1_egarch11_std_dummy)

cat("\n===== GARCH(1,1) =====\n")
print(test1_garch11_std_dummy)

cat("\n===== GJR-GARCH(1,1) =====\n")
print(test1_gjr11_std_dummy)
sink()

var1_backtest_comparison <- data.frame(
  Model = c("EGARCH(1,1)",
            "GARCH(1,1)",
            "GJR-GARCH(1,1)"),
  Expected_Exceed = c(test1_egarch11_std_dummy$expected.exceed,
                      test1_garch11_std_dummy$expected.exceed,
                      test1_gjr11_std_dummy$expected.exceed),
  Actual_Exceed = c(test1_egarch11_std_dummy$actual.exceed,
                    test1_garch11_std_dummy$actual.exceed,
                    test1_gjr11_std_dummy$actual.exceed),
  UC_LRstat = c(test1_egarch11_std_dummy$uc.LRstat,
                test1_garch11_std_dummy$uc.LRstat,
                test1_gjr11_std_dummy$uc.LRstat),
  UC_pvalue = c(test1_egarch11_std_dummy$uc.LRp,
                test1_garch11_std_dummy$uc.LRp,
                test1_gjr11_std_dummy$uc.LRp),
  UC_Decision = c(test1_egarch11_std_dummy$uc.Decision,
                  test1_garch11_std_dummy$uc.Decision,
                  test1_gjr11_std_dummy$uc.Decision),
  CC_LRstat = c(test1_egarch11_std_dummy$cc.LRstat,
                test1_garch11_std_dummy$cc.LRstat,
                test1_gjr11_std_dummy$cc.LRstat),
  CC_pvalue = c(test1_egarch11_std_dummy$cc.LRp,
                test1_garch11_std_dummy$cc.LRp,
                test1_gjr11_std_dummy$cc.LRp),
  CC_Decision = c(test1_egarch11_std_dummy$cc.Decision,
                  test1_garch11_std_dummy$cc.Decision,
                  test1_gjr11_std_dummy$cc.Decision)
)

write.xlsx(
  var1_backtest_comparison,
  file = file.path(tables_dir, "VaR_1pct_backtest_comparison_std.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "VaR_1pct_backtest_comparison_std.txt"))
cat("COMPARATIVE 1% VaR BACKTEST RESULTS - Student-t\n\n")
print(var1_backtest_comparison)
sink()

# 15. VaR backtesting at the 5% level
test5_egarch11_std_dummy <- VaRTest(
  alpha = 0.05,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var5_egarch11_std_dummy)
)

test5_garch11_std_dummy <- VaRTest(
  alpha = 0.05,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var5_garch11_std_dummy)
)

test5_gjr11_std_dummy <- VaRTest(
  alpha = 0.05,
  actual = as.numeric(realizedret),
  VaR = as.numeric(var5_gjr11_std_dummy)
)

sink(file.path(output_dir, "VaR_tests_5pct_all_models_std.txt"))
cat("5% VaR Backtest Results - Student-t Models\n\n")

cat("===== EGARCH(1,1) =====\n")
print(test5_egarch11_std_dummy)

cat("\n===== GARCH(1,1) =====\n")
print(test5_garch11_std_dummy)

cat("\n===== GJR-GARCH(1,1) =====\n")
print(test5_gjr11_std_dummy)
sink()

var5_backtest_comparison <- data.frame(
  Model = c("EGARCH(1,1)",
            "GARCH(1,1)",
            "GJR-GARCH(1,1)"),
  Expected_Exceed = c(test5_egarch11_std_dummy$expected.exceed,
                      test5_garch11_std_dummy$expected.exceed,
                      test5_gjr11_std_dummy$expected.exceed),
  Actual_Exceed = c(test5_egarch11_std_dummy$actual.exceed,
                    test5_garch11_std_dummy$actual.exceed,
                    test5_gjr11_std_dummy$actual.exceed),
  UC_LRstat = c(test5_egarch11_std_dummy$uc.LRstat,
                test5_garch11_std_dummy$uc.LRstat,
                test5_gjr11_std_dummy$uc.LRstat),
  UC_pvalue = c(test5_egarch11_std_dummy$uc.LRp,
                test5_garch11_std_dummy$uc.LRp,
                test5_gjr11_std_dummy$uc.LRp),
  UC_Decision = c(test5_egarch11_std_dummy$uc.Decision,
                  test5_garch11_std_dummy$uc.Decision,
                  test5_gjr11_std_dummy$uc.Decision),
  CC_LRstat = c(test5_egarch11_std_dummy$cc.LRstat,
                test5_garch11_std_dummy$cc.LRstat,
                test5_gjr11_std_dummy$cc.LRstat),
  CC_pvalue = c(test5_egarch11_std_dummy$cc.LRp,
                test5_garch11_std_dummy$cc.LRp,
                test5_gjr11_std_dummy$cc.LRp),
  CC_Decision = c(test5_egarch11_std_dummy$cc.Decision,
                  test5_garch11_std_dummy$cc.Decision,
                  test5_gjr11_std_dummy$cc.Decision)
)

write.xlsx(
  var5_backtest_comparison,
  file = file.path(tables_dir, "VaR_5pct_backtest_comparison_std.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "VaR_5pct_backtest_comparison_std.txt"))
cat("COMPARATIVE 5% VaR BACKTEST RESULTS - Student-t\n\n")
print(var5_backtest_comparison)
sink()

# 16. Summary table of VaR backtests
backtest_summary_all <- data.frame(
  Model = rep(c("EGARCH(1,1)",
                "GARCH(1,1)",
                "GJR-GARCH(1,1)"), each = 2),
  VaR_Level = rep(c("1%", "5%"), times = 3),
  Expected_Exceed = c(
    test1_egarch11_std_dummy$expected.exceed, test5_egarch11_std_dummy$expected.exceed,
    test1_garch11_std_dummy$expected.exceed,  test5_garch11_std_dummy$expected.exceed,
    test1_gjr11_std_dummy$expected.exceed,    test5_gjr11_std_dummy$expected.exceed
  ),
  Actual_Exceed = c(
    test1_egarch11_std_dummy$actual.exceed, test5_egarch11_std_dummy$actual.exceed,
    test1_garch11_std_dummy$actual.exceed,  test5_garch11_std_dummy$actual.exceed,
    test1_gjr11_std_dummy$actual.exceed,    test5_gjr11_std_dummy$actual.exceed
  ),
  UC_LRstat = c(
    test1_egarch11_std_dummy$uc.LRstat, test5_egarch11_std_dummy$uc.LRstat,
    test1_garch11_std_dummy$uc.LRstat,  test5_garch11_std_dummy$uc.LRstat,
    test1_gjr11_std_dummy$uc.LRstat,    test5_gjr11_std_dummy$uc.LRstat
  ),
  UC_pvalue = c(
    test1_egarch11_std_dummy$uc.LRp, test5_egarch11_std_dummy$uc.LRp,
    test1_garch11_std_dummy$uc.LRp,  test5_garch11_std_dummy$uc.LRp,
    test1_gjr11_std_dummy$uc.LRp,    test5_gjr11_std_dummy$uc.LRp
  ),
  UC_Decision = c(
    test1_egarch11_std_dummy$uc.Decision, test5_egarch11_std_dummy$uc.Decision,
    test1_garch11_std_dummy$uc.Decision,  test5_garch11_std_dummy$uc.Decision,
    test1_gjr11_std_dummy$uc.Decision,    test5_gjr11_std_dummy$uc.Decision
  ),
  CC_LRstat = c(
    test1_egarch11_std_dummy$cc.LRstat, test5_egarch11_std_dummy$cc.LRstat,
    test1_garch11_std_dummy$cc.LRstat,  test5_garch11_std_dummy$cc.LRstat,
    test1_gjr11_std_dummy$cc.LRstat,    test5_gjr11_std_dummy$cc.LRstat
  ),
  CC_pvalue = c(
    test1_egarch11_std_dummy$cc.LRp, test5_egarch11_std_dummy$cc.LRp,
    test1_garch11_std_dummy$cc.LRp,  test5_garch11_std_dummy$cc.LRp,
    test1_gjr11_std_dummy$cc.LRp,    test5_gjr11_std_dummy$cc.LRp
  ),
  CC_Decision = c(
    test1_egarch11_std_dummy$cc.Decision, test5_egarch11_std_dummy$cc.Decision,
    test1_garch11_std_dummy$cc.Decision,  test5_garch11_std_dummy$cc.Decision,
    test1_gjr11_std_dummy$cc.Decision,    test5_gjr11_std_dummy$cc.Decision
  )
)

write.xlsx(
  backtest_summary_all,
  file = file.path(tables_dir, "VaR_backtests_summary_all_models_std.xlsx"),
  rowNames = FALSE
)

sink(file.path(output_dir, "VaR_backtests_summary_all_models_std.txt"))
cat("FINAL SUMMARY OF VaR BACKTESTS FOR ALL Student-t MODELS\n\n")
print(backtest_summary_all)
sink()

# 16A. Overall model comparison and final ranking
# Overall ranking combines:
# 1. In-sample fit.
# 2. Volatility forecasting performance measured by QLIKE.
# 3. Risk forecasting performance measured by VaR backtests.
# Lower final weighted score indicates a better overall model.
# Component weights:
# 20% In-sample fit.
# 40% Volatility forecasting, measured by QLIKE.
# 40% Risk forecasting, measured by VaR backtests.

# Helper ranking functions.
rank_asc <- function(x) rank(x, ties.method = "min", na.last = "keep")
rank_desc <- function(x) rank(-x, ties.method = "min", na.last = "keep")

# Build the raw comparison table.
insample_metrics <- model_comparison[, c(
  "Model", "LogLikelihood", "AIC", "BIC", "Shibata", "HannanQuinn"
)]

qlike_metrics <- qlike_summary[, c(
  "Model", "Mean_QLIKE", "Median_QLIKE", "SD_QLIKE", "Min_QLIKE", "Max_QLIKE"
)]

var1_metrics <- data.frame(
  Model = var1_backtest_comparison$Model,
  VaR1_Expected_Exceed = var1_backtest_comparison$Expected_Exceed,
  VaR1_Actual_Exceed = var1_backtest_comparison$Actual_Exceed,
  VaR1_Exceed_Error = abs(var1_backtest_comparison$Actual_Exceed - var1_backtest_comparison$Expected_Exceed),
  VaR1_UC_pvalue = var1_backtest_comparison$UC_pvalue,
  VaR1_CC_pvalue = var1_backtest_comparison$CC_pvalue,
  VaR1_UC_Pass = ifelse(var1_backtest_comparison$UC_pvalue >= 0.05, 1, 0),
  VaR1_CC_Pass = ifelse(var1_backtest_comparison$CC_pvalue >= 0.05, 1, 0)
)

var5_metrics <- data.frame(
  Model = var5_backtest_comparison$Model,
  VaR5_Expected_Exceed = var5_backtest_comparison$Expected_Exceed,
  VaR5_Actual_Exceed = var5_backtest_comparison$Actual_Exceed,
  VaR5_Exceed_Error = abs(var5_backtest_comparison$Actual_Exceed - var5_backtest_comparison$Expected_Exceed),
  VaR5_UC_pvalue = var5_backtest_comparison$UC_pvalue,
  VaR5_CC_pvalue = var5_backtest_comparison$CC_pvalue,
  VaR5_UC_Pass = ifelse(var5_backtest_comparison$UC_pvalue >= 0.05, 1, 0),
  VaR5_CC_Pass = ifelse(var5_backtest_comparison$CC_pvalue >= 0.05, 1, 0)
)

overall_raw <- Reduce(
  function(x, y) merge(x, y, by = "Model", all = TRUE),
  list(insample_metrics, qlike_metrics, var1_metrics, var5_metrics)
)

overall_raw$Total_VaR_Passes <- with(
  overall_raw,
  VaR1_UC_Pass + VaR1_CC_Pass + VaR5_UC_Pass + VaR5_CC_Pass
)

# Rank each model by individual metric.
overall_ranks <- overall_raw

overall_ranks$Rank_LogLikelihood <- rank_desc(overall_ranks$LogLikelihood)
overall_ranks$Rank_AIC           <- rank_asc(overall_ranks$AIC)
overall_ranks$Rank_BIC           <- rank_asc(overall_ranks$BIC)
overall_ranks$Rank_Shibata       <- rank_asc(overall_ranks$Shibata)
overall_ranks$Rank_HannanQuinn   <- rank_asc(overall_ranks$HannanQuinn)

overall_ranks$Rank_Mean_QLIKE <- rank_asc(overall_ranks$Mean_QLIKE)

overall_ranks$Rank_Total_VaR_Passes <- rank_desc(overall_ranks$Total_VaR_Passes)
overall_ranks$Rank_VaR1_Exceed_Error <- rank_asc(overall_ranks$VaR1_Exceed_Error)
overall_ranks$Rank_VaR5_Exceed_Error <- rank_asc(overall_ranks$VaR5_Exceed_Error)
overall_ranks$Rank_VaR1_UC_pvalue    <- rank_desc(overall_ranks$VaR1_UC_pvalue)
overall_ranks$Rank_VaR1_CC_pvalue    <- rank_desc(overall_ranks$VaR1_CC_pvalue)
overall_ranks$Rank_VaR5_UC_pvalue    <- rank_desc(overall_ranks$VaR5_UC_pvalue)
overall_ranks$Rank_VaR5_CC_pvalue    <- rank_desc(overall_ranks$VaR5_CC_pvalue)

# Create component scores.
overall_ranks$InSample_Score <- rowMeans(cbind(
  overall_ranks$Rank_LogLikelihood,
  overall_ranks$Rank_AIC,
  overall_ranks$Rank_BIC,
  overall_ranks$Rank_Shibata,
  overall_ranks$Rank_HannanQuinn
), na.rm = TRUE)

overall_ranks$Volatility_Score <- overall_ranks$Rank_Mean_QLIKE

overall_ranks$Risk_Score <- rowMeans(cbind(
  overall_ranks$Rank_Total_VaR_Passes,
  overall_ranks$Rank_VaR1_Exceed_Error,
  overall_ranks$Rank_VaR5_Exceed_Error,
  overall_ranks$Rank_VaR1_UC_pvalue,
  overall_ranks$Rank_VaR1_CC_pvalue,
  overall_ranks$Rank_VaR5_UC_pvalue,
  overall_ranks$Rank_VaR5_CC_pvalue
), na.rm = TRUE)

# Calculate the final weighted score.
overall_ranks$Final_Weighted_Score <- 0.20 * overall_ranks$InSample_Score +
  0.40 * overall_ranks$Volatility_Score +
  0.40 * overall_ranks$Risk_Score

overall_ranks$Overall_Rank <- rank_asc(overall_ranks$Final_Weighted_Score)

overall_ranks <- overall_ranks[order(overall_ranks$Overall_Rank), ]

# Create the final compact ranking table.
final_model_ranking <- overall_ranks[, c(
  "Overall_Rank",
  "Model",
  "Final_Weighted_Score",
  "InSample_Score",
  "Volatility_Score",
  "Risk_Score",
  "LogLikelihood",
  "AIC",
  "BIC",
  "Shibata",
  "HannanQuinn",
  "Mean_QLIKE",
  "Total_VaR_Passes",
  "VaR1_Exceed_Error",
  "VaR5_Exceed_Error",
  "VaR1_UC_pvalue",
  "VaR1_CC_pvalue",
  "VaR5_UC_pvalue",
  "VaR5_CC_pvalue"
)]

final_model_ranking$Comment <- c(
  "Best overall model",
  "Second-best overall model",
  "Third-best overall model"
)

# Create separate rankings by evaluation purpose.
purpose_ranking <- data.frame(
  Model = overall_ranks$Model,
  InSample_Rank = rank_asc(overall_ranks$InSample_Score),
  Volatility_Forecasting_Rank = rank_asc(overall_ranks$Volatility_Score),
  Risk_Forecasting_Rank = rank_asc(overall_ranks$Risk_Score),
  Overall_Rank = overall_ranks$Overall_Rank
)

purpose_ranking <- purpose_ranking[order(purpose_ranking$Overall_Rank), ]

# Save final comparison reports.
sink(file.path(output_dir, "overall_model_ranking_std.txt"))
cat("OVERALL MODEL COMPARISON AND FINAL RANKING - Student-t Models\n")
cat("Lower Final_Weighted_Score = better overall model\n\n")
cat("Weights used:\n")
cat("- 20% In-sample fit\n")
cat("- 40% Volatility forecasting (QLIKE)\n")
cat("- 40% Risk forecasting (VaR backtests)\n\n")

cat("FINAL RANKING TABLE\n\n")
print(final_model_ranking)

cat("\n\nPURPOSE-SPECIFIC RANKING\n\n")
print(purpose_ranking)

cat("\n\nBEST OVERALL MODEL:\n")
cat(as.character(final_model_ranking$Model[1]), "\n")
sink()

# Save the final comparison workbook.
write.xlsx(
  list(
    "Overall_Ranking" = final_model_ranking,
    "Purpose_Ranking" = purpose_ranking,
    "Raw_Metrics" = overall_raw,
    "Metric_Ranks" = overall_ranks
  ),
  file = file.path(tables_dir, "overall_model_ranking_std.xlsx"),
  rowNames = FALSE
)

cat("\n=====================================\n")
cat("FINAL MODEL RANKING COMPLETED\n")
cat("Best overall model:\n")
cat(final_model_ranking$Model[1], "\n")
cat("Excel file saved at:\n")
cat(file.path(tables_dir, "overall_model_ranking_std.xlsx"), "\n")
cat("=====================================\n")

# 17. Save workspace
save.image(file = file.path(output_dir, "workspace_SP500_45_with_VIX_dummy_std_1pct_5pct.RData"))

cat("Script execution completed successfully.\n")
cat("VaR was calculated and backtested at the 1% and 5% levels.\n")
cat("QLIKE and the overall model ranking were also calculated.\n")
cat("All results were saved to:\n")
cat(output_dir, "\n")

convergence(egarch11_std_dummy)
