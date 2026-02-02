
# Reproduzierbares Setup
rm(list = ls()); cat("\014")

# Pakete für Datenaufbereitung, Zeitreihen und VAR-Schätzung
library(readr)
library(lubridate)
library(zoo)
library(vars)
library(ggplot2)
library(dplyr)

# Pfad zur Quartalsdaten-Datei (Deutschland
panel_path <- "C:/Bachelor/ErgebnisseNeu/quarterly_germany_VAR_with_taxes(1).csv"

# Einlesen der Daten und Umwandlung des Quartalsindexes
df <- read_csv(panel_path, show_col_types = FALSE)
out_dir <- "C:/Bachelor/ErgebnisseNeu/Output"
if (!dir.exists(out_dir)) dir.create(out_dir)
OUT <- function(fname) file.path(out_dir, fname)


df <- read_csv(panel_path, show_col_types = FALSE) %>%
  mutate(quarter = as.yearqtr(quarter, format = "%YQ%q")) %>%
  arrange(quarter)

# Zusammenstellung des VAR-Systems
# Reihenfolge entspricht der Cholesky-Identifikation
Y <- df %>%
  select(quarter,
         gov_g_q_dlog_pct,
         tax_g_q_dlog_pct,
         gdp_g_q_dlog_pct,
         infl_q_dlog_pct,
         euribor3m_q_avg) %>%
  rename(
    gov_g = gov_g_q_dlog_pct,
    tax_g = tax_g_q_dlog_pct,
    gdp_g = gdp_g_q_dlog_pct,
    infl  = infl_q_dlog_pct,
    i     = euribor3m_q_avg
  ) %>%
  arrange(quarter)


# Definition des Stichprobenzeitraums
library(dplyr)
start_q <- Y$quarter[1]
end_q   <- Y$quarter[nrow(Y)]
start_year <- as.integer(format(start_q, "%Y"))
start_quarter <- as.integer(cycle(start_q))  # 1..4


# Umwandlung in ts-Objekt (Quartalsfrequenz)
Y_ts <- ts(
  Y %>% select(gov_g, tax_g, gdp_g, infl, i),
  start = c(start_year, start_quarter),
  frequency = 4
)

cat("Sample:", format(start_q), "to", format(end_q), "N =", nrow(Y), "\n")
cat("Variables (Cholesky order): gov_g, tax_g, gdp_g, infl, i\n")


# Lag-Längen-Selektion anhand gängiger Informationskriterien
lag_max <- 8
lag_sel <- VARselect(Y_ts, lag.max = lag_max, type = "const")

crit <- data.frame(
  lag  = 1:lag_max,
  AIC  = lag_sel$criteria["AIC(n)", ],
  BIC  = lag_sel$criteria["SC(n)", ],
  HQIC = lag_sel$criteria["HQ(n)", ],
  FPE  = lag_sel$criteria["FPE(n)", ]
)


# Übersicht der Informationskriterien
print(lag_sel$selection)

# Übersicht der Informationskriterien
print(lag_sel$selection)
cat("\nLag criteria (head):\n")
print(head(crit, 8))

write.csv(crit, OUT("var_lag_criteria_quarterly_with_taxes.csv"), row.names = FALSE)


# Schätzung des reduzierten VAR-Modells
p <- 4
var_fit <- VAR(Y_ts, p = p, type = "const")


print(summary(var_fit))
eqs <- var_fit$varresult   

names(eqs) 

coef_all <- do.call(rbind, lapply(names(eqs), function(eq) {
  fit <- eqs[[eq]]                
  
  ct <- coef(summary(fit))         
  ci <- confint(fit, level = 0.95) 
  
  data.frame(
    equation = eq,
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"],
    t = ct[, "t value"],
    p = ct[, "Pr(>|t|)"],
    ci_low = ci[, 1],
    ci_high = ci[, 2],
    row.names = NULL
  )
}))

View(coef_all)
write.csv(coef_all, OUT("var_coeff_fulltable_all_equations.csv"), row.names = FALSE)

# Stabilitätsprüfung (Eigenwerte des Companion-Matrix)
roots_mod <- Mod(roots(var_fit))
is_stable <- all(roots_mod < 1)


# Residuen
whiteness_12 <- serial.test(var_fit, lags.pt = 12, type = "PT.asymptotic")
whiteness_16 <- serial.test(var_fit, lags.pt = 16, type = "PT.asymptotic")
norm_test <- normality.test(var_fit)
arch_res  <- arch.test(var_fit, lags.multi = 12)


cat("Stable (all roots modulus > 1):", is_stable, "\n")
cat("Roots modulus (first 10):\n")
print(head(roots_mod, 10))
cat("\nWhiteness (nlags=12):\n"); print(whiteness_12)
cat("\nWhiteness (nlags=16):\n"); print(whiteness_16)
cat("\nNormality:\n"); print(norm_test)
cat("\nARCH (optional):\n"); print(arch_res)



h <- 20
nboot <- 2000


# Orthogonalisierte Impuls-Antwort-Funktionen mit Bootstrap-Konfidenzintervall
irf_obj <- irf(
  var_fit,
  n.ahead = h,
  ortho = TRUE,
  boot = TRUE,
  runs = nboot,
  ci = 0.95
)

plot(irf_obj)

make_irf_long <- function(irf_obj, vars_order, h) {
  out <- list()
  for (shock in vars_order) {
    irf_mat <- irf_obj$irf[[shock]]
    lo_mat  <- irf_obj$Lower[[shock]]
    hi_mat  <- irf_obj$Upper[[shock]]
    horizons <- 0:h
    for (resp in colnames(irf_mat)) {
      out[[length(out) + 1]] <- data.frame(
        horizon_q = horizons,
        response = resp,
        shock = shock,
        irf = as.numeric(irf_mat[, resp]),
        lower = as.numeric(lo_mat[, resp]),
        upper = as.numeric(hi_mat[, resp])
      )
    }
  }
  bind_rows(out)
}

vars_order <- c("gov_g","tax_g","gdp_g","infl","i")
irf_long <- make_irf_long(irf_obj, vars_order, h)

# Speicherung der IRF-Ergebnisse zur weiteren Auswertung
write.csv(irf_long, OUT("var_quarterly_with_taxes_irf_orth_values.csv"), row.names = FALSE)

plot_irf_one <- function(irf_long, response, shock, file, title) {
  sub <- irf_long %>%
    filter(response == !!response, shock == !!shock) %>%
    arrange(horizon_q)
  
  p <- ggplot(sub, aes(x = horizon_q, y = irf)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_hline(yintercept = 0) +
    labs(x = "Horizon (quarters)", y = "Response", title = title) +
    theme_minimal()
  
  ggsave(file, p, width = 7.5, height = 4.5, dpi = 200)
}

fc <- predict(var_fit, n.ahead = 8, ci = 0.95)  
fc

# Grafische Darstellung ausgewählter Impuls-Antworten
plot_irf_one(irf_long, "gdp_g", "gov_g", OUT("irf_gdp_to_govshock.png"),
             "IRF: GDP growth response to a government spending shock (orthogonalised)")
plot_irf_one(irf_long, "infl",  "gov_g", OUT("irf_infl_to_govshock.png"),
             "IRF: Inflation response to a government spending shock (orthogonalised)")

plot_irf_one(irf_long, "gdp_g", "tax_g", OUT("irf_gdp_to_taxshock.png"),
             "IRF: GDP growth response to a taxes shock (orthogonalised)")
plot_irf_one(irf_long, "infl",  "tax_g", OUT("irf_infl_to_taxshock.png"),
             "IRF: Inflation response to a taxes shock (orthogonalised)")

plot_irf_one(irf_long, "gdp_g", "i", OUT("irf_gdp_to_rateshock.png"),
             "IRF: GDP growth response to an interest-rate shock (orthogonalised)")
plot_irf_one(irf_long, "infl",  "i", OUT("irf_infl_to_rateshock.png"),
             "IRF: Inflation response to an interest-rate shock (orthogonalised)")
plot_irf_one(irf_long, "i",     "i", OUT("irf_rate_to_rateshock.png"),
             "IRF: Interest rate response to an interest-rate shock (orthogonalised)")

cat("DONE. Outputs saved in:", normalizePath(out_dir), "\n")
