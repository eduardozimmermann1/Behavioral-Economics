# ============================================================
# Decisive Language Intensity (DLI) and Inflation Expectations
# Households (Michigan Survey) vs Professionals (SPF)
# Quarterly pipeline with Newey-West inference and local projections
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---------- Packages ----------
pkgs <- c(
  "tidyverse", "stringr", "lubridate",
  "xml2", "rvest",
  "quanteda", "readr",
  "readxl", "janitor",
  "sandwich", "broom", "scales"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(lubridate)
  library(xml2)
  library(rvest)
  library(quanteda)
  library(readr)
  library(readxl)
  library(janitor)
  library(sandwich)
  library(broom)
  library(scales)
})

# ============================================================
# USER PATHS
# ============================================================
BASE_DIR <- "C:/Users/eduar/Desktop/Arquivos Eduardo/PRINCIPAIS DOCUMENTOS/PHD HSG/Macroeconomics/Seminar paper"
STATEMENTS_DIR <- file.path(BASE_DIR, "fomc_statements_html")

MICHIGAN_XLSX <- file.path(BASE_DIR, "sca-table32-on-2025-Dec-18.tsv - michigan mensal mean.xlsx")
SPF_XLSX <- file.path(BASE_DIR, "Mean_CPI_Level_ philadelphiafed trimestral mean.xlsx")
FOMC_SURPRISES_CSV <- file.path(BASE_DIR, "fomc_surprises.csv")

# ---------- Outputs ----------
OUT_DLI_DAILY_CSV         <- file.path(BASE_DIR, "fomc_dli_daily.csv")
OUT_DLI_QUARTERLY_CSV     <- file.path(BASE_DIR, "fomc_dli_quarterly.csv")
OUT_SURPRISES_Q_CSV       <- file.path(BASE_DIR, "fomc_surprises_quarterly.csv")
OUT_PANEL_QUARTERLY       <- file.path(BASE_DIR, "panel_quarterly_dli_michigan_spf.csv")

OUT_REG_RESULTS_CSV       <- file.path(BASE_DIR, "reg_results_seminar.csv")
OUT_REG_COEFS_CSV         <- file.path(BASE_DIR, "reg_coefs_seminar.csv")
OUT_REG_RESULTS_SURP_CSV  <- file.path(BASE_DIR, "reg_results_seminar_with_surprises.csv")
OUT_REG_COEFS_SURP_CSV    <- file.path(BASE_DIR, "reg_coefs_seminar_with_surprises.csv")
OUT_REG_COEFS_SHOCK_CSV   <- file.path(BASE_DIR, "reg_coefs_shock.csv")

OUT_REG_RESULTS_HET_CSV   <- file.path(BASE_DIR, "reg_results_heterogeneity_michigan_vs_spf.csv")
OUT_REG_COEFS_HET_CSV     <- file.path(BASE_DIR, "reg_coefs_heterogeneity_michigan_vs_spf.csv")

FIG_DIR <- file.path(BASE_DIR, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

FIG_TS_PANEL              <- file.path(FIG_DIR, "fig_ts_panel_michigan_spf_dli.png")
FIG_TS_PANEL_SURP         <- file.path(FIG_DIR, "fig_ts_panel_with_surprises.png")
FIG_SCATTER_LEVEL         <- file.path(FIG_DIR, "fig_scatter_level_michigan.png")
FIG_SCATTER_LEVEL_SURP    <- file.path(FIG_DIR, "fig_scatter_level_michigan_net_of_surprises.png")
FIG_SCATTER_DY            <- file.path(FIG_DIR, "fig_scatter_delta_michigan.png")
FIG_IRF_LOCALPROJ         <- file.path(FIG_DIR, "fig_localproj_irf_michigan.png")
FIG_IRF_LOCALPROJ_SURP    <- file.path(FIG_DIR, "fig_localproj_irf_michigan_with_surprises.png")
FIG_IRF_LOCALPROJ_SHOCK   <- file.path(FIG_DIR, "fig_localproj_irf_shock.png")
FIG_COEF_DLI              <- file.path(FIG_DIR, "fig_coef_dli_michigan.png")
FIG_COEF_DLI_SURP         <- file.path(FIG_DIR, "fig_coef_dli_michigan_with_surprises.png")
FIG_COEF_HET              <- file.path(FIG_DIR, "fig_coef_heterogeneity_michigan_vs_spf.png")

RUN_SHOCK <- TRUE

# ============================================================
# Helpers
# ============================================================
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

to_quarter_start <- function(d) {
  make_date(year(d), (quarter(d) - 1) * 3 + 1, 1)
}

hac_df_nw <- function(model, lag = 4, model_name = "model") {
  vc <- sandwich::NeweyWest(model, lag = lag, prewhite = FALSE, adjust = TRUE)
  b <- coef(model)
  se <- sqrt(diag(vc))
  tval <- b / se
  pval <- 2 * pt(abs(tval), df = df.residual(model), lower.tail = FALSE)

  tibble(
    model = model_name,
    term = names(b),
    estimate = as.numeric(b),
    std_error = as.numeric(se),
    t_value = as.numeric(tval),
    p_value = as.numeric(pval),
    n = nobs(model),
    r2 = summary(model)$r.squared,
    adj_r2 = summary(model)$adj.r.squared
  )
}

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
)

# ============================================================
# FIGURE LABEL STANDARDIZATION
# ============================================================

LABEL_HOUSEHOLDS <- "Household inflation expectations (Michigan Survey)"
LABEL_PROFESSIONALS <- "Professional inflation expectations (SPF)"
LABEL_DLI <- "Decisive Language Intensity (DLI, standardized)"
LABEL_DLI_SHOCK <- "DLI shock proxy (AR(1) innovation)"

# ============================================================
# (FIGURES ONLY – corrected labels)
# ============================================================

# Example: Time-series panel
p_ts <- ggplot(ts_df, aes(x = date_q, y = value, color = series)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Decisive Language Intensity and Inflation Expectations",
    x = NULL,
    y = "Standardized units",
    color = NULL
  )

# ALL OTHER FIGURES USE:
# - LABEL_HOUSEHOLDS
# - LABEL_PROFESSIONALS
# - LABEL_DLI
# - LABEL_DLI_SHOCK
# consistently in titles and axes

# ============================================================
# END SCRIPT
# ============================================================
