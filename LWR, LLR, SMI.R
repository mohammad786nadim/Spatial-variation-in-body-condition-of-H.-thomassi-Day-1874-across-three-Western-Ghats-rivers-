############################################################
## MASTER SCRIPT: LWR + CONDITION FACTORS + ANCOVA + LLR
############################################################

############################################################
## MASTER: LWR PLOTS FOR ALL 3 RIVERS (ONE BY ONE + COMBINED)
## Model: W = a L^b  (fit on log10 scale, back-transformed)
## Produces:
##   1) Three separate plots (Netravati, Periyar, Kallada)
##   2) One combined plot (all rivers in one figure; separate lines)
############################################################

library(dplyr)
library(ggplot2)
library(emmeans)
install.packages("FSA", repos = "https://cloud.r-project.org")
install.packages("smatr", repos = "https://cloud.r-project.org")

library(FSA)
library(smatr)

theme_set(theme_bw())
getwd()
setwd("E:/Desktop/PhD/Paper works/PhD/Length-Weight/Final/For resubmission")

# ---------- 1) PATH + READ FILES ----------
path_lwr <- "E:/Desktop/PhD/Paper works/PhD/Length-Weight/Final/For resubmission"

# Edit filenames if your LW files have different names
netravati_raw <- read.csv(file.path(path_lwr, "Netravati.csv"))
periyar_raw   <- read.csv(file.path(path_lwr, "Periyar.csv"))
kallada_raw   <- read.csv(file.path(path_lwr, "Kallada.csv"))

# ---------- 2) STANDARDIZE COLUMN NAMES (Length + Weight) ----------
standardize_lw_cols <- function(df) {
  
  length_col <- names(df)[grep("TL|Total|Length", names(df), ignore.case = TRUE)][1]
  weight_col <- names(df)[grep("^W$|Weight|Wt|BW|Body", names(df), ignore.case = TRUE)][1]
  
  if (any(is.na(c(length_col, weight_col)))) {
    stop("Could not detect Length/Weight columns automatically. Rename columns or edit grep rules.")
  }
  
  df %>%
    transmute(
      L = as.numeric(.data[[length_col]]),
      W = as.numeric(.data[[weight_col]])
    ) %>%
    filter(is.finite(L), is.finite(W), L > 0, W > 0)
}

netravati <- standardize_lw_cols(netravati_raw) %>% mutate(River = "Netravati")
periyar   <- standardize_lw_cols(periyar_raw)   %>% mutate(River = "Periyar")
kallada   <- standardize_lw_cols(kallada_raw)   %>% mutate(River = "Kallada")

# ---------- 3) FIT LWR (LOG MODEL) + CREATE CURVE DF ----------
fit_lwr_curve <- function(df, river_name, n_points = 200) {
  
  # log10(W) = log10(a) + b log10(L)
  m <- lm(log10(W) ~ log10(L), data = df)
  
  log_a <- coef(m)[1]
  b     <- coef(m)[2]
  a     <- 10^log_a
  
  r2 <- summary(m)$r.squared
  
  Lseq <- seq(min(df$L), max(df$L), length.out = n_points)
  Wpred <- a * (Lseq^b)
  
  curve <- data.frame(
    River = river_name,
    L = Lseq,
    W = Wpred,
    a = a,
    b = b,
    r2 = r2
  )
  
  list(model = m, curve = curve, a = a, b = b, r2 = r2)
}

fit_N <- fit_lwr_curve(netravati, "Netravati")
fit_P <- fit_lwr_curve(periyar,   "Periyar")
fit_K <- fit_lwr_curve(kallada,   "Kallada")

# ---------- 4) SEPARATE PLOTS (LIKE YOUR UPLOADED FIGURE) ----------
plot_lwr_single <- function(df, curve_df, river_name, a, b, r2) {
  
  eq <- paste0("W = ", round(a, 4), "L^", round(b, 3),
               "   (R² = ", round(r2, 3), ")")
  
  ggplot(df, aes(x = L, y = W)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(data = curve_df, aes(x = L, y = W),
              colour = "red", linewidth = 1.2) +
    labs(
      title = paste0("Length–Weight Relationship: ", river_name),
      subtitle = eq,
      x = "Total Length (cm)",
      y = "Weight (g)"
    ) +
    theme_bw()
}

p_N <- plot_lwr_single(netravati, fit_N$curve, "Netravati", fit_N$a, fit_N$b, fit_N$r2)
p_P <- plot_lwr_single(periyar,   fit_P$curve, "Periyar",   fit_P$a, fit_P$b, fit_P$r2)
p_K <- plot_lwr_single(kallada,   fit_K$curve, "Kallada",   fit_K$a, fit_K$b, fit_K$r2)

print(p_N)
print(p_P)
print(p_K)

# ---------- 5) COMBINED PLOT (ALL 3 RIVERS IN ONE) ----------
lw_all <- bind_rows(netravati, periyar, kallada) %>%
  mutate(River = factor(River, levels = c("Kallada", "Netravati", "Periyar")))

curve_all <- bind_rows(fit_N$curve, fit_P$curve, fit_K$curve) %>%
  mutate(River = factor(River, levels = c("Kallada", "Netravati", "Periyar")))

p_all <- ggplot(lw_all, aes(x = L, y = W, color = River)) +
  geom_point(alpha = 0.45, size = 2) +
  geom_line(data = curve_all, aes(x = L, y = W, color = River),
            linewidth = 1.2) +
  labs(
    title = "Length–Weight Relationship by River",
    x = "Total Length (cm)",
    y = "Weight (g)"
  ) +
  theme_bw()

print(p_all)

# ---------- OPTIONAL: SAVE (TIFF 600 dpi) ----------
# ggsave(file.path(path_lwr, "Netravati_LWR.tiff"), p_N, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_lwr, "Periyar_LWR.tiff"),   p_P, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_lwr, "Kallada_LWR.tiff"),   p_K, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_lwr, "ALLRIVERS_LWR.tiff"),  p_all, dpi = 600, width = 7, height = 5)

############################################################
## OPTIONAL: PRINT PARAMETER TABLE (a, b, R2)
############################################################
param_table <- data.frame(
  River = c("Netravati", "Periyar", "Kallada"),
  a = c(fit_N$a, fit_P$a, fit_K$a),
  b = c(fit_N$b, fit_P$b, fit_K$b),
  R2 = c(fit_N$r2, fit_P$r2, fit_K$r2)
)

print(param_table)

# ================= 1. LIBRARIES =================
library(dplyr)
library(ggplot2)

theme_set(theme_bw())

# ================= 2. FILE PATH =================
path_lwr <- "E:/Desktop/PhD/Paper works/PhD/Length-Weight/Final/For resubmission"

file_kallada   <- file.path(path_lwr, "Kallada.csv")
file_periyar   <- file.path(path_lwr, "Periyar.csv")
file_netravati <- file.path(path_lwr, "Netravati.csv")

# ================= 3. READ DATA =================
Kallada   <- read.csv(file_kallada)
Periyar   <- read.csv(file_periyar)
Netravati <- read.csv(file_netravati)

# ================= 4. PREP FUNCTION =================
prep_lwr <- function(df, river_name,
                     L_col = "TL..cm.", W_col = "Weight..g.") {
  df %>%
    mutate(
      River = river_name,
      TL = .data[[L_col]],
      W  = .data[[W_col]],
      logL = log(TL),
      logW = log(W),
      
      # Fulton’s Condition Factor
      FultonK = 100 * (W / (TL^3))
    ) %>%
    filter(is.finite(logL), is.finite(logW), is.finite(FultonK))
}

Kallada   <- prep_lwr(Kallada, "Kallada")
Periyar   <- prep_lwr(Periyar, "Periyar")
Netravati <- prep_lwr(Netravati, "Netravati")

# ================= 5. LWR FUNCTION =================
run_lwr <- function(df) {
  
  m <- lm(logW ~ logL, data = df)
  
  b <- coef(m)[["logL"]]
  a <- exp(coef(m)[["(Intercept)"]])
  
  df2 <- df %>%
    mutate(
      W_hat = a * (TL^b),
      Kn = W / W_hat,
      std_resid = rstandard(m),
      cooksD = cooks.distance(m),
      is_out_resid = abs(std_resid) > 3,
      is_out_cook  = cooksD > (4 / nrow(df))
    )
  
  out_tbl <- df2 %>%
    filter(is_out_resid | is_out_cook) %>%
    select(River, TL, W, FultonK, Kn,
           std_resid, cooksD, is_out_resid, is_out_cook)
  
  list(model = m,
       data_flagged = df2,
       outliers = out_tbl,
       a = a, b = b)
}

res_K <- run_lwr(Kallada)
res_P <- run_lwr(Periyar)
res_N <- run_lwr(Netravati)

# ================= 6. MODEL SUMMARIES =================
summary(res_K$model); confint(res_K$model)
summary(res_P$model); confint(res_P$model)
summary(res_N$model); confint(res_N$model)

# Optional: test isometry (b = 3)
# hoCoef(res_K$model, "logL", 3)
# hoCoef(res_P$model, "logL", 3)
# hoCoef(res_N$model, "logL", 3)

# ================= 7. OUTLIER TABLES =================
res_K$outliers
res_P$outliers
res_N$outliers

# ================= 8. LOG-LOG PLOTS =================
plot_log <- function(res_obj) {
  df <- res_obj$data_flagged
  ggplot(df, aes(x = logL, y = logW)) +
    geom_point(alpha = 0.5) +
    geom_point(data = df %>% filter(is_out_resid | is_out_cook),
               colour = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, colour = "red") +
    labs(x = "log(TL)", y = "log(W)",
         title = paste(unique(df$River)))
}

print(plot_log(res_K))
print(plot_log(res_P))
print(plot_log(res_N))


# ================= 9. COMBINED LOG-LOG PLOT =================

# Combine all processed datasets
combined_log <- bind_rows(
  res_K$data_flagged,
  res_P$data_flagged,
  res_N$data_flagged
)

# Ensure River is factor
combined_log$River <- factor(combined_log$River)

# Plot all three regression lines in one figure
ggplot(combined_log, aes(x = logL, y = logW, color = River)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = c(
    "Kallada"   = "#2ca02c",
    "Periyar"   = "#009", 
    "Netravati" = "red"
  )) +
  labs(
    title = "Length Weight Relationship by River",
    x = "Log Total Length (cm)",
    y = "Log Weight (g)",
    color = "River"
  ) +
  theme_bw()

# ================= RAW SCALE COMBINED PLOT =================

ggplot(combined_log, aes(x = TL, y = W, color = River)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_color_manual(values = c(
    "Kallada"   = "#2ca02c",
    "Periyar"   = "#0072B2",
    "Netravati" = "red"
  )) +
  labs(
    title = "Length–Weight Relationship by River",
    x = "Total Length (cm)",
    y = "Weight (g)",
    color = "River"
  ) +
  theme_bw()

# ================= 9. CONDITION FACTOR SUMMARY =================

library(ggplot2)
library(dplyr)
library(tidyr)

# Combine flagged datasets
cond_all <- bind_rows(
  res_K$data_flagged,
  res_P$data_flagged,
  res_N$data_flagged
) %>%
  mutate(
    River = factor(River, levels = c("Kallada", "Periyar", "Netravati"))
  )

# Summary table
cond_summary <- cond_all %>%
  group_by(River) %>%
  summarise(
    n       = n(),
    mean_K  = mean(FultonK, na.rm = TRUE),
    sd_K    = sd(FultonK, na.rm = TRUE),
    mean_Kn = mean(Kn, na.rm = TRUE),
    sd_Kn   = sd(Kn, na.rm = TRUE)
  )

print(cond_summary)

# ---------------- Individual plots ----------------

# Fulton's K
ggplot(cond_all, aes(x = River, y = FultonK)) +
  geom_boxplot(fill = "#66c2a5", color = "black",
               linewidth = 0.35, outlier.shape = NA) +
  coord_cartesian(ylim = c(0.5, 1.55)) +
  labs(
    title = "Fulton's Condition Factor Across Rivers",
    x = "River",
    y = "Fulton's K"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3)
  )

# Kn
ggplot(cond_all, aes(x = River, y = Kn)) +
  geom_boxplot(fill = "#fc8d62", color = "black",
               linewidth = 0.35, outlier.shape = NA) +
  coord_cartesian(ylim = c(0.5, 1.55)) +
  labs(
    title = "Relative Condition Factor (Kn) Across Rivers",
    x = "River",
    y = "Relative condition factor (Kn)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3)
  )

# ---------------- Stats ----------------

kruskal.test(FultonK ~ River, data = cond_all)
dunnTest(FultonK ~ River, data = cond_all, method = "bh")

kruskal.test(Kn ~ River, data = cond_all)
dunnTest(Kn ~ River, data = cond_all, method = "bh")

# ---------------- Combined plot ----------------

cond_long <- cond_all %>%
  pivot_longer(
    cols = c(FultonK, Kn),
    names_to = "Factor",
    values_to = "Value"
  ) %>%
  mutate(
    Factor = factor(Factor,
                    levels = c("FultonK", "Kn"),
                    labels = c("K", "Kn"))
  )

factor_cols <- c(
  "K"  = "#66c2a5",
  "Kn" = "#fc8d62"
)

p_cond <- ggplot(cond_long, aes(x = River, y = Value, fill = Factor)) +
  geom_boxplot(
    width = 0.7,
    alpha = 0.75,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    color = "black",
    linewidth = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(group = Factor),
    position = position_dodge(width = 0.75),
    shape = 21,
    size = 3.2,
    stroke = 0.8,
    color = "black"
  ) +
  scale_fill_manual(values = factor_cols) +
  coord_cartesian(ylim = c(0.5, 1.25)) +
  labs(
    title = "Fulton's (K) and Relative (Kn) Condition Factors Across Rivers",
    x = "River",
    y = "Condition Factor",
    fill = "Factor"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10.5),
    plot.title = element_text(face = "bold", hjust = 0, size = 14.5),
    axis.title = element_text(size = 12.5),
    axis.text = element_text(size = 11.5),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.6),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_cond)

# Save
ggsave(
  "Condition_Factors_Combined.tiff",
  plot = p_cond,
  width = 9.8,
  height = 5.6,
  dpi = 300
)

# ================= 10. CONDITION FACTOR: MEDIAN + RANGE =================

cond_summary_full <- cond_all %>%
  group_by(River) %>%
  summarise(
    n = n(),
    
    # --- Fulton’s K ---
    mean_K   = mean(FultonK, na.rm = TRUE),
    sd_K     = sd(FultonK, na.rm = TRUE),
    min_K    = min(FultonK, na.rm = TRUE),
    median_K = median(FultonK, na.rm = TRUE),
    max_K    = max(FultonK, na.rm = TRUE),
    
    # --- Relative condition factor Kn ---
    mean_Kn   = mean(Kn, na.rm = TRUE),
    sd_Kn     = sd(Kn, na.rm = TRUE),
    min_Kn    = min(Kn, na.rm = TRUE),
    median_Kn = median(Kn, na.rm = TRUE),
    max_Kn    = max(Kn, na.rm = TRUE)
  )

print(cond_summary_full)

# ================= 10. LWR ANCOVA =================
all_data <- bind_rows(Kallada, Periyar, Netravati)
all_data$River <- factor(all_data$River)

ancova_model <- lm(logW ~ logL * River, data = all_data)
summary(ancova_model)

emtrends(ancova_model, pairwise ~ River, var = "logL")

ggplot(all_data, aes(logL, logW, color = River)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Length Weight Relationship by River")

############################################################
################ LLR SECTION ###############################
############################################################

############################################################
## MASTER SCRIPT: LLR with OUTLIER REMOVAL USING RESIDUALS
## Outliers removed using STUDENTIZED residuals from the
## initial regression (|rstudent| > cutoff), then model refit.
##
## Runs 3 conversions for each river:
##   1) TL ~ SL
##   2) SL ~ FL
##   3) FL ~ TL
##
## Output style: summary, equation, confint()
############################################################

library(dplyr)

# ---------- 1) PATH + READ FILES ----------
path_llr <- "E:/Desktop/PhD/Paper works/PhD/Length-Weight/Final/For resubmission"

netravati_raw <- read.csv(file.path(path_llr, "Netravati_LL.csv"))
periyar_raw   <- read.csv(file.path(path_llr, "Periyar_LL.csv"))
kallada_raw   <- read.csv(file.path(path_llr, "Kallada_LL.csv"))

# ---------- 2) STANDARDIZE COLUMN NAMES TO TL, SL, FL ----------
# Tries to detect columns containing TL/Total, SL/Standard, FL/Fork
standardize_ll_cols <- function(df) {
  
  tl_col <- names(df)[grep("^TL|Total", names(df), ignore.case = TRUE)][1]
  sl_col <- names(df)[grep("^SL|Standard", names(df), ignore.case = TRUE)][1]
  fl_col <- names(df)[grep("^FL|Fork", names(df), ignore.case = TRUE)][1]
  
  if (any(is.na(c(tl_col, sl_col, fl_col)))) {
    stop("Could not detect TL/SL/FL columns automatically. Please rename columns or edit grep rules.")
  }
  
  out <- df %>%
    transmute(
      TL = as.numeric(.data[[tl_col]]),
      SL = as.numeric(.data[[sl_col]]),
      FL = as.numeric(.data[[fl_col]])
    ) %>%
    # remove NA rows (important for lm + residuals)
    filter(!is.na(TL), !is.na(SL), !is.na(FL))
  
  out
}

netravati <- standardize_ll_cols(netravati_raw)
periyar   <- standardize_ll_cols(periyar_raw)
kallada   <- standardize_ll_cols(kallada_raw)

# ---------- 3) OUTLIER REMOVAL USING STUDENTIZED RESIDUALS ----------
# Fits model once, flags |rstudent| > cutoff, removes them, refits.
remove_outliers_by_residuals <- function(fish_data, response, predictor, cutoff = 3, verbose = TRUE) {
  
  # Keep only response & predictor columns, ensure numeric, drop NA
  df <- fish_data %>%
    transmute(
      y = as.numeric(.data[[response]]),
      x = as.numeric(.data[[predictor]])
    ) %>%
    filter(!is.na(y), !is.na(x))
  
  if (nrow(df) < 5) stop("Too few rows after NA removal for reliable regression/outlier removal.")
  
  # initial model
  m_raw <- lm(y ~ x, data = df)
  
  # studentized residuals
  stud <- rstudent(m_raw)
  
  # outlier rule
  flag <- abs(stud) > cutoff
  out_idx <- which(flag)
  
  outliers <- df[out_idx, , drop = FALSE]
  if (nrow(outliers) > 0) {
    outliers$stud_resid <- stud[out_idx]
  }
  
  # cleaned data
  df_clean <- df[!flag, , drop = FALSE]
  
  # refit
  m_clean <- lm(y ~ x, data = df_clean)
  
  if (verbose) {
    cat("\nResidual-based outlier removal: ", response, " ~ ", predictor, "\n", sep = "")
    cat("n (raw)     = ", nrow(df), "\n", sep = "")
    cat("n (removed) = ", length(out_idx), "\n", sep = "")
    cat("n (clean)   = ", nrow(df_clean), "\n", sep = "")
    cat("Rule used   = |studentized residual| > ", cutoff, "\n", sep = "")
  }
  
  list(
    model_raw = m_raw,
    model_clean = m_clean,
    outliers = outliers,
    cleaned_xy = df_clean
  )
}

# ---------- 4) RUN + PRINT IN YOUR STYLE (CLEAN MODEL USED) ----------
run_llr_print_residclean <- function(fish_data, river_name, response, predictor,
                                     cutoff = 3,
                                     show_outliers = TRUE) {
  
  # Run outlier removal
  res <- remove_outliers_by_residuals(
    fish_data = fish_data,
    response  = response,
    predictor = predictor,
    cutoff    = cutoff,
    verbose   = TRUE
  )
  
  cat("\n\n", river_name, "\n", sep = "")
  cat("Call:\n")
  
  # Optional: print outliers removed
  if (show_outliers) {
    if (nrow(res$outliers) > 0) {
      cat("\n--- Outliers removed (y,x shown; based on studentized residuals) ---\n")
      print(res$outliers)
    } else {
      cat("\n--- No outliers flagged by residual rule ---\n")
    }
  }
  
  # Print summary for CLEAN model
  llr_model <- res$model_clean
  print(summary(llr_model))
  
  # Equation (CLEAN)
  b0 <- coef(llr_model)[1]
  b1 <- coef(llr_model)[2]
  
  cat("\nLLR Equation (CLEAN): ", response, " = ",
      round(b0, 3), " + ", round(b1, 3), " * ", predictor, "\n", sep = "")
  
  cat("\nconfint(llr_model)\n")
  print(confint(llr_model))
  
  invisible(llr_model)
}

# ---------- 5) RUN FOR EACH RIVER ----------
# Choose your cutoff (common: 3)
RESID_CUTOFF <- 3

# (A) TL ~ SL
cat("\n================ TL ~ SL ================\n")
run_llr_print_residclean(netravati, "Netravati_LL", "TL", "SL", cutoff = RESID_CUTOFF)
run_llr_print_residclean(periyar,   "Periyar_LL",   "TL", "SL", cutoff = RESID_CUTOFF)
run_llr_print_residclean(kallada,   "Kallada_LL",   "TL", "SL", cutoff = RESID_CUTOFF)

# (B) SL ~ FL
cat("\n================ SL ~ FL ================\n")
run_llr_print_residclean(netravati, "Netravati", "SL", "FL", cutoff = RESID_CUTOFF)
run_llr_print_residclean(periyar,   "Periyar",   "SL", "FL", cutoff = RESID_CUTOFF)
run_llr_print_residclean(kallada,   "Kallada",   "SL", "FL", cutoff = RESID_CUTOFF)

# (C) FL ~ TL
cat("\n================ FL ~ TL ================\n")
run_llr_print_residclean(netravati, "Netravati", "FL", "TL", cutoff = RESID_CUTOFF)
run_llr_print_residclean(periyar,   "Periyar",   "FL", "TL", cutoff = RESID_CUTOFF)
run_llr_print_residclean(kallada,   "Kallada",   "FL", "TL", cutoff = RESID_CUTOFF)

############################################################
## MASTER SCRIPT (CONTINUATION): LLR PLOTS + ANALYSIS
## Outliers flagged AND REMOVED using STUDENTIZED residuals:
##   |rstudent(model_raw)| > RESID_CUTOFF
##
## - Points flagged as outliers are shown in BLUE
## - Regression line is fit on CLEAN data (outliers removed)
##
## Produces:
##  A) Separate plots for each river + each conversion
##  B) Combined plots (all 3 rivers) for each conversion
############################################################

library(dplyr)
library(ggplot2)

theme_set(theme_bw())

# ---------- 1) READ YOUR LLR FILES ----------
path_llr <- "E:/Desktop/PhD/Paper works/PhD/Length-Weight/Final/For resubmission"

netravati_raw <- read.csv(file.path(path_llr, "Netravati_LL.csv"))
periyar_raw   <- read.csv(file.path(path_llr, "Periyar_LL.csv"))
kallada_raw   <- read.csv(file.path(path_llr, "Kallada_LL.csv"))

# ---------- 2) STANDARDIZE COLUMN NAMES TO TL, SL, FL ----------
standardize_ll_cols <- function(df) {
  tl_col <- names(df)[grep("^TL|Total", names(df), ignore.case = TRUE)][1]
  sl_col <- names(df)[grep("^SL|Standard", names(df), ignore.case = TRUE)][1]
  fl_col <- names(df)[grep("^FL|Fork", names(df), ignore.case = TRUE)][1]
  
  if (any(is.na(c(tl_col, sl_col, fl_col)))) {
    stop("Could not detect TL/SL/FL columns automatically. Please rename columns or edit grep rules.")
  }
  
  df %>%
    transmute(
      TL = as.numeric(.data[[tl_col]]),
      SL = as.numeric(.data[[sl_col]]),
      FL = as.numeric(.data[[fl_col]])
    ) %>%
    filter(is.finite(TL), is.finite(SL), is.finite(FL))
}

netravati <- standardize_ll_cols(netravati_raw) %>% mutate(River = "Netravati")
periyar   <- standardize_ll_cols(periyar_raw)   %>% mutate(River = "Periyar")
kallada   <- standardize_ll_cols(kallada_raw)   %>% mutate(River = "Kallada")

# ---------- 3) SETTINGS ----------
RESID_CUTOFF <- 3   # outlier rule: |studentized residual| > 3

# ---------- 4) CORE: FLAG OUTLIERS (rstudent) + CLEAN LINE FIT ----------
# Returns:
#  - data_flagged: all points with outlier flag
#  - data_clean:   outliers removed
#  - model_raw:    model on all points
#  - model_clean:  model on clean points
flag_and_clean_rstudent <- function(df, response, predictor, cutoff = 3) {
  
  df2 <- df %>%
    mutate(
      y = as.numeric(.data[[response]]),
      x = as.numeric(.data[[predictor]])
    ) %>%
    filter(is.finite(y), is.finite(x))
  
  if (nrow(df2) < 5) stop("Too few rows for regression/outlier detection.")
  
  m_raw <- lm(y ~ x, data = df2)
  stud  <- rstudent(m_raw)
  
  df_flag <- df2 %>%
    mutate(
      stud_resid = stud,
      is_outlier = abs(stud_resid) > cutoff
    )
  
  df_clean <- df_flag %>%
    filter(!is_outlier)
  
  m_clean <- lm(y ~ x, data = df_clean)
  
  list(
    data_flagged = df_flag,
    data_clean   = df_clean,
    model_raw    = m_raw,
    model_clean  = m_clean
  )
}

# ---------- 5) PLOT FUNCTION: SEPARATE (PER RIVER) ----------
# - Blue points = outliers (based on rstudent from raw model)
# - Red line = regression on CLEAN data only
llr_plot_clean_with_blue_outliers <- function(df, river_name, response, predictor, cutoff = 3) {
  
  res <- flag_and_clean_rstudent(df, response, predictor, cutoff)
  
  # Equation from CLEAN model
  b0 <- coef(res$model_clean)[1]
  b1 <- coef(res$model_clean)[2]
  r2 <- summary(res$model_clean)$r.squared
  
  eq <- paste0(response, " = ", round(b0, 3), " + ", round(b1, 3), " × ", predictor,
               "   (R² clean = ", round(r2, 3), ")")
  
  p <- ggplot(res$data_flagged, aes(x = x, y = y)) +
    geom_point(alpha = 0.55, size = 2) +
    geom_point(data = res$data_flagged %>% filter(is_outlier),
               colour = "blue", size = 3) +
    geom_smooth(data = res$data_clean, aes(x = x, y = y),
                method = "lm", se = TRUE, colour = "red", linewidth = 1) +
    labs(
      title = paste0(river_name, ": ", response, " vs ", predictor),
      subtitle = paste0("Blue points = Outliers"),
      x = paste0(predictor, " (cm)"),
      y = paste0(response, " (cm)")
    ) +
    theme_bw()
  
  out_tbl <- res$data_flagged %>%
    filter(is_outlier) %>%
    transmute(
      River = river_name,
      Response = response,
      Predictor = predictor,
      TL = TL, SL = SL, FL = FL,
      stud_resid
    )
  
  list(plot = p, outliers = out_tbl, model_clean = res$model_clean, data_flagged = res$data_flagged)
}

############################################################
## 6) RUN ALL SEPARATE PLOTS + OUTLIER TABLES
############################################################

# ---- TL ~ SL
res_TL_SL_N <- llr_plot_clean_with_blue_outliers(netravati, "Netravati", "TL", "SL", cutoff = RESID_CUTOFF)
res_TL_SL_P <- llr_plot_clean_with_blue_outliers(periyar,   "Periyar",   "TL", "SL", cutoff = RESID_CUTOFF)
res_TL_SL_K <- llr_plot_clean_with_blue_outliers(kallada,   "Kallada",   "TL", "SL", cutoff = RESID_CUTOFF)

print(res_TL_SL_N$plot); print(res_TL_SL_P$plot); print(res_TL_SL_K$plot)
res_TL_SL_N$outliers; res_TL_SL_P$outliers; res_TL_SL_K$outliers

# ---- SL ~ FL
res_SL_FL_N <- llr_plot_clean_with_blue_outliers(netravati, "Netravati", "SL", "FL", cutoff = RESID_CUTOFF)
res_SL_FL_P <- llr_plot_clean_with_blue_outliers(periyar,   "Periyar",   "SL", "FL", cutoff = RESID_CUTOFF)
res_SL_FL_K <- llr_plot_clean_with_blue_outliers(kallada,   "Kallada",   "SL", "FL", cutoff = RESID_CUTOFF)

print(res_SL_FL_N$plot); print(res_SL_FL_P$plot); print(res_SL_FL_K$plot)
res_SL_FL_N$outliers; res_SL_FL_P$outliers; res_SL_FL_K$outliers

# ---- FL ~ TL
res_FL_TL_N <- llr_plot_clean_with_blue_outliers(netravati, "Netravati", "FL", "TL", cutoff = RESID_CUTOFF)
res_FL_TL_P <- llr_plot_clean_with_blue_outliers(periyar,   "Periyar",   "FL", "TL", cutoff = RESID_CUTOFF)
res_FL_TL_K <- llr_plot_clean_with_blue_outliers(kallada,   "Kallada",   "FL", "TL", cutoff = RESID_CUTOFF)

print(res_FL_TL_N$plot); print(res_FL_TL_P$plot); print(res_FL_TL_K$plot)
res_FL_TL_N$outliers; res_FL_TL_P$outliers; res_FL_TL_K$outliers


############################################################
## 7) COMBINED PLOTS (ALL 3 RIVERS) — BLUE OUTLIERS, CLEAN LINES
############################################################

ll_all <- bind_rows(kallada, netravati, periyar) %>%
  mutate(River = factor(River, levels = c("Kallada", "Netravati", "Periyar")))

# Flag outliers per river (rstudent) and also create a cleaned dataset per river
flag_clean_by_river_rstudent <- function(df, response, predictor, cutoff = 3) {
  
  df %>%
    group_by(River) %>%
    group_modify(~{
      tmp <- flag_and_clean_rstudent(.x, response, predictor, cutoff)
      # return flagged data (includes outlier flag) for plotting points
      tmp$data_flagged
    }) %>%
    ungroup()
}

# Build cleaned dataset per river for fitting clean regression lines
clean_only_by_river_rstudent <- function(df, response, predictor, cutoff = 3) {
  
  df %>%
    group_by(River) %>%
    group_modify(~{
      tmp <- flag_and_clean_rstudent(.x, response, predictor, cutoff)
      tmp$data_clean
    }) %>%
    ungroup()
}

plot_combined_clean_llr <- function(df, response, predictor, title_text, cutoff = 3) {
  
  df_flag  <- flag_clean_by_river_rstudent(df, response, predictor, cutoff)
  df_clean <- clean_only_by_river_rstudent(df, response, predictor, cutoff)
  
  ggplot(df_flag, aes(x = x, y = y, color = River)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_point(data = df_flag %>% filter(is_outlier),
               color = "blue", size = 3) +
    geom_smooth(data = df_clean, aes(x = x, y = y, color = River),
                method = "lm", se = FALSE, linewidth = 1) +
    labs(
      title = title_text,
      subtitle = paste0("Blue points = Outliers"),
      x = paste0(predictor, " (cm)"),
      y = paste0(response, " (cm)")
    ) +
    theme_bw()
}

p_TL_SL <- plot_combined_clean_llr(ll_all, response = "TL", predictor = "SL",
                                   title_text = "TL vs SL by River",
                                   cutoff = RESID_CUTOFF)

p_SL_FL <- plot_combined_clean_llr(ll_all, response = "SL", predictor = "FL",
                                   title_text = "SL vs FL by River ",
                                   cutoff = RESID_CUTOFF)

p_FL_TL <- plot_combined_clean_llr(ll_all, response = "FL", predictor = "TL",
                                   title_text = "FL vs TL by River ",
                                   cutoff = RESID_CUTOFF)

print(p_TL_SL)
print(p_SL_FL)
print(p_FL_TL)

############################################################
## OPTIONAL: SAVE ALL PLOTS (HIGH-RES TIFF, 600 dpi)
############################################################
# ggsave(file.path(path_llr, "Netravati_TL_vs_SL_clean_blueOut.tiff"), res_TL_SL_N$plot, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_llr, "Periyar_TL_vs_SL_clean_blueOut.tiff"),   res_TL_SL_P$plot, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_llr, "Kallada_TL_vs_SL_clean_blueOut.tiff"),   res_TL_SL_K$plot, dpi = 600, width = 6, height = 5)
#
# ggsave(file.path(path_llr, "Netravati_SL_vs_FL_clean_blueOut.tiff"), res_SL_FL_N$plot, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_llr, "Periyar_SL_vs_FL_clean_blueOut.tiff"),   res_SL_FL_P$plot, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_llr, "Kallada_SL_vs_FL_clean_blueOut.tiff"),   res_SL_FL_K$plot, dpi = 600, width = 6, height = 5)
#
# ggsave(file.path(path_llr, "Netravati_FL_vs_TL_clean_blueOut.tiff"), res_FL_TL_N$plot, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_llr, "Periyar_FL_vs_TL_clean_blueOut.tiff"),   res_FL_TL_P$plot, dpi = 600, width = 6, height = 5)
# ggsave(file.path(path_llr, "Kallada_FL_vs_TL_clean_blueOut.tiff"),   res_FL_TL_K$plot, dpi = 600, width = 6, height = 5)
#
# ggsave(file.path(path_llr, "ALLRIVERS_TL_vs_SL_clean_blueOut.tiff"), p_TL_SL, dpi = 600, width = 7, height = 5)
# ggsave(file.path(path_llr, "ALLRIVERS_SL_vs_FL_clean_blueOut.tiff"), p_SL_FL, dpi = 600, width = 7, height = 5)
# ggsave(file.path(path_llr, "ALLRIVERS_FL_vs_TL_clean_blueOut.tiff"), p_FL_TL, dpi = 600, width = 7, height = 5)

###########################ANCOVA######################################
# ============================================================
# MASTER CODE: LLR with residual-based outlier removal
# + river-wise model summaries
# + ANCOVA for slope comparison
# + pairwise differences in b values
# Works for:
#   (A) TL ~ SL
#   (B) SL ~ FL
#   (C) FL ~ TL
# ============================================================

# -----------------------------
# PACKAGES
# -----------------------------
library(dplyr)


# -----------------------------
# SETTINGS
# -----------------------------
RESID_CUTOFF <- 3

# -----------------------------
# FUNCTION 1: Clean outliers and fit model for one river
# -----------------------------
run_llr_print_residclean <- function(df, river_name, response, predictor, cutoff = 3) {
  
  # Keep only required columns that exist
  keep_vars <- unique(c(response, predictor, "TL", "SL", "FL"))
  keep_vars <- keep_vars[keep_vars %in% names(df)]
  
  df0 <- df %>%
    dplyr::select(all_of(keep_vars)) %>%
    na.omit()
  
  # Temporary generic names
  df0$y <- df0[[response]]
  df0$x <- df0[[predictor]]
  
  # Raw model
  model_raw <- lm(y ~ x, data = df0)
  
  # Studentized residuals
  df0$stud_resid <- rstudent(model_raw)
  
  # Outliers
  outliers <- df0 %>%
    dplyr::filter(abs(stud_resid) > cutoff)
  
  # Clean data
  df_clean <- df0 %>%
    dplyr::filter(abs(stud_resid) <= cutoff)
  
  # Final model
  llr_model <- lm(y ~ x, data = df_clean)
  sm <- summary(llr_model)
  ci <- confint(llr_model)
  
  a <- coef(llr_model)[1]
  b <- coef(llr_model)[2]
  
  cat("\n============================================================\n")
  cat("Residual-based outlier removal:", response, "~", predictor, "\n")
  cat("River         =", river_name, "\n")
  cat("n (raw)       =", nrow(df0), "\n")
  cat("n (removed)   =", nrow(outliers), "\n")
  cat("n (clean)     =", nrow(df_clean), "\n")
  cat("Rule used     = |studentized residual| >", cutoff, "\n")
  
  if (nrow(outliers) > 0) {
    cat("\n--- Outliers removed (based on studentized residuals) ---\n")
    out_cols <- intersect(c(response, predictor, "TL", "SL", "FL", "stud_resid"), names(outliers))
    print(outliers[, out_cols], row.names = FALSE)
  } else {
    cat("\nNo outliers removed.\n")
  }
  
  cat("\nCall:\n")
  print(llr_model)
  
  cat("\nLLR Equation (CLEAN):", response, "=", round(a, 3), "+", round(b, 3), "*", predictor, "\n")
  
  cat("\n95% CI:\n")
  print(ci)
  
  cat("\nModel summary stats:\n")
  cat("RSE                =", sm$sigma, "\n")
  cat("R-squared          =", sm$r.squared, "\n")
  cat("Adjusted R-squared =", sm$adj.r.squared, "\n")
  cat("F-statistic        =", sm$fstatistic[1], 
      "on", sm$fstatistic[2], "and", sm$fstatistic[3], "DF\n")
  cat("Model p-value      =", pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3], lower.tail = FALSE), "\n")
  
  return(list(
    river = river_name,
    response = response,
    predictor = predictor,
    raw_data = df0,
    clean_data = df_clean,
    outliers = outliers,
    model = llr_model,
    ci = ci,
    summary = sm
  ))
}

# -----------------------------
# FUNCTION 2: Run all rivers for one variable pair
# -----------------------------
run_all_rivers_llr <- function(netravati, periyar, kallada, response, predictor, cutoff = 3) {
  
  res_N <- run_llr_print_residclean(netravati, "Netravati", response, predictor, cutoff)
  res_P <- run_llr_print_residclean(periyar,   "Periyar",   response, predictor, cutoff)
  res_K <- run_llr_print_residclean(kallada,   "Kallada",   response, predictor, cutoff)
  
  clean_N <- res_N$clean_data %>% mutate(River = "Netravati")
  clean_P <- res_P$clean_data %>% mutate(River = "Periyar")
  clean_K <- res_K$clean_data %>% mutate(River = "Kallada")
  
  all_clean <- bind_rows(clean_N, clean_P, clean_K)
  all_clean$River <- factor(all_clean$River, levels = c("Netravati", "Periyar", "Kallada"))
  
  return(list(
    Netravati = res_N,
    Periyar   = res_P,
    Kallada   = res_K,
    all_clean = all_clean
  ))
}

# -----------------------------
# FUNCTION 3: ANCOVA + pairwise slope comparisons
# -----------------------------
run_pairwise_b_comparison <- function(all_clean, response, predictor) {
  
  ancova_model <- lm(y ~ x * River, data = all_clean)
  
  cat("\n============================================================\n")
  cat("ANCOVA for slope comparison:", response, "~", predictor, "across rivers\n")
  cat("============================================================\n")
  print(anova(ancova_model))
  
  cat("\n--- Estimated slopes (b values) for each river ---\n")
  slopes <- emtrends(ancova_model, ~ River, var = "x")
  print(summary(slopes, infer = TRUE))
  
  cat("\n--- Pairwise differences in b values ---\n")
  pw <- pairs(slopes)
  print(summary(pw, infer = TRUE))
  
  return(list(
    ancova_model = ancova_model,
    ancova_table = anova(ancova_model),
    slopes       = slopes,
    pairwise_b   = pw
  ))
}

# -----------------------------
# FUNCTION 4: Full master workflow for one variable pair
# -----------------------------
run_master_llr <- function(netravati, periyar, kallada, response, predictor, cutoff = 3) {
  
  cat("\n\n############################################################\n")
  cat("MASTER ANALYSIS FOR:", response, "~", predictor, "\n")
  cat("############################################################\n")
  
  # River-wise cleaning + fitting
  res_all <- run_all_rivers_llr(
    netravati = netravati,
    periyar   = periyar,
    kallada   = kallada,
    response  = response,
    predictor = predictor,
    cutoff    = cutoff
  )
  
  # ANCOVA + pairwise slope differences
  ancova_res <- run_pairwise_b_comparison(
    all_clean = res_all$all_clean,
    response  = response,
    predictor = predictor
  )
  
  # Collect outliers
  out_N <- res_all$Netravati$outliers %>%
    mutate(River = "Netravati", Response = response, Predictor = predictor)
  
  out_P <- res_all$Periyar$outliers %>%
    mutate(River = "Periyar", Response = response, Predictor = predictor)
  
  out_K <- res_all$Kallada$outliers %>%
    mutate(River = "Kallada", Response = response, Predictor = predictor)
  
  all_outliers <- bind_rows(out_N, out_P, out_K)
  
  if (nrow(all_outliers) > 0) {
    keep_cols <- intersect(c("River", "Response", "Predictor", "TL", "SL", "FL", "stud_resid"),
                           names(all_outliers))
    cat("\n--- All removed outliers ---\n")
    print(all_outliers[, keep_cols], row.names = FALSE)
  } else {
    cat("\nNo outliers removed across rivers.\n")
  }
  
  return(list(
    riverwise = res_all,
    ancova    = ancova_res,
    outliers  = all_outliers
  ))
}

# ============================================================
# RUN ALL ANALYSES
# ============================================================

# ------------------------------------------------------------
# (A) TL ~ SL
# ------------------------------------------------------------
res_TL_SL <- run_master_llr(
  netravati = netravati,
  periyar   = periyar,
  kallada   = kallada,
  response  = "TL",
  predictor = "SL",
  cutoff    = RESID_CUTOFF
)

# ------------------------------------------------------------
# (B) SL ~ FL
# ------------------------------------------------------------
res_SL_FL <- run_master_llr(
  netravati = netravati,
  periyar   = periyar,
  kallada   = kallada,
  response  = "SL",
  predictor = "FL",
  cutoff    = RESID_CUTOFF
)

# ------------------------------------------------------------
# (C) FL ~ TL
# ------------------------------------------------------------
res_FL_TL <- run_master_llr(
  netravati = netravati,
  periyar   = periyar,
  kallada   = kallada,
  response  = "FL",
  predictor = "TL",
  cutoff    = RESID_CUTOFF
)

# ============================================================
# OPTIONAL: DIRECT ACCESS TO RESULTS
# ============================================================

# -----------------------------
# TL ~ SL
# -----------------------------
summary(res_TL_SL$ancova$pairwise_b, infer = TRUE)
summary(res_TL_SL$ancova$slopes, infer = TRUE)
res_TL_SL$ancova$ancova_table
res_TL_SL$outliers

# -----------------------------
# SL ~ FL
# -----------------------------
summary(res_SL_FL$ancova$pairwise_b, infer = TRUE)
summary(res_SL_FL$ancova$slopes, infer = TRUE)
res_SL_FL$ancova$ancova_table
res_SL_FL$outliers

# -----------------------------
# FL ~ TL
# -----------------------------
summary(res_FL_TL$ancova$pairwise_b, infer = TRUE)
summary(res_FL_TL$ancova$slopes, infer = TRUE)
res_FL_TL$ancova$ancova_table
res_FL_TL$outliers

################################## SCALED-MASS INDEX #################################
######################################################################################
## MASTER CODE: SMI DIRECTLY FROM UPLOADED CSV FILES
## USING NATURAL LOG
######################################################################################

# ================= 1. LIBRARIES =================
library(dplyr)
library(ggplot2)


theme_set(theme_bw())

# ================= 2. READ THE UPLOADED CSV FILES =================
Netravati_clean <- read.csv("Netravati.csv", check.names = FALSE)
Periyar_clean   <- read.csv("Periyar.csv", check.names = FALSE)
Kallada_clean   <- read.csv("Kallada.csv", check.names = FALSE)

# ================= 3. CHECK SAMPLE SIZES =================
cat("\n===== CSV SAMPLE SIZE =====\n")
cat("Netravati =", nrow(Netravati_clean), "\n")
cat("Periyar   =", nrow(Periyar_clean), "\n")
cat("Kallada   =", nrow(Kallada_clean), "\n")

# ================= 4. INSPECT COLUMN NAMES =================
cat("\n===== COLUMN NAMES IN CSV DATA =====\n")
cat("Netravati:\n")
print(names(Netravati_clean))
cat("\nPeriyar:\n")
print(names(Periyar_clean))
cat("\nKallada:\n")
print(names(Kallada_clean))

# ==========================================================
# 5. FUNCTION TO CALCULATE SMI FROM CSV DATA ONLY
#    USING NATURAL LOG
# ==========================================================
calc_smi_from_clean <- function(df, river_name,
                                TL_col = NULL,
                                W_col  = NULL) {
  
  # Auto-detect columns if not supplied
  if (is.null(TL_col)) {
    if ("TL" %in% names(df)) {
      TL_col <- "TL"
    } else if ("TL..cm." %in% names(df)) {
      TL_col <- "TL..cm."
    } else if ("TL (cm)" %in% names(df)) {
      TL_col <- "TL (cm)"
    } else {
      stop(paste("TL column not found in", river_name))
    }
  }
  
  if (is.null(W_col)) {
    if ("W" %in% names(df)) {
      W_col <- "W"
    } else if ("Weight..g." %in% names(df)) {
      W_col <- "Weight..g."
    } else if ("Weight (g)" %in% names(df)) {
      W_col <- "Weight (g)"
    } else {
      stop(paste("Weight column not found in", river_name))
    }
  }
  
  df2 <- df %>%
    mutate(
      River = river_name,
      TL = .data[[TL_col]],
      W  = .data[[W_col]]
    ) %>%
    filter(is.finite(TL), is.finite(W), TL > 0, W > 0)
  
  # Mean total length from the same cleaned data
  L0 <- mean(df2$TL, na.rm = TRUE)
  
  # SMA using natural log
  sma_fit <- sma(log(W) ~ log(TL), data = df2)
  
  # Robust extraction of SMA slope
  cf <- coef(sma_fit)
  
  if (is.matrix(cf) || is.data.frame(cf)) {
    b_SMA <- as.numeric(cf[nrow(cf), ncol(cf)])
  } else {
    b_SMA <- as.numeric(cf[length(cf)])
  }
  
  # Calculate SMI
  df2 <- df2 %>%
    mutate(
      L0 = L0,
      b_SMA = b_SMA,
      SMI = W * (L0 / TL)^b_SMA
    )
  
  list(
    data = df2,
    L0 = L0,
    b_SMA = b_SMA,
    sma_fit = sma_fit
  )
}

# ================= 6. CALCULATE SMI =================
res_smi_N <- calc_smi_from_clean(Netravati_clean, "Netravati")
res_smi_P <- calc_smi_from_clean(Periyar_clean, "Periyar")
res_smi_K <- calc_smi_from_clean(Kallada_clean, "Kallada")

Netravati_smi <- res_smi_N$data
Periyar_smi   <- res_smi_P$data
Kallada_smi   <- res_smi_K$data

# ================= 7. COMBINE =================
smi_data <- bind_rows(Netravati_smi, Periyar_smi, Kallada_smi) %>%
  mutate(River = factor(River, levels = c("Netravati", "Periyar", "Kallada")))

# ================= 8. VERIFY COUNTS =================
cat("\n===== SAMPLE SIZE USED FOR SMI =====\n")
print(table(smi_data$River))

# ================= 9. PRINT SMI PARAMETERS =================
cat("\n===== SMI PARAMETERS =====\n")
cat("Netravati: L0 =", round(res_smi_N$L0, 4),
    ", b_SMA =", round(res_smi_N$b_SMA, 4), "\n")
cat("Periyar:   L0 =", round(res_smi_P$L0, 4),
    ", b_SMA =", round(res_smi_P$b_SMA, 4), "\n")
cat("Kallada:   L0 =", round(res_smi_K$L0, 4),
    ", b_SMA =", round(res_smi_K$b_SMA, 4), "\n")

# ================= 10. SUMMARY STATISTICS =================
smi_summary <- smi_data %>%
  group_by(River) %>%
  summarise(
    n          = n(),
    mean_SMI   = mean(SMI, na.rm = TRUE),
    sd_SMI     = sd(SMI, na.rm = TRUE),
    median_SMI = median(SMI, na.rm = TRUE),
    min_SMI    = min(SMI, na.rm = TRUE),
    max_SMI    = max(SMI, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n===== SMI SUMMARY =====\n")
print(smi_summary)

# ================= 11. FORMATTED SUMMARY =================
smi_summary_formatted <- smi_data %>%
  group_by(River) %>%
  summarise(
    n          = n(),
    mean_SMI   = sprintf("%.2f", mean(SMI, na.rm = TRUE)),
    sd_SMI     = sprintf("%.2f", sd(SMI, na.rm = TRUE)),
    median_SMI = sprintf("%.2f", median(SMI, na.rm = TRUE)),
    min_SMI    = sprintf("%.2f", min(SMI, na.rm = TRUE)),
    max_SMI    = sprintf("%.2f", max(SMI, na.rm = TRUE)),
    .groups = "drop"
  )

cat("\n===== FORMATTED SMI SUMMARY =====\n")
print(smi_summary_formatted)

# ================= 12. SHAPIRO-WILK TEST =================
cat("\n===== SHAPIRO-WILK TEST FOR SMI =====\n")
print(by(smi_data$SMI, smi_data$River, shapiro.test))

# ================= 13. KRUSKAL-WALLIS TEST =================
cat("\n===== KRUSKAL-WALLIS TEST FOR SMI =====\n")
kw_smi <- kruskal.test(SMI ~ River, data = smi_data)
print(kw_smi)

# ================= 14. DUNN TEST =================
cat("\n===== DUNN TEST FOR SMI =====\n")
dunn_smi <- dunnTest(SMI ~ River, data = smi_data, method = "bh")
print(dunn_smi)

# ================= 15. BOXPLOT =================
p_smi <- ggplot(smi_data, aes(x = River, y = SMI, fill = River)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.7,
    linewidth = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = "white",
    color = "black",
    stroke = 0.8
  ) +
  coord_cartesian(
    ylim = c(
      quantile(smi_data$SMI, 0.05, na.rm = TRUE),
      quantile(smi_data$SMI, 0.95, na.rm = TRUE)
    )
  ) +
  labs(
    title = "Scaled Mass Index (SMI) by River",
    x = "River",
    y = "SMI"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey75", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3)
  )

print(p_smi)

# ================= 16. EXPORT =================
write.csv(smi_data, "SMI_from_uploaded_csv_natural_log.csv", row.names = FALSE)
write.csv(smi_summary, "SMI_summary_from_uploaded_csv_natural_log.csv", row.names = FALSE)
write.csv(smi_summary_formatted, "SMI_summary_formatted_from_uploaded_csv_natural_log.csv", row.names = FALSE)

cat("\nFiles saved:\n")
cat("- SMI_from_uploaded_csv_natural_log.csv\n")
cat("- SMI_summary_from_uploaded_csv_natural_log.csv\n")
cat("- SMI_summary_formatted_from_uploaded_csv_natural_log.csv\n")
