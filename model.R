####################################
# Coagulation-flocculation model   #
# Jose A. D. Lopez-Coronado        #
# Agriculture Victoria - 2026      #
####################################

# Libraries
library(dplyr)
library(lme4)

# Dataset
ch3 <- read.csv(
  "https://raw.githubusercontent.com/jlopezco-wm/tanflocpfsmodel/main/Lopez_coagulation_JWREUSE.csv",
  stringsAsFactors = TRUE
)

# -------------------------------
# FIT MODELS
# -------------------------------

tanfloc_rem <- lmer(
  removal ~ dose * turb_ini + I(dose^2) + I(turb_ini^2) + (1 | round),
  data = subset(ch3, coagulant == "tanfloc")
)

pfs_rem <- lmer(
  removal ~ dose * turb_ini + I(dose^2) + I(turb_ini^2) + (1 | round),
  data = subset(ch3, coagulant == "pfs")
)

# Tag models
attr(tanfloc_rem, "coagulant") <- "tanfloc"
attr(pfs_rem, "coagulant") <- "pfs"

# Save models
saveRDS(tanfloc_rem, "tanfloc_rem.rds")
saveRDS(pfs_rem, "pfs_rem.rds")

# -------------------------------
# DOSE FUNCTION
# -------------------------------

dose_mod <- function(model,
                     turbidity,
                     target = 0,
                     max_dose = 10000,
                     effluent,
                     concentration,
                     nsim = 1000) {
  
  # ---- MODEL CHECK ----
  if (!inherits(model, "lmerMod")) {
    stop("model must be a fitted lmer model")
  }
  
  # ---- INITIAL CONDITIONS ----
  turb_ini <- turbidity * (1 - 0.681)
  doses <- 0:max_dose
  
  newdat <- data.frame(
    dose = doses,
    turb_ini = turb_ini
  )
  
  # ---- HELPER FUNCTION ----
  get_dose <- function(preds) {
    if (length(preds) < 2 || all(is.na(preds))) return(NA_real_)
    
    crossing <- target >= turb_ini - preds
    idx <- which(crossing)[1]
    
    if (is.na(idx)) return(NA_real_)
    
    idx - 1
  }
  
  # ---- POINT ESTIMATE ----
  preds <- predict(model, newdata = newdat, re.form = NA)
  dose_hat <- get_dose(preds)
  
  # ---- BOOTSTRAP ----
  sim_fun <- function(fit) {
    preds_sim <- predict(fit, newdata = newdat, re.form = NA)
    get_dose(preds_sim)
  }
  
  sims <- lme4::bootMer(
    model,
    FUN = sim_fun,
    nsim = nsim,
    use.u = FALSE,
    type = "parametric"
  )
  
  dose_sims <- sims$t[,1]
  dose_sims <- dose_sims[!is.na(dose_sims)]
  
  # ---- CHECK ----
  if (is.na(dose_hat) || length(dose_sims) < 10) {
    return("No reliable dose estimate")
  }
  
  # ---- CONVERT TO mL ----
  dose_ml_hat <- dose_hat * effluent / concentration
  dose_ml <- dose_sims * effluent / concentration
  
  # ---- CI ----
  dose_ml_ci <- quantile(
    dose_ml,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  # ---- OUTPUT ----
  paste(
    paste("Dose:", round(dose_hat, 1), "mg/L"),
    paste("Volume:", round(dose_ml_hat, 1), "mL"),
    paste(
      "Acceptable range:",
      round(dose_ml_ci[1], 1), "-",
      round(dose_ml_ci[3], 1), "mL"
    ),
    sep = "\n"
  )
}

# Save function
saveRDS(dose_mod, "dose_mod.rds")

# -------------------------------
# TESTS
# -------------------------------

# Tanfloc (auto uses 207)
dose_mod(
  model = tanfloc_rem,
  turbidity = 2856.5,
  effluent = 80
)

# PFS (auto uses 100)
dose_mod(
  model = pfs_rem,
  turbidity = 2856.5,
  effluent = 80
)

# Compare model coefficients
fixef(tanfloc_rem)
fixef(pfs_rem)
