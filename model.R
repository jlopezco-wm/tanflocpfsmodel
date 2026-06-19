#Library
library(dplyr)
library(lme4)

#Dataset

ch3 <- read.csv(
  "https://raw.githubusercontent.com/jlopezco-wm/coagulation/main/Lopez_coagulation_JWREUSE.csv",
  stringsAsFactors = TRUE
)

#Model tanfloc
tanfloc_rem<-lmer(removal ~ dose * turb_ini + I(dose^2) + I(turb_ini^2) + (1 | round), data = subset(ch3,coagulant == "tanfloc"))
saveRDS(tanfloc_rem, "tanfloc_rem.rds") #this RDS contains the model formula with coefficients, no data is stored

#Model PFS
pfs_rem<-lmer(removal ~ dose * turb_ini + I(dose^2) + I(turb_ini^2) + (1 | round), data = subset(ch3,coagulant == "pfs"))
saveRDS(pfs_rem, "pfs_rem.rds") #this RDS contains the model formula with coefficients, no data is stored


# predict(tanfloc_rem, newdata = data.frame(dose = 1337, turb_ini = 2856.5*(1 - 0.681)), re.form = NA)   test

# Max clarification model -------------------------------------------------

#turbidity is in NTU
#effluent in in litres
#concentration is in mg/L

# Function  1------------------------------------------------------------
dose_mod <- function(model,
                     turbidity,
                     target = 0 ,
                     max_dose = 10000,
                     effluent,
                     concentration,
                     nsim = 1000) {
  
  # ---- MODEL CHECK ----
  if (!inherits(model, "lmerMod")) {
    stop("model must be a fitted lmer model")
  }
  
  # ---- DEFAULT CONCENTRATION LOGIC ----
  if (missing(concentration)) {
    
    model_name <- deparse(substitute(model))
    
    if (model_name == "tanfloc_rem") {
      concentration <- 207
    } else if (model_name == "pfs_rem") {
      concentration <- 100
    } else {
      stop("Model not recognised for default concentration")
    }
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
  
  # ---- BOOTSTRAP SIMULATION ----
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
  
  dose_sims <- sims$t[, 1]
  dose_sims <- dose_sims[!is.na(dose_sims)]
  
  # ---- CHECK RESULTS ----
  if (is.na(dose_hat) || length(dose_sims) < 10) {
    return("No reliable dose estimate")
  }
  
  # ---- CONVERT TO mL ----
  dose_ml_hat <- dose_hat * effluent / concentration
  dose_ml <- dose_sims * effluent / concentration
  
  # ---- CONFIDENCE INTERVAL ----
  dose_ml_ci <- quantile(
    dose_ml,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  # ---- FORMATTED OUTPUT ----
  paste(
    paste("Dose:", round(dose_hat, 1), "mg/L"),
    paste("Volume:",
          round(dose_ml_hat, 1),
          "mL"),
    paste(
      "Acceptable range:",
      round(dose_ml_ci[1], 1), "-",
      round(dose_ml_ci[3], 1),
      "mL"
    ),
    sep = "\n"
  )
}
  
  # ---- rest of your function ----
#dose_mod(model = pfs_rem,turbidity = 2856.5,concentration = 207, effluent = 80)
saveRDS(dose_mod, "dose_mod.rds")

