pacman::p_load(data.table, openxlsx, dplyr, magrittr, survival, plyr)

endpoints <- c("ACD", "AD", "VaD")
omics <- c("Metabolome", "Proteome", "Clinical")
stages <- c("stage1", "stage2")

# Load outcome and covariate data
outcome <- fread("Dementia_outcome_allpeople.csv") %>% 
  as.data.frame()
load("BloodDict")
load("Covar_imputed")
load("dat_sta")

# Extract multiomics IDs
multiomics_id <- dat_sta$eid[dat_sta$Multiomics == 1]

# Define covariates
c1 <- c("age", "sex", "College", "apoe4", "ethn")  # Basic model
c2 <- c("tdi", "smk", "bmi", "diab_hst", "cvd_hst", "med_bp", "med_tc")
Covar_imputed <- Covar_imputed %>% 
  select(all_of(c("eid", c1, c2)))

# Define function for incident analysis
run_incident_analysis <- function(data, variables, time_var, status_var, covariates) {
  data_var <- data
  results <- lapply(variables, function(var) {
    print(paste("Variable:", var))
    if (nrow(data_var) == 0) return(NULL)
    formula_str <- paste0("Surv(", time_var, ",", status_var, " == 1) ~ ", var, " + ", paste(covariates, collapse = " + "))
    formula <- as.formula(formula_str)
    fit <- coxph(formula, data = data_var)
    sfit <- summary(fit)
    beta <- sfit$coefficients[1, "coef"]
    se <- sfit$coefficients[1, "se(coef)"]
    hr <- exp(beta)
    pvalue <- sfit$coefficients[1, "Pr(>|z|)"]
    lci <- exp(beta - 1.96 * se)
    uci <- exp(beta + 1.96 * se)
    data.frame(
      Characteristics = var,
      Outcome = status_var,
      N_case = fit$nevent,
      N = fit$n,
      Pvalue = pvalue,
      Beta = beta,
      Se = se,
      HR = hr,
      LCI = lci,
      UCI = uci
    )
  })
  results_df <- do.call(rbind, results)
  return(results_df)
}

# Loop through endpoints, omics, and stages
for (endpoint in endpoints) {
  for (omic in omics) {
    for (stage in stages) {
      cat("Endpoint:", endpoint, "\n")
      cat("Omic:", omic, "\n")
      cat("Stage:", stage, "\n")
      
      # Load and preprocess omic data
      omic_data <- get(load(paste0(omic, "_imputed"))) %>% 
        dplyr::select(eid, BloodDict$Omics_feature[which(BloodDict$Omics_group == omic)]) %>% 
        mutate_at(vars(BloodDict$Omics_feature[which(BloodDict$Omics_group == omic)]), scale)
      
      # Define outcome variables
      outcome_status <- paste0(endpoint, "_status")
      outcome_years <- paste0(endpoint, "_years")
      
      # Filter ACD baseline cases
      ACD <- outcome %>% 
        filter(ACD_status == 1 & ACD_years < 0)
      
      # Merge data and filter based on stage
      data_merged <- omic_data %>% 
        inner_join(Covar_imputed, by = "eid") %>% 
        inner_join(outcome, by = "eid") %>% 
        filter(!eid %in% ACD$eid)
      
      if (stage == "stage1") {
        data_filtered <- data_merged %>% 
          filter(!eid %in% multiomics_id)
      } else {
        data_filtered <- data_merged %>% 
          filter(eid %in% multiomics_id)
      }
      
      # Prepare data for incident analysis
      controls <- data_filtered %>% 
        filter(.data[[outcome_status]] == 0 & .data[[outcome_years]] > 0)
      incident_cases <- data_filtered %>% 
        filter(.data[[outcome_status]] == 1 & .data[[outcome_years]] > 0)
      incident_data <- bind_rows(incident_cases, controls)
      
      # Define variables and covariates for analysis
      variables <- BloodDict$Omics_feature[which(BloodDict$Omics_group == omic)]
      time_var <- outcome_years
      status_var <- outcome_status
      covariates <- c1
      
      # Run incident analysis and save results
      incident_results <- run_incident_analysis(incident_data, variables, time_var, status_var, covariates)
      fwrite(
        incident_results,
        file = paste0(stage, "_", endpoint, "_", omic, "_incident.csv"),
        row.names = FALSE,
        quote = FALSE
      )
    }
  }
}