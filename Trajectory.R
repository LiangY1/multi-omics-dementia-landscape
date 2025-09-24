pacman::p_load(data.table,dplyr,plyr,magrittr,MatchIt,cobalt,stringr)
rm(list=ls())

endpoints <- c("ACD", "AD", "VaD")
omics <- c("Metabolome", "Proteome", "Clinical")
span <- 0.8  # Span values for LOESS smoothing

# Load outcome and covariate data
outcome <- fread("Dementia_outcome_allpeople.csv") %>% 
  as.data.frame()
load("BloodDict")
load("Covar_imputed")
load("dat_sta")

dat_sta <- dat_sta %>% 
  dplyr::mutate(Trajectory = ifelse(Metabolome + Proteome + Clinical == 3, 1, 0))
Trajectory_id <- dat_sta$eid[dat_sta$Trajectory == 1]

# Filter ACD baseline cases
ACD <- outcome %>% 
  filter(ACD_status == 1 & ACD_years < 0)

Neurod <- fread("NeuroD_processed.csv") %>% as.data.frame() %>% 
  filter(ONeuroD_y == 1 & follow_up < 0)  # Neurod baseline

# Define covariates
c1 <- c("age", "sex", "College", "apoe4", "ethn")  # Basic model
c2 <- c("tdi", "smk", "bmi", "diab_hst", "cvd_hst", "med_bp", "med_tc")
Covar_imputed <- Covar_imputed %>% 
  dplyr::select(all_of(c("eid", c1, c2)))

for (endpoint in endpoints) {
  for (omic in omics) {
    cat("Endpoint:", endpoint, "\n")
    cat("Omic:", omic, "\n")
    cat("Span:", span, "\n")
    
    omic_data <- get(load(paste0(omic, "_imputed"))) %>% 
      dplyr::select(eid, BloodDict$Omics_feature[which(BloodDict$Omics_group == omic)]) 
    data_merged <- omic_data %>% 
      inner_join(Covar_imputed, by = "eid") %>% 
      inner_join(outcome, by = "eid") %>% 
      filter(eid %in% Trajectory_id)
    
    # Define outcome variables
    outcome_status <- paste0(endpoint, "_status")
    outcome_years <- paste0(endpoint, "_years")
    
    # Define case and control groups
    case <- data_merged %>% 
      filter(.data[[outcome_status]] == 1)
    
    control <- data_merged %>% 
      filter(!eid %in% c(ACD$eid, Neurod$eid), .data[[outcome_status]] == 0)
    match_final <- rbind(case, control)
    
    # Perform matching
    mt_out1 <- matchit(
      as.formula(paste(outcome_status, "~ age + College + ethn + bmi")),
      method = "nearest",
      distance = "mahalanobis",
      ratio = 10,
      link = "logit",
      exact = ~sex,
      data = match_final
    )
    mt_data <- match.data(mt_out1)
    
    # Assign same time to case and control within matched groups
    mt_final <- mt_data %>% 
      arrange(subclass)
    time_scale <- mt_final[[outcome_years]][which(mt_final[[outcome_status]] == 1)]
    mt_final <- mt_final %>% 
      mutate(time_scale = rep(time_scale, each = 11) * sign(-1))
    
    # Residuals
    resid_mt <- mt_final %>% magrittr::set_rownames(.$eid)
    for (i in names(omic_data)[2:ncol(omic_data)]) {
      print(i)
      FML <- paste0(i, " ~ ", paste(c(c1, c2), collapse = " + "))
      sub_data <- resid_mt[, c("eid", i, c1, c2)] %>% 
        na.omit() %>% 
        magrittr::set_rownames(.$eid)
      sub_data[, i] <- resid(lm(FML, data = sub_data))
      resid_mt[, i] <- NA
      resid_mt[rownames(sub_data), i] <- as.numeric(scale(sub_data[, i]))
    }
    save(resid_mt, file = paste0("AC_resid_", endpoint, "_", omic))
    
    # Calculate Z-scores    
    df_group <- resid_mt %>% 
      ungroup()
    df_group_case <- df_group %>% 
      filter(.data[[outcome_status]] == 1)
    df_group_control <- df_group %>% 
      filter(.data[[outcome_status]] == 0)
    
    stats_control <- df_group_control %>%
      summarise_at(
        vars(all_of(variable.names)),
        list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))
      ) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = c("Omics_feature", ".value"),
        names_pattern = "(.*)_(mean|sd)"
      ) %>% 
      as.data.frame() %>% 
      set_rownames(.$Omics_feature)
    
    df_z_scores <- df_group_case
    for (var_name in variable.names) {
      df_z_scores[[var_name]] <- (df_group_case[[var_name]] - stats_control[var_name, "mean"]) / 
        stats_control[var_name, "sd"]
    }
    t <- df_z_scores %>% 
      dplyr::select(eid, subclass, Year = time_scale, all_of(variable.names))
    save(t, file = paste0("AC_trajectories_Zscore_", endpoint, "_", omic))
    
    # LOESS
    df_t <- t %>% 
      tidyr::pivot_longer(cols = variable.names, names_to = "Characteristics", values_to = "Estimate")
    df_loess <- do.call(rbind, lapply(variable.names, function(var_name) {
      df2 <- df_t %>% 
        filter(Characteristics == var_name) %>% 
        na.omit()
      loess_fit <- loess(Estimate ~ Year, data = df2, span = span)
      gap <- seq(-15, 0, 0.5)
      data.frame(Year = gap, Characteristics = rep(var_name, length(gap))) %>% 
        dplyr::mutate(Estimate_loess = predict(loess_fit, .))
    }))
    save(df_loess, file = paste0("AC_trajectories_loess", span, "_", endpoint, "_", omic))
  }
}

pacman::p_load(Mfuzz, corrplot, cowplot, ggplot2)
omics_color <- c(
  Metabolome = "#BC80BD",
  Proteome = "#8DD3C7",
  Clinical = "#FDB462"
)

endpoints <- c("ACD", "AD", "VaD")
omics <- c("Metabolome", "Proteome", "Clinical")
span <- 0.8  # Span value for LOESS smoothing

load("BloodDict")

for (endpoint in endpoints) {
  
  sig_features_combined <- fread(paste0("sig_features_combined_", endpoint, ".csv")) %>% 
    as.data.frame() %>% pull(Omics_feature)
  
  loess_omics <- list()
  for (omic in omics) {
    load(paste0("AC_trajectories_loess", span, "_", endpoint, "_", omic))
    df_loess <- df_loess %>% 
      filter(Characteristics %in% BloodDict$Omics_feature[which(BloodDict$Omics_group == omic)])
    loess_omics[[omic]] <- df_loess
  }
  loess_omics <- plyr::rbind.fill(loess_omics)
  
  df_wide <- loess_omics %>% 
    filter(Characteristics %in% sig_features_combined) %>% 
    tidyr::pivot_wider(names_from = Characteristics, values_from = Estimate_loess) %>% 
    tibble::column_to_rownames("Year") %>% 
    t() %>% 
    as.matrix()
  mfuzz_class <- new("ExpressionSet", exprs = df_wide)
  m1 <- mestimate(mfuzz_class)
  
  # Select optimal number of clusters and export plot data
  set.seed(123)
  dir.create(endpoint)
  pdf(paste0(endpoint, "/", endpoint, ".pdf"), width = 5, height = 5)
  min_distances <- Dmin(mfuzz_class, m = m1, crange = 2:15, repeats = 3, visu = TRUE)
  dev.off()
  
  # Perform fuzzy clustering with optimal number of clusters
  optimal_clusters <- 4  # ACD:6, AD:4, VaD:4
  fcm_final <- mfuzz(mfuzz_class, c = optimal_clusters, m = m1)
  center <- get_mfuzz_center(data = mfuzz_class, c = fcm_final, membership_cutoff = membership_cutoff)
  rownames(center) <- paste("Cluster", rownames(center), sep = " ")
  cluster_info <- data.frame(
    Omics_feature = names(fcm_final$cluster),
    fcm_final$membership,
    cluster = fcm_final$cluster,
    stringsAsFactors = FALSE
  ) %>% 
    arrange(cluster)
  
  idx <- 1
  cluster_specific <- cluster_center <- plot_idx <- list()
  for (idx in 1:optimal_clusters) {
    cat("idx: ", idx, "\n")
    
    cluster_data <- cluster_info %>% 
      dplyr::filter(cluster == idx) %>% 
      dplyr::select(1, 1 + idx, cluster)
    colnames(cluster_data)[2] <- "membership"
    cluster_data <- cluster_data %>% 
      dplyr::filter(membership > membership_cutoff)
    
    cluster_center[[idx]] <- center[idx, , drop = TRUE] %>% 
      unlist() %>%
      data.frame(
        cluster = idx,
        time = names(.),
        value = .,
        stringsAsFactors = FALSE
      ) %>%
      dplyr::mutate(time = as.numeric(time))
    
    cluster_specific[[idx]] <- mfuzz_class_scaled[cluster_data$Omics_feature, ] %>%
      data.frame(
        cluster = cluster_data$cluster,
        membership = cluster_data$membership,
        .,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ) %>%
      tibble::rownames_to_column(var = "Omics_feature") %>%
      tidyr::pivot_longer(
        cols = -c(cluster, Omics_feature, membership),
        names_to = "time",
        values_to = "value"
      ) %>%
      dplyr::mutate(time = as.numeric(time)) %>%
      dplyr::left_join(BloodDict[, c("Omics_feature", "Omics_group")], by = "Omics_feature") %>%
      dplyr::arrange(membership, Omics_feature) %>%
      dplyr::arrange(desc(Omics_group)) %>%
      dplyr::mutate(Omics_feature = factor(Omics_feature, levels = unique(Omics_feature)))
    
    plot_idx[[idx]] <- cluster_specific[[idx]] %>%
      ggplot(aes(time, value, group = Omics_feature)) +
      geom_line(aes(color = Omics_group), alpha = 0.7) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)
      ) +
      labs(
        x = "",
        y = "Z-score",
        title = paste("Cluster ", idx, " (", nrow(cluster_data), " molecules)", sep = "")
      ) +
      geom_line(mapping = aes(time, value, group = 1), data = cluster_center[[idx]], size = 2) +
      geom_hline(yintercept = 0) +
      scale_color_manual(values = omics_color)
  }
  
  plot_cluster <- do.call(plot_grid, c(plot_idx, align = "hv", ncol = 3))
  
  single_width <- 3
  single_height <- 3
  ncol <- 3
  nrow <- ceiling(optimal_clusters / ncol)
  total_width <- single_width * ncol
  total_height <- single_height * nrow
  ggsave(
    filename = paste0(endpoint, "/Cluster_", optimal_clusters, ".pdf"),
    plot = plot_cluster,
    width = total_width,
    height = total_height,
    units = "in"
  )
  
  cluster_info <- unique(cluster_info$cluster) %>% 
    purrr::map(function(x) {
      temp <- cluster_info %>% 
        dplyr::select(Omics_feature, paste0("X", x), cluster)
      colnames(temp)[2] <- "membership"
      temp <- temp %>% 
        dplyr::filter(membership >= membership_cutoff)
      temp
    }) %>% 
    dplyr::bind_rows() %>% 
    as.data.frame()
  cluster_info %>% 
    dplyr::count(cluster)
  save(
    cluster_info,
    cluster_specific,
    cluster_center,
    file = paste0(endpoint, "/Cluster_", optimal_clusters)
  )
}


pacman::p_load(DEswan)

endpoints <- c("ACD", "AD", "VaD")
omics <- c("Metabolome", "Proteome", "Clinical")
span <- 0.8  # Span value for LOESS smoothing

deswan_c <- c("age", "sex", "ethn")  # Covariates for DEswan

outcome <- fread("Dementia_outcome_allpeople.csv") %>% 
  as.data.frame()
load("BloodDict")
load("Covar_imputed")
load("dat_sta")

dat_sta <- dat_sta %>% 
  dplyr::mutate(Trajectory = ifelse(Metabolome + Proteome + Clinical == 3, 1, 0))
Trajectory_id <- dat_sta$eid[dat_sta$Trajectory == 1]

Covar_imputed <- Covar_imputed %>% select(all_of(c("eid", deswan_c)))

remove_outliers <- function(x) {
  median_x <- median(x, na.rm = TRUE)
  iqr_x <- IQR(x, na.rm = TRUE)
  x[which(abs(x - median_x) > 4 * iqr_x)] <- NA
  return(x)
}

for (endpoint in endpoints) {
  cat("Endpoint:", endpoint, "\n")
  outcome_status <- paste0(endpoint, "_status")
  outcome_years <- paste0(endpoint, "_years")
  
  # Load original omics data
  ukb_met <- fread("Metabolomics.tsv.gz") %>% 
    as.data.frame() %>% 
    filter(eid %in% Trajectory_id) %>% 
    dplyr::select(eid, BloodDict$Omics_feature[which(BloodDict$Omics_group == "Metabolome")]) %>% 
    mutate_at(vars(BloodDict$Omics_feature[which(BloodDict$Omics_group == "Metabolome")]), remove_outliers)
  
  ukb_pro <- fread("Proteomics.tsv.gz") %>% 
    as.data.frame() %>% 
    filter(eid %in% Trajectory_id) %>% 
    dplyr::select(eid, BloodDict$Omics_feature[which(BloodDict$Omics_group == "Proteome")]) %>% 
    mutate_at(vars(BloodDict$Omics_feature[which(BloodDict$Omics_group == "Proteome")]), remove_outliers)
  
  ukb_clin <- fread("Clinical.tsv.gz") %>% 
    as.data.frame() %>% 
    filter(eid %in% Trajectory_id) %>% 
    dplyr::select(eid, BloodDict$Omics_feature[which(BloodDict$Omics_group == "Clinical")]) %>% 
    mutate_at(vars(BloodDict$Omics_feature[which(BloodDict$Omics_group == "Clinical")]), remove_outliers)

  omics_dat <- Reduce(function(x, y) merge(x, y, by = "eid"), list(ukb_met, ukb_clin, ukb_pro)) %>% 
    as.data.frame() %>%
    inner_join(Covar_imputed, by = "eid") %>% 
    inner_join(outcome, by = "eid") %>% 
    filter(.data[[outcome_status]] == 1) %>%
    dplyr::select(eid, Year = .data[[outcome_years]], BloodDict$Omics_feature, deswan_c) %>% 
    dplyr::mutate(Year = -Year)

  res_DEswan <- res_p <- list()
  for (parcel_width in 4:9) {
    character_parcel_width <- as.character(parcel_width)
    res_DEswan[[character_parcel_width]] <- DEswan(
      data.df = omics_dat[, BloodDict$Omics_feature],
      qt = omics_dat[, "Year"],
      covariates = omics_dat[, deswan_c],
      window.center = seq(-13, -3, 1),
      buckets.size = parcel_width
    )
    res_p[[character_parcel_width]] <- res_DEswan[[character_parcel_width]] %>% 
      reshape.DEswan(parameter = 1, factor = "qt")
  }
  save(res_DEswan, res_p, file = paste0(endpoint))
}