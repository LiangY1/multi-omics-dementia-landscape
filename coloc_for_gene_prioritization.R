pacman::p_load(data.table, openxlsx, dplyr, tidyr, stringr, coloc, plyr)

window <- 1e6  # 1Mb window for analysis
pp_threshold <- 0.5  # Posterior probability threshold for first filtering

wgs_file_normed <- fread("wgs_file.tsv") %>% as.data.frame()
lead_snvs <- fread("lead_snvs.tsv") %>% as.data.frame()

lead_snvs_data <- wgs_file_normed %>% dplyr::filter(ID %in% lead_snvs$ID)

# Filter lead SNPs
filter_lead_snvs <- function(lead_snvs_data, window) {
  lead_snvs_data <- lead_snvs_data[order(lead_snvs_data$Pvalue),]
  filtered_lead_snvs <- data.frame()
  for (j in 1:nrow(lead_snvs_data)) {
    current_snv <- lead_snvs_data[j,]
    if (nrow(filtered_lead_snvs) == 0) {
      filtered_lead_snvs <- rbind(filtered_lead_snvs, current_snv)
      next
    }
    is_far_enough <- all(abs(current_snv$GENPOS - filtered_lead_snvs$GENPOS) > window | current_snv$CHR != filtered_lead_snvs$CHR)
    if (is_far_enough) filtered_lead_snvs <- rbind(filtered_lead_snvs, current_snv)
  }
  return(filtered_lead_snvs)
}

filtered_lead_snvs <- filter_lead_snvs(lead_snvs_data, window) %>% 
  mutate(identifier = paste(CHROM, GENPOS, sep = ":")) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    region_start = min(wgs_file_normed$POS[wgs_file_normed$CHROM == CHROM & wgs_file_normed$POS >= (POS - (window / 2))], na.rm = TRUE),
    region_end = max(wgs_file_normed$POS[wgs_file_normed$CHROM == CHROM & wgs_file_normed$POS <= (POS + (window / 2))], na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  as.data.frame()
  
snv_list <- filtered_lead_snvs$ID

xQTL_data_normed <- fread("xQTL_file.tsv") %>% 
  as.data.frame() %>% 
  dplyr::select(SNP, CHROM, POS, Allele1, Allele2, BETA, SE, Pvalue, ENSG, SYMBOL)

# N_xQTL = sample size in xQTL studies

# Define function for colocalization analysis
run_coloc_analysis <- function(xQTL_data, wgs_file, snv_list, N_xQTL) {
  coloc_results <- list()
  for (snv in snv_list) { 
    lead_snv <- wgs_file_normed[which(wgs_file_normed$ID == snv), ]
    if (nrow(lead_snv) == 0) next
    lead_snv_position <- lead_snv$POS
    lead_snv_chromosome <- lead_snv$CHR
    
    # Extract xQTL and GWAS data within the window
    wgs_sub <- wgs_file_normed %>% 
      filter(
        CHROM == lead_snv_chromosome, 
        GENPOS >= (lead_snv_position - (window / 2)), 
        GENPOS <= (lead_snv_position + (window / 2)),
        is.finite(BETA), is.finite(SE), !is.na(BETA), !is.na(SE)
      ) %>% 
      mutate(identifier = paste(CHROM, GENPOS, sep = ":"))
    xQTL_sub <- xQTL_data %>% 
      dplyr::filter(
        CHROM == lead_snv_chromosome, 
        POS >= (lead_snv_position - (window / 2)), 
        POS <= (lead_snv_position + (window / 2)),
        is.finite(BETA), is.finite(SE), !is.na(BETA), !is.na(SE)
      ) %>% 
      mutate(identifier = paste(CHROM, POS, sep = ":"))
    
    common_snps <- intersect(wgs_sub$identifier, xQTL_sub$identifier)
    if (length(common_snps) == 0) next
    cat("SNP:", snv, "; Number of common SNPs:", length(common_snps), "\n")
    
    # Filter and align data for common SNPs
    xQTL_snvs_subset_clean <- xQTL_sub %>% 
      filter(identifier %in% common_snps) %>% 
      arrange(Pvalue) %>% 
      distinct(identifier, .keep_all = TRUE) %>% 
      arrange(identifier)
    wgs_snvs_subset_clean <- wgs_sub %>% 
      filter(identifier %in% common_snps) %>% 
      arrange(Pvalue) %>% 
      distinct(identifier, .keep_all = TRUE) %>% 
      arrange(identifier)
    
    # Align alleles between datasets
    keep_rows <- rep(TRUE, nrow(xQTL_snvs_subset_clean))
    for (r in 1:nrow(xQTL_snvs_subset_clean)) {
      if (xQTL_snvs_subset_clean$Allele1[r] == wgs_snvs_subset_clean$Allele1[r] && 
          xQTL_snvs_subset_clean$Allele2[r] == wgs_snvs_subset_clean$Allele2[r]) {
        next
      } else if (xQTL_snvs_subset_clean$Allele1[r] == wgs_snvs_subset_clean$Allele2[r] && 
                 xQTL_snvs_subset_clean$Allele2[r] == wgs_snvs_subset_clean$Allele1[r]) {
        xQTL_snvs_subset_clean$Allele1[r] <- wgs_snvs_subset_clean$Allele1[r]
        xQTL_snvs_subset_clean$Allele2[r] <- wgs_snvs_subset_clean$Allele2[r]
        xQTL_snvs_subset_clean$BETA[r] <- -xQTL_snvs_subset_clean$BETA[r]
      } else {
        keep_rows[r] <- FALSE
      }
    }
    
    xQTL_snvs_subset_clean <- xQTL_snvs_subset_clean[keep_rows, ]
    wgs_snvs_subset_clean <- wgs_snvs_subset_clean[keep_rows, ]
    if (nrow(xQTL_snvs_subset_clean) == 0 || nrow(wgs_snvs_subset_clean) == 0) next
    
    # Perform colocalization analysis
    coloc_result <- coloc.abf(
      dataset1 = list(
        snp = wgs_snvs_subset_clean$ID, 
        beta = wgs_snvs_subset_clean$BETA, 
        varbeta = (wgs_snvs_subset_clean$SE)^2, 
        p = wgs_snvs_subset_clean$Pvalue, 
        N = wgs_snvs_subset_clean$N[1], 
        type = "cc"
      ),
      dataset2 = list(
        snp = wgs_snvs_subset_clean$ID, # This is only valid under the premise that the IDs of WGS data and QTL data are aligned.
        beta = xQTL_snvs_subset_clean$BETA, 
        varbeta = (xQTL_snvs_subset_clean$SE)^2, 
        N = as.numeric(N_xQTL), 
        type = "quant"
      ),
      MAF = wgs_snvs_subset_clean$A1FREQ, 
      p1 = 1e-4, 
      p2 = 1e-4, 
      p12 = 1e-5
    )
    coloc_results[[snv]] <- coloc_result 
  }
  return(coloc_results)
}

# Process colocalization results for each candidate gene
coloc_results_ensg <- list()
final_probelist <- fread("final_probelist.tsv") %>% as.data.frame() %>% pull(ENSG) %>% unique() # candidate gene
for (i in 1:length(final_probelist)) {
  Gene_ensg <- final_probelist[i]
  cat("Processing Gene:", Gene_ensg, "\n")
  
  xQTL_Gene <- xQTL_data_normed %>% dplyr::filter(ENSG == Gene_ensg)
  if (nrow(xQTL_Gene) == 0) next
  coloc_results <- run_coloc_analysis(xQTL_Gene, wgs_file_normed, snv_list, N_xQTL)
  
  if (length(coloc_results) == 0) {
    cat(Gene_ensg, "No common SNPs in xQTL data.\n")
    next
  }
  
  coloc_results_ensg[[Gene_ensg]] <- data.frame()
  for (snp_name in names(coloc_results)) {
    results_df <- coloc_results[[snp_name]]$results
    summary_df <- coloc_results[[snp_name]]$summary
    snp_filtered <- data.frame(
      sig_SNPs = results_df$snp, 
      PPH4_abf = summary_df["PP.H4.abf"], 
      PPH4_SNP = results_df$SNP.PP.H4, 
      row.names = NULL
    ) %>% 
      dplyr::filter(PPH4_abf > pp_threshold) %>% 
      dplyr::slice_max(PPH4_SNP, n = 1)
    if (nrow(snp_filtered) == 0) {
      cat("filtered_lead_SNPs,", snp_name, ", no significant SNPs with PPH4 >", pp_threshold, "\n")
      next
    }
    snp_filtered$filtered_lead_SNPs <- snp_name
    coloc_results_ensg[[Gene_ensg]] <- rbind(coloc_results_ensg[[Gene_ensg]], snp_filtered)
  }
  
  if (nrow(coloc_results_ensg[[Gene_ensg]]) == 0) next 
  
  coloc_results_ensg[[Gene_ensg]] <- coloc_results_ensg[[Gene_ensg]] %>% 
    tidyr::separate(
      sep = ":",
      col = "sig_SNPs",
      into = c("CHROM", "POS"),
      remove = FALSE,
      extra = "drop",
      convert = TRUE
    ) %>% 
    mutate(ENSG = Gene_ensg)
}

coloc_results_PMID <- plyr::rbind.fill(coloc_results_ensg) # then save the combined results
