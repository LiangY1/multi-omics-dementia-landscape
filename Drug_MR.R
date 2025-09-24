#Metabolome/Clinical lab tests-Disease MR#
library(TwoSampleMR)
library(data.table)
F11 <- c("F5_DEMENTIA","F5_VASCDEM","G6_AD_WIDE")
ins <- format_data(IV,type = "exposure",header = T,phenotype_col = "PHENO",snp_col = "ID",beta_col = "BETA",se_col = "SE",eaf_col = "A1_FREQ",effect_allele_col = "A1",other_allele_col = "AX",pval_col = "P")
for (i in 1:length(F11)) {
  pheno <- as.data.frame(fread(paste0("Finngen_GWAS_summary/R11/finngen_R11_",F11[i]),header = T))
  pheno$pheno <- F11[i]
  out <- format_data(pheno,type = "outcome",header = T,snps = ins$SNP,phenotype_col = "pheno",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",effect_allele_col ="alt",other_allele_col = "ref",eaf_col = "af_alt",pval_col = "pval")
  harmo <- harmonise_data(exposure_dat = ins,outcome_dat = out)
  mr_res <- mr(harmo,method_list = c("mr_wald_ratio","mr_ivw"))
  mr_pleio <- mr_pleiotropy_test(harmo)
  mr_hetero <- mr_heterogeneity(harmo)
  results_OR <- generate_odds_ratios(mr_res)
  results_OR <- results_OR[,-c(1,2)]
  write.csv(results_OR,paste0(F11[i],".csv"),row.names = F)
}

#Proteome-Metabolome cis-MR#
library(TwoSampleMR)
library(data.table)
ins <- format_data(cis_IV,type = "exposure",header = T,phenotype_col = "Pro_code",snp_col = "rsid",beta_col = "BETA",se_col = "SE",eaf_col = "A1FREQ",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",pval_col = "P")
for (i in 1:length(list)) {
  pheno <- data.frame()
  for (j in 1:22) {
    dt <- as.data.frame(fread(paste0("UKB_Met_ex_Pro_GWAS/raw/chr",j,".",list[i],".glm.linear"),header = T))
    pheno <- rbind(pheno,dt)
  }
  pheno$PHENO <- list[i]
  out <- format_data(pheno,type = "outcome",header = T,snps = ins$SNP,phenotype_col = "PHENO",snp_col = "ID",beta_col = "BETA",se_col = "SE",effect_allele_col ="A1",other_allele_col = "AX",eaf_col = "A1_FREQ",pval_col = "P")
  harmo <- harmonise_data(exposure_dat = ins,outcome_dat = out)
  mr_res <- mr(harmo,method_list = c("mr_wald_ratio","mr_ivw"))
  mr_pleio <- mr_pleiotropy_test(harmo)
  mr_hetero <- mr_heterogeneity(harmo)
  results_OR <- generate_odds_ratios(mr_res)
  results_OR <- results_OR[,-c(1,2)]
  write.csv(results_OR,paste0(list[i],".csv"),row.names = F)
}

#Proteome-Clinical lab tests cis-MR#
library(TwoSampleMR)
library(data.table)
ins <- format_data(cis_IV,type = "exposure",header = T,phenotype_col = "Pro_code",snp_col = "rsid",beta_col = "BETA",se_col = "SE",eaf_col = "A1FREQ",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",pval_col = "P")
for (i in 1:length(list)) {
  pheno <- data.frame()
  for (j in 1:22) {
    dt <- as.data.frame(fread(paste0("UKB_C_ex_Pro_GWAS/raw/",list[i],"/chr",j,".",list[i],".glm.linear"),header = T))
    pheno <- rbind(pheno,dt)
  }
  pheno$PHENO <- list[i]
  out <- format_data(pheno,type = "outcome",header = T,snps = ins$SNP,phenotype_col = "PHENO",snp_col = "ID",beta_col = "BETA",se_col = "SE",effect_allele_col ="A1",other_allele_col = "AX",eaf_col = "A1_FREQ",pval_col = "P")
  harmo <- harmonise_data(exposure_dat = ins,outcome_dat = out)
  mr_res <- mr(harmo,method_list = c("mr_wald_ratio","mr_ivw"))
  mr_pleio <- mr_pleiotropy_test(harmo)
  mr_hetero <- mr_heterogeneity(harmo)
  results_OR <- generate_odds_ratios(mr_res)
  results_OR <- results_OR[,-c(1,2)]
  write.csv(results_OR,paste0(list[i],".csv"),row.names = F)
}

#Proteome-Disease cis-MR#
library(TwoSampleMR)
library(data.table)
F11 <- c("F5_DEMENTIA","F5_VASCDEM","G6_AD_WIDE")
ins <- format_data(cis_IV,type = "exposure",header = T,phenotype_col = "Pro_code",snp_col = "rsid",beta_col = "BETA",se_col = "SE",eaf_col = "A1FREQ",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",pval_col = "P")
for (i in 1:length(F11)) {
  pheno <- as.data.frame(fread(paste0("Finngen_GWAS_summary/R11/finngen_R11_",F11[i]),header = T))
  pheno$pheno <- F11[i]
  out <- format_data(pheno,type = "outcome",header = T,snps = ins$SNP,phenotype_col = "pheno",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",effect_allele_col ="alt",other_allele_col = "ref",eaf_col = "af_alt",pval_col = "pval")
  harmo <- harmonise_data(exposure_dat = ins,outcome_dat = out)
  mr_res <- mr(harmo,method_list = c("mr_wald_ratio","mr_ivw"))
  mr_pleio <- mr_pleiotropy_test(harmo)
  mr_hetero <- mr_heterogeneity(harmo)
  results_OR <- generate_odds_ratios(mr_res)
  results_OR <- results_OR[,-c(1,2)]
  write.csv(results_OR,paste0(F11[i],".csv"),row.names = F)
}