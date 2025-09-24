# Activate conda environment for regenie
source conda.sh
conda activate regenie_env

# Single-variant Analysis ======================================================

# Logistic Regression - Step 1
# Description:  (m0: Stage 1, m1： Stage 2)
for part in m0 m1; do
    for phenotype in ACD AD VaD; do
        regenie \
            --step 1 \
            --bed ukb_cal_allChrs_hg38 \
            --extract qc_pass.snplist \
            --phenoFile ${part}_${phenotype}_pheno_step1 \
            --covarFile ${part}_${phenotype}_cov_step1 \
            --catCovarList batch,X54,sex \
            --maxCatLevels 27 \
            --bt \
            --bsize 1000 \
            --threads 64 \
            --lowmem \
            --lowmem-prefix ${part}_${phenotype}_tmp_preds \
            --out ${part}_${phenotype}_step1
    done
done

# Logistic Regression - Step 2
for part in m0 m1; do
    for phenotype in ACD AD VaD; do
        for chr in {1..22}; do
            regenie \
                --step 2 \
                --chr ${chr} \
                --bed Q0_unre_Caucasian_c${chr} \
                --phenoFile ${part}_${phenotype}_pheno_step2 \
                --covarFile ${part}_${phenotype}_cov_step2 \
                --catCovarList batch,X54,sex \
                --maxCatLevels 27 \
                --pred ${part}_${phenotype}_step1_pred.list \
                --bt \
                --firth --approx \
                --firth-se \
                --pThresh 0.05 \
                --minMAC 5 \
                --bsize 1000 \
                --threads 128 \
                --write-samples \
                --print-pheno \
                --out ${part}_${phenotype}_chr${chr}
        done
    done
done

# Survival Analysis - Step 1
for part in m0 m1; do
    for phenotype in ACD AD VaD; do
        regenie \
            --step 1 \
            --t2e \
            --bed ukb_cal_allChrs_hg38 \
            --extract qc_pass.snplist \
            --phenoFile ${part}_${phenotype}_pheno_survival \
            --phenoColList ${phenotype}_age \
            --eventColList ${phenotype}_status \
            --covarFile ${part}_${phenotype}_cov_survival \
            --catCovarList batch,X54,sex \
            --maxCatLevels 27 \
            --bt \
            --bsize 1000 \
            --threads 100 \
            --lowmem \
            --lowmem-prefix surv_${part}_${phenotype}_tmp_preds \
            --out surv_${part}_${phenotype}
    done
done

# Survival Analysis - Step 2
for part in m0 m1; do
    for phenotype in ACD AD VaD; do
        for chr in {1..22}; do
            regenie \
                --step 2 \
                --t2e \
                --chr ${chr} \
                --bed Q0_unre_Caucasian_c${chr} \
                --phenoFile ${part}_${phenotype}_pheno_survival \
                --phenoColList ${phenotype}_age \
                --eventColList ${phenotype}_status \
                --covarFile ${part}_${phenotype}_cov_survival \
                --catCovarList batch,X54,sex \
                --maxCatLevels 27 \
                --pred surv_${part}_${phenotype}_pred.list \
                --bt \
                --firth --approx \
                --firth-se \
                --pThresh 0.05 \
                --minMAC 5 \
                --bsize 1000 \
                --threads 80 \
                --write-samples \
                --print-pheno \
                --out ${part}_${phenotype}_chr${chr}
        done
    done
done

# Gene-level Analysis =========================================================

# Logistic Regression
for part in m0 m1; do
    for phenotype in ACD AD VaD; do
        for genebased in PTV Missense Splice UTR_3 UTR_5; do
            regenie \
                --step 2 \
                --chr ${chr} \
                --bed Q0_unre_Caucasian_c${chr} \
                --phenoFile ${part}_${phenotype}_pheno_step2 \
                --covarFile ${part}_${phenotype}_cov_step2 \
                --catCovarList batch,X54,sex \
                --maxCatLevels 27 \
                --pred ${part}_${phenotype}_step1_pred.list \
                --anno-file ${genebased}_chr${chr}.txt \
                --set-list chr${chr}_${genebased}.setlist \
                --mask-def Mask_${genebased}.txt \
                --aaf-bins 0.01,0.001,0.0001 \
                --vc-tests skato,acato-full \
                --vc-maxAAF 0.01 \
                --joint minp,acat,sbat \
                --rgc-gene-p \
                --check-burden-files \
                --write-mask \
                --bt \
                --firth --approx \
                --firth-se \
                --pThresh 0.05 \
                --bsize 200 \
                --threads 100 \
                --write-samples \
                --write-mask-snplist \
                --print-pheno \
                --out ${part}_${phenotype}_chr${chr}_${genebased}
        done
    done
done

# Survival Analysis
for part in m0 m1; do
    for phenotype in ACD AD VaD; do
        for genebased in PTV Missense Splice UTR_3 UTR_5; do
            regenie \
                --step 2 \
                --t2e \
                --chr ${chr} \
                --bed Q0_unre_Caucasian_c${chr} \
                --phenoFile ${part}_${phenotype}_pheno_survival \
                --phenoColList ${phenotype}_age \
                --eventColList ${phenotype}_status \
                --covarFile ${part}_${phenotype}_cov_survival \
                --catCovarList batch,X54,sex \
                --maxCatLevels 27 \
                --pred surv_${part}_${phenotype}_pred.list \
                --anno-file ${genebased}_chr${chr}.txt \
                --set-list chr${chr}_${genebased}.setlist \
                --mask-def Mask_${genebased}.txt \
                --aaf-bins 0.01,0.001,0.0001 \
                --vc-maxAAF 0.01 \
                --check-burden-files \
                --write-mask \
                --bt \
                --firth --approx \
                --firth-se \
                --pThresh 0.05 \
                --bsize 200 \
                --threads 100 \
                --write-samples \
                --write-mask-snplist \
                --print-pheno \
                --out ${part}_${phenotype}_chr${chr}_${genebased}
        done
    done
done
