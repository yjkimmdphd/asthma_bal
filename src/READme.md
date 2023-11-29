nasal_BAL_rnaseq_v2.R sources phenotype_cleanup_nasal_BAL_rnaseq.R to get phenotype information
nasal_BAL_rnaseq_v3.R sources phenotype_cleanup_nasal_BAL_rnaseq_v2.R to get phenotype information
phenotype_cleanup_nasal_BAL_rnaseq.R sets gene cutoff as lcpm = 0 
phenotype_cleanup_nasal_BAL_rnaseq_v2.R sets gene cutoff as lcpm = log2(10/M+2/L), # M is median. L is mean. library size. M = 48.24072. L = 51.48902
