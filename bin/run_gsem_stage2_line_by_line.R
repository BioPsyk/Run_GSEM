#!/usr/bin/env Rscript

require(GenomicSEM)
require(dplyr)

args           = commandArgs(trailingOnly = TRUE)
sum_stats      = read.table(args[1], header = T)
ldsc_out       = args[2]
reference      = args[3]
info_threshold = args[4]
maf_threshold  = args[5]
out_prefix     = args[6]
chromosome     = args[7]

colnames(sum_stats) = c("TRAIT", 
                        "FILE_PATH", 
                        "N_EFF", 
                        "SAMPLE_PREV", 
                        "POP_PREV")
LDSCoutput = readRDS(ldsc_out)

# ----- STEP 4: PREPARE SUMMARY STATS FOR GWAS -----

sum_stats = sum_stats %>% 
    mutate(CHR_FILE = paste0(TRAIT, "_CHR", chromosome, ".split.assoc"))

p_sumstats = sumstats(files = sum_stats$CHR_FILE,
                      ref = reference,
                      trait.names = sum_stats$TRAIT,
                      se.logit = rep(T, nrow(sum_stats)),
                      OLS = NULL,
                      linprob = NULL,
                      N = sum_stats$N_EFF,
                      betas = NULL,
                      info.filter = info_threshold,
                      maf.filter = maf_threshold,
                      keep.indel = FALSE,
                      parallel = FALSE,
                      cores=NULL)

# ----- STEP 5: RUN COMMON FACTOR GWAS, WRITE OUTPUT LINE BY LINE -----

sum_stats_chunk = p_sumstats[1, ]

pfactor = commonfactorGWAS(covstruc = LDSCoutput, SNPs = sum_stats_chunk)
write.table(pfactor, 
            paste0(out_prefix, 
                   "_", 
                   "COMMON_FACTOR_GWAS_CHR", 
                   chromosome, 
                   ".txt"), 
            row.names = FALSE, 
            sep = " ", 
            quote = FALSE,
            append = FALSE)

for(i in 2:nrow(p_sumstats)) {
    sum_stats_chunk = p_sumstats[i,]
    pfactor = commonfactorGWAS(covstruc = LDSCoutput, SNPs = sum_stats_chunk)

    write.table(pfactor, 
                paste0(out_prefix, 
                       "_", 
                       "COMMON_FACTOR_GWAS_CHR", 
                       chromosome, 
                       ".txt"), 
                row.names = FALSE, 
                col.names = FALSE,
                sep = " ", 
                quote = FALSE,
                append = TRUE)
}