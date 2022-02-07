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
                        "N", 
                        "SAMPLE_PREV", 
                        "POP_PREV")
LDSCoutput = readRDS(ldsc_out)

# ----- STEP 4: PREPARE SUMMARY STATS FOR GWAS -----

sum_stats = sum_stats %>% 
    mutate(CHR_FILE = paste0(TRAIT, "_chr", chromosome, ".split.assoc"))

p_sumstats = sumstats(files = sum_stats$CHR_FILE,
                      ref = reference,
                      trait.names = sum_stats$TRAIT,
                      se.logit = rep(T, nrow(sum_stats)),
                      OLS = NULL,
                      linprob = NULL,
                      N = sum_stats$N,
                      betas = NULL,
                      info.filter = info_threshold,
                      maf.filter = maf_threshold,
                      keep.indel = FALSE,
                      parallel = FALSE,
                      cores=NULL)

# ----- STEP 5: RUN COMMON FACTOR GWAS, WRITE OUTPUT -----

model_2fm = "NA * ADHD + BIP + MDD + NonSczPsyc
F2 =~ NA * ANO + AUT + SCZ + SczAFF"

pfactor = userGWAS(covstruc = LDSCoutput, 
                   model = model_2fm,
                   SNPs = p_sumstats,
                   smooth_check = TRUE)

write.table(pfactor, 
            paste0(out_prefix, "_", "UserGWAS_2FM_chr", chromosome, ".txt"), 
            row.names = F, 
            sep = " ", 
            quote = F)