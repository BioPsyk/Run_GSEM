#!/usr/bin/env Rscript

require(GenomicSEM)
require(dplyr)
require(data.table)

args           = commandArgs(trailingOnly = TRUE)
sum_stats      = read.table(args[1])
hm3_snps_path  = args[2]
ld_folder      = args[3]
reference      = args[4]
info_threshold = args[5]
maf_threshold  = args[6]
out_prefix     = args[7]
chromosome     = args[8]

colnames(sum_stats) = c("TRAIT", 
                        "FILE_PATH", 
                        "N_EFF", 
                        "SAMPLE_PREV", 
                        "POP_PREV")

# ----- STEP 1: MUNGE SUMMARY STATS ----- 

munge(files = sum_stats$FILE_PATH, 
      hm3 = hm3_snps_path,
      trait.names = sum_stats$TRAIT, 
      N = sum_stats$N_EFF, 
      info.filter = info_threshold, 
      maf.filter = maf_threshold)


# ----- STEP 2: RUN MULTI-VARIABLE LDSC -----

sum_stats       = sum_stats %>% 
    mutate(MUNGED_FILE = paste0(TRAIT, ".sumstats.gz"))
LDSCoutput      = ldsc(traits = sum_stats$MUNGED_FILE, 
                       sample.prev = sum_stats$SAMPLE_PREV, 
                       population.prev = sum_stats$POP_PREV, 
                       ld = ld_folder, 
                       wld = ld_folder, 
                       trait.names = sum_stats$TRAIT,
                       ldsc.log = paste0(out_prefix, ".log"),
                       stand = TRUE)

# ----- STEP 3: PREPARE SUMMARY STATS FOR GWAS -----

for ( i in 1:length(sum_stats)) {
    file = fread(sum_stats$FILE_PATH[i])
    file_chr = file %>% filter(CHR == chromosome)
    write.table(file_chr, 
                paste0(sum_stats$TRAIT, "_CHR", chromosome, ".split.assoc"),
                row.names = F,
                quote = F,
                sep = " ")
}

sum_stats = sum_stats %>% 
    mutate(CHR_FILE = paste0(TRAIT, "_CHR", chromosome, ".split.assoc"))

p_sumstats <-sumstats(files = sum_stats$CHR_FILE,
                      ref = reference,
                      trait.names = sum_stats$TRAIT,
                      se.logit = rep(T, length(sum_stats)),
                      OLS = NULL,
                      linprob = NULL,
                      N = sum_stats$N_EFF,
                      betas = NULL,
                      info.filter = info_threshold,
                      maf.filter = maf_threshold,
                      keep.indel = FALSE,
                      parallel = FALSE,
                      cores=NULL)

# ----- STEP 4: RUN COMMON FACTOR GWAS, WRITE OUTPUT -----

pfactor = commonfactorGWAS(covstruc = LDSCoutput, SNPs = p_sumstats)

write.table(pfactor, 
            paste0(out_prefix, "_", "COMMON_FACTOR_GWAS_CHR", chromosome, ".txt"), 
            row.names = F, sep = " ", quote = F)