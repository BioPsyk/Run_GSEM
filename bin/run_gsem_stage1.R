#!/usr/bin/env Rscript

require(GenomicSEM)
require(dplyr)
require(data.table)

args           = commandArgs(trailingOnly = TRUE)
sum_stats      = read.table(args[1], header = T)
hm3_snps_path  = args[2]
ld_folder      = args[3]
reference      = args[4]
info_threshold = args[5]
maf_threshold  = args[6]
out_prefix     = args[7]

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

sum_stats = sum_stats %>% mutate(MUNGED_FILE = paste0(TRAIT, ".sumstats.gz"))
LDSCoutput = ldsc(traits = sum_stats$MUNGED_FILE,
                  sample.prev = sum_stats$SAMPLE_PREV,
                  population.prev = sum_stats$POP_PREV,
                  ld = ld_folder,
                  wld = ld_folder,
                  trait.names = sum_stats$TRAIT,
                  ldsc.log = out_prefix,
                  stand = TRUE)

saveRDS(LDSCoutput, paste0(out_prefix, "_LDSC.rds"))

# ----- STEP 3: SPLIT SUMMARY STATS FOR GWAS -----

for (i in 1:nrow(sum_stats)) {
  dir_name = paste0(out_prefix, "_", i)
  dir.create(dir_name, showWarnings = TRUE)
  file = fread(sum_stats$FILE_PATH[i], data.table = FALSE, header = T)
  for (chromosome in c(1:22)) {
    file_chr = file %>% filter(CHR == chromosome)
    write.table(file_chr, 
                paste0(dir_name,
                       "/",
                       sum_stats$TRAIT[i],
                       "_chr",
                       chromosome,
                       ".split.assoc"),
                row.names = F,
                quote = F,
                sep = " ")
    }
}