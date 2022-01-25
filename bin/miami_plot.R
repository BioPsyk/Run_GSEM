#!/usr/bin/env Rscript

# - CODE TO MAKE A MIAMI PLOT for GENOMIC SEM OUTPUT HIGHLIGHTING ASSOCIATION 
# - TEST P-VALUE AND HETEROGENEITY TEST Q_VALUE ---------------------------- #

require(dplyr, quietly = TRUE)
require(data.table, quietly = TRUE)
require(qqman, quietly = TRUE)

args       = commandArgs(trailingOnly = TRUE)
assoc      = fread(args[1], header = TRUE)
assoc      = assoc %>% 
    filter(fail == 0) %>% 
    select(CHR, SNP, BP, Pval_Estimate, Q_pval)
out_prefix = args[2]
assoc_P    = assoc %>% 
    select(CHR, SNP, BP, Pval_Estimate) %>%
    rename(P = Pval_Estimate)
assoc_Q    = assoc %>% 
    select(CHR, SNP, BP, Q_pval) %>%
    rename(P = Q_pval)

# --------------------------- MIAMI PLOT ---------------------------#

png(paste0(out_prefix, "_MiamiPlot.png"), 
    res = 300, width = 10, height = 8, units = "in")

par(mfrow = c(2, 1))
par(mar = c(0, 5, 3, 3))
manhattan(assoc_P,
          ylim = c(0, max(assoc_P$P)), 
          annotatePval = 5e-8)
par(mar = c(5, 5, 3, 3))
manhattan(assoc_Q, 
          ylim = c(max(assoc_Q$P), 0),  
          xlab = "", 
          annotatePval = 5e-8)
dev.off()

# --------------------------- QQ PLOT ---------------------------

png(paste0(out_prefix, "_QQ.png"), 
    res = 300, width = 8, height = 8, units = "in")
par(mfrow = c(2, 1))
par(mar = c(0, 5, 3, 3))
qq(assoc_P$P)
par(mar = c(5, 5, 3, 3))
qq(assoc_Q$P)

dev.off()