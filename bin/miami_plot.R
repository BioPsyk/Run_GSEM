#!/usr/bin/env Rscript

#- CODE TO MAKE A MIAMI PLOT for GENOMIC SEM OUTPUT HIGHLIGHTING ASSOCIATION 
#- TEST P-VALUE AND HETEROGENEITY TEST Q_VALUE -----------------------------

require(dplyr, quietly = TRUE)
require(data.table, quietly = TRUE)
require(ggplot2, quietly = TRUE)

args       = commandArgs(trailingOnly = TRUE)
assoc      = fread(args[1], header = TRUE)
assoc      = assoc %>% 
    filter(fail == 0) %>% 
    select(CHR, BP, Pval_Estimate, Q_pval)
out_prefix = args[2]
assoc_P    = assoc %>% select(CHR, BP, Pval_Estimate) %>% 
    mutate(Test = "Association") %>% 
    mutate(Highlight = ifelse(Pval_Estimate <= 5e-8, 1, 0)) %>% 
    rename(P = Pval_Estimate)
assoc_Q    = assoc %>% select(CHR, BP, Q_pval) %>% 
    mutate(Test = "Heterogeneity") %>% 
    mutate(Highlight = ifelse(Q_pval <= 5e-8, 1, 0)) %>%
    rename(P = Q_pval)
assoc      = rbind(assoc_P, assoc_Q)

png(paste0(out_prefix, "_MiamiPlot.png"), 
    res = 300, width = 10, height = 8, units = "in")

ggplot(assoc, aes(x = BP, y = -1 * log10(P), color = Highlight)) + 
    geom_point() + 
    theme_bw() + 
    facet_grid(CHR ~ Test) + 
    theme(axis.text.x = element_blank(), legend.position = "none") +
    xlab("") +
    ylab("-log10(P)") + 
    geom_hline(yintercept = -1 * log10(5e-8), lty = 2, color = "blue") +
    geom_hline(yintercept = -1 * log10(1e-6), lty = 2, color = "green") + 
    scale_color_manual(values = c("black", "red"))
    
dev.off()