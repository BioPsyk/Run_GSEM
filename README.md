# Run Genomic SEM common factor GWAS

These two scripts, together implement the common factor GWAS using GenomicSEM (Grotzinger et. al 2019) on a set of summary statistics for which one may wish to investigate the presence of a latent common underlying factor.

Read more about GenomicSEM common factor GWAS [here](https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS) and their Nature Human Behavior paper [here](https://www.nature.com/articles/s41562-019-0566-x)

## Script1 USAGE

    Rscript run_gsem_stage1.R <file_of_sumstats>
        <path_to_hapmap3_snps_list]
        <path_to_ldsc_ld_folder> 
        <path_to_1KG_reference_dataset> 
        <info_score_threshold> 
        <maf_threshold> 
        <output_prefix>

## Script 1 OUTPUTs

    <out_prefix>_LDSC.rds LDSC Covariance Matrix between traits in the summary stats files
    <out_prefix>_ldsc.log 
    <out_prefix>_1/*_chr1.split.assoc 
    <out_prefix>_2/*_chr2.split.assoc ... 
    <out_prefix>_22/_chr22.split.assoc : Folders with summary stats split by chromosome

## Script 2 USAGE

    Rscript run_gsem_stage2.R <file_of_sumstats>
        <path_to_ldsc_covariance_from_stage1_out>
        <path_to_1KG_reference_dataset> 
        <info_score_threshold> 
        <maf_threshold> 
        <output_prefix>
        <chromosome>

## Script 2 OUTPUTs

    <out_prefix>_COMMON_FACTOR_GWAS_chr<chromosome>.txt: Summary stats from common factor GWAS output for specified chromosome

 Output format described [here](https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS) in detail

### FILE FORMATS

#### The format for the file of sumstats must have the following columns

    TRAIT : Name of the trait
    FILE_PATH  : Path to corresponding sumstats file
    N : Sample size
    SAMPLE_PREV : N_CASES / N
    POPULATION_PREV : Population based prevalence of the trait

#### Remaining files, including HapMap3 SNP lists, 1000 genomes reference etc can be downloaded from the author's page [here](https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v)

### NOTE 1: Before running the GWAS in stage 2, it is a good idea to interactively fit the model without SNP effects as shown [here,](https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects)to make sure it fits the data

### NOTE 2: Stage 2 requires high memory, I've had success with 180g and a wall time of 4-5 days

Results can be plotted using existing code and libraries for Manhattan and QQ plots
