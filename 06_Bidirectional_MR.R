#!/usr/bin/env Rscript
# Generic bidirectional two-sample MR workflow accompanying manuscript Methods.
# INPUT: two harmonised-format GWAS tables per phenotype, NOT participant-level data.
# Inputs must be GRCh37/hg19, allele-oriented and prepared with valid frequencies.
# Run: Rscript 06_Bidirectional_MR.R
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# ========================= USER SETTINGS =========================
# Each GWAS file needs: SNP, chr, pos, ea, oa, beta, se, p, eaf, n.
# Binary traits additionally need ncase, ncontrol (columns, or add them upstream),
# and population prevalence (NOT the case fraction in the sampled GWAS).
# eaf is EFFECT-allele frequency (not MAF).
traits <- list(
  TraitA = list(file='data/trait_a.tsv', type='quant', prevalence=NA_real_),
  TraitB = list(file='data/trait_b.tsv', type='cc', prevalence=0.01)
)
# All ordered directions below, e.g. A -> B and B -> A.
pairs <- list(c('TraitA','TraitB'), c('TraitB','TraitA'))
plink_bin <- '/path/to/plink'               # PLINK 1.9 executable
ld_ref_bfile <- '/path/to/1000G_EUR'        # PLINK .bed/.bim/.fam prefix
output_dir <- 'results/mr'
p_threshold <- 5e-8
clump_kb <- 10000
clump_r2 <- 0.001
min_F <- 10
exclude_mhc <- TRUE                      # match analysis plan
mhc_start <- 25000000L; mhc_end <- 35000000L
run_presso <- TRUE
# ================================================================

stopifnot(length(traits) >= 2L)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
needed <- c('SNP','chr','pos','ea','oa','beta','se','p','eaf','n')
load_gwas <- function(name) {
  info <- traits[[name]]
  if(is.null(info) || !info$type %in% c('quant','cc')) stop('Unknown trait/type: ',name)
  d <- fread(info$file)
  if(length(setdiff(needed, names(d)))) stop(name, ': missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
  if(info$type=='cc' && !all(c('ncase','ncontrol') %in% names(d)))
    stop(name, ': binary GWAS must include ncase and ncontrol; supply study-correct metadata')
  if(info$type=='cc' && (!is.finite(info$prevalence) || info$prevalence<=0 || info$prevalence>=1))
    stop(name, ': supply justified population prevalence for Steiger testing')
  if(anyDuplicated(d$SNP)) stop(name, ': duplicated SNPs; resolve before running')
  d <- d[!is.na(SNP) & is.finite(beta) & is.finite(se) & se>0 &
         is.finite(p) & p>=0 & p<=1 & is.finite(eaf) & eaf>0 & eaf<1 &
         is.finite(n) & n>0]
  if(info$type=='cc') d <- d[is.finite(ncase)&is.finite(ncontrol)&ncase>0&ncontrol>0]
  d[, `:=`(ea=toupper(ea),oa=toupper(oa))]
  d
}

as_tsmr <- function(d, name, role) {
  isexp <- identical(role,'exposure')
  x <- data.table(
    SNP=d$SNP,
    beta=d$beta,
    se=d$se,
    effect_allele=d$ea,
    other_allele=d$oa,
    eaf=d$eaf,
    pval=d$p,
    samplesize=d$n
  )
  x <- as.data.frame(x)
  names(x)[-1] <- paste0(names(x)[-1], if(isexp) '.exposure' else '.outcome')
  if(isexp) {
    x$exposure <- name; x$id.exposure <- name
  } else {
    x$outcome <- name; x$id.outcome <- name
  }
  ty <- traits[[name]]$type
  x[[paste0('units.',role)]] <- if(ty=='cc') 'log odds' else 'SD/trait units'
  x[[paste0('type.',role)]] <- ty
  if(ty=='cc') {
    x[[paste0('ncase.',role)]] <- d$ncase
    x[[paste0('ncontrol.',role)]] <- d$ncontrol
    x[[paste0('prevalence.',role)]] <- traits[[name]]$prevalence
  }
  x
}

run_pair <- function(exp_name,out_name,raw) {
  prefix <- file.path(output_dir,paste0(exp_name,'_to_',out_name))
  e <- copy(raw[[exp_name]])
  if(exclude_mhc) e <- e[!(chr==6 & pos>=mhc_start & pos<=mhc_end)]
  e <- e[p<p_threshold]
  if(!nrow(e)) {warning('No GWS instruments: ',exp_name);return(NULL)}
  e <- as_tsmr(e,exp_name,'exposure')
  # TwoSampleMR clump_data requires SNP position fields and exposure identifier.
  pos <- raw[[exp_name]][match(e$SNP,raw[[exp_name]]$SNP),.(chr,pos)]
  e$chr_name <- pos$chr; e$chrom_start <- pos$pos
  e <- clump_data(e, clump_kb=clump_kb,clump_r2=clump_r2,
                  clump_p1=p_threshold, bfile=ld_ref_bfile,plink_bin=plink_bin)
  if(!nrow(e)) return(NULL)
  o <- raw[[out_name]][SNP %in% e$SNP]
  if(!nrow(o)) {warning('No overlapping outcome SNPs: ', exp_name,' -> ',out_name);return(NULL)}
  h <- harmonise_data(e,as_tsmr(o,out_name,'outcome'),action=2)
  h <- h[h$mr_keep %in% TRUE,,drop=FALSE]
  if(!nrow(h)) return(NULL)
  h$F_stat <- (h$beta.exposure / h$se.exposure)^2
  h <- h[is.finite(h$F_stat)&h$F_stat>min_F,,drop=FALSE]
  if(!nrow(h)) return(NULL)
  # Steiger is a per-variant FILTER here, not merely a diagnostic report.
  # Requires reliable allele frequencies, sample size, trait type and
  # population prevalence for case-control GWAS. Do not assume unknown data.
  h <- tryCatch(steiger_filtering(h),error=function(err)
    stop('Steiger failed for ',exp_name,' -> ',out_name,': ',conditionMessage(err),
         '. Check GWAS metadata and units.'))
  fwrite(as.data.table(h),paste0(prefix,'_harmonised_steiger.tsv'),sep='\t')
  h <- h[h$steiger_dir %in% TRUE & h$mr_keep %in% TRUE,,drop=FALSE]
  if(!nrow(h)) {warning('No SNPs passed Steiger: ',exp_name,' -> ',out_name);return(NULL)}
  fwrite(as.data.table(h),paste0(prefix,'_instruments.tsv'),sep='\t')

  methods <- c('mr_ivw','mr_egger_regression','mr_weighted_median','mr_simple_mode','mr_weighted_mode')
  # Estimators that require more instruments are skipped by TwoSampleMR.
  estimate <- mr(h,method_list=methods)
  fwrite(as.data.table(estimate),paste0(prefix,'_estimates.tsv'),sep='\t')
  if(nrow(h)>=2) {
    q <- tryCatch(mr_heterogeneity(h),error=function(e)NULL)
    if(!is.null(q)) fwrite(as.data.table(q),paste0(prefix,'_heterogeneity.tsv'),sep='\t')
  }
  if(nrow(h)>=3) {
    pleio <- tryCatch(mr_pleiotropy_test(h),error=function(e)NULL)
    if(!is.null(pleio)) fwrite(as.data.table(pleio),paste0(prefix,'_egger_intercept.tsv'),sep='\t')
    loo <- tryCatch(mr_leaveoneout(h),error=function(e)NULL)
    if(!is.null(loo)) fwrite(as.data.table(loo),paste0(prefix,'_leave_one_out.tsv'),sep='\t')
  }
  if(run_presso && nrow(h)>=4) {
    if(!requireNamespace('MRPRESSO',quietly=TRUE))
      warning('MRPRESSO unavailable; install separately to run this sensitivity analysis')
    else {
      pres <- tryCatch(MRPRESSO::mr_presso(BetaOutcome='beta.outcome',BetaExposure='beta.exposure',
        SdOutcome='se.outcome',SdExposure='se.exposure',OUTLIERtest=TRUE,DISTORTIONtest=TRUE,
        data=h,NbDistribution=1000,SignifThreshold=0.05),error=function(e)e)
      capture.output(print(pres),file=paste0(prefix,'_MRPRESSO.txt'))
    }
  }
  estimate
}

raw <- lapply(names(traits),load_gwas); names(raw) <- names(traits)
all_res <- lapply(pairs,function(p) run_pair(p[1],p[2],raw))
all_res <- rbindlist(lapply(all_res,as.data.table),fill=TRUE)
if(nrow(all_res)) {
  # Multiple testing on PRIMARY IVW estimates across all analysed directions.
  all_res[,p_fdr:=NA_real_]
  all_res[method=='Inverse variance weighted',p_fdr:=p.adjust(pval,method='BH')]
  fwrite(all_res,file.path(output_dir,'all_MR_estimates_FDR.tsv'),sep='\t')
}
cat('MR analysis completed. Outputs in: ',output_dir,'\n',sep='')
