#!/usr/bin/env Rscript
# Summarise GWAS-PW blocks and, when extracted SNP posteriors are supplied,
# construct nominal 90% Model-3 credible sets.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
if (length(a) < 3L) stop('Usage: Rscript 02_pw_gwas_postprocess.R BLOCKS.segbfs.gz EXTRACTED.bfs.gz OUTPUT_PREFIX [EXCLUDE_APOE=FALSE]')
segfile <- a[1]; snpfile <- a[2]; out <- a[3]
exclude_apoe <- if (length(a)>=4) as.logical(a[4]) else FALSE
if (is.na(exclude_apoe)) stop('EXCLUDE_APOE must be TRUE or FALSE')
seg <- fread(segfile)
cols <- c('chunk','chr','st','sp',paste0('PPA_',1:4))
if (length(setdiff(cols,names(seg)))) stop('Missing segment columns: ',paste(setdiff(cols,names(seg)),collapse=', '))
s <- seg[,.(block=chunk,chr,st,sp,P1=PPA_1,P2=PPA_2,P3=PPA_3,P4=PPA_4)]
s[,P0:=pmax(0,1-P1-P2-P3-P4)]
cols <- paste0('P',0:4)
s[,MAP_model:=max.col(.SD,ties.method='first')-1L,.SDcols=cols]
s[,MAP_prob:=do.call(pmax,c(.SD,list(na.rm=TRUE))),.SDcols=cols]
setcolorder(s,c('block','chr','st','sp',cols,'MAP_model','MAP_prob'))
fwrite(s,paste0(out,'_blocks_all.tsv'),sep='\t')
keep <- !(s$chr==6 & s$st<35000000 & s$sp>25000000)
if(exclude_apoe) keep <- keep & !(s$chr==19 & s$st<45900000 & s$sp>44400000)
b <- s[keep]
fwrite(b,paste0(out,'_blocks_filtered.tsv'),sep='\t')
for (threshold in c(0.50,0.80,0.90)) {
 fwrite(b[P3>=threshold],sprintf('%s_P3_ge_%0.2f.tsv',out,threshold),sep='\t')
}
fwrite(b[P4>=0.90],paste0(out,'_P4_ge_0.90.tsv'),sep='\t')
cat('Eligible blocks:',nrow(b),'| P3>=0.80:',nrow(b[P3>=0.80]),'| P4>=0.90:',nrow(b[P4>=0.90]),'\n')
if (snpfile == 'NONE') quit(save='no',status=0)
if (!file.exists(snpfile)) stop('SNP posterior file not found: ',snpfile)
x <- fread(snpfile)
needed <- c('chr','pos','chunk','PPA_3')
if (length(setdiff(needed,names(x)))) stop('Missing SNP columns: ',paste(setdiff(needed,names(x)),collapse=', '))
selected <- unique(b[P3>=0.80,block]); x <- x[chunk %in% selected]
sets <- list(); summaries <- list()
for (id in selected) {
 d <- copy(x[chunk==id]); d <- d[is.finite(PPA_3)&PPA_3>=0]
 d <- unique(d,by=c('chr','pos'))
 if(!nrow(d) || sum(d$PPA_3)<=0) next
 setorder(d,-PPA_3)
 d[,posterior_weight:=PPA_3/sum(PPA_3)]
 d[,cumulative_weight:=cumsum(posterior_weight)]
 n <- which(d$cumulative_weight>=0.90)[1L]
 cs <- d[seq_len(n)]; cs[,block:=id]
 sets[[as.character(id)]] <- cs
 summaries[[as.character(id)]] <- data.table(block=id,lead_chr=cs$chr[1],lead_pos=cs$pos[1],lead_variant_PPA3=cs$PPA_3[1],CS_size=n,CS_start=min(cs$pos),CS_end=max(cs$pos),mass=sum(cs$posterior_weight))
}
if (!length(sets)) {
 warning('No eligible model-3 credible sets: writing empty outputs')
 fwrite(data.table(chr=integer(),pos=integer(),block=integer(),PPA_3=numeric()),paste0(out,'_CS90_variants.tsv'),sep='\t')
 fwrite(data.table(block=integer(),CS_size=integer()),paste0(out,'_CS90_summary.tsv'),sep='\t')
} else {
 fwrite(rbindlist(sets,fill=TRUE),paste0(out,'_CS90_variants.tsv'),sep='\t')
 fwrite(rbindlist(summaries),paste0(out,'_CS90_summary.tsv'),sep='\t')
}
cat('Nominal 90% credible sets:',length(sets),'\n')
# Assumption: extracted PPA_3 values can be normalised within block to
# represent relative shared-variant posterior weights; validate against
# GWAS-PW output semantics before treating as calibrated credible sets.
