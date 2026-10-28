#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(coloc)})
# Run: Rscript scripts/05_colocalisation.R config/coloc_config.tsv
cfgargs<-commandArgs(TRUE)
if(length(cfgargs)!=1L) stop('Usage: Rscript 05_colocalisation.R CONFIG.tsv')
cfg<-fread(cfgargs[1]); required<-c('trait1_file','trait2_file','blocks_file','trait1_type','trait2_type','trait1_sdY','trait2_sdY','trait1_N','trait2_N','trait1_case_fraction','trait2_case_fraction','output_prefix','p12')
if(length(setdiff(required,names(cfg))) || nrow(cfg)!=1L) stop('Config requires one row and: ',paste(required,collapse=', '))
cfg<-as.list(cfg[1]); A<-fread(cfg$trait1_file); B<-fread(cfg$trait2_file); blocks<-fread(cfg$blocks_file)
# Input standard: chr,pos,snp,EA,OA,beta,se,p,MAF for each GWAS.
needed<-c('chr','pos','snp','EA','OA','beta','se','p','MAF')
for(d in list(A,B)) if(length(setdiff(needed,names(d)))) stop('GWAS missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
if(length(setdiff(c('block','chr','st','sp','P3','P4'),names(blocks)))) stop('Block table missing fields')
if(anyDuplicated(A$snp)||anyDuplicated(B$snp)) stop('Duplicate SNP IDs: resolve before analysis')
if(!all(c(cfg$trait1_type,cfg$trait2_type)%in%c('quant','cc'))) stop('type must be quant or cc')
make_dataset<-function(d,typ,sdY,N,frac){
 z<-list(beta=d$beta,varbeta=d$se^2,snp=d$snp,type=typ)
 if(typ=='quant'){
  if(is.finite(sdY) && sdY>0) z$sdY<-sdY
  else if(is.finite(N) && N>0){z$N<-N;z$MAF<-d$MAF}
  else stop('Quant trait needs documented sdY, or N and MAF for sdY estimation')
 } else {
  if(!is.finite(N)||N<=0||!is.finite(frac)||frac<=0||frac>=1) stop('Case-control trait requires N and case fraction')
  z$N<-N;z$s<-frac
 }
 z
}
calc<-function(bb,signal){
 w<-blocks[block==bb][1L]
 a<-A[chr==w$chr & pos>=max(1,w$st-250000) & pos<=w$sp+250000]
 b<-B[chr==w$chr & pos>=max(1,w$st-250000) & pos<=w$sp+250000]
 m<-merge(a,b,by='snp',suffixes=c('.a','.b'))
 if(!nrow(m)) return(data.table(block=bb,signal=signal,n_snps=0,note='no overlap'))
 m<-m[chr.a==chr.b & pos.a==pos.b]
 for(nm in c('EA.a','OA.a','EA.b','OA.b')) set(m,j=nm,value=toupper(m[[nm]]))
 # Exact strand-matched alleles only; palindromes removed even when aligned.
 ok<-m$EA.a==m$EA.b & m$OA.a==m$OA.b
 flip<-m$EA.a==m$OA.b & m$OA.a==m$EA.b
 amb<-paste0(m$EA.a,m$OA.a)%in%c('AT','TA','CG','GC') | paste0(m$EA.b,m$OA.b)%in%c('AT','TA','CG','GC')
 m<-m[(ok|flip)&!amb]
 if(!nrow(m)) return(data.table(block=bb,signal=signal,n_snps=0,note='no aligned SNPs'))
 m[EA.a==OA.b & OA.a==EA.b,beta.b:=-beta.b]
 m<-m[is.finite(beta.a)&is.finite(beta.b)&is.finite(se.a)&is.finite(se.b)&se.a>0&se.b>0&is.finite(MAF.a)&is.finite(MAF.b)&MAF.a>=0.01&MAF.b>=0.01&MAF.a<=0.5&MAF.b<=0.5&is.finite(p.a)&p.a<0.05]
 if(nrow(m)<100) return(data.table(block=bb,signal=signal,n_snps=nrow(m),note='fewer than 100 QC SNPs'))
 d1<-data.table(snp=m$snp,beta=m$beta.a,se=m$se.a,MAF=m$MAF.a)
 d2<-data.table(snp=m$snp,beta=m$beta.b,se=m$se.b,MAF=m$MAF.b)
 f<-tryCatch(coloc.abf(make_dataset(d1,cfg$trait1_type,as.numeric(cfg$trait1_sdY),as.numeric(cfg$trait1_N),as.numeric(cfg$trait1_case_fraction)),make_dataset(d2,cfg$trait2_type,as.numeric(cfg$trait2_sdY),as.numeric(cfg$trait2_N),as.numeric(cfg$trait2_case_fraction)),p1=1e-4,p2=1e-4,p12=as.numeric(cfg$p12)),error=function(e)e)
 if(inherits(f,'error'))return(data.table(block=bb,signal=signal,n_snps=nrow(m),note=conditionMessage(f)))
 s<-f$summary
 data.table(block=bb,signal=signal,n_snps=nrow(m),PP.H0=unname(s['PP.H0.abf']),PP.H1=unname(s['PP.H1.abf']),PP.H2=unname(s['PP.H2.abf']),PP.H3=unname(s['PP.H3.abf']),PP.H4=unname(s['PP.H4.abf']),note='ok')
}
# MHC chr6:25-35 Mb removed in this generic implementation.
blocks<-blocks[!(chr==6 & st<35000000 & sp>25000000)]
shared<-unique(blocks[P3>=0.90,block]); distinct<-unique(blocks[P4>=0.90,block])
res<-rbindlist(c(lapply(shared,calc,signal='P3_shared'),lapply(distinct,calc,signal='P4_distinct')),fill=TRUE)
if(!nrow(res)) stop('No blocks meeting the >=0.90 criteria')
fwrite(res,paste0(cfg$output_prefix,'_coloc.tsv'),sep='\t')
cat('Wrote ',cfg$output_prefix,'_coloc.tsv\n',sep='')
