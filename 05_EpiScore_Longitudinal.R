#!/usr/bin/env Rscript
# Generic longitudinal mixed-effects FI and protein trajectories.
# Usage: Rscript 06_EpiScore_Longitudinal.R
suppressPackageStartupMessages({library(data.table);library(lme4)})

# ================== USER SETTINGS ==================
# Edit these settings. One row per participant/wave; FI waves 1–5, biomarker waves 1–4.
input_file <- 'data/longitudinal.tsv'
output_dir <- 'results/episcore_longitudinal'
id_col <- 'participant_id'; wave_col <- 'wave'
age_col <- 'age'; sex_col <- 'sex'; fi_col <- 'FI'
protein_cols <- c('Epi_A','Epi_B')
min_people <- 50
max_protein_wave <- 4
max_fi_wave <- 5
min_slope_sd <- 0.01
age_scaling <- 1
# ===================================================

required<-c('input_file','output_dir','id_col','age_col','sex_col','wave_col','fi_col','protein_cols')
if(length(setdiff(required,ls())))stop('Missing config fields: ',paste(setdiff(required,ls()),collapse=', '))
if(!exists('min_slope_sd')) min_slope_sd<-0.01
if(!exists('min_people')) min_people<-50L
if(!exists('max_protein_wave')) max_protein_wave<-4
if(!exists('max_fi_wave')) max_fi_wave<-5
if(!exists('age_scaling')) age_scaling<-1  # 1 = per year; set explicitly for other units
if(!is.numeric(age_scaling)||age_scaling<=0)stop('age_scaling must be positive')
d<-fread(input_file);dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
cols<-unique(c(id_col,age_col,sex_col,wave_col,fi_col,protein_cols))
if(length(setdiff(cols,names(d))))stop('Missing input columns: ',paste(setdiff(cols,names(d)),collapse=', '))
setnames(d,c(id_col,age_col,sex_col,wave_col,fi_col),c('ID','age','sex','wave','FI'))
if(anyDuplicated(d[,.(ID,wave)]))stop('Duplicate ID/wave records: resolve before modelling')
# Center age for more interpretable intercepts; maintain slopes in units per year.
ref_age<-mean(d$age,na.rm=TRUE);d[,age_c:=(age-ref_age)/age_scaling]
# FI outcome on original scale by default; optional scale supplied explicitly upstream.
dfi<-d[wave<=max_fi_wave & complete.cases(d[,.(ID,age_c,sex,FI)])]
if(uniqueN(dfi$ID)<min_people)stop('Too few people for FI model')
fi_fit<-lmer(FI~age_c+sex+(1+age_c|ID),data=dfi,REML=TRUE)
fi_re<-as.data.table(ranef(fi_fit)$ID,keep.rownames='ID')
setnames(fi_re,c('(Intercept)','age_c'),c('FI_intercept_random','FI_slope_random'))
fi_re[,FI_slope:=unname(fixef(fi_fit)['age_c'])+FI_slope_random]
# Person baseline age and sex.
base<-d[order(age),.SD[1],by='ID'][,.(ID,baseline_age=age,sex)]
fi_re<-merge(fi_re,base,by='ID',all.x=TRUE)
fwrite(fi_re,file.path(output_dir,'FI_person_slopes.tsv'),sep='\t')
res<-list();outslopes<-list();idx<-0L
for(protein in protein_cols){
  dp<-d[wave<=max_protein_wave & complete.cases(d[,c('ID','age_c','sex',protein),with=FALSE]),
        c('ID','age_c','sex',protein),with=FALSE]
  if(uniqueN(dp$ID)<min_people)next
  setnames(dp,protein,'protein_value')
  if(sd(dp$protein_value)==0)next
  fit<-tryCatch(lmer(protein_value~age_c+sex+(1+age_c|ID),data=dp,REML=TRUE),error=function(e)NULL)
  if(is.null(fit))next
  re<-as.data.table(ranef(fit)$ID,keep.rownames='ID')
  setnames(re,c('(Intercept)','age_c'),c('protein_intercept_random','protein_slope_random'))
  re[,protein_slope:=unname(fixef(fit)['age_c'])+protein_slope_random]
  re[,protein_intercept:=unname(fixef(fit)['(Intercept)'])+protein_intercept_random]
  slope_sd<-sd(re$protein_slope,na.rm=TRUE)
  if(!is.finite(slope_sd)||slope_sd<min_slope_sd)next
  z<-merge(re,fi_re[,.(ID,FI_slope,baseline_age,sex)],by='ID')
  if(nrow(z)<min_people)next
  z[,protein_slope_z:=as.numeric(scale(protein_slope))]
  # Model mirrors manuscript: FI slope ~ protein slope + protein intercept + baseline age + sex.
  lmfit<-tryCatch(lm(FI_slope~protein_slope_z+protein_intercept+baseline_age+sex,data=z),error=function(e)NULL)
  if(is.null(lmfit))next
  co<-coef(summary(lmfit));if(!'protein_slope_z'%in%rownames(co))next
  idx<-idx+1L
  res[[idx]]<-data.table(protein=protein,n_people=nrow(z),protein_slope_SD=slope_sd,
                         beta=unname(co['protein_slope_z',1]),se=unname(co['protein_slope_z',2]),
                         p=unname(co['protein_slope_z',4]),singular=isSingular(fit))
  z[,protein:=protein];outslopes[[idx]]<-z
}
results<-rbindlist(res,fill=TRUE)
if(nrow(results)){
 results[,q:=p.adjust(p,'BH')]
 results[,classification:=fifelse(q<.05 & beta>0,'harmful',
                           fifelse(q<.05 & beta<0,'protective','not_significant'))]
}
fwrite(results,file.path(output_dir,'protein_FI_trajectory_associations.tsv'),sep='\t')
if(length(outslopes))fwrite(rbindlist(outslopes,fill=TRUE),file.path(output_dir,'protein_person_slopes.tsv'),sep='\t')
cat('Longitudinal results saved to',output_dir,'\n')
