#!/usr/bin/env Rscript
# Generic cross-sectional PRS/protein/EpiScore and network analysis.
# Usage: Rscript 05_EpiScore_PRS_Network.R
suppressPackageStartupMessages({library(data.table); library(igraph)})

# ================== USER SETTINGS ==================
# Edit these settings to match your dataset. One row per person per wave.
input_file <- 'data/cross_sectional.tsv'
output_dir <- 'results/episcore_cross_sectional'
id_col <- 'participant_id'; cohort_col <- 'cohort'; wave_col <- 'wave'
fi_col <- 'FI'; age_col <- 'age'; sex_col <- 'sex'
pc_cols <- c('PC1','PC2','PC3','PC4')
prs_cols <- c('PRS_RA','PRS_CP')
episcore_cols <- c('Epi_A','Epi_B')
protein_cols <- c('protein_A','protein_B') # optional measured proteins
min_n <- 50
fdr_cutoff <- 0.05
correlation_threshold <- 0.30
# Set TRUE only if inputs have NOT already been technically adjusted.
adjust_episcores <- FALSE
technical_cols <- c('set','array','hybridisation_date')
technical_mode <- 'fixed'
# ===================================================

required <- c('input_file','output_dir','id_col','cohort_col','wave_col','fi_col',
              'age_col','sex_col','pc_cols','prs_cols','episcore_cols')
if(length(setdiff(required,ls()))) stop('Config missing: ', paste(setdiff(required,ls()),collapse=', '))
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
d <- fread(input_file)
reqcols <- unique(c(id_col,cohort_col,wave_col,fi_col,age_col,sex_col,pc_cols,prs_cols,episcore_cols))
if(!all(reqcols %in% names(d))) stop('Missing columns: ',paste(setdiff(reqcols,names(d)),collapse=', '))
if(!exists('min_n')) min_n <- 50L
if(!exists('fdr_cutoff')) fdr_cutoff <- 0.05
if(!exists('correlation_threshold')) correlation_threshold <- 0.30
if(!exists('protein_cols')) protein_cols <- character()
if(!exists('technical_cols')) technical_cols <- character()
if(!exists('adjust_episcores')) adjust_episcores <- FALSE
if(!exists('technical_mode')) technical_mode <- 'fixed'
if(!technical_mode %in% c('fixed')) stop('Template supports fixed-effect technical residualisation only; see README')
if(length(setdiff(c(protein_cols,technical_cols), names(d)))) stop('Missing protein/technical columns')
# Standardise PRSs separately in each cohort and wave.
for(p in prs_cols) d[, (p):=if(sum(is.finite(get(p)))>1 && sd(get(p),na.rm=TRUE)>0)
  as.numeric(scale(get(p))) else rep(NA_real_,.N), by=c(cohort_col,wave_col)]
# OPTIONAL: technical residualisation separately per cohort/wave (requires verification
# against the exact technical adjustment used to generate the paper's results).
if(adjust_episcores){
  if(!length(technical_cols)) stop('Provide technical_cols when adjust_episcores=TRUE')
  for(e in episcore_cols){
    for(g in unique(d[,paste(get(cohort_col),get(wave_col),sep='::')])){
      idx <- which(paste(d[[cohort_col]],d[[wave_col]],sep='::')==g)
      sub <- d[idx]; cols <- c(e,technical_cols); sub <- sub[complete.cases(sub[,..cols])]
      if(nrow(sub)<min_n) {d[idx,(e):=NA_real_];next}
      form <- reformulate(technical_cols,response=e)
      fit <- lm(form,data=sub)
      residuals_all <- rep(NA_real_,length(idx))
      ok <- complete.cases(d[idx,..cols]); residuals_all[ok] <- residuals(fit)
      set(d,i=idx,j=e,value=residuals_all)
    }
  }
}
# Each result uses one exposure at a time (not a combined multi-PRS model).
fit_one <- function(sub, outcome, exposure, covs, direction){
  cols <- unique(c(outcome,exposure,covs)); x <- sub[complete.cases(sub[,..cols])]
  if(nrow(x)<min_n || sd(x[[exposure]])==0 || sd(x[[outcome]])==0) return(NULL)
  fit <- tryCatch(lm(reformulate(c(exposure,covs), response=outcome),data=x),error=function(e)NULL)
  if(is.null(fit)) return(NULL)
  cf <- coef(summary(fit)); if(!exposure %in% rownames(cf)) return(NULL)
  data.table(direction=direction, predictor=exposure, outcome=outcome, n=nrow(x),
             beta=unname(cf[exposure,1]),se=unname(cf[exposure,2]),p=unname(cf[exposure,4]))
}
prs_res <- list(); fi_res <- list(); j<-0L;k<-0L
for(coh in unique(d[[cohort_col]])) for(w in unique(d[[wave_col]])){
  sub<-d[get(cohort_col)==coh & get(wave_col)==w]
  for(p in prs_cols) for(e in unique(c(episcore_cols,protein_cols))){
    ans<-fit_one(sub,e,p,c(age_col,sex_col,pc_cols),'PRS_to_biomarker')
    if(!is.null(ans)){j<-j+1L;ans[,`:=`(cohort=as.character(coh),wave=as.character(w))];prs_res[[j]]<-ans}
  }
  for(e in episcore_cols){
    ans<-fit_one(sub,fi_col,e,c(age_col,sex_col),'EpiScore_to_FI')
    if(!is.null(ans)){k<-k+1L;ans[,`:=`(cohort=as.character(coh),wave=as.character(w))];fi_res[[k]]<-ans}
  }
}
prs_res<-rbindlist(prs_res,fill=TRUE);fi_res<-rbindlist(fi_res,fill=TRUE)
# FDR family: all biomarkers for each PRS within cohort and wave.
if(nrow(prs_res)) prs_res[,q:=p.adjust(p,'BH'),by=.(cohort,wave,predictor)]
if(nrow(fi_res)) {
  fi_res[,q:=p.adjust(p,'BH'),by=.(cohort,wave)]
  fi_res[,classification:=fifelse(q<fdr_cutoff & beta>0,'harmful',
                             fifelse(q<fdr_cutoff & beta<0,'protective','not_significant'))]
}
fwrite(prs_res,file.path(output_dir,'PRS_biomarker_associations.tsv'),sep='\t')
fwrite(fi_res,file.path(output_dir,'EpiScore_FI_associations.tsv'),sep='\t')
# Directed two-stage network: PRS -> EpiScore -> FI.
# Do not interpret graph direction as causal proof.
if(nrow(prs_res) && nrow(fi_res)){
  edges<-list(); z<-0L
  for(coh in unique(prs_res$cohort)) for(w in unique(prs_res$wave)){
    a<-prs_res[cohort==coh & wave==w & q<fdr_cutoff & outcome %in% episcore_cols]
    b<-fi_res[cohort==coh & wave==w & q<fdr_cutoff]
    intermediates<-intersect(a$outcome,b$predictor)
    if(!length(intermediates))next
    e1<-a[outcome %in% intermediates,.(from=paste0('PRS:',predictor),to=paste0('Epi:',outcome),weight=beta)]
    e2<-b[predictor %in% intermediates,.(from=paste0('Epi:',predictor),to='FI',weight=beta)]
    ed<-rbind(e1,e2);ed[,`:=`(cohort=coh,wave=w)];z<-z+1L;edges[[z]]<-ed
    graph<-graph_from_data_frame(ed[,.(from,to,weight)],directed=TRUE)
    pdf(file.path(output_dir,sprintf('PRS_EpiScore_FI_network_%s_%s.pdf',coh,w)),width=10,height=7)
    plot(graph,layout=layout_with_fr(graph),vertex.size=6,vertex.label.cex=.6,
         edge.width=1+pmin(5,abs(E(graph)$weight)*4))
    dev.off()
  }
  if(length(edges))fwrite(rbindlist(edges),file.path(output_dir,'network_edges.tsv'),sep='\t')
}
# Optional correlation network (not restricted to significant regression edges).
if(length(protein_cols)>1) for(coh in unique(d[[cohort_col]])) for(w in unique(d[[wave_col]])){
  sub<-d[get(cohort_col)==coh & get(wave_col)==w,..protein_cols]
  C<-cor(as.matrix(sub),use='pairwise.complete.obs');diag(C)<-0
  C[!is.finite(C) | abs(C)<=correlation_threshold]<-0
  if(!any(C!=0))next
  graph<-graph_from_adjacency_matrix(C,mode='undirected',weighted=TRUE,diag=FALSE)
  pdf(file.path(output_dir,sprintf('protein_correlation_network_%s_%s.pdf',coh,w)),width=10,height=7)
  plot(graph,layout=layout_with_fr(graph),vertex.size=6,vertex.label.cex=.6)
  dev.off()
}
cat('Cross-sectional analyses saved to',output_dir,'\n')
