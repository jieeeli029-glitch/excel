############################################################
# ACTL90007 Assignment 2, Q3: Total population, 1961–2005, ages 40–90
# Adapted from the lecturer's LC_APC.R and Unit 4 pp.50–58.
# Method retained: StMoMo Poisson LC and Renshaw–Haberman (RH).
# No forecasts, log-rate SVD replacement, data trimming or model selection search.
############################################################
# Install once in R/RStudio if necessary:
# install.packages(c("StMoMo", "reshape2", "forecast"))
suppressPackageStartupMessages({library(StMoMo);library(reshape2);library(forecast)})
args_all <- commandArgs(trailingOnly=FALSE)
script_arg <- grep("^--file=",args_all,value=TRUE)
if(length(script_arg)) setwd(dirname(normalizePath(sub("^--file=","",script_arg[1]))))
# In RStudio, set the working directory to this script's folder before sourcing it.
dir.create("output",showWarnings=FALSE)
set.seed(90007) # reproducible numerical initialisation; model unchanged
ages.fit <- 40:90
years.fit <- 1961:2005

read_population <- function(country) {
  deaths_raw <- read.table(file.path("data",paste0(country,"_Deaths_1x1.txt")),skip=2,header=TRUE,stringsAsFactors=FALSE,na.strings=".")
  exposures_raw <- read.table(file.path("data",paste0(country,"_Exposures_1x1.txt")),skip=2,header=TRUE,stringsAsFactors=FALSE,na.strings=".")
  clean_age <- function(age) as.numeric(gsub("\\+","",age))
  deaths_raw$Age_num <- clean_age(deaths_raw$Age)
  exposures_raw$Age_num <- clean_age(exposures_raw$Age)
  deaths <- subset(deaths_raw,Year %in% years.fit & Age_num %in% ages.fit)
  exposures <- subset(exposures_raw,Year %in% years.fit & Age_num %in% ages.fit)
  stopifnot(nrow(deaths)==51*45,nrow(exposures)==51*45,
            !anyDuplicated(deaths[c("Year","Age_num")]),
            !anyDuplicated(exposures[c("Year","Age_num")]))
  # Total (both sexes), not Male as in the classroom demonstration.
  Dxt <- acast(deaths,Age_num~Year,value.var="Total")
  Ext <- acast(exposures,Age_num~Year,value.var="Total")
  Dxt <- as.matrix(Dxt[as.character(ages.fit),as.character(years.fit)])
  Ext <- as.matrix(Ext[as.character(ages.fit),as.character(years.fit)])
  stopifnot(identical(dim(Dxt),c(51L,45L)),all(is.finite(Dxt)),
            all(is.finite(Ext)),all(Dxt>0),all(Ext>0))
  list(Dxt=Dxt,Ext=Ext)
}

# Unscaled Poisson deviance residuals on the same scale for all models.
# Using scale=FALSE avoids changing the colour scale through model-specific
# dispersion estimates. The sign of each residual is the sign of D - fitted D.
get_residuals <- function(model) residuals(model,scale=FALSE)$residuals

export_fit <- function(model,label) {
  stopifnot(isTRUE(model$conv))
  write.csv(data.frame(age=ages.fit,ax=model$ax,bx=model$bx[,1]),
            file.path("output",paste0(label,"_age_parameters.csv")),row.names=FALSE)
  write.csv(data.frame(year=years.fit,kt=as.numeric(model$kt[1,])),
            file.path("output",paste0(label,"_period_parameters.csv")),row.names=FALSE)
  write.csv(fitted(model,type="rates"),file.path("output",paste0(label,"_fitted_rates.csv")))
  write.csv(get_residuals(model),file.path("output",paste0(label,"_deviance_residuals.csv")))
  saveRDS(model,file.path("output",paste0(label,"_fit.rds")))
}

# Heatmap coordinates: year on horizontal axis, age on vertical axis.
# Same-cohort observations lie along diagonals t-x=constant.
# A common symmetric colour limit is used within each comparison.
heat <- function(mat,title,limit,sign_only=FALSE) {
  if(sign_only) {
    z <- sign(mat);br <- c(-1.5,-.5,.5,1.5);cols <- c("#2166AC","#F7F7F7","#B2182B")
  } else {
    z <- mat; z[z < -limit] <- -limit; z[z > limit] <- limit; br<-seq(-limit,limit,length.out=102)
    cols<-colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101)
  }
  image(years.fit,ages.fit,t(z),breaks=br,col=cols,xlab="Calendar year t",ylab="Age x",main=title,useRaster=TRUE)
  for(cohort in seq(1880,1960,20)) {
    xx <- years.fit-cohort;ok<-xx>=40 & xx<=90
    lines(years.fit[ok],xx[ok],col="#555555",lty=3,lwd=.6)
  }
  box()
}
heat_pair <- function(mats,titles,file,sign_only=FALSE) {
  lim <- max(abs(unlist(mats)),na.rm=TRUE)
  draw <- function() {
    par(mfrow=c(1,2),mar=c(4,4,3,1),oma=c(2.5,0,0,0))
    for(j in 1:2)heat(mats[[j]],titles[j],lim,sign_only)
    mtext(if(sign_only) "Blue: observed deaths below fitted; red: above fitted. Dotted lines: constant birth cohorts."
          else sprintf("Common colour scale: blue = negative, white = 0, red = positive; range +/- %.2f. Dotted lines: constant cohorts.",lim),side=1,outer=TRUE,cex=.75)
  }
  pdf(file.path("output",paste0(file,".pdf")),width=12,height=5.8);draw();dev.off()
  png(file.path("output",paste0(file,".png")),width=1800,height=870,res=150);draw();dev.off()
}

pop <- lapply(c("AUS","UK"),read_population);names(pop)<-c("AUS","UK")
LC <- lc(link="log")
LC.fits <- lapply(pop,function(dat) fit(LC,Dxt=dat$Dxt,Ext=dat$Ext,
                                       ages=ages.fit,years=years.fit,verbose=TRUE))
for(country in names(pop))export_fit(LC.fits[[country]],paste0(country,"_LC"))

# Q3(a): parameter plots, same axis scales within each parameter column.
plot_parameters <- function() {
  par(mfrow=c(2,3),mar=c(4,4,3,1))
  for(country in names(pop)) {
    fit<-LC.fits[[country]];label<-if(country=="AUS")"Australia" else "United Kingdom"
    plot(ages.fit,fit$ax,type="l",col="#17365D",xlab="Age",ylab=expression(a[x]),main=paste(label,"LC: a[x]"),ylim=range(unlist(lapply(LC.fits,`[[`,"ax"))))
    plot(ages.fit,fit$bx[,1],type="l",col="#17365D",xlab="Age",ylab=expression(b[x]),main=paste(label,"LC: b[x]"),ylim=range(unlist(lapply(LC.fits,`[[`,"bx"))))
    plot(years.fit,fit$kt[1,],type="l",col="#17365D",xlab="Year",ylab=expression(k[t]),main=paste(label,"LC: k[t]"),ylim=range(unlist(lapply(LC.fits,`[[`,"kt"))))
  }
}
pdf("output/Q3a_LC_parameters.pdf",width=12,height=7);plot_parameters();dev.off()
png("output/Q3a_LC_parameters.png",width=1800,height=1050,res=150);plot_parameters();dev.off()
LC.res <- lapply(LC.fits,get_residuals)
heat_pair(LC.res,c("Australia: LC residuals","United Kingdom: LC residuals"),"Q3a_LC_residual_heatmaps")
heat_pair(LC.res,c("Australia: LC residual signs","United Kingdom: LC residual signs"),"Q3a_LC_residual_signs",TRUE)
cohort_means <- function(mat) {
  cohort<-outer(ages.fit,years.fit,function(x,t)t-x)
  aggregate(as.vector(mat),list(cohort=as.vector(cohort)),function(z)c(mean=mean(z),n=length(z)))
}
for(country in names(pop)) {
  zz<-cohort_means(LC.res[[country]])
  write.csv(data.frame(cohort=zz$cohort,mean_residual=zz$x[,1],n=zz$x[,2]),paste0("output/",country,"_LC_cohort_means.csv"),row.names=FALSE)
}
if("--lc-only" %in% commandArgs(trailingOnly=TRUE))quit(save="no",status=0)

# Q3(b): UK chosen after examining the part-(a) diagonal residual patterns.
# Retain the lecturer's age-modulated RH model, not an unrelated APC variant.
selected_country <- "UK"
RH <- rh(link="log",cohortAgeFun="NP",approxConst=TRUE)
RH.fit <- fit(RH,Dxt=pop[[selected_country]]$Dxt,Ext=pop[[selected_country]]$Ext,
              ages=ages.fit,years=years.fit,verbose=TRUE)
export_fit(RH.fit,paste0(selected_country,"_RH"))
write.csv(data.frame(age=ages.fit,cx=as.numeric(RH.fit$b0x)),"output/UK_RH_cohort_age_sensitivity.csv",row.names=FALSE)
write.csv(data.frame(cohort=as.numeric(names(RH.fit$gc)),hc=as.numeric(RH.fit$gc)),"output/UK_RH_cohort_factor.csv",row.names=FALSE)
RH.res <- get_residuals(RH.fit)
heat_pair(list(LC.res[[selected_country]],RH.res),c("United Kingdom: LC residuals","United Kingdom: RH residuals"),"Q3b_UK_LC_RH_residual_heatmaps")
heat_pair(list(LC.res[[selected_country]],RH.res),c("United Kingdom: LC residual signs","United Kingdom: RH residual signs"),"Q3b_UK_LC_RH_residual_signs",TRUE)
plot_RH <- function() {
  par(mfrow=c(1,2),mar=c(4,4,3,1))
  plot(ages.fit,RH.fit$b0x,type="l",xlab="Age",ylab=expression(c[x]),main="RH cohort age sensitivity")
  plot(as.numeric(names(RH.fit$gc)),RH.fit$gc,type="l",xlab="Birth cohort t-x",ylab=expression(h[c]),main="RH cohort factor")
}
pdf("output/Q3b_UK_RH_cohort_parameters.pdf",width=11,height=4.8);plot_RH();dev.off()
png("output/Q3b_UK_RH_cohort_parameters.png",width=1650,height=720,res=150);plot_RH();dev.off()
models<-list(AUS_LC=LC.fits$AUS,UK_LC=LC.fits$UK,UK_RH=RH.fit)
diagnostics<-do.call(rbind,lapply(names(models),function(n) {
 f<-models[[n]];r<-get_residuals(f)
 data.frame(model=n,nobs=f$nobs,npar=f$npar,deviance=f$deviance,
            logLik=f$loglik,RMSE_log_rate=sqrt(mean((log(f$Dxt/f$Ext)-log(fitted(f,type="rates")))^2)),
            max_abs_deviance_residual=max(abs(r)),converged=f$conv)
}))
write.csv(diagnostics,"output/model_diagnostics.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"output/sessionInfo.txt")
writeLines(capture.output(lapply(models,summary)),"output/model_summaries.txt")
print(diagnostics)
cat("Completed Q3. No future mortality forecasts are required or produced.\n")
