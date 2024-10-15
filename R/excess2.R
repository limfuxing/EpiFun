#' Estimate Excess Number of Events in Treatment Group With Condition-Specific Exposure
#'
#' This function uses Poisson or Negative Binomial regression to estimate the Excess Number of Events (ENE) per 100,000 person-years. Agegroup-specific and Overall ENE statistic is calculated for each condition.
#'
#' @param file path to an Excel file containing the data. The file must contain multiple sheets with the following sheet naming convention: YYYY_Control_BD containing data for control population in year YYYY and YYYY_Treat_BD containing data for population under treatment. On top of these, the file must also contain the matching exposure that will act as denominator in the regression model. The naming format for these sheets containing exposures are YYYY_Control_AD and YYYY_Treat_AD. Inside each sheet, the condition name/ID must be specified in the first column and the first row contains variable names. Each row contains data for a condition across different agegroups (from column 2 onwards) with the last row contains the population size to be used as offset in the regression model. 
#' @param reg.model Regression model used to estimate the statistic. The default is Poisson regression with Negative Binomial regression as alternative.
#' @param wt.type Which population to use for calculating the 'Overall' Excess estimate as weighted average of the the age-specific Excess estimate. The default is to weigh using the treated (exposed) population. This gives a realistic assessment of the disease burden and spread for the exposed population, and also allows comparisons between the exposed and general populations, accounting for differences in age, so can be interpreted as the effect of exposures if the general population had the same age structure as the population with exposure.  
#' @param wt.var weighting variable to estimate overall effect combining age-specific effect (must be ordered according to age-category). If specified, it will override wt.type specification.
#' @param save whether to automatically save the output as a CSV file.  
#' @param save.plot whether to produce volcano plot and save it as a plotly object in an HTML file. 

#' @return A table containing the ENE statistic (per 100,000 person-years) for each condition and agegroup. 

#' @export


ENE2 <- function(file,reg.model=c('poisson','NB'),wt.type='treat',wt.var=NULL,save=TRUE,save.plot=TRUE) {

# process and add each file onto data frame and prepare data for volc plot
volc.df <- NULL
TEvent   <- TEvent.ctl <- TEvent.trt <- NULL
adj.est <- NULL

sheet <- sort(openxlsx::getSheetNames(file=file))

# control sheet names
sheet.ctrl <- sheet[grep('Control',sheet)]
sheet.ctrl.AD <- sheet.ctrl[grep('AD',sheet.ctrl)]
sheet.ctrl.BD <- sheet.ctrl[grep('BD',sheet.ctrl)]

# Treat sheet names
sheet.diab <- sheet[-grep('Control',sheet)]
sheet.diab.AD <- sheet.diab[grep('AD',sheet.diab)]
sheet.diab.BD <- sheet.diab[grep('BD',sheet.diab)]

# read control sheets
out.AD <- NULL
for(SH in sheet.ctrl.AD) {
  tmp=openxlsx::read.xlsx(file,sheet=SH,rowNames=TRUE)
  Year <- as.numeric(unlist(strsplit(SH,"_"))[1])
  out.AD<-rbind(out.AD,data.frame(Category=file,SubCat=rep(substring(rownames(tmp),1,1000),ncol(tmp)),agegrp=rep(colnames(tmp),
		rep(nrow(tmp),ncol(tmp))),count=c(unlist(tmp)),status=unlist(strsplit(SH,"_"))[2],Year=rep(Year,ncol(tmp))))
}

out.BD <- NULL
for(SH in sheet.ctrl.BD) {
  tmp=openxlsx::read.xlsx(file,sheet=SH,rowNames=TRUE)
  Year <- as.numeric(unlist(strsplit(SH,"_"))[1])
  out.BD <-  rbind(out.BD,data.frame(Category=file,SubCat=rep(substring(rownames(tmp),1,1000),ncol(tmp)),agegrp=rep(colnames(tmp),
		rep(nrow(tmp),ncol(tmp))),count=c(unlist(tmp)),status=unlist(strsplit(SH,"_"))[2],Year=rep(Year,ncol(tmp))))
}
  
# read treat sheets
out2.AD <- NULL
for(SH in sheet.diab.AD) {
  tmp=openxlsx::read.xlsx(file,sheet=SH,rowNames=TRUE)
  Year <- as.numeric(unlist(strsplit(SH,"_"))[1])
  out2.AD <-  rbind(out2.AD,data.frame(Category=file,SubCat=rep(substring(rownames(tmp),1,1000),ncol(tmp)),agegrp=rep(colnames(tmp),
		rep(nrow(tmp),ncol(tmp))),count=c(unlist(tmp)),status=unlist(strsplit(SH,"_"))[2],Year=rep(Year,ncol(tmp))))
}

out2.BD <- NULL
for(SH in sheet.diab.BD) {
  tmp=openxlsx::read.xlsx(file,sheet=SH,rowNames=TRUE)
  Year <- as.numeric(unlist(strsplit(SH,"_"))[1])
  out2.BD <-  rbind(out2.BD,data.frame(Category=file,SubCat=rep(substring(rownames(tmp),1,1000),ncol(tmp)),agegrp=rep(colnames(tmp),
		rep(nrow(tmp),ncol(tmp))),count=c(unlist(tmp)),status=unlist(strsplit(SH,"_"))[2],Year=rep(Year,ncol(tmp))))
}

# make sure rows are aligned
subcat = intersect(intersect(unique(out.AD$SubCat),unique(out.BD$SubCat)),intersect(unique(out2.AD$SubCat),unique(out2.BD$SubCat)))
out2.AD=out2.AD[out2.AD$SubCat %in% subcat,]
out2.BD=out2.BD[out2.BD$SubCat %in% subcat,]
out.AD=out.AD[out.AD$SubCat %in% subcat,]
out.BD=out.BD[out.BD$SubCat %in% subcat,]

# combine cases (DM) and controls (non-DM)
data.AD <- rbind(out.AD,out2.AD)
data.BD <- rbind(out.BD,out2.BD)
data.BD <- subset(data.BD, data.AD$count>0)
data.AD <- subset(data.AD, data.AD$count>0)


# check missing data
row.withmiss.BD <- which(is.na(data.BD$count) | is.infinite(data.BD$count))
df.toprint <- data.BD[,c('SubCat','agegrp','count','status', 'Year')]
row.withmiss.AD <- which(is.na(data.AD$count) | is.infinite(data.AD$count))
df2.toprint <- data.AD[,c('SubCat','agegrp','count','status', 'Year')]

if(length(row.withmiss.BD)>0) {
  print(paste0('The following data may have missing values. Please CHECK your input file!'))
  rownames(df.toprint) <- NULL
  print(df.toprint[row.withmiss.BD,])
}

if(length(row.withmiss.AD)>0) {
  print(paste0('The following EXPOSURE data may have missing values. Please CHECK your input file!'))
  rownames(df2.toprint) <- NULL
  print(df2.toprint[row.withmiss.AD,])
}

# fit model for each condition and calculate statistic
subc <- unique(data.BD$SubCat)
agegroup <- unique(data.BD$agegrp)
for(i in 1:length(subc)) {
    data.sub.AD <- subset(data.AD,SubCat==subc[i])
    data.sub.BD <- subset(data.BD,SubCat==subc[i])
    data.sub.BD$Year <- factor(data.sub.BD$Year)
    data.sub.BD$agegrp <- factor(data.sub.BD$agegrp)
    contrasts(data.sub.BD$Year) <- contr.sum
    PY <- data.sub.AD$count
    # model for agegroup-specific estimate
    if(reg.model[1]=='poisson')
     model.adj   <- glm(count~status*agegrp+Year+offset(log(PY)),family='poisson',data=data.sub.BD)
    if(reg.model[1]=='NB') {
     model.adj   <- tryCatch({MASS::glm.nb(count~status*agegrp+Year+offset(log(PY)),data=data.sub.BD)}, error = function(e) 
			{glm(count~status*agegrp+Year+offset(log(PY)),family='poisson',data=data.sub.BD)})
    }
    # model for overall estimate
    contrasts(data.sub.BD$agegrp) <- contr.sum
    X  <- model.matrix(~status*agegrp+Year,contr.arg=list(agegrp='contr.sum',Year='contr.sum'),data=data.sub.BD)
    if(reg.model[1]=='poisson')
     model2.adj   <- glm(count~X[,-1]+offset(log(PY)),family='poisson',data=data.sub.BD)
    if(reg.model[1]=='NB') {
     model2.adj   <- tryCatch({MASS::glm.nb(count~X[,-1]+offset(log(PY)),data=data.sub.BD)}, error = function(e) 
			{glm(count~X[,-1]+offset(log(PY)),family='poisson',data=data.sub.BD)})
    }
    # start calculating agegroup-specific statistics
    lev <- levels(data.sub.BD$agegrp)
    nlev <- length(lev)
    ncoef<- length(coef(model.adj))
    for(agegroup in lev) {
      beta0 = model.adj$coef[1] + (model.adj$coef[1+match(agegroup,lev)])*(match(agegroup,lev)>1)
      beta = model.adj$coef[2]  + (model.adj$coef[ncoef-nlev+match(agegroup,lev)])*(match(agegroup,lev)>1)
      if(ncol(vcov(model.adj)) >= (ncoef-nlev+match(agegroup,lev)) ) {
        vbeta0=vcov(model.adj)[1,1]+
  	      (vcov(model.adj)[1+match(agegroup,lev),1+match(agegroup,lev)] + 2*vcov(model.adj)[1,1+match(agegroup,lev)])*(match(agegroup,lev)>1)
        vbeta = diag(vcov(model.adj))[2] +
  	      (diag(vcov(model.adj))[ncoef-nlev+match(agegroup,lev)] + 2*vcov(model.adj)[2,ncoef-nlev+match(agegroup,lev)])*(match(agegroup,lev)>1) 
        covbbeta0 = vcov(model.adj)[1,2] + (vcov(model.adj)[2,1+match(agegroup,lev)] + vcov(model.adj)[1,ncoef-nlev+match(agegroup,lev)] +
			vcov(model.adj)[ncoef-nlev+match(agegroup,lev),1+match(agegroup,lev)])*(match(agegroup,lev)>1)
      }

      if(ncol(vcov(model.adj)) < (ncoef-nlev+match(agegroup,lev)) ) {
       vbeta0 <- vbeta <- covbbeta0 <- NA
      }

      Hosp.ctl <- exp(beta0)*100000 
      Hosp.trt <- exp(beta+beta0)*100000 
      SE.log.Hosp.ctl <- sqrt(vbeta0)
      tmp <- vbeta0+vbeta+2*covbbeta0
      SE.log.Hosp.trt <- sqrt(tmp)
      SE.Hosp.ctl     <- Hosp.ctl*SE.log.Hosp.ctl
      SE.Hosp.trt     <- Hosp.trt*SE.log.Hosp.trt
      ExHosp <- (exp(beta)-1)*exp(beta0)*100000
      var.ExHosp <- (exp(beta+beta0)^2*tmp + exp(beta0)^2*vbeta0 - 2*exp(beta+beta0)*exp(beta0)*(covbbeta0 + vbeta0))*100000^2
      SE.ExHosp  <- sqrt(var.ExHosp)
      
            adj.est <- rbind(adj.est,data.frame(SubCat=subc[i],agegrp=agegroup,
				logPR=beta,SE=sqrt(vbeta),logPR_ci_lower=beta-1.96*sqrt(vbeta),
				logPR_ci_upper=beta+1.96*sqrt(vbeta),ExEvent=ExHosp,ExEvent_ci_lower=ExHosp-1.96*SE.ExHosp,
				ExEvent_ci_upper=ExHosp+1.96*SE.ExHosp,Event.ctl=Hosp.ctl,Event.ctl_ci_lower=exp(log(Hosp.ctl)-1.96*SE.log.Hosp.ctl),
				Event.ctl_ci_upper=exp(log(Hosp.ctl)+1.96*SE.log.Hosp.ctl),Event.trt=Hosp.trt,
				Event.trt_ci_lower=exp(log(Hosp.trt)-1.96*SE.log.Hosp.trt),Event.trt_ci_upper=exp(log(Hosp.trt)+1.96*SE.log.Hosp.trt)))
    }
    # fixed CI with missing/infinite values
    adj.est$Event.ctl_ci_lower[(is.na(adj.est$Event.ctl_ci_lower) | is.infinite(adj.est$Event.ctl_ci_lower)) & adj.est$SubCat==subc[i]] <- 0
    adj.est$Event.ctl_ci_upper[(is.na(adj.est$Event.ctl_ci_upper) | is.infinite(adj.est$Event.ctl_ci_upper)) & adj.est$SubCat==subc[i]] <- 1
    adj.est$Event.trt_ci_lower[(is.na(adj.est$Event.trt_ci_lower) | is.infinite(adj.est$Event.trt_ci_lower)) & adj.est$SubCat==subc[i]] <- 0
    adj.est$Event.trt_ci_upper[(is.na(adj.est$Event.trt_ci_upper) | is.infinite(adj.est$Event.trt_ci_upper)) & adj.est$SubCat==subc[i]] <- 1
    # those with very small event rates, fixed upper CI at 1
    adj.est$Event.ctl_ci_upper[adj.est$Event.ctl<0.1 & adj.est$SubCat==subc[i]] <- 1
    adj.est$Event.trt_ci_upper[adj.est$Event.trt<0.1 & adj.est$SubCat==subc[i]] <- 1

    # get overall statistic
    if(is.null(wt.var)) {
     if(wt.type=='treat')
       data.sub.ctl = subset(data.sub.AD,status!='Control')
     if(wt.type!='treat')
       data.sub.ctl = subset(data.sub.AD,status=='Control')
    
    virtual.pop <- aggregate(data.sub.ctl$count,by=list(agegroup=data.sub.ctl$agegrp),sum)
    virtual.pop$x <- virtual.pop$x/sum(virtual.pop$x)
    virtual.pop <- virtual.pop[match(lev,virtual.pop$agegroup),]
    }
    
    if(!is.null(wt.var)) {
      virtual.pop <- data.frame(x=wt.var[match(lev,names(wt.var))])
    }
 
    interval   <- adj.est$Event.ctl_ci_upper[adj.est$SubCat==subc[i]]-adj.est$Event.ctl_ci_lower[adj.est$SubCat==subc[i]]
    interval   <- ifelse(interval>sqrt(.Machine$double.xmax),sqrt(.Machine$double.xmax),interval)
    SE.Hosp.ctl<-sqrt(sum(virtual.pop$x^2*(interval/3.92)^2))
    SE.Hosp.ctl<-ifelse(SE.Hosp.ctl> sqrt(.Machine$double.xmax), sqrt(.Machine$double.xmax), SE.Hosp.ctl)
    Hosp.ctl <- exp(weighted.mean(log(adj.est$Event.ctl[adj.est$SubCat==subc[i]]),w=virtual.pop$x))
    SE.log.Hosp.ctl <- SE.Hosp.ctl/Hosp.ctl
    interval   <- adj.est$Event.trt_ci_upper[adj.est$SubCat==subc[i]]-adj.est$Event.trt_ci_lower[adj.est$SubCat==subc[i]]
    SE.Hosp.trt<-sqrt(sum(virtual.pop$x^2*(interval/3.92)^2))
    SE.Hosp.trt<-ifelse(SE.Hosp.trt> sqrt(.Machine$double.xmax), sqrt(.Machine$double.xmax), SE.Hosp.trt)
    Hosp.trt <- exp(weighted.mean(log(adj.est$Event.trt[adj.est$SubCat==subc[i]]),w=virtual.pop$x))
    SE.log.Hosp.trt <- SE.Hosp.trt/Hosp.trt

    beta   <- log(Hosp.trt/Hosp.ctl)
    vbeta  <- SE.log.Hosp.trt^2 + SE.log.Hosp.ctl^2
    ExHosp <- Hosp.trt-Hosp.ctl
    SE.ExHosp<- sqrt(SE.Hosp.trt^2 + SE.Hosp.ctl^2)

    adj.est <- rbind(adj.est,data.frame(SubCat=subc[i],agegrp="Overall",
				logPR=beta,
                                SE=sqrt(vbeta),logPR_ci_lower=beta-1.96*sqrt(vbeta),
				logPR_ci_upper=beta+1.96*sqrt(vbeta),ExEvent=ExHosp,ExEvent_ci_lower=ExHosp-1.96*SE.ExHosp,
				ExEvent_ci_upper=ExHosp+1.96*SE.ExHosp,Event.ctl=Hosp.ctl,Event.ctl_ci_lower=exp(log(Hosp.ctl)-1.96*SE.log.Hosp.ctl),
				Event.ctl_ci_upper=exp(log(Hosp.ctl)+1.96*SE.log.Hosp.ctl),Event.trt=Hosp.trt,
				Event.trt_ci_lower=exp(log(Hosp.trt)-1.96*SE.log.Hosp.trt),Event.trt_ci_upper=exp(log(Hosp.trt)+1.96*SE.log.Hosp.trt)))
  }
  # ensure all CIs are finite
  adj.est$Event.ctl_ci_upper[is.infinite(adj.est$Event.ctl_ci_upper)] <- sqrt(.Machine$double.xmax)
  adj.est$Event.trt_ci_upper[is.infinite(adj.est$Event.trt_ci_upper)] <- sqrt(.Machine$double.xmax)

TEvent       <- aggregate(data.BD$count,by=list(SubCat=data.BD$SubCat,agegrp=data.BD$agegrp),sum)
TEvent.ctl   <- aggregate(out.BD$count,by=list(SubCat=out.BD$SubCat,agegrp=out.BD$agegrp),sum)
TEvent.trt   <- aggregate(out2.BD$count,by=list(SubCat=out2.BD$SubCat,agegrp=out2.BD$agegrp),sum)
# Overall
TEvent2       <- aggregate(data.BD$count,by=list(SubCat=data.BD$SubCat),sum)
TEvent2       <- data.frame(TEvent2[,-ncol(TEvent2)],agegrp = 'Overall', TEvent2[,ncol(TEvent2)])
TEvent2.ctl   <- aggregate(out.BD$count,by=list(SubCat=out.BD$SubCat),sum)
TEvent2.ctl   <- data.frame(TEvent2.ctl[,-ncol(TEvent2.ctl)],agegrp = 'Overall', TEvent2.ctl[,ncol(TEvent2.ctl)])
TEvent2.trt   <- aggregate(out2.BD$count,by=list(SubCat=out2.BD$SubCat),sum)
TEvent2.trt   <- data.frame(TEvent2.trt[,-ncol(TEvent2.trt)],agegrp = 'Overall', TEvent2.trt[,ncol(TEvent2.trt)])

# merge
colnames(TEvent2) <- colnames(TEvent2.ctl) <- colnames(TEvent2.trt) <- colnames(TEvent)
TEvent        <- rbind(TEvent,TEvent2)
TEvent.ctl    <- rbind(TEvent.ctl,TEvent2.ctl)
TEvent.trt    <- rbind(TEvent.trt,TEvent2.trt)

volc.df <- adj.est
# get Total, Ctl and Treat Events by SubCat and agegrp. This total to be used as size of the dots in volc plot
volc.df$TEvent <- TEvent$x[match(paste0(volc.df$SubCat,volc.df$agegrp),paste0(TEvent$SubCat,TEvent$agegrp))]
volc.df$TEvent.ctl  <- TEvent.ctl$x[match(paste0(volc.df$SubCat,volc.df$agegrp),paste0(TEvent.ctl$SubCat,TEvent.ctl$agegrp))]
volc.df$TEvent.trt  <- TEvent.trt$x[match(paste0(volc.df$SubCat,volc.df$agegrp),paste0(TEvent.trt$SubCat,TEvent.trt$agegrp))]

# get log pval
volc.df$logP <- pnorm(abs(volc.df$logPR/volc.df$SE),lower.tail=FALSE,log=TRUE) + log(2)
# get FDR (adjusted p-value)
volc.df$FDR  <- p.adjust(exp(volc.df$logP),method='BY')
# sort based on p-value (smallest to largest)
volc.df <- volc.df[order(volc.df$logP,decreasing=FALSE),]

# save table
if(save) 
  write.csv(volc.df,file='analysis_output.csv', row.names=FALSE, quote=FALSE)

# volcano plot (optional)
if(save.plot) {
 require(ggplot2)
 require(plotly)
 p <- ggplot(volc.df,aes(x=logPR,y=log(-logP),shape=agegrp,col=SubCat,size=TEvent)) + geom_point() + 
	ggtitle('Volcano plots: ALL subcategories') + xlab('log RR') + ylab('log(-log pvalue)') + guides(size="none")
 # save as plotly object in HTML format
 p <- ggplotly(p)
 htmlwidgets::saveWidget(p, "volcanoplot.html")
}
volc.df
}


