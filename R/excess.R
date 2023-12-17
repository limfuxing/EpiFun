#' Estimate Excess Number of Events in Treatment Group
#'
#' This function uses Poisson or Negative Binomial regression to estimate the Excess Number of Events (ENE) per 100,000 person-years. Agegroup-specific and Overall ENE statistic is calculated for each condition.
#'
#' @param file path to an Excel file containing the data. The file mus contain multiple sheets with the following sheet naming convention: YYYY_Control containing data for control population in year YYYY and YYYY_Treat containing data for population under treatment. Inside each sheet, the condition name/ID must be specified in the first column and the first row contains variable names. Each row contains data for a condition across different agegroups (from column 2 onwards) with the last row contains the population size to be used as offset in the regression model. 
#' @param reg.model Regression model used to estimate the statistic. The default is Poisson regression with Negative Binomial regression as alternative.
#' @param save whether to automatically save the output as a CSV file.  
#' @param save.plot whether to produce volcano plot and save it as a plotly object in an HTML file. 

#' @return A table containing the ENE statistic (per 100,000 person-years) for each condition and agegroup. 

#' @export


ENE <- function(file,reg.model=c('poisson','NB'),save=TRUE,save.plot=TRUE) {

# process and add each file onto data frame and prepare data for volc plot
volc.df <- NULL
TEvent   <- TEvent.ctl <- TEvent.trt <- NULL
adj.est <- NULL

sheet <- sort(openxlsx::getSheetNames(file=file))
# control sheet names
sheet.ctrl <- sheet[grep('Control',sheet)]
# treat sheet names
sheet.trt <- sheet[-grep('Control',sheet)]
# read control sheets
  out <- NULL
  for(SH in sheet.ctrl) {
    tmp=openxlsx::read.xlsx(file,sheet=SH,rowNames=TRUE)
    Year <- as.numeric(unlist(strsplit(SH,"_"))[1])
    # get PY 
    py  = c(tmp[nrow(tmp),])
    # remove last row (PY)
    tmp = tmp[-nrow(tmp),]
    out <-  rbind(out,data.frame(SubCat=rep(substring(rownames(tmp),1,1000),ncol(tmp)),
		agegrp=rep(colnames(tmp),rep(nrow(tmp),ncol(tmp))),count=c(unlist(tmp)),status=unlist(strsplit(SH,"_"))[2],
		Year=rep(Year,ncol(tmp)),PY=unlist(rep(py,rep(nrow(tmp),ncol(tmp))))))
  }
  
  # read treat sheets
  out2 <- NULL
  for(SH in sheet.trt) {
    tmp=openxlsx::read.xlsx(file,sheet=SH,rowNames=TRUE)
    Year <- as.numeric(unlist(strsplit(SH,"_"))[1])
    # get PY 
    py  = c(tmp[nrow(tmp),])
    # remove last row (PY)
    tmp = tmp[-nrow(tmp),]
    out2 <- rbind(out2,data.frame(SubCat=rep(substring(rownames(tmp),1,1000),ncol(tmp)),
		agegrp=rep(colnames(tmp),rep(nrow(tmp),ncol(tmp))),count=c(unlist(tmp)),status=unlist(strsplit(SH,"_"))[2],
		Year=rep(Year,ncol(tmp)),PY=unlist(rep(py,rep(nrow(tmp),ncol(tmp))))))
  }
  
  # make sure rows are aligned
  subcat = intersect(unique(out$SubCat),unique(out2$SubCat))
  out2=out2[out2$SubCat %in% subcat,]
  out=out[out$SubCat %in% subcat,]
  # combine treat and controls data
  data.all <- rbind(out,out2)
  # fit model for each condition and calculate statistic
  subc <- unique(data.all$SubCat)
  agegroup <- unique(data.all$agegrp)
  for(i in 1:length(subc)) {
    data.sub  <- data.all[data.all$SubCat==subc[i],]
    data.sub$Year <- factor(data.sub$Year)
    data.sub$agegrp <- factor(data.sub$agegrp)
    contrasts(data.sub$Year) <- contr.sum
    # model for agegroup-specific estimate
    if(reg.model[1]=='poisson')
     model.adj   <- glm(count~status*agegrp+Year+offset(log(PY)),family='poisson',data=data.sub)
    if(reg.model[1]=='NB')
     model.adj   <- MASS::glm.nb(count~status*agegrp+Year+offset(log(PY)),data=data.sub)
    # model for overall estimate
    contrasts(data.sub$agegrp) <- contr.sum
    X  <- model.matrix(~status*agegrp+Year,contr.arg=list(agegrp='contr.sum',Year='contr.sum'),data=data.sub)
    if(reg.model[1]=='poisson')
     model2.adj   <- glm(count~X[,-1]+offset(log(PY)),family='poisson',data=data.sub)
    if(reg.model[1]=='NB')
     model2.adj   <- MASS::glm.nb(count~X[,-1]+offset(log(PY)),data=data.sub)
    # start calculating agegroup-specific statistics
    lev <- levels(data.sub$agegrp)
    nlev <- length(lev)
    ncoef<- length(coef(model.adj))
    for(agegroup in lev) {
      beta0 = summary(model.adj)$coef[1,1] + (summary(model.adj)$coef[1+match(agegroup,lev),1])*(match(agegroup,lev)>1)
      vbeta0=vcov(model.adj)[1,1]+
  	      (vcov(model.adj)[1+match(agegroup,lev),1+match(agegroup,lev)] + 2*vcov(model.adj)[1,1+match(agegroup,lev)])*(match(agegroup,lev)>1)
      beta = summary(model.adj)$coef[2,1]  + (summary(model.adj)$coef[ncoef-nlev+match(agegroup,lev),1])*(match(agegroup,lev)>1)
      vbeta = diag(vcov(model.adj))[2] +
  	      (diag(vcov(model.adj))[ncoef-nlev+match(agegroup,lev)])*(match(agegroup,lev)>1)
      covbbeta0 = vcov(model.adj)[1,2] + (vcov(model.adj)[1,1+match(agegroup,lev)] + vcov(model.adj)[2,ncoef-nlev+match(agegroup,lev)] +
			vcov(model.adj)[ncoef-nlev+match(agegroup,lev),1+match(agegroup,lev)])*(match(agegroup,lev)>1)
      ExHosp <- (exp(beta)-1)*exp(vbeta0)*100000
      SE.log.Exhosp=sqrt((exp(beta)/(exp(beta)-1))^2 * (vbeta+2*covbbeta0*(exp(beta)-1)/exp(beta)) + vbeta0)
            adj.est <- rbind(adj.est,data.frame(SubCat=subc[i],agegrp=agegroup,
					               logPR=beta,SE=sqrt(vbeta),logPR_ci_lower=beta-1.96*sqrt(vbeta),
					logPR_ci_upper=beta+1.96*sqrt(vbeta),ExHosp=ExHosp,ExHosp_ci_lower=ExHosp*exp(-1.96*SE.log.Exhosp),
					ExHosp_ci_upper=ExHosp*exp(1.96*SE.log.Exhosp)))
    }
    # get overall statistic
    beta0 = summary(model2.adj)$coef[1,1] 
    vbeta0=vcov(model2.adj)[1,1]
    beta = summary(model2.adj)$coef[2,1] 
    vbeta = diag(vcov(model2.adj))[2] 
    covbbeta0 = vcov(model2.adj)[1,2]
    ExHosp <- (exp(beta)-1)*exp(vbeta0)*100000
    SE.log.Exhosp=sqrt((exp(beta)/(exp(beta)-1))^2 * (vbeta+2*covbbeta0*(exp(beta)-1)/exp(beta)) + vbeta0)
    adj.est <- rbind(adj.est,data.frame(SubCat=subc[i],agegrp="Overall",
					logPR=beta,
                                        SE=sqrt(vbeta),logPR_ci_lower=beta-1.96*sqrt(vbeta),
					logPR_ci_upper=beta+1.96*sqrt(vbeta),ExHosp=ExHosp,ExHosp_ci_lower=ExHosp*exp(-1.96*SE.log.Exhosp),
					ExHosp_ci_upper=ExHosp*exp(1.96*SE.log.Exhosp)))

  }
TEvent       <- rbind(TEvent,aggregate(data.all$count,by=list(SubCat=data.all$SubCat,agegrp=data.all$agegrp),sum))
TEvent.ctl   <- rbind(TEvent.ctl,aggregate(out$count,by=list(SubCat=out$SubCat,agegrp=out$agegrp),sum))
TEvent.trt   <- rbind(TEvent.trt,aggregate(out2$count,by=list(SubCat=out2$SubCat,agegrp=out2$agegrp),sum))
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

