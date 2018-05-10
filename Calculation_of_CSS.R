

#calculation of DSS scores for single drugs and CSS scores for direct and reverse drug combos


library( dplyr )
library( tidyr )
library( plyr )
library( grid )
library( tibble )
library( data.table )
library( reshape2 )
library( readr )
library( stringr )
library( readxl )
library( ggplot2 )

#DSS function
dss<-function(ic50,slope,max,min.conc.tested,max.conc.tested,y=10,DSS.type=2,concn_scale=1e-9){
  #rdata should be in in format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration
  a=as.numeric(unname(max))
  
  b=as.numeric(unname(slope))
  d=0 # min response
  ic50 = as.numeric(unname(ic50))
  min.conc.tested = as.numeric(unname(min.conc.tested))
  max.conc.tested = as.numeric(unname(max.conc.tested))
  Min.Conc<- log10(min.conc.tested*concn_scale) #
  Max.Conc<- max.conc.tested
  x2<-log10(Max.Conc*concn_scale)  
  
  
  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){
    dss<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    dss<-0
  }
  else if(isTRUE(b==0)){
    dss<-0
  }
  else{
    if(a>100){ a<-100  }
    if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
        else if(isTRUE(x1 > x2)){x1<-x2}
      }
      else {x1<-Min.Conc}
      
      # This is a logistic function used in Dotmatics.com
      # y = d+(a-d)/(1+10^(b*(c-x)))
      #inverse function
      # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
      int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
      
      total_area<-(x2-Min.Conc)*(100-y)
      
      if(DSS.type==1){
        norm_area<-((int_y/total_area)*100)#DSS1
      }
      if(DSS.type==2){
        #	if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
      }
      if(DSS.type==3){
        #	if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3 #AUC5
      }
      if(isTRUE(norm_area < 0|norm_area > 100)){
        dss<-0
      }else{
        dss<-round(norm_area,digits=4)}
    } else {dss<-0} 
  } 
  return (dss)
}

#functions to save DSS values and other statistics
get_param<-function(tbl){
  
  
  sel_drug <- unique(tbl$Drug)
  sel_patient <- unique(tbl$Sample)
  
  mat_tbl <- data.frame( inhibition=tbl$inhibition, dose=tbl$dose, logconc = log10( tbl$dose ) )
  mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]	
  
  
  estimate_param <- tryCatch({drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drc::drmc(errorm = FALSE))}, 
                             warning=function(w){drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                             error=function(e){drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
  # (extract and name coefficients)
  coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
  coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1 
  
  # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
  coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=TRUE), max(mat_tbl$dose,na.rm=TRUE),coef_estim["IC50"])
  # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=TRUE),coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=TRUE),coef_estim["IC50"])
  # similar to previous step but now compare log10(IC50) with log(min. conc.).
  coef_estim["IC50"] <- log10(coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
  # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
  coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=TRUE),coef_estim["IC50"])
  #(Trying to fix curves that need outlier kickout)
  coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=TRUE)
  #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.	
  min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=TRUE) > 0,min(mat_tbl$inhibition,na.rm=TRUE),0)
  min_lower <- ifelse(min_lower >= 100,99,min_lower)
  #similar to previous step but for MAX
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
  #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
  max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=TRUE)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=TRUE))
  max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=TRUE))
  max_lower <- ifelse(max_lower < 0,0,max_lower)
  max_lower <- ifelse(max_lower > 100,100,max_lower)
  #(Fix upper maximum for negative slopes)
  run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
  max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
  max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
  max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
  max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
  max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
  # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen. 
  mean_inh_last = mean(mat_tbl$inhibition[4:5],na.rm=TRUE)
  if(mean_inh_last < 60) 
  {if(mean_inh_last > 25)
    coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=TRUE)
  else if(mean_inh_last < 25)
    coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=TRUE)}
  if(mean(mat_tbl$inhibition[1:3],na.rm=TRUE)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=TRUE)
  #add a bit of positive noise to MAX if it is the same as MIN. 
  if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
  
  
  #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
  nls_result_ic50 <- tryCatch({
    nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", 
        start=list(SLOPE=1,MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
        lower=list(SLOPE=0,MIN=0,MAX=max_lower,IC50=min(mat_tbl$logconc)),
        upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
        control=list(warnOnly=TRUE,minFactor = 1/2048))
  }, error = function(e) {
    
    minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                      start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                      lower=c(SLOPE=0.5, MIN=min_lower,MAX=100,	IC50=min(mat_tbl$logconc)),
                      upper=c(SLOPE=2.5, MIN=coef_estim["MIN"],MAX=100, IC50=max(mat_tbl$logconc)))
  })
  
  #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
  if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
  {
    if(mean_inh_last > 60)
      coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=TRUE)
    nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",
                           start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                           lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),
                           upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
                           control=list(warnOnly=TRUE,minFactor = 1/2048))
  }
  max_signal <- max(mat_tbl$dose,na.rm=TRUE)
  min_signal <- min(mat_tbl$dose,na.rm=TRUE)
  #############################  
  #############   Final modification & STD error
  
  #prepare final data and convert IC50 back from log scale (inverse)
  coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; 
  
  coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
  
  #(Fix ic50 for curves in wrong direction)
  coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
  
  
  #(Fix based on MAX)
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
  coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
  
  #(Fix over sensitive drugs)
  coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=TRUE),min(mat_tbl$inhibition,na.rm=TRUE))>50),min_signal,coef_ic50["IC50"])
  
  
  #Calculate the standard error scores
  sumIC50 = summary(nls_result_ic50); 
  ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]
  ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),2)
  
  
  dss_score <- dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=2)
  
  #x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=200)
  #y <- predict(nls_result_ic50, data.frame(logconc=x))
  #icpl <-  ggplot2::ggplot(mat_tbl, aes(logconc, inhibition)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
  #  geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
  #   geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8)  + theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
  #    geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T)
  #icpl <- icpl + ggtitle( paste( "Drug:",sel_drug,"Sample=", sel_patient ) )
  # plot(icpl)
  #plot_list[[i]] = icpl
  
  i <- i+1
  # icpl
  # y1<- coef_ic50["MIN"] + (coef_ic50["MAX"] - coef_ic50["MIN"])/(1 + 10^(- coef_ic50["SLOPE"]*(x - log10(coef_ic50["IC50"]  ))))
  #plot( x,y1 , type = 'l', ylim=c(-2,120))
  
  
  #points(mat_tbl$logconc,mat_tbl$inhibition, col= 'red' )
  
  res = data.frame(Drug = sel_drug, Sample = sel_patient, min_val =coef_ic50["MIN"], max_val = coef_ic50["MAX"], Slope= coef_ic50["SLOPE"], IC50 = coef_ic50["IC50"],min_conc = min_signal,  DSS = dss_score,  min_log_conc= min(mat_tbl$logconc), max_log_conc = max(mat_tbl$logconc),max_conc = max_signal  )
}
get_param1<-function(tbl){
  
  
  sel_drug <- unique(tbl$Drug)
  sel_patient <- unique(tbl$Sample)
  
  mat_tbl <- data.frame( inhibition=tbl$inhibition, dose=tbl$dose, logconc = log10( tbl$dose ) )
  mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]	
  
  
  estimate_param <- tryCatch({drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drc::drmc(errorm = FALSE))}, 
                             warning=function(w){drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                             error=function(e){drc::drm(inhibition ~ logconc, data = mat_tbl, fct = drc::L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
  # (extract and name coefficients)
  coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
  coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1 
  
  # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
  coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=TRUE), max(mat_tbl$dose,na.rm=TRUE),coef_estim["IC50"])
  # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=TRUE),coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=TRUE),coef_estim["IC50"])
  # similar to previous step but now compare log10(IC50) with log(min. conc.).
  coef_estim["IC50"] <- log10(coef_estim["IC50"])
  coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
  # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
  coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=TRUE),coef_estim["IC50"])
  #(Trying to fix curves that need outlier kickout)
  coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=TRUE)
  #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.	
  min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=TRUE) > 0,min(mat_tbl$inhibition,na.rm=TRUE),0)
  min_lower <- ifelse(min_lower >= 100,99,min_lower)
  #similar to previous step but for MAX
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
  coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
  #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
  max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=TRUE)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=TRUE))
  max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=TRUE))
  max_lower <- ifelse(max_lower < 0,0,max_lower)
  max_lower <- ifelse(max_lower > 100,100,max_lower)
  #(Fix upper maximum for negative slopes)
  run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
  max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
  max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
  max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
  max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
  max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
  # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen. 
  mean_inh_last = mean(mat_tbl$inhibition[4:5],na.rm=TRUE)
  if(mean_inh_last < 60) 
  {if(mean_inh_last > 25)
    coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=TRUE)
  else if(mean_inh_last < 25)
    coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=TRUE)}
  if(mean(mat_tbl$inhibition[1:3],na.rm=TRUE)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=TRUE)
  #add a bit of positive noise to MAX if it is the same as MIN. 
  if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
  
  
  #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
  nls_result_ic50 <- tryCatch({
    nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", 
        start=list(SLOPE=1,MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
        lower=list(SLOPE=0,MIN=0,MAX=max_lower,IC50=min(mat_tbl$logconc)),
        upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
        control=list(warnOnly=TRUE,minFactor = 1/2048))
  }, error = function(e) {
    
    minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                      start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                      lower=c(SLOPE=0.5, MIN=min_lower,MAX=100,	IC50=min(mat_tbl$logconc)),
                      upper=c(SLOPE=2.5, MIN=coef_estim["MIN"],MAX=100, IC50=max(mat_tbl$logconc)))
  })
  
  #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
  if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
  {
    if(mean_inh_last > 60)
      coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=TRUE)
    nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",
                           start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                           lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),
                           upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),
                           control=list(warnOnly=TRUE,minFactor = 1/2048))
  }
  max_signal <- max(mat_tbl$dose,na.rm=TRUE)
  min_signal <- min(mat_tbl$dose,na.rm=TRUE)
  #############################  
  #############   Final modification & STD error
  
  #prepare final data and convert IC50 back from log scale (inverse)
  coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; 
  
  coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
  
  #(Fix ic50 for curves in wrong direction)
  coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
  
  
  #(Fix based on MAX)
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
  coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
  coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
  
  #(Fix over sensitive drugs)
  coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=TRUE),min(mat_tbl$inhibition,na.rm=TRUE))>50),min_signal,coef_ic50["IC50"])
  
  
  #Calculate the standard error scores
  sumIC50 = summary(nls_result_ic50); 
  ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]
  ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),2)
  
  
  dss_score <- dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=2)
  
  #x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=200)
  #y <- predict(nls_result_ic50, data.frame(logconc=x))
  #icpl <-  ggplot2::ggplot(mat_tbl, aes(logconc, inhibition)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
  #  geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
  #   geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8)  + theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
  #    geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T)
  #icpl <- icpl + ggtitle( paste( "Drug:",sel_drug,"Sample=", sel_patient ) )
  # plot(icpl)
  #plot_list[[i]] = icpl
  
  i <- i+1
  # icpl
  # y1<- coef_ic50["MIN"] + (coef_ic50["MAX"] - coef_ic50["MIN"])/(1 + 10^(- coef_ic50["SLOPE"]*(x - log10(coef_ic50["IC50"]  ))))
  #plot( x,y1 , type = 'l', ylim=c(-2,120))
  
  
  #points(mat_tbl$logconc,mat_tbl$inhibition, col= 'red' )
  
  res = data.frame( Num = tbl$num, Drug = sel_drug, Sample = sel_patient, min_val =coef_ic50["MIN"], max_val = coef_ic50["MAX"], Slope= coef_ic50["SLOPE"], IC50 = coef_ic50["IC50"],min_conc = min_signal,  DSS = dss_score,  min_log_conc= min(mat_tbl$logconc), max_log_conc = max(mat_tbl$logconc),max_conc = max_signal  )
}

setwd( '~/Desktop/oneil drug combo' )



# read data for single drugs
single_drug <- read_excel("~/Desktop/oneil drug combo/156849_1_supp_0_w2lh45.xlsx")
single_drug[, 8:10] <- as.numeric( unlist( single_drug[, 8:10] ) )

#calculate mean response over replicates
ind <- vector()
for ( i in 1:nrow( single_drug ) ) {
  
  r <- single_drug[i,5:10]
  ind<-  which(  r != 'NA' )

  single_drug$viability[i] <- mean( as.matrix( r[ind] ) )
  
} 

single_drug <- single_drug[, c( 2, 3, 4, 13 )]
colnames( single_drug )[3] <- 'dose'


single_drug<- 
  single_drug %>%
  group_by( cell_line, drug_name, dose ) %>%
  dplyr::summarise( viability = mean( viability))

#calculate inhibition
single_drug$inhibition <- 1 - single_drug$viability

cl <- unique( single_drug$cell_line )
dr <- unique( single_drug$drug_name )






#read data for combo
combo_drug <- read_excel("~/Desktop/oneil drug combo/156849_1_supp_1_w2lrww.xlsx")
colnames( combo_drug )[c( 4,6 )] <- c( 'drugA_conc', 'drugB_conc' )
combo_drug[, 11] <- as.numeric( unlist(combo_drug[, 11]  ) )


#calculate mean response over replicates
for ( i in 1:nrow( combo_drug ) ) {
  
  r <- combo_drug[i,8:11]
  ind<-  which(  r != 'NA' )
  print(i)
  combo_drug$viability[i] <- mean( as.matrix( r[ind] ) )
  
} 

combo_drug<- 
  combo_drug %>%
  group_by( cell_line, drugA_name, drugA_conc, drugB_name, drugB_conc) %>%
  dplyr::summarise( viability = mean( viability))

#calculate inhibition
combo_drug$inhibition <- 1 - combo_drug$viability
combo_drug$Sample <- paste( combo_drug$cell_line, paste( combo_drug$drugB_name, paste( combo_drug$drugA_name,'_background_',combo_drug$drugA_conc, sep ='' ), sep= '_'), sep = '_' )
colnames( combo_drug )[c( 2,3,4,5 )] <- c( 'Background_drug', 'Background_dose', 'Drug', 'dose' )
combo_drug$inhibition <- combo_drug$inhibition*100
#save( combo_drug, file = "oneil_combo_inh.RData" )





#call DSS function for single
single_drug$inhibition <- single_drug$inhibition*100
single_drug$Sample <- paste( single_drug$cell_line, single_drug$drug_name,sep = '_' )
colnames( single_drug )[2] <- 'Drug' 
attr(single_drug, "vars") <- NULL
attr(single_drug, "drop") <- NULL

library(drc)

get_stat <- function( tbl ) {
 
  tbl <-
    tbl %>%
    dplyr::select( Drug, dose, inhibition, Sample )
 res <- get_param( tbl )
 
  
}

datalist = list()
cl <- unique( single_drug$cell_line )


for ( i in 1:length( cl ) ) {
  
  print(i)
  datalist[[i]] <-
    single_drug %>%
    filter( cell_line == cl[i] )%>%
    group_by( Drug ) %>%
    do( { get_stat(.) } ) 
  
  dss_cell_line <- dplyr::bind_rows( datalist )
  
}

detach("package:drc", unload=TRUE)

dss_cell <-dss_cell_line

#dss for single is saved in oneil_single.RData
#save( dss_cell, file = "oneil_single.RData" )




#call DSS function for combos:direct and reverse
#load('oneil_combo_inh.RData')
attr( combo_drug, "vars") <- NULL
attr( combo_drug, "drop") <- NULL
colnames( combo_drug )[2:5] <- c( 'Drug', 'dose', 'Background_drug', 'Background_dose' ) ### FOR REVERSE CASE
#colnames( combo_drug )[2:5] <- c( 'Background_drug', 'Background_dose', 'Drug', 'dose' ) #### FOR DIRECT CASE


library(drc)

get_stat1 <- function( tbl ) {
  
  tbl <-
    tbl %>%
    dplyr::select( Drug, dose, inhibition, Sample )
  
  tbl$num <- nrow( tbl )
  res <- get_param1( tbl )
  
  
}

datalist = list()
cl <- unique( combo_drug$cell_line )


for ( i in 1:length( cl ) ) {
  
  print(i)
  datalist[[i]] <-
    combo_drug %>%
    filter( cell_line == cl[i] )%>%
    group_by( Background_drug, Background_dose, Drug ) %>%
    do( { get_stat1(.) } ) 
  
  combo_cell_line <- dplyr::bind_rows( datalist )
  
}

detach("package:drc", unload=TRUE)

combo_dss_cell <-
  combo_dss_cell %>%
  distinct()
combo_dss_cell <- combo_cell_line


#save( combo_dss_cell, file = "oneil_combo.RData" ) ### DIRECT COMBO CASE
#save( combo_dss_cell, file = "oneil_reverse_combo.RData" ) ### REVERSE COMBO CASE





#Just check 
############################################################################################

#check whether IC50 values are within the range of tested concentrations in the combo experiment
#seperately for selected cell line ( too many cell lines to plot all together )



for ( i in 1:length( cl ) ) {
dd <- combo_dss_cell
dd$cl <-word(dd$Sample,1,sep = "\\_")


selected_cell_line <- cl[i]
res <-
  dd %>%
  ungroup() %>%
  dplyr::select( cl, Background_drug, Background_dose ) %>%
  distinct() %>%
  filter( cl == selected_cell_line )

res$cond <- 'tested_in_combo'

f <- dss_cell_line
f$cl <-word( f$Sample,1,sep = "\\_")
f <- 
  f %>%
  dplyr::select( cl,Drug, IC50 )%>%
  filter( cl == selected_cell_line )

colnames( f )[c(2, 3)] <-  c('Background_drug','Background_dose' )
f$cond <- 'IC50'
f <- as.data.frame(f)

res <- rbind( res, f )



icpl <-   ggplot2::ggplot( res, aes( Background_drug, Background_dose, color = cond ) ) +  geom_point( size = 1.8, alpha = 0.5 ) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
icpl <- icpl + ggtitle( paste0('Cell line:', cl[i] ) )
plot(icpl)


}

# it seems that for every cell line the IC50 for every background drug is within the range of the combo tested concentrations


