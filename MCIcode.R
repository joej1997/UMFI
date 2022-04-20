library(rje)
library(ranger)
library(xgboost)
library(iml)
library(reshape2)
library(ggplot2)
library(gridExtra)
rm(list = ls())
gc()

######################################################################################3
##############################   start experiments   #################################
#######################################################################################

#import data
camels_clim <- read.csv("C:/Users/joeja/Desktop/researchMasters/Snow_Paper/camels_attributes_v2.0/camels_attributes_v2.0/camels_clim.txt", sep=";")
camels_geol <- read.csv("C:/Users/joeja/Desktop/researchMasters/Snow_Paper/camels_attributes_v2.0/camels_attributes_v2.0/camels_geol.txt", sep=";")
camels_soil <- read.csv("C:/Users/joeja/Desktop/researchMasters/Snow_Paper/camels_attributes_v2.0/camels_attributes_v2.0/camels_soil.txt", sep=";")
camels_topo <- read.csv("C:/Users/joeja/Desktop/researchMasters/Snow_Paper/camels_attributes_v2.0/camels_attributes_v2.0/camels_topo.txt", sep=";")
camels_vege <- read.csv("C:/Users/joeja/Desktop/researchMasters/Snow_Paper/camels_attributes_v2.0/camels_attributes_v2.0/camels_vege.txt", sep=";")
camels_hydro <- read.csv("C:/Users/joeja/Desktop/researchMasters/Snow_Paper/camels_attributes_v2.0/camels_attributes_v2.0/camels_hydro.txt", sep=";")

# combine dataset
dat<-camels_clim
dat<-merge(dat,camels_geol,by="gauge_id",all.x=T)
dat<-merge(dat,camels_soil,by="gauge_id",all.x=T)
dat<-merge(dat,camels_topo,by="gauge_id",all.x=T)
dat<-merge(dat,camels_vege,by="gauge_id",all.x=T)
dat<-merge(dat,camels_hydro,by="gauge_id",all.x=T)


rm(camels_clim,camels_geol,camels_soil,camels_topo,camels_vege,camels_hydro)

# remove character variables
dat$gauge_lat<-NULL
dat$gauge_lon<-NULL

dat$high_prec_timing<-NULL
dat$geol_1st_class<-NULL
dat$geol_2nd_class<-NULL
dat$dom_land_cover<-NULL
dat$dom_land_cover_frac<-NULL
dat$glim_1st_class_frac<-NULL
dat$glim_2nd_class_frac<-NULL
dat$dom_land_cover_frac<-NULL
dat$area_geospa_fabric<-NULL
dat$low_prec_timing<-NULL
dat$gauge_id<-NULL
dat$area_gages2<-NULL
dat$water_frac<-NULL
dat$organic_frac<-NULL
dat$other_frac<-NULL
dat$zero_q_freq<-NULL
dat$runoff_ratio<-NULL
dat$stream_elas<-NULL
dat$high_q_freq<-NULL
dat$high_q_dur<-NULL
dat$low_q_dur<-NULL
dat$low_q_freq<-NULL

dat<-dat[complete.cases(dat),]

responses<-dat[,(ncol(dat)-5):ncol(dat)]
dat[,(ncol(dat)-5):ncol(dat)]<-NULL

#import BRCA data instead
Bigdat<- read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCA.csv")
Bigdat$Sample.ID<-NULL
responses<-as.factor(Bigdat$BRCA_Subtype_PAM50)
Bigdat$BRCA_Subtype_PAM50<-NULL

set.seed(1)
dat<-Bigdat[,sample(1:ncol(Bigdat),15)]


#modify data functions
modifty_linreg<-function(dat,protect){
  #remove dependedence via linear regression
  modifiedDAT<-dat
  tomodify<-setdiff(1:ncol(dat),protect)
  for(i in tomodify){
    mod<-lm(dat[,i]~dat[,protect])
    if(summary(mod)$coefficients[2,4]<0.01) modifiedDAT[,i]<- mod$residuals
    if(var(modifiedDAT[,i])==0) modifiedDAT[,i]<-rnorm(nrow(modifiedDAT))
  }
  modifiedDAT
}


#remove dependence via pairwise optimal transport
modifty_otpw_quantiles_lin<-function(dat,protect){
  modifiedDAT<-dat #the new dataframe that will be returned
  tomodify<-setdiff(1:ncol(dat),protect) #the columns in dat to modify
  n<-nrow(dat) #number of rows in dat
  z=dat[,protect] #the protected attribute
  n_quan=ceiling(nrow(dat)/100) #number of quantiles to use (20 points per regression)
  quans<-(seq(from=0,to=1,length.out = n_quan)) #quantiles of interest
  quans<-quantile(z,quans) #quantiles of z
  #loop through each feature we need to modify
  for(j in tomodify){
    x=dat[,j] #feature we will modifty
    newx<-x
    orderedCONDF<-sort(x) #sorted x
    for(quan in 2:n_quan){
      cur_obs<- (z<=quans[quan] & z>=quans[quan-1])
      x_curquan=x[cur_obs]
      z_curquan=z[cur_obs]
      if(sd(x_curquan)<1e-6) x_curquan<-x_curquan+rnorm(length(x_curquan),sd=sd(x)/length(x))
      if(sd(z_curquan)<1e-6) z_curquan<-z_curquan+rnorm(length(z_curquan),sd=sd(z)/length(z))
      condF<-x_curquan
      mod<-lm(x_curquan~z_curquan)
      rv<-x_curquan
      if(summary(mod)$coefficients[2,4]<1) rv<-as.numeric(mod$residuals)
      #condF<-pnorm(rv/var(rv))
      for(k in 1:length(x_curquan)){
        condF[k]<- mean(rv<= rv[k])
        #newx[cur_obs][k]<-orderedCONDF[max(sum((1:n)/n<=condF[k]),1)]
        newx[cur_obs][k]<-condF[k]
      }
    }
    modifiedDAT[,j]<-newx
  }
  modifiedDAT
}



calc_dependence<-function(dat,tt,method,mod_meth){
  
  to_return<-rep(0,ncol(dat))
  
  for(j in 1:length(to_return)){
    if(mod_meth=="nothing"){
      dat2<-dat
    }
    
    if(mod_meth=="lin"){
      dat2<-modifty_linreg(dat,j)
    }
    
    if(mod_meth=="nonlin"){
      dat2<-modifty_nonlinreg(dat,j)
    }
    
    if(mod_meth=="ot_pw"){
      dat2<-modifty_otpw(dat,j)
    }
    
    if(mod_meth=="ot_pw_quan"){
      dat2<-modifty_otpw_quantiles_lin(dat,j)
    }
    
    train<-dat2[tt<0.8,]
    test<-dat2[tt>0.8,]
    
    if(method=="lin"){
      y<-train[,j]
      x<-train[,-c(j)]
      mod<-lm(y~.,data=x)
      preds<-as.numeric(predict(mod,newdata = test[,-c(j)]))
      to_return[j]<-max(1-sum((test[,j]-preds)^2)/sum((test[,j]-mean(test[,j]))^2),0)
    }
    
    if(method=="rf"){
      rfmod<-ranger(y=train[,j],x=train[,-c(j)])
      preds<-predict(rfmod,test[,-c(j)])$predictions
      to_return[j]<-max(1-sum((test[,j]-preds)^2)/sum((test[,j]-mean(test[,j]))^2),0)
    }
  }
  to_return
}

#get plot of dependence of protected attribute with rest of dataset
tt<-runif(nrow(dat))

#lin_none<-calc_dependence(dat=dat,tt=tt,method="lin",mod_meth = "nothing")
rf_none<-calc_dependence(dat=dat,tt=tt,method="rf",mod_meth = "nothing")
rf_resid<-calc_dependence(dat=dat,tt=tt,method="rf",mod_meth = "lin")
#lin_resid<-calc_dependence(dat=dat,tt=tt,method="lin",mod_meth = "ot_pw_quan")
rf_ot_pw<-calc_dependence(dat=dat,tt=tt,method="rf",mod_meth = "ot_pw_quan")


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/MutualInfoR2.jpeg"),width = 12, height = 9, units = 'in',res = 600)
par(mar = c(4, 6, 1, 1))
plot(rf_none,col="black",xlab="Feature",ylab = "Mutual Information (R^2)",type="l",ylim=c(0,0.75),lwd=5, cex.axis=1.6,cex.lab=1.8)
lines(rf_ot_pw,col="blue",lwd=5)
lines(rf_resid,col="red",lwd=5)
legend("topleft",col=c("black","red", "blue"),lty = 1,lwd=5,legend = c("Raw Data","Linear Regression","Pairwise Optimal Transport"),cex=1.8)
dev.off()



#test how much we are modifying original dataset
library(minerva)
library(corrplot)

modification_ot<-matrix(0,nrow = ncol(dat),ncol = ncol(dat))
modification_lin<-matrix(0,nrow = ncol(dat),ncol = ncol(dat))
for(i in 1:nrow(modification_ot)){
  dat2<-modifty_otpw_quantiles_lin(dat,i)
  dat3<-modifty_linreg(dat,i)
  for(j in 1:ncol(modification_ot)){
    modification_ot[i,j]<-mine(dat[,j],dat2[,j],est = "mic_e")[[1]]
    modification_lin[i,j]<-mine(dat[,j],dat3[,j],est = "mic_e")[[1]]
  }
}

modification_ot<-round(modification_ot,2)
modification_lin<-round(modification_lin,2)

colnames(modification_ot)<-colnames(dat)
colnames(modification_lin)<-colnames(dat)
rownames(modification_ot)<-colnames(dat)
rownames(modification_lin)<-colnames(dat)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/MutualInfoMIC_both.jpeg"),width = 12, height = 6, units = 'in',res = 600)
par(mar = c(1, 1, 1, 1), mfrow=c(1,2))
corrplot(modification_ot, method="color", col=col(200), order="original", is.corr = F,
         addCoef.col = "black", cl.cex = 1.2, cl.length = 3,cl.lim = c(0, 1), # Add coefficient of correlation
         tl.col="black", tl.srt=45, number.cex = 0.8,number.digits = 2,tl.cex = 0.8,addgrid.col = "grey" #Text label color and rotation
)
mtext("Optimal Transport", at=8.5, line=-1.5, cex=1.5)

corrplot(modification_lin, method="color", col=col(200), order="original", is.corr = F,
         addCoef.col = "black", cl.cex = 1.2, cl.length = 3,cl.lim = c(0, 1), # Add coefficient of correlation
         tl.col="black", tl.srt=45, number.cex = 0.8,number.digits = 2,tl.cex = 0.8,addgrid.col = "grey" #Text label color and rotation
)
mtext("Linear Regression", at=8.5, line=-1.5, cex=1.5)


dev.off()

# jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/MutualInfoMIC_lin.jpeg"),width = 12, height = 12, units = 'in',res = 600)
# corrplot(modification_lin, method="color", col=col(200), order="original", is.corr = F,
#          addCoef.col = "black", cl.cex = 1.5, cl.length = 3,cl.lim = c(0, 1), # Add coefficient of correlation
#          tl.col="black", tl.srt=45, number.cex = 1.5,number.digits = 2,tl.cex = 1.5,addgrid.col = "grey" #Text label color and rotation
# )
# dev.off()



#################################################################################################################
###############################   feature importance tests
#################################################################################################################

UMFI<- function(X,y,mod_meth){
  fi<-rep(0,ncol(X))
  for(i in 1:length(fi)){
    if(mod_meth=="otpw") newX<-modifty_otpw_quantiles_lin(X,i)
    if(mod_meth=="lin") newX<-modifty_linreg(X,i)
    #newX[,apply(newX,2,sd)<1e-6]<-matrix(rnorm(nrow(X)*sum(apply(newX,2,sd)<1e-6)),nrow = nrow(X))
    #newX<-cbind(newX,rnorm(nrow(newX)),rnorm(nrow(newX)))
    rfwith<-ranger(x=newX,y=y,num.trees = 100)
    rfwithout<-ranger(x=newX[,-c(i)],y=y,num.trees = 100)
    if(is.numeric(y)) fi[i]<-max(rfwith$r.squared,0)-max(rfwithout$r.squared,0)
    if(is.factor(y)) fi[i]<- max(1-rfwith$prediction.error,0.5)-max(1-rfwithout$prediction.error,0.5)
  }
  fi[fi<0]<-0
  fi
}

#Variable importance functions
MCI<-function(X,y){
  colvec<-1:ncol(X)
  if(ncol(X)<=15) CompleteSet<-powerSet(colvec)
  if(ncol(X)>15) CompleteSet<-powerSet(colvec,m=3)
  CompleteSetErrors<-rep(0,length(CompleteSet))
  
  for(e in 1:length(CompleteSetErrors)){
    if(length(CompleteSet[[e]])>0){
      rfmod<-ranger(y=y,x=as.data.frame(X[,CompleteSet[[e]]]),num.trees = 100)
      if(is.numeric(y)) CompleteSetErrors[e]<-rfmod$r.squared
      if(is.factor(y)) CompleteSetErrors[e]<- 1- rfmod$prediction.error
    }
  }
  
  if(is.numeric(y)) CompleteSetErrors[CompleteSetErrors<0]<-0
  if(is.factor(y)) CompleteSetErrors[CompleteSetErrors<0.5]<- 0.5
  
  OUTPUT<-rep(0,ncol(X))
  for(j in 1:ncol(X)){
    jsHERE<-unlist(lapply(CompleteSet, is.element,el=j))
    jSET<-CompleteSet[jsHERE]
    
    NOjSET<-lapply(jSET, setdiff,y=j)
    NOjSET<-intersect(NOjSET,CompleteSet)
    jSET<-lapply(NOjSET, c,j)
    jSET<-lapply(jSET, sort)
    
    charlistjSET<-unlist(lapply(jSET,paste,collapse=" "))
    charlistNOjSET<-unlist(lapply(NOjSET,paste,collapse=" "))
    charlistCompleteSet<-unlist(lapply(CompleteSet,paste,collapse=" "))
    errorWITH<-CompleteSetErrors[order(match(charlistCompleteSet, charlistjSET),na.last = NA)]
    errorWITHOUT<-CompleteSetErrors[order(match(charlistCompleteSet, charlistNOjSET),na.last = NA)]
    
    OUTPUT[j]<-max(errorWITH- errorWITHOUT)
  }
  OUTPUT
}



getResultsPlot<-function(simTest,nobs,niter,nX){
  
  Imp<-list(MCI=matrix(0,nrow = niter,ncol = nX),UMFI_LI=matrix(0,nrow = niter,ncol=nX),UMFI_OT=matrix(0,nrow = niter,ncol=nX))
  for(i in 1:niter){
    if(simTest=="Correlated Interaction"){
      A<-rnorm(nobs,mean = 0,sd=1)
      B<-rnorm(nobs,mean = 0,sd=1)
      C<-rnorm(nobs,mean = 0,sd=1)
      D<-rnorm(nobs,mean = 0,sd=1)
      E<-rnorm(nobs,mean = 0,sd=1)
      G<-rnorm(nobs,mean = 0,sd=1)
      Boston=data.frame(x1=A+B,x2=B+C,x3=D+E,x4=E+G)
      Boston$y<-Boston$x1+Boston$x2+sign(Boston$x1*Boston$x2)+Boston$x3+Boston$x4
    }
    
    if(simTest=="Nonlinearity"){
      Boston=data.frame(x1=runif(nobs,-3,3),x2=runif(nobs,-3,3),x3=runif(nobs,-3,3),x4=runif(nobs,-3,3))
      X1<-.3*(Boston$x1^3-9*Boston$x1)
      X2<-.47*exp(Boston$x2)
      X3<-rep(0,nobs)
      X3[Boston$x3<1 & Boston$x3> -1]<-4.75
      Boston$y<-X1+X2+X3
    }
    
    if(simTest=="Interaction"){
      Boston=data.frame(x1=rnorm(nobs,mean=0,sd=1),x2=rnorm(nobs,mean=0,sd=1),x3=rnorm(nobs,mean=0,sd=1),x4=rnorm(nobs,mean=0,sd=1))
      Boston$y<-sign(Boston$x1*Boston$x2)*rexp(nobs,1/sqrt(2))
    }
    
    if(simTest=="Noise"){
      Boston=data.frame(x1=rnorm(nobs,mean=0,sd=1),x2=rnorm(nobs,mean=0,sd=1),x3=rnorm(nobs,mean=0,sd=1),x4=rnorm(nobs,mean=0,sd=1))
      Boston$y<-2*Boston$x1+Boston$x2+5*rnorm(nobs)
    }
    
    if(simTest=="Correlation"){
      Boston=data.frame(x1=rnorm(nobs,mean=0,sd=1),x2=rnorm(nobs,mean=0,sd=1))
      Boston$x3<-Boston$x1 +rnorm(nobs,mean=0,sd=0.2)
      Boston$x4<-rnorm(nobs,mean=0,sd=1)
      Boston$y<-Boston$x1+Boston$x2
    }
    
    if(simTest=="Dependence"){
      Boston=data.frame(x1=rnorm(nobs,mean=0,sd=5),x2=rnorm(nobs,mean=0,sd=5))
      Boston$x3<-abs(Boston$x1)
      Boston$x4<-abs(Boston$x1)
      Boston$y<-Boston$x1+Boston$x2
    }
    
    Imp$MCI[i,]<-MCI(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y)
    
    Imp$UMFI_LI[i,]<-UMFI(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y,mod_meth = "lin")
    Imp$UMFI_OT[i,]<-UMFI(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y,mod_meth = "otpw")
  }
  
  
  Imp2<-Imp
  for(i in 1:length(Imp)){
    curdat<-Imp[[i]]
    Imp2[[i]]<-t(apply(curdat,1,function(x) 100*x/sum(x)))
  }
  
  df <- lapply(names(Imp2),function(x)cbind(name=x,as.data.frame(Imp2[[x]])))
  df <- melt(do.call(rbind,df),id="name")
  levels(df$variable)<-colnames(Boston[,1:(ncol(Boston)-1)])
  df
}

simtests<-c("Correlated Interaction","Interaction","Correlation")

for(i in 1:length(simtests)){
  Composed<-getResultsPlot(simTest=simtests[i], nobs=500,niter=100,nX=4)
  
  jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/",simtests[i],".jpg"),width = 7, height = 5, units = 'in',res = 600)
  print(ggplot(Composed, aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + ggtitle(simtests[i])+ xlab("Variable")
        + ylab("Variable Importance (%)") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)))
  dev.off()
  print(i)
}


#downsides plot
# Composed<-getResultsPlot(simTest="Correlated Interaction", nobs=500,niter=100,nX=4)
# jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/downsides.jpg"),width = 7, height = 5, units = 'in',res = 600)
# print(ggplot(Composed[Composed$name=="MCI",], aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + 
#         xlab("Variable") + ylab("Variable Importance (%)") + 
#         theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")))
# dev.off()
# 
# 
# what<-UMFI(X=dat,y=responses$q_mean,mod_meth = "otpw")
# 
# what1<-UMFI(X=dat,y=responses$q_mean,mod_meth = "lin")


#calculate computational time of MCI vs UMCI
tests=5:15
time_UMFI<-rep(0,length(tests))
time_MCI<-rep(0,length(tests))
for(i in 1:length(tests)){
  start_time <- Sys.time()
  fi<-MCI(y=responses,X=dat[,1:(tests[i])])
  end_time <-Sys.time()
  time_MCI[i]<-as.numeric(end_time- start_time,units="mins")
  
  start_time <- Sys.time()
  fi<-UMFI(y=responses,X=dat[,1:(tests[i])],mod_meth="otpw")
  end_time <-Sys.time()
  time_UMFI[i]<-as.numeric(end_time- start_time,units="mins")
  print(tests[i])
}


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/time_complexity.jpg"),width = 7, height = 5, units = 'in',res = 600)
par(mar = c(4, 6, 4, 4))
plot(tests,time_MCI,type="l",col="red",xlab="# of Features", ylab="Time Complexity (Minutes)",lwd=5, cex.axis=1.6,cex.lab=1.8)
lines(tests,time_UMFI,col="blue",lwd=5)
legend("topleft",col=c("red", "blue"),lty = 1,lwd=5,legend = c("MCI","UMFI"),cex=1.8)
dev.off()


#calculate importance for the BRCA dataset
BRCA_otpw<-UMFI(X=Bigdat,y=responses,mod_meth="otpw")
BRCA_lin<-UMFI(X=Bigdat,y=responses,mod_meth="lin")
BRCA_MCI<-MCI(X=Bigdat,y=responses)

getcols<-function(nam){
  importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
  colnum<-rep(1,length(nam))
  colnum[nam %in% importantGenes]<-2
  colcol<-c("grey","darkblue")[colnum]
  colcol
}

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCA.jpg"),width = 9, height = 9, units = 'in',res = 600)
par(mar = c(6, 5, 3, 3), mfrow=c(3,1))
barplot(BRCA_MCI[order(- BRCA_MCI)], ylab = "MCI Importance", cex.axis=1.1,cex.lab=1.5,main="(a)",
        names.arg = colnames(Bigdat)[order(- BRCA_MCI)], col = getcols(colnames(Bigdat)[order(- BRCA_MCI)]),las=2)
legend("topright",col=c("darkblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

barplot(BRCA_lin[order(- BRCA_lin)], ylab = "UMFI_LR Importance",cex.axis=1.1,cex.lab=1.5,main="(b)",
        names.arg = colnames(Bigdat)[order(- BRCA_lin)], col = getcols(colnames(Bigdat)[order(- BRCA_lin)]),las=2)
legend("topright",col=c("darkblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

barplot(BRCA_otpw[order(- BRCA_otpw)], ylab = "UMFI_OT Importance", cex.axis=1.1,cex.lab=1.5, main="(c)",
        names.arg = colnames(Bigdat)[order(- BRCA_otpw)], col = getcols(colnames(Bigdat)[order(- BRCA_otpw)]),las=2)
legend("topright",col=c("darkblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)
dev.off()


#calculate importance for qmean and baseflow index
qmean_Imp<-UMFI(X=dat,y=responses$q_mean,mod_meth="otpw")
baseflow_Imp<-UMFI(X=dat,y=responses$baseflow_index,mod_meth="otpw")
q95_Imp<-UMFI(X=dat,y=responses$q95,mod_meth="otpw")


####################################################################################################################################3
######################################33     example code    ##########################################################################
#####################################################################################################################################

#test time of UMFI with more features
tests=5:50
time_UMFI<-rep(0,length(tests))
for(i in 1:length(tests)){
  start_time <- Sys.time()
  fi<-UMFI(y=responses,X=Bigdat[,1:(tests[i])],mod_meth="otpw")
  end_time <-Sys.time()
  time_UMFI[i]<-as.numeric(end_time- start_time,units="mins")
  print(tests[i])
}

#test ranger with classification
rfmod<-ranger(x=Bigdat,y=responses,importance = "permutation")
1-rfmod$prediction.error
rfmod$variable.importance
rfmod<-ranger(x=Bigdat[,-3],y=responses,importance = "permutation")
1-rfmod$prediction.error


#run Johndrow example 
z<-sample(c(0,1),replace = T,1000)
x<-z
x[z==0]<-rnorm(sum(z==0),sd=1,mean=4)
x[z==1]<-rnorm(sum(z==1),sd=1,mean=4+1)

newx<-rep(0,1000)
mod<-lm(x~z)
rv<-as.numeric(mod$residuals)
condF<-rv
for(i in 1:length(rv)){
  condF[i]<- mean(rv<= rv[i])
}
newx<-rep(0,1000)
orderedCONDF<-sort(x)
for(i in 1:length(x)){
  newx[i]<-orderedCONDF[sum((1:length(x))/length(x)<=condF[i])]
}
cor(newx,x)
cor(newx,z)
cor(x,z)
plot(newx,z)


#run johndrow example with quantiles of z
rm(list = ls())
z<-sample(c(0,1),replace = T,1000)
x<-z
x[z==0]<-rnorm(sum(z==0),sd=1,mean=4)
x[z==1]<-rnorm(sum(z==1),sd=1,mean=4+1)

#divide into quantiles
quans<-(seq(from=0,to=1,by=.1))
n_quan=length(quans)
x_quantiles<-list()
for(quan in 2:n_quan){
  x_quantiles[[quan-1]]<-x[z<=quans[quan] & z>=quans[quan-1]]
}
condF_quantiles=x_quantiles
newx_quantiles=x_quantiles

for(quan in 2:n_quan){
  rv=x_quantiles[[quan-1]]
  for(i in 1:length(rv)){
    condF_quantiles[[quan-1]][i]<- mean(rv<= rv[i])
  }
}

newx<-x
#put quantiles back into vector form
for(quan in 2:n_quan){
  newx[z<=quans[quan] & z>=quans[quan-1]]<-condF_quantiles[[quan-1]]
}

newx2<-rep(0,1000)
orderedCONDF<-sort(x)
for(i in 1:length(x)){
  newx2[i]<-orderedCONDF[sum((1:length(x))/length(x)<=newx[i])]
}

cor(newx2,x)
cor(newx2,z)
cor(x,z)
plot(newx2,z)




#run nonlinear example
x<-dat$carbonate_rocks_frac
z<-dat$p_mean


n<-length(x) #number of rows in dat
n_quan=12 #number of quantiles to use
quans<-(seq(from=0,to=1,length.out = n_quan)) #quantiles of interest
quans<-quantile(z,quans) #quantiles of z

newx<-x
orderedCONDF<-sort(x) #sorted x
for(quan in 2:n_quan){
  cur_obs<- (z<=quans[quan] & z>=quans[quan-1])
  x_curquan=x[cur_obs]
  z_curquan=z[cur_obs]
  if(sd(x_curquan)==0) x_curquan<-rnorm(length(x_curquan))
  condF<-x_curquan
  
  mod<-lm(x_curquan~z_curquan)
  rv<-mod$residuals
  for(k in 1:length(x_curquan)){
    condF[k]<- mean(rv<= mod$residuals[k])
    #newx[cur_obs][k]<-orderedCONDF[sum((1:n)/n<=condF[k])]
    newx[cur_obs][k]<-condF[k]
  }
}

cor(newx2,x)
cor(newx2,z)
cor(x,z)
plot(newx2,z)
plot(x,z)




#checking feature importance after dependence removal
dat2<-modifty_otpw_quantiles_lin(dat,1)
mod<-ranger(p_mean~.,data = dat2,importance="permutation")
mod$variable.importance
mod$r.squared

#what if some features are kept the same
dat2<-modifty_otpw_quantiles_lin(dat,1)
dat2$frac_snow<-dat$frac_snow
dat2$carbonate_rocks_frac<-dat$carbonate_rocks_frac
mod<-ranger(y=dat2$p_mean,x=dat2[,2:ncol(dat2)],importance="permutation")
mod$variable.importance
mod$r.squared





dat2<-modifty_otpw_quantiles(dat,1)
mod<-ranger(y=dat2$p_mean,x=dat2[,2:ncol(dat2)],importance="permutation")
mod$variable.importance
mod$r.squared

mod<-ranger(y=dat$p_mean,x=dat[,2:ncol(dat)],importance="permutation")
mod$variable.importance
mod$r.squared



#testing optimal transport
x<-rnorm(1000,sd=5)
z<-runif(1000,min=0,max=1)+3*sin(x)
newx<-rep(0,1000)

mod<-lm(x~z)

for(i in 1:length(x)){
  condF<-pnorm((x[i]-mod$fitted.values[i])/var(mod$residuals))
  newx[i]<-qnorm(condF,mean = mean(x),sd=sd(x))
}

plot(newx,z)
plot(newx,x)


#testing conditional feature importance
library("party")
library("permimp")
set.seed(542863)
airq <- subset(airquality, !(is.na(Ozone) | is.na(Solar.R)))
airq$Wind2<-airq$Wind
airq$Wind3<-airq$Wind
cfAirq50 <- cforest(Ozone ~ ., data = airq,
                    control = cforest_unbiased(mtry = 2, ntree = 50,
                                               minbucket = 5, 
                                               minsplit = 10))

CPI_permimp <- permimp(cfAirq50, conditional = TRUE, progressBar = FALSE)
CPI_varimp <- varimp(cfAirq50, conditional = TRUE)

CPI_permimp
CPI_varimp



#testing optimality of linear regression
library(energy)
nobs=50000
A<-rnorm(nobs,mean = 0,sd=1)
B<-rnorm(nobs,mean = 0,sd=1)
C<-rnorm(nobs,mean = 0,sd=1)
x1=A+B
x2=B+C

mod<-lm(x1~x2)
resid<-mod$residuals
summary(mod)



data <- data.frame(x1 = x1,x2 = x2)
mvnorm.etest(data, R=100)





