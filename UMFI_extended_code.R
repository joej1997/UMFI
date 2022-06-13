#clear all
rm(list = ls())
gc()


#import packages
library(ranger)
library(caret)
library(rje)
library(iml)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doParallel)
library(minerva)
library(corrplot)

#Import data

Bigdat<- read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCA.csv")
Bigdat$Sample.ID<-NULL
responses<-as.factor(Bigdat$BRCA_Subtype_PAM50)
Bigdat$BRCA_Subtype_PAM50<-NULL

rfmod<-ranger(y=responses,x=Bigdat,num.trees = 100) # accuracy of random forest
1-rfmod$prediction.error

#set.seed(1)
#dat<-Bigdat[,sample(1:ncol(Bigdat),15)] #get subset for analysis of linear regression vs ot


#modify data functions
modifty_linreg<-function(dat,protect){
  #remove dependedence via linear regression
  modifiedDAT<-dat
  tomodify<-setdiff(1:ncol(dat),protect)
  if(sd(dat[,protect])<1e-6) dat[,protect]<-dat[,protect]+rnorm(nrow(dat))
  for(i in tomodify){
    if(sd(dat[,i])<1e-6) dat[,i]<-dat[,i]+rnorm(nrow(dat))
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
  z=dat[,protect] #the protected attribute
  n_quan=ceiling(nrow(dat)/150) #number of quantiles to use (20 points per regression)
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
      mod<-lm(x_curquan~z_curquan)
      rv<-as.numeric(mod$residuals)
      condF<-rank(rv)/length(rv)
      newx[cur_obs]<-condF
    }
    modifiedDAT[,j]<-newx
  }
  modifiedDAT
}


calc_dependence<-function(dat,method,mod_meth){
  
  to_return<-rep(0,ncol(dat))
  
  for(j in 1:length(to_return)){
    if(mod_meth=="nothing"){
      dat2<-dat
    }
    
    if(mod_meth=="lin"){
      dat2<-modifty_linreg(dat,j)
    }
    
    if(mod_meth=="ot_pw_quan"){
      dat2<-modifty_otpw_quantiles_lin(dat,j)
    }
    
    rfmod<-ranger(y=dat2[,j],x=dat2[,-c(j)])
    to_return[j]<-max(rfmod$r.squared,0)
  }
  to_return
}

#get plot of dependence of protected attribute with rest of dataset
rf_none<-calc_dependence(dat=Bigdat,mod_meth = "nothing")
rf_resid<-calc_dependence(dat=Bigdat,mod_meth = "lin")
rf_ot_pw<-calc_dependence(dat=Bigdat,mod_meth = "ot_pw_quan")


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/MutualInfoR2.jpeg"),width = 12, height = 9, units = 'in',res = 600)
par(mar = c(4, 6, 1, 1))
plot(rf_none,col="black",xlab="Feature",ylab = "Mutual Information (R^2)",type="l",ylim=c(0,0.8),lwd=5, cex.axis=1.6,cex.lab=1.8)
lines(rf_ot_pw,col="blue",lwd=5)
lines(rf_resid,col="red",lwd=5)
legend("topleft",col=c("black","red", "blue"),lty = 1,lwd=5,legend = c("Raw Data","Linear Regression","Pairwise Optimal Transport"),cex=1.8)
dev.off()



#test how much we are modifying original dataset
modification_ot<-matrix(0,nrow = ncol(Bigdat),ncol = ncol(Bigdat))
modification_lin<-matrix(0,nrow = ncol(Bigdat),ncol = ncol(Bigdat))
for(i in 1:nrow(modification_ot)){
  dat2<-modifty_otpw_quantiles_lin(Bigdat,i)
  dat3<-modifty_linreg(Bigdat,i)
  for(j in 1:ncol(modification_ot)){
    modification_ot[i,j]<-mine(Bigdat[,j],dat2[,j],est = "mic_e")[[1]]
    modification_lin[i,j]<-mine(Bigdat[,j],dat3[,j],est = "mic_e")[[1]]
  }
}

modification_ot<-round(modification_ot,2)
modification_lin<-round(modification_lin,2)

colnames(modification_ot)<-colnames(Bigdat)
colnames(modification_lin)<-colnames(Bigdat)
rownames(modification_ot)<-colnames(Bigdat)
rownames(modification_lin)<-colnames(Bigdat)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/MutualInfoMIC_both.jpeg"),width = 12, height = 6, units = 'in',res = 600)
par(mar = c(1, 1, 1, 1), mfrow=c(1,2))
corrplot(modification_ot[1:15,1:15], method="color", col=col(200), order="original", is.corr = F,
         addCoef.col = "black", cl.cex = 1.2, cl.length = 3,cl.lim = c(0, 1), # Add coefficient of correlation
         tl.col="black", tl.srt=45, number.cex = 0.8,number.digits = 2,tl.cex = 0.8,addgrid.col = "grey" #Text label color and rotation
)
mtext("Optimal Transport", at=8.5, line=-1.5, cex=1.5)

corrplot(modification_lin[1:15,1:15], method="color", col=col(200), order="original", is.corr = F,
         addCoef.col = "black", cl.cex = 1.2, cl.length = 3,cl.lim = c(0, 1), # Add coefficient of correlation
         tl.col="black", tl.srt=45, number.cex = 0.8,number.digits = 2,tl.cex = 0.8,addgrid.col = "grey" #Text label color and rotation
)
mtext("Linear Regression", at=8.5, line=-1.5, cex=1.5)


dev.off()



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


UMFI_par<- function(X,y,mod_meth){
  fi<-foreach(i=1:ncol(X),  .inorder = FALSE, .export = c("modifty_otpw_quantiles_lin","modifty_linreg"),
              .packages = c("ranger", "doParallel"),.combine = 'c')%dopar%{
                if(mod_meth=="otpw") newX<-modifty_otpw_quantiles_lin(X,i)
                if(mod_meth=="lin") newX<-modifty_linreg(X,i)
                #newX[,apply(newX,2,sd)<1e-6]<-matrix(rnorm(nrow(X)*sum(apply(newX,2,sd)<1e-6)),nrow = nrow(X))
                #newX<-cbind(newX,rnorm(nrow(newX)),rnorm(nrow(newX)))
                rfwith<-ranger(x=newX,y=y,num.trees = 100)
                rfwithout<-ranger(x=newX[,-c(i)],y=y,num.trees = 100)
                if(is.numeric(y)) return(max(rfwith$r.squared,0)-max(rfwithout$r.squared,0))
                if(is.factor(y)) return(max(1-rfwith$prediction.error,0.5)-max(1-rfwithout$prediction.error,0.5))
              }
  fi[fi<0]<-0
  fi
}

#Variable importance functions
MCI<-function(X,y,k){
  colvec<-1:ncol(X)
  CompleteSet<-powerSet(colvec,m=k)
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


MCI_par<-function(X,y,k){
  colvec<-1:ncol(X)
  CompleteSet<-powerSet(colvec,m=k)
  
  CompleteSetErrors<-foreach(e=1:length(CompleteSet),  .inorder = FALSE,
                             .packages = c("ranger", "doParallel"),.combine = 'c')%dopar%{
                               if(length(CompleteSet[[e]])>0){
                                 rfmod<-ranger(y=y,x=as.data.frame(X[,CompleteSet[[e]]]),num.trees = 100)
                                 if(is.numeric(y)) return(rfmod$r.squared)
                                 if(is.factor(y)) return(1- rfmod$prediction.error)
                               }
                             }
  
  CompleteSetErrors<-c(0,CompleteSetErrors) #add accuracy for no features
  
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
  
  Imp<-list(MCI=matrix(0,nrow = niter,ncol = nX),UMFI_LR=matrix(0,nrow = niter,ncol=nX),UMFI_OT=matrix(0,nrow = niter,ncol=nX))
  for(i in 1:niter){
    if(simTest=="Correlated_Interaction"){
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
    
    Imp$MCI[i,]<-MCI_par(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y)
    
    Imp$UMFI_LR[i,]<-UMFI_par(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y,mod_meth = "lin")
    Imp$UMFI_OT[i,]<-UMFI_par(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y,mod_meth = "otpw")
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

registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

simtests<-c("Correlated_Interaction","Interaction","Correlation")

for(i in 1:length(simtests)){
  Composed<-getResultsPlot(simTest=simtests[i], nobs=500,niter=100,nX=4)
  
  jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/",simtests[i],".jpg"),width = 7, height = 5, units = 'in',res = 600)
  print(ggplot(Composed, aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + xlab("Variable")
        + ylab("Variable Importance (%)") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)))
  dev.off()
  print(i)
}

stopCluster(cl)


#get combined plot of interaction and correlation for paper
registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

set.seed(123)
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/InteractionCorrelation.jpg"),width = 12, height = 5, units = 'in',res = 600)
Composed<-getResultsPlot(simTest="Correlation", nobs=500,niter=100,nX=4)
plot1<-ggplot(Composed, aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + ggtitle('(a)') + xlab("Variable") +
      ylab("Variable Importance (%)") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14),
                                                axis.text.y = element_text(size=14))

Composed<-getResultsPlot(simTest="Interaction", nobs=500,niter=100,nX=4)
plot2<-ggplot(Composed, aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + ggtitle('(b)') + xlab("Variable") +
      ylab("Variable Importance (%)") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14),
                                                axis.text.y = element_text(size=14))

grid.arrange(plot1,plot2,ncol=2)
dev.off()

stopCluster(cl)




#evaluate real data with uncertainty
importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
Bigdatrand<-Bigdat

n_change<- as.numeric(colnames(Bigdat) %in% importantGenes)
for(j in 1:ncol(Bigdat)){
  if(n_change[j] ==0){
    set.seed(j)
    Bigdatrand[,j]<-Bigdat[sample(1:nrow(Bigdat),nrow(Bigdat)),j]
    colnames(Bigdatrand)[j]<-paste0(colnames(Bigdatrand)[j],"_p")
  }
}

rfmod<-ranger(y=responses,x=Bigdatrand,num.trees = 100) # accuracy of random forest
1-rfmod$prediction.error

registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

niterMCI<-200
niterUMFI<-200
BRCA_OTPW_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niterUMFI)
BRCA_LIN_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niterUMFI)
BRCA_MCI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niterMCI)
for(i in 1:niterUMFI){
  # n_change<- as.numeric(colnames(Bigdat) %in% importantGenes)
  # for(j in 1:ncol(Bigdat)){
  #   if(n_change[j] ==0){
  #     Bigdatrand[,j]<-Bigdat[sample(1:nrow(Bigdat),nrow(Bigdat)),j]
  #   }
  # }
  set.seed(i)
  samp<-sample(1:nrow(Bigdatrand),500)
  BRCA_OTPW_mat[i,]<-UMFI_par(X=Bigdatrand[samp,],y=responses[samp],mod_meth="otpw")
  BRCA_LIN_mat[i,]<-UMFI_par(X=Bigdatrand[samp,],y=responses[samp],mod_meth="lin")
  if(i<=niterMCI) BRCA_MCI_mat[i,]<-MCI_par(X=Bigdatrand[samp,],y=responses[samp],k=3)
  print(i)
}

BRCA_MCI<-apply(BRCA_MCI_mat, 2, median)
BRCA_lin<-apply(BRCA_LIN_mat, 2, median)
BRCA_otpw<-apply(BRCA_OTPW_mat, 2, median)

BRCA_MCI_u<-apply(BRCA_MCI_mat, 2, quantile,0.75)
BRCA_lin_u<-apply(BRCA_LIN_mat, 2, quantile,0.75)
BRCA_otpw_u<-apply(BRCA_OTPW_mat, 2, quantile,0.75)

BRCA_MCI_l<-apply(BRCA_MCI_mat, 2, quantile,0.25)
BRCA_lin_l<-apply(BRCA_LIN_mat, 2, quantile,0.25)
BRCA_otpw_l<-apply(BRCA_OTPW_mat, 2, quantile,0.25)


getcols<-function(nam){
  importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
  colnum<-rep(1,length(nam))
  colnum[nam %in% importantGenes]<-2
  colcol<-c("grey","lightblue")[colnum]
  colcol
}
cluster<- as.numeric(colnames(Bigdatrand) %in% importantGenes)

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCArand.jpg"),width = 9, height = 9, units = 'in',res = 600)
par(mar = c(6, 5, 3, 3), mfrow=c(3,1))
med<-barplot(BRCA_MCI[order(- BRCA_MCI)], ylab = "MCI Importance", cex.axis=1.1,cex.lab=1.5,main="(a)",ylim = c(0,max(BRCA_MCI_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_MCI)], col = getcols(colnames(Bigdatrand)[order(- BRCA_MCI)]),las=2)
segments(med, BRCA_MCI_l[order(- BRCA_MCI)], med, BRCA_MCI_u[order(- BRCA_MCI)], lwd = 1.5)
arrows(med, BRCA_MCI_l[order(- BRCA_MCI)], med,BRCA_MCI_u[order(- BRCA_MCI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_lin[order(- BRCA_lin)], ylab = "UMFI_LR Importance",cex.axis=1.1,cex.lab=1.5,main="(b)",ylim = c(0,max(BRCA_MCI_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_lin)], col = getcols(colnames(Bigdatrand)[order(- BRCA_lin)]),las=2)
segments(med, BRCA_lin_l[order(- BRCA_lin)], med, BRCA_lin_u[order(- BRCA_lin)], lwd = 1.5)
arrows(med, BRCA_lin_l[order(- BRCA_lin)], med,BRCA_lin_u[order(- BRCA_lin)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_otpw[order(- BRCA_otpw)], ylab = "UMFI_OT Importance", cex.axis=1.1,cex.lab=1.5, main="(c)",ylim = c(0,max(BRCA_MCI_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_otpw)], col = getcols(colnames(Bigdatrand)[order(- BRCA_otpw)]),las=2)
segments(med, BRCA_otpw_l[order(- BRCA_otpw)], med, BRCA_otpw_u[order(- BRCA_otpw)], lwd = 1.5)
arrows(med, BRCA_otpw_l[order(- BRCA_otpw)], med,BRCA_otpw_u[order(- BRCA_otpw)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)
dev.off()

stopCluster(cl)

class_mci<-as.numeric(BRCA_MCI>0)
class_lin<-as.numeric(BRCA_lin>0)
class_otpw<-as.numeric(BRCA_otpw>0)


conf_matrix_mci<-table(class_mci,cluster)
conf_matrix_lin<-table(class_lin,cluster)
conf_matrix_ot<-table(class_otpw,cluster)

sensitivity(as.factor(class_mci),as.factor(cluster))
specificity(as.factor(class_mci),as.factor(cluster))

sensitivity(as.factor(class_lin),as.factor(cluster))
specificity(as.factor(class_lin),as.factor(cluster))

sensitivity(as.factor(class_otpw),as.factor(cluster))
specificity(as.factor(class_otpw),as.factor(cluster))


mean(class_mci==cluster)
mean(class_lin==cluster)
mean(class_otpw==cluster)



# show that as we increase the number of iterations, more zeros appear for duds
registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

niterUMFI<-10000
BRCA_OTPW_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niterUMFI)
BRCA_LIN_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niterUMFI)
for(i in 1:niterUMFI){
  set.seed(i)
  samp<-sample(1:nrow(Bigdatrand),500)
  BRCA_OTPW_mat[i,]<-UMFI_par(X=Bigdatrand[samp,],y=responses[samp],mod_meth="otpw")
  BRCA_LIN_mat[i,]<-UMFI_par(X=Bigdatrand[samp,],y=responses[samp],mod_meth="lin")
  print(i)
}

BRCA_lin<-apply(BRCA_LIN_mat, 2, median)
BRCA_otpw<-apply(BRCA_OTPW_mat, 2, median)

BRCA_lin_u<-apply(BRCA_LIN_mat, 2, quantile,0.75)
BRCA_otpw_u<-apply(BRCA_OTPW_mat, 2, quantile,0.75)

BRCA_lin_l<-apply(BRCA_LIN_mat, 2, quantile,0.25)
BRCA_otpw_l<-apply(BRCA_OTPW_mat, 2, quantile,0.25)


getcols<-function(nam){
  importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
  colnum<-rep(1,length(nam))
  colnum[nam %in% importantGenes]<-2
  colcol<-c("grey","lightblue")[colnum]
  colcol
}
cluster<- as.numeric(colnames(Bigdatrand) %in% importantGenes)

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCArand_5000UMFI.jpg"),width = 12, height = 9, units = 'in',res = 600)
par(mar = c(6, 5, 3, 3), mfrow=c(2,1))

med<-barplot(BRCA_lin[order(- BRCA_lin)], ylab = "UMFI_LR Importance",cex.axis=1,cex.lab=1.2,main="(a)",ylim = c(0,max(BRCA_otpw_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_lin)], col = getcols(colnames(Bigdatrand)[order(- BRCA_lin)]),las=2)
segments(med, BRCA_lin_l[order(- BRCA_lin)], med, BRCA_lin_u[order(- BRCA_lin)], lwd = 1.5)
arrows(med, BRCA_lin_l[order(- BRCA_lin)], med,BRCA_lin_u[order(- BRCA_lin)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.2)

med<-barplot(BRCA_otpw[order(- BRCA_otpw)], ylab = "UMFI_OT Importance", cex.axis=1,cex.lab=1.2, main="(b)",ylim = c(0,max(BRCA_otpw_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_otpw)], col = getcols(colnames(Bigdatrand)[order(- BRCA_otpw)]),las=2)
segments(med, BRCA_otpw_l[order(- BRCA_otpw)], med, BRCA_otpw_u[order(- BRCA_otpw)], lwd = 1.5)
arrows(med, BRCA_otpw_l[order(- BRCA_otpw)], med,BRCA_otpw_u[order(- BRCA_otpw)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.2)
dev.off()

stopCluster(cl)

class_lin<-as.numeric(BRCA_lin>0)
class_otpw<-as.numeric(BRCA_otpw>0)


conf_matrix_lin<-table(class_lin,cluster)
conf_matrix_ot<-table(class_otpw,cluster)


sensitivity(as.factor(class_lin),as.factor(cluster))
specificity(as.factor(class_lin),as.factor(cluster))

sensitivity(as.factor(class_otpw),as.factor(cluster))
specificity(as.factor(class_otpw),as.factor(cluster))

mean(class_lin==cluster)
mean(class_otpw==cluster)




#calculate the most amount of features that can be done within an hour
extracol=10000
BigBigdat<-cbind(Bigdat,matrix(rnorm(nrow(Bigdat)*extracol),nrow=nrow(Bigdat),ncol = extracol))

registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

time_MCI<-list()
p<-list()
curtime<-0
i=1
while(curtime<1){
  p[[i]]<-5+(i-1)*2
  start_time <- Sys.time()
  fi<-MCI_par(y=responses,X=BigBigdat[,1:p[[i]]],k=p[[i]])
  end_time <-Sys.time()
  curtime<-as.numeric(end_time- start_time,units="hours")
  time_MCI[[i]]<-curtime
  print(paste(p[[i]],curtime))
  i=i+1
}

time_MCI<-as.numeric(time_MCI)
p_MCI<-as.numeric(p)
time_dat<-data.frame(p_MCI=p_MCI,time_MCI=time_MCI)
write.csv(time_dat,"C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_MCI.csv",row.names = F)



time_MCI_approx<-list()
p<-list()
curtime<-0
i=1
while(curtime<1){
  p[[i]]<-5+(i-1)*10
  start_time <- Sys.time()
  fi<-MCI_par(y=responses,X=BigBigdat[,1:p[[i]]],k=3)
  end_time <-Sys.time()
  curtime<-as.numeric(end_time- start_time,units="hours")
  time_MCI_approx[[i]]<-curtime
  print(paste(p[[i]],curtime))
  i=i+1
}

time_MCI_approx<-as.numeric(time_MCI_approx)
p_MCI_approx<-as.numeric(p)

time_dat<-data.frame(p_MCI_approx=p_MCI_approx,time_MCI_approx=time_MCI_approx)
write.csv(time_dat,"C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_MCI_approx.csv",row.names = F)


time_UMFI<-list()
p<-list()
curtime<-0
i=1
while(curtime<1){
  p[[i]]<-5+(i-1)*250
  start_time <- Sys.time()
  fi<-UMFI_par(y=responses,X=BigBigdat[,1:p[[i]]],mod_meth = "otpw")
  end_time <-Sys.time()
  curtime<-as.numeric(end_time- start_time,units="hours")
  time_UMFI[[i]]<-curtime
  print(paste(p[[i]],curtime))
  i=i+1
}

time_UMFI<-as.numeric(time_UMFI)
p_UMFI<-as.numeric(p)

time_dat<-data.frame(p_UMFI=p_UMFI,time_MCI=time_UMFI)
write.csv(time_dat,"C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_UMFI.csv",row.names = F)


time_UMFI<-list()
p<-list()
curtime<-0
i=1
while(curtime<1){
  p[[i]]<-5+(i-1)*500
  start_time <- Sys.time()
  fi<-UMFI_par(y=responses,X=BigBigdat[,1:p[[i]]],mod_meth = "lin")
  end_time <-Sys.time()
  curtime<-as.numeric(end_time- start_time,units="hours")
  time_UMFI[[i]]<-curtime
  print(paste(p[[i]],curtime))
  i=i+1
}

time_UMFI<-as.numeric(time_UMFI)
p_UMFI<-as.numeric(p)

time_dat<-data.frame(p_UMFI=p_UMFI,time_MCI=time_UMFI)
write.csv(time_dat,"C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_UMFI_lin.csv",row.names = F)
stopCluster(cl)


MCI_time<-read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_MCI.csv")
MCI_approx_time<-read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_MCI_approx.csv")
UMFI_time<-read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_UMFI.csv")
UMFI_lin_time<-read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/timebug_UMFI_lin.csv")


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/time_budget.jpg"),width = 10, height = 7, units = 'in',res = 600)
par(mar = c(4, 6, 4, 4))
plot(UMFI_lin_time$p_UMFI,UMFI_lin_time$time_MCI,type="l",col="darkblue",xlab="# of Features", 
     ylab="Time (Hours)",lwd=5, cex.axis=1.6,cex.lab=1.8, log="",ylim=c(0.01,1))
lines(MCI_approx_time$p_MCI_approx,MCI_approx_time$time_MCI_approx,col="pink",lwd=5)
lines(MCI_time$p_MCI,MCI_time$time_MCI,col="darkred",lwd=5)
lines(UMFI_time$p_UMFI,UMFI_time$time_MCI,col="lightblue",lwd=5)
legend("bottomright",col=c("darkred","pink","lightblue","darkblue"),lty = 1,lwd=5,legend = c("MCI","Approximate MCI","UMFI_OT","UMFI_LR"),cex=1.8)
dev.off()








###############################################################################################################################
##################################################3    supplement experiments    #############################################
##############################################################################################################################
rm(list=ls())
gc()
library(randomForest)
library(permimp)
library(rje)
library(iml)
library(foreach)
library(doParallel)

#modify data functions
modifty_linreg<-function(dat,protect){
  #remove dependedence via linear regression
  modifiedDAT<-dat
  tomodify<-setdiff(1:ncol(dat),protect)
  if(sd(dat[,protect])<1e-6) dat[,protect]<-dat[,protect]+rnorm(nrow(dat))
  for(i in tomodify){
    if(sd(dat[,i])<1e-6) dat[,i]<-dat[,i]+rnorm(nrow(dat))
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
  z=dat[,protect] #the protected attribute
  n_quan=ceiling(nrow(dat)/150) #number of quantiles to use (20 points per regression)
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
      mod<-lm(x_curquan~z_curquan)
      rv<-as.numeric(mod$residuals)
      condF<-rank(rv)/length(rv)
      newx[cur_obs]<-condF
    }
    modifiedDAT[,j]<-newx
  }
  modifiedDAT
}

Ablation<- function(X,y){
  rfwith <- ranger(x=X,y=y,num.trees = 100)
  OUTPUT<-rep(0,ncol(X))
  
  for(j in 1:length(OUTPUT)){
    rfwithout <- ranger(x=X[,-c(j)],y=y,num.trees = 100)
    if(is.numeric(y)) OUTPUT[j]<-max(rfwith$r.squared,0)-max(rfwithout$r.squared,0)
    if(is.factor(y)) OUTPUT[j]<- max(1-rfwith$prediction.error,0.5)-max(1-rfwithout$prediction.error,0.5)
    
  }
  OUTPUT[OUTPUT<0]<-0
  OUTPUT
}

PI<- function(X,y){
  rf <- ranger(x=X,y=y,num.trees = 100,importance = "permutation")
  OUTPUT<-as.numeric(rf$variable.importance)
  OUTPUT[OUTPUT<0]<-0
  OUTPUT
}


CPI<-function(X,y){
  readingSkills.rf <- randomForest::randomForest(x=X,y=y,ntree=100, importance = TRUE,
                                                 keep.forest = TRUE, keep.inbag = TRUE)
  permimp_rf <-permimp(readingSkills.rf, conditional = TRUE,do_check = FALSE,progressBar = FALSE)
  OUTPUT<-as.numeric(permimp_rf$values)
  OUTPUT[OUTPUT<0]<-0
  OUTPUT
}

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


UMFI_par<- function(X,y,mod_meth){
  fi<-foreach(i=1:ncol(X),  .inorder = FALSE, .export = c("modifty_otpw_quantiles_lin","modifty_linreg"),
              .packages = c("ranger", "doParallel"),.combine = 'c')%dopar%{
                if(mod_meth=="otpw") newX<-modifty_otpw_quantiles_lin(X,i)
                if(mod_meth=="lin") newX<-modifty_linreg(X,i)
                #newX[,apply(newX,2,sd)<1e-6]<-matrix(rnorm(nrow(X)*sum(apply(newX,2,sd)<1e-6)),nrow = nrow(X))
                #newX<-cbind(newX,rnorm(nrow(newX)),rnorm(nrow(newX)))
                rfwith<-ranger(x=newX,y=y,num.trees = 100)
                rfwithout<-ranger(x=newX[,-c(i)],y=y,num.trees = 100)
                if(is.numeric(y)) return(max(rfwith$r.squared,0)-max(rfwithout$r.squared,0))
                if(is.factor(y)) return(max(1-rfwith$prediction.error,0.5)-max(1-rfwithout$prediction.error,0.5))
              }
  fi[fi<0]<-0
  fi
}

#Variable importance functions
MCI<-function(X,y,k){
  colvec<-1:ncol(X)
  CompleteSet<-powerSet(colvec,m=k)
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


MCI_par<-function(X,y,k){
  colvec<-1:ncol(X)
  CompleteSet<-powerSet(colvec,m=k)
  
  CompleteSetErrors<-foreach(e=1:length(CompleteSet),  .inorder = FALSE,
                             .packages = c("ranger", "doParallel"),.combine = 'c')%dopar%{
                               if(length(CompleteSet[[e]])>0){
                                 rfmod<-ranger(y=y,x=as.data.frame(X[,CompleteSet[[e]]]),num.trees = 100)
                                 if(is.numeric(y)) return(rfmod$r.squared)
                                 if(is.factor(y)) return(1- rfmod$prediction.error)
                               }
                             }
  
  CompleteSetErrors<-c(0,CompleteSetErrors) #add accuracy for no features
  
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


Bigdat<- read.csv("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCA.csv")
Bigdat$Sample.ID<-NULL
responses<-as.factor(Bigdat$BRCA_Subtype_PAM50)
Bigdat$BRCA_Subtype_PAM50<-NULL

importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
Bigdatrand<-Bigdat

n_change<- as.numeric(colnames(Bigdat) %in% importantGenes)
for(j in 1:ncol(Bigdat)){
  if(n_change[j] ==0){
    set.seed(j)
    Bigdatrand[,j]<-Bigdat[sample(1:nrow(Bigdat),nrow(Bigdat)),j]
    colnames(Bigdatrand)[j]<-paste0(colnames(Bigdatrand)[j],"_p")
  }
}

niter<-200
BRCA_Ablation_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_PI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_CPI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
for(i in 1:niter){
  set.seed(i)
  samp<-sample(1:nrow(Bigdatrand),500)
  BRCA_Ablation_mat[i,]<-Ablation(X=Bigdatrand[samp,],y=responses[samp])
  BRCA_PI_mat[i,]<-PI(X=Bigdatrand[samp,],y=responses[samp])
  BRCA_CPI_mat[i,]<-CPI(X=Bigdatrand[samp,],y=responses[samp])
  print(i)
}

BRCA_Ablation_mat[BRCA_Ablation_mat<0]<-0
BRCA_PI_mat[BRCA_PI_mat<0]<-0
BRCA_CPI_mat[BRCA_CPI_mat<0]<-0

BRCA_Ablation<-apply(BRCA_Ablation_mat, 2, median)
BRCA_PI<-apply(BRCA_PI_mat, 2, median)
BRCA_CPI<-apply(BRCA_CPI_mat, 2, median)

BRCA_Ablation_u<-apply(BRCA_Ablation_mat, 2, quantile,0.75)
BRCA_PI_u<-apply(BRCA_PI_mat, 2, quantile,0.75)
BRCA_CPI_u<-apply(BRCA_CPI_mat, 2, quantile,0.75)

BRCA_Ablation_l<-apply(BRCA_Ablation_mat, 2, quantile,0.25)
BRCA_PI_l<-apply(BRCA_PI_mat, 2, quantile,0.25)
BRCA_CPI_l<-apply(BRCA_CPI_mat, 2, quantile,0.25)


getcols<-function(nam){
  importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
  colnum<-rep(1,length(nam))
  colnum[nam %in% importantGenes]<-2
  colcol<-c("grey","lightblue")[colnum]
  colcol
}
cluster<- as.numeric(colnames(Bigdatrand) %in% importantGenes)

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCArand_OTHERmethods.jpg"),width = 9, height = 9, units = 'in',res = 600)
par(mar = c(6, 5, 3, 3), mfrow=c(3,1))
med<-barplot(BRCA_Ablation[order(- BRCA_Ablation)], ylab = "Ablation Importance", cex.axis=1.1,cex.lab=1.5,main="(a)",ylim = c(0,max(BRCA_Ablation_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_Ablation)], col = getcols(colnames(Bigdatrand)[order(- BRCA_Ablation)]),las=2)
segments(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med, BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5)
arrows(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med,BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_PI[order(- BRCA_PI)], ylab = "PI",cex.axis=1.1,cex.lab=1.5,main="(b)",ylim = c(0,max(BRCA_PI_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_PI)], col = getcols(colnames(Bigdatrand)[order(- BRCA_PI)]),las=2)
segments(med, BRCA_PI_l[order(- BRCA_PI)], med, BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5)
arrows(med, BRCA_PI_l[order(- BRCA_PI)], med,BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_CPI[order(- BRCA_CPI)], ylab = "CPI", cex.axis=1.1,cex.lab=1.5, main="(c)",ylim = c(0,max(BRCA_CPI_u)),
             names.arg = colnames(Bigdatrand)[order(- BRCA_CPI)], col = getcols(colnames(Bigdatrand)[order(- BRCA_CPI)]),las=2)
segments(med, BRCA_CPI_l[order(- BRCA_CPI)], med, BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5)
arrows(med, BRCA_CPI_l[order(- BRCA_CPI)], med,BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)
dev.off()


class_Ablation<-as.numeric(BRCA_Ablation>0)
class_PI<-as.numeric(BRCA_PI>0)
class_CPI<-as.numeric(BRCA_CPI>0)

sensitivity(as.factor(class_Ablation),as.factor(cluster))
specificity(as.factor(class_Ablation),as.factor(cluster))

sensitivity(as.factor(class_PI),as.factor(cluster))
specificity(as.factor(class_PI),as.factor(cluster))

sensitivity(as.factor(class_CPI),as.factor(cluster))
specificity(as.factor(class_CPI),as.factor(cluster))


mean(class_Ablation==cluster)
mean(class_PI==cluster)
mean(class_CPI==cluster)





getResultsPlot<-function(simTest,nobs,niter,nX){
  
  Imp<-list(Ablation=matrix(0,nrow = niter,ncol = nX),PI=matrix(0,nrow = niter,ncol=nX),CPI=matrix(0,nrow = niter,ncol=nX))
  for(i in 1:niter){
    if(simTest=="Correlated_Interaction"){
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
    
    Imp$Ablation[i,]<-Ablation(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y)
    
    Imp$PI[i,]<-PI(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y)
    Imp$CPI[i,]<-CPI(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y)
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

simtests<-c("Correlated_Interaction","Interaction","Correlation")

for(i in 1:length(simtests)){
  Composed<-getResultsPlot(simTest=simtests[i], nobs=500,niter=100,nX=4)
  
  jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/",simtests[i],"_OTHERmethods",".jpg"),width = 7, height = 5, units = 'in',res = 600)
  print(ggplot(Composed, aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + xlab("Variable")
        + ylab("Variable Importance (%)") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)))
  dev.off()
  print(i)
}




registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)


tests=seq(from=5,by=5,to=50)
niter=10
time_UMFI<-matrix(0,nrow=length(tests),ncol = niter)
time_MCI<-matrix(0,nrow=length(tests),ncol = niter)
time_UMFI_par<-matrix(0,nrow=length(tests),ncol = niter)
time_MCI_par<-matrix(0,nrow=length(tests),ncol = niter)
time_Ablation<-matrix(0,nrow=length(tests),ncol = niter)
time_PI<-matrix(0,nrow=length(tests),ncol = niter)
time_CPI<-matrix(0,nrow=length(tests),ncol = niter)
for(i in 1:length(tests)){
  for(n in 1:niter){
    start_time <- Sys.time()
    fi<-MCI(y=responses,X=Bigdat[,1:(tests[i])],k=3)
    end_time <-Sys.time()
    time_MCI[i,n]<-as.numeric(end_time- start_time,units="mins")
    
    start_time <- Sys.time()
    fi<-UMFI(y=responses,X=Bigdat[,1:(tests[i])],mod_meth="otpw")
    end_time <-Sys.time()
    time_UMFI[i,n]<-as.numeric(end_time- start_time,units="mins")
    
    start_time <- Sys.time()
    fi<-MCI_par(y=responses,X=Bigdat[,1:(tests[i])],k=3)
    end_time <-Sys.time()
    time_MCI_par[i,n]<-as.numeric(end_time- start_time,units="mins")
    
    start_time <- Sys.time()
    fi<-UMFI_par(y=responses,X=Bigdat[,1:(tests[i])],mod_meth="otpw")
    end_time <-Sys.time()
    time_UMFI_par[i,n]<-as.numeric(end_time- start_time,units="mins")
    
    start_time <- Sys.time()
    fi<-Ablation(y=responses,X=Bigdat[,1:(tests[i])])
    end_time <-Sys.time()
    time_Ablation[i,n]<-as.numeric(end_time- start_time,units="mins")
    
    start_time <- Sys.time()
    fi<-PI(y=responses,X=Bigdat[,1:(tests[i])])
    end_time <-Sys.time()
    time_PI[i,n]<-as.numeric(end_time- start_time,units="mins")
    
    start_time <- Sys.time()
    fi<-CPI(y=responses,X=Bigdat[,1:(tests[i])])
    end_time <-Sys.time()
    time_CPI[i,n]<-as.numeric(end_time- start_time,units="mins")
  }
  print(tests[i])
}

stopCluster(cl)

time_MCI_avg<-apply(time_MCI, 1, mean)
time_UMFI_avg<-apply(time_UMFI, 1, mean)
time_Ablation_avg<-apply(time_Ablation, 1, mean)
time_PI_avg<-apply(time_PI, 1, mean)
time_CPI_avg<-apply(time_CPI, 1, mean)
time_MCI_par_avg<-apply(time_MCI_par, 1, mean)
time_UMFI_par_avg<-apply(time_UMFI_par, 1, mean)


options(scipen=5)
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/time_complexity_OTHERmethods.jpg"),width = 12, height = 8, units = 'in',res = 600)
par(mar = c(4, 6, 4, 4))
plot(tests,time_MCI_avg,type="l",col="red",xlab="# of Features", ylab="Time Complexity (Minutes)",lwd=5, cex.axis=1.6,cex.lab=1.8,ylim = c(0.000,1))
lines(tests,time_UMFI_avg,col="pink",lwd=5)
lines(tests,time_UMFI_par_avg,col="pink",lwd=5,lty=2)
lines(tests,time_MCI_par_avg,col="red",lwd=5,lty=2)
lines(tests,time_Ablation_avg,col="darkblue",lwd=5)
lines(tests,time_PI_avg,col="blue",lwd=5)
lines(tests,time_CPI_avg,col="lightblue",lwd=5)
legend("topleft",col=c("red","red","pink", "pink","darkblue","blue","lightblue"),lty = c(1,2,1,2,1,1,1),lwd=5,
       legend = c("MCI","MCI_par","UMFI","UMFI_par","Ablation","PI","CPI"),cex=1.8)
dev.off()


#feature importance on the unperturbed data
niter<-100
BRCA_Ablation_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_PI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_CPI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
for(i in 1:niter){
  set.seed(i)
  samp<-sample(1:nrow(Bigdat),500)
  BRCA_Ablation_mat[i,]<-Ablation(X=Bigdat[samp,],y=responses[samp])
  BRCA_PI_mat[i,]<-PI(X=Bigdat[samp,],y=responses[samp])
  BRCA_CPI_mat[i,]<-CPI(X=Bigdat[samp,],y=responses[samp])
  print(i)
}

BRCA_Ablation_mat[BRCA_Ablation_mat<0]<-0
BRCA_PI_mat[BRCA_PI_mat<0]<-0
BRCA_CPI_mat[BRCA_CPI_mat<0]<-0

BRCA_Ablation<-apply(BRCA_Ablation_mat, 2, median)
BRCA_PI<-apply(BRCA_PI_mat, 2, median)
BRCA_CPI<-apply(BRCA_CPI_mat, 2, median)

BRCA_Ablation_u<-apply(BRCA_Ablation_mat, 2, quantile,0.75)
BRCA_PI_u<-apply(BRCA_PI_mat, 2, quantile,0.75)
BRCA_CPI_u<-apply(BRCA_CPI_mat, 2, quantile,0.75)

BRCA_Ablation_l<-apply(BRCA_Ablation_mat, 2, quantile,0.25)
BRCA_PI_l<-apply(BRCA_PI_mat, 2, quantile,0.25)
BRCA_CPI_l<-apply(BRCA_CPI_mat, 2, quantile,0.25)


getcols<-function(nam){
  importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
  colnum<-rep(1,length(nam))
  colnum[nam %in% importantGenes]<-2
  colcol<-c("grey","lightblue")[colnum]
  colcol
}
cluster<- as.numeric(colnames(Bigdat) %in% importantGenes)

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCA_OTHERmethods.jpg"),width = 9, height = 9, units = 'in',res = 600)
par(mar = c(6, 5, 3, 3), mfrow=c(3,1))
med<-barplot(BRCA_Ablation[order(- BRCA_Ablation)], ylab = "Ablation Importance", cex.axis=1.1,cex.lab=1.5,main="(a)",ylim = c(0,max(BRCA_Ablation_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_Ablation)], col = getcols(colnames(Bigdat)[order(- BRCA_Ablation)]),las=2)
segments(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med, BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5)
arrows(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med,BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_PI[order(- BRCA_PI)], ylab = "PI",cex.axis=1.1,cex.lab=1.5,main="(b)",ylim = c(0,max(BRCA_PI_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_PI)], col = getcols(colnames(Bigdat)[order(- BRCA_PI)]),las=2)
segments(med, BRCA_PI_l[order(- BRCA_PI)], med, BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5)
arrows(med, BRCA_PI_l[order(- BRCA_PI)], med,BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_CPI[order(- BRCA_CPI)], ylab = "CPI", cex.axis=1.1,cex.lab=1.5, main="(c)",ylim = c(0,max(BRCA_CPI_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_CPI)], col = getcols(colnames(Bigdat)[order(- BRCA_CPI)]),las=2)
segments(med, BRCA_CPI_l[order(- BRCA_CPI)], med, BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5)
arrows(med, BRCA_CPI_l[order(- BRCA_CPI)], med,BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)
dev.off()



#feature importance on the unperturbed data
registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

niter<-100
BRCA_Ablation_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_PI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_CPI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
for(i in 1:niter){
  set.seed(i)
  samp<-sample(1:nrow(Bigdat),500)
  BRCA_Ablation_mat[i,]<-MCI_par(X=Bigdat[samp,],y=responses[samp],k=3)
  BRCA_PI_mat[i,]<-UMFI_par(X=Bigdat[samp,],y=responses[samp],mod_meth = "lin")
  BRCA_CPI_mat[i,]<-UMFI_par(X=Bigdat[samp,],y=responses[samp],mod_meth = "otpw")
  print(i)
}

stopCluster(cl)

BRCA_Ablation_mat[BRCA_Ablation_mat<0]<-0
BRCA_PI_mat[BRCA_PI_mat<0]<-0
BRCA_CPI_mat[BRCA_CPI_mat<0]<-0

BRCA_Ablation<-apply(BRCA_Ablation_mat, 2, median)
BRCA_PI<-apply(BRCA_PI_mat, 2, median)
BRCA_CPI<-apply(BRCA_CPI_mat, 2, median)

BRCA_Ablation_u<-apply(BRCA_Ablation_mat, 2, quantile,0.75)
BRCA_PI_u<-apply(BRCA_PI_mat, 2, quantile,0.75)
BRCA_CPI_u<-apply(BRCA_CPI_mat, 2, quantile,0.75)

BRCA_Ablation_l<-apply(BRCA_Ablation_mat, 2, quantile,0.25)
BRCA_PI_l<-apply(BRCA_PI_mat, 2, quantile,0.25)
BRCA_CPI_l<-apply(BRCA_CPI_mat, 2, quantile,0.25)


getcols<-function(nam){
  importantGenes<-c("BCL11A","EZH2","IGF1R","LFNG","BRCA1","SLC22A5","CDK6","BRCA2","TEX14","CCND1")
  colnum<-rep(1,length(nam))
  colnum[nam %in% importantGenes]<-2
  colcol<-c("grey","lightblue")[colnum]
  colcol
}
cluster<- as.numeric(colnames(Bigdat) %in% importantGenes)

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/BRCA.jpg"),width = 9, height = 9, units = 'in',res = 600)
par(mar = c(6, 5, 3, 3), mfrow=c(3,1))
med<-barplot(BRCA_Ablation[order(- BRCA_Ablation)], ylab = "MCI", cex.axis=1.1,cex.lab=1.5,main="(a)",ylim = c(0,max(BRCA_Ablation_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_Ablation)], col = getcols(colnames(Bigdat)[order(- BRCA_Ablation)]),las=2)
segments(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med, BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5)
arrows(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med,BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_PI[order(- BRCA_PI)], ylab = "UMFI_LR",cex.axis=1.1,cex.lab=1.5,main="(b)",ylim = c(0,max(BRCA_PI_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_PI)], col = getcols(colnames(Bigdat)[order(- BRCA_PI)]),las=2)
segments(med, BRCA_PI_l[order(- BRCA_PI)], med, BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5)
arrows(med, BRCA_PI_l[order(- BRCA_PI)], med,BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)

med<-barplot(BRCA_CPI[order(- BRCA_CPI)], ylab = "UMFI_OT", cex.axis=1.1,cex.lab=1.5, main="(c)",ylim = c(0,max(BRCA_CPI_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_CPI)], col = getcols(colnames(Bigdat)[order(- BRCA_CPI)]),las=2)
segments(med, BRCA_CPI_l[order(- BRCA_CPI)], med, BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5)
arrows(med, BRCA_CPI_l[order(- BRCA_CPI)], med,BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
legend("topright",col=c("lightblue","grey"),lty = 1,lwd=10,legend = c("BRCA Associated","Not BRCA Associated"),cex=1.5)
dev.off()




######################################################################################################################################
# experiments with extra trees on hydrology dataset
#####################################################################################################################################
modifty_otpw_quantiles_lin<-function(dat,protect){
  modifiedDAT<-dat #the new dataframe that will be returned
  tomodify<-setdiff(1:ncol(dat),protect) #the columns in dat to modify
  z=dat[,protect] #the protected attribute
  n_quan=ceiling(nrow(dat)/30) #number of quantiles to use (20 points per regression)
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
      mod<-lm(x_curquan~z_curquan)
      rv<-as.numeric(mod$residuals)
      condF<-rank(rv)/length(rv)
      newx[cur_obs]<-condF
      #for(k in 1:length(x_curquan)){
      #  condF[k]<- mean(rv<= rv[k])
      #  newx[cur_obs][k]<-condF[k]
      #}
    }
    modifiedDAT[,j]<-newx
  }
  modifiedDAT
}

UMFI_par<- function(X,y,mod_meth){
  fi<-foreach(i=1:ncol(X),  .inorder = FALSE, .export = c("modifty_otpw_quantiles_lin","modifty_linreg"),
              .packages = c("ranger", "doParallel"),.combine = 'c')%dopar%{
                if(mod_meth=="otpw") newX<-modifty_otpw_quantiles_lin(X,i)
                if(mod_meth=="lin") newX<-modifty_linreg(X,i)
                #newX[,apply(newX,2,sd)<1e-6]<-matrix(rnorm(nrow(X)*sum(apply(newX,2,sd)<1e-6)),nrow = nrow(X))
                #newX<-cbind(newX,rnorm(nrow(newX)),rnorm(nrow(newX)))
                rfwith<-ranger(x=newX,y=y,num.trees = 100,splitrule = "extratrees")
                rfwithout<-ranger(x=newX[,-c(i)],y=y,num.trees = 100,splitrule = "extratrees")
                if(is.numeric(y)) return(max(rfwith$r.squared,0)-max(rfwithout$r.squared,0))
                if(is.factor(y)) return(max(1-rfwith$prediction.error,0.5)-max(1-rfwithout$prediction.error,0.5))
              }
  fi[fi<0]<-0
  fi
}

MCI_par<-function(X,y,k){
  colvec<-1:ncol(X)
  CompleteSet<-powerSet(colvec,m=k)
  
  CompleteSetErrors<-foreach(e=1:length(CompleteSet),  .inorder = FALSE,
                             .packages = c("ranger", "doParallel"),.combine = 'c')%dopar%{
                               if(length(CompleteSet[[e]])>0){
                                 rfmod<-ranger(y=y,x=as.data.frame(X[,CompleteSet[[e]]]),num.trees = 100,splitrule = "extratrees")
                                 if(is.numeric(y)) return(rfmod$r.squared)
                                 if(is.factor(y)) return(1- rfmod$prediction.error)
                               }
                             }
  
  CompleteSetErrors<-c(0,CompleteSetErrors) #add accuracy for no features
  
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

Bigdat<-dat[complete.cases(dat),]

responses<-Bigdat[,(ncol(Bigdat)-5):ncol(Bigdat)]
responses<-responses$q_mean
Bigdat[,(ncol(Bigdat)-5):ncol(Bigdat)]<-NULL


rfmod<-ranger(y=responses,x=Bigdat,num.trees = 100,splitrule = "extratrees") # accuracy of extra trees 
rfmod$r.squared
cor(rfmod$predictions,responses)^2


rf_none<-calc_dependence(dat=Bigdat,mod_meth = "nothing")
rf_resid<-calc_dependence(dat=Bigdat,mod_meth = "lin")
rf_ot_pw<-calc_dependence(dat=Bigdat,mod_meth = "ot_pw_quan")


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/HydrologyMutualInfoR2.jpeg"),width = 12, height = 9, units = 'in',res = 600)
par(mar = c(4, 6, 1, 1))
plot(rf_none,col="black",xlab="Feature",ylab = "Mutual Information (R^2)",type="l",ylim=c(-0.1,1),lwd=5, cex.axis=1.6,cex.lab=1.8)
lines(rf_ot_pw,col="blue",lwd=5)
lines(rf_resid,col="red",lwd=5)
legend("bottomleft",col=c("black","red", "blue"),lty = 1,lwd=5,legend = c("Raw Data","Linear Regression","Pairwise Optimal Transport"),cex=1.8)
dev.off()





registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

niter<-100
BRCA_Ablation_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_PI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
BRCA_CPI_mat<-matrix(0,ncol = ncol(Bigdat),nrow = niter)
for(i in 1:niter){
  set.seed(i)
  samp<-sample(1:nrow(Bigdat),500)
  BRCA_Ablation_mat[i,]<-MCI_par(X=Bigdat[samp,],y=responses[samp],k=3)
  BRCA_PI_mat[i,]<-UMFI_par(X=Bigdat[samp,],y=responses[samp],mod_meth = "lin")
  BRCA_CPI_mat[i,]<-UMFI_par(X=Bigdat[samp,],y=responses[samp],mod_meth = "otpw")
  print(i)
}

stopCluster(cl)

BRCA_Ablation_mat[BRCA_Ablation_mat<0]<-0
BRCA_PI_mat[BRCA_PI_mat<0]<-0
BRCA_CPI_mat[BRCA_CPI_mat<0]<-0

BRCA_Ablation<-apply(BRCA_Ablation_mat, 2, median)
BRCA_PI<-apply(BRCA_PI_mat, 2, median)
BRCA_CPI<-apply(BRCA_CPI_mat, 2, median)

BRCA_Ablation_u<-apply(BRCA_Ablation_mat, 2, quantile,0.75)
BRCA_PI_u<-apply(BRCA_PI_mat, 2, quantile,0.75)
BRCA_CPI_u<-apply(BRCA_CPI_mat, 2, quantile,0.75)

BRCA_Ablation_l<-apply(BRCA_Ablation_mat, 2, quantile,0.25)
BRCA_PI_l<-apply(BRCA_PI_mat, 2, quantile,0.25)
BRCA_CPI_l<-apply(BRCA_CPI_mat, 2, quantile,0.25)


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/Hydrology.jpg"),width = 9, height = 11, units = 'in',res = 600)
par(mar = c(10, 5, 3, 3), mfrow=c(3,1))
med<-barplot(BRCA_Ablation[order(- BRCA_Ablation)], ylab = "MCI", cex.axis=1.1,cex.lab=1.5,main="(a)",ylim = c(0,max(BRCA_Ablation_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_Ablation)], col = "lightblue",las=2)
segments(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med, BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5)
arrows(med, BRCA_Ablation_l[order(- BRCA_Ablation)], med,BRCA_Ablation_u[order(- BRCA_Ablation)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)

med<-barplot(BRCA_PI[order(- BRCA_PI)], ylab = "UMFI_LR",cex.axis=1.1,cex.lab=1.5,main="(b)",ylim = c(0,max(BRCA_PI_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_PI)], col = "lightblue",las=2)
segments(med, BRCA_PI_l[order(- BRCA_PI)], med, BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5)
arrows(med, BRCA_PI_l[order(- BRCA_PI)], med,BRCA_PI_u[order(- BRCA_PI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)

med<-barplot(BRCA_CPI[order(- BRCA_CPI)], ylab = "UMFI_OT", cex.axis=1.1,cex.lab=1.5, main="(c)",ylim = c(0,max(BRCA_CPI_u)),
             names.arg = colnames(Bigdat)[order(- BRCA_CPI)], col = "lightblue",las=2)
segments(med, BRCA_CPI_l[order(- BRCA_CPI)], med, BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5)
arrows(med, BRCA_CPI_l[order(- BRCA_CPI)], med,BRCA_CPI_u[order(- BRCA_CPI)], lwd = 1.5, 
       angle = 90,code = 3, length = 0.05)
dev.off()



getResultsPlot<-function(simTest,nobs,niter,nX){
  
  Imp<-list(MCI=matrix(0,nrow = niter,ncol = nX),UMFI_LR=matrix(0,nrow = niter,ncol=nX),UMFI_OT=matrix(0,nrow = niter,ncol=nX))
  for(i in 1:niter){
    if(simTest=="Correlated_Interaction"){
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
    
    Imp$MCI[i,]<-MCI_par(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y)
    
    Imp$UMFI_LR[i,]<-UMFI_par(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y,mod_meth = "lin")
    Imp$UMFI_OT[i,]<-UMFI_par(X=Boston[,1:(ncol(Boston)-1)],y=Boston$y,mod_meth = "otpw")
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

registerDoParallel(cores = 12)
cl <- makeCluster(12)
registerDoParallel(cl)

simtests<-c("Correlated_Interaction","Interaction","Correlation")

for(i in 1:length(simtests)){
  Composed<-getResultsPlot(simTest=simtests[i], nobs=500,niter=100,nX=4)
  
  jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH605D/Project/",simtests[i],"ExtraTree.jpg"),width = 7, height = 5, units = 'in',res = 600)
  print(ggplot(Composed, aes(x=variable, y=value)) + geom_boxplot() + facet_grid(~name) + xlab("Variable")
        + ylab("Variable Importance (%)") + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)))
  dev.off()
  print(i)
}

stopCluster(cl)




