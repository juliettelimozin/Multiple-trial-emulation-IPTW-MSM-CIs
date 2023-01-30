
## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 


DATA_GEN_censored_reduced<-function(ns, nv, conf = 0.5, treat_prev = 0, 
                                    outcome_prev = -3.8, all_treat = FALSE, 
                                    all_control = FALSE, censor = TRUE){   
  # ns= number of subjects, nv=no of visits including baseline visit
  
  
  nvisit<-nv+1
  
  X2<-rep(0,nvisit*ns)          ## place holders for time-varying covariates
  Z2<-rnorm(nvisit*ns,0,1)
  X4<-rep(rnorm(ns,0,1),each=nvisit) # baseline continuous covariate
  
  A<-rep(0,nvisit*ns) ##place holders for current  treatments
  Ap<-rep(0,nvisit*ns) ##place holders for  previous treatments
  
  CAp<-rep(0,nvisit*ns)   ##place holders for sum of previous treatment A
  
  
  
  Y<-rep(0,nvisit*ns)     ##place holders for outcome vector
  Yp<-rep(0,nvisit*ns)    ##place holders for previous outcome vector
  
  ##Fill in initial values
  seq1<-seq(1,nvisit*ns-nv,nvisit)
  
  X2[seq1]<-0
  
  P0<-list() ##list of treatment probabilities 
  P0[[1]]<-rep(0, ns) 
  seqlist<-list()                              
  seqlist[[1]]<-seq1
  CAp[seq1]<-rep(0, ns)
  
  for (k in 2:nvisit){  
    ## update covariates
    
    
    seqlist[[k]]<-seqlist[[k-1]]+1
    Ap[seqlist[[k]]]<-A[seqlist[[k-1]]]
    CAp[seqlist[[k]]]<-CAp[seqlist[[k-1]]]+Ap[seqlist[[k]]]
    
    X2[seqlist[[k]]]<-Z2[seqlist[[k]]]-0.3*Ap[seqlist[[k]]] ## continuous time-varying confounder 
    
    ## update treatment
    
    ########### Old formula: lpp<- as.numeric(treat_prev) + Ap[seqlist[[k]]]+0.5*X1[seqlist[[k]]]+as.numeric(conf)*X2[seqlist[[k]]]
    ###########                    -0.2*X3[seqlist[[k]]]+X4[seqlist[[k]]]-0.3*(age[seqlist[[k]]]-35)/12
    lpp<- as.numeric(treat_prev) + 0.05*Ap[seqlist[[k]]] + as.numeric(conf)*X2[seqlist[[k]]] + 0.2*X4[seqlist[[k]]]
    P0[[k]]<-1/(1+exp(-lpp))
    
    if (all_treat == TRUE){
      A[seqlist[[k]]]<- 1.0
    } else{ if (all_control == TRUE){
      A[seqlist[[k]]]<- 0.0
    } else{
      A[seqlist[[k]]]<-rbinom(ns,1,P0[[k]]) ##Generate treatment at current visit based on  covariates, previous treatment
    }
    }
    ##Generate outcome
    
    ##### Old formula: intercept was -7 -3.7
    lp<- as.numeric(outcome_prev) -0.5*A[seqlist[[k]]]+as.numeric(conf)*X2[seqlist[[k]]]+X4[seqlist[[k]]]
    
    Yp[seqlist[[k]]]<-Y[seqlist[[k-1]]]
    Y[seqlist[[k]]]<-(rbinom(ns,1,1/(1+exp(-lp))))*as.numeric(Yp[seqlist[[k]]]==0)+as.numeric(Yp[seqlist[[k]]]==1)
    
  }
  
  ##Make data frame
  
  ID<-rep(1:ns,each=nv)
  
  ##Align data by removing values 
  NSEQ<-seq1
  
  X2<-X2[-NSEQ]
  X4<-X4[-NSEQ]

  A<-A[-NSEQ]
  Ap<-Ap[-NSEQ]
  CAp<-CAp[-NSEQ]
  Y<-Y[-seq(1,nvisit*ns-nv,nvisit)]
  Yp<-Yp[-seq(1,nvisit*ns-nv,nvisit)]
  
  ##Create data frame
  
  DATA<-data.frame(ID,t=rep(c(0:(nv-1)),ns),A,Ap,CAp,X2,X4,Y,Yp)
  DATA$eligible<-as.numeric(CAp==0 & Yp==0)  ## eligibility criteria: age>=18, had no treatment so far, no event so far
  
  ##censoring
  if (censor == T){
    Dprob<-1/(1+exp(2.5 + Ap-0.5*X2-0.2*X4)) ##Probability of dropout
    
    DATA$C<-rbinom(nv*ns,1,Dprob) ##C=0 is remain in the study
    
    
    
    indfun<-function(n){
      if (sum(n)==0) {rep(0,nv)}
      else{k<-min(which(n==1))
      c(rep(0,k),rep(1,nv-k))}}
    
    RL<-ave(DATA$C,DATA$ID,FUN=indfun)
    
    
    eligCum<-ave(DATA$eligible,DATA$ID,FUN=cumsum)
    
    DATA[RL==0 & DATA$Yp==0 & eligCum>0,] #remove observations after event occurrence and censoring, and not eligible
  } else {
    DATA$C <- 0
    eligCum<-ave(DATA$eligible,DATA$ID,FUN=cumsum)
    
    DATA[DATA$Yp==0 & eligCum>0,] 
  }
  
}