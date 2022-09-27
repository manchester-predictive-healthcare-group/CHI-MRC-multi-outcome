#2020-11-20
###########################################################
#Anamaria Savu, PhD                                                                                         #
#Canadian VIGOUR Centre                                                                                #
#University of Alberta                                                                                        #
# savu@ualberta.ca                                                                                            #
#November 2020                                                                                                #
###########################################################
#PDI – polytomous discrimination index                                                        #
###########################################################

###############################################################
# data = input dataset of data.frame type
#             required
#              
# nbs = number of bootstrap samples to perform of numeric type
#           required
###############################################################
#  Usage notes
#  data must have one of its components named $outcome for storing 
# the outcome values. The values of the outcome must cover all integer 
# values between a lower and an upper bound. The remaining components 
# of the data frame # ($p1, $p2, …) must record the estimated probabilities
# in a specific order: $p1 probabilities for the lowest value of the outcome, 
# $p2 probabilities for the second lowest value of the outcome, and so on.
###############################################################

#Estimates of PDI and its components
pdiest<-function(data){
  
  y<-data$outcome
  ymin<-min(y)
  ymax<-max(y)
  noutcome<-ymax-ymin
  p<-prod(table(y))
  pdi<-c()
  
  for (i in 1:(noutcome+1)){
    
    predprob<-data[,(i+1)]  #READ predicted probabilities for level i
    t0<-table(predprob,y)   #CALCULATE frequencies of predicted probabilities for level i by outcome
    
    dim1<-dim(t0)[1]
    dim2<-dim(t0)[2]
    t<-cbind(t0[,i],t0[,-i]) #REORDER columns
    restrictt<- if (noutcome == 1){matrix(t[,2:(noutcome+1)],ncol=1)} else {t[,2:(noutcome+1)] } #REMOVE first column of t
    
    c<-apply(restrictt,2,cumsum) #CALCULATE cumulative frequencies of predicted probabilities for level i by outcome
    cnew<- if (noutcome == 1) {rbind(rep(0,noutcome),matrix(c[1:(dim(c)[1]-1),],ncol=))} else {rbind(rep(0,noutcome),c[1:(dim(c)[1]-1),])} #INTRODUCE a row of zeros at the begining of c
    
    mat<-c()                     #MATRIX of 0s and 1s of dimension 2^(noutcome) x noutcome
    for (j in 1:noutcome){
      mat0<-cbind(mat,0)
      mat1<-cbind(mat,1)
      mat<-rbind(mat0,mat1)}
    
    r<-0
    for (k in 1:dim(mat)[1]){
      dt<-t(apply(restrictt, 1, function(x) mat[k,]*x))
      dcnew<-t(apply(cnew, 1, function(x) (1-mat[k,])*x))
      dfinal<-if (noutcome == 1) {cbind(t[,1],t(dt+dcnew))} else {cbind(t[,1],dt+dcnew)} #TAKE all combinations of frequencies and cumulative frequencies
      r<-r+sum(apply(dfinal,1,prod))/(1+sum(mat[k,]))}                                   #MULTIPLYIES across rows
    
    r<-r/p     #PDI component for outcome i
    pdi<-rbind(pdi,r)
  }
  pdi<-rbind(mean(pdi),pdi)
  pdi}

#Estimates and bootstrap 95% confidence intervals for PDI and its components
pdifunction<-function(data,nbs){
  #PDI estimate
  estimate<-pdiest(data)
  #BOOTSTRAP
  samplesize<-dim(data)[1]
  for (i in 1:nbs)
  {vec<-sample.int(samplesize,size=samplesize, replace=TRUE)
  mydatabs<-data[vec,]
  if (i<2) {pdibs<-pdiest(mydatabs)
  } else {
    pdibs<-cbind(pdibs,pdiest(mydatabs))
  }
  }
  
  stderr <- sqrt(apply(pdibs, 1, var))
  lowerci<- pmax(0, estimate - 1.96*stderr)
  upperci<- pmin(1, estimate + 1.96*stderr)
  
  estci<-cbind(estimate, lowerci, upperci)
  estci
}
