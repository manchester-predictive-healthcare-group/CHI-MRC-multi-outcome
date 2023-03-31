## This file will calculate the outcomes for cohort A and B
## Hopefully it will do it in a fraction of the time!

setwd("/mnt/bmh01-rds/mrc-multi-outcome/data_Aurum_65plus/data_intermediate/")

# install.packages("parallel", repos="http://cran.rstudio.com/", lib="/mnt/ja01-home01/mbrxsap3/phd_risk/R/packages")
# install.packages("dplyr", repos="http://cran.rstudio.com/")


library(parallel)
library(dplyr)

# Create algotihm to look for two entrys within 90 days, both below 60 
# This will search each dataset for all CKD occurences

check_all_CKD<-function(dat){
  # Create output dataset
  output<-data.frame(person_id=numeric(),EntryDate=factor(),
                     CodeValue=numeric(),
                     num=numeric(),stringsAsFactors=FALSE)
  
  # Loop through each patid using the 'num' variable
  # First check that what i'm looping through is correct
  # print(dat[1,"num"])
  # print(dat[dim(dat)[1],"num"])
  # Now loop through ID = ...
  for(ID in dat[1,"num"]:dat[dim(dat)[1],"num"]){
    d=dat[as.character(dat$num)==as.character(ID),]
    control=FALSE
    i=1
    index=0
    # date.CKD=NULL
    while((!control)&(i<dim(d)[1])){
      if((d[i,]$CodeValue<60)){
        j=1
        while((!control)&(j<=(dim(d)[1]-i))){
          if((d[i+j,]$CodeValue<60)){
            diff=(as.Date(d[i,]$EntryDate)-as.Date(d[i+j,]$EntryDate))
            if(diff>=90){
              control=TRUE
              index=i
              #date.CKD=as.character(d[index,]$EntryDate)
            }
          }else{
            j=dim(d)[1]-i
          }
          j=j+1
        }
      }
      i=i+1
    }
    # If index is > 0, add the entry to the output
    if (index > 0){
      output<-rbind(output,d[index,])
    }
  }
  # At the end of the for loop, output the dataset
  output
}

##########################################################
################ Do cohort A first #######################
##########################################################

# Read in GFR data
mydata <- read.table("/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primhist_gfr_scores_A.txt", header=TRUE)
head(mydata)

# Change relevant variable names
names(mydata)[names(mydata)=="measurement_date"]<- "EntryDate"
names(mydata)[names(mydata)=="gfr"]<- "CodeValue"

## Create datasets to go into lapply
maxnumA<-mydata[dim(mydata)[1],"num"]

n.dat<-20
for (i in 1:n.dat){
  assign(paste("mydata",i,sep=""), filter(mydata,(num<=maxnumA*i/n.dat)&(maxnumA*(i-1)/n.dat<num)))
}


# Now do the lapply stuff in the cluster
cl <- makeCluster(20)
clusterExport(cl, "check_all_CKD")

start_time_clusA20 <- Sys.time()
CKD_stage345.l<-parLapply(cl, list(mydata1,mydata2,mydata3,mydata4,mydata5,mydata6,mydata7,mydata8,mydata9,mydata10,
                   mydata11,mydata12,mydata13,mydata14,mydata15,mydata16,mydata17,mydata18,mydata19,mydata20), 
          function(x) check_all_CKD(x))
end_time_clusA20 <- Sys.time()
time_clusA20<-start_time_clusA20-end_time_clusA20
stopCluster(cl)
time_clusA20

# str(CKD_stage345.l)

## Combine results into one dataset
CKD_from_test_A <-bind_rows(CKD_stage345.l[[1]],CKD_stage345.l[[2]],CKD_stage345.l[[3]],CKD_stage345.l[[4]],
                          CKD_stage345.l[[5]],CKD_stage345.l[[6]],CKD_stage345.l[[7]],CKD_stage345.l[[8]],
                          CKD_stage345.l[[9]],CKD_stage345.l[[10]],CKD_stage345.l[[11]],CKD_stage345.l[[12]],
                          CKD_stage345.l[[13]],CKD_stage345.l[[14]],CKD_stage345.l[[15]],CKD_stage345.l[[16]],
                          CKD_stage345.l[[17]],CKD_stage345.l[[18]],CKD_stage345.l[[19]],CKD_stage345.l[[20]])

## Export to csv file
write.csv(CKD_from_test_A,file="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primhist_from_test_A.csv", row.names= FALSE)


#########################################################
################ Do cohort B next #######################
#########################################################

# Read in GFR data
mydata <- read.table("/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primhist_gfr_scores_B.txt", header=TRUE)
head(mydata)

# Change relevant variable names
names(mydata)[names(mydata)=="measurement_date"]<- "EntryDate"
names(mydata)[names(mydata)=="gfr"]<- "CodeValue"

## Create datasets to go into lapply
maxnumA<-mydata[dim(mydata)[1],"num"]

n.dat<-20
for (i in 1:n.dat){
  assign(paste("mydata",i,sep=""), filter(mydata,(num<=maxnumA*i/n.dat)&(maxnumA*(i-1)/n.dat<num)))
}


# Now do the lapply stuff in the cluster
cl <- makeCluster(20)
clusterExport(cl, "check_all_CKD")

start_time_clusA20 <- Sys.time()
CKD_stage345.l<-parLapply(cl, list(mydata1,mydata2,mydata3,mydata4,mydata5,mydata6,mydata7,mydata8,mydata9,mydata10,
                                   mydata11,mydata12,mydata13,mydata14,mydata15,mydata16,mydata17,mydata18,mydata19,mydata20), 
                          function(x) check_all_CKD(x))
end_time_clusA20 <- Sys.time()
time_clusA20<-start_time_clusA20-end_time_clusA20
stopCluster(cl)
time_clusA20

str(CKD_stage345.l)

## Combine results
CKD_from_test_B <-bind_rows(CKD_stage345.l[[1]],CKD_stage345.l[[2]],CKD_stage345.l[[3]],CKD_stage345.l[[4]],
                            CKD_stage345.l[[5]],CKD_stage345.l[[6]],CKD_stage345.l[[7]],CKD_stage345.l[[8]],
                            CKD_stage345.l[[9]],CKD_stage345.l[[10]],CKD_stage345.l[[11]],CKD_stage345.l[[12]],
                            CKD_stage345.l[[13]],CKD_stage345.l[[14]],CKD_stage345.l[[15]],CKD_stage345.l[[16]],
                            CKD_stage345.l[[17]],CKD_stage345.l[[18]],CKD_stage345.l[[19]],CKD_stage345.l[[20]])

## Write to csv file
write.csv(CKD_from_test_B,file="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primhist_from_test_B.csv", row.names= FALSE)

rm(list=setdiff(ls(),list("CKD_from_test_A", "CKD_from_test_B")))

## Save image
save.image("varCKD_from_test_histprim_image.RData")