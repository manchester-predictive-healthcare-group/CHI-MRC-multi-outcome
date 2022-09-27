### This program will summarise the results from the simulation
rm(list=ls())

## Load image
load("R_out/multinom simulation scen1to6 500000val.RData")


### Create lists of dataframes of the simulations output
scen1.dataframes <- vector("list", 5)
for (i in 1:5){scen1.dataframes[[i]] <- data.frame(scenario1.res[[(i)]])}

scen2.dataframes <- vector("list", 5)
for (i in 1:5){scen2.dataframes[[i]] <- data.frame(scenario2.res[[(i)]])}

scen3.dataframes <- vector("list", 5)
for (i in 1:5){scen3.dataframes[[i]] <- data.frame(scenario3.res[[(i)]])}

scen4.dataframes <- vector("list", 5)
for (i in 1:5){scen4.dataframes[[i]] <- data.frame(scenario4.res[[(i)]])}

scen5.dataframes <- vector("list", 5)
for (i in 1:5){scen5.dataframes[[i]] <- data.frame(scenario5.res[[(i)]])}

scen6.dataframes <- vector("list", 5)
for (i in 1:5){scen6.dataframes[[i]] <- data.frame(scenario6.res[[(i)]])}


### Give sample sizes as names to list elements
sample.sizes1 <- c(250,500,1000,541,576)
sample.sizes2 <- c(250,500,1000,569,628)
sample.sizes3 <- c(250,500,1000,566,582)
sample.sizes4 <- c(250,500,1000,648,901)
sample.sizes5 <- c(250,500,1000,706,1458)
sample.sizes6 <- c(250,500,1000,1558,1616)

names(scen1.dataframes) <- sample.sizes1
names(scen2.dataframes) <- sample.sizes2
names(scen3.dataframes) <- sample.sizes3
names(scen4.dataframes) <- sample.sizes4
names(scen5.dataframes) <- sample.sizes5
names(scen6.dataframes) <- sample.sizes6


## Want to calculate the heuristic shrinkage factor for each model
calc.S_VH <- function(data.in){
  data.out <- data.in
  for (i in 1:5){
      data.out[[i]]$S_VH_multinom <- 1 - 10/(data.out[[i]]$LR.multinom)
      data.out[[i]]$S_VH_dislog1 <- 1 - 5/(data.out[[i]]$LR.dislog.1)
      data.out[[i]]$S_VH_dislog2 <- 1 - 5/(data.out[[i]]$LR.dislog.2)
    }
  return(data.out)
}


scen1.dataframes.S <- calc.S_VH(scen1.dataframes)
scen2.dataframes.S <- calc.S_VH(scen2.dataframes)
scen3.dataframes.S <- calc.S_VH(scen3.dataframes)
scen4.dataframes.S <- calc.S_VH(scen4.dataframes)
scen5.dataframes.S <- calc.S_VH(scen5.dataframes)
scen6.dataframes.S <- calc.S_VH(scen6.dataframes)


### Calculate R2_CS_APP and R2_CS_ADJ for each model
### Going to just write a separate function for each as they all have a different set of sample sizes
for (i in 1:5){scen1.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen1.dataframes.S[[i]]$LR.multinom/sample.sizes1[i])
               scen1.dataframes.S[[i]]$R2_CS_ADJ <- scen1.dataframes.S[[i]]$R2_CS_APP*scen1.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:5){scen2.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen2.dataframes.S[[i]]$LR.multinom/sample.sizes2[i])
               scen2.dataframes.S[[i]]$R2_CS_ADJ <- scen2.dataframes.S[[i]]$R2_CS_APP*scen2.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:5){scen3.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen3.dataframes.S[[i]]$LR.multinom/sample.sizes3[i])
               scen3.dataframes.S[[i]]$R2_CS_ADJ <- scen3.dataframes.S[[i]]$R2_CS_APP*scen3.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:5){scen4.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen4.dataframes.S[[i]]$LR.multinom/sample.sizes4[i])
               scen4.dataframes.S[[i]]$R2_CS_ADJ <- scen4.dataframes.S[[i]]$R2_CS_APP*scen4.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:5){scen5.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen5.dataframes.S[[i]]$LR.multinom/sample.sizes5[i])
               scen5.dataframes.S[[i]]$R2_CS_ADJ <- scen5.dataframes.S[[i]]$R2_CS_APP*scen5.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:5){scen6.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen6.dataframes.S[[i]]$LR.multinom/sample.sizes6[i])
               scen6.dataframes.S[[i]]$R2_CS_ADJ <- scen6.dataframes.S[[i]]$R2_CS_APP*scen6.dataframes.S[[i]]$S_VH_multinom}


### Calculate median of each entity of interest and put into a table
get.output.tables.med <- function(data.in){
  tables.out <- rbind(
    # First row
    c(capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_dislog2),3), ")", sep = ""))),
    # Second row
    c(capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_dislog2),3), ")", sep = ""))),
    # Third row
    c(capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_dislog2),3), ")", sep = ""))),
    # Fourth row
    c(capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_dislog2),3), ")", sep = ""))),
    # Fifth row
    c(capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_dislog2),3), ")", sep = "")))
    )
    return(tables.out)}
      


### And apply it to each set of results
output.scen1.med <- get.output.tables.med(scen1.dataframes.S)
output.scen2.med <- get.output.tables.med(scen2.dataframes.S)
output.scen3.med <- get.output.tables.med(scen3.dataframes.S)
output.scen4.med <- get.output.tables.med(scen4.dataframes.S)
output.scen5.med <- get.output.tables.med(scen5.dataframes.S)
output.scen6.med <- get.output.tables.med(scen6.dataframes.S)

rm(list=setdiff(ls(), list("output.scen1.med","output.scen2.med",
                           "output.scen3.med","output.scen4.med",
                           "output.scen5.med","output.scen6.med")))

save.image("R_out/multinom results Table S2.RData")
