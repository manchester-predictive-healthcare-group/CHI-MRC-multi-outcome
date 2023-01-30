### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("Project 4/code/sim_function_load_packages.R")
### Load functions
source("Project 4/code/sim_function_clinical_example.R")
### Load data
variables.vec <- c("Age", "gender", "Smoking", "SBP", "Cholhdl_ratio", "IMD", "BMI", "Ethnicity6")
model.type <- "C"
source("Project 4/code/p4_clinical_example_load_data_200000.R")


##################
##################
### FIT MODELS ###
##################
##################

###
### Product
###
print("fit product")
Sys.time()
predrisk.product <- calc.predrisk.product(data.devel = data.devel, 
                                          data.valid = data.valid, 
                                          t.eval = t.eval, 
                                          variables.vec = variables.vec)

###
### Dual
###
print("fit dual")
Sys.time()
predrisk.dual <- calc.predrisk.dual(data.devel = data.devel, 
                                    data.valid = data.valid, 
                                    t.eval = t.eval, 
                                    variables.vec = variables.vec)

###
### Clayton
###
print("fit clayton")
Sys.time()
predrisk.clayton.list <- vector("list", 4)
names(predrisk.clayton.list) <- c("r0", "r90", "r180", "r270")

print("r0")
predrisk.clayton.list[["r0"]] <- calc.predrisk.copula(data.devel = data.devel, 
                                                      data.valid = data.valid, 
                                                      t.eval = t.eval, 
                                                      variables.vec = variables.vec,
                                                      copula.in = "clayton", 
                                                      rotate.in = 0)

# print("r90")
# predrisk.clayton.list[["r90"]] <- calc.predrisk.copula(data.devel = data.devel, 
#                                                        data.valid = data.valid, 
#                                                        t.eval = t.eval, 
#                                                        variables.vec = variables.vec,
#                                                        copula.in = "clayton", 
#                                                        rotate.in = 90)

print("r180")
predrisk.clayton.list[["r180"]] <- calc.predrisk.copula(data.devel = data.devel, 
                                                        data.valid = data.valid, 
                                                        t.eval = t.eval, 
                                                        variables.vec = variables.vec,
                                                        copula.in = "clayton", 
                                                        rotate.in = 180)

# print("r270")
# predrisk.clayton.list[["r270"]] <- calc.predrisk.copula(data.devel = data.devel, 
#                                                         data.valid = data.valid, 
#                                                         t.eval = t.eval, 
#                                                         variables.vec = variables.vec,
#                                                         copula.in = "clayton", 
#                                                         rotate.in = 270)


###
### Clayton
###
print("fit gumbel")
Sys.time()
predrisk.gumbel.list <- vector("list", 4)
names(predrisk.gumbel.list) <- c("r0", "r90", "r180", "r270")

print("r0")
predrisk.gumbel.list[["r0"]] <- calc.predrisk.copula(data.devel = data.devel, 
                                                     data.valid = data.valid, 
                                                     t.eval = t.eval, 
                                                     variables.vec = variables.vec,
                                                     copula.in = "gumbel", 
                                                     rotate.in = 0)

# print("r90")
# predrisk.gumbel.list[["r90"]] <- calc.predrisk.copula(data.devel = data.devel, 
#                                                        data.valid = data.valid, 
#                                                        t.eval = t.eval, 
#                                                        variables.vec = variables.vec,
#                                                        copula.in = "gumbel", 
#                                                        rotate.in = 90)

# print("r180")
# predrisk.gumbel.list[["r180"]] <- calc.predrisk.copula(data.devel = data.devel, 
#                                                        data.valid = data.valid, 
#                                                        t.eval = t.eval, 
#                                                        variables.vec = variables.vec,
#                                                        copula.in = "gumbel", 
#                                                        rotate.in = 180)

# print("r270")
# predrisk.gumbel.list[["r270"]] <- calc.predrisk.copula(data.devel = data.devel, 
#                                                         data.valid = data.valid, 
#                                                         t.eval = t.eval, 
#                                                         variables.vec = variables.vec,
#                                                         copula.in = "gumbel", 
#                                                         rotate.in = 270)


###
### Frank
###
print("fit frank r0")
Sys.time()
predrisk.frank <- calc.predrisk.copula(data.devel = data.devel, 
                                       data.valid = data.valid, 
                                       t.eval = t.eval, 
                                       variables.vec = variables.vec,
                                       copula.in = "frank", 
                                       rotate.in = 0)

save.image(paste("Project 4/data/clinical_example_fit_models", model.type, "_n", data.size,".RData", sep = ""))


