####### All covariates for each scenario in one place, so I have easy access to them

#######################
### DGM MULTINOMIAL ###
#######################

### Scenario 1, K = 5
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))
beta0.sim <- c(-1, -1, -1, -1)*0.35

### Scenario 2, K = 5
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(1, 1, 1, 1)*0.1

### Scenario 3, K = 5
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))
beta0.sim <- c(-1, -1, -1, -1)*0.35

### Scenario 4, K = 5
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(-1, -1, -1, -1)*0.35

### Scenario 1, K = 3
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))
beta0.sim <- c(-1, -1)*0.35

### Scenario 2, K = 3
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(1, 1)*0.1

### Scenario 3, K = 3
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))
beta0.sim <- c(-1, -1)*0.35

### Scenario 4, K = 3
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(-1, -1)*0.35


##################
### DGM SEQLOG ###
##################

### Scenario 1, K = 5
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))
beta0.sim <- c(-0.75, -0.5, -0.25, 0)*2.25

### Scenario 2, K = 5
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(-0.75, -0.5, -0.25, 0)*2.25

### Scenario 3, K = 5
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))
beta0.sim <- c(-0.75, -0.5, -0.25, 0)*2.25

### Scenario 4, K = 5
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(-0.75, -0.5, -0.25, 0)*2.25

### Scenario 1, K = 3
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))
beta0.sim <- c(-0.5, 0)*2.25

### Scenario 2, K = 3
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(-0.5, 0)*2.25

### Scenario 3, K = 3
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))
beta0.sim <- c(-0.5, 0)*2.25

### Scenario 4, K = 3
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(-0.5, 0)*2.25