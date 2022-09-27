# Claas Heuer, September 2015
#
# Multinomial class probabilities given pairwise probabilities
# from bianry calssifiers.
# Reference: Probability Estimates for Multi-class Classification by Pairwise Coupling, Wu et al. (2003)

# number of classes
K <- 4

# some random probabilities from binary contrasts
probs <- array(rbeta(K,1,1), dim=c(1,(K*(K-1)/2)))
#probs <- array(c(0.6,0.8,0.2, 0.5, 0.9, 0.1), dim=c(1,(K*(K-1)/2)))
#colnames(probs) <- c("1 vs 2", "1 vs 3", "1 vs 4", "2 vs 3", "2 vs 4", "3 vs 4")

# library(xtable)
#print(xtable(probs), include.rownames=FALSE)
str(probs)
# the Q matrix from Wu et al. (2003)
Q <- matrix(0,K,K)
Q[lower.tri(Q)] <- 1 - probs
Qt <- t(Q)
Q[upper.tri(Q)] <- 1 - Qt[upper.tri(Qt)]
diag(Q) <- rowSums(Q)
Q <- Q / (K-1)


#x=xtable(Q,align=rep("",ncol(Q)+1))
#print(x, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)


# initial vector
p <- rbeta(K,1,1)
p <- p/sum(p)

# updating the prob vector until equilibrium is reached
for(i in 1:1000) p <- Q%*%p

round(t(p), digits=2)


https://gist.github.com/cheuerde/7c649749892c8623eee2#file-pairwise_coupling-r