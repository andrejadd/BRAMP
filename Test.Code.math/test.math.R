
source("../Code/ginv.R")
source("../Code/mvrnorm.R")


n.obs = 5

alpha.sigma = 1
beta.sigma = 0.1

delta.snr = 0.5

##
## This if for testing to see if my code does the same as Marcos
##

total.parents = 11  # 9 + 2 additional

parents.vec = c(1,2,3)

## init a [parents,parents] matrix
#Cov_mat = matrix(c(1:parents))


X = matrix( c(1.0, 1.1, 2.1,
              1.0, 1.2, 2.2,
              1.0, 1.3, 2.3,
              1.0, 1.4, 2.4,
              1.0, 1.5, 2.5), nrow=3, ncol=5)

y = c(2.1,2.2,2.3,2.4,2.5)

mu.vec = c(1,11,21)

mu.tilde = t(X) %*% mu.vec

Cov.mat = matrix( c(1,1.1,21, 2, 1.8, 27.8, 3, 1.3, 23), nrow=3, ncol=3)
  
inv.Sigma.tilde = diag(n.obs) - t(X) %*% ginv(ginv(delta.snr * Cov.mat) + (X %*% t(X)) ) %*% X

mahalanobis.dist = t(y - mu.tilde) %*% inv.Sigma.tilde %*% (y - mu.tilde)


alpha.sigma = alpha.sigma + (n.obs/2)
        
beta.sigma = beta.sigma + (mahalanobis.dist/2)


inv.var.all = rgamma(1,shape=alpha.sigma, scale=(1/beta.sigma))
variance = 1/inv.var.all


## mean and covariance for weight
Sigma.inv.star = ginv(delta.snr * Cov.mat) + X %*% t(X)
mu.star = ginv(Sigma.inv.star) %*% (ginv(delta.snr * Cov.mat) %*% mu.vec + X %*% y)

browser()
weights = mvrnorm(1, mu=mu.star, Sigma= (variance * ginv(Sigma.inv.star)))

