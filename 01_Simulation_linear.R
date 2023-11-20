#############################################################################
###########This code is used for continuous outcome variable#################
#############################################################################
#n=sample siz
#alpha = the effect of exposure on mediator
#eta1 = the effect of instrument of exposure on exposure
#eta2 = the effect of instrument of mediator on mediator
#Tot = the total effect of exposure on outcome
#proportion= the percentage of the mediation effect in relation to the overall effect
#theta = the effect of interaction between exposure and mediator
#NIE = natural indirect effect
#NDE = natural direct effect
#X = exposure
#M = mediator
#Y = outcome

library(foreach);library(doParallel);library(parallel)

#assign specific values to the parameters
n <- 50000
alpha <- 0.2
R_square <- 0.1
Tot <- c(0.2, 0.5, 1)
proportion <- c(0.05, 0.25, 0.75)
theta <- c(-0.5, -0.2, 0, 0.2, 0.5)

#derive eta1 and eta2 based on the values of alpha and R-squared.
eta1 <- sqrt(2*R_square/(1-R_square))
eta2 <- sqrt(R_square*(alpha^2*eta1^2+2*alpha^2+2*alpha+2)/(1-R_square))
F.statistics <- (n-1-1)/1*(R_square)/(1-R_square)

#set the scenarios
senario <- matrix(nrow = 45,ncol = 11,dimnames = list(paste('senario_',1:45,sep=''),c('Total','Proportion','alpha','theta','beta','gamma','eta1','eta2','NIE','NDE','n')))
i <- 1  
for (j in Tot) {for (k in proportion) {for (m in theta) {
        senario[i,1:4] <- c(j,k,alpha,m)
        senario[i,5] <- j*(1-k)
        senario[i,6] <- (j-m*alpha-senario[i,5])/alpha
        senario[i,7:8] <- c(eta1,eta2)
        senario[i,9] <- j*k
        senario[i,10] <- j*(1-k)
        senario[i,11] <- n
        i <- i+1
}}}
sim.fun <- function(n=n,alpha=alpha,beta=beta,gamma=gamma,theta=theta){
  #simulate the data
  U <- rnorm(mean = 0,sd = 1,n = n)
  V <- rnorm(mean = 0,sd = 1,n = n)
  E <- rnorm(mean = 0,sd = 1,n = n)
  C <- rnorm(mean = 0,sd = 1,n = n)
  Z1 <- rnorm(mean = 0,sd = 1,n = n)
  Z2 <- rnorm(mean = 0,sd = 1,n = n)
  X <- C+eta1*Z1+V 
  M <- C+E+eta2*Z2+alpha*X
  Y <- C+U+beta*X+gamma*M+theta*X*M
  #regression models
  exposure_error <- residuals(lm(X~Z1))
  mediator_error <- residuals(lm(M~Z2+X+exposure_error))
  exposure_mediator_error <- residuals(lm(X*M~Z1*Z2+I(Z1^2)))
  exposure_mediator_model <- lm(M~X+exposure_error)
  model.iv <- lm(Y~X*M+exposure_error+mediator_error+exposure_mediator_error)
  exposure_mediator.vand.model <- lm(M~X+C)
  van.model <- lm(Y~X*M+C)
  res <- data.frame('NIE'=gamma*alpha+theta*alpha,'NIE.iv'=coef(exposure_mediator_model)[2]*(coef(model.iv)[3]+coef(model.iv)[7]),
                    'NIE.vand'=coef(exposure_mediator.vand.model)[2]*(coef(van.model)[3]+coef(van.model)[5]),
                    'alpha'=alpha,'alpha.iv'=coef(exposure_mediator_model)[2],
                    'beta'=beta,'beta.iv'=coef(model.iv)[2],
                    'gamma'=gamma,'gamma.iv'=coef(model.iv)[3],
                    'theta'=theta,'theta.iv'=coef(model.iv)[7])
  return(res)
}
set.seed(1200)
clus <- makeCluster(20)
registerDoParallel(clus)
sim.res.linear <- foreach(i=1:nrow(senario)) %:%
  foreach(1:5000,.combine='rbind') %dopar% 
  sim.fun(n=senario[i,'n'],alpha=senario[i,'alpha'],beta=senario[i,'beta'],gamma=senario[i,'gamma'],theta=senario[i,'theta'])
stopCluster(clus)
names(sim.res.linear) <- paste('senario_',1:45,sep='')
save(sim.res.linear,file='Simulation_linear.RData')