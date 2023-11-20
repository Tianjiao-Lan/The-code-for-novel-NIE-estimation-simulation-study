#############################################################################
####This code is used for dichotomy outcome variable with weak instrument####
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

library(foreach);library(doParallel);library(parallel);library(boot)

#assign specific values to the parameters
n <- 50000
alpha <- 0.2
F.statistics <- 5
Tot <- exp(c(0.2, 0.5, 1))
proportion <- c(0.05, 0.25, 0.75)
theta <- c(-0.5, -0.2, 0, 0.2, 0.5)

#determine R-squared based on the F statistics
R_square <- F.statistics*2/(F.statistics*2+n-2-1)
#determine eta1 and eta2 based on the alpha and R-squared
eta1 <- sqrt(2*R_square/(1-R_square))
eta2 <- sqrt(R_square*(alpha^2*eta1^2+2*alpha^2+2*alpha+2)/(1-R_square))

#set the scenarios
senario <- matrix(nrow = 45,ncol = 11,dimnames = list(paste('senario_',1:45,sep=''),c('Total','Proportion','alpha','theta','beta','gamma','eta1','eta2','NIE','NDE','n')))
i <- 1  
for (j in Tot) {for (k in proportion) {for (m in theta) {
  senario[i,1:4] <- c(j,k,alpha,m)
  senario[i,10] <- j-k*(j-1) #OR_NDE=OR_Total-(Proportion)*(OR_Total-1)
  senario[i,9] <- (senario[i,10]-k)/(1-k)/(senario[i,10]) #OR_NIE=(OR_NDE-Proportion)/((1-Proportion)*OR_NDE)
  senario[i,6] <- (log(senario[i,9])-m*alpha)/alpha #gamma=(log(OR_NIE)-theta*alpha)/alpha
  senario[i,5] <- log(senario[i,10])-m*alpha-m*senario[i,6]-0.5*m^2 #beta=log(OR_NDE)-theta*alpha-theta*gamma-0.5*theta^2
  senario[i,7:8] <- c(eta1,eta2)
  senario[i,11] <- n
  i <- i+1
}}}

sim.fun <- function(n=n,alpha=alpha,beta=beta,gamma=gamma,theta=theta){
  E <- rnorm(mean = 0,sd = 1,n = n)
  V <- rnorm(mean = 0,sd = 1,n = n)
  C <- rnorm(mean = 0,sd = 1,n = n)
  Z1 <- rnorm(mean = 0,sd = 1,n = n)
  Z2 <- rnorm(mean = 0,sd = 1,n = n)
  X <- C+eta1*Z1+V
  M <- C+E+eta2*Z2+alpha*X
  P_Y.indi <- boot::inv.logit(C+beta*X+gamma*M+theta*X*M)
  #due to computational errors arising from computer floating-point numbers, restricting the values of P to 1 and 0 facilitates the identification of the intercept beta_0.
  P_Y.indi[which(P_Y.indi==1)] <- 0.99^20
  P_Y.indi[which(P_Y.indi==0)] <- 0.01^20
  rate <- 0.1
  #find the intercept beta_0 such that it satisfies the occurrence rate of Y as set
  beta0.fun <- function(x){return(mean(P_Y.indi*exp(x)/(P_Y.indi*exp(x)+1-P_Y.indi))-rate)}
  beta_0 <- uniroot(beta0.fun,interval = c(-500,500))$root
  P_Y <- boot::inv.logit(beta_0+C+beta*X+gamma*M+theta*X*M)
  Y <- rbinom(n = n,size = 1,prob = P_Y)
  #regression models
  exposure_error <- residuals(lm(X~Z1))
  mediator_error <- residuals(lm(M~Z2+X+exposure_error))
  exposure_mediator_error <- residuals(lm(X*M~Z1*Z2+I(Z1^2)))
  exposure_mediator_model <- lm(M~X+exposure_error)
  model.iv <- glm(Y~X*M+exposure_error+mediator_error+exposure_mediator_error,family=binomial(logit))
  exposure_mediator.vand.model <- lm(M~X+C)
  van.model <- glm(Y~X*M+C,family=binomial(logit))
  alpha.iv <- coef(exposure_mediator_model)[2]
  gamma.iv <- coef(model.iv)[3]
  theta.iv <- coef(model.iv)[7]
  res <- data.frame('P.y'=mean(Y),'NIE'=exp(gamma*alpha+theta*alpha),'NIE.iv'=exp(alpha.iv*(gamma.iv+theta.iv)),
                    'NIE.vand'=exp(coef(exposure_mediator.vand.model)[2]*(coef(van.model)[3]+coef(van.model)[5])),
                    'alpha'=alpha,'alpha.iv'=alpha.iv,
                    'beta'=beta,'beta.iv'=coef(model.iv)[2],
                    'gamma'=gamma,'gamma.iv'=gamma.iv,
                    'theta'=theta,'theta.iv'=theta.iv)
  return(res)
}
set.seed(1200)
clus <- makeCluster(20)
registerDoParallel(clus)
sim.res.nonlinear.weakiv <- foreach(i=1:nrow(senario)) %:%
  foreach(1:5000,.combine='rbind') %dopar% 
  sim.fun(n=senario[i,'n'],alpha=senario[i,'alpha'],beta=senario[i,'beta'],gamma=senario[i,'gamma'],theta=senario[i,'theta'])
stopCluster(clus)
names(sim.res.nonlinear.weakiv) <- paste('senario_',1:45,sep='')
save(sim.res.nonlinear.weakiv,file='Simulation_nonlinear_weakIV.RData')