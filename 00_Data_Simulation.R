#########################################################################################
#===============================SYNTHETIC DATA FUNCTIONS================================#
#########################################################################################



#########################################################################################
#=================================SIMPLE DATA FUNCTIONS=================================#
#########################################################################################
# The output of each function contains:
# [[1]] <- simulated data
# [[2]] <- true values
# [[3]] <- regression parameters obtained with the full data
# [[4]] <- standard deviation of the regression parameters estimated with the full dataset
#==========================SIMULATE DATA FOR LINEAR REGRESSION==========================

# LINEAR RELATIONSHIP
# Y ~ Normal( mean = Intercept + theta_1*x1 + theta_2*x2, variance = sigma.e)
# N: length of the data set
# Intercept/theta_0
# theta_1: main effect of x1 on Y
# theta_2: main effect of x2 on Y
# x1 and x2  are continuous

sim_LinR<-function(N, theta.0, theta.1, theta.2, sigma.e){

  x1<-runif(n = N, -1, 1)
  x2<-runif(n = N, -1, 1)
  y<-rnorm(N , mean = theta.0 + theta.1*x1 + theta.2*x2, sd = sigma.e)
  
  df<-data.frame(y , x1, x2)
  
  theta<- c(theta.0, theta.1, theta.2)
  
  mod<- lm(y ~ . ,data = df)
  
  theta.hat<-coefficients(mod)
  theta.hat.SE<- sqrt(diag(vcov(mod)))
  
  sim_out<-list(
    sim.data = df,
    theta.true = theta,
    theta.hat = theta.hat,
    theta.hat.SE = theta.hat.SE
  )
  
  return(sim_out)
}

###########Simulation LinR###########

set.seed(123)
df_LinR<-sim_LinR(N = 100000, theta.0 = 1, theta.1 = 2, theta.2 = 3, sigma.e = 1)
#=========================================================================================#
#==========================SIMULATE DATA FOR LOGISTIC REGRESSION ==========================

# BINOMIAL 
# Y ~ Binom(N, p = 1/(1+exp(-(Intercept + theta_1*x1 + theta_2*x2)) )
# N: length of the data set
# Intercept/theta_0 is the plogis of the mean probability
# theta_1: main effect of x1 on Y
# theta_2: main effect of x2 on Y
# x1 is Continuous 
# x2 is Discrete with 3 leves

sim_LogR<-function(N, theta.1, theta.2B, theta.2C, mean.prob) {

  
  theta.0<- plogis(mean.prob)
  categories.prob<-c(0.2,0.5,0.2)
  
  categories<- sample(c("A","B","C"), size = N, replace=TRUE, prob=categories.prob)
  theta.2<-0
  for (i in 1:N) {
    theta.2[i]<-ifelse(categories[i]=="A", theta.0,
                       ifelse(categories[i]=="B", theta.2B+theta.0, theta.2C+theta.0))
  }
  
  x1<-runif(n = N, -1, 1)
  
  logit.p<- theta.1*x1 +theta.2
  p<- plogis(logit.p)
  
  y<-rbinom(N,1,p)

  
  df<-data.frame(y , x1, x2 = categories)
  
  theta<- c(theta.0, theta.1, theta.2B, theta.2C)
  
  mod<- glm(y ~ x1 + x2 ,data = df, family = binomial(link = "logit"))

  
  
  theta.hat<-coefficients(mod)
  theta.hat.SE<-sqrt(diag(vcov(mod)))
  
  sim_out<-list(
    sim.data = df,
    theta.true = theta,
    theta.hat = theta.hat,
    theta.hat.SE = theta.hat.SE
  )
  
  return(sim_out)
  
}

###########Simulation LogR###########

set.seed(456)
df_LogR<-sim_LogR(N = 100000, mean.prob = 0.65, theta.1 = -0.3, theta.2B = -0.4, theta.2C = 0.7)
#==========================SIMULATE DATA FOR NON LINEAR REGRESSION==========================


# NON LINEAR RELATIONSHIP
# Michaelis-Menten equation for species richness
# Y = (theta.1*x / theta.2 + x) + error
# error ~ Normal(mean = 0, sd = sigma.e)
# N: length of the data set
# theta_1: main effect of x on Y
# theta_2: main effect of 1 on Y
# max.x: max number for the independent variable x
# x is discrete and positive

sim_NlinR<-function(N, theta.1, theta.2, max.x, sigma.e){


  
  x<-round(runif(n = N, 0, max.x))
  error<-rnorm(N, mean = 0, sd = sigma.e)

  y<- ((theta.1*x)/(theta.2 + x)) + error
  
  
  df<-data.frame(y , x)
  
  theta<- c(theta.1, theta.2)
  
  mod<- nls(y ~ ((a*x)/(b + x)), start = c(a = 1, b = 1),data = df)
  

  
  theta.hat<-coefficients(mod)
  theta.hat.SE<-sqrt(diag(vcov(mod)))
  
  sim_out<-list(
    sim.data = df,
    theta.true = theta,
    theta.hat = theta.hat,
    theta.hat.SE = theta.hat.SE
  )
  
  return(sim_out)
}

###########Simulation NLinR###########

set.seed(789)

df_NLinR<-sim_NlinR(N = 100000, theta.1 = 50, theta.2 = 2, max.x = 10, sigma.e = 1)


#########################################################################################
#================================GROUPED DATA FUNCTIONS=================================#
library(lme4)
#########################################################################################
# The output of each function contains:
# [[1]] <- simulated data
# [[2]] <- true values
# [[3]] <- regression parameters obtained with the full data
# [[4]] <- standard deviation of the regression parameters estimated with the full dataset
# [[5]] <-
# [[6]] <-
#=================SIMULATE DATA FOR NON LINEAR MIX-EFFECT MODEL BALANCED=========== ======

# LINEAR RELATIONSHIP WITH ONE LEVEL OF GROUPING
# Y = x1 + x2 +(1|grp) 
# error ~ Normal(mean = 0, sd = sigma.e)
# re ~ Normal(mean = 0, sd = ranef.sd)
# N: length of the data set N.grp*N.Obs
# N.grp : Number of groups
# N.Obs : Number of observations per group
# theta_1: main effect of x1 on Y
# theta_2: main effect of x2 on Y
# x1 and x2 are continuous

sim_LMM1<-function(N.grp, N.Obs, theta.0, theta.1, theta.2, sigma.e ,ranef.sd){
  
  N<-N.grp*N.Obs

  theta = c(theta.0, theta.1, theta.2)

  x1<-runif(n = N, -1, 1)
  x2<-runif(n = N, -1, 1)
  
  error<-rnorm(N , mean = 0, sd = sigma.e) 
  re<-rnorm(N.grp, mean = 0, sd = ranef.sd)
  
  grp<- rep(LETTERS[1:N.grp], each = N.Obs)
  
  grp.eff<-rep(re, each = N.Obs)
  
  y <- theta.0 + theta.1*x1 + theta.2*x2 + error + grp.eff
  
  
  df<- data.frame(y, x1, x2, grp)


  mod<- lmer(y ~ x1+ x2 + (1|grp), data = df)
  
  
  ranef.hat<-ranef(mod)[[1]]
  
  theta.hat<-fixef(mod)
  theta.hat.SE<-sqrt(diag(vcov(mod)))
  
  var.matrix<-as.data.frame(VarCorr(mod))
  
  sim_out<- list(sim.data = df,
                 theta.true = theta,
                 theta.hat = theta.hat,
                 theta.hat.SE = theta.hat.SE,
                 ranef.true = re,
                 ranef.hat = ranef.hat,
                 var.re = var.matrix[,ncol(var.matrix)])
  
  
  return(sim_out)
}

###########Simulation LMM1###########

set.seed(321)
df_LMM1<- sim_LMM1(N.grp = 10, N.Obs = 10000, theta.0 = 1, theta.1 = 2, theta.2 = 3, sigma.e =  1, ranef.sd = 2)
#==========================SIMULATE DATA FOR NON LINEAR MIX-EFFECT MODEL UNBALANCED==========================

# LINEAR RELATIONSHIP WITH ONE LEVEL OF GROUPING
# Y = x1 + x2 +(1|grp) 
# error ~ Normal(mean = 0, sd = sigma.e)
# re ~ Normal(mean = 0, sd = ranef.sd)
# N: length of the data set N.grp*N.Obs
# N.grp : Number of groups
# N.Obs : Number of observations per group
# theta_1: main effect of x1 on Y
# theta_2: main effect of x2 on Y
# x1 and x2 are continuous

sim_LMM2<-function(N , N.grp, N.Obs, theta.0, theta.1, theta.2, sigma.e, ranef.sd ){

  
  theta = c(theta.0, theta.1, theta.2)
  
  x1<-runif(n = N.grp*N.Obs, -1, 1)
  x2<-runif(n = N.grp*N.Obs, -1, 1)
  
  error<-rnorm(N.grp*N.Obs , mean = 0, sd = sigma.e ) 
  re<-rnorm(N.grp, mean = 0, sd = ranef.sd)
  
  grp<- rep(LETTERS[1:N.grp], each = N.Obs)
  grp.eff<- rep(re, each = N.Obs)
  
  y <- theta.0 + theta.1*x1 + theta.2*x2 + error + grp.eff
  
  
  df.full<- data.frame(y, x1, x2, grp)
  
  length(seq(from = 0, to = 0.5, by = 0.05)[-1])
  
  #SAMPLE A DIFFERENT NUMBER OF GROUPS
  grp.prob<-c(rep(seq(from = 0, to = 0.5, by = 0.05)[-1], each = N.Obs))
  
  df<-df.full[sample(N.grp*N.Obs , N , replace = F, prob = grp.prob),]
  
  mod<- lmer(y ~ x1+ x2 + (1|grp), data = df)
  
  ranef.hat<-ranef(mod)[[1]]
  
  theta.hat<-fixef(mod)
  theta.hat.SE<-sqrt(diag(vcov(mod)))
  
  var.matrix<-as.data.frame(VarCorr(mod))
  
  sim_out<- list(sim.data = df,
                 theta.true = theta,
                 theta.hat = theta.hat,
                 theta.hat.SE = theta.hat.SE,
                 ranef.true = re,
                 ranef.hat = ranef.hat,
                 var.re = var.matrix[,ncol(var.matrix)])
  
  
  return(sim_out)
}


###########Simulation LMM1###########
set.seed(654)

df_LMM2<-sim_LMM2(N = 100000, N.grp = 10, N.Obs = 20000, theta.0 = 1, theta.1 = 2, theta.2 = 3, sigma.e = 1, ranef.sd = 2)
#==========================SIMULATE DATA FOR NON LINEAR MIX-EFFECT MODEL RANDOM SLOPE==========================

# LINEAR RELATIONSHIP WITH ONE LEVEL OF GROUPING
# Y = x1 + x2 +(0+x2|grp) 
# error ~ Normal(mean = 0, sd = sigma.e)
# re ~ Normal(mean = 0, sd = ranef.sd)
# N: length of the data set N.grp*N.Obs
# N.grp : Number of groups
# N.Obs : Number of observations per group
# theta_1: main effect of x1 on Y
# theta_2: main effect of x2 on Y
# x1 and x2 are continuous


sim_LMM3<-function(N.grp, N.Obs, theta.0, theta.1, theta.2, sigma.e ,ranef.sd){

  N<-N.Obs*N.grp
  
  theta = c(theta.0, theta.1, theta.2)
  
  
  x1<-runif(n = N, -1, 1)
  x2<-runif(n = N, -1, 1)
  
  grp<-rep(LETTERS[1:N.grp], each = N.Obs)
  
  
  error<-rnorm(N , mean = 0, sd = sigma.e) 
  
  re<-rnorm(N.grp, mean = 0, sd = ranef.sd)
  
  grp.eff<- rep(re, each =  N.Obs)
  
  
  y <- theta.0 + theta.1*x1 +(theta.2 + grp.eff)*x2 + error
  
  df<- data.frame(y, x1, x2, grp)
  
  
  mod<- lmer(y ~ x1+ x2 + (0+x2|grp), data = df)
  
  
  ranef.hat<-ranef(mod)[[1]]
  
  theta.hat<-fixef(mod)
  theta.hat.SE<-sqrt(diag(vcov(mod)))
  
  var.matrix<-as.data.frame(VarCorr(mod))
  
  
  sim_out<- list(sim.data = df,
                 theta.true = theta,
                 theta.hat = theta.hat,
                 theta.hat.SE = theta.hat.SE,
                 ranef.true = re,
                 ranef.hat = ranef.hat,
                 var.re = var.matrix[,ncol(var.matrix)])
  
  
  return(sim_out)
}
###########Simulation LMM1###########
set.seed(987)

df_LMM3<- sim_LMM3(N.grp = 10, N.Obs = 10000, theta.0 = 1, theta.1 = 2, theta.2 = 3, sigma.e = 1, ranef.sd = 2.5)



#save(df_LinR, df_LogR, df_NLinR , df_LMM1, df_LMM2, df_LMM3 , file = "synthetic_data.RData")
