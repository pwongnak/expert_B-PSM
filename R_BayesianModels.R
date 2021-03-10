##########################################################################################
#Bayesian parametric survival model
#Bayesian models in JAGS 1) Without expert; 2) Mixing experts; 3) Ageraging experts
##########################################################################################
library(rjags)
#Tick survival data###################################################################
survData <- read.table("table__survival_data_Milne1949_NEW.csv", sep = ",", header = T)
survData$RH_n <- survData$RH_n/100
uTemp <- unique(survData$temp)
uRH <- unique(survData$RH_n)
survData$condition <- as.factor(paste(survData$temp, survData$RH_n, sep = "-"))
condData <- survData[,-1]
conData2 <-unique(condData)
J <- match(survData$condition, conData2$condition) 
Q <- conData2$temp
U <- conData2$RH_n
I <- nrow(survData) 
ncond <- nrow(conData2)

#Model 1 Without experts#############################################################
M1 <- "model{
 #Likelihood##################################################################################
   for(i in 1:I){
     Time[i] ~ dweib(p ,lambda[J[i]])
   }
   
   for(j in 1:ncond){
      logLambda[j] <- (beta0 + (beta1 * (U[j])^k) + (beta2 * (Q[j])) + (beta3 * (A[j]) * (U[j])^k))
      lambda[j] <- exp(logLambda[j])
   }
   
   beta0 ~ dunif(-50, 50)
   beta1 ~ dunif(-50, 50)
   beta2 ~ dunif(-50, 50)
   beta3 ~ dunif(-50, 50)
 k ~ dunif(1, 6)
 p ~ dunif(0.01, 5)
}"
data_M1 <- list(I = I, Time = survData$survivalTime, ncond = ncond, Q = Q, U = U, J = J)
params_M1 <- c("beta0", "beta1", "beta2", "beta3", "p", "k")
model1 <- jags.model(textConnection(M1), data = data_M1, n.chains = 3)
update(model1, 5000*10)
mcmc1 <- coda.samples(model1, variable.names = params_M1, n.iter = 5000*700 , thin = 700)
plot(mcmc1)
post_M1 <- as.data.frame(as.matrix(mcmc1))
summary(post_M1)
traceplot(mcmc1)
gelman.plot(mcmc1)
autocorr.plot(mcmc1)
write.csv(post_M1, "post_mod1new.csv", row.names = F)
#Expert data#########################################################################
ExpertList <- c("ExpertData_1.csv", "ExpertData_2.csv", "ExpertData_3.csv", 
                "ExpertData_4.csv", "ExpertData_5.csv", "ExpertData_6.csv")
N <- length(ExpertList)
rep1Ex <- rep(1, N)
pi <- rep(1/N, N)
expData <- NULL
for(i in 1:N){
  data <- read.csv(ExpertList[i])
  sdata <- c(data$Tp[1], data$U[1], data$meanLog[1], data$sdLog[1],
             data$Tp[2], data$U[2], data$meanLog[2], data$sdLog[2],
             data$Tp[3], data$U[3], data$meanLog[3], data$sdLog[3],
             data$Tp[4], data$U[4], data$meanLog[4], data$sdLog[4])
  expData <- rbind(sdata, expData)
}

#M2: Pooling expert opinions##############################################################################################################################
M2 <- 'model{
 #Likelihood##################################################################################
   for(i in 1:I){
     Time[i] ~ dweib(p ,lambda[J[i]])
   }
   
   for(j in 1:ncond){
      logLambda[j] <- (beta0 + (beta1 * (U[j])^k) + (beta2 * (Q[j])) + (beta3 * (Q[j]) * (U[j])^k))
      lambda[j] <- exp(logLambda[j])
   }
 
 for(ne in 1:N){
  #Define X matrix
  X[1:4,1, ne] <- c(1,1,1,1)
  X[1:4,2, ne] <- expData[ne, c(2,6,10,14)]^k                             #U^k
  X[1:4,3, ne] <- expData[ne, c(1,5,9,13)]                                #Tp
  X[1:4,4, ne] <- expData[ne, c(1,5,9,13)] * expData[ne, c(2,6,10,14)]^k  #interaction
 
 #Determinant for 4 x 4 matrix##########################################################
  #S1#########################################
  pS1[ne] <- X[4,1,ne]
  S1[1:3,1:3, ne] <- X[1:3, 2:4,ne]
  DetS1[ne] <- (S1[1,1,ne] * S1[2,2,ne] * S1[3,3,ne]) + (S1[1,2,ne] * S1[2,3,ne] * S1[3,1,ne]) + (S1[1,3,ne] * S1[2,1,ne] * S1[3,2,ne]) -
           (S1[1,3,ne] * S1[2,2,ne] * S1[3,1,ne]) - (S1[1,1,ne] * S1[2,3,ne] * S1[3,2,ne]) - (S1[1,2,ne] * S1[2,1,ne] * S1[3,3,ne])
  #S2#########################################
  pS2[ne] <- X[4,2,ne]
  S2[1:3,1:3, ne] <- X[1:3, c(1,3,4),ne]
  DetS2[ne] <- (S2[1,1,ne] * S2[2,2,ne] * S2[3,3,ne]) + (S2[1,2,ne] * S2[2,3,ne] * S2[3,1,ne]) + (S2[1,3,ne] * S2[2,1,ne] * S2[3,2,ne]) -
           (S2[1,3,ne] * S2[2,2,ne] * S2[3,1,ne]) - (S2[1,1,ne] * S2[2,3,ne] * S2[3,2,ne]) - (S2[1,2,ne] * S2[2,1,ne] * S2[3,3,ne])
  #S3#########################################
  pS3[ne] <- X[4,3,ne]
  S3[1:3,1:3, ne] <- X[1:3, c(1,2,4),ne]
  DetS3[ne] <- (S3[1,1,ne] * S3[2,2,ne] * S3[3,3,ne]) + (S3[1,2,ne] * S3[2,3,ne] * S3[3,1,ne]) + (S3[1,3,ne] * S3[2,1,ne] * S3[3,2,ne]) -
           (S3[1,3,ne] * S3[2,2,ne] * S3[3,1,ne]) - (S3[1,1,ne] * S3[2,3,ne] * S3[3,2,ne]) - (S3[1,2,ne] * S3[2,1,ne] * S3[3,3,ne])
  #S4#########################################
  pS4[ne] <- X[4,4,ne]
  S4[1:3,1:3, ne] <- X[1:3, c(1,2,3),ne]
  DetS4[ne] <- (S4[1,1,ne] * S4[2,2,ne] * S4[3,3,ne]) + (S4[1,2,ne] * S4[2,3,ne] * S4[3,1,ne]) + (S4[1,3,ne] * S4[2,1,ne] * S4[3,2,ne]) -
           (S4[1,3,ne] * S4[2,2,ne] * S4[3,1,ne]) - (S4[1,1,ne] * S4[2,3,ne] * S4[3,2,ne]) - (S4[1,2,ne] * S4[2,1,ne] * S4[3,3,ne])

  DetX[ne] <-  (-1 * pS1[ne] * DetS1[ne]) + (pS2[ne] * DetS2[ne]) - (pS3[ne] * DetS3[ne]) + (pS4[ne] * DetS4[ne])
 
  #Adjugate#########################################################################
  A[1,1,ne] <- (X[2,2,ne]*X[3,3,ne]*X[4,4,ne]) + (X[2,3,ne]*X[3,4,ne]*X[4,2,ne]) + (X[2,4,ne]*X[3,2,ne]*X[4,3,ne]) -
            (X[2,4,ne]*X[3,3,ne]*X[4,2,ne]) - (X[2,3,ne]*X[3,2,ne]*X[4,4,ne]) - (X[2,2,ne]*X[3,4,ne]*X[4,3,ne])
  A[1,2,ne] <- (-1*X[1,2,ne]*X[3,3,ne]*X[4,4,ne]) - (X[1,3,ne]*X[3,4,ne]*X[4,2,ne]) - (X[1,4,ne]*X[3,2,ne]*X[4,3,ne]) +
            (X[1,4,ne]*X[3,3,ne]*X[4,2,ne]) + (X[1,3,ne]*X[3,2,ne]*X[4,4,ne]) + (X[1,2,ne]*X[3,4,ne]*X[4,3,ne])
  A[1,3,ne] <- (X[1,2,ne]*X[2,3,ne]*X[4,4,ne]) + (X[1,3,ne]*X[2,4,ne]*X[4,2,ne]) + (X[1,4,ne]*X[2,2,ne]*X[4,3,ne]) -
            (X[1,4,ne]*X[2,3,ne]*X[4,2,ne]) - (X[1,3,ne]*X[2,2,ne]*X[4,4,ne]) - (X[1,2,ne]*X[2,4,ne]*X[4,3,ne])
  A[1,4,ne] <- (-1*X[1,2,ne]*X[2,3,ne]*X[3,4,ne]) - (X[1,3,ne]*X[2,4,ne]*X[3,2,ne]) - (X[1,4,ne]*X[2,2,ne]*X[3,3,ne]) +
            (X[1,4,ne]*X[2,3,ne]*X[3,2,ne]) + (X[1,3,ne]*X[2,2,ne]*X[3,4,ne]) + (X[1,2,ne]*X[2,4,ne]*X[3,3,ne])

  A[2,1,ne] <- (-1*X[2,1,ne]*X[3,3,ne]*X[4,4,ne]) - (X[2,3,ne]*X[3,4,ne]*X[4,1,ne]) - (X[2,4,ne]*X[3,1,ne]*X[4,3,ne]) +
            (X[2,4,ne]*X[3,3,ne]*X[4,1,ne]) + (X[2,3,ne]*X[3,1,ne]*X[4,4,ne]) + (X[2,1,ne]*X[3,4,ne]*X[4,3,ne])
  A[2,2,ne] <- (X[1,1,ne]*X[3,3,ne]*X[4,4,ne]) + (X[1,3,ne]*X[3,4,ne]*X[4,1,ne]) + (X[1,4,ne]*X[3,1,ne]*X[4,3,ne]) -
            (X[1,4,ne]*X[3,3,ne]*X[4,1,ne]) - (X[1,3,ne]*X[3,1,ne]*X[4,4,ne]) - (X[1,1,ne]*X[3,4,ne]*X[4,3,ne])
  A[2,3,ne] <- (-1*X[1,1,ne]*X[2,3,ne]*X[4,4,ne]) - (X[1,3,ne]*X[2,4,ne]*X[4,1,ne]) - (X[1,4,ne]*X[2,1,ne]*X[4,3,ne]) +
            (X[1,4,ne]*X[2,3,ne]*X[4,1,ne]) + (X[1,3,ne]*X[2,1,ne]*X[4,4,ne]) + (X[1,1,ne]*X[2,4,ne]*X[4,3,ne])
  A[2,4,ne] <- (X[1,1,ne]*X[2,3,ne]*X[3,4,ne]) + (X[1,3,ne]*X[2,4,ne]*X[3,1,ne]) + (X[1,4,ne]*X[2,1,ne]*X[3,3,ne]) -
            (X[1,4,ne]*X[2,3,ne]*X[3,1,ne]) - (X[1,3,ne]*X[2,1,ne]*X[3,4,ne]) - (X[1,1,ne]*X[2,4,ne]*X[3,3,ne])

  A[3,1,ne] <- (X[2,1,ne]*X[3,2,ne]*X[4,4,ne]) + (X[2,2,ne]*X[3,4,ne]*X[4,1,ne]) + (X[2,4,ne]*X[3,1,ne]*X[4,2,ne]) -
            (X[2,4,ne]*X[3,2,ne]*X[4,1,ne]) - (X[2,2,ne]*X[3,1,ne]*X[4,4,ne]) - (X[2,1,ne]*X[3,4,ne]*X[4,2,ne])
  A[3,2,ne] <- (-1*X[1,1,ne]*X[3,2,ne]*X[4,4,ne]) - (X[1,2,ne]*X[3,4,ne]*X[4,1,ne]) - (X[1,4,ne]*X[3,1,ne]*X[4,2,ne]) +
            (X[1,4,ne]*X[3,2,ne]*X[4,1,ne]) + (X[1,2,ne]*X[3,1,ne]*X[4,4,ne]) + (X[1,1,ne]*X[3,4,ne]*X[4,2,ne])
  A[3,3,ne] <- (X[1,1,ne]*X[2,2,ne]*X[4,4,ne]) + (X[1,2,ne]*X[2,4,ne]*X[4,1,ne]) + (X[1,4,ne]*X[2,1,ne]*X[4,2,ne]) -
            (X[1,4,ne]*X[2,2,ne]*X[4,1,ne]) - (X[1,2,ne]*X[2,1,ne]*X[4,4,ne]) - (X[1,1,ne]*X[2,4,ne]*X[4,2,ne])
  A[3,4,ne] <- (-1*X[1,1,ne]*X[2,2,ne]*X[3,4,ne]) - (X[1,2,ne]*X[2,4,ne]*X[3,1,ne]) - (X[1,4,ne]*X[2,1,ne]*X[3,2,ne]) +
            (X[1,4,ne]*X[2,2,ne]*X[3,1,ne]) + (X[1,2,ne]*X[2,1,ne]*X[3,4,ne]) + (X[1,1,ne]*X[2,4,ne]*X[3,2,ne])

  A[4,1,ne] <- (-1*X[2,1,ne]*X[3,2,ne]*X[4,3,ne]) - (X[2,2,ne]*X[3,3,ne]*X[4,1,ne]) - (X[2,3,ne]*X[3,1,ne]*X[4,2,ne]) +
            (X[2,3,ne]*X[3,2,ne]*X[4,1,ne]) + (X[2,2,ne]*X[3,1,ne]*X[4,3,ne]) + (X[2,1,ne]*X[3,3,ne]*X[4,2,ne])
  A[4,2,ne] <- (X[1,1,ne]*X[3,2,ne]*X[4,3,ne]) + (X[1,2,ne]*X[3,3,ne]*X[4,1,ne]) + (X[1,3,ne]*X[3,1,ne]*X[4,2,ne]) -
            (X[1,3,ne]*X[3,2,ne]*X[4,1,ne]) - (X[1,2,ne]*X[3,1,ne]*X[4,3,ne]) - (X[1,1,ne]*X[3,3,ne]*X[4,2,ne])
  A[4,3,ne] <- (-1*X[1,1,ne]*X[2,2,ne]*X[4,3,ne]) - (X[1,2,ne]*X[2,3,ne]*X[4,1,ne]) - (X[1,3,ne]*X[2,1,ne]*X[4,2,ne]) +
            (X[1,3,ne]*X[2,2,ne]*X[4,1,ne]) + (X[1,2,ne]*X[2,1,ne]*X[4,3,ne]) + (X[1,1,ne]*X[2,3,ne]*X[4,2,ne])
  A[4,4,ne] <- (X[1,1,ne]*X[2,2,ne]*X[3,3,ne]) + (X[1,2,ne]*X[2,3,ne]*X[3,1,ne]) + (X[1,3,ne]*X[2,1,ne]*X[3,2,ne]) -
            (X[1,3,ne]*X[2,2,ne]*X[3,1,ne]) - (X[1,2,ne]*X[2,1,ne]*X[3,3,ne]) - (X[1,1,ne]*X[2,3,ne]*X[3,2,ne])

  X_inv[1:4, 1:4, ne] <- A[1:4, 1:4, ne]/DetX[ne]
  #################################################################################################################
  #Defining Lambda from Tbar and p
  L1[ne] <- -p * (log(Tbar1[ne]) - loggam(1 + (1/p)))
  L2[ne] <- -p * (log(Tbar2[ne]) - loggam(1 + (1/p)))
  L3[ne] <- -p * (log(Tbar3[ne]) - loggam(1 + (1/p)))
  L4[ne] <- -p * (log(Tbar4[ne]) - loggam(1 + (1/p)))
  
  Tbar1[ne] ~ dlnorm(expData[ne,3], 1/(expData[ne,4]^2))
  Tbar2[ne] ~ dlnorm(expData[ne,7], 1/(expData[ne,8]^2))
  Tbar3[ne] ~ dlnorm(expData[ne,11], 1/(expData[ne,12]^2))
  Tbar4[ne] ~ dlnorm(expData[ne,15], 1/(expData[ne,16]^2))
  
  Beta0[ne] <- X_inv[1,1,ne]*L1[ne] + X_inv[1,2,ne]*L2[ne] + X_inv[1,3,ne]*L3[ne] + X_inv[1,4,ne]*L4[ne]
  Beta1[ne] <- X_inv[2,1,ne]*L1[ne] + X_inv[2,2,ne]*L2[ne] + X_inv[2,3,ne]*L3[ne] + X_inv[2,4,ne]*L4[ne]
  Beta2[ne] <- X_inv[3,1,ne]*L1[ne] + X_inv[3,2,ne]*L2[ne] + X_inv[3,3,ne]*L3[ne] + X_inv[3,4,ne]*L4[ne]
  Beta3[ne] <- X_inv[4,1,ne]*L1[ne] + X_inv[4,2,ne]*L2[ne] + X_inv[4,3,ne]*L3[ne] + X_inv[4,4,ne]*L4[ne]
  
 }

 epsilon ~ dcat(pi)

  beta0 <- Beta0[epsilon] 
  beta1 <- Beta1[epsilon]  
  beta2 <- Beta2[epsilon]  
  beta3 <- Beta3[epsilon]   

 k ~ dunif(1, 6)
 p ~ dunif(0.01, 5)
}'

data_M2 <- list(I = I, Time = survData$survivalTime, ncond = ncond, Q = Q, U = U, J = J,
                N = N, expData = expData, pi = pi)
params_M2 <- c("beta0", "beta1", "beta2", "beta3", "p", "k", "epsilon", "Beta0", "Beta1", "Beta2", "Beta3")
model2 <- jags.model(textConnection(M2), data = data_M2, n.chains = 3)
update(model2, 5000*10)
mcmc2 <- coda.samples(model2, variable.names = params_M2, n.iter = 5000*700 , thin = 700)
plot(mcmc2)
post_M2 <- as.data.frame(as.matrix(mcmc2))
summary(post_M2)
traceplot(mcmc2)
gelman.plot(mcmc2)
autocorr.plot(mcmc2)
write.csv(post_M2, "post_mod2new.csv", row.names = F)

#Model 3 Averaging expert opinions############################################
M3 <- 'model{
 #Likelihood##################################################################################
   for(i in 1:I){
     Time[i] ~ dweib(p ,lambda[J[i]])
   }
   
   for(j in 1:ncond){
      logLambda[j] <- (beta0 + (beta1 * (U[j])^k) + (beta2 * (Q[j])) + (beta3 * (Q[j]) * (U[j])^k))
      lambda[j] <- exp(logLambda[j])
   }
 
 for(ne in 1:N){
  #Define X matrix
  X[1:4,1, ne] <- c(1,1,1,1)
  X[1:4,2, ne] <- expData[ne, c(2,6,10,14)]^k                             #U^k
  X[1:4,3, ne] <- expData[ne, c(1,5,9,13)]                                #Tp
  X[1:4,4, ne] <- expData[ne, c(1,5,9,13)] * expData[ne, c(2,6,10,14)]^k  #interaction
 
 #det4*4##########################################################
  #S1#########################################
  pS1[ne] <- X[4,1,ne]
  S1[1:3,1:3, ne] <- X[1:3, 2:4,ne]
  DetS1[ne] <- (S1[1,1,ne] * S1[2,2,ne] * S1[3,3,ne]) + (S1[1,2,ne] * S1[2,3,ne] * S1[3,1,ne]) + (S1[1,3,ne] * S1[2,1,ne] * S1[3,2,ne]) -
           (S1[1,3,ne] * S1[2,2,ne] * S1[3,1,ne]) - (S1[1,1,ne] * S1[2,3,ne] * S1[3,2,ne]) - (S1[1,2,ne] * S1[2,1,ne] * S1[3,3,ne])
  #S2#########################################
  pS2[ne] <- X[4,2,ne]
  S2[1:3,1:3, ne] <- X[1:3, c(1,3,4),ne]
  DetS2[ne] <- (S2[1,1,ne] * S2[2,2,ne] * S2[3,3,ne]) + (S2[1,2,ne] * S2[2,3,ne] * S2[3,1,ne]) + (S2[1,3,ne] * S2[2,1,ne] * S2[3,2,ne]) -
           (S2[1,3,ne] * S2[2,2,ne] * S2[3,1,ne]) - (S2[1,1,ne] * S2[2,3,ne] * S2[3,2,ne]) - (S2[1,2,ne] * S2[2,1,ne] * S2[3,3,ne])
  #S3#########################################
  pS3[ne] <- X[4,3,ne]
  S3[1:3,1:3, ne] <- X[1:3, c(1,2,4),ne]
  DetS3[ne] <- (S3[1,1,ne] * S3[2,2,ne] * S3[3,3,ne]) + (S3[1,2,ne] * S3[2,3,ne] * S3[3,1,ne]) + (S3[1,3,ne] * S3[2,1,ne] * S3[3,2,ne]) -
           (S3[1,3,ne] * S3[2,2,ne] * S3[3,1,ne]) - (S3[1,1,ne] * S3[2,3,ne] * S3[3,2,ne]) - (S3[1,2,ne] * S3[2,1,ne] * S3[3,3,ne])
  #S4#########################################
  pS4[ne] <- X[4,4,ne]
  S4[1:3,1:3, ne] <- X[1:3, c(1,2,3),ne]
  DetS4[ne] <- (S4[1,1,ne] * S4[2,2,ne] * S4[3,3,ne]) + (S4[1,2,ne] * S4[2,3,ne] * S4[3,1,ne]) + (S4[1,3,ne] * S4[2,1,ne] * S4[3,2,ne]) -
           (S4[1,3,ne] * S4[2,2,ne] * S4[3,1,ne]) - (S4[1,1,ne] * S4[2,3,ne] * S4[3,2,ne]) - (S4[1,2,ne] * S4[2,1,ne] * S4[3,3,ne])

  DetX[ne] <-  (-1 * pS1[ne] * DetS1[ne]) + (pS2[ne] * DetS2[ne]) - (pS3[ne] * DetS3[ne]) + (pS4[ne] * DetS4[ne])
 
  #Adjugate#########################################################################
  A[1,1,ne] <- (X[2,2,ne]*X[3,3,ne]*X[4,4,ne]) + (X[2,3,ne]*X[3,4,ne]*X[4,2,ne]) + (X[2,4,ne]*X[3,2,ne]*X[4,3,ne]) -
            (X[2,4,ne]*X[3,3,ne]*X[4,2,ne]) - (X[2,3,ne]*X[3,2,ne]*X[4,4,ne]) - (X[2,2,ne]*X[3,4,ne]*X[4,3,ne])
  A[1,2,ne] <- (-1*X[1,2,ne]*X[3,3,ne]*X[4,4,ne]) - (X[1,3,ne]*X[3,4,ne]*X[4,2,ne]) - (X[1,4,ne]*X[3,2,ne]*X[4,3,ne]) +
            (X[1,4,ne]*X[3,3,ne]*X[4,2,ne]) + (X[1,3,ne]*X[3,2,ne]*X[4,4,ne]) + (X[1,2,ne]*X[3,4,ne]*X[4,3,ne])
  A[1,3,ne] <- (X[1,2,ne]*X[2,3,ne]*X[4,4,ne]) + (X[1,3,ne]*X[2,4,ne]*X[4,2,ne]) + (X[1,4,ne]*X[2,2,ne]*X[4,3,ne]) -
            (X[1,4,ne]*X[2,3,ne]*X[4,2,ne]) - (X[1,3,ne]*X[2,2,ne]*X[4,4,ne]) - (X[1,2,ne]*X[2,4,ne]*X[4,3,ne])
  A[1,4,ne] <- (-1*X[1,2,ne]*X[2,3,ne]*X[3,4,ne]) - (X[1,3,ne]*X[2,4,ne]*X[3,2,ne]) - (X[1,4,ne]*X[2,2,ne]*X[3,3,ne]) +
            (X[1,4,ne]*X[2,3,ne]*X[3,2,ne]) + (X[1,3,ne]*X[2,2,ne]*X[3,4,ne]) + (X[1,2,ne]*X[2,4,ne]*X[3,3,ne])

  A[2,1,ne] <- (-1*X[2,1,ne]*X[3,3,ne]*X[4,4,ne]) - (X[2,3,ne]*X[3,4,ne]*X[4,1,ne]) - (X[2,4,ne]*X[3,1,ne]*X[4,3,ne]) +
            (X[2,4,ne]*X[3,3,ne]*X[4,1,ne]) + (X[2,3,ne]*X[3,1,ne]*X[4,4,ne]) + (X[2,1,ne]*X[3,4,ne]*X[4,3,ne])
  A[2,2,ne] <- (X[1,1,ne]*X[3,3,ne]*X[4,4,ne]) + (X[1,3,ne]*X[3,4,ne]*X[4,1,ne]) + (X[1,4,ne]*X[3,1,ne]*X[4,3,ne]) -
            (X[1,4,ne]*X[3,3,ne]*X[4,1,ne]) - (X[1,3,ne]*X[3,1,ne]*X[4,4,ne]) - (X[1,1,ne]*X[3,4,ne]*X[4,3,ne])
  A[2,3,ne] <- (-1*X[1,1,ne]*X[2,3,ne]*X[4,4,ne]) - (X[1,3,ne]*X[2,4,ne]*X[4,1,ne]) - (X[1,4,ne]*X[2,1,ne]*X[4,3,ne]) +
            (X[1,4,ne]*X[2,3,ne]*X[4,1,ne]) + (X[1,3,ne]*X[2,1,ne]*X[4,4,ne]) + (X[1,1,ne]*X[2,4,ne]*X[4,3,ne])
  A[2,4,ne] <- (X[1,1,ne]*X[2,3,ne]*X[3,4,ne]) + (X[1,3,ne]*X[2,4,ne]*X[3,1,ne]) + (X[1,4,ne]*X[2,1,ne]*X[3,3,ne]) -
            (X[1,4,ne]*X[2,3,ne]*X[3,1,ne]) - (X[1,3,ne]*X[2,1,ne]*X[3,4,ne]) - (X[1,1,ne]*X[2,4,ne]*X[3,3,ne])

  A[3,1,ne] <- (X[2,1,ne]*X[3,2,ne]*X[4,4,ne]) + (X[2,2,ne]*X[3,4,ne]*X[4,1,ne]) + (X[2,4,ne]*X[3,1,ne]*X[4,2,ne]) -
            (X[2,4,ne]*X[3,2,ne]*X[4,1,ne]) - (X[2,2,ne]*X[3,1,ne]*X[4,4,ne]) - (X[2,1,ne]*X[3,4,ne]*X[4,2,ne])
  A[3,2,ne] <- (-1*X[1,1,ne]*X[3,2,ne]*X[4,4,ne]) - (X[1,2,ne]*X[3,4,ne]*X[4,1,ne]) - (X[1,4,ne]*X[3,1,ne]*X[4,2,ne]) +
            (X[1,4,ne]*X[3,2,ne]*X[4,1,ne]) + (X[1,2,ne]*X[3,1,ne]*X[4,4,ne]) + (X[1,1,ne]*X[3,4,ne]*X[4,2,ne])
  A[3,3,ne] <- (X[1,1,ne]*X[2,2,ne]*X[4,4,ne]) + (X[1,2,ne]*X[2,4,ne]*X[4,1,ne]) + (X[1,4,ne]*X[2,1,ne]*X[4,2,ne]) -
            (X[1,4,ne]*X[2,2,ne]*X[4,1,ne]) - (X[1,2,ne]*X[2,1,ne]*X[4,4,ne]) - (X[1,1,ne]*X[2,4,ne]*X[4,2,ne])
  A[3,4,ne] <- (-1*X[1,1,ne]*X[2,2,ne]*X[3,4,ne]) - (X[1,2,ne]*X[2,4,ne]*X[3,1,ne]) - (X[1,4,ne]*X[2,1,ne]*X[3,2,ne]) +
            (X[1,4,ne]*X[2,2,ne]*X[3,1,ne]) + (X[1,2,ne]*X[2,1,ne]*X[3,4,ne]) + (X[1,1,ne]*X[2,4,ne]*X[3,2,ne])

  A[4,1,ne] <- (-1*X[2,1,ne]*X[3,2,ne]*X[4,3,ne]) - (X[2,2,ne]*X[3,3,ne]*X[4,1,ne]) - (X[2,3,ne]*X[3,1,ne]*X[4,2,ne]) +
            (X[2,3,ne]*X[3,2,ne]*X[4,1,ne]) + (X[2,2,ne]*X[3,1,ne]*X[4,3,ne]) + (X[2,1,ne]*X[3,3,ne]*X[4,2,ne])
  A[4,2,ne] <- (X[1,1,ne]*X[3,2,ne]*X[4,3,ne]) + (X[1,2,ne]*X[3,3,ne]*X[4,1,ne]) + (X[1,3,ne]*X[3,1,ne]*X[4,2,ne]) -
            (X[1,3,ne]*X[3,2,ne]*X[4,1,ne]) - (X[1,2,ne]*X[3,1,ne]*X[4,3,ne]) - (X[1,1,ne]*X[3,3,ne]*X[4,2,ne])
  A[4,3,ne] <- (-1*X[1,1,ne]*X[2,2,ne]*X[4,3,ne]) - (X[1,2,ne]*X[2,3,ne]*X[4,1,ne]) - (X[1,3,ne]*X[2,1,ne]*X[4,2,ne]) +
            (X[1,3,ne]*X[2,2,ne]*X[4,1,ne]) + (X[1,2,ne]*X[2,1,ne]*X[4,3,ne]) + (X[1,1,ne]*X[2,3,ne]*X[4,2,ne])
  A[4,4,ne] <- (X[1,1,ne]*X[2,2,ne]*X[3,3,ne]) + (X[1,2,ne]*X[2,3,ne]*X[3,1,ne]) + (X[1,3,ne]*X[2,1,ne]*X[3,2,ne]) -
            (X[1,3,ne]*X[2,2,ne]*X[3,1,ne]) - (X[1,2,ne]*X[2,1,ne]*X[3,3,ne]) - (X[1,1,ne]*X[2,3,ne]*X[3,2,ne])

  X_inv[1:4, 1:4, ne] <- A[1:4, 1:4, ne]/DetX[ne]
  #################################################################################################################
  #Defining Lambda from Tbar and p
  L1[ne] <- -p * (log(Tbar1[ne]) - loggam(1 + (1/p)))
  L2[ne] <- -p * (log(Tbar2[ne]) - loggam(1 + (1/p)))
  L3[ne] <- -p * (log(Tbar3[ne]) - loggam(1 + (1/p)))
  L4[ne] <- -p * (log(Tbar4[ne]) - loggam(1 + (1/p)))
  
  Tbar1[ne] ~ dlnorm(expData[ne,3], 1/(expData[ne,4]^2))
  Tbar2[ne] ~ dlnorm(expData[ne,7], 1/(expData[ne,8]^2))
  Tbar3[ne] ~ dlnorm(expData[ne,11], 1/(expData[ne,12]^2))
  Tbar4[ne] ~ dlnorm(expData[ne,15], 1/(expData[ne,16]^2))
  
  Beta0[ne] <- X_inv[1,1,ne]*L1[ne] + X_inv[1,2,ne]*L2[ne] + X_inv[1,3,ne]*L3[ne] + X_inv[1,4,ne]*L4[ne]
  Beta1[ne] <- X_inv[2,1,ne]*L1[ne] + X_inv[2,2,ne]*L2[ne] + X_inv[2,3,ne]*L3[ne] + X_inv[2,4,ne]*L4[ne]
  Beta2[ne] <- X_inv[3,1,ne]*L1[ne] + X_inv[3,2,ne]*L2[ne] + X_inv[3,3,ne]*L3[ne] + X_inv[3,4,ne]*L4[ne]
  Beta3[ne] <- X_inv[4,1,ne]*L1[ne] + X_inv[4,2,ne]*L2[ne] + X_inv[4,3,ne]*L3[ne] + X_inv[4,4,ne]*L4[ne]
  
 }

  beta0 <- mean(Beta0[1:N])  
  beta1 <- mean(Beta1[1:N])  
  beta2 <- mean(Beta2[1:N])  
  beta3 <- mean(Beta3[1:N])  


 k ~ dunif(1, 6)
 p ~ dunif(0.01, 5)
}'

data_M3 <- list(I = I, Time = survData$survivalTime, ncond = ncond, Q = Q, U = U, J = J,
                N = N, expData = expData)
params_M3 <- c("beta0", "beta1", "beta2", "beta3", "p", "k", "Beta0", "Beta1", "Beta2", "Beta3")
model3 <- jags.model(textConnection(M3), data = data_M3, n.chains = 3)
update(model3, 5000*10)
mcmc3 <- coda.samples(model3, variable.names = params_M3, n.iter = 5000*700 , thin = 700)
plot(mcmc3)
post_M3 <- as.data.frame(as.matrix(mcmc3))
summary(post_M3)
traceplot(mcmc3)
gelman.plot(mcmc3)
autocorr.plot(mcmc3)
write.csv(post_M3, "post_mod3new.csv", row.names = F)
