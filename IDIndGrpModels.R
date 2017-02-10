# Preliminary Modeling ----
rm(list=ls())
library(rstan)
library(rethinking)

load("IDIndGrpData.Rda")
d.all.s <-d.all[,c("ans.bio", "AgeMean", "ParticipantPopulation")]        #less data for stan (NB! stupid mean age: unknown adult ages currently=20)
d.all.s$ParticipantIDLiberal <-d.all$ParticipantID                        #assume each answer in groups w/ unknown id's (ParticipantPopulation)>6) is different individual
d.all.s$ParticipantIDLiberal[as.numeric(d.all.s$ParticipantPopulation)>6] <- seq(1:length(d.all.s$ParticipantIDLiberal[as.numeric(d.all.s$ParticipantPopulation)>6])) 
d.all.s$ParticipantIDLiberal <- factor(d.all.s$ParticipantIDLiberal)
n.pop<-length(table(d.all.s$ParticipantPopulation))
Pop.Estimates <-data.frame("ParticipantPopulation"=levels(d.all.s$ParticipantPopulation))         #for storing estimates by pop


# 1-L models

hum.1s <- map2stan(                      #no predictors, no MLM
  list(
    ans.bio ~ dbinom(1,p) ,
    logit(p) ~ a ,
    a~dnorm(0,1)
  ), 
  data=d.all.s, start=list(a=0), iter=200, warmup = 100)
precis(hum.1s)

hum.2s <- map2stan(                   #age predictor, no MLM
  list(
    ans.bio ~ dbinom(1,p) ,
    logit(p) ~ a +b*AgeMean ,
    a~dnorm(0,1) ,
    b~dnorm(0,1)
  ), 
  data=d.all.s, start=list(a=0, b=0), iter=200, warmup = 100)
precis(hum.2s)

post <-extract.samples(hum.1s, replace==T)
mcmcpairs(post)

# 2-L models
hum.1s.ml <- map2stan(                      #no predictors, Pop as level
  alist(
    bio ~ dbinom(1,p) ,
    logit(p) ~ a +aj ,
    aj[PartPop] ~ dnorm(0, sigma_pop) ,
    a~dnorm(0,1) ,
    sigma_pop ~dcauchy(0, 1)
  ), 
  data=list(
    bio=d.all.s$ans.bio,
    PartPop=as.integer(d.all.s$ParticipantPopulation)
  ),
  start=list(a=0, aj=rep(0,n.pop), sigma_pop=1) ,
  iter=200, warmup = 100)
precis(hum.1s.ml)
post<-extract.samples(hum.1s.ml)
est.h.1s.ml <- sapply(1:n.pop, function(j)
  logistic( median( post$a + post$aj[,j] ) ) ) 
plot(est.h.1s.ml)

hum.2s.ml <- map2stan(                      #age predictor, Pop as level
  alist(
    bio ~ dbinom(1,p) ,
    logit(p) ~ a +aj +b*age ,
    aj[PartPop] ~ dnorm(0, sigma_pop) ,
    c(a, b)~dnorm(0,1) ,
    sigma_pop ~dcauchy(0, 1)
  ), 
  data=list(
    bio=d.all.s$ans.bio,
    PartPop=as.integer(d.all.s$ParticipantPopulation) ,
    age=d.all.s$AgeMean
  ),
  start=list(a=0, aj=rep(0,n.pop), b=0, sigma_pop=1) ,
  iter=200, warmup = 100)

precis(hum.2s.ml)
post <-extract.samples(hum.2s.ml)
Pop.Estimates$bio.est <-sapply(1:n.pop, function(j)
  logistic(median(post$a + post$aj[,j] + post$b*20)) )    #NB: shitty estimate given age 
par(las=2, oma=c(10,0,0,0))
plot(Pop.Estimates$ParticipantPopulation, Pop.Estimates$bio.est, ylim=c(0,1), ylab="estimated pr. birth parent ID")

hum.3s.ml <- map2stan(                      #age predictor- varying by Pop, Pop as level
  alist(
    bio ~ dbinom(1,p) ,
    logit(p) ~ a +aj +(b+bj)*age ,
    c(a, b)~dnorm(0,1) ,
    c(aj, bj)[PartPop] ~ dmvnorm2(0, sigma_pop, rho_pop) ,
    sigma_pop ~dcauchy(0, 1) ,
    rho_pop ~dlkjcorr(2)
  ), 
  data=list(
    bio=d.all.s$ans.bio,
    age=d.all.s$AgeMean ,
    PartPop=as.integer(d.all.s$ParticipantPopulation) 
  ),
  start=list(a=0, aj=rep(0,n.pop), b=0, bj=rep(0,n.pop), sigma_pop=c(1,1), rho_pop=diag(2)) ,
  iter=200, warmup = 100)

precis(hum.3s.ml, depth = 2)
compare(hum.1s.ml, hum.2s.ml, hum.3s.ml) #no age, constant age, slope-varying age (always with pop RE)

# 3-L models
hum.1s.ml3 <- map2stan(                      #no predictors, Pop & (liberally defined) Individuals as levels, w/ stan
  alist(
    bio ~ dbinom(1,p) ,
    logit(p) ~ a +aj +ak,
    aj[PartPop] ~ dnorm(0, sigma_pop) ,
    ak[Partic] ~ dnorm(0, sigma_partic) ,
    a~dnorm(0,1) ,
    sigma_pop ~dcauchy(0, 1) ,
    sigma_partic ~dcauchy(0, 1)
  ), 
  data=list(
    bio=d.all.s$ans.bio ,
    PartPop=as.integer(d.all.s$ParticipantPopulation) ,
    Partic=as.integer(d.all.s$ParticipantIDLiberal)
  ),
  start=list(a=0, aj=rep(0,length(table(d.all.s$ParticipantPopulation))), 
             ak=rep(0,length(table(d.all.s$ParticipantIDLiberal)) ), sigma_pop=1, sigma_partic=1) ,
  iter=200, warmup = 100)

precis(hum.1s.ml3)
