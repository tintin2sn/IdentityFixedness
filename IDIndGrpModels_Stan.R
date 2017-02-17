# Preliminary Modeling ----
rm(list=ls())
library(rstan)
library(rethinking)

load("IDIndGrpData20170217.Rda")

# cbind( d.all$ParticipantPopulation , d.all$Age , d.all$Adult , d.all$AgeMin , d.all$AgeMax , d.all$AgeSD , d.all$AgeMean )

# looks like 5 types in data:
# (1) individual age observed: 1, 2, 3, 4, 5, 6
# (2) min/max/mean known: 7, 9, 11, 13, 17, 18, 20, 21, 22, 25
# (3) mean/sd known: 8, 16
# (4) only adult known: 9, 11, 12, 23
# (5) only min/max known: 14, 19, 20, 21, 24
# pops 9, 11 are strange, as 1 year ranges and mean at min of range

# prep age impute meta data
N_missing_ages <- sum(is.na(d.all$Age))
miss_idx <- which(is.na(d.all$Age))
age_miss_min <- d.all$AgeMin[miss_idx]
age_miss_max <- d.all$AgeMax[miss_idx]
for ( i in 1:N_missing_ages ) {
    if ( is.na(age_miss_min[i]) ) {
        if ( !is.na(d.all$Adult[miss_idx[i]]) ) {
            if ( d.all$Adult[miss_idx[i]]==1 ) {
                # adult
                age_miss_min[i] <- 18
                age_miss_max[i] <- 100
            } else {
                # not adult
                age_miss_min[i] <- 0
                age_miss_max[i] <- 18
            }
        } else {
            # adult status missing, so set maximal limits
            age_miss_min[i] <- 0
            age_miss_max[i] <- 100
        }
    }
}
miss_idx2 <- rep(0,nrow(d.all))
miss_idx2[miss_idx] <- 1:N_missing_ages

NA2neg <- function(x) ifelse(is.na(x),-1,x)
d.all$pop_id <- as.integer(d.all$ParticipantPopulation)
age_meta <- as.matrix(cbind( d.all$Adult , d.all$AgeMin , d.all$AgeMax , d.all$AgeMean , d.all$AgeSD ))
age_meta <- NA2neg(age_meta)

dat_list <- list(
    y = d.all$ans.bio,
    N = nrow(d.all),
    N_pops = max(d.all$pop_id),
    pop_id = d.all$pop_id,
    age = NA2neg(d.all$Age),
    age_meta = age_meta,
    N_missing_ages = N_missing_ages,
    age_miss_min = age_miss_min,
    age_miss_max = age_miss_max,
    miss_idx = miss_idx2
)

m <- stan( file="model_impute.stan" , data=dat_list , 
    chains=1 , control=list(adapt_delta=0.99) )

stanergy(m,binwidth=NULL)

post <- extract.samples(m)
plot(precis(m,2,pars="age_impute"))

