#Erin Bohaboy simulated population code
#modified by Eva Schemmel 6_1_2020
#reviewed by Jon Brodziak 6_21_2020
#adding POS sampling to simulation
#Erin added selectivity parameters to the pop simulation

#load packages
library(reshape)
library(dplyr)
library(ggplot2)
library(magrittr)


#Onaga sampling example

#  1. SIMULATE A POPULATION --------

simulate_population_harvest <- function(Linf,Linf_sd, M, F, mincat, catsd, maxcat, maxcatsd, L0, L0_sd, k, Amax, N){
  # create Amax + 1 cohorts of N Age zero fish, put them all into a list object
  all_cohorts <- list()
  all_harvest <- list()
  hr = 1-exp(-F)
  mr = 1-exp(-M)
  for (i in 1:(Amax+1)) {
    cohort <- matrix(nrow=N,ncol=Amax+3)
    cohort[,3] = rnorm(N, mean = L0, sd = L0_sd)
    # assign an Linf
    cohort[,2] = rnorm(N, mean = Linf, sd = Linf_sd)
    # we gave length at A0, calculate t0 for each fish to use in von Bert
    cohort[,1]= (log(1-(cohort[,3]/cohort[,2])))/k
    
    # ----- step A. grow the fish, all survive
    for (A in 1:Amax) {
      cohort[,A+3] = cohort[,2]*(1-exp(-k*(A-cohort[,1])))
    }
    # ----  step B. apply natural mortality 
    # natural mortality
    mr_cohort = matrix(mr, nrow=N,ncol=Amax+1)    
    surv_mr = matrix(rbinom(N*(Amax+1), size = 1, prob = 1-(mr_cohort)),nrow=N, ncol=Amax+1)
    # all fish survive to age 0
    surv_mr[,1] = 1
    # modify mr into a cumulative cross-product
    for (A in 1:(Amax)) {
      surv_mr[,A+1]=surv_mr[,A]*surv_mr[,A+1]
    }
    # lengths for fish that survive natural mortality, 0s for those that don't
    mr_surv_lengths = cohort[,3:(Amax+3)]*surv_mr
    
    # ----- step C. apply fishing mortality and collect catch
    # will depend on fishery selectivity at length
    selex_cohort = pmin((1-pnorm(mincat,cohort[,3:(Amax+3)],catsd)),pnorm(maxcat ,cohort[,3:(Amax+3)],maxcatsd))
    hr_cohort = hr*selex_cohort
    catch_cohort = matrix(rbinom(N*(Amax+1), size = 1, prob = hr_cohort),nrow=N, ncol=Amax+1)
    
    #  each fish only gets caught once
    # modify catch_cohort into additive vector over ages so once a fish gets "caught" the second time, catch_cohort > 1
    for (A in 1:(Amax)) {
      catch_cohort[,A+1]=catch_cohort[,A]+catch_cohort[,A+1]
    }
    # this is goofy- take another additive product so that the age a fish gets caught first = 1, subsequent ages > 1
    for (A in 1:(Amax)) {
      catch_cohort[,A+1]=catch_cohort[,A]+catch_cohort[,A+1]
    }
    # the age a fish gets caught first = 1, all others = 0
    catch_once_cohort <- replace(catch_cohort, catch_cohort > 1, 0)
    
    # lengths for fish that survive natural mortality AND GET CAUGHT, 0s for those that don't
    catch_lengths = mr_surv_lengths*catch_once_cohort
    #  save it in the all_harvest list
    all_harvest[[i]] <- catch_lengths
    
    # apply fishing mortality by taking the fish that got caught (identified above) out of the population of fish that survived natural mortality
    hr_surv <- 1-catch_once_cohort
    # turn this into a cumulative cross-product
    for (A in 1:(Amax)) {
      hr_surv[,A+1]=hr_surv[,A]*hr_surv[,A+1]
    }
    # lengths for fish that survive, 0s for those that don't
    cohort_lengths = mr_surv_lengths*hr_surv
    #  save it
    all_cohorts[[i]] <- cohort_lengths
  }
  
  # now pull one age from each cohort to build the overall population
  # create a population object using Age0
  population <- data.frame(age=0, length = all_cohorts[[1]][,1])
  #drop dead fish
  population <- subset(population , length != 0)
  
  for (i in 1:(Amax)) {
    each_age <- data.frame(age=i, length = all_cohorts[[(i+1)]][,(i+1)])
    #drop dead fish
    each_age <- subset(each_age , length != 0)
    #save
    population <- rbind(population, each_age)
  }
  #  we have a population
  
  #  repeat for harvest 
  harvest <- data.frame(age=0, length = all_harvest[[1]][,1])
  #drop dead fish
  harvest <- subset(harvest , length != 0)
  for (i in 1:(Amax)) {
    each_h_age <- data.frame(age=i, length = all_harvest[[(i+1)]][,(i+1)])
    #drop dead fish
    each_h_age <- subset(each_h_age , length != 0)
    #save
    harvest <- rbind(harvest, each_h_age)
  }
  #  we have a harvest
  
  
  # add a warning message if L0_sd was large compared to L0 such that we have negative length fish
  if (min(population$length)<0) { print("Warning: negative length fish, L0_sd too large compared to L0") }
  
  return(list(population = population, harvest = harvest))
  
}# end function


pop_harv <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.17, F=0.17, mincat=0, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.4, Amax=40, N=10000)
population<-pop_harv$harvest[]
population_true<-pop_harv$population[]



#comparing mincat 
pop_harv10 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.0, F=0.17, mincat=10, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.13, Amax=40, N=10000)
population10<-pop_harv10$harvest[]
population10$mincat<-10

pop_harv20 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.0, F=0.17, mincat=20, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.13, Amax=40, N=10000)
population20<-pop_harv20$harvest[]
population20$mincat<-20

pop_harv30 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.0, F=0.17, mincat=30, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.13, Amax=40, N=10000)
population30<-pop_harv30$harvest[]
population30$mincat<-30

pop_harv40 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.0, F=0.17, mincat=40, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.13, Amax=40, N=10000)
population40<-pop_harv40$harvest[]
population40$mincat<-40

pop_harv60 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.0, F=0.17, mincat=60, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.13, Amax=40, N=10000)
population60<-pop_harv60$harvest[]
population60$mincat<-60

population<-rbind(population10, population20, population30,population40, population60)

#comparing k with moderate to high mortality
pop_harv1 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.17, F=0.17, mincat=00, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.1, Amax=40, N=10000)
population1<-pop_harv1$harvest[]
population1$k<-.1

pop_harv2 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.17, F=0.17, mincat=00, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.2, Amax=40, N=10000)
population2<-pop_harv2$harvest[]
population2$k<-.20

pop_harv3 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.17, F=0.17, mincat=00, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.3, Amax=40, N=10000)
population3<-pop_harv3$harvest[]
population3$k<-.30

pop_harv4 <-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.17, F=0.17, mincat=00, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.4, Amax=40, N=10000)
population4<-pop_harv4$harvest[]
population4$k<-.40

pop_harv5<-simulate_population_harvest(Linf=100,Linf_sd=5, M=0.17, F=0.17, mincat=00, catsd=0,maxcat = 200, maxcatsd = 0, L0=10, L0_sd=2.5, k=0.5, Amax=40, N=10000)
population5<-pop_harv5$harvest[]
population5$k<-.50

population<-rbind(population1, population2, population3,population4, population5)



#  ------------------------------------------------------------------------
#	1B. Add function plot_Fselex to visualize fishery selectivity

# for dome shaped, chose inflection (maxcat) and slope (maxcatsd) for righthand side
# for logistic, set maxcat >> Linf and maxcatsd = 0

plot_Fselex <- function(Linf, mincat, catsd, maxcat, maxcatsd) {
  xmax = round((Linf + 0.20*Linf),0)
  plot(seq(0,xmax ,1),(pmin((1-pnorm(mincat,seq(0,xmax ,1),catsd)),pnorm(maxcat ,seq(0,xmax ,1),maxcatsd ))),
       type="l",xlab="Fish Length",ylab="Fishery Selectivity")
}


# examples
# plot_Fselex(Linf=100, mincat=40, catsd=10, maxcat=200, maxcatsd=0) 
# plot_Fselex(Linf=101, mincat=0, catsd=0, maxcat=90, maxcatsd=5) 

#plot_Fselex(Linf=100, mincat=25, catsd=5, maxcat=200, maxcatsd=5) 

#   -----------------------------------------------------------------------
#	2. SAMPLE THE POPULATION

#  draw a sample from the population
#	in this instance, we sample entire population

#  this is where selectivity of the survey would be accounted for.


# -------
# a truly random sample:
random_sample <- function(population, sample_size) {
  return(population[sample(nrow(population),sample_size),])
}
sample <- random_sample(population, 300)

# -------
# an FOS sample- 

# read in the function
FOS_sample <- function(population, FOS_plan) {
  rm(grow_samp)
  grow_samp = data.frame(age=numeric(), length=numeric())
  for (i in 1:nrow(FOS_plan)) {
    pop_bin <- subset(population, length > FOS_plan$bin_lower[i] & length <= FOS_plan$bin_upper[i])
    
    if(nrow(pop_bin) >=  FOS_plan$n_samps[i]) {
      bin_samp <- pop_bin[sample(nrow(pop_bin),FOS_plan$n_samps[i]),]
      grow_samp <- rbind(grow_samp,bin_samp)
    } else {
      grow_samp <- rbind(grow_samp,pop_bin)
    }
  }
  return(grow_samp)
}

# set up a "sampling plan", can streamline this later

# chose number of bins and sample size so samples/bin is an integer-
# i.e. 11 bins, 300 samples = 27 samps in most bins, 28 in 3 of them...

FOS_plan <- data.frame(bin_lower = seq(0,100,10),
                       bin_upper = seq(10,110,10),
                       n_samps = c(27, 27, 28, 28, 27, 27, 27, 27, 28, 27, 27))

#remove juveniles <40
FOS_plan <- data.frame(bin_lower = seq(0,100,10),
                       bin_upper = seq(10,110,10),
                       n_samps = c(0, 0, 0, 0, 42, 42, 42, 42, 42, 42, 42))


#remove juveniles <20
FOS_plan <- data.frame(bin_lower = seq(0,100,10),
                       bin_upper = seq(10,110,10),
                       n_samps = c(0, 0, 33, 33, 33, 33, 33, 33, 33, 33, 33))

sample <- FOS_sample(population, FOS_plan)

#generate weights for POS based on length frequency
brks <- seq(0,120,10)
population<-population %>% subset(length>=0)
population$binl <- cut(population[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
library(magrittr)
library(dplyr)
population<-population %>% 
  group_by(binl) %>%
  mutate(prop=n()) %>%
  group_by() %>%
  mutate(n=sum(n())) %>%
  mutate(prop1=prop/n)
head(population)

POS<-sample_n(population, 300 ,replace=FALSE, weight=prop1) #13 bins (bins * n per bin = sample size)
POS_n<-POS %>%
  group_by(binl) %>%
  summarize(n=length(binl))

hist(population$length)
par(mfrow=c(1,2))
hist(population$age, xlab= "Age", ylab="Pop N at Age", main="", cex=1.5, cex.lab=1.5, ,cex.axis=1.2)
plot(population$age, population$length, col="black", ylab="Length (cm)", xlab="Age", ylim=c(0, 120), xlim=c(0,40), cex=1.5, cex.lab=1.5, ,cex.axis=1.2)



#sample <-  POS


#   -----------------------------------------------------------------------
#	3. FIT A VON BERTALANFFY GROWTH CURVE FOR INDIVIDUAL SAMPLING RUNS AND RESULTS (SEE BOOTSTRAP WRAPPER BELOW)

#	this is just sum of squares, does not estimate L0_sd or Linf_sd
#	add that later
#	I put a normal error structure in, so should still be true? Probably needs to be modified.


#  fit a von bert to the population data
#	when making the von Bert func, be cautious of reusing variable names

# initialize parameter vector
theta <- c(100,0,0.2)
Amax=40
ssr = function(theta){	
  LINF = theta[1]
  A0 = theta[2]
  K = theta[3]
  L_pred = LINF*(1-exp(-K*(sample$age-A0)))
  ssr = sum((sample$length-L_pred)^2)
  return(ssr)
}

fit = optim(theta,ssr)
fit$par


# make a plot
#looking at mincat on random sampling
population10<-population %>% subset(mincat==10)

random_sample <- function(population10, sample_size) {
  return(population[sample(nrow(population10),sample_size),])
}
sample <- random_sample(population10, 300)

#30
population30<-population %>% subset(mincat==30)

random_sample <- function(population30, sample_size) {
  return(population[sample(nrow(population30),sample_size),])
}
sample <- random_sample(population30, 300)

#60
population60<-population %>% subset(mincat==60)

random_sample <- function(population60, sample_size) {
  return(population[sample(nrow(population60),sample_size),])
}
sample <- random_sample(population60, 300)

plot(sample$age, sample$length, ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", col="white" ,cex.axis=1.2, cex.lab=1.5)

fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths, lwd=2,col="blue") #change color of line for each sample
lines(plot_ages,pred_lengths, lwd=2,col="purple")
lines(plot_ages,pred_lengths, lwd=2,col="grey")

#sample=FOS
#remove juveniles <60
FOS_plan <- data.frame(bin_lower = seq(0,100,10),
                       bin_upper = seq(10,110,10),
                       n_samps = c(0, 0, 0, 0, 0, 0, 60, 60, 60, 60, 60))
sample <- FOS_sample(population60, FOS_plan)

#remove juveniles <30
FOS_plan <- data.frame(bin_lower = seq(0,100,10),
                       bin_upper = seq(10,110,10),
                       n_samps = c(0, 0, 0, 37, 37 , 37, 37, 37, 37, 37, 37))
sample <- FOS_sample(population30, FOS_plan)
#remove juveniles <10
FOS_plan <- data.frame(bin_lower = seq(0,100,10),
                       bin_upper = seq(10,110,10),
                       n_samps = c(0, 30,30,30,30,30,30,30,30,30,30))
sample <- FOS_sample(population10, FOS_plan)

plot(sample$age, sample$length, col="white", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", cex.axis=1.2, cex.lab=1.5)
fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths, lwd=2,col="grey") #change color of line for each sample
lines(plot_ages,pred_lengths, lwd=2,col="purple") 
lines(plot_ages,pred_lengths, lwd=2,col="blue") 

#POS juveniles #set population to population10, population30, population60 and rerun
population60$binl <- cut(population60[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
population60<-population60 %>% 
  group_by(binl) %>%
  mutate(prop=n()) %>%
  group_by() %>%
  mutate(n=sum(n())) %>%
  mutate(prop1=prop/n)

POS<-sample_n(population60, 300 ,replace=FALSE, weight=prop1) #13 bins (bins * n per bin = sample size)
sample <- POS

population30$binl <- cut(population30[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
population30<-population30 %>% 
  group_by(binl) %>%
  mutate(prop=n()) %>%
  group_by() %>%
  mutate(n=sum(n())) %>%
  mutate(prop1=prop/n)

POS<-sample_n(population30, 300 ,replace=FALSE, weight=prop1) #13 bins (bins * n per bin = sample size)
sample <- POS

population10$binl <- cut(population10[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
population10<-population10 %>% 
  group_by(binl) %>%
  mutate(prop=n()) %>%
  group_by() %>%
  mutate(n=sum(n())) %>%
  mutate(prop1=prop/n)

POS<-sample_n(population10, 300 ,replace=FALSE, weight=prop1) #13 bins (bins * n per bin = sample size)
sample <- POS


plot(sample$age, sample$length, col="white", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", cex.axis=1.2, cex.lab=1.5)

fit = optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths,lwd=2, col="grey") #change color of line for each sample
lines(plot_ages,pred_lengths,lwd=2, col="purple") 
lines(plot_ages,pred_lengths,lwd=2, col="blue") 

#####look at k on sampling with moderate to high mortality 
#looking at mincat on random sampling
par=mfrow=c(1,3)
population1<-population %>% subset(k==.10)

random_sample <- function(population1, sample_size) {
  return(population1[sample(nrow(population1),sample_size),])
}
sample <- random_sample(population1, 300)

#k=.3
population3<-population %>% subset(k==.30)

random_sample <- function(population3, sample_size) {
  return(population3[sample(nrow(population3),sample_size),])
}
sample <- random_sample(population3, 300)

#k=.50
population5<-population %>% subset(k==.50)

random_sample <- function(population5, sample_size) {
  return(population5[sample(nrow(population5),sample_size),])
}
sample <- random_sample(population5, 300)

plot(sample$age, sample$length, ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", col="white" ,cex.axis=1.2, cex.lab=1.5)

fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths,lwd=2, col="blue") #change color of line for each sample
lines(plot_ages,pred_lengths, lwd=2,col="purple")
lines(plot_ages,pred_lengths, lwd=2,col="grey")

#sample=FOS
sample <- FOS_sample(population1, FOS_plan)
sample <- FOS_sample(population3, FOS_plan)
sample <- FOS_sample(population5, FOS_plan)
plot(sample$age, sample$length, col="white", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", cex.axis=1.2, cex.lab=1.5)
fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))
lines(plot_ages,pred_lengths, lwd=2,col="blue")
lines(plot_ages,pred_lengths,lwd=2, col="purple")
lines(plot_ages,pred_lengths,lwd=2, col="grey") #change color of line for each sample


#POS juveniles #set population to population10, population30, population60 and rerun
population1$binl <- cut(population1[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
population1<-population1 %>% 
  group_by(binl) %>%
  mutate(prop=n()) %>%
  group_by() %>%
  mutate(n=sum(n())) %>%
  mutate(prop1=prop/n)


POS<-sample_n(population1, 300 ,replace=FALSE, weight=prop1) #13 bins (bins * n per bin = sample size)

sample <- POS
plot(sample$age, sample$length, col="white", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", cex.axis=1.2, cex.lab=1.5)
fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths,lwd=2, col="blue")
lines(plot_ages,pred_lengths,lwd=2, col="purple")
lines(plot_ages,pred_lengths,lwd=2, col="grey") #change color of line for each sample


### plot different sampling designs - random, FOS, POS #####
#Random
par(mfrow=c(1,3))
sample <- random_sample(population, 300) 


plot(sample$age, sample$length, col="black", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age",  cex.axis=1.2, cex.lab=1.5)

fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths, col="blue")
text(30, 80, substitute(L[inf]==LINF, list(LINF=LINF)))
text(30, 70, substitute(A[0]==A0, list(A0=A0)))
text(30, 60, substitute("K"==K, list(K=K)))



#estimate Linf from von B
#Random
modrand<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
              data=sample, 
              start = list (Linf = 100,
                            K1 = 0.2,
                            t0=0.2))
summary(modrand)
confint(modrand)
AIC(modrand)


#FOS
sample <- FOS_sample(population, FOS_plan)
#sample=FOS

plot(sample$age, sample$length, col="black", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", cex.axis=1.2, cex.lab=1.5)
fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths, col="blue")
text(30, 80, substitute(L[inf]==LINF, list(LINF=LINF)))
text(30, 70, substitute(A[0]==A0, list(A0=A0)))
text(30, 60, substitute("K"==K, list(K=K)))

#estimate Linf from von B
modFOS<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)
                                     ))),
             data=sample, 
             start = list (Linf = 100,
                           K1 = 0.2,
                           t0=0.2))
summary(modFOS)
confint(modFOS)
AIC(modFOS)



#POS
sample <- POS
plot(sample$age, sample$length, col="black", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age", cex.axis=1.2, cex.lab=1.5)
fit = optim(theta,ssr)
fit$par

fit= optim(theta,ssr)
fit$par

# add fit to the sample plot
LINF = round(fit$par[1],2)
A0 = round(fit$par[2],2)
K = round(fit$par[3],2)
plot_ages = seq(0,Amax, 0.1)
pred_lengths = LINF*(1-exp(-K*(plot_ages-A0)))

lines(plot_ages,pred_lengths, col="blue")
text(30, 80, substitute(L[inf]==LINF, list(LINF=LINF)))
text(30, 70, substitute(A[0]==A0, list(A0=A0)))
text(30, 60, substitute("K"==K, list(K=K)))

#estimate Linf from von B
modFOS<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
             data=sample, 
             start = list (Linf = 100,
                           K1 = 0.2,
                           t0=0.2))
summary(modFOS)
confint(modFOS)
AIC(modFOS)





###########bootstrap wrapper ############
#packages
library(reshape2)


#set up data for resampling
#generate weights for POS based on length frequency
brks <- seq(0,120,10)
population$binl <- cut(population[,"length"], breaks=brks, labels=brks[1:(length(brks)-1)],right=F)
population<-population %>% 
  group_by(binl) %>%
  mutate(prop=n()) %>%
  group_by() %>%
  mutate(n=sum(n())) %>%
  mutate(prop1=prop/n)
head(population)

POS_n<-population %>%
  group_by(binl) %>%
  summarize(n=length(binl))


#POS
df_total = data.frame()
age<-seq(0,40, by =1)
mod<-list()
POS<-sample_n(population, 300, replace=FALSE, weight=prop1)
output <- data.frame(age=numeric(),L_lower95=numeric(),L_upper95=numeric())
for(x in 1:100){
  POS[[x]]<-sample_n(population, 300, replace=FALSE, weight=prop1)
  mod[[x]]<- nls(length ~ Linf* (1-exp(-(K*(age-t0)))),
                 data=POS[[x]], 
                 start = list (Linf = 100,
                               K = 0.13,
                               t0=0.1))
  print(summary(mod[[x]]))
  pred.length_POS <- coef(mod[[x]])[1]*(1 - exp(-coef(mod[[x]])[2]*(age -0)))
  age= age
  length= pred.length_POS
  df=data.frame(age, length)
  df$run<-x
  df_total <- rbind(df_total,df)
  
}

POSdata<-df_total

bootstrap_output <- data.frame(age=numeric(),L_lower95=numeric(),L_upper95=numeric(), meanl=numeric())
Amax=40
for (i in 0:Amax) {
  sub <- subset(POSdata, age==i)
  lower95 <- as.numeric(quantile(sub$length,0.025,FALSE))
  upper95 <- as.numeric(quantile(sub$length,0.975,FALSE))
  mean<-as.numeric(mean(sub$length))
  bootstrap_output[i+1,1] <- i
  bootstrap_output[i+1,2] <- lower95
  bootstrap_output[i+1,3] <- upper95
  bootstrap_output[i+1,4] <- mean
}

POSCI<-bootstrap_output

#for non-contrained models

POS_coef<-(lapply(mod,coef))
test2<-melt(data=POS_coef, "Linf", "K1", "t0")

t0_test <- test2[c(rep(FALSE,2),TRUE), ]
Linf_test <- test2[c(TRUE,rep(FALSE,2)), ]
test0<-test2[-1,]
K_test<-test0[c(TRUE,rep(FALSE,2)), ]

mean(Linf_test$value)
quantile(Linf_test$value, prob=0.025)
quantile(Linf_test$value, prob=0.975) 

mean(K_test$value)
quantile(K_test$value, prob=0.025)
quantile(K_test$value, prob=0.975) 



#FOS 
df_total = data.frame()
#add sample size per bin for FOS
n <- 42 #max length/10cm size bins  for juveniles removed <40cm
n<-33 ##max length/10cm size bins  for juveniles removed <20cm
population$n<-n #needed to add a row to the table to use sample_n with bins that have less than prefered sample size

FOS <- population %>%dplyr::group_by(binl)%>% sample_n(pmin(n(), n))
FOS_n<-FOS %>%
  group_by(binl) %>%
  summarize(n=length(binl))
sum(FOS_n$n)
age<-seq(0,40, by =1)
mod<-list()
for(x in 1:100){
    FOS[[x]] <- population %>%dplyr::group_by(binl)%>% sample_n(pmin(n(), n)) #takes the n specified above or if there are not enough in the bin than the total number of samples
    mod[[x]]<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
                   data=FOS[[x]], 
                   start = list (Linf = 100,
                                 K1 = 0.13,
                                 t0=0.1))
    print(summary(mod[[x]]))
    pred.length_POS <- coef(mod[[x]])[1]*(1 - exp(-coef(mod[[x]])[2]*(age -0)))
    age= age
    length= pred.length_POS
    df=data.frame(age, length)
    df$run<-x
    df_total <- rbind(df_total,df)
    
  }

LOSdata <- df_total
bootstrap_output <- data.frame(age=numeric(),L_lower95=numeric(),L_upper95=numeric(), meanl=numeric())
Amax=40
for (i in 0:Amax) {
  sub <- subset(LOSdata, age==i)
  lower95 <- as.numeric(quantile(sub$length,0.025,FALSE))
  upper95 <- as.numeric(quantile(sub$length,0.975,FALSE))
  mean<-as.numeric(mean(sub$length))
  bootstrap_output[i+1,1] <- i
  bootstrap_output[i+1,2] <- lower95
  bootstrap_output[i+1,3] <- upper95
  bootstrap_output[i+1,4] <- mean
}


LOSCI<-bootstrap_output

LOS_coef<-(lapply(mod,coef))

#for non-contrained models
test2<-melt(data=LOS_coef, "Linf", "K1", "t0")

t0_test <- test2[c(rep(FALSE,2),TRUE), ]
Linf_test <- test2[c(TRUE,rep(FALSE,2)), ]
test0<-test2[-1,]
K_test<-test0[c(TRUE,rep(FALSE,2)), ]

mean(Linf_test$value)
quantile(Linf_test$value, prob=0.025)
quantile(Linf_test$value, prob=0.975) 

mean(K_test$value)
quantile(K_test$value, prob=0.025)
quantile(K_test$value, prob=0.975) 

##random
ran <- random_sample(population, 300)
age<-seq(0,40, by =1)
mod<-list()
for(x in 1:100){
  ran[[x]] <- random_sample(population, 300)
  mod[[x]]<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
                 data=ran[[x]], 
                 start = list (Linf = 100,
                               K1 = 0.13,
                               t0=0.1))
  print(summary(mod[[x]]))
  pred.length_POS <- coef(mod[[x]])[1]*(1 - exp(-coef(mod[[x]])[2]*(age -0)))
  age= age
  length= pred.length_POS
  df=data.frame(age, length)
  df$run<-x
  df_total <- rbind(df_total,df)
  
}

Randata <- df_total

bootstrap_output <- data.frame(age=numeric(),L_lower95=numeric(),L_upper95=numeric(), meanl=numeric())
Amax=40
for (i in 0:Amax) {
  sub <- subset(Randata, age==i)
  lower95 <- as.numeric(quantile(sub$length,0.025,FALSE))
  upper95 <- as.numeric(quantile(sub$length,0.975,FALSE))
  mean<-as.numeric(mean(sub$length))
  bootstrap_output[i+1,1] <- i
  bootstrap_output[i+1,2] <- lower95
  bootstrap_output[i+1,3] <- upper95
  bootstrap_output[i+1,4] <- mean
}


RanCI<-bootstrap_output

Ran_coef<-(lapply(mod,coef))

#for non-contrained models
test2<-melt(data=LOS_coef, "Linf", "K1", "t0")

t0_test <- test2[c(rep(FALSE,2),TRUE), ]
Linf_test <- test2[c(TRUE,rep(FALSE,2)), ]
test0<-test2[-1,]
K_test<-test0[c(TRUE,rep(FALSE,2)), ]

mean(Linf_test$value)
quantile(Linf_test$value, prob=0.025)
quantile(Linf_test$value, prob=0.975) 

mean(K_test$value)
quantile(K_test$value, prob=0.025)
quantile(K_test$value, prob=0.975) 


###sim. pop bootstrapped von b
# df_total = data.frame()
# age<-seq(0,40, by =1)
# mod<-list()
# 
# for (i in 1:100) {
#   mod[[i]]<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
#                  data=population_true, 
#                  start = list (Linf = 100,
#                                K1 = 0.13,
#                                t0=0.1))
#   print(summary(mod[[i]]))
#   pred.length_pop <- coef(mod[[i]])[1]*(1 - exp(-coef(mod[[i]])[2]*(age -0)))
#   age= age
#   length= pred.length_pop
#   df=data.frame(age, length)
#   df$run<-i
#   df_total <- rbind(df_total,df)
# }


#CI for simulated population #true_population

output <- data.frame(age=numeric(),L_lower95=numeric(),L_upper95=numeric(),meanl=numeric())

for (i in 0:Amax) {
  sub <- subset(population_true, age==i)
  lower95 <- as.numeric(quantile(sub$length,0.025,FALSE))
  upper95 <- as.numeric(quantile(sub$length,0.975,FALSE))
  mean<-as.numeric(mean(sub$length))
  output[i+1,1] <- i
  output[i+1,2] <- lower95
  output[i+1,3] <- upper95
  output[i+1,4] <- mean
}

# plot(population_true$age, population_true$length, col="black", ylim=c(0, 120), xlim=c(0,40), ylab="Length (cm)", xlab="Age")
plot(output$age,output$L_lower95,type="l",lty=2,ylim=c(0,120))
lines(output$age,output$L_upper95, lty=2)


POPdata <- df_total


POP_coef<-(lapply(mod,coef))

#######looking at CI for entire population#########


mod <- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
          data=population_true, 
          start = list (Linf = 100,
                        K1 = 0.13,
                        t0=0.1))

## Form data plot and smooth line for the predictions
plot(length ~ age, data = population_true, bty="l", col = "white",
     main = "")
#lines(age, predict(mod, list(age = age)))
lines(output$age, output$meanl, lty=1, lwd=2)
lines(output$age,output$L_lower95, lty=2, lwd=2)
lines(output$age,output$L_upper95, lty=2, lwd=2)
#POS
lines(POSCI$age, POSCI$meanl, col="tomato3",lty=1,lwd=2)
lines(POSCI$age,POSCI$L_lower95,col="tomato3", lty=2, lwd=2)
lines(POSCI$age,POSCI$L_upper95,col="tomato3", lty=2, lwd=2)
lines(LOSCI$age, LOSCI$meanl, col="steelblue1",lty=1, lwd=2)
lines(LOSCI$age,LOSCI$L_lower95,col="steelblue1", lty=2, lwd=2)
lines(LOSCI$age,LOSCI$L_upper95,col="steelblue1", lty=2, lwd=2)
#lines(RanCI$age, RanCI$meanl, col="darkgreen",lty=1, lwd=2)
#lines(RanCI$age,RanCI$L_lower95,col="darkgreen", lty=2)
#lines(RanCI$age,RanCI$L_upper95,col="darkgreen", lty=2)




library(nlstools) 
bootpop <- nlsBoot(mod)
bootpop$coefboot[]
bootpop$bootCI[]
bootpop$estiboot[]



###plot it by run
library(ggplot2)
##plot mean and CI for 1000 runs
sim_boot<-ggplot() + 
  geom_line(data=POSmean, aes(x=age, y=length),size=1, color='tomato3')+
  geom_line(data=LOSmean, aes(x=age, y=length), size=1,color='slateblue1') + 
  geom_line(data=Ranmean, aes(x=age, y=mlength), size=1,color="darkgreen") +
  geom_line(data=POPmean, aes(x=age, y=length), size=1,color="grey") +
  xlim(0, 40)+
  ylim(0, 120)+
  labs(x = "Age (years)", y = "Length (cm)", size=15) +
  theme_minimal(base_size = 15)+
  theme(legend.title = element_blank()) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=1, size=15), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=1, size=15), 
        axis.title.x = element_text( size=15),
        axis.title.y = element_text( size=15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank(),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_manual(values=c("slateblue1","tomoato3", "grey"))+
  theme(legend.position="bottomright")

simboot<-sim_boot+
  geom_ribbon(data=POSmean, aes(x=age,ymin=lower, ymax=upper) ,fill="tomato3", alpha=0.2)+
  geom_ribbon(data=LOSmean, aes(x=age,ymin=lower, ymax=upper) ,fill="slateblue1", alpha=0.2)+
  geom_ribbon(data=Ranmean, aes(x=age,ymin=lower, ymax=upper) ,fill="darkgreen", alpha=0.2)


##see each individual run
LOSdata$run<-as.factor(LOSdata$run)
POSdata$run<-as.factor(POSdata$run)
POPdata$run<-as.factor(POPdata$run)
randata$run<-as.factor(POPdata$run)

#plot them together

ggsave(sim_boot, file="sim_boot_unconstrained_0.1Fm_dome_onaga.png")


#####constained BOOT############
#POS
df_total = data.frame()
t0=0
age<-seq(0,40, by =1)
mod<-list()
for(x in 1:100){
    POS[[x]]<-sample_n(population, 300, replace=FALSE, weight=prop1)
    mod[[x]]<- nls(length ~ Linf* (1-exp(-(K*(age-t0)))),
                   data=POS[[x]], 
                   start = list (Linf = 100,
                                 K = 0.13
                            ))
    print(summary(mod[[x]]))
    pred.length_POS <- coef(mod[[x]])[1]*(1 - exp(-coef(mod[[x]])[2]*(age -0)))
    age= age
    length= pred.length_POS
    df=data.frame(age, length)
    df$run<-x
    df_total <- rbind(df_total,df)
    
  }

POSdata<-df_total

#constrained model
POS_coef<-(lapply(mod,coef))
test2<-melt(data=POS_coef, "Linf", "K1")


even_indexes<-seq(2,200,2)
odd_indexes<-seq(1,199,2)

#get means and (95%CI)
Linf_test <- data.frame(x=test2[odd_indexes,1])
K_test <- data.frame(x=test2[even_indexes,1])
mean(Linf_test$x)
quantile(Linf_test$x, prob=0.025)
quantile(Linf_test$x, prob=0.975) 

mean(K_test$x)
quantile(K_test$x, prob=0.025)
quantile(K_test$x, prob=0.975) 






#FOS (LOS)
t0=0
df_total = data.frame()
age<-seq(0,40, by =1)
mod<-list()
for(x in 1:100){
    FOS[[x]] <- population %>%dplyr::group_by(binl)%>% sample_n(pmin(n(), n)) #takes the n specified above or if there are not enough in the bin than the total number of samples
    mod[[x]]<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
                   data=FOS[[x]], 
                   start = list (Linf = 100,
                                 K1 = 0.13))
    print(summary(mod[[x]]))
    pred.length_POS <- coef(mod[[x]])[1]*(1 - exp(-coef(mod[[x]])[2]*(age -0)))
    age= age
    length= pred.length_POS
    df=data.frame(age, length)
    df$run<-x
    df_total <- rbind(df_total,df)
    
  }
LOSdata <- df_total
LOS_coef<-(lapply(mod,coef))



test2<-melt(data=LOS_coef, "Linf", "K1")

even_indexes<-seq(2,200,2)
odd_indexes<-seq(1,199,2)

#get means and (95%CI)
Linf_test <- data.frame(x=test2[odd_indexes,1])
K_test <- data.frame(x=test2[even_indexes,1])
mean(Linf_test$x)
quantile(Linf_test$x, prob=0.025)
quantile(Linf_test$x, prob=0.975) 

mean(K_test$x)
quantile(K_test$x, prob=0.025)
quantile(K_test$x, prob=0.975) 





###sim. pop bootstrapped von b
##for population with differing mincat
str(population)
subset

df_total = data.frame()
age<-seq(0,40, by =1)
mod<-list()
for (i in 1:100) {
  mod[[i]]<- nls(length ~ Linf* (1-exp(-(K1*(age-t0)))),
                 data=population_true, 
                 start = list (Linf = 100,
                               K1 = 0.13))
  print(summary(mod[[i]]))
  pred.length_pop <- coef(mod[[i]])[1]*(1 - exp(-coef(mod[[i]])[2]*(age -0)))
  age= age
  length= pred.length_pop
  df=data.frame(age, length)
  df$run<-i
  df_total <- rbind(df_total,df)
  
}
POPdata <- df_total

POP_coef<-(lapply(mod,coef))

test2<-melt(data=POP_coef, "Linf", "K1")

even_indexes<-seq(2,200,2)
odd_indexes<-seq(1,199,2)


#get means and (95%CI)
Linf_test <- data.frame(x=test2[odd_indexes,1])
K_test <- data.frame(x=test2[even_indexes,1])
mean(Linf_test$x)
quantile(Linf_test$x, prob=0.025)
quantile(Linf_test$x, prob=0.975) 

mean(K_test$x)
quantile(K_test$x, prob=0.025)
quantile(K_test$x, prob=0.975) 


###plot it 
LOSdata$run<-as.factor(LOSdata$run)
POSdata$run<-as.factor(POSdata$run)
POPdata$run<-as.factor(POPdata$run)
Randata$run<-as.factor(Randata$run)


#plot them together
sim_boot<-ggplot() + 
  geom_line(data=POSdata, aes(x=age, y=length, group=run), color='tomato3')+
  geom_line(data=LOSdata, aes(x=age, y=length, group=run), color='slateblue1') + 
  geom_line(data=POPdata, aes(x=age, y=length), size=1.5,color="grey") +

  xlim(0, 40)+
  ylim(0, 120)+
  labs(x = "Age (years)", y = "Fork Length (cm)", size=15) +
  theme_minimal(base_size = 15)+
  theme(legend.title = element_blank()) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=1, size=15), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=1, size=15), 
        axis.title.x = element_text( size=15),
        axis.title.y = element_text( size=15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank(),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_manual(values=c("slateblue1","tomoato3", "grey"))+
  theme(legend.position="bottomright")


ggsave(sim_boot, file="sim_boot_fish_constratined.png")
