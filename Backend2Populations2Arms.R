# Loading libraries
library(mvtnorm)

# This R file creates the necessary backend files for the optimizer to call
# There are three key functions 
# 1) construct.test.statistics.joint.distribution creates mean and covariance 
# matrices associated with the vector of statistics
# 2) get.eff.bound calculates the efficacy boundaries for the design
# 3) design.evaluate for different vectors of test statistics calculates
# which hypothesis are rejected and at which stage.

# Throughout the sequence of test statistics is given by blocks of stages and within a block 
# of a stage k the vector of test statistics is given by 
# (Z_{1,1,k}, Z_{1,2,k}, Z_{2,1,k}, Z_{2,2,k}). Here, the first subscript
# indicates treatment, the second sub-populaiton and the third stage.



### is.scalar ##################################################################
# Description: determine if object is scalar.
# Input:
#   x (object): object to test
# Output:
#   (logical): is object a scalar
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

### is.whole.number ############################################################
# Description: determine if numeric value is a whole number 
#   (see ?is.integer) for more information.
# Input:
#   x (numeric): numeric value
# Output:
#   (logical): is whole number to a given tolerance
is.whole.number <-
function(x, tol=.Machine$double.eps^0.5) {
  stopifnot(is.numeric(x))
  abs(x - round(x)) < tol
}

### in.open.interval ###########################################################
# Description: determine if numeric value is in an open interval (lb, ub)
# Input:
#   x (numeric): numeric value
# Output:
#   (logical): is whole number in the interval
in.open.interval <-
function(x, lb=0, ub=1) {
  stopifnot(is.numeric(x), is.numeric(lb), is.numeric(ub), ub>lb)
  (x > lb) & (x < ub)
}

### reals.to.probability #######################################################
# Description: map matrix of reals to a vector of probabilities by a logistic
#   transformation and normalization.
# Input:
#   x (numeric matrix): real values to be mapped
# Output:
#   (numeric matrix): non-negative values in (0, 1)
reals.to.probability <- function(x) plogis(x)/sum(plogis(x))

### squash #####################################################################
# Description: map matrix of reals to an interval (min, max) using a linear
#   threshold (method="linear"), a logistic function (method="logit"),
#   or other (tbd) methods.
# Input:
#   x (numeric matrix): real values to be mapped to (min, max)
# Output:
#   (numeric matrix): values in (min, max)
squash <- function(x, x.min=0, x.max=1,
 method=c("linear", "logit") [1]) {
  stopifnot(is.numeric(x), 
    x.max>x.min)
  switch(method,
   linear=pmin(pmax(x.min, x), x.max), # Linear Interpolation 
   logit=x.min + (x.max-x.min)*plogis(x, location=(x.max+x.min)/2,
    scale=(x.max-x.min)/6)
   
   )
}

# This function calculates the covariate matrix for a binary and a continuous outcome
# It assumes two treatments and 
# a common control, two sub-populations, and arbitrary number of stages
# prop.samp.vec.pop.1 
# Inputs: var.vec.pop.1: variance vector for population 1
#  var.vec.pop.2: variance vector for population 2
#  prop.samp.vec.pop.1: the proportion of total number of subjects in sub-population one
#  which are enrolled at each stage. E.g. c(0.5, 0.5) means 50 \% of the obs are enrolled 
#  at stage one and 50% at stage 2.
#  prop.samp.vec.pop.1: the proportion of total number of subjects in sub-population one
#  which are enrolled at each stage. E.g. c(0.5, 0.5) means 50 \% of the obs are enrolled 
#  at stage one and 50% at stage 2.

# Output: The covariance matrix associated with the vector of test statistics

cov.mat.cont.bin = function(var.vec.pop.1, var.vec.pop.2,
  prop.samp.vec.pop.1, prop.samp.vec.pop.2){
  # K is the number of stages
  K = length(prop.samp.vec.pop.1)
  # Proportion of sample size enrolled up to and until stage k in each population
  cumsum.prop.samp.pop.1 = cumsum(prop.samp.vec.pop.1)
  cumsum.prop.samp.pop.2 = cumsum(prop.samp.vec.pop.2)
  cov.mat = matrix(0, nrow = 2 * 2 * K, ncol = 2 * 2 * K)
  
  # Filling in the covariance matrix
  
  # Filling in by blocks
  for(i in 1:K){
    for(j in 1:K){
      min.max.term.pop.1 =
      sqrt(min(cumsum.prop.samp.pop.1[i], cumsum.prop.samp.pop.1[j])/
       max(cumsum.prop.samp.pop.1[i], cumsum.prop.samp.pop.1[j]))
      min.max.term.pop.2 =
      sqrt(min(cumsum.prop.samp.pop.2[i], cumsum.prop.samp.pop.2[j])/
       max(cumsum.prop.samp.pop.2[i], cumsum.prop.samp.pop.2[j]))
      sigma.term.pop.1 =
      var.vec.pop.1[1]/sqrt((var.vec.pop.1[1] + var.vec.pop.1[2]) *
        (var.vec.pop.1[1] + var.vec.pop.1[3]))
      sigma.term.pop.2 =
      var.vec.pop.2[1]/sqrt((var.vec.pop.2[1] + var.vec.pop.2[2]) *
        (var.vec.pop.2[1] + var.vec.pop.2[3]))
      cov.mat[(i-1)*4 + 1, (j-1)*4 + 1] = min.max.term.pop.1
      cov.mat[(i-1)*4 + 1, (j-1)*4 + 3] = min.max.term.pop.1 * sigma.term.pop.1
      cov.mat[(i-1)*4 + 2, (j-1)*4 + 2] = min.max.term.pop.2
      cov.mat[(i-1)*4 + 2, (j-1)*4 + 4] = min.max.term.pop.2 * sigma.term.pop.2
      cov.mat[(i-1)*4 + 3, (j-1)*4 + 1] = min.max.term.pop.1 * sigma.term.pop.1
      cov.mat[(i-1)*4 + 3, (j-1)*4 + 3] = min.max.term.pop.1
      cov.mat[(i-1)*4 + 4, (j-1)*4 + 2] = min.max.term.pop.2 * sigma.term.pop.2
      cov.mat[(i-1)*4 + 4, (j-1)*4 + 4] = min.max.term.pop.2
    }
  }
  return(cov.mat)
}


# This function calculates the covariance matrix associatec with 
# a survival outcome. 
# Input d.l.j, l = 0,1,2, and j = 1,2. Here, d.l.j is a vector of length
# K where element k is the expected number of deaths at or before analysis
# k in subpopulation j and treatment l.

# Output covariance matrix associate with vector of test statistics

cov.mat.surv = function(d.0.1, d.1.1, d.2.1, d.0.2, d.1.2, d.2.2){

  K = length(d.0.1)
  cov.mat = matrix(0, nrow = 2 * 2 * K, ncol = 2 * 2 * K)

# Filling in by blocks
  for(i in 1:K){
   for(j in 1:K){
    min.max.term.pop.1.treatment.1 = sqrt((d.0.1[min(i,j)] + d.1.1[min(i,j)])/(d.0.1[max(i,j)] + d.1.1[max(i,j)]))
    min.max.term.pop.1.treatment.2 = sqrt((d.0.1[min(i,j)] + d.2.1[min(i,j)])/(d.0.1[max(i,j)] + d.2.1[max(i,j)]))
    min.max.term.pop.2.treatment.1 = sqrt((d.0.2[min(i,j)] + d.1.2[min(i,j)])/(d.0.2[max(i,j)] + d.1.2[max(i,j)]))
    min.max.term.pop.2.treatment.2 = sqrt((d.0.2[min(i,j)] + d.2.2[min(i,j)])/(d.0.2[max(i,j)] + d.2.2[max(i,j)]))
      # Both treatment 1 and subpopulation 1 different stages
    cov.mat[(i-1) * 4 + 1, (j-1) * 4 + 1] = min.max.term.pop.1.treatment.1
      # Different treatments and same subpopulation 1 different stages
    cov.mat[(i-1) * 4 + 1, (j-1) * 4 + 3] = mean(min.max.term.pop.1.treatment.1, min.max.term.pop.1.treatment.2) * 0.5
      # Both treatment 1 and subpopulation 2 different stages
    cov.mat[(i-1) * 4 + 2, (j-1) * 4 + 2] = min.max.term.pop.2.treatment.1
      # Different treatments and subpopulation 2 different stages      
    cov.mat[(i-1) * 4 + 2, (j-1) * 4 + 4] = mean(min.max.term.pop.2.treatment.1, min.max.term.pop.2.treatment.2) * 0.5
      # Different treatments and subpopulation 1 different stages      
    cov.mat[(i-1) * 4 + 3, (j-1) * 4 + 1] = mean(min.max.term.pop.1.treatment.1, min.max.term.pop.1.treatment.2) * 0.5
      # Same treatment 2 and subpopulation 1 different stages      
    cov.mat[(i-1) * 4 + 3, (j-1) * 4 + 3] = min.max.term.pop.1.treatment.2
      # Different treatments 1 and subpopulation 2 different stages      
    cov.mat[(i-1) * 4 + 4, (j-1) * 4 + 2] =  mean(min.max.term.pop.2.treatment.1, min.max.term.pop.2.treatment.2) * 0.5
      # Same treatment 2 and subpopulation 2 different stages      
    cov.mat[(i-1) * 4 + 4, (j-1) * 4 + 4] = min.max.term.pop.2.treatment.2
  }
}
# Throughout the sequence of test statistics is given by blocks of stages 
# and within a block  of a stage k the vector of test statistics is given by 
# (Z_{1,1,k}, Z_{1,2,k}, Z_{2,1,k}, Z_{2,2,k}), where the first subscript
# indicates treatment, the second sub-populaiton and the third stage.

return(cov.mat)
}


# This function creates the covariance matrix and mean vector 
# associated with the test statistic
# Inputs: #   analytic.n.per.stage - [K x J(L+1)] matrix: patients with primary outcome
#             at each interim analysis.
#             stage 1: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#             stage 2: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#                ...
#             stage k: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#         outcome.type: type of outcome, one of continuous, binary or survival
#         mean.sub.pop.1: the assumed means assocated with each treatment in sub-population 1
#         the mean vector is input in the order (control, treatment 1, treatment 2) 
#         mean.sub.pop.2: the assumed means assocated with each treatment in sub-population 1
#         the mean vector is input in the order (control, treatment 1, treatment 2) 
#         var.vec.pop.1: the variance vector associated with each treatment in sub-population one
#         var.vec.pop.2: the variance vector associated with each treatment in sub-population two 
#         prop.pop.1: The proportion of subjects in population one. Assumed known.
#         max.follow: For survival outcome, how long each participant is followed up
#         enrollment.period: For survival outcome, the maximum time participants are enrolled
#         hazard.rate.pop.1: For survival outcome, hazard rate for subpopulation 1 in 
#         the order (control, treatment 1, treatment 2)
#         hazard.rate.pop.2: For survival outcome, hazard rate for subpopulation 2 in 
#         the order (control, treatment 1, treatment 2)
#         time: for time-to-event outcome, time is the timing of all analysis.
#         relative.efficiency: ratio of the asymptotic variance of unadjusted estimator to 
#         asymptotic variance of adjusted estimator.
#         censoring.rate: For a survival outcome only. It is the proportion of participants that are not
#         administratively censored which drop out of the study. For example, if 100 events are expected
#         without any dropout then setting censoring rate to 0.5 means that 100*0.5 events are expected



# Output: A list with three elements:
#         cov.mat.used: Covariance matrix associated with test statistic.
#         non.centrality.parameter.vec = The mean vector associated with each test statistic.
#         information.vector =  for a given stage k, elements [(1+ (k-1) * 4):(4 + (k-1) * 4)] are 
#         (var(\beta_{1,1,k}, var(\beta_{1,2,k}), var(\beta_{2,1,k}), var(beta_{2,2,k}) where the
#         first subscript indicates treatment, the second sub-population and the third stage. 
#         beta is the estimator of the treatment effect

construct.test.statistics.joint.distribution <- function(analytic.n.per.stage,
  mean.sub.pop.1=NULL,
  mean.sub.pop.2=NULL,
  var.vec.pop.1=NULL,
  var.vec.pop.2=NULL,
  outcome.type,
  prop.pop.1,
  max.follow = NULL,
  enrollment.period = NULL,
  hazard.rate.pop.1 = NULL,
  hazard.rate.pop.2 = NULL,
  time = NULL,
  censoring.rate = NULL,
  relative.efficiency = NULL){
  # Number of stages
  K <- nrow(analytic.n.per.stage)
  
  # Calculating the total number of subjects at each analysis for both sub-populations
  # Note: We assume that an equal number is enrolled to treatment and control.
  n.pop.1 = analytic.n.per.stage[, 1]
  n.pop.2 = analytic.n.per.stage[, 2]
  
  # Calculating the proportion of observation sampled at each stage for the
  # Two treatments
  prop.samp.vec.pop.1 = diff(c(0, n.pop.1))/n.pop.1[K]
  prop.samp.vec.pop.2 = diff(c(0, n.pop.2))/n.pop.2[K]
  
  
  # Creating storage space for mean vector
  mean.vec = rep(NA, 4 * K)
  
  # Do the calculations seperately depending on the type of outcome
  if(outcome.type == "continuous"){
    for(i in 1:K){
      mean.vec[((i-1)*4+1):((i-1)*4+4)] = 
      c(sqrt(n.pop.1[i])*(mean.sub.pop.1[2]-mean.sub.pop.1[1])/
        sqrt(var.vec.pop.1[2]+var.vec.pop.1[1]), 
        sqrt(n.pop.2[i])*(mean.sub.pop.2[2] - mean.sub.pop.2[1])/
        sqrt(var.vec.pop.2[2]+var.vec.pop.2[1]),
        sqrt(n.pop.1[i])*(mean.sub.pop.1[3] - mean.sub.pop.1[1])/
        sqrt(var.vec.pop.1[3]+var.vec.pop.1[1]),
        sqrt(n.pop.2[i])*(mean.sub.pop.2[3] - mean.sub.pop.2[1])/
        sqrt(var.vec.pop.2[3]+var.vec.pop.2[1]))
    }
    cov.mat.used = cov.mat.cont.bin(var.vec.pop.1,
      var.vec.pop.2,
      prop.samp.vec.pop.1,
      prop.samp.vec.pop.2)
  }
  
  if(outcome.type == "binary"){
    var.vec.pop.1 = mean.sub.pop.1*(1 - mean.sub.pop.1)
    var.vec.pop.2 = mean.sub.pop.2*(1 - mean.sub.pop.2)
    for(i in 1:K){
      mean.vec[((i-1)*4+1):((i-1)*4+4)] =
      c(sqrt(n.pop.1[i])*(mean.sub.pop.1[2] - mean.sub.pop.1[1])/
        sqrt(var.vec.pop.1[2]+var.vec.pop.1[1]),
        sqrt(n.pop.2[i])*(mean.sub.pop.2[2] - mean.sub.pop.2[1])/
        sqrt(var.vec.pop.2[2]+var.vec.pop.2[1]),
        sqrt(n.pop.1[i])*(mean.sub.pop.1[3] - mean.sub.pop.1[1])/
        sqrt(var.vec.pop.1[3]+var.vec.pop.1[1]),
        sqrt(n.pop.2[i])*(mean.sub.pop.2[3] - mean.sub.pop.2[1])/
        sqrt(var.vec.pop.2[3]+var.vec.pop.2[1]))
    }
    cov.mat.used = cov.mat.cont.bin(var.vec.pop.1,
      var.vec.pop.2,
      prop.samp.vec.pop.1,
      prop.samp.vec.pop.2)
  }
  

  if(outcome.type == "survival"){

# The ratios of log-rank tests
    theta = -c(log(hazard.rate.pop.1[2]/hazard.rate.pop.1[1]), log(hazard.rate.pop.2[2]/hazard.rate.pop.2[1]), log(hazard.rate.pop.1[3]/hazard.rate.pop.1[1]), log(hazard.rate.pop.2[3]/hazard.rate.pop.2[1]))

# information vector in same order as test statistics mentioned above.
    mean.vec = rep(NA, 4 * K)

    d.0.1 = rep(0,K) # number of deaths in control group pop 1
    d.1.1 = rep(0,K) # number of deaths in treatment group 1 pop 1
    d.2.1 = rep(0,K) # number of deaths in treatment group 2 pop 1
    d.0.2 = rep(0,K) # number of deaths in control group pop 2
    d.1.2 = rep(0,K) # number of deaths in treatment group 1 pop 2
    d.2.2 = rep(0,K) # number of deaths in treatment group 2 pop 2

    for(i in 1:K){

    # Calculating the expected number of deaths for each treatment + sub-population combination 
    # at interim analys i

    # We cycle through the 6 different cases
      if(enrollment.period >= time[i] & time[i] >= max.follow){
       d.0.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[1] * max.follow))/(max.follow * hazard.rate.pop.1[1]))
       d.1.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[2] * max.follow))/(max.follow * hazard.rate.pop.1[2]))
       d.2.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[3] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[3] * max.follow))/(max.follow * hazard.rate.pop.1[3]))
       d.0.2[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.2[1] * max.follow))/(max.follow * hazard.rate.pop.2[1]))
       d.1.2[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.2[2] * max.follow))/(max.follow * hazard.rate.pop.2[2]))
       d.2.2[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.2[3] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.2[3] * max.follow))/(max.follow * hazard.rate.pop.2[3]))
     }
     
     if(enrollment.period >= max.follow & max.follow >= time[i]){
       d.0.1[i] = 1 - (1-exp(-hazard.rate.pop.1[1] * time[i]))/(time[i] * hazard.rate.pop.1[1])
       d.1.1[i] = 1 - (1-exp(-hazard.rate.pop.1[2] * time[i]))/(time[i] * hazard.rate.pop.1[2])
       d.2.1[i] = 1 - (1-exp(-hazard.rate.pop.1[3] * time[i]))/(time[i] * hazard.rate.pop.1[3])
       d.0.2[i] = 1 - (1-exp(-hazard.rate.pop.2[1] * time[i]))/(time[i] * hazard.rate.pop.2[1])
       d.1.2[i] = 1 - (1-exp(-hazard.rate.pop.2[2] * time[i]))/(time[i] * hazard.rate.pop.2[2])
       d.2.2[i] = 1 - (1-exp(-hazard.rate.pop.2[3] * time[i]))/(time[i] * hazard.rate.pop.2[3])
     }

     if(time[i] >= enrollment.period & enrollment.period >= max.follow){
       k = time[i] - enrollment.period
       d.0.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[1] * k) - exp(-hazard.rate.pop.1[1] * max.follow))/((max.follow - k) * hazard.rate.pop.1[1]))
       d.1.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[2] * k) - exp(-hazard.rate.pop.1[2] * max.follow))/((max.follow - k) * hazard.rate.pop.1[2]))
       d.2.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[3] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[3] * k) - exp(-hazard.rate.pop.1[3] * max.follow))/((max.follow - k) * hazard.rate.pop.1[3]))
       d.0.2[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[1] * k) - exp(-hazard.rate.pop.2[1] * max.follow))/((max.follow - k) * hazard.rate.pop.2[1]))
       d.1.2[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[2] * k) - exp(-hazard.rate.pop.2[2] * max.follow))/((max.follow - k) * hazard.rate.pop.2[2]))
       d.2.2[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[3] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[3] * k) - exp(-hazard.rate.pop.2[3] * max.follow))/((max.follow - k) * hazard.rate.pop.2[3]))
     }

     if(max.follow >= enrollment.period & enrollment.period >= time[i]){
       d.0.1[i] = 1 - (1-exp(-hazard.rate.pop.1[1] * time[i]))/(time[i] * hazard.rate.pop.1[1])
       d.1.1[i] = 1 - (1-exp(-hazard.rate.pop.1[2] * time[i]))/(time[i] * hazard.rate.pop.1[2])
       d.2.1[i] = 1 - (1-exp(-hazard.rate.pop.1[3] * time[i]))/(time[i] * hazard.rate.pop.1[3])
       d.0.2[i] = 1 - (1-exp(-hazard.rate.pop.2[1] * time[i]))/(time[i] * hazard.rate.pop.2[1])
       d.1.2[i] = 1 - (1-exp(-hazard.rate.pop.2[2] * time[i]))/(time[i] * hazard.rate.pop.2[2])
       d.2.2[i] = 1 - (1-exp(-hazard.rate.pop.2[3] * time[i]))/(time[i] * hazard.rate.pop.2[3])
     }

     if(max.follow >= time[i] & time[i] >= enrollment.period){
       d.0.1[i] = 1 - (exp(-hazard.rate.pop.1[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[1] * time[i]))/(enrollment.period * hazard.rate.pop.1[1])
       d.1.1[i] = 1 - (exp(-hazard.rate.pop.1[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[2] * time[i]))/(enrollment.period * hazard.rate.pop.1[2])
       d.2.1[i] = 1 - (exp(-hazard.rate.pop.1[3] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[3] * time[i]))/(enrollment.period * hazard.rate.pop.1[3])
       d.0.2[i] = 1 - (exp(-hazard.rate.pop.2[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[1] * time[i]))/(enrollment.period * hazard.rate.pop.2[1])
       d.1.2[i] = 1 - (exp(-hazard.rate.pop.2[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[2] * time[i]))/(enrollment.period * hazard.rate.pop.2[2])
       d.2.2[i] = 1 - (exp(-hazard.rate.pop.2[3] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[3] * time[i]))/(enrollment.period * hazard.rate.pop.2[3])
     }

     if(time[i] >= max.follow & max.follow >= enrollment.period){
       d.0.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[1] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[1]))
       d.1.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[2] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[2]))
       d.2.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[3] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[3] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[3] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[3]))
       d.0.2[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[1] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.2[1]))
       d.1.2[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[2] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.2[2]))
       d.2.2[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[3] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[3] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[3] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.2[3]))
     }

     d.0.1[i] = d.0.1[i] * n.pop.1[i] * (1 - censoring.rate)
     d.1.1[i] = d.1.1[i] * n.pop.1[i] * (1 - censoring.rate)
     d.2.1[i] = d.2.1[i] * n.pop.1[i] * (1 - censoring.rate)
     d.0.2[i] = d.0.2[i] * n.pop.2[i] * (1 - censoring.rate)
     d.1.2[i] = d.1.2[i] * n.pop.2[i] * (1 - censoring.rate)
     d.2.2[i] = d.2.2[i] * n.pop.2[i] * (1 - censoring.rate)

    # Calculating the information and covariance matrix
     mean.vec[((i-1) * 4 + 1):((i-1) * 4 + 4)] = - theta * sqrt(c((d.0.1[i] + d.1.1[i])/4, (d.0.2[i] + d.1.2[i])/4, (d.0.1[i] + d.2.1[i])/4, (d.0.2[i] + d.2.2[i])/4))
   }
   
   cov.mat.used = cov.mat.surv(d.0.1, d.1.1, d.2.1, d.0.2, d.1.2, d.2.2)
 }


# Create information vector on the extimator scale
# for a given stage k elements [(1+ (k-1) * 4):(4 + (k-1) * 4)] are 
# (var(\beta_{1,1,k}, var(\beta_{1,2,k}), var(\beta_{2,1,k}), var(beta_{2,2,k})
#, where the first subscript indicates treatment, the second sub-populaiton and the third stage.

# For a continuous outcome 
 if(outcome.type == "continuous"){
# Initialize the vector
  information.vector.inv = rep(NA, 4 * K)

  for(i in 1:K){
    information.vector.inv[((i-1)*4+1):((i-1)*4+4)] = 
    c(1/n.pop.1[i] * (var.vec.pop.1[2]+var.vec.pop.1[1]), 
      1/n.pop.2[i] * (var.vec.pop.2[2]+var.vec.pop.2[1]),
      1/n.pop.1[i] * (var.vec.pop.1[3]+var.vec.pop.1[1]),
      1/n.pop.2[i] * (var.vec.pop.2[3]+var.vec.pop.2[1]))
  }
} # End if outcome.type is continuous

# For a binary outcome 
if(outcome.type == "binary"){
# Initialize the vector
  information.vector.inv = rep(NA, 4 * K)

  var.vec.pop.1 = mean.sub.pop.1*(1 - mean.sub.pop.1)
  var.vec.pop.2 = mean.sub.pop.2*(1 - mean.sub.pop.2)

  for(i in 1:K){
    information.vector.inv[((i-1)*4+1):((i-1)*4+4)] = 
    c(1/n.pop.1[i] * (var.vec.pop.1[2]+var.vec.pop.1[1]), 
      1/n.pop.2[i] * (var.vec.pop.2[2]+var.vec.pop.2[1]),
      1/n.pop.1[i] * (var.vec.pop.1[3]+var.vec.pop.1[1]),
      1/n.pop.2[i] * (var.vec.pop.2[3]+var.vec.pop.2[1]))
  }
} # End if outcome.type is binary

if(outcome.type == "survival"){

  # Initialize the vector
  information.vector.inv = rep(NA, 4 * K)

  for(i in 1:K){
    information.vector.inv[((i-1)*4+1):((i-1)*4+4)] = 
    c(4/((d.0.1[i] + d.1.1[i])), 
      4/((d.0.2[i] + d.1.2[i])),
      4/((d.0.1[i] + d.2.1[i])),
      4/((d.0.2[i] + d.2.2[i])))
  }

}

if(!is.null(relative.efficiency)){
  mean.vec = mean.vec * sqrt(relative.efficiency)
  information.vector.inv = information.vector.inv/relative.efficiency
}

information.vector.inv.matrix = matrix(information.vector.inv, nrow = K, byrow = TRUE)



return(list(cov.mat.used=cov.mat.used,
  non.centrality.parameter.vec = mean.vec,
  information.vector = 1/information.vector.inv.matrix))
}



# This function calculates the efficacy boundaries
# Input: alpha.alloc = Vector of alpha allocations:
#        The alpha allocation vector is of length 2 * number of stages
#        the first two elements are the alpha allocations at stage one to each subpopulation
#        the next two elements are the alpha allocations at the next stage and so on
#        cov.mat.used: covariance matrix under the scenario
#        err.tol = how precise is the binary search

# Output:eff.bound: The vector u_{j,k} with blocks corresponding to stages and (u_{1,k}, u_{2,k}) within stage
#        eff.bound.z: The vector z_{j,k} with blocks corresponding to stages and (z_{1,k}, z_{2,k}) within stage
#        eff.bound.alpha: K blocks where each block is (\tilde u_{j,K}, \tilde u_{j,K}) where block
#        k corresponds to alpha reallocated if both treatments in other sub-pop are rejected at stage k
#        eff.bound.z.alpha: K blocks where each block is (\tilde z_{j,K}, \tilde z_{j,K}) where block
#        k corresponds to alpha reallocated if both treatments in other sub-pop are rejected at stage k


get.eff.bound = function(alpha.allocation, cov.mat.used, err.tol = 10^-3){
  
  # Number of stages
  K = length(alpha.allocation)/2
  
  # Getting index corresponding to which sup-population is being used in each 
  # alpha allocation
  index.sub.pop = rep(c(1, 2), K)
  
  # eff.bound is the vector of efficacy boundaries with the elements corresponding to the
  # same stage and subpopulation combinations as in the alpha.allocation vector
  eff.bound = rep(NA, 2 * K)
  eff.bound.z = rep(NA, 2 * K)

  # Cumulative alpha allocation with subpopulation 1 and 2
  cum.alpha.1 = cumsum(alpha.allocation[which(index.sub.pop == 1)])
  cum.alpha.2 = cumsum(alpha.allocation[which(index.sub.pop == 2)])
  
  # Calculating the first elements of the efficacy boundary u_{j,1} corresponding to each 
  # subpopulation
  temp.1 = rmvnorm(10^6, mean = rep(0, 2), sigma = cov.mat.used[1:4, 1:4][c(1,3), c(1,3)])
  temp.2 = rmvnorm(10^6, mean = rep(0, 2), sigma = cov.mat.used[1:4, 1:4][c(2,4), c(2,4)])
  temp.1.max = apply(temp.1, 1, max)
  temp.2.max = apply(temp.2, 1, max)
  eff.bound[1] = quantile(temp.1.max, 1-alpha.allocation[1])
  eff.bound[2] = quantile(temp.2.max, 1-alpha.allocation[2])
  
  # Calculating the first elements of the efficacy boundary z_{j,1} 
  eff.bound.z[1] = qnorm(1-alpha.allocation[1])
  eff.bound.z[2] = qnorm(1-alpha.allocation[2]) 


  # A function that calculates the cumulative type one error corresponding 
  # to a sub-population j
  # effecacy boundaries is the efficacy boundary
  # cov.mat.used is the covariance matrix
  # j is the subpopulation
  sign.lev = function(eff.bound.used, cov.mat.used, sub.pop.numb){
    
    numb.stages = length(eff.bound.used)
    
    # Index which test-statistic belongs to population and treatment, respectivly
    index.sub.pop.eff = rep(c(1, 2), 2 * numb.stages)
    index.treatment.eff = rep(c(1,1, 2,2), numb.stages)
    
    
    cov.mat.eff.bound =
    cov.mat.used[1:(4* numb.stages), 1:(4* numb.stages)][which(index.sub.pop.eff == sub.pop.numb), which(index.sub.pop.eff == sub.pop.numb)]
    # Calculating the overall type one error under the global null
    type.1.err = 1- pmvnorm(mean = rep(0, 2 * numb.stages), sigma= cov.mat.eff.bound,lower = rep(-Inf, 2 * numb.stages), upper= rep(eff.bound.used, each = 2), algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))[1]
    
    return(type.1.err)
  }

  sign.lev.z = function(eff.bound.used, cov.mat.used, sub.pop.numb){
    
    numb.stages = length(eff.bound.used)
    
    # Index which test-statistic belongs to population and treatment, respectivly
    index.sub.pop.eff = rep(c(1, 2), 2 * numb.stages)

    # Getting the treatment assignment
    index.treatment.eff = rep(c(1,1,2,2), numb.stages)
    
    # Only use treatment 1 and subpopulation of interest
    index.used = which(index.treatment.eff == 1 & index.sub.pop.eff == sub.pop.numb)
    
    cov.mat.eff.bound = cov.mat.used[index.used, index.used]

   # Calculating the overall type one error 
    type.1.err = 1- pmvnorm(mean = rep(0, numb.stages), sigma= cov.mat.eff.bound,lower = rep(-Inf, numb.stages), upper= eff.bound.used, algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))[1]
    
    return(type.1.err)
  }

  
  
  # Start calculating the efficacy boundaries associated with population 1
  # Cycling through the stages after and calculating the efficacy boundary u_{1,k} for k = 2, \ldots, K
  if(K > 1){
    for(i in 2:K){

      # Start by doing binary search for u_{1, k}
      # upper and lower values of interval
      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound[index.sub.pop == 1]
      
      while(length.int > err.tol){
        eff.bound.j[i] = upper.bound.term
        alpha.used = sign.lev(eff.bound.j[1:i], cov.mat.used, 1)
        
        if(alpha.used < cum.alpha.1[i]){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= cum.alpha.1[i]){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound[(i-1)*2 + 1] = upper.bound.term

      # Now do binary search for z_{1,k}
      # upper and lower values of interval

      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound.z[index.sub.pop == 1]
      
      while(length.int > err.tol){
        eff.bound.j[i] = upper.bound.term
        alpha.used = sign.lev.z(eff.bound.j[1:i], cov.mat.used, 1)
        
        if(alpha.used < cum.alpha.1[i]){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= cum.alpha.1[i]){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound.z[(i-1)*2 + 1] = upper.bound.term
    }
    
    # Calculate the efficacy boundaries associated with population 2
    # Cycling through the stages and calculating the efficacy boundary 
    # u_{2,k}, z_{2,k} for k = 1, \ldots, K
    for(i in 2:K){
      # upper and lower values of interval
      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound[index.sub.pop == 2]
      
      while(length.int > err.tol){
        eff.bound.j[i] = upper.bound.term
        alpha.used = sign.lev(eff.bound.j[1:i], cov.mat.used, 2)
        
        if(alpha.used < cum.alpha.2[i]){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= cum.alpha.2[i]){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound[(i-1)*2 + 2] = upper.bound.term

      # Now do binary search for z_{2,k}
      # upper and lower values of interval

      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound.z[index.sub.pop == 2]
      
      while(length.int > err.tol){
        eff.bound.j[i] = upper.bound.term
        alpha.used = sign.lev.z(eff.bound.j[1:i], cov.mat.used, 2)
        
        if(alpha.used < cum.alpha.2[i]){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= cum.alpha.2[i]){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound.z[(i-1)*2 + 2] = upper.bound.term
    }
  } # End if K >1 statement

  # efficacy boundaries associated with alphar reallocation
  # eff.bound.alpha[2*(k-1) +1] is \tilde u_{1,K} if both H_0 in subpopulation two 
  # are rejected at stage k
  # eff.bound.alpha[2*k] is \tilde u_{2,K} if both H_0 in subpopulation one 
  # are rejected at stage k
  # eff.bound.z.alpha[2*(k-1) +1] is \tilde z_{1,K} if both H_0 in subpopulation two 
  # are rejected at stage k
  # eff.bound.alpha[2*k] is \tilde z_{2,K} if both H_0 in subpopulation one 
  # are rejected at stage k

  eff.bound.alpha = rep(NA, 2 * K)
  eff.bound.z.alpha = rep(NA, 2 * K)
  
  for(k in 1:K){
  # Start with sub-population 1
  # Now we calculate the efficacy boundaries for pop 1 if both null hypothesis corresponding to
  # pop 1 are rejected at stage k

      # Start a binary search for \tilde u_{1,K}
      # upper and lower values of interval
    upper = 20
    lower = -20
    length.int = upper - lower
    
      # Initial guess
    upper.bound.term = mean(c(upper, lower))
    eff.bound.j = eff.bound[index.sub.pop == 1]
      # The cumulative alpha level that the last stage is allowed to test at (note not \alpha_{1,K})
      # \sum_{j=1}^K \alpha_{1,j} + \sum_{j=k}^K \alpha_{2,j}
    alpha.allowed = cum.alpha.1[K] + (cum.alpha.2[K] - c(0,cum.alpha.2)[k])
    
    while(length.int > err.tol){
      eff.bound.j[K] = upper.bound.term
      alpha.used = sign.lev(eff.bound.j, cov.mat.used, 1)

      if(alpha.used < alpha.allowed){
        upper = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, lower))
      }
      if(alpha.used >= alpha.allowed){
        lower = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, upper))
      }
      length.int = upper - lower
    }
    
      # "Rounding" up to preserve type 1 error
    upper.bound.term = upper.bound.term + length.int
    eff.bound.alpha[(k-1) * 2 + 1] = upper.bound.term

    # Start a binary search for \tilde z
      # upper and lower values of interval
    upper = 20
    lower = -20
    length.int = upper - lower
    
      # Initial guess
    upper.bound.term = mean(c(upper, lower))
    eff.bound.j = eff.bound.z[index.sub.pop == 1]
      # The alpha level that the last stage is allowed to test at
      # \alpha_{1,K} + \sum_{j=k}^K \alpha_{2,j}
    alpha.allowed = cum.alpha.1[K] + (cum.alpha.2[K] - c(0,cum.alpha.2)[k])
    
    while(length.int > err.tol){
      eff.bound.j[K] = upper.bound.term
      alpha.used = sign.lev.z(eff.bound.j, cov.mat.used, 1)

      if(alpha.used < alpha.allowed){
        upper = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, lower))
      }
      if(alpha.used >= alpha.allowed){
        lower = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, upper))
      }
      length.int = upper - lower
    }
    
      # "Rounding" up to preserve type 1 error
    upper.bound.term = upper.bound.term + length.int
    eff.bound.z.alpha[(k-1)*2 + 1] = upper.bound.term

  # Now sub-population 2
  # Calculate the efficacy boundaries for pop 2 if both null hypothesis corresponding to
  # pop 1 are rejected at stage k

      # Start a binary search for \tilde u
      # upper and lower values of interval
    upper = 20
    lower = -20
    length.int = upper - lower
    
      # Initial guess
    upper.bound.term = mean(c(upper, lower))
    eff.bound.j = eff.bound[index.sub.pop == 2]
      # The alpha level that the last stage is allowed to test at
      # \alpha_{1,K} + \sum_{j=k}^K \alpha_{2,j}
    alpha.allowed = cum.alpha.2[K] + (cum.alpha.1[K] - c(0,cum.alpha.1)[k])
    
    while(length.int > err.tol){
      eff.bound.j[K] = upper.bound.term
      alpha.used = sign.lev(eff.bound.j, cov.mat.used, 2)

      if(alpha.used < alpha.allowed){
        upper = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, lower))
      }
      if(alpha.used >= alpha.allowed){
        lower = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, upper))
      }
      length.int = upper - lower
    }
    
      # "Rounding" up to preserve type 1 error
    upper.bound.term = upper.bound.term + length.int
    eff.bound.alpha[k*2] = upper.bound.term

    # Start a binary search for \tilde z
      # upper and lower values of interval
    upper = 20
    lower = -20
    length.int = upper - lower
    
      # Initial guess
    upper.bound.term = mean(c(upper, lower))
    eff.bound.j = eff.bound.z[index.sub.pop == 2]
      # The alpha level that the last stage is allowed to test at
      # \alpha_{1,K} + \sum_{j=k}^K \alpha_{2,j}
    alpha.allowed = cum.alpha.2[K] + (cum.alpha.1[K] - c(0,cum.alpha.1)[k])
    
    while(length.int > err.tol){
      eff.bound.j[K] = upper.bound.term
      alpha.used = sign.lev.z(eff.bound.j, cov.mat.used, 2)

      if(alpha.used < alpha.allowed){
        upper = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, lower))
      }
      if(alpha.used >= alpha.allowed){
        lower = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, upper))
      }
      length.int = upper - lower
    }
    
      # "Rounding" up to preserve type 1 error
    upper.bound.term = upper.bound.term + length.int
    eff.bound.z.alpha[k*2] = upper.bound.term
  } # end for k loop
  
  # Make sure that eff.bound >= eff.bound.z
  eff.bound = pmax(eff.bound.z, eff.bound)
  eff.bound.alpha = pmax(eff.bound.z.alpha, eff.bound.alpha)


  return(list(eff.bound = eff.bound, eff.bound.z = eff.bound.z, eff.bound.alpha = eff.bound.alpha, eff.bound.z.alpha = eff.bound.z.alpha))
}



# This function evaluates the performance of a given design
# Inputs: test.statistics: A matrix of test statistics where each row is a vector of test statistics
#         efficacy.boundary: all the different effacacy boundaries outputted from get.eff.bound
#         futility.boundary: A vector of futility boundaries with blocks corresponding to stages
#         and within a block the null hypothesis are ordered as in cov.mat.bin
#         alpha.allocation: The alpha allocation to each stage

# Output: A list consisting of three elements. The first one is a matrix of dim
#         n.sim times 4 where each column corresponds to if H_0 is rejected where the 
#         hypothesis are in the same order as in the covariance matrix. One means rejected 
#         and zero means not rejected.
#         The second element is of the same nature as the first element where
#         each column indicates at what stage the decision to reject or not reject
#         the corresponding hypothesis is made.
#         The third vector is the list of efficacy boundaries as outputted by get.eff.bound
design.evaluate <- function(test.statistics,
  efficacy.boundary,
  futility.boundary,
  alpha.allocation){
  
  # Number of stages
  K = length(alpha.allocation)/2
  
  # Number of MC evaluations
  n.sim <- nrow(test.statistics)
  
  # u_{j,k} boundaries
  est.eff.bound= efficacy.boundary$eff.bound
  # The alpha re-allocated u_{j,k}
  est.eff.bound.alpha = efficacy.boundary$eff.bound.alpha
  # z_{j,k} bounaries
  est.eff.bound.z= efficacy.boundary$eff.bound.z
  # The alpha re-allocated u_{j,k}
  est.eff.bound.alpha.z = efficacy.boundary$eff.bound.z.alpha
  
  # Going through each hypothesis being tested finding out if rejected or not
  # and when rejected/stopped
  
  reject.hyp = matrix(0, nrow = n.sim, ncol = 4)
  stage.decision = matrix(NA, nrow = n.sim, ncol = 4)
  
  for(i in 1:n.sim){
    # Start by looking at subpopulation 1
    index.sub.pop = rep(c(1, 2), 2 * K)
    index.stage = rep(1:K, each = 4)
    is.rejected = FALSE
    
    # Looking at sub-population one at stage one
    index.used = which(index.sub.pop == 1 & index.stage == 1)
    if(max(test.statistics[i, index.used]) > est.eff.bound[index.used[1]]){
      reject.hyp[i, c(1,3)[which.max(test.statistics[i, index.used])]] = 1
      is.rejected = TRUE
      stage.decision[i, c(1,3)[which.max(test.statistics[i, index.used])]] = 1
    }
    
    # If the max test is rejected continue onto the min test
    if(is.rejected & min(test.statistics[i, index.used]) > est.eff.bound.z[1]){
      reject.hyp[i, c(1,3)[which.min(test.statistics[i, index.used])]] = 1
      stage.decision[i, c(1,3)[which.min(test.statistics[i, index.used])]] = 1
    }
    
    # Testing for futility at stage one in sub-population one
    # Finding the treatments that should be stopped for futility
    stage.decision[i, index.used[which(test.statistics[i, index.used] <= futility.boundary[index.used] & reject.hyp[i, index.used] != 1)]] <- 1
    
    
    # Sub-population 2
    index.used = which(index.sub.pop == 2 & index.stage == 1)
    is.rejected = FALSE
    
    if(max(test.statistics[i, index.used]) > est.eff.bound[index.used[1]]){
      reject.hyp[i, index.used[which.max(test.statistics[i, index.used])] ] = 1
      is.rejected = TRUE
      stage.decision[i, index.used[which.max(test.statistics[i, index.used])]] = 1
    }
    
    # If the max test is rejected continue onto the next stage
    if(is.rejected & min(test.statistics[i, index.used]) > est.eff.bound.z[2]){
      reject.hyp[i, index.used[which.min(test.statistics[i, index.used])]] = 1
      stage.decision[i, index.used[which.min(test.statistics[i, index.used])]] = 1
    }
    
    # Testing for futility at stage one in sub-population two
    # Finding the treatments that should be stopped for futility
    stage.decision[i, index.used[which(test.statistics[i, index.used] <= futility.boundary[index.used] & reject.hyp[i, index.used] != 1)]] <- 1
    
    
    # Cycling through the stages
    if(K>1){
      for(k in 2:K){
        
        # Start by sub-population one
        # Look if already stopped at last stage
        index.used = which(index.sub.pop == 1 & index.stage == k)
        
        # Finding which treatments in sub-population one continued onto
        # stage k 
        which.at.stage = which(is.na(stage.decision[i, c(1,3)]))
        is.rejected = FALSE
        # If both continue onto second stage 
        if(length(which.at.stage) == 2){
          if(max(test.statistics[i, index.used]) > est.eff.bound[2 * (k-1) + 1]){
            reject.hyp[i, c(1,3)[which.max(test.statistics[i, index.used])]] = 1
            is.rejected = TRUE
            stage.decision[i, c(1,3)[which.max(test.statistics[i, index.used])]] = k
          }
          
          # If the max test is rejected continue onto test the the other hypothesis
          if(is.rejected & min(test.statistics[i, index.used]) > est.eff.bound.z[2 * (k-1) + 1]){
            reject.hyp[i, c(1,3)[which.min(test.statistics[i, index.used])]] = 1
            stage.decision[i, c(1,3)[which.min(test.statistics[i, index.used])]] = k
          }
        }
        
        # If only one continues onto second stage 
        if(length(which.at.stage) == 1){
          # Finding the other one
          not.at.stage = setdiff(1:2, which.at.stage)
          
          # If the other one stopped for futility at earlier stages
          if(reject.hyp[i, c(1,3)[not.at.stage]] == 0)
            if(test.statistics[i, index.used[which.at.stage]] > est.eff.bound[2 * (k-1) + 1]){
              reject.hyp[i, c(1,3)[which.at.stage]] = 1
              stage.decision[i, c(1,3)[which.at.stage]] = k
            }
            
          # If the other one stopped for efficacy at earlier stages use the step down procedure
            if(reject.hyp[i, c(1,3)[not.at.stage]] == 1)
              if(test.statistics[i, index.used[which.at.stage]] >  est.eff.bound.z[2 * (k-1) + 1]){
                reject.hyp[i, c(1,3)[which.at.stage]] = 1
                stage.decision[i, c(1,3)[which.at.stage]] = k
              }
            }
            
        # Sub-population 2
        # Look if already stopped at last stage
            index.used = which(index.sub.pop == 2 & index.stage == k)
            is.rejected = FALSE
            
        # Finding which treatments in sub-population one continued onto
        # stage k 
            which.at.stage = which(is.na(stage.decision[i, c(2,4)]))
            
        # If both continue onto second stage 
            if(length(which.at.stage) == 2){
              if(max(test.statistics[i, index.used]) > est.eff.bound[2 * (k-1) + 2]){
                reject.hyp[i, c(2,4)[which.max(test.statistics[i, index.used])]] = 1
                is.rejected = TRUE
                stage.decision[i,  c(2,4)[which.max(test.statistics[i, index.used])]] = k
              }
              
          # If the max test is rejected continue onto the next test
              if(is.rejected & min(test.statistics[i, index.used]) > est.eff.bound.z[2 * (k-1) + 2]){
                reject.hyp[i,  c(2,4)[which.min(test.statistics[i, index.used])]] = 1
                stage.decision[i,  c(2,4)[which.min(test.statistics[i, index.used])]] = k
              }
            }
            
        # If only one continues onto second stage 
            if(length(which.at.stage) == 1){
          # Finding the other one
              not.at.stage = setdiff(1:2, which.at.stage)
              
          # If the other one stopped for futility at earlier stages
              if(reject.hyp[i, c(2,4)[not.at.stage]] == 0)
                if(test.statistics[i, index.used[which.at.stage]] > est.eff.bound[2 * (k-1) + 2]){
                  reject.hyp[i, c(2,4)[which.at.stage]] = 1
                  stage.decision[i, c(2,4)[which.at.stage]] = k
                }
                
          # If the other one stopped for efficacy at earlier stages use the step down procedure
                if(reject.hyp[c(2,4)[not.at.stage]] == 1)
                  if(test.statistics[i, index.used[which.at.stage]] >  est.eff.bound.z[2 * (k-1) + 2]){
                    reject.hyp[i, c(2,4)[which.at.stage]] = 1
                    stage.decision[i, c(2,4)[which.at.stage]] = k
                  }
                }
                
              } # end k loop
            } else {
              k=1
            }
            
    # All H_0 not already stopped for efficacy or futility for futility at stage K
            stage.decision[i, which(is.na(stage.decision[i, ]))] = K
          }

  # Do the alpha reallocation

          for(i in 1:n.sim){
  # Cycle through stages
            for(k in 1:K){


  # Start with subpopulation 1
    # if both subpop 2 are rejected at or before stage k and at least one at k
              if(stage.decision[i, 2] <= k & stage.decision[i, 4] <= k & reject.hyp[i, 2] == 1 & reject.hyp[i, 4] == 1 & (stage.decision[i, 2] == k | stage.decision[i, 4] == k)){

                index.used = which(index.sub.pop == 1 & index.stage == K)
        # If both treatment arms are enrolled in subpopulation 1 at stage K
                if(stage.decision[i, 1] == K &  stage.decision[i, 3] == K){
            # Doing the maximum test
                  if(max(test.statistics[i, index.used]) > est.eff.bound.alpha[1]){
                # Reject the larger test statistic
                    reject.hyp[i, c(1,3)[which.max(test.statistics[i, index.used])]] = 1
                # Do the minimum test
                    if(min(test.statistics[i, index.used]) > est.eff.bound.alpha.z[1]){
                      reject.hyp[i, c(1,3)[which.min(test.statistics[i, index.used])]] = 1
                    }
                  }
                }

       # If only one treatment is enrolled at stage K
                if(sum(stage.decision[i, c(1,3)] == K) == 1){
                  which.stopped = which(stage.decision[i, c(1,3)] != K)
                  which.enrolled = which(stage.decision[i, c(1,3)] == K)
          # if the other treatment was stopped for efficacy
                  if(reject.hyp[i, c(1,3)[which.stopped]] ==  1){
                    index.used = which(index.sub.pop == 1 & index.stage == K)[which.enrolled]
              # Doing test for efficacy 
                    if(test.statistics[i, index.used] > est.eff.bound.alpha.z[1]){
                      reject.hyp[i, c(1,3)[which.enrolled]] = 1
                    }
                  }

          # if the other treatment was stopped for futility
                  if(reject.hyp[i, c(1,3)[which.stopped]] ==  0){
                    index.used = which(index.sub.pop == 1 & index.stage == K)[which.enrolled]
              # Doing test for efficacy 
                    if(test.statistics[i, index.used] > est.eff.bound.alpha[1]){
                      reject.hyp[i, c(1,3)[which.enrolled]] = 1
                    }
                  }
                }

              } # end if loop if both subpop 2 are rejected at or before stage k and at least one at k


  # Now with subpopulation 2
    # if both subpop 1 are rejected at or before stage k and at least one at k
              if(stage.decision[i, 1] <= k & stage.decision[i, 3] <= k & reject.hyp[i, 1] == 1 & reject.hyp[i, 3] == 1 & (stage.decision[i, 1] == k | stage.decision[i, 3] == k)){
                index.used = which(index.sub.pop == 2 & index.stage == K)
        # If both treatment arms are enrolled in subpopulation 1 at stage K
                if(stage.decision[i, 2] == K &  stage.decision[i, 4] == K){
            # Doing the maximum test
                  if(max(test.statistics[i, index.used]) > est.eff.bound.alpha[2]){
                # Reject the larger test statistic
                    reject.hyp[i, c(2,4)[which.max(test.statistics[i, index.used])]] = 1
                # Do the minimum test
                    if(min(test.statistics[i, index.used]) > est.eff.bound.alpha.z[2]){
                      reject.hyp[i, c(2,4)[which.min(test.statistics[i, index.used])]] = 1
                    }
                  }
                }

       # If only one treatment is enrolled at stage K
                if(sum(stage.decision[i, c(2,4)] == K) == 1){
                  which.stopped = which(stage.decision[i, c(2,4)] != K)
                  which.enrolled = which(stage.decision[i, c(2,4)] == K)
          # if the other treatment was stopped for efficacy
                  if(reject.hyp[i, c(2,4)[which.stopped]] ==  1){
                    index.used = which(index.sub.pop == 2 & index.stage == K)[which.enrolled]
              # Doing test for efficacy 
                    if(test.statistics[i, index.used] > est.eff.bound.alpha.z[2]){
                      reject.hyp[i, c(2,4)[which.enrolled]] = 1
                    }
                  }

          # if the other treatment was stopped for futility
                  if(reject.hyp[i, c(2,4)[which.stopped]] ==  0){
                    index.used = which(index.sub.pop == 2 & index.stage == K)[which.enrolled]
              # Doing test for efficacy 
                    if(test.statistics[i, index.used] > est.eff.bound.alpha[2]){
                      reject.hyp[i, c(2,4)[which.enrolled]] = 1
                    }
                  }
                }

              } # end if loop if both subpop 1 are rejected at or before stage k and at least one at k


            } # end for k loop
          } # end for i loop
          
          colnames(reject.hyp) = c("A1", "A2", "B1", "B2")
          colnames(stage.decision) = c("A1", "A2", "B1", "B2")

          return(list(rejection.matrix = reject.hyp,
            stage.decision = stage.decision,
            eff.boundaries = efficacy.boundary))
        }


# Test code
test.code = function(){
# Creating an example
          alpha.allocation = rep(0.025/4, 4)
          time = c(1,2)
          accrural.rate = 20
          delay = 0
          outcome.type = "continuous"
          prop.pop.1 = 0.5
          var.vec.pop.1 = c(1,1,1)
          var.vec.pop.2 = c(1,1,1)
          mean.sub.pop.1 = c(0.1, 0.1, 0.1)
          mean.sub.pop.2 = c(0.1, 0.1, 0.1)
          futility.boundary = c(-Inf, -Inf, -Inf, -Inf, Inf, Inf, Inf, Inf)
          n.sim = 10^4

# Create input into the 
          analytic.n.per.stage = get.sample.sizes.per.stage(n.arms = 3, n.per.arm = 100, subpopulation.sizes = c(0.5, 0.5),
           interim.info.times = c(0.3,1), accrual.rate = 20, delay = 0)$analytic.n.per.stage

# Create joint distribution
          temp = construct.test.statistics.joint.distribution(analytic.n.per.stage,
            mean.sub.pop.1=mean.sub.pop.1,
            mean.sub.pop.2=mean.sub.pop.2,
            var.vec.pop.1=var.vec.pop.1,
            var.vec.pop.2=var.vec.pop.2,
            outcome.type = outcome.type,
            prop.pop.1 = prop.pop.1,
            max.follow = NULL,
            enrollment.period = NULL,
            hazard.rate.pop.1 = NULL,
            hazard.rate.pop.2 = NULL,
            time = NULL)

          non.centrality.parameter.vec = temp$non.centrality.parameter.vec
          cov.mat.used = temp$cov.mat.used
          information.vector = temp$information.vector
          
          n.mc = 10^6
          test.statistics <- rmvnorm(n.mc, mean=non.centrality.parameter.vec, sigma=cov.mat.used)
          eff.bound.used = get.eff.bound(alpha.allocation, cov.mat.used, err.tol = 10^-3)

          eval.design = design.evaluate(test.statistics,                                  
            efficacy.boundary = eff.bound.used,
            futility.boundary,
            alpha.allocation)

          which.false = c(1,2,3,4)
          mean(apply(eval.design$rejection.matrix, 1, function(x){(sum(x[which.false] == 1) > 0)}))

          tianchen = FALSE
          if(tianchen == TRUE){
            true.estimand.value = c(mean.sub.pop.1[2] - mean.sub.pop.1[1], mean.sub.pop.1[3] - mean.sub.pop.1[1], mean.sub.pop.2[2] - mean.sub.pop.2[1], mean.sub.pop.2[3] - mean.sub.pop.2[1])
            stage.decision = eval.design$stage.decision
            information.level = temp$information.vector

            compute.performance(test.statistics, cov.mat.used, stage.decision, information.level, true.estimand.value)
            compute.performance(test.statistics, cov.mat.used, stage.decision, information.level, true.estimand.value, ci.level = 0.80)
          }

# Survival outcome
          max.follow = 2
          hazard.rate.pop.1 = c(1,1,1)
          hazard.rate.pop.2 = c(1,1,1)
          time = c(1,2)
          enrollment.period = 3 

          accrural.rate = 100
          prop.pop.1 = 0.5
          outcome.type = "survival"
          futility.boundary = c(-Inf, -Inf, -Inf, -Inf, Inf, Inf, Inf, Inf)
          n.sim = 10^4
          alpha.allocation = rep(0.025/4, 4)
          censoring.rate = 0.1


# Create input into the 
          analytic.n.per.stage = get.sample.sizes.per.stage(n.arms = 3, n.per.arm = 100, subpopulation.sizes = c(0.5, 0.5),
           interim.info.times = c(0.3,1), accrual.rate = 20, delay = 0)$analytic.n.per.stage

# Create joint distribution
          temp = construct.test.statistics.joint.distribution(analytic.n.per.stage,
            mean.sub.pop.1=NULL,
            mean.sub.pop.2=NULL,
            var.vec.pop.1=NULL,
            var.vec.pop.2=NULL,
            outcome.type = outcome.type,
            prop.pop.1 = prop.pop.1,
            max.follow = max.follow,
            enrollment.period = enrollment.period,
            hazard.rate.pop.1 = hazard.rate.pop.1,
            hazard.rate.pop.2 = hazard.rate.pop.2,
            time = time,
            relative.efficiency = 1.1,
            censoring.rate = censoring.rate)
          
          non.centrality.parameter.vec = temp$non.centrality.parameter.vec
          cov.mat.used = temp$cov.mat.used

          n.mc = 10^6
          test.statistics <- rmvnorm(n.mc, mean=non.centrality.parameter.vec, sigma=cov.mat.used)
          eff.bound.used = get.eff.bound(alpha.allocation, cov.mat.used, err.tol = 10^-3)

          eval.design = design.evaluate(test.statistics,                                  
            efficacy.boundary = eff.bound.used,
            futility.boundary,
            alpha.allocation)

          which.false = c(1,2,3,4)
          mean(apply(eval.design$rejection.matrix, 1, function(x){(sum(x[which.false] == 1) > 0)}))
        }

test.smart.av = function(){

# Creating an example
alpha.allocation = c(0.239, 0.251, 0.255, 0.255) * 0.025
time = c(0.409,1)
accrual.rate = 663
delay = 0
outcome.type = "continuous"
prop.pop.1 = 0.49
var.vec.pop.1 = c(60^2,60^2,60^2)
var.vec.pop.2 = c(60^2,60^2,60^2)
mean.sub.pop.1 = c(0, 15, 15)
mean.sub.pop.2 = c(0, 15, 15)
futility.boundary = c(-2.40,-5.69,0.71,0.83, Inf, Inf, Inf, Inf)
delay = 0.5
n.sim = 10^4

# Create input into the 

# Create input into the 
analytic.n.per.stage = get.sample.sizes.per.stage(n.arms = 3, n.per.arm = accrual.rate * max(time)/3, subpopulation.sizes = c(prop.pop.1, 1- prop.pop.1),
           interim.info.times = time/max(time), accrual.rate = accrual.rate, delay = delay)$analytic.n.per.stage

# Create joint distribution
          temp = construct.test.statistics.joint.distribution(analytic.n.per.stage,
            mean.sub.pop.1=mean.sub.pop.1,
            mean.sub.pop.2=mean.sub.pop.2,
            var.vec.pop.1=var.vec.pop.1,
            var.vec.pop.2=var.vec.pop.2,
            outcome.type = outcome.type,
            prop.pop.1 = prop.pop.1,
            max.follow = NULL,
            enrollment.period = NULL,
            hazard.rate.pop.1 = NULL,
            hazard.rate.pop.2 = NULL,
            time = NULL)

          non.centrality.parameter.vec = temp$non.centrality.parameter.vec
          cov.mat.used = temp$cov.mat.used
          
          n.mc = 10^6
          test.statistics <- rmvnorm(n.mc, mean=non.centrality.parameter.vec, sigma=cov.mat.used)
          eff.bound.used = get.eff.bound(alpha.allocation, cov.mat.used, err.tol = 10^-3)

          eval.design = design.evaluate(test.statistics,                                  
            efficacy.boundary = eff.bound.used,
            futility.boundary,
            alpha.allocation)

          which.false = c(1,2,3,4)
          mean(apply(eval.design$rejection.matrix, 1, function(x){(sum(x[which.false] == 1) > 0)}))
 
}
