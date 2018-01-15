# Note to self: 2017-08-15 19:52:
# We can use names(formals(transform.parameters)) to get the set of formal
# arguments to a function: this might be a useful way of making sure all
# required parameters are passed to avoid downstream errors.


### Wrappers to User-Defined Functions for the 
### User's Class of Designs to Be Optimized Over
#   Must contain the following 3 functions described below: 
#   construct.test.statistics.joint.distribution
#   get.eff.bound
#   design.evaluate
source("Backend2Populations2Arms.R")
# Computes distribution of test statistics in a given scenario, using canonical joint distribution
get.z.distribution <- 
  function(...){
    construct.test.statistics.joint.distribution(...)
  }
# Computes efficacy stopping boundaries
get.efficacy.dunnett <-
  function(...){
    get.eff.bound(...)
  }
# Evaluates performance of simulated trials 
evaluate.design.dunnett <- 
  function(...){
    design.evaluate(...)
  }

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
squash <- 
  function(x,
           x.min=0,
           x.max=1,
           method=c("linear", "logit") [1]) {
  stopifnot(is.numeric(x), 
            x.max>x.min)
  switch(method,
         linear=pmin(pmax(x.min, x), x.max), # Linear Interpolation 
         logit=x.min + (x.max-x.min)*plogis(x, location=(x.max+x.min)/2,
                                            scale=(x.max-x.min)/6)
  )
}

### get.sample.sizes.per.stage #################################################
# Description: Calculate the number of patients per arm in multi-stage clinical
#   trials with one or more treatments compared with controls in one population,
#   potentially divided into disjoint subgroups. Single subpopulations and 
#   single-stage trials are also supported.
#
# Input:
#   n.arms (numeric scalar - whole): number of total arms (L+1)
#   n.per.arm (numeric scalar - whole: number of patients per arm
#   subpopulation.sizes (numeric vector - vector of probabilities): the 
#     proportion of the population in each disjoint subgroup. Must sum to 1.
#     Using subpopulation.sizes = NULL or 1 gives a trial in a single
#     population.
#   interim.info.times (numeric vector - sorted proportions, ending in 1): The 
#     proportion of patients enrolled by each interim analysis. If no interim
#     analyses are desired, set to NULL or 1. 
#   accrual.rate (numeric scalar): number of patients enrolled per year
#   delay - numeric : time in years from randomization to observation of primary
#     outcome
#   warn.last.accrual.before.interim (logical scalar): warn if an interim 
#     analysis occurs after enrollment has stopped: such an interim analysis can
#     not reduce sample size, but may reduce costs of follow-up visits.
#
# Output: 
#   analytic.n.per.stage ([K x J(L+1)] matrix), where K is the number of stages,
#     J is the number of subpopulations, and L is the number of treatments
#     compared to the control arm. For continuous or binary outcomes, this 
#     represents the patients with primary outcome at each interim analysis; for
#     survival outcome, this represents number of patients enrolled. Each row is
#     one stage of the study, with each column indicating the number of outcomes
#     available in that (treatment x subpopulation) stratum. Columns are
#     arranged by subpopulation (S1, ..., SJ) within treatments (T0, ..., TL):
#
#     stage 1: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#     stage 2: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#       ...
#     stage k: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#
#   total.n.per.stage ([K x J(L+1)] matrix): all enrolled patients at each
#     interim analysis, irrespective of when the primary outcome is observed.
#     stage 1: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#     stage 2: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#       ...
#     stage k: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#
#   accrual.years (K x 1 numeric vector): time to accrue the patients in each
#     stage of a trial, assuming no interim stopping.
#
#   outcome.years (K x 1 numeric vector): time to observe the outcomes of the
#     patients accrued in each of K stages of a trial: may exceed the length
#     of enrollment, depending on enrollment rate and delay. Again, this assumes
#     no interim stopping.
#
# NOTE: When delay = 0, total.n.per.stage and analytic.n.per.stage will be
#   equal when sample sizes are integers: when not integers, the number of
#   observed outcomes is rounded down using floor() and the number of enrolled
#   participants is rounded up using ceiling() to be conservative. For survival
#   outcomes, delay should be set to 0.
get.sample.sizes.per.stage <-
  function(n.arms, n.per.arm, subpopulation.sizes,
           interim.info.times, accrual.rate, delay,
           warn.last.accrual.before.interim=TRUE, restrict.enrollment = FALSE){

    max.arms <- 8 # This can be increased, but not beyond 26 (=dim(letters))

    interim.info.times = sort(interim.info.times)
    
    # Input Checks
    stopifnot(is.whole.number(n.arms) & n.arms > 0 & n.arms <= max.arms,
              is.finite(n.arms) & is.scalar(n.arms),
              is.whole.number(n.per.arm) & n.per.arm > 0,
              is.finite(n.per.arm) & is.scalar(n.per.arm),
              is.numeric(accrual.rate) & accrual.rate > 0,
              is.finite(accrual.rate) & is.scalar(accrual.rate),
              is.numeric(delay) & delay >=0,
              is.finite(delay) & is.scalar(delay),
              is.logical(warn.last.accrual.before.interim) &
                is.scalar(warn.last.accrual.before.interim))
    
    # Check interim.info.times and subpopulation.sizes for coherence
    if(!is.null(subpopulation.sizes)){
      stopifnot(is.numeric(subpopulation.sizes),
                length(subpopulation.sizes) >= 1,
                all.equal(sum(subpopulation.sizes), 1))
      if(length(subpopulation.sizes)>1) {
        stopifnot(in.open.interval(subpopulation.sizes, lb=0, ub=1))
      }
    } else {
      subpopulation.sizes <- 1
    }
    
    if(!is.null(interim.info.times)){
      # Numeric, Not empty, strictly increasing
#      stopifnot(is.numeric(interim.info.times),
#                length(interim.info.times) >= 1, 
#                !is.unsorted(interim.info.times, strictly=TRUE))
      if(length(interim.info.times)>=1) {
        # Positive, less than or equal to 1, last value == 1
        stopifnot(sum(interim.info.times<0)==0,
                  mean(interim.info.times<=1)==1)
        if(tail(interim.info.times, 1)!=1) interim.info.times <- 
            c(interim.info.times, 1)
      }
    } else {
      interim.info.times <- 1
    }
    
    # Enforce n.per.arm limit
    if(length(subpopulation.sizes)>1){
      max.per.subpop <- round(head(subpopulation.sizes, -1)*n.per.arm)
      max.per.subpop <- c(max.per.subpop, n.per.arm-sum(max.per.subpop))
    } else max.per.subpop <- n.per.arm
    
    # See if a subpopulation is empty or near empty
    if(sum(round(subpopulation.sizes*n.per.arm)<1)){
      warning(paste0("Empty subgroup: Check n.per.arm (", n.per.arm,
                     ") and subpopulation.sizes (", 
                     paste(subpopulation.sizes, collapse=", "), ")."))
    }
    
    # Calculate time to accrue n.per.arm
    stages <- length(interim.info.times) + 1*(tail(interim.info.times, 1)!=1)
    n.subgroups <- length(subpopulation.sizes)
    accrual.years <- n.arms*interim.info.times*n.per.arm/accrual.rate
    outcome.years <- accrual.years + delay
    
    # Calculate total patients accrued at interim analyses
    accrued.patients <-
      floor(kronecker(outcome.years*accrual.rate/n.arms,
                        t(subpopulation.sizes)))
    
    # Cap max enrollment - can overshoot with long delay
    for(i in 1:nrow(accrued.patients)) {
      if(sum(accrued.patients[i,] > max.per.subpop)>0) {
        accrued.patients[i,] <- max.per.subpop
      }
    }
    
    # Count the number of outcomes observed at each information time
    accrued.outcomes <-  
      floor(kronecker(interim.info.times,
                      t(subpopulation.sizes))*n.per.arm)
    
    # Enforce n.per.arm limit
    accrued.outcomes[nrow(accrued.outcomes),] <- tail(accrued.patients, 1)
    
    # Apply names to each arm, subgroup; C for control is column 1 in each block
    arm.names <- c(LETTERS[3], LETTERS[1:max.arms][-3])[1:n.arms]
    analytic.n.per.stage <- kronecker(t(rep(1, n.arms)), accrued.outcomes)
    
    colnames(analytic.n.per.stage) <-
      paste0(rep(arm.names, each=n.subgroups), rep(1:n.subgroups, n.arms))
    
    total.n.per.stage <- kronecker(t(rep(1, n.arms)), accrued.patients)
    colnames(total.n.per.stage) <-
      paste0(rep(arm.names, each=n.subgroups), rep(1:n.subgroups, n.arms))

    # Optionally warn if accrual completes before interim analysis - 
    # If all patients are enrolled before an interim analysis, sample size
    # can't be reduced by efficacy/futility stopping.
    if(warn.last.accrual.before.interim &
       sum(duplicated(total.n.per.stage)>0)) {
      warning("Accrual completed before interim analyses.")
    }
    
    if(sum(analytic.n.per.stage==0)>0){
      warning("Some stage did not accrue a patient in one subpopulation.")
    }

   if(restrict.enrollment == TRUE){
      for(k in 1:n.arms){
          analytic.n.per.stage[, 2*k] = (analytic.n.per.stage[, 2*k] - analytic.n.per.stage[1, 2*k])
          total.n.per.stage[, 2*k] = (total.n.per.stage[, 2*k] - total.n.per.stage[1, 2*k])
      }
    }
    
    return(list(analytic.n.per.stage=analytic.n.per.stage,
                total.n.per.stage=total.n.per.stage,
                accrual.years=accrual.years,
                outcome.years=outcome.years))
  }

### get.sample.sizes.per.stage example #########################################
# get.sample.sizes.per.stage(n.arms=1, n.per.arm=100,
#                            subpopulation.sizes=NULL, interim.info.times=NULL,
#                            accrual.rate=100, delay=100)
# get.sample.sizes.per.stage(n.arms=1, n.per.arm=100,
#                            subpopulation.sizes=c(0.1, 0.9),
#                            interim.info.times=NULL,
#                            accrual.rate=100, delay=100)
# get.sample.sizes.per.stage(n.arms=2, n.per.arm=100,
#                            subpopulation.sizes=c(0.1, 0.9),
#                            interim.info.times=NULL,
#                            accrual.rate=100, delay=100)
# get.sample.sizes.per.stage(n.arms=3, n.per.arm=500,
#                            subpopulation.sizes=c(0.1, 0.9),
#                            interim.info.times=c(.25, .5, .75, 1),
#                            accrual.rate=100, delay=0.1)


### list.potential.trials ######################################################
# Description: Create a look-up table that gives the sample size and duration
# of each potential outcome of a clinical trial. This can be used to calculate
# ESS and duration from a matrix of simulated trial results. If each row is a
# simulated trial, and each column specifies the stage at which each (arm x 
# subgroup) strata ended the trial, simply tabulate the number of outcomes of
# each type, merge with the look-up table, and then use weighted.mean() on
# sample size and duration.

list.potential.trials <-
  function(total.n.per.stage,
           accrual.rate,
           subpopulation.sizes){
    
    # Get number of arms, subgroups, and stages from total.n.per.stage
    n.arms <- 
      length(unique(substring(colnames(total.n.per.stage), first=1, last=1)))
    n.subgroups <-
      length(unique(substring(colnames(total.n.per.stage), first=2)))
    stages <- nrow(total.n.per.stage)
    
    # Get the increment in sample size at each stage
    incremental.n <- diff(rbind(0, total.n.per.stage))
    
    # Create a list of columns for each subgroup;
    # Vector of columns for each control arm
    subgroup.col <- list()
    for(j in 1:n.subgroups) subgroup.col[[j]] <-
      which(substring(colnames(total.n.per.stage), first=2)==j)
    control.col <-
      match(paste0("C", 1:n.subgroups), colnames(total.n.per.stage))
    
    ### Calculate possible sample sizes, durations
    # Create all potential outcomes: Arm x Subopulation x Stage combinations
    all.trials <- matrix(1:stages, nrow=stages, ncol=n.subgroups*n.arms)
    colnames(all.trials) <- colnames(total.n.per.stage)
    all.trials <- expand.grid(data.frame(all.trials))
    
    if(nrow(all.trials)>1){
      
      # Restrict such that active treatments can't enroll without controls:
      # Keep only trials where the last stage reached by any treatment is equal
      # to the last stage reached by control
      for(j in 1:n.subgroups){
        tx.cols <- setdiff(subgroup.col[[j]], control.col[j])
        if(length(tx.cols)>0){
          all.trials <-
            all.trials[which(apply(as.matrix(all.trials[,tx.cols]), 1, max) == 
                               all.trials[,control.col[j]]),]
        }
      }
      
      stage.duration <- matrix(0, nrow=nrow(all.trials), ncol=n.subgroups)
      trial.duration <- matrix(0, nrow=nrow(all.trials))
      
      for (k in 1:stages){
        patients.stage.k <- 
          (all.trials >= k)* # If Treatment x Subpopulation strata in stage k
          kronecker(matrix(1, nrow=nrow(all.trials)), 
                    t(incremental.n[k,])) # Incremental sample size by strata
        
        # Count arms enrolling; Get accrual rate in subpopulation by adjusting
        # subpopulation accrual rate by arms enrolling in current stage
        for (j in 1:n.subgroups){
          arms.enrolling <- 
            rowSums(as.matrix((all.trials >= k)[, subgroup.col[[j]]]))
          randomization.rate <- 
            (accrual.rate)*subpopulation.sizes[j]*n.arms/arms.enrolling
          stage.duration[,j] <- 
            rowSums(as.matrix(patients.stage.k[, subgroup.col[[j]]]))/
            randomization.rate
        }
        trial.duration <- trial.duration + apply(stage.duration, 1, max)
      }
      
      # Calculate sample size of each trial
      potential.trials <-
        data.frame(all.trials,
                   sample.size=apply(t(all.trials), 2, function(x) 
                     sum(total.n.per.stage[cbind(x, 1:ncol(all.trials))])),
                   duration=trial.duration)

    } else {
      potential.trials <-
        data.frame(all.trials,
                   sample.size=sum(total.n.per.stage),
                   duration=sum(total.n.per.stage)/accrual.rate)
    }
    return(potential.trials)
}

per.stage.sample.sizes <-
  get.sample.sizes.per.stage(n.arms=3, n.per.arm=500,
                             subpopulation.sizes=c(0.1, 0.9),
                             interim.info.times=c(.5, .75, 1),
                             accrual.rate=100, delay=0.1)
# NOTE: this does not match because total.n.per.stage is rounded, resulting in
# over-estimation of duration.
list.potential.trials(total.n.per.stage=
                        per.stage.sample.sizes$total.n.per.stage,
                      accrual.rate=100,
                      subpopulation.sizes=c(0.1, 0.9))




### get.trial.criteria #########################################################
# Description: determine the maximum sample size, expected sample size, expected
#   duration, power, type I error, and familywise type I error rate (FWER) for
#   a clinical trial based on a trial.simulation object, a matrix of hypothesis
#   rejections (decisions.rejections), a matrix of stages at which decisions
#   were made (decisions.stages), the columns corresponding to null hypotheses
#   (null.hypotheses), a look-up table of potential trials with their sample
#   sizes and durations calculated (potential.trials)
#
#
#   if no hypotheses are null, specify null.hypotheses=NULL
get.trial.criteria <-
  function(decisions.rejections,
           decisions.stages,
           null.hypotheses=NULL,
           potential.trials){
    
    # Add warning if names of decisions.stages don't match potential.trials
    
    # Get number of arms, subgroups, and stages from total.n.per.stage
    n.arms <- 
      length(unique(substring(colnames(decisions.stages), first=1, last=1)))
    n.subgroups <- length(unique(substring(colnames(decisions.stages), first=2)))
    n.simulations <- nrow(decisions.stages)
    
    # Create a list of columns for each subgroup;
    # Vector of columns for each control arm
    subgroup.col <- list()
    for(j in 1:n.subgroups) subgroup.col[[j]] <-
      which(substring(colnames(decisions.stages), first=2)==j)
    
    conj.power.n.subgroups <-
      matrix(NA, nrow=n.simulations, ncol=n.subgroups)
    disj.power.n.subgroups <-
      matrix(NA, nrow=n.simulations, ncol=n.subgroups)
    
    # Copy decisions.rejections
    empirical.alpha <- empirical.power <- decisions.rejections
    if(length(null.hypotheses) > 0){
      empirical.power[, null.hypotheses] <- NA
      empirical.alpha[, -null.hypotheses] <- NA
    } else {
      empirical.alpha <- empirical.alpha*NA
    }
    
    for(j in 1:n.subgroups){
      non.null.sub.cols <- setdiff(subgroup.col[[j]], null.hypotheses)
      if(length(non.null.sub.cols)>0){
        sub.cols <- as.matrix(empirical.power[,non.null.sub.cols])
        conj.power.n.subgroups[,j] <- mean(rowMeans(sub.cols)==1)
        disj.power.n.subgroups[,j] <- mean(rowMeans(sub.cols)>0)
      }
    }
    
    conj.power.n.subgroups <- colMeans(conj.power.n.subgroups)
    disj.power.n.subgroups <- colMeans(disj.power.n.subgroups)
    names(disj.power.n.subgroups) <- paste0("DIS_Power_", 1:n.subgroups)
    names(conj.power.n.subgroups) <- paste0("CON_Power_", 1:n.subgroups)

    conj.power <- mean(rowMeans(empirical.power, na.rm=TRUE)==1)
    disj.power <- mean(rowMeans(empirical.power, na.rm=TRUE)>0)
    empirical.power <- colMeans(empirical.power==1)
    
    if(length(null.hypotheses)>0){
      fwer <- mean(rowMeans(empirical.alpha, na.rm=TRUE)>0)
      type.1.error <- colMeans(empirical.alpha==1)
    } else {
      fwer <- NA
      type.1.error <- head(decisions.stages, 1)*NA
    }
    
    distribution.of.trials <- count(decisions.stages)
    names(distribution.of.trials) <- c(colnames(decisions.stages), "frequency")
    distribution.of.trials$proportion <- 
      distribution.of.trials$frequency/n.simulations
    
    distribution.of.trials <- merge(potential.trials,
                                    distribution.of.trials, all=TRUE)
    distribution.of.trials <-
      distribution.of.trials[union(names(potential.trials),
                                   names(distribution.of.trials))]
    
    return(list(distribution.of.trials=distribution.of.trials,
                empirical.power=empirical.power,
                conj.power=conj.power,
                disj.power=disj.power,
                conj.power.n.subgroups=conj.power.n.subgroups,
                disj.power.n.subgroups=disj.power.n.subgroups,
                type.1.error=type.1.error,
                fwer=fwer))
  }


### linear.threshold.penalty ###################################################
# Description: determine if a matrix of values meets or exceeds a matrix of
#   thresholds. A linear penalty is applied to values below the thresholds.
#   A loss matrix is returned of the same dimension as the power matrix and 
#   thresholds.
#
# If any scenario x hypothesis combinations do not contribute to the penalty,
#   they should be replaced with NA values before using power.penalty().
#
# Input:
#   x (numeric matrix): a (RxC) matrix of values
#   threshold.matrix (numeric matrix): a (RxC) matrix of thresholds
#   linear.penalty (numeric scalar - positive): a positive linear penalty to
#     scale the absolute difference between x and threshold.matrix
#
# Output: 
#   loss (numeric matrix): a (MxH) matrix of loss values.
linear.threshold.penalty <-
  function(x,
           threshold.matrix,
           linear.penalty) {
    # Check to see if x are probability scale and equal dimensions to
    # power thresholds, and thresholds are feasible.
    stopifnot(is.numeric(linear.penalty) & is.scalar(linear.penalty),
              linear.penalty>=0,
              identical(dim(x), dim(threshold.matrix)))
    
    return(linear.penalty*(x<threshold.matrix)*abs(x-threshold.matrix))
  }

### dunnett.means ##############################################################
# Description: For continuous or binary outcomes, evaluate trial design performance,
#   simulate trials under each user-defined scenario and return trial 
#   i.e., performance including expected sample size, power, Type I error, and duration.
#   Intermediate steps include generating the sample size per stage; then, for each scenario:                           
#   the multivariate normal distribution of test statistics is constructed by a call to
#   construct.test.statistics.joint.distribution, the stopping times and null hypotheses                           
#   rejected are computed for each simulated trial by a call to design.evaluate
#   and then performance metrics are computed. 
#                   
# Input:
#   n.arms (numeric scalar - whole): number of total arms (currently 
#     restricted to 2)
#   n.per.arm (numeric scalar - whole: number of patients per arm
#   accrual.rate (positive numeric scalar): number of patients per year accrued
#     in the combined population
#   delay (numeric scalar): time in years from randomization to observation of
#     primary outcome
#   subpopulation.sizes (numeric vector - vector of probabilities): the 
#     proportion of the population in each disjoint subgroup. Must sum to 1.
#     Using subpopulation.sizes = NULL or 1 gives a trial in a single
#     population.
#   interim.info.times (numeric vector - sorted proportions, ending in 1): The 
#     proportion of patients enrolled by each interim analysis. If no interim
#     analyses are desired, set to NULL or 1. (currently restricted to length 
#     up to 2)
#   total.alpha (numeric scalar): total allowable familywise type I error rate.
#     Defaults to 0.025.
#   alpha.allocation (numeric vector): proportion of error rate spent on each
#     analysis. Must sum to 1.
#   outcome.type (string): outcome distribution: either "continuous" or "binary"
#   outcome.mean (M x JT numeric matrix): mean outcomes by scenario.
#   outcome.sd (M x JT numeric matrix): mean SDs by scenario.
#   mcid (numeric scalar): the minimum clinically important difference - on the
#     same scale as outcome.mean
#   futility.boundaries (numeric vector): futility boundaries
#   n.simulations (numeric scalar): Number of simulated trials to estimate trial
#     performance. Defaults to 10000.

dunnett.means <- 
  function(n.arms,
           n.per.arm,
           accrual.rate,
           delay=0,
           subpopulation.sizes,
           interim.info.times,
           total.alpha=0.025,
           alpha.allocation,
           outcome.type,
           outcome.mean,
           outcome.sd=NULL,
           mcid=0,
           futility.boundaries=NULL,
           n.simulations=50000){
    
    n.subgroups <- length(subpopulation.sizes)
    n.stages <- length(interim.info.times) + 1*(tail(interim.info.times, 1)!=1)
    n.scenarios <- nrow(outcome.mean)
    
    if(outcome.type=="binary"){
      outcome.sd <- outcome.mean*(1-outcome.mean)
    }
    
    efficacy.boundary.function <- get.efficacy.dunnett
    efficacy.boundary.parameters <- list(err.tol=10^-3)
    test.procedure.function <- evaluate.design.dunnett
    
    per.stage.sample.sizes <-
      get.sample.sizes.per.stage(n.arms=n.arms, 
                                 n.per.arm=n.per.arm,
                                 subpopulation.sizes=subpopulation.sizes,
                                 interim.info.times=interim.info.times,
                                 accrual.rate=accrual.rate,
                                 delay=delay)
    
    analytic.n.per.stage <- per.stage.sample.sizes$analytic.n.per.stage
    total.n.per.stage <- per.stage.sample.sizes$total.n.per.stage
    
    alpha.allocation <- alpha.allocation*total.alpha
    trial.criteria <- trial.criteria.by.scenario <- list()
    
    # List potential trial outcomes
    potential.trials <-
      list.potential.trials(total.n.per.stage=total.n.per.stage,
                            accrual.rate=accrual.rate,
                            subpopulation.sizes=subpopulation.sizes)
    
    # Determine which hypotheses are null, which are above MCID
    mean.differences <-
      with(data.frame(outcome.mean), data.frame(A1=A1-C1, A2=A2-C2, B1=B1-C1, B2=B2-C2)) # JON CHANGED
    null.hypotheses <- 1*(abs(mean.differences)<1e-6)
    above.mcid <- 1*(abs(mean.differences>=mcid))
    
    # In case of MCID=0:
    above.mcid[which(null.hypotheses==1, arr.ind=TRUE)] <- 0
    
    
    for(m in 1:n.scenarios){ 
      ### SPECIFIC TO 2-ARM and 3-ARM DESIGNS ###
      prop.pop.1 <- subpopulation.sizes[1]
      
      if(n.arms==2){
        mean.sub.pop.1 <- outcome.mean[m, c("C1", "A1")]
        mean.sub.pop.2 <- outcome.mean[m, c("C2", "A2")]
        var.vec.pop.1 <- (outcome.sd[m, c("C1", "A1")])^2
        var.vec.pop.2 <- (outcome.sd[m, c("C2", "A2")])^2
      } else if(n.arms==3){
        mean.sub.pop.1 <- outcome.mean[m, c("C1", "A1", "B1")]
        mean.sub.pop.2 <- outcome.mean[m, c("C2", "A2", "B2")]
        var.vec.pop.1 <- (outcome.sd[m, c("C1", "A1", "B1")])^2
        var.vec.pop.2 <- (outcome.sd[m, c("C2", "A2", "B2")])^2
      } else {stop()}
      
      z.parameters <-
        get.z.distribution(analytic.n.per.stage,
                           mean.sub.pop.1=mean.sub.pop.1,
                           mean.sub.pop.2=mean.sub.pop.2,
                           var.vec.pop.1=var.vec.pop.1,
                           var.vec.pop.2=var.vec.pop.2,
                           outcome.type=outcome.type,
                           prop.pop.1=prop.pop.1)
      
      # Get efficacy boundaries
      efficacy.boundaries <-
        do.call(what=efficacy.boundary.function,
                args=c(list(cov.mat.used=z.parameters$cov.mat.used,
                            alpha.allocation=alpha.allocation),
                       efficacy.boundary.parameters))
      
      # Sample test statistics 
      test.statistics <- 
        with(z.parameters,
             rmvnorm(n=n.simulations,
                     mean=non.centrality.parameter.vec,
                     sigma=cov.mat.used))
      
      # Determine which hypotheses are rejected at each stage
      test.decisions <-
        do.call(what=test.procedure.function,
                args=list(test.statistics=test.statistics,
                          efficacy.boundary=efficacy.boundaries,
                          ### SPECIFIC TO 2-ARM AND 3-ARM DESIGNS ###
                          futility.boundary=c(futility.boundaries, rep(-Inf,4)),
                          alpha.allocation=alpha.allocation))
      
      # Get ESS, Duration, Power, Type I Error, FWER
      trial.criteria.by.scenario[[m]] <-
        get.trial.criteria(decisions.rejections=test.decisions$rejection.matrix,
                           decisions.stages=test.decisions$stage.decision,
                           null.hypotheses=which(null.hypotheses[m,]==1),
                           potential.trials=potential.trials)
      # Label distribution.of.trials by scenario
      trial.criteria.by.scenario[[m]]$distribution.of.trials <-
        data.frame(scenario=m,
                   trial.criteria.by.scenario[[m]]$distribution.of.trials)
    }
    
    # Unpack list of trial criteria - Label with scenario number
    criteria.names <- unique(unlist(lapply(trial.criteria.by.scenario, names)))
    for(i in 1:length(criteria.names)){
      criteria.result <- 
        do.call(what=rbind,
                args=lapply(trial.criteria.by.scenario, function(x)
                  get(criteria.names[i], x)))
      if(is.null(colnames(criteria.result))){
        colnames(criteria.result) <- criteria.names[i]
      }
      trial.criteria[[i]] <- criteria.result
    }
    
    names(trial.criteria) <- criteria.names
    
    trial.criteria$ess <-
      aggregate(sample.size*proportion~scenario,
                FUN=sum,
                data=trial.criteria$distribution.of.trials)[,2]
    
    trial.criteria$duration <-
      aggregate(duration*proportion~scenario,
                FUN=sum,
                data=trial.criteria$distribution.of.trials)[,2]
    
    trial.criteria$null.hypotheses <- null.hypotheses
    trial.criteria$above.mcid <- above.mcid
    
    return(trial.criteria)
  }


### dunnett.survival ###########################################################
# Description: For survival outcomes, evaluate trial design performance, i.e., 
#   simulate trials under each user-defined scenario and return trial 
#   performance including expected sample size, power, Type I error, and duration.
#   Intermediate steps include generating the sample size per stage; then, for each scenario:                           
#   the multivariate normal distribution of test statistics is constructed by a call to
#   construct.test.statistics.joint.distribution, the stopping times and null hypotheses                           
#   rejected are computed for each simulated trial by a call to design.evaluate
#   and then performance metrics are computed.                            
#
# Input:
#   n.arms (numeric scalar - whole): number of total arms (currently 
#     restricted to 2)
#   accrual.rate (positive numeric scalar): number of patients per year accrued
#     in the combined population
#   subpopulation.sizes (numeric vector - vector of probabilities): the 
#     proportion of the population in each disjoint subgroup. Must sum to 1.
#     Using subpopulation.sizes = NULL or 1 gives a trial in a single
#     population.
#   hazard.rate (M x JT numeric matrix): hazard rates by scenario.
#   total.alpha (numeric scalar): total allowable familywise type I error rate.
#     Defaults to 0.025.
#   alpha.allocation (numeric vector): proportion of error rate spent on each
#     analysis. Must sum to 1.
#   time (numeric vector - sorted): The years at which analyses occur.
#   max.follow (numeric scalar):
#   censoring rate (numeric scalar):
#   ni.margin (numeric scalar):
#   non.inferiority (logical scalar):
#   mcid (numeric scalar): the minimum clinically important difference - on the
#     same scale as outcome.mean
#   futility.boundaries (numeric vector): futility boundaries
#   n.simulations (numeric scalar): Number of simulated trials to estimate trial
#     performance. Defaults to 10000.
dunnett.survival <-
  function(n.arms,
           accrual.rate,
           enrollment.period,
           subpopulation.sizes,
           hazard.rate,
           total.alpha=0.025,
           alpha.allocation,
           time,
           max.follow,
           censoring.rate,
           futility.boundaries,
           ni.margin=NULL,
           non.inferiority=FALSE,
           mcid=0,
           restrict.enrollment=FALSE,
           boundary.to.enroll,
           n.simulations=1000){
    
    n.subgroups <- length(subpopulation.sizes)
    n.stages <- length(time)
    
    hazard.ratios <-
      with(data.frame(hazard.rate), 
           data.frame(A1=A1/C1, A2=A2/C2))
    
    prop.pop.1 <- subpopulation.sizes[1]
    n.scenarios <- nrow(hazard.rate)
    
    n.per.arm <- ceiling((accrual.rate*enrollment.period)/n.arms)
    
    # Functions for Dunnett Enrollment Modification & Efficacy Boundaries
    test.procedure.function=evaluate.design.dunnett
    efficacy.boundary.function=get.efficacy.dunnett
    efficacy.boundary.parameters=list(err.tol=10^-3)
    
    per.stage.sample.sizes <-
      get.sample.sizes.per.stage(n.arms=n.arms, 
                                 n.per.arm=n.per.arm,
                                 subpopulation.sizes=subpopulation.sizes,
                                 interim.info.times=time/max(time),
                                 accrual.rate=accrual.rate,
                                 delay=0,
                                 restrict.enrollment=restrict.enrollment)
    
    analytic.n.per.stage <- per.stage.sample.sizes$analytic.n.per.stage
    total.n.per.stage <- per.stage.sample.sizes$total.n.per.stage
    
    alpha.allocation <- alpha.allocation*total.alpha
    trial.criteria <- trial.criteria.by.scenario <- list()
    
    # List potential trial outcomes
    potential.trials <-
      list.potential.trials(total.n.per.stage=total.n.per.stage,
                            accrual.rate=accrual.rate,
                            subpopulation.sizes=subpopulation.sizes)
    
    
    # Determine which hypotheses are truly null
    if(non.inferiority){
      # In non-inferiority trial, null is treatment effect is inferior
      null.hypotheses <- 1*(hazard.ratios>=ni.margin)
    } else { # Superiority
      # In superiority trial
      null.hypotheses <- 1*(abs(hazard.ratios-1)<1e-6)
      above.mcid <- 1*(hazard.ratios>=mcid)
      above.mcid[which(null.hypotheses==1, arr.ind=TRUE)] <- 0
    }
    
    for(m in 1:n.scenarios){
      ### SPECIFIC TO 2-ARM AND 3-ARM DESIGNS ###
      if(n.arms==2){
        hazard.rate.pop.1 <- hazard.rate[m, c(1,3)]
        hazard.rate.pop.2 <- hazard.rate[m, c(2,4)]
      } else if(n.arms==3){
        hazard.rate.pop.1 <- hazard.rate[m, c(1,3,5)]
        hazard.rate.pop.2 <- hazard.rate[m, c(2,4,6)]
      } else {stop()}
      
      z.parameters <-
        do.call(what=get.z.distribution,
                args=list(analytic.n.per.stage=analytic.n.per.stage, 
                          outcome.type="survival",
                          prop.pop.1 = prop.pop.1,
                          NI=non.inferiority,
                          ni.margin=ni.margin,
                          max.follow=max.follow,
                          enrollment.period=enrollment.period,
                          hazard.rate.pop.1=hazard.rate.pop.1,
                          hazard.rate.pop.2=hazard.rate.pop.2,
                          time=time,
                          restrict.enrollment=restrict.enrollment,
                          censoring.rate=censoring.rate))
      
      # Get efficacy boundaries
      efficacy.boundaries <-
        do.call(what=efficacy.boundary.function,
                args=c(list(cov.mat.used=z.parameters$cov.mat.used,
                            alpha.allocation=alpha.allocation,restrict.enrollment=restrict.enrollment),
                       efficacy.boundary.parameters))
      
      # Sample test statistics 
      test.statistics <- 
        with(z.parameters,
             rmvnorm(n=n.simulations,
                     mean=non.centrality.parameter.vec,
                     sigma=cov.mat.used))
      
      # Determine which hypotheses are rejected at each stage
      test.decisions <-
        do.call(what=test.procedure.function,
                args=list(test.statistics=test.statistics,
                          efficacy.boundary=efficacy.boundaries,
                          ### SPECIFIC TO 2-ARM AND 3-ARM DESIGNS ###
                          futility.boundary=c(futility.boundaries, rep(-Inf,n.arms)),
                          alpha.allocation=alpha.allocation,
                          restrict.enrollment=restrict.enrollment,
                          boundary.to.enroll=boundary.to.enroll))
      
      # Get ESS, Duration, Power, Type I Error, FWER
      trial.criteria.by.scenario[[m]] <-
        get.trial.criteria(decisions.rejections=test.decisions$rejection.matrix,
                           decisions.stages=test.decisions$stage.decision,
                           null.hypotheses=which(null.hypotheses[m,]==1),
                           potential.trials=potential.trials)
      # Label distribution.of.trials by scenario
      trial.criteria.by.scenario[[m]]$distribution.of.trials <-
        data.frame(scenario=m,
                   trial.criteria.by.scenario[[m]]$distribution.of.trials)
    }
    
    # Unpack list of trial criteria - Label with scenario number
    criteria.names <- unique(unlist(lapply(trial.criteria.by.scenario, names)))
    for(i in 1:length(criteria.names)){
      criteria.result <- 
        do.call(what=rbind,
                args=lapply(trial.criteria.by.scenario, function(x)
                  get(criteria.names[i], x)))
      if(is.null(colnames(criteria.result))){
        colnames(criteria.result) <- criteria.names[i]
      }
      trial.criteria[[i]] <- criteria.result
    }
    
    # Unpack list of trial criteria - Label with scenario number
    criteria.names <- unique(unlist(lapply(trial.criteria.by.scenario, names)))
    for(i in 1:length(criteria.names)){
      criteria.result <- 
        do.call(what=rbind,
                args=lapply(trial.criteria.by.scenario, function(x)
                  get(criteria.names[i], x)))
      if(is.null(colnames(criteria.result))){
        colnames(criteria.result) <- criteria.names[i]
      }
      trial.criteria[[i]] <- criteria.result
    }
    
    names(trial.criteria) <- criteria.names
    
    trial.criteria$ess <-
      aggregate(sample.size*proportion~scenario,
                FUN=sum,
                data=trial.criteria$distribution.of.trials)[,2]
    
    trial.criteria$duration <-
      aggregate(duration*proportion~scenario,
                FUN=sum,
                data=trial.criteria$distribution.of.trials)[,2]
    
    trial.criteria$null.hypotheses <- null.hypotheses
    
    if(!non.inferiority){
      trial.criteria$above.mcid <- above.mcid
    }
    
    return(trial.criteria)
  }



### dunnett.wrapper ############################################################
# Description: a generic function for calling dunnett.survival for survival
#   outcomes or dunnett.means for continuous/binary outcomes.
# 
# Input:
#   outcome.type (string): outcome distribution: either "continuous", "binary",
#     or "survival"
#   ... : other arguments passed on to either dunnett.survival or dunnett.
dunnett.wrapper <- function(outcome.type, ...){
  parameters=list(...)
  if(outcome.type=="survival"){
    n.scenarios <- nrow(parameters$hazard.rate)
    n.arms <- 
      length(unique(substring(colnames(parameters$hazard.rate),
                              first=1, last=1)))
    n.subgroups <- length(unique(substring(colnames(parameters$hazard.rate),
                                           first=2)))
    n.stages <- length(parameters$time)
    restrict.enrollment <- parameters$restrict.enrollment
    if(restrict.enrollment==FALSE){
      with(parameters,
         stopifnot(length(futility.boundaries)==
                     (n.stages-1)*(n.arms-1)*n.subgroups,
                   length(alpha.allocation)==
                     n.stages*(n.arms-1)*n.subgroups,
                   n.arms==2,
                   n.subgroups==2,
                   n.stages<=2))
    } else {
      with(parameters,
           stopifnot(length(futility.boundaries)==
                       (n.stages-1)*(n.arms-1)*n.subgroups-1,
                     length(alpha.allocation)==
                       n.stages*(n.arms-1)*n.subgroups-1,
                     n.arms==2,
                     n.subgroups==2,
                     n.stages<=2))
    }
    
    dunnett.survival(...)
  } else if (outcome.type %in% c("continuous", "binary")){
    n.scenarios <- nrow(parameters$outcome.mean)
    n.arms <- 
      length(unique(substring(colnames(parameters$outcome.mean),
                              first=1, last=1)))
    n.subgroups <- length(unique(substring(colnames(parameters$outcome.mean),
                                           first=2)))
    n.stages <- length(parameters$interim.info.times) + 
      1*(tail(parameters$interim.info.times, 1)!=1)
 #   with(parameters,
 #        stopifnot(length(futility.boundaries)==
 #                    (n.stages-1)*(n.arms-1)*n.subgroups,
 #                  length(alpha.allocation)==
 #                    n.stages*(n.arms-2)*n.subgroups,
 #                  n.arms==2,
 #                  n.subgroups==2,
 #                  n.stages<=2))
    
    dunnett.means(outcome.type=outcome.type, ...)
  } else {
    stop("Supported distributions are 'survival', 'continuous', and 'binary'.")
  }
}

### dunnett.means example - continuous outcome two treatment arms #################################
ui.outcome.mean <-
  matrix(c(0.00, 0.00, 0.00, 0.00, 0, 0,
           0.00, 0, 15, 0, 0, 0, 
           0.00, 0,  15, 15, 0, 0,
          0.00, 0,  15, 0, 15, 0,
          0.00, 0,  15, 15, 15, 0,
           0.00, 0, 15, 15, 15, 15),
         ncol=6, byrow=TRUE)
colnames(ui.outcome.mean) <- c("C1", "C2", "A1", "A2", "B1", "B2")
ui.outcome.sd <- matrix(60, ncol=ncol(ui.outcome.mean),
                        nrow=nrow(ui.outcome.mean))
colnames(ui.outcome.sd) <- c("C1", "C2", "A1", "A2", "B1", "B2")

dunnett.wrapper(n.arms= 3,
                n.per.arm= 700,
                subpopulation.sizes=c(0.49, 1-0.49),
                interim.info.times= 1,
                accrual.rate=20,
                outcome.type = "continuous",
                outcome.mean=ui.outcome.mean,
                outcome.sd=ui.outcome.sd,
                delay=0.5,
                mcid=0,
                futility.boundaries=NULL,
                alpha.allocation=c(0.5, 0.5),
                total.alpha=0.025)

### transform.parameters #######################################################
# Description: transform parameters to restricted parameters spaces, such as
#   non-negative reals, positive reals, probabilities, normalized vectors, etc.
#
#
# Input:
#   parameters (list): a named list of parameters which are assumed to take any
#     real value.
#
#   transforms (list): a named list of functions used to transform the 
#     parameters from real values to another set. This can be used to handle
#     parameters that are strictly positive, require normalization, and so on,
#     such as alpha (re)allocation and information times. The names of the list
#     should correspond to the parameters intended for the transform.
#
# Output:
#   tf.params (list): a named list of transformed parameters
#
transform.parameters <- 
  function(parameters,
           transforms) {
  
  # 1. Make sure parameters and transforms are named
  if(sum(nchar(names(parameters))==0) > 0) {
    stop("parameters must be named to avoid ambiguity")
  }
  
  if(sum(nchar(names(transforms))==0) > 0) {
    stop("transforms must be named to avoid ambiguity")
  }
  
  if(!is.null(transforms)){
    
    # 2. Make sure transforms map uniquely to parameters
    which.tf <- match(names(transforms), names(parameters))
    
    if(length(unique(which.tf)) != length(transforms)){
      if(sum(duplicated(which.tf))>0) {
        stop("transforms did not uniquely match to parameters")
      }
      if(sum(is.na(which.tf))>0){
        stop(paste0("parameters did not contain: ",
                    names(transforms)[which(is.na(which.tf))]))
      }
    }
  }
  
  tf.params <- parameters
  
  for(i in 1:length(which.tf)) {
    tf.params[[which.tf[i]]] <- 
      transforms[[i]](parameters[[which.tf[i]]])
  }
  
  return(tf.params)
}

### transform.parameters example ###############################################
# search.parameters <- list(n.per.arm=100,
#                        total.alpha=0.025,
#                        alpha.allocation=c(0, 0, 0, 0))
# search.transforms <- list(alpha.allocation=reals.to.probability)
# 
# transform.parameters(search.parameters,
#                      search.transforms)

### power.penalized.weighted.ess ###############################################
# Description: an objective function that encorporates a weighted average of
#   expected sample size (ESS) across outcome scenarios, with an added penalty
#   for each violation of the power constraints. 
#
# Input:
#   trial.performance (trial performance object): trial performance to evaluate.
#     This should have the following named elements: mss (numeric scalar),
#     ess (numeric vector), duration (numeric vector), power (numeric matrix),
#     fwer (numeric vector), and type.1.error (numeric matrix).
#   scenario.weights (numeric vector): a numeric vector of weights for
#     calculating a weighted mean of ESS over outcome scenarios. If no weight is
#     supplied, the weights are assumed equal.
#   power.constraints (numeric matrix): a numeric matrix of power constraints
#     matching the dimensions of trial.performance$power. If no power constraint
#     is supplied, 80% power is taken as the constraint.
#   power.penalty (numeric scalar): a penalty applied to the difference between
#     power achieved and power desired.
#
#
# Output: an objective function value.
power.penalized.weighted.ess <- 
  function(trial.performance,
           scenario.weights=NULL,
           power.constraints=NULL,
           power.penalty=100,
           objective.scale=1) {
    
    stopifnot("ess" %in% names(trial.performance),
              "empirical.power" %in% names(trial.performance),
              is.numeric(objective.scale) & is.scalar(objective.scale),
              is.finite(objective.scale) & objective.scale > 0)
    
    trial.power <- trial.performance$empirical.power
    
    # For superiority trials, don't penalize for power in alternatives < MCID
    if(exists("above.mcid", where=trial.performance)){
      above.mcid <- trial.performance$mcid
      trial.power[which(above.mcid==0, arr.ind=TRUE)] <- NA
    }
    
    if(is.null(power.constraints)){
      power.constraints <-
        matrix(0.8, nrow=nrow(trial.power), ncol=ncol(trial.power))
    }
    
    if(is.null(scenario.weights)){
      scenario.weights <- rep(1, length(trial.performance$ess))
    }
    
    objective.value <- 
      weighted.mean(trial.performance$ess, w=scenario.weights) +
      sum(linear.threshold.penalty(trial.power,
                                   threshold.matrix=power.constraints,
                                   linear.penalty=power.penalty), na.rm=TRUE)
    
    stopifnot(is.scalar(objective.value))
    return(objective.value/objective.scale)
  }


### get.optim.trajectory ######################################################
# Description: by default, optim() only returns the value of the objective
# function at convergence. By passing options to optim, the iteration history
# can be directed to console, and capture.output() can be used to turn this
# into a character vector. This function takes in such a character vector,
# extracts the iteration history, and places it in a data frame.
#
# Input: the output of a call to capture.output(optim(...))
#
# Output: a data frame with iteration (the iteration number) and obj.fx
# (the objective function at that iteration).

get.optim.trajectory <- 
  function(optim.text) {
    optim.text <- optim.text[grep("iter\\s*\\d*\\s*value\\s*", optim.text)]
    optim.text <- data.frame(iteration=NA, obj.fx=NA, optim.text=optim.text)
    optim.text$iteration <-
      as.numeric(gsub("(iter)|(value.*)", "", optim.text$optim.text))
    optim.text$obj.fx <- 
      as.numeric(gsub("iter\\s*\\d*\\s*value", "", optim.text$optim.text))
    optim.text[c("iteration", "obj.fx")]
  }

### sa.optimize ################################################################
# Description: this is a wrapper for optimizing objects, which requires two
#   functions: a function from the parameters to create an object, and a 
#   function that assigns a value to the object (a loss function that takes
#   in an object and produces a numeric scalar value). This allows optimization
#   of objects using simulated annealing by attempting to minimize the loss for
#   a given set of parameters.
#
#   create.object() is a function that takes the lists search.parameters and 
#   fixed.parameters, and produces an object. Optional transforms specified
#   by search.transforms can constrain this parameter space. The function 
#   evaluate.object() maps the object to a numeric scalar, using the optional
#   parameters specified by ellipsis (...).
#
#   search.parameters = {p1, p2, ..., pd}
#   search.transforms = {f1(), f2(), ... , fd()}
#   transformed.parameters := {x1=f1(p1), x2=f2(p2), ..., xd=fd(pd)
#   object.to.optimize <- create.object(transformed.parameters)
#   objective.value = evaluate.object(object.to.optimize, ...)
#
#   If G is create.object and F is evaluate.object, SA will try to find: 
#     xmin = argmin(F(G(x1, x2, ..., xd)))
#
# Input:
#   search.parameters (list): a named list of all parameters that are to be
#     searched over by simulated annealing, with starting values: passed to
#     create.object()
#
#   search.transforms (list): a named list of functions used to transform the
#     search parameters from real values to another set. This can be used to
#     handle parameters that are strictly positive, require normalization, and
#     so on. If NULL, parameters are assumed to be able to assume any real
#     value. See the documentation of transform.parameters() for more
#     information.
#
#   fixed.parameters (list): a named list of all parameters that are required to
#     fully specify the creation of the object, along with their values: passed
#     to create.object()
#
#   create.object (function): a function which accepts the arguments 
#     search.parameters (after transformation by search.transforms) and
#     fixed.parameters, and returns an object.
#
#   evaluate.object (function): a function that takes an object and assigns it
#     a scalar numeric value, with smaller values being preferred.
#
#   max.iterations (positive numeric scalar): number of iterations of simulated
#     annealing. See ?optim() for more details.
#
#   temperature (positive numeric scalar): 'temperature' parameter for simulated
#     annealing. See ?optim() for more details.
#
#   evals.per.temp (positive numeric scalar): number of evaluations per 
#     temperature. See ?optim() for more details.
#
#   report.iteration (positive numeric scalar): interval for reporting objective
#     function - report every x iterations. See ?optim() for more details.
#
# Output:
# parameters: the list of parameters for the optimized object, after any 
#   transformations were applied.
# optimization.trajectory: a data frame of the objective function values by
#   iteration. See get.optim.trajectory().
# optim.result: the object returned by the call to optim().
sa.optimize <-
  function(search.parameters,
           search.transforms=NULL,
           fixed.parameters=NULL,
           create.object,
           evaluate.object,
           max.iterations,
           temperature,
           evals.per.temp,
           report.iteration,
           function.scale=1,
           parameter.scale=1,
           ...){

    
    # optim() can only take vectors: need to unlist/relist to retain structure.
    initial.param <- as.relistable(search.parameters)
    
    optim.trajectory <-
      capture.output(
        sa.result <- 
          optim(par=unlist(initial.param),
                fn=function(varying,
                            transforms=search.transforms,
                            fixed=fixed.parameters,
                            objective.fx=evaluate.object){ 
                  # Re-structure parameters as list
                  varying <- relist(varying, skeleton=search.parameters)
                  
                  # Call user-specified create.object function: capture any text
                  discard.text <-
                    capture.output(
                      sa.object <-
                        do.call(what=create.object,
                                args=c(list(),
                                       transform.parameters(varying,
                                                            transforms),
                                       fixed))
                    )
                      
                  # Call user-specified objective function: capture any text
                  discard.text <-
                    capture.output(
                      obj.fx <- evaluate.object(sa.object, ...)
                      )
                  obj.fx
                },
                method="SANN",
                control=list(trace=6,
                             maxit=max.iterations,
                             temp=temperature,
                             tmax=evals.per.temp,
                             REPORT=report.iteration,
                             fnscale=function.scale,
                             parscale=parameter.scale))
      )
        
    # Re-structure result as list
    search.parameters <- relist(sa.result$par, skeleton=search.parameters)
    
    # Transform results and return
    parameters <- 
      c(transform.parameters(search.parameters, search.transforms),
        fixed.parameters)
    
    optima <- do.call(what=create.object, args=parameters)

    return(list(parameters=parameters,
                optima=optima,
#                optimization.trajectory=get.optim.trajectory(optim.trajectory),
                optim.result=relist(sa.result$par, skeleton=search.parameters)))
  }

  ### sa.optimize example ########################################################


  ### dunnett.means example - continuous outcome two treatment arms #################################
  ui.outcome.mean <-
    matrix(c(0.00, 0.00, 0.00, 0.00, 0, 0,
             0.00, 0, 15, 0, 0, 0, 
             0.00, 0,  15, 15, 0, 0,
            0.0, 0,  15, 0, 15, 0,
            0, 0, 15, 15, 15, 0,
             0.00, 0, 15, 15, 15, 15),
           ncol=6, byrow=TRUE)
  colnames(ui.outcome.mean) <- c("C1", "C2", "A1", "A2", "B1", "B2")
  ui.outcome.sd <- matrix(60, ncol=ncol(ui.outcome.mean),
                          nrow=nrow(ui.outcome.mean))
  colnames(ui.outcome.sd) <- c("C1", "C2", "A1", "A2", "B1", "B2")

# The following code snippets run every combination of 1-4 stages and optimized vs fixed design. 
# Need to uncomment for the code to run the corresponding optimization.

########################################################################################################################################
#################.  A fixed two stage design 
########################################################################################################################################

  smart.av.trial.two.stage.fixed <-
    sa.optimize(search.parameters=                
                list(n.per.arm=650),
                search.transforms=
                  list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4))),
                fixed.parameters=list(n.arms=3,
                                      interim.info.times = 0.5,
                                      subpopulation.sizes=c(0.49, 0.51),
                                      alpha.allocation=c(0.25, 0.25, 0.25, 0.25),
                                      accrual.rate=20 * 12,
                                      outcome.type="continuous",
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      delay=0.5,
                                      mcid=0,
                                      futility.boundaries= c(0,0,0,0),
                                      total.alpha=0.05),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=1,
                parameter.scale=c(10),
                max.iterations=500,
                temperature=10,
                evals.per.temp=10,
                report.iteration=1,
                power.penalty=10^6)


########################################################################################################################################
#################.  A one stage design with equal alpha allocation
########################################################################################################################################

  smart.av.trial.one.stage.equal.alpha <-
    sa.optimize(search.parameters=                
                list(n.per.arm=750),
                search.transforms=
                  list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4))),
                fixed.parameters=list(n.arms=3,
                                      interim.info.times = NULL,
                                      subpopulation.sizes=c(0.49, 0.51),
                                      alpha.allocation=c(0.5, 0.5),
                                      accrual.rate=20 * 12,
                                      outcome.type="continuous",
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      delay=0.5,
                                      mcid=0,
                                      futility.boundaries= NULL,
                                      total.alpha=0.05),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=1,
                parameter.scale=c(10),
                max.iterations=500,
                temperature=10,
                evals.per.temp=10,
                report.iteration=1,
                power.penalty=10^6)

########################################################################################################################################
#################.  A one stage design with optimized alpha allocation
########################################################################################################################################


  smart.av.trial.one.stage.optimized.alpha <-
    sa.optimize(search.parameters=                
                list(n.per.arm=630,
                       alpha.allocation=c(0, 0)),
                search.transforms=
                  list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4)),                                               
                       alpha.allocation=reals.to.probability),
                fixed.parameters=list(n.arms=3,
                                      interim.info.times = NULL,
                                      futility.boundaries = NULL,
                                     subpopulation.sizes=c(0.49, 0.51),
                                      accrual.rate=20 * 12,
                                      outcome.type="continuous",
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      delay=0.5,
                                      mcid=0,
                                      total.alpha=0.05),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=1,
                parameter.scale=c(10, rep(1, 2)),
                max.iterations=500,
                temperature=10,
                evals.per.temp=10,
                report.iteration=1,
                power.penalty=10^6)



########################################################################################################################################
#################.  A two stage design that optimizes alpha.allocation, interim analysis timing and futility boundaries
########################################################################################################################################

  smart.av.trial.optimized.adaptive <-
    sa.optimize(search.parameters=                
                list(n.per.arm=650,
                       interim.info.times=0,
                       alpha.allocation=c(0, 0, 0, 0),
                       futility.boundaries=c(-1, -1, -1, -1)),
                search.transforms=
                  list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4)),
                       interim.info.times=function(x)
                         squash(plogis(x),0.1, 0.9),
                       alpha.allocation=reals.to.probability),
                fixed.parameters=list(n.arms=3,
                                      subpopulation.sizes=c(0.49, 0.51),
                                      accrual.rate=20 * 12,
                                      outcome.type="continuous",
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      delay=0.5,
                                      mcid=0,
                                      total.alpha=0.05),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=1,
                parameter.scale=c(10, 1, rep(1, 4), rep(1, 4)),
                max.iterations=500,
                temperature=10,
                evals.per.temp=10,
                report.iteration=1,
                power.penalty=10^6)


########################################################################################################################################
#################.  A three stage design that optimizes alpha.allocation, interim analysis timing and futility boundaries
########################################################################################################################################

  smart.av.trial.optimized.three.stage.adaptive <-
        sa.optimize(search.parameters=
                                    list(n.per.arm=650,
                                                                interim.info.times=c(0.3, 0.5),
                                                                alpha.allocation=c(0, 0, 0, 0, 0, 0),
                                                                futility.boundaries=c(-1, -1, -1, -1, -1, -1, -1, -1)),
                                    search.transforms=
                                      list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4)),
                                                                  interim.info.times=function(x)
                                                                    squash(sort(plogis(x)),0.1, 0.9),
                                                                  alpha.allocation=reals.to.probability),
                                    fixed.parameters=list(n.arms=3,
                                                                              subpopulation.sizes=c(0.49, 0.51),
                                                                              accrual.rate=20 * 12,
                                                                              outcome.type="continuous",
                                                                              outcome.mean=ui.outcome.mean,
                                                                              outcome.sd=ui.outcome.sd,
                                                                              delay=0.5,
                                                                              mcid=0,
                                                                              total.alpha=0.05),
                                    create.object=dunnett.wrapper,
                                    evaluate.object=power.penalized.weighted.ess,
                                    function.scale=1,
                                    parameter.scale=c(10, rep(1,2), rep(1, 6), rep(1, 8)),
                                    max.iterations=500,
                                    temperature=10,
                                    evals.per.temp=10,
                                    report.iteration=1,
                                    power.penalty=10^6)

########################################################################################################################################
#################.  A fixed three stage design with interim analysis at times 1/3 and 2/3
########################################################################################################################################

  smart.av.trial.three.stage.fixed <-
        sa.optimize(search.parameters=
                                    list(n.per.arm=650),
                                    search.transforms=
                                      list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4))),
                                    fixed.parameters=list(n.arms=3,
                                                                              interim.info.times = c(1/3, 2/3),
                                                                              subpopulation.sizes=c(0.49, 0.51),
                                                                              alpha.allocation=rep(1/6, 6),
                                                                              accrual.rate=20 * 12,
                                                                              outcome.type="continuous",
                                                                              outcome.mean=ui.outcome.mean,
                                                                              outcome.sd=ui.outcome.sd,
                                                                              delay=0.5,
                                                                              mcid=0,
                                                                              futility.boundaries= c(0,0,0,0,0,0,0,0),
                                                                              total.alpha=0.05),
                                    create.object=dunnett.wrapper,
                                    evaluate.object=power.penalized.weighted.ess,
                                    function.scale=1,
                                    parameter.scale=c(10),
                                    max.iterations=500,
                                    temperature=10,
                                    evals.per.temp=10,
                                    report.iteration=1,
                                    power.penalty=10^6)

########################################################################################################################################
#################.  An optimized four stage design
######################################################################################################################################## 

  smart.av.trial.optimized.four.stage.adaptive <-
        sa.optimize(search.parameters=
                                    list(n.per.arm=650,
                                                                interim.info.times=c(0, 1/3,2/3),
                                                                alpha.allocation=c(0, 0, 0, 0, 0, 0, 0, 0),
                                                                futility.boundaries=c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)),
                                    search.transforms=
                                      list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4)),
                                                                  interim.info.times=function(x)
                                                                    squash(sort(plogis(x)),0.1, 0.9),
                                                                  alpha.allocation=reals.to.probability),
                                    fixed.parameters=list(n.arms=3,
                                                                              subpopulation.sizes=c(0.49, 0.51),
                                                                              accrual.rate=20 * 12,
                                                                              outcome.type="continuous",
                                                                              outcome.mean=ui.outcome.mean,
                                                                              outcome.sd=ui.outcome.sd,
                                                                              delay=0.5,
                                                                              mcid=0,
                                                                              total.alpha=0.05),
                                    create.object=dunnett.wrapper,
                                    evaluate.object=power.penalized.weighted.ess,
                                    function.scale=1,
                                    parameter.scale=c(10, rep(1,3), rep(1, 8), rep(1, 12)),
                                    max.iterations=500,
                                    temperature=10,
                                    evals.per.temp=10,
                                    report.iteration=1,
                                    power.penalty=10^6)

########################################################################################################################################
#################.  A fixed four stage design with interim analysis at times 1/4 and 2/4, 3/4
########################################################################################################################################

  smart.av.trial.four.stage.fixed <-
        sa.optimize(search.parameters=
                                    list(n.per.arm=650),
                                    search.transforms=
                                      list(n.per.arm=function(x) ceiling(squash(x, 50, 1e4))),
                                    fixed.parameters=list(n.arms=3,
                                                                              interim.info.times = c(1/4, 2/4, 3/4),
                                                                              subpopulation.sizes=c(0.49, 0.51),
                                                                              alpha.allocation=rep(1/8, 8),
                                                                              accrual.rate=20 * 12,
                                                                              outcome.type="continuous",
                                                                              outcome.mean=ui.outcome.mean,
                                                                              outcome.sd=ui.outcome.sd,
                                                                              delay=0.5,
                                                                              mcid=0,
                                                                              futility.boundaries= c(0,0,0,0,0,0,0,0,0,0,0,0),
                                                                              total.alpha=0.05),
                                    create.object=dunnett.wrapper,
                                    evaluate.object=power.penalized.weighted.ess,
                                    function.scale=1,
                                    parameter.scale=c(10),
                                    max.iterations=500,
                                    temperature=10,
                                    evals.per.temp=10,
                                    report.iteration=1,
                                    power.penalty=10^6)
