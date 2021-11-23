# 05/15/2019
# utility functions for inference on incomplete data

library(Matrix)
#library(igraph)
#library(ggplot2)

source("./Truncated_Exponential.R")
 

# 1. utility functions for evaluating likelihood
# 1.1 trace system status to a given time point
# default model: SIR
obtain_system <- function(G0, I0, events, t.obs, 
                          model="SIR", quarantine=T){
  tmax =  max(events$time)
  if(t.obs >  tmax){
    cat("Required time out of range, use tmax =",tmax, "instead.\n")
    t.obs = tmax
  }
  rowmax = min(which(events$time >= t.obs))
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; infected =  I0; epid[infected] = 1
  
  for(r in 1:rowmax){
    z =  events[r,]$event
    if (z==1){
      # infection
      p1 = events[r,]$per1
      epid[p1] = 1
    }else if (z==2){
      # recovery
      p1 = events[r,]$per1
      if(model=="SIR"){
        epid[p1] = -1
      }else{
        epid[p1] = 0
      }
    }else{
      # some edge stuff
      p1 = events[r,]$per1
      p2 = events[r,]$per2
      if(quarantine){
        if(z %in% c(3:5)){
          # reconnection
          adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        }else{
          # disconnection
          adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        }
      }else{
        if(z==3){
          adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        }else{
          adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        }
      }
    }
    # report progress
    cat("Processing...", r/rowmax,"\r")
  }
  
  # also return useful summary statistics
  N.I = sum(epid==1); N.S = sum(epid==0)
  # cat("At time", t.obs, ":\n",
  #     "epid vector =", epid, "\n")
  if(quarantine){
    susceptible = which(epid==0)
    infected = which(epid==1)
    Mmax.SS = N.S*(N.S-1)/2
    Mmax.SI = N.S * N.I
    Mmax.II = N.I*(N.I-1)/2
    if(Mmax.SS==0){
      M.SS = 0
    }else{
      M.SS = nnzero(adjmat[susceptible,susceptible])/2
    }
    if(Mmax.SI==0){
      M.SI = 0
    }else{
      M.SI = nnzero(adjmat[susceptible,infected])
    }
    if(Mmax.II==0){
      M.II = 0
    }else{
      M.II = nnzero(adjmat[infected,infected])/2
    }
    
    # stats = c(N.S, N.I, 
    #           Mmax.SS, Mmax.SI, Mmax.SS,
    #           M.SS, M.SI, M.II)
    # names(stats) =  c("N.S","N.I","Mmax.SS","Mmax.SI","Mmax.SS",
    #                   "M.SS","M.SI","M.II")
    
    stats = data.frame(N.S = N.S, N.I = N.I, 
                       Mmax.SS=Mmax.SS, Mmax.SI=Mmax.SI, Mmax.II=Mmax.II,
                       M.SS=M.SS, M.SI=M.SI, M.II=M.II)
    
  }else{
    Mmax = N*(N-1)/2
    M =  nnzero(adjmat)
    
    if(N.S > 0 & N.I > 0){
      N.SI = nnzero(adjmat[susceptible,infected])
    }else{
      N.SI = 0
    }
    
    # stats  = c(N.S, N.I, Mmax, M)
    # names(stats) = c("N.S","N.I","Mmax", "M")
    
    stats  = data.frame(N.S = N.S, N.I = N.I, N.SI = N.SI, 
                        Mmax = Mmax, M = M)
  }
  
  return(list(adjmat = adjmat, epid = epid,
              stats = stats))
}


# 1.2 evaluate the log likelihood of events in a given interval [st,en)
## use log-likelihood for better computation accuracy
## only deal with the coupled case for now
## 05/09/2019
## augmentation: calculate important stats (for inference) as well
eval_interval_loglik <- function(dat, st, en, model="SIR", quarantine=T,
                                 bet = 0.03, gam = 0.15, 
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05),
                                 cal_stats = T){
  
  params = c(bet, gam, alpha.r, alpha.d)
  if(quarantine & length(params) < 8){
    stop("Type-dependent edge rates should be specified! What are you using instead??\n")
  }
  if(!quarantine & length(params) > 4){
    stop("Type-independent edge rates should be specified! What are you using instead??\n")
  }
  
  G0 = dat$G0; I0 = dat$I0; events = dat$events
  tmax =  max(events$time)
  if(st >= en){
    stop("Invalid time interval.\n")
  }
  if(en >  tmax){
    cat("Required ending time out of range, use tmax =",tmax, "instead.\n")
    en = tmax
  }
  
  rowmin = min(which(events$time > st))
  rowmax = max(which(events$time <= en))
  
  if(rowmin > 1){
    sys = obtain_system(G0, I0, events, st, model, quarantine)
    adjmat = sys$adjmat; epid = sys$epid; stats = sys$stats
    
    t.cur = events$time[rowmin-1]
  }else{
    adjmat = G0; epid = rep(0,nrow(G0)); epid[I0] = 1
    susceptible = which(epid==0); infected = which(epid==1)
    N.S = sum(epid==0); N.I = sum(epid==1)
    stats = data.frame(N.S = N.S, N.I = N.I,
                       Mmax.SS = N.S*(N.S-1)/2, 
                       Mmax.SI = N.S * N.I, 
                       Mmax.II = N.I*(N.I-1)/2,
                       M.SS = nnzero(adjmat[susceptible,susceptible])/2, 
                       M.SI = nnzero(adjmat[susceptible,infected]), 
                       M.II = nnzero(adjmat[infected,infected])/2)
    
    t.cur = 0
  }
  
  N.S = stats$N.S; N.I = stats$N.I
  if(quarantine){
    # basic edge stats
    N.SI = stats$M.SI
    M = c(stats$M.SS, stats$M.SI, stats$M.II)
    Mmax= c(stats$Mmax.SS, stats$Mmax.SI, stats$Mmax.II)
    #M.d = Mmax - M
    cat("At the start...\n",
        "N.S:",N.S, "N.I:", N.I,"\n",
        "M:",M,"\n",
        "Mmax:",Mmax,"\n",
        "Md:",Mmax-M,"\n")
  }else{
    stop("This function doesn't deal with the decoupled case for now.\n")
  }
  
  n.E = sum(events$event[rowmin:rowmax] == 1)
  n.R = sum(events$event[rowmin:rowmax] == 2)
  C = rep(0,3); D = rep(0,3)
  
  logsum = 0; expo = 0
  if(cal_stats){
    # data storage for "big.sums"
    big.sums = numeric(length(params))
  }else{
    big.sums = 0
  }
  
  for(r in rowmin:rowmax){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    
    # update "big.sums"
    if(cal_stats){
      big.sums = big.sums + params * c(N.SI, N.I, Mmax-M, M) * del.t
    }
    # update the exponent part (excluding "-")
    expo = expo + 
      sum(params * c(N.SI, N.I, Mmax-M, M)) * del.t
    
    z =  events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # update logsum term
      logsum = logsum + log(I.p1)
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      M[1] = M[1] - S.p1
      M[2] = M[2] + S.p1 - I.p1
      M[3] = M[3] + I.p1
      N.SI = M[2]
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        epid[p1] = -1
        M[1] = M[1]
        M[2] = M[2] - S.p1
        M[3] = M[3] - I.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
      }
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      N.SI = M[2]
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]

      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          logsum = logsum + log(Mmax[1]-M[1])
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          logsum = logsum + log(Mmax[3]-M[3])
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else{
          # S-I type
          logsum = logsum + log(Mmax[2]-M[2])
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = M[2]
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          logsum = logsum + log(M[1])
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          logsum = logsum + log(M[3])
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else{
          # S-I type
          logsum = logsum + log(M[2])
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = M[2]
        }
      }
      
    }
    # report progress
    cat("Processing...", 
        (r-rowmin)/(rowmax-rowmin),"\r")
    
    t.cur = t.next
    
  }
  cat("\nDone.\n")
  cat("Event counts:",c(n.E, n.R, C, D),"\n",
      "logsum:",logsum,"\n",
      "expo:",expo,"\n")
  
  # calculate log-lik
  ll = logsum + 
    sum(log(params)*c(n.E, n.R, C, D))-
    expo
  
  res = list(ll=ll, event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}


# 2.3 a function that augments the observed data w/ imputed recovery times
# and outputs event counts & all the summations
# assume a coupled process
# **deprecated version**
# it disregards all the edge changes related to R people
parse_augment <- function(G0, I0, events, recovers, model="SIR"){
  # dats: a list consisting of everything observed
  # recovers: a dataset w/ variables `time` & `per1`, all imputed recovery times

  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; infected =  I0; epid[infected] = 1
  
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2) + nrow(recovers)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the case of "quarantine" - coupled case
  big.sums = numeric(8) 
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  events = events[,-c(5:6)]
  recovers$per2 = NA
  recovers$event = 2
  recovers = recovers[,names(events)]
  events = rbind(events, recovers)
  events = events[order(events$time),]
  
  # parse through the augmented data
  t.cur = 0
  for(r in 1:nrow(events)){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t

    z =  events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      M[1] = M[1] - S.p1
      M[2] = M[2] + S.p1 - I.p1
      M[3] = M[3] + I.p1
      N.SI = M[2]
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        epid[p1] = -1
        M[1] = M[1]
        M[2] = M[2] - S.p1
        M[3] = M[3] - I.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
      }
      Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
      N.SI = M[2]
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
     
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = M[2]
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]==0 & epid[p2]==0){
          # S-S type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = M[2]
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next

  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}


# 05/23/2019
# the actually used version
# modify the `parse_augment` function to include R-related edge changes
# 06/02/2019
# debugged
parse_augment2 <- function(G0, I0, events, recovers, model="SIR"){
  # dats: a list consisting of everything observed
  # recovers: a dataset w/ variables `time` & `per1`, all imputed recovery times
  
  N = nrow(G0); epid = rep(0,N)
  adjmat =  G0; epid[I0] = 1
  
  susceptible = which(epid==0); infected = which(epid==1)
  N.S = sum(epid==0); N.I = sum(epid==1)
  if(model == "SIR"){ N.R = 0 }
  Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
  M = c(nnzero(adjmat[susceptible,susceptible])/2, 
        nnzero(adjmat[susceptible,infected]), 
        nnzero(adjmat[infected,infected])/2)
  N.SI = M[2]
  
  # data storage
  n.E = sum(events$event == 1)
  n.R = sum(events$event == 2) + nrow(recovers)
  C = rep(0,3); D = rep(0,3)
  
  # FOR NOW
  # assume the coupled case 
  big.sums = numeric(8) 
  
  # combine the two datasets to make a "complete" dataset
  # and order by ascending time
  events = events[,-c(5:6)]
  recovers$per2 = NA
  recovers$event = 2
  recovers = recovers[,names(events)]
  events = rbind(events, recovers)
  events = events[order(events$time),]
  
  # parse through the augmented data
  t.cur = 0
  for(r in 1:nrow(events)){
    t.next = events$time[r]
    del.t = t.next - t.cur
    
    # update "big.sums"
    big.sums = big.sums + c(N.SI, N.I, Mmax-M, M) * del.t
    
    z = events$event[r]
    if (z==1){
      # infection
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      epid[p1] = 1
      N.S = N.S - 1
      N.I = N.I + 1
      # update counts
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1 - R.p1
        M[2] = M[2] + S.p1 + R.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = N.SI + S.p1 - I.p1
      }else{
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        M[1] = M[1] - S.p1
        M[2] = M[2] + S.p1 - I.p1
        M[3] = M[3] + I.p1
        N.SI = M[2]
      }
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      # figure out neighborhood
      I.p1 = nnzero(adjmat[p1, which(epid==1)])
      S.p1 = nnzero(adjmat[p1, which(epid==0)])
      # bookkeeping
      N.I = N.I - 1
      if(model=="SIR"){
        R.p1 = nnzero(adjmat[p1, which(epid==-1)])
        epid[p1] = -1
        N.R = N.R + 1
        M[1] = M[1] + S.p1 + R.p1
        M[2] = M[2] - S.p1 - R.p1 + I.p1
        M[3] = M[3] - I.p1
        #Mmax = c((N.S+N.R)*(N.S+N.R-1)/2, (N.S+N.R) * N.I, N.I*(N.I-1)/2)
        Mmax = c((N-N.I)*(N-N.I-1)/2, (N-N.I) * N.I, N.I*(N.I-1)/2)
        N.SI = N.SI - S.p1
      }else{
        epid[p1] = 0
        N.S = N.S + 1
        M[1] = M[1] + S.p1
        M[2] = M[2] - S.p1 + I.p1
        M[3] = M[3] - I.p1
        Mmax = c(N.S*(N.S-1)/2, N.S * N.I, N.I*(N.I-1)/2)
        N.SI = M[2]
      }
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
  
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          C[1] = C[1] + 1
          M[1] = M[1] + 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          C[3] = C[3] + 1
          M[3] = M[3] + 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type: neither is recovered
          C[2] = C[2] + 1
          M[2] = M[2] + 1
          N.SI = N.SI + 1
        }else{
          # R-I type
          C[2] = C[2] + 1
          M[2] = M[2] + 1
        }
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        # bookkeeping
        if (epid[p1]!=1 & epid[p2]!=1){
          # S-S/S-R/R-R type
          D[1] = D[1] + 1
          M[1] = M[1] - 1
        }else if (epid[p1]==1 & epid[p2]==1){
          # I-I type
          D[3] = D[3] + 1
          M[3] = M[3] - 1
        }else if (epid[p1]!=-1 & epid[p2]!=-1){
          # S-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
          N.SI = N.SI - 1
        }else{
          # R-I type
          D[2] = D[2] + 1
          M[2] = M[2] - 1
        }
      }
      
    }
    # report progress
    # cat("Processing augmented dataset...", 
    #     r/nrow(events),"\r")
    
    t.cur = t.next
    
  }
  #cat("\nDone.\n")
  
  res = list(event.counts = c(n.E, n.R, C, D), big.sums = big.sums)
  
  return(res)
}



# 2. Impute missing recovery times on intervals with missingness
# 2.1 a function to obtain time intervals to do imputation on, and involved individuals
## return a table of interval boundaries, (st, en)'s, and a list of vectors of individual labels
get_miss_recov <- function(report, times, events){
  # data storage
  lbs = NULL; ubs = NULL
  miss_recov = list()
  # go through all intervals and check
  nt = length(times)
  for(ix in 2:nt){
    lb = times[ix-1]; ub = times[ix]
    epid.change = report[ix,] - report[ix-1,]
    recovered = which(epid.change < 0)
    
    cat("Interval:",lb,ub,"\n")
    cat("Recovered:",recovered,"\n")
    
    if(length(recovered) > 0){
      # check if they are included in the event log already
      st = min(which(events$time > lb)); en = max(which(events$time <= ub))
      events.sel = events[st:en,]
      events.sel = events[events$event == 2,]
      if(nrow(events.sel) > 0){
        exact.recovered = events.sel$per1
        recovered = recovered[! recovered %in% exact.recovered]
      }
      # record if 'recovered' non-empty
      if(length(recovered) > 0){
        lbs = c(lbs,lb); ubs = c(ubs, ub)
        miss_recov = append(miss_recov, list(recovered))
      }
    }
  }
  return(list(intervals = data.frame(lb=lbs, ub=ubs), recover = miss_recov))
}


# 2.2 a function to obtain the INFECTED neighbors AT TIME OF INFECTION for each individual
# and return a list for permanent storage
# "Infected": up to the latest status report + new infections since then
# option: only return results for a subset of individuals
get_nei_infection <- function(G0, events, report, times, subset=NULL){
  adjmat = G0
  events = events[events$event != 2,]
  nei_infec = list()
  new_infecs = NULL
  lb = 0
 
  # updat adjmat on the fly
  # and record whenever a POI gets infected
  for(r in 1:nrow(events)){
    z =  events$event[r]
    if(z==1){
      # infection: record neighborhood
      p1 = events$per1[r]
      t = events$time[r]
      # keep track of new infections since latest status report
      lb.t = max(which(times <= t))
      if(lb.t == lb){
        new_infecs = c(new_infecs, p1)
      }else{
        new_infecs = p1
        lb = lb.t
      }
      if(is.null(subset) | (!is.null(subset) & p1 %in% subset)){
        nei = which(adjmat[p1,]==1)
        ix = max(which(times <= t))
        infected.pre = which(report[ix,]==1)
        infected = union(infected.pre,new_infecs)
        nei_infec[[as.character(p1)]] = intersect(nei, infected)
      }
      # if(!is.null(subset) & p1 %in% subset){
      #   nei_infec[[as.character(p1)]] = which(adjmat[p1,]==1)
      # }
      # if(is.null(subset)){
      #   nei_infec[[as.character(p1)]] = which(adjmat[p1,]==1)
      # }
    }else{
      # edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
      # assume type-dependent edge evolution
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
      }
    }
  }
  return(nei_infec)
}



# 2.3 propose recovery times for a given time interval
# 2.3.a "REJECT": propose times from TE and keep rejecting if not consistent with infection trajectory
propose_recov_rej <- function(lb, ub, recovers, events, nei_infec, gam=0.2){
  
  st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  
  # pull up infection log in this interval 
  events.infec = events[st:en,]
  events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # propose initial candicate times
  cands = sam_texp(length(recovers), gam, 0, ub-lb) + lb
  
  #cat("Interval: [",lb,",",ub,"]\n")
  #cat("To recover:",recovers,"\n")
  #cat("infection log:\n")
  #print(events.infec)
  #cat("Proposed recovery times:\n",cands,"\n")
  
  # a sub-function for consistency-checking
  check_consist <- function(cands){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      # check self: can't recover before infection
      if(p1 %in% recovers){
        t = events.infec$time[r]
        if(cands[recovers == p1] <= t){
          return(FALSE)
        }
      }
      # then check neighborhood
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        t = events.infec$time[r]
        poi.t = cands[recovers %in% poi]
        if(all(poi.t <= t)){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  
  # if no infection cases to accommodate at all, directly return the proposal
  if(nrow(events.infec)==0){
    return(cands)
  }
  
  # otherwise, rejection sampling
  while (!check_consist(cands)) {
    #cat("Rejected!\n")
    cands = sam_texp(length(recovers), gam, 0, ub-lb) + lb
    #cat("Proposed recovery times:\n",cands,"\n")
  }
  #cat("Accepted!\n")
  return(cands)
}


# 2.3.b "MH": propose times from TE and keep the previous sample if not consistent with infection trajectory
# return something if it is viable, otherwise return a vector of NA's
propose_recov_MH <- function(lb, ub, recovers, events, nei_infec, gam=0.2){
  
  st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  
  # pull up infection log in this interval 
  events.infec = events[st:en,]
  events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # propose candicate times
  cands = sam_texp(length(recovers), gam, 0, ub-lb) + lb
  
  #cat("Interval: [",lb,",",ub,"]\n")
  #cat("To recover:",recovers,"\n")
  #cat("infection log:\n")
  #print(events.infec)
  #cat("Proposed recovery times:\n",cands,"\n")
  
  # a sub-function for consistency-checking
  check_consist <- function(cands){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      # check self: can't recover before infection
      if(p1 %in% recovers){
        t = events.infec$time[r]
        if(cands[recovers == p1] <= t){
          return(FALSE)
        }
      }
      # then check neighborhood
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        t = events.infec$time[r]
        poi.t = cands[recovers %in% poi]
        if(all(poi.t <= t)){
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }
  
  # if no infection cases to accommodate at all, directly return the proposal
  if(nrow(events.infec)==0){
    return(cands)
  }
  
  # otherwise, return cands only if they are consistent with observed data
  if(check_consist(cands)){
    return(cands)
  }else{
    return(rep(NA, length(recovers)))
  }
}

# 2.3.c "CHEWBACCA": filter through each infected person's neighborhood to ensure consistency
propose_recov_filter <- function(lb, ub, recovers, events, nei_infec, gam=0.2){
  
  st = min(which(events$time > lb)); en = max(which(events$time <= ub))
  
  # pull up infection log in this interval 
  events.infec = events[st:en,]
  events.infec = events.infec[events.infec$event == 1,c("time","per1")]
  
  # if events.infec non-empty:
  # obtain adjusted feasible sampling lower bounds for each person in `recovers`
  # if about-to-recover people are the only sick neighbors of someone infected at t, 
  # then randomly select one of them to mandately recover after t
  bounds = rep(lb, length(recovers))
  if(nrow(events.infec) > 0){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      if(p1 %in% recovers){
        t = events.infec$time[r]
        bounds[recovers==p1] = max(bounds[recovers==p1],t)
      }
      nei = nei_infec[[as.character(p1)]]
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        #cat("POIs for individual",p1,":\n",poi,"\n")
        t = events.infec$time[r]
        if(length(poi)==1){
          p = poi
        }else{
          p = sample(poi, 1)
        }
        bounds[recovers==p] = max(bounds[recovers==p],t)
      }
    }
  }
  
  # sample recovery times under the adjusted bounds
  cands = sam_texp(length(recovers), gam, bounds-lb, ub-lb) + lb
  
  # cat("Interval: [",lb,",",ub,"]\n")
  # cat("To recover:",recovers,"\n")
  # cat("Feasible lower bounds:\n",bounds,"\n")
  # cat("Proposed recovery times:\n",cands,"\n")
  
  return(cands)
}



# ## bench mark the two functions (REJECT v.s. CHEWBACCA)
# bp_res =
# bench::press(
#   ix = c(1:length(recovers)),
#   {
#     bench::mark(
#       length(propose_recov_rej(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
#                                   events = miss1$events, nei_infec = nei_infec_miss)),
#       length(propose_recov_filter(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
#                            events = miss1$events, nei_infec = nei_infec_miss))
#     )
#   }
# )
# 
# bp_tab = bp_res %>% dplyr::select(ix, min, median)
# bp_tab$num_recov = unlist(bp_res$result)
# bp_tab$method = rep(c("rejection","domain"),5)
# bp_tab$min = as.character(bp_tab$min)
# bp_tab$median = as.character(bp_tab$median)
# library(xtable)
# xtable(bp_tab)
# ## not much difference;
# ## but when #(imputation) is big or contraints are complex,
# ## the 'filter' function is more efficient
