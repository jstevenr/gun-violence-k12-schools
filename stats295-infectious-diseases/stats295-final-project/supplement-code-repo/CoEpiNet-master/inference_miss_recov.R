# 05/16/2019
# Inference on partially observed data
# --Missing (some) recovery times 

# 05/22/2019
# add another function that does the entire process:
# 1) read in a complete dataset
# 2) generate missingness in recovery times
# 3) do inference (and make plots) in 3 different imputation methods
# 4) calculate ESS and conduct Geweke diagnostics

# preparations: loading stuff

source("./inference_util.R")
source("./pre_process.R")

# try some other prior
# just to shake things up
pr = data.frame(count = rep(1,4), avg = c(0.02, 0.1, 0.004, 0.06))


# 1. function to do Bayesian inference on data w/ missing recovery times
# output/plot results every `output.sams` recorded samples
# priors: data frame of vars `count` & `avg`
# assume "quarantine" case!

# 05/21/2019
# add a third method for proposing recovery times: ``..._MH'': 
# keep the previous sample if new proposal is the new proposal is no good

# 05/22/2019
# pull true params from the dataset, if such info is available

# 06/04/2019
# use `parse_augment2` function 
infer_miss_recov <- function(dats, priors, init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dats$G0; I0 = dats$I0; events = dats$events
  reports = dats$report; report.times = dats$report.times
  
  if("truth" %in% names(dats)){
    true.params = dats$truth
  }else{
    true.params= c(.03, .15, .005, 0.001, .005, .05, 0.1, .05)
  }
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, reports, report.times)
  
  # process the priors and initialize param values
  if(nrow(priors)==4){
    a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
    avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
    b.pr = a.pr/avg.pr
  }else{
    a.pr = priors$count
    b.pr = a.pr/priors$avg
  }
  if(is.null(init.params)){
    # draw from priors
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # a "use new proposal" flag for the MH method
  # use.new=T -> need to recompute suff.stats.
  if(impute=="MH"){
    use.new = T
  }
  
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    if(impute == "MH"){
      if(it > 1){
        prev.times = times
      }else{
        prev.times = NULL
      }
    }
    times = NULL
    for(ix in 1:length(MR$recover)){
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      if (impute == "filter"){
        imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }else if (impute == "MH"){
        # use filter to get a viable imputation in the 1st interation
        if(it == 1){
          imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }else{
          imputed = propose_recov_MH(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
        }
      }else{
        imputed = propose_recov_rej(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }
      times = c(times, imputed)
    }
    if(impute == "MH" & it > 1){
      # update the flag
      use.new = !all(is.na(times))
      # use the previous proposals if new proposals are not good
      # for each particular interval!
      times[is.na(times)] = prev.times[is.na(times)]
    }
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    #print(recover.dat)
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    # (2) compute event counts and summations
    #PA = parse_augment(G0, I0, events, recover.dat)
    PA = parse_augment2(G0, I0, events, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if(it > burn & (it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  # traceplot
  # modification: add in 95% posterior credible intervals
  for(ix in 1:8){
    v = vars[ix]
    sams = params[,ix]
    if(!is.null(true.params)){
      truth = true.params[ix]
      yl = range(sams, truth)
      plot(sams~c(1:samples), main=paste("Traceplot for",v), 
           xlab = "sample", ylab = v, type="l", ylim = yl)
      abline(h=truth,col="red",lwd=2)
    }else{
      yl = range(sams)
      plot(sams~c(1:samples), main=paste("Traceplot for",v), 
           xlab = "sample", ylab = v, type="l", ylim = yl)
    }
    bounds = quantile(sams, c(.025, .975))
    abline(h=bounds[1],col="gray",lty=2,lwd=2)
    abline(h=bounds[2],col="gray",lty=2,lwd=2)
  }
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  return(params)
}




# 2. the entire process function

# 2.1 a diagnostic function for inference results 
# (a matrix with each column being a chain for a particular parameter)
infer.diag <- function(res, method, plot=F){
  vars = colnames(res)
  
  ess = coda::effectiveSize(res)#/nrow(res)
  zscore = coda::geweke.diag(res)$z
  pvals = sapply(zscore, function(z) ifelse(z<0, pnorm(z), 1-pnorm(z))) * 2
  if(plot){
    par(mfrow=c(2,2))
    res = coda::as.mcmc(res)
    coda::geweke.plot(res, auto.layout = F)
  }
  res.dat = as.data.frame(rbind(ess, zscore, pvals))
  names(res.dat) = vars
  res.dat$method = as.character(method)
  row.names(res.dat) = c("ESS","z.score","p.values")
  return(res.dat)
}

# 06/03/2019
# modified: able to set pdfnames
# (but all plots are saved to a pdf file)
pipeline_miss_recov <- function(datname, fpath = "~/Documents/",
                                interval = 7, miss_prop = 1, miss_model = "SIR",
                                doMH = T, save_miss = T, 
                                pdfname = NULL, ...){
  dats = readRDS(paste0(fpath,datname,".rds"))
  miss_dats = miss_recovery(dats, interval, miss_prop, miss_model)
  if(save_miss){
    misspath = paste0(fpath,datname,"_miss",miss_prop*100,".rds")
    saveRDS(miss_dats, misspath)
  }
  
  diag.dat = NULL
  
  par(mfrow=c(2,2))
  if(is.null(pdfname)){
    pdf(paste0(fpath,datname,".pdf"), width = 9, height = 6)
  }else{
    pdf(paste0(fpath,pdfname,".pdf"), width = 9, height = 6)
  }
  
  res.fil = infer_miss_recov(dats = miss_dats, impute = "filter", ...)
  diag.dat = rbind(diag.dat, infer.diag(res.fil, "filter"))
  res.rej = infer_miss_recov(dats = miss_dats, impute = "rej", ...)
  diag.dat = rbind(diag.dat, infer.diag(res.rej, "reject"))
  if(doMH){
    res.MH = infer_miss_recov(dats = miss_dats, impute = "MH", ...)
    diag.dat = rbind(diag.dat, infer.diag(res.MH, "MH"))
  }
  
  dev.off()
  
  return(diag.dat)
}

