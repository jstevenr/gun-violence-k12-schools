
# CoEpiNet
Codes for _Likelihood-based Inference for Partially Observed Epidemics on Dynamic Networks_

## Things included
1. simulation codes (both SIR- and SIS-type) epidemic processes on adaptive networks
2. inference codes on complete data
3. inference codes on partially observed data w/ missing recovery times
4. codes for generating visualizations of toy examples
5. example datasets from simulations

## Things **not** included
1. real data (proprietary)
2. codes for inference on real data

## How to run things

### Demo and visualization

One of our motivations is that social intervention ("isolation") can help contain disease spread. For a particular contagion, if we enforce a simple non-pharmaceutical intervention by, say, "encouraging" people to avoid contact with their sick friends, or, asking sick people to stay home as much as possible, then usually we see a decrease in overall infection counts (peak infection cases and/or the final epidemic size). 

Run `toy_example_plots.R` to see some example epidemic processes on a small community with `20` individuals.

### Simulate complete event data and estimate parameters

The function `stochastic_coevolve_infer2` in `sim_inference.R` simulates one realization of a temporal network epidemic process (it can be coupled or decoupled) and carried out maximum likelihood estimation or Bayesian estimation for the parameters.

Another function `rep_stochastic_coevolve_infer` defined in the same file does the "simulate+infer" procedure repeatedly.

Some example commands of running simulations and complete data inference are included in `run_sim_inference.R`:

 - Simulations with N=100 people, starting with an ER(100, 0.1) network
 - Simulations with N=100 people, starting with a "hubnet" (one person is connected to all, and the rest form an ER(N-1, 0.1) network)
 - Simulations with N=500 people, starting with an ER(500, 0.1) network

An example of simulated complete event data is in `ex3dat_5.rds`.

### Inference from partially observed epidemic events

The function `infer_miss_recov` in `inference_miss_recov.R` carries out inference from event data with missing recovery times. Example datasets are `coupled_6_miss50.rds` (50% missingness) and `coupled_6_miss100.rds` (100% missingness).

To run the inference algorithm:
```r
source("./inference_miss_recov.R")
pr = data.frame(count = rep(1,4), avg = c(0.05, 0.1, 0.005, 0.05)) # or some other prior settings
miss_dats = readRDS("./coupled_6_miss50.rds") # or use "./coupled_6_miss100.rds"

inf.fil = infer_miss_recov(miss_dats, priors = pr, output.sams = 100, 
                           samples = 1000, burn = 500, thin = 2, impute = "filter")
```

And to compare our proposed data augmentation sampler with rejection sampling and Metropolis-Hastings sampling:
```r
inf.rej = infer_miss_recov(miss_dats, priors = pr, output.sams = 100, 
                           samples = 1000, burn = 500, thin = 2, impute = "reject")
inf.MH =  infer_miss_recov(miss_dats, priors = pr, output.sams = 100, 
                           samples = 1000, burn = 500, thin = 2, impute = "MH")
                           
# run diagonostics
infer.diag(inf.fil, method="filter", plot=T)
infer.diag(inf.rej, method="reject", plot=T)
infer.diag(inf.MH, method="MH",  plot=T)
```

Also, to compare the run times of the proposed algorithm and rejection sampling:
```r
source("./inference_util.R")
bp_res =
bench::press(
  ix = c(1:length(recovers)),
  {
    bench::mark(
      length(propose_recov_rej(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
                                  events = miss1$events, nei_infec = nei_infec_miss)),
      length(propose_recov_filter(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
                           events = miss1$events, nei_infec = nei_infec_miss))
    )
  }
 )
}
```

Another function `pipeline_miss_recov` in the same file goes through the following pipeline: generate missingness, conduct inference, and run diagonostics. For example, to generate 100% missing from the data in "ex3dat_5.rds" and then conduct inference using all three data augmentation samplers (proposed, rejection, MH) as well as compare them:
```r
source("./inference_miss_recov.R")
res.pipe = pipeline_miss_recov("ex3dat_5, fpath = "./", interval = 7, miss_prop = 1, miss_model = "SIR",
                                doMH = T, save_miss = F, pdfname = "pipeline_plots")
```
And all the plots will be saved to a pdf file named "pipeline_plots.pdf".
