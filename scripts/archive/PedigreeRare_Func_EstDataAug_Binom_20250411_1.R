
## Function for binomial-based pooled-sex data augmentation estimator

est.data_aug_binom = function(ped2, Aug_Mult = 8)
{
  AllPars = c(unique(ped2$mom), unique(ped2$dad))
  nDet = length(AllPars)
  
  ParDet.mat = matrix(0, ncol = nDet, nrow = 1)
  
  for(j in seq(nDet[1]))
  {
    # ParDet.mat[1,j] = sum(ped2_nofsib$mom == ParsDetected[j], na.rm = T) +
    #   sum(ped2_nofsib$dad == ParsDetected[j], na.rm = T)
    ParDet.mat[1,j] = sum(ped2$mom %in% AllPars[j] | ped2$dad %in% AllPars[j] , na.rm = T)
      
  }

  
  ###########
  Aug_Mult = Aug_Mult
  
  Detects = ParDet.mat
  
  nAug = round(nDet * Aug_Mult,0) - nDet
  nDet_aug = nDet + nAug
  
  Detects_aug = matrix(0, nrow = dim(Detects)[1], ncol = max(nDet_aug))
  
  Detects_aug[1:dim(Detects)[1], 1:dim(Detects)[2]] = Detects
  
  jags_data = list(
    nObs = rowSums(Detects),
    nDet = nDet_aug,
    Detects = as.vector(Detects_aug))
  
  # Write model
  jags_file = ".\\RareEst_Binom_All.txt"
  ### Don't write model for looped processing (use preexisting)
  #write_model(Rare_est, jags_file)
  
  # Set inits
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        omega = runif(length(nDet),0, 1),
        p = runif(1,0, 1),
        Pres = sample(c(0,1),sum(nDet_aug), replace = T)
        # prop_f = runif(length(nObs),0, 1)
        #IsFemale = sample(c(0,1), nDet_aug, replace = T)
      )
    }
    
    return(inits)
  }
  
  
  #------------- STEP 5: SET MCMC DIMENSIONS ---------------#
  jags_dims = c(
    ni = 2000,  # number of post-burn-in samples per chain
    na = 100, # number of adaptive samples
    nb = 500,  # number of burn-in samples
    nt = 20,     # thinning rate
    nc = 3      # number of chains
  )
  
  jags_params = c("omega", "Npar_all", "p")
  n_params = length(jags_params) #used to autofill dataframe later
  
  MCMC.settings <- paste0("thin", jags_dims[names(jags_dims) == "nt"], "_draw", jags_dims[names(jags_dims) == "ni"], "_burn", jags_dims[names(jags_dims) == "nb"])
  
  #---------------- STEP 6: RUN JAGrS ---------------#
  post = jagsUI::jags(data = jags_data, inits = jags_inits(jags_dims["nc"]),
                      parameters.to.save = jags_params, jags_file, 
                      n.adapt = jags_dims["na"],
                      n.burnin = jags_dims["nb"], n.chains = jags_dims["nc"], 
                      n.iter = sum(jags_dims[c("ni", "nb")]), parallel = F,  verbose = T)
  
  post
  
  est_dataaug = post$q50$Npar_all
  est_dataaug.sd = post$sd$Npar_all
  est_dataaug.p = post$mean$mu_p
  est_dataaug.omega = post$mean$omega
  est_dataaug.ndet = nDet_aug
  est_dataaug.omega = post$mean$omega
  est_dataaug.lci = post$q2.5$Npar_all
  est_dataaug.uci = post$q97.5$Npar_all
  
  post$summary
  
  
  return(data.frame(
    est_dataaug_binom = est_dataaug, est_dataaug_binom.sd = est_dataaug.sd, 
    est_dataaug_binom.lci = est_dataaug.lci, est_dataaug_binom.uci = est_dataaug.uci,
    est_dataaug_binom.rhat = post$Rhat$Npar_all))
}





