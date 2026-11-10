
## Function for multinomial-based sex-specific data augmentation estimator

est.data_aug_multinom = function(ped2, Aug_Mult = 8)
{
  
  nDet = c(length(unique(ped2$mom)), length(unique(ped2$dad)))
  
  ParDet.mat = matrix(0, ncol = max(nDet), nrow = 2)
  
  for(j in seq(nDet[1]))
  {
    # ParDet.mat[1,j] = sum(ped2_nofsib$mom == ParsDetected[j], na.rm = T) +
    #   sum(ped2_nofsib$dad == ParsDetected[j], na.rm = T)
    ParDet.mat[1,j] = sum(ped2$mom == unique(ped2$mom)[j], na.rm = T)

  }
  
  for(j in seq(nDet[2]))
  {
    # ParDet.mat[1,j] = sum(ped2_nofsib$mom == ParsDetected[j], na.rm = T) +
    #   sum(ped2_nofsib$dad == ParsDetected[j], na.rm = T)
    ParDet.mat[2,j] = sum(ped2$dad == unique(ped2$dad)[j], na.rm = T)
    
  }
  
  ###########
  Aug_Mult = Aug_Mult

  Detects = ParDet.mat
  
  nAug = round(max(nDet) * Aug_Mult,0) - max(nDet)
  nDet_aug = max(nDet) + nAug
  
  Detects_aug = matrix(0, nrow = dim(Detects)[1], ncol = nDet_aug)
  
  Detects_aug[1:dim(Detects)[1], 1:dim(Detects)[2]] = Detects
  
  jags_data = list(
    nObs = rowSums(Detects),
    nDet = ncol(Detects_aug),
    Detects = Detects_aug)
  
  # Write model
  jags_file = ".\\RareEst_Multinom_BySex.txt"
  ### Don't write model for looped processing (use preexisting)
  #write_model(Rare_est, jags_file)
  
  # Set inits
  jags_inits = function(nc) {
    inits = list()
    for(c in 1:nc){
      inits[[c]] = list(
        omega = runif(length(nDet),0, 1),
        p = runif(2,0, 1),
        Pres = matrix(nrow = 2, sample(c(0,1),nDet_aug * 2, replace = T))
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
  
  jags_params = c("omega","Npar", "Npar_all", "p")
  n_params = length(jags_params) #used to autofill dataframe later
  
  MCMC.settings <- paste0("thin", jags_dims[names(jags_dims) == "nt"], "_draw", jags_dims[names(jags_dims) == "ni"], "_burn", jags_dims[names(jags_dims) == "nb"])
  
  #---------------- STEP 6: RUN JAGrS ---------------#
  post = jagsUI::jags(data = jags_data, inits = jags_inits(jags_dims["nc"]),
                      parameters.to.save = jags_params, jags_file, 
                      n.adapt = jags_dims["na"],
                      n.burnin = jags_dims["nb"], n.chains = jags_dims["nc"], 
                      n.iter = sum(jags_dims[c("ni", "nb")]), parallel = F,  verbose = T)
  
  post

  est_dataaug = post$q50$Npar_all #Mean or q50? Trying latter due to potential lopsided dist
  est_dataaug.sd = post$sd$Npar_all
  est_dataaug.p = post$mean$mu_p
  est_dataaug.omega = post$mean$omega
  est_dataaug.ndet = nDet_aug
  est_dataaug.omega = post$mean$omega
  est_dataaug.lci = post$q2.5$Npar_all
  est_dataaug.uci = post$q97.5$Npar_all
  
  post$summary
  
  
  return(data.frame(
    est_dataaug_multin = est_dataaug, est_dataaug_multin.sd = est_dataaug.sd, 
    est_dataaug_multin.lci = est_dataaug.lci, est_dataaug_multin.uci = est_dataaug.uci,
    est_dataaug_multin.rhat = post$Rhat$Npar_all))
}





