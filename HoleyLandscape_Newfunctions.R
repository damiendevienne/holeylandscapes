### April 2024
### D de Vienne/J Joseph
### Additional functions for answer to Dochtermann et al. (2023) on holey landscapes

### Gaussian_sim2: Function to perform Gaussian simulations Similar to Dochtermann et al. with thez following changes: 
### - The function performs a single replicate
### - The ISS matrix can either be simulated as in Dochtermann et al. (2023) OR be input directly. 
### - The function returns the ISS matrix in addition to the final G-matrix
### - The number of traits (nbtraits) can be changed when calling the function. 

Gaussian_sim<-function(nbgenerations, ISSmat=NULL, nbtraits=10) {
	if (is.null(ISSmat)) n.traits<-nbtraits
	else n.traits<-nrow(ISSmat)

	if (is.null(ISSmat)) {
		Wright.LS<-rlkjcorr(1, n.traits , eta=1 )
		ISS.var<-rep(1,n.traits) #this currently makes no difference because I'm using truncation selection
		Wright.ISS<-cor2cov(ISS.var,Wright.LS)
	}
	else {
		Wright.ISS<-ISSmat
	}
    optimum=matrix(rep(0,n.traits),ncol=1)

    #starting values
    mu.G=rep(0,n.traits)
    mu.E=rep(0,n.traits)
    Sigma.G=diag(n.traits)
    Sigma.E=diag(n.traits)
    eigen.stor<-matrix(NA,nbgenerations,n.traits)
    #fit.stor <- matrix(NA,nbgenerations,2)
    #main model
    for (i in 1:nbgenerations){
   	  # print(i)
      G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
      eigen.stor[i,]<-c(eigen(cov(G))$values)
      E<-mvrnorm(N,mu=mu.E,Sigma=Sigma.E)
      P<-h2*G+(1-h2)*E
      fit<-apply(P,1,fitness.AL.Lande1980,omega=Wright.ISS,optimum=optimum)
      #fit.stor[i,1] <- mean(fit)
      #fit.stor[i,2] <- sd(fit)
      pop<-cbind(fit,P,G,E)
      surv.pop<-trunc.select(pop,prop=prop)	      
      Sigma.G<-cov(surv.pop[,c((n.traits+2):(n.traits*2+1))])
      mu.G<-colMeans(surv.pop[,c((n.traits+2):(n.traits*2+1))])	      
    }
    COVG<-cov(G)
    return(list(covSTART=Wright.ISS, covEND=COVG))
}


Holey_sim<-function(nbsimulations, nbgenerations, nbtraits=10) { 
  n.traits<-nbtraits

  COVG<-list()

  #### IVC Model (Holey LS) P = 0.5 ####
  p=0.5 #holey landscape proportion parameter

  #Fitness landscape/ISS
  HLS.phen<-permutations(2,n.traits,c(0,1),repeats.allowed = T)
  phen.combs<-dim(HLS.phen)[1]
  HLS.t<-c(t(HLS.phen))

  for(j in 1:nbsimulations){
    print(j)
    HLS.fit<-sample(c(0,1),phen.combs,prob=c(1-p,p),replace=T)
    #starting values
    mu.G=rep(0,n.traits)
    mu.E=rep(0,n.traits)
    Sigma.G=diag(n.traits)
    Sigma.E=diag(n.traits)
    #main model
    for (i in 1:nbgenerations){
      G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
      E<-mvrnorm(N,mu=mu.E,Sigma=Sigma.E)
      P<-h2*G+(1-h2)*E
      P2<-make.binary(P)#convert P to 0 or 1
      fit<-apply(P2,1,fitness.HL,mat=HLS.t,HLS.fit=HLS.fit, n.traits=n.traits) #n.traits argument was added to function fitness.HL for ease.
      pop<-cbind(fit,P,G,E)
      surv.pop<-pop[which(pop[,1]==1),]     
      Sigma.G<-cov(surv.pop[,c((n.traits+2):(n.traits*2+1))])
      mu.G<-colMeans(surv.pop[,c((n.traits+2):(n.traits*2+1))])
    }
  COVG[[j]]<-cov(G)
  print(G)
  }
  return(COVG)
}


computel2l1<-function(mat) {
  eigenvalues<-eigen(mat)$values
  l2l1<-eigenvalues[2]/eigenvalues[1]
  return(l2l1)
}