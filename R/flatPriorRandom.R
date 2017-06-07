flatPriorRandom <-
function(shrinkOverall=NULL,shrinkSpeed=NULL,shrinkLocal=NULL,family="probit"){

	if(family=="probit"){
		### Controls the overall shrinkage level
		if(is.null(shrinkOverall)){
			shrinkOverall<-c(50,1)
		}

		### Controls how fast the shrinkage increases as the factor number increases
		if(is.null(shrinkSpeed)){
			shrinkSpeed<-c(50,1)
		}

		### Hyperparameter for the local shrinkage parameter (nu in Bhattacharya and Dunson 2011)
		if(is.null(shrinkLocal)){
			shrinkLocal<-3
		}
		print("The priors for the latent variables should be OK for probit models but not necessarily for other models, be careful")
	}

	if(family=="poisson" | family=="overPoisson"){
		### Controls the overall shrinkage level
		if(is.null(shrinkOverall)){
			shrinkOverall<-c(50,1)
		}

		### Controls how fast the shrinkage increases as the factor number increases
		if(is.null(shrinkSpeed)){
			shrinkSpeed<-c(50,1)
		}

		### Hyperparameter for the local shrinkage parameter (nu in Bhattacharya and Dunson 2011)
		if(is.null(shrinkLocal)){
			shrinkLocal<-3
		}
#		print("The priors for the latent variables should be OK for Poisson models but not necessarily for other models, be careful")
	}

	if(family=="gaussian"){
		### Controls the overall shrinkage level
		if(is.null(shrinkOverall)){
			shrinkOverall<-c(2,1)
		}

		### Controls how fast the shrinkage increases as the factor number increases
		if(is.null(shrinkSpeed)){
			shrinkSpeed<-c(2,1)
		}

		### Hyperparameter for the local shrinkage parameter (nu in Bhattacharya and Dunson 2011)
		if(is.null(shrinkLocal)){
			shrinkLocal<-3
		}
		print("The priors for the latent variables should be OK for Gaussian models but not necessarily for other models, be careful")
	}

	priors<-list(shrinkOverall=shrinkOverall,
				 shrinkSpeed=shrinkSpeed,
				 shrinkLocal=shrinkLocal)

	return(priors)
}
