hmscProbit <-
function(data,param=NULL,priors=NULL,niter=12000,nburn=2000,thin=100,verbose=TRUE){
#### F. Guillaume Blanchet - July 2014
##########################################################################################
	outlierSp=FALSE ### There are conceptual issues with outlierSp, so for now it will be set to OFF until this is solved
	call<-match.call()
	
	### General checks
	if(niter < nburn){
		stop("'niter' should be equal or larger than 'burning'")
	}
	### A few basic objects
	nparamX<-ncol(data$X) # *** This assumes a linear model
	nsp<-ncol(data$Y)
	nsite<-nrow(data$Y)

	### Transform each data into a matrix
	Y<-as.matrix(data$Y)
	X<-as.matrix(data$X)
	
	if(!is.null(data$Tr)){
		Tr<-as.matrix(data$Tr)
	}else{
		Tr<-matrix(1,nrow=1,ncol=nsp)
	}
	### Some more basic objects
	nTr<-nrow(Tr)
	if(!is.null(data$Random)){
		nRandom<-ncol(data$Random) #nr
		nRandomLev<-mapply(nlevels,data$Random) #np
	}
	
	#====================================================
	### Initiate prior values if they have not been given
	#====================================================
	if(is.null(priors)){
		#-----------------------------------------------------------------------------------
		### Controls the overall shrinkage level (ad1 = a_1 in Bhattacharya and Dunson 2011)
		### c(ad1, bd1)
		### ad1 = Shape parameter of a gamma distribution
		### bd1 = Scale parameter of a gamma distribution
		###
		### * ad1 will likely need to be tweaked
		#-----------------------------------------------------------------------------------
		shrinkOverall.hypGamma<-c(10,1) # Shape parameter of a gamma distribution        # ad1, bd1

		#-----------------------------------------------------------------------------------------------------------------------
		### Controls how fast the shrinkage increases as the factor number increases (ad2 = a_2 in Bhattacharya and Dunson 2011)
		### c(ad2, bd2)
		### ad2 = Shape parameter of a gamma distribution
		### bd2 = Scale parameter of a gamma distribution
		###
		### * ad2 will likely need to be tweaked
		#-----------------------------------------------------------------------------------------------------------------------
		shrinkSpeed.hypGamma<-c(15,1) # Shape parameter of a gamma distribution        # ad2, bd2
		
		### Hyperparameter for the local shrinkage parameter (nu in Bhattacharya and Dunson 2011)
		shrinkLocal.hypDf<-3 # df
		
		### Prior for varX
		varXDf.hypInvWish<-nparamX+1 # f0
		varXScaleMat.hypInvWish<-diag(nparamX) # G0
		
		### Prior for paramTr
		paramTr.prior<-matrix(0,nrow=nparamX,ncol=nTr) # U0
		varparamTr.prior<-diag(nTr) # V0
		
		### Prior for the outlier species estimation
		outlierSp.prior<-4 # nu2
		
		priors<-list(paramTr=paramTr.prior,
					 varTr=varparamTr.prior,
					 varXDf=varXDf.hypInvWish,
					 varXScaleMat=varXScaleMat.hypInvWish,
					 shrinkOverall=shrinkOverall.hypGamma,
					 shrinkSpeed=shrinkSpeed.hypGamma,
					 shrinkLocal=shrinkLocal.hypDf,
					 outlierSp=outlierSp.prior)
	}
	
	#=================================================================
	### Initiate starting parameter values if they have not been given
	#=================================================================
	if(is.null(param)){
		#------------------
		### Initiate paramX
		#------------------
		options(warn=-1)
		paramX<-matrix(NA,nrow=nsp,ncol=nparamX)
		for(i in 1:nsp){
			paramX[i,]<-coef(glm(data$Y[,i]~-1+.,data=as.data.frame(data$X),family=binomial(link = "probit")))
		}
		options(warn=0)
		
		### Correct for very extreme values
		paramX[paramX > 3]<- 3 # * much faster than ifelse()
		paramX[paramX < -3]<--3
		
		#-------------------
		### Initiate paramTr
		#-------------------
		varTr<-solve(tcrossprod(Tr)+solve(priors$varTr)) # Vn
		paramTr<-t(varTr%*%(Tr%*%paramX+tcrossprod(solve(priors$varTr),priors$paramTr))) # Un
		
		#-----------------------
		### Initiate meansparamX
		#-----------------------
		meansparamX<-t(paramTr%*%Tr)
		
		#----------------
		### Initiate varX
		#----------------
		varX<-cov(paramX)+diag(0.1,nparamX)
		
		if(any(is.na(varX))){
			precX <- rWishart(1,nparamX+1,diag(nparamX))[,,1]
			varX <- solve(precX)
		}else{
			precX<-solve(varX)
		}
		
		param<-list(paramX=paramX,
					varX=varX,
					precX=precX,
					meansparamX=meansparamX,
					paramTr=paramTr,
					varTr=varTr)
	}
	
	### Degrees of freedom for Wishart distribution to update precX and varX (removed from the for loop)
	varXDf<-priors$varXDf+nsp
	
	#------------------
	### Initiate latent
	#------------------
	if(!is.null(data$Random)){
		### Initial number of latent variables
		nLatent<-rep(max(2,floor(log(nsp)*3)),nRandom) # k
		
		### Parameters for the adaptation when calculating the number and importance of latent variables
		adapt<-c(1,0.0005) # c(b0,b1)

		### redund[1] (prop) : Proportion of redundant elements within factor loadings
		### redund [2] (epsilon) : Proportion of redundant elements within factor loadings
		redund<-c(1,0.001) # c(prop,epsilon)
		
		### initiate latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
		latent<-vector("list",length=nRandom) # eta
		paramLatent<-vector("list",length=nRandom) # Lambda
		shrinkLocal<- vector("list",length=nRandom) # psijh
		paramShrinkGlobal<- vector("list",length=nRandom) # delta
		shrinkGlobal<- vector("list",length=nRandom) # tauh
		shrink<- vector("list",length=nRandom) # Plam
		for(i in 1:nRandom){
			latent[[i]]<-matrix(rnorm(nRandomLev[i]*nLatent[i]),nrow=nRandomLev[i],ncol=nLatent[i])
			paramLatent[[i]]<-matrix(rnorm(nsp*nLatent[i]),nrow=nsp,ncol=nLatent[i])
			shrinkLocal[[i]]<-matrix(rgamma(nsp*nLatent[i],shape=priors$shrinkLocal/2,rate=priors$shrinkLocal/2),nrow=nsp,ncol=nLatent[i])
			paramShrinkGlobal[[i]]<-c(rgamma(1,shape=priors$shrinkOverall[1],rate=priors$shrinkOverall[2]),rgamma(nLatent[i]-1,shape=priors$shrinkSpeed[1],rate=priors$shrinkSpeed[2]))
			shrinkGlobal[[i]]<-cumprod(paramShrinkGlobal[[i]])
			shrink[[i]]<-t(t(shrinkLocal[[i]])*shrinkGlobal[[i]])
		}
	}
	
	#======================
	### Initiate a latent Y
	#======================
	EstModel<-tcrossprod(X,param$paramX)
	ProbModel<-pnorm(0,EstModel,1) # cut
	unifSmpl<-matrix(runif(nsite*nsp),nrow=nsite,ncol=nsp) # u
	unifSmpl[Y==0]<-unifSmpl[Y==0]*ProbModel[Y==0]
	unifSmpl[Y==1]<-ProbModel[Y==1] + unifSmpl[Y==1]*(1-ProbModel[Y==1])
	Ylatent<-qnorm(unifSmpl,EstModel,1) # z
	
	### Correct for extreme values
	Ylatent[Ylatent > 20] <- 20
	Ylatent[Ylatent < -20] <- -20
	
	### Initiate outlier species parameter
	if(outlierSp){
		outlierSp<-rep(1,nsp)
	}
	
	#=========================================================================================
	### Construct the model using the developments proposed in Bhattacharya and Dunson (2011)
	#=========================================================================================
	#----------------
	### Store results
	#----------------
	saveIter<-seq(0,niter,by=thin)[-1]
	BurnCount<-which(saveIter <= nburn)
	nBurnValues<-length(BurnCount)
	BurnValues<-saveIter[BurnCount]
	IterCount<-which(saveIter > nburn)
	nIterValues<-length(IterCount)
	IterValues<-saveIter[IterCount]
	
	#__________
	### Burning
	#__________
	paramX.burn<-array(dim=c(nsp,nparamX,nBurnValues))
	paramTr.burn<-array(dim=c(nparamX,nTr,nBurnValues))
	varX.burn<-array(dim=c(nparamX,nparamX,nBurnValues))
	
	if(!is.null(data$Random)){
		latent.burn<-vector("list",length=nBurnValues)
		paramLatent.burn<-vector("list",length=nBurnValues)
	}
	
	if(!is.logical(outlierSp)){
		outlierSp.burn<-matrix(NA,nrow=nBurnValues,ncol=nsp)
	}
	
	burnCounter<-1
	#_____________________________
	### Estimation (after burning)
	#_____________________________
	paramX.res<-array(dim=c(nsp,nparamX,nIterValues))
	paramTr.res<-array(dim=c(nparamX,nTr,nIterValues))
	varX.res<-array(dim=c(nparamX,nparamX,nIterValues))
	
	if(!is.null(data$Random)){
		latent.res<-vector("list",length=nIterValues)
		paramLatent.res<-vector("list",length=nIterValues)
	}

	if(!is.logical(outlierSp)){
		outlierSp.res<-matrix(NA,nrow=nIterValues,ncol=nsp)
	}
	
	interCounter<-1
	#-----------------------
	### Start Gibbs sampling
	#-----------------------
	for(i in 1:niter){
		if(!is.null(data$Random)){
			#_____________________________________________________________________________
			### Update paramLatent (following the procedure proposed by Rue and Held 2005)
			#_____________________________________________________________________________
			for(j in 1:nRandom){
				Yresid<-Ylatent-tcrossprod(X,param$paramX)
				for(k in 1:nRandom){
					if(j != k){
						### Remove the effect of the other factor(s)
						Yresid<-Yresid-tcrossprod(latent[[k]][data$Random[,k],],paramLatent[[k]])
					}
				}
				
				latentCross<-crossprod(latent[[j]][data$Random[,j],])
				### Update for each species
				for(k in 1:nsp){
					shrinkLatent<-diag(shrink[[j]][k,])+latentCross
					latentSp<-crossprod(latent[[j]][data$Random[,j],],Yresid[,k])
					cholShrinkLatent<-chol(shrinkLatent) ## Llam
					meanparamLatentSp<-solve(cholShrinkLatent,solve(t(cholShrinkLatent),latentSp))
					smplparamLatentSp<-rnorm(nLatent[j])
					StructLatentSp<-solve(cholShrinkLatent,smplparamLatentSp)
					paramLatent[[j]][k,]<-StructLatentSp+meanparamLatentSp
				}
			}
			#________________
			### Update latent
			#________________
			for(j in 1:nRandom){
				Yresid<-Ylatent-tcrossprod(X,param$paramX)
				for(k in 1:nRandom){
					if(j != k){
						### Remove the effect of the other factor(s)
						Yresid<-Yresid-tcrossprod(latent[[k]][data$Random[,k],],paramLatent[[k]])
					}
				}
				
				#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
				### If the random effect is at the sampling unit level (faster updater... to check)
				#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
				if(nRandomLev[j]==nsite){
					### Precision matrix for the latent variables
					precLatent<-diag(nLatent[j])+crossprod(paramLatent[[j]])
					### Cholesky decomposition of the precision matrix
					cholprecLatent<-chol(precLatent)
					### Extract the R matrix from a QR decomposition of the Cholesky decomposition
					RprecLatent <- qr.R(qr(cholprecLatent))
					### Inverse the R matrix
					invRprecLatent <- solve(RprecLatent)
					### Calculate the variance-covariance matrix to sample the of the inverted R matrix
					varLatent<-tcrossprod(invRprecLatent)
					### Calculate the means for the new latent variables
					meansLatent<-Yresid%*%paramLatent[[j]]%*%varLatent
					### Calculate the new latent variables
					latent[[j]]<-meansLatent+matrix(rnorm(nRandomLev[j]*nLatent[j]),nrow=nRandomLev[j],ncol=nLatent[j]) %*% t(invRprecLatent)
				#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
				### If the random effect is more general
				#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
				}else{
					RandomLev<-levels(data$Random[,j])
					for(k in 1:nRandomLev[j]){
						lev<-data$Random[,j]==RandomLev[k]
						### Weight to add to consider the importance of each levels
						wRandomLev<-sum(lev)
						### Precision matrix for the latent variables
						precLatent<-diag(nLatent[j])+crossprod(paramLatent[[j]])*wRandomLev
						### Cholesky decomposition of the precision matrix
						cholprecLatent<-chol(precLatent)
						### Extract the R matrix from a QR decomposition of the Cholesky decomposition
						RprecLatent <- qr.R(qr(cholprecLatent))
						### Inverse the R matrix
						invRprecLatent <- solve(RprecLatent)
						### Calculate the variance-covariance matrix to sample the of the inverted R matrix
						varLatent<-tcrossprod(invRprecLatent)
						### Calculate the means for the new latent variables
						meansLatent<-matrix(rep(1,wRandomLev),nrow=1)%*%Yresid[lev,]%*%paramLatent[[j]]%*%varLatent
						### Calculate the new latent variables
						latent[[j]][k,]<-meansLatent+matrix(rnorm(nLatent[j]),nrow=1,ncol=nLatent[j]) %*% t(invRprecLatent)
					}
				}
			}
		}
		#__________________
		### Impute latent Y
		#__________________
		EstModel<-tcrossprod(X,param$paramX)
		if(!is.null(data$Random)){
			for(j in 1:nRandom){
				### Remove the effect of the other factor(s)
				EstModel<-EstModel+tcrossprod(latent[[j]][data$Random[,j],],paramLatent[[j]])
			}
		}
		
		ProbModel<-pnorm(0,EstModel,1) # cut
		unifSmpl<-matrix(runif(nsite*nsp),nrow=nsite,ncol=nsp) # u
		unifSmpl[Y==0]<-unifSmpl[Y==0]*ProbModel[Y==0]
		unifSmpl[Y==1]<-ProbModel[Y==1] + unifSmpl[Y==1]*(1-ProbModel[Y==1])
		Ylatent<-qnorm(unifSmpl,EstModel,1) # z
		
		### Correct for extreme values
		Ylatent[Ylatent > 20] <- 20
		Ylatent[Ylatent < -20] <- -20
		
		#________________
		### Update paramX
		#________________
		Yresid<-Ylatent
		if(!is.null(data$Random)){
			for(j in 1:nRandom){
				### Remove the effect of the other factor(s)
				Yresid<-Yresid-tcrossprod(latent[[j]][data$Random[,j],],paramLatent[[j]])
			}
		}
		### When considering outlier species
		if(!is.logical(outlierSp)){
			for(j in 1:nsp){
				varXEst<-solve(outlierSp[j]*param$precX+crossprod(X))
				meanparamXEst<-varXEst%*%(outlierSp[j]*param$precX%*%param$meansparamX[j,]+crossprod(X,Yresid[,j]))
				param$paramX[j,]<-mvrnorm(1,meanparamXEst,varXEst)
			}
		### Without considering outlier species
		}else{
			for(j in 1:nsp){
				varXEst<-solve(param$precX+crossprod(X))
				meanparamXEst<-varXEst%*%(param$precX%*%param$meansparamX[j,]+crossprod(X,Yresid[,j]))
				param$paramX[j,]<-mvrnorm(1,meanparamXEst,varXEst)
			}
		}
		#___________________
		### Update outlierSp
		#___________________
		if(!is.logical(outlierSp)){
			shapeOutlierSp<-(priors$outlierSp+nparamX)/2 ### Peut etre deplacer a l'exterieur de la boucle for()
			for(j in 1:nsp){
				paramXct<-param$paramX[j,]-param$meansparamX[j,]
				rateOutlierSp<-0.5*(priors$outlierSp+paramXct%*%param$precX%*%paramXct)
				outlierSp[j]<-rgamma(1,shape=shapeOutlierSp,rate=rateOutlierSp)
			}
		}
		#_________________________________
		### Update paramTr, varX and precX
		#_________________________________
		### Scaling of Tr and of paramX to account for outlierSp (Note : the scaling leads to scaledparamX ~ N(scaledTr %*% paramTr, V))
		if(!is.logical(outlierSp)){
			### Scaling
			outlierSpSqrt<-sqrt(outlierSp)
			scaledparamX<-param$paramX*matrix(rep(outlierSpSqrt,nparamX),nrow=nsp,ncol=nparamX)
			scaledTr<-Tr*matrix(rep(outlierSpSqrt,each=nTr),nrow=nTr,ncol=nsp)
			
			### Recalculate varTr and paramTr
			varTr<-solve(tcrossprod(scaledTr)+solve(priors$varTr))
			paramTr<-t(varTr%*%(scaledTr%*%scaledparamX+tcrossprod(solve(priors$varTr),priors$paramTr)))
			
			### Calculate the scale matrix to be used for sampling varX from the inverse Wishart distribution
			varXScaleMat<-solve(priors$varXScaleMat+crossprod(scaledparamX)+priors$paramTr%*%tcrossprod(solve(priors$varTr),priors$paramTr)-paramTr%*%tcrossprod(solve(varTr),paramTr))
		}else{
			### --- No scaling of Tr or paramX ---
			### Recalculate paramTr
			paramTr<-t(param$varTr%*%(Tr%*%param$paramX+tcrossprod(solve(priors$varTr),priors$paramTr)))
			
			### Calculate the scale matrix to be used for sampling varX from the inverse Wishart distribution
			varXScaleMat<-solve(priors$varXScaleMat+crossprod(param$paramX)+priors$paramTr%*%tcrossprod(solve(priors$varTr),priors$paramTr)-paramTr%*%tcrossprod(solve(param$varTr),paramTr))
		}

		### Update precX and varX
		param$precX<-rWishart(1,varXDf,varXScaleMat)[,,1]
		param$varX<-solve(param$precX)

		### Update paramTr
		paramTrsmpl<-matrix(rnorm(nTr*nparamX),nrow=nTr,ncol=nparamX)
		if(!is.logical(outlierSp)){
			param$paramTr<-paramTr+t(crossprod(chol(varTr),paramTrsmpl)%*%chol(param$varX))
		}else{
			param$paramTr<-paramTr+t(crossprod(chol(param$varTr),paramTrsmpl)%*%chol(param$varX))
		}
		
		### Recalculate meansparamX
		param$meansparamX<-t(paramTr%*%Tr)
		#__________________________________________________________________
		### Update shrinkLocal, paramShrinkGlobal, shrinkGlobal, and shrink
		#__________________________________________________________________
		if(!is.null(data$Random)){
			for(j in 1:nRandom){
				#_+_+_+_+_+_+_+_+_+_+_
				### Update shrinkLocal
				#_+_+_+_+_+_+_+_+_+_+_
				rateShrinkLocal<-as.vector(t(t(paramLatent[[j]]^2)*shrinkGlobal[[j]])*0.5+priors$shrinkLocal/2)
				shrinkLocal[[j]]<-matrix(rgamma(nsp*nLatent[j],shape=priors$shrinkLocal/2+0.5,rate=rateShrinkLocal),nrow=nsp,ncol=nLatent[j])
				
				#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
				### Update paramShrinkGlobal and shrinkGlobal
				#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
				paramLatentw<-shrinkLocal[[j]]*paramLatent[[j]]^2
				
				### For the overall shrinkage level
				shapeParamShrinkGlobal<-priors$shrinkOverall[1]+0.5*nsp*nLatent[j]
				rateParamShrinkGlobal<-priors$shrinkOverall[2]+0.5*(1/paramShrinkGlobal[[j]][1])*sum(shrinkGlobal[[j]]*colSums(paramLatentw)) ### *** should "(1/paramShrinkGlobal[[j]][1])" be left in the model
				paramShrinkGlobal[[j]][1]<-rgamma(1,shape=shapeParamShrinkGlobal,rate=rateParamShrinkGlobal)
				shrinkGlobal[[j]]<-cumprod(paramShrinkGlobal[[j]])
				
				### For the shrinkage speed
				for(h in 2:nLatent[j]){
					shapeParamShrinkGlobal<-priors$shrinkSpeed[1]+0.5*nsp*(nLatent[j]-h+1)
					rateParamShrinkGlobal<-priors$shrinkSpeed[2]+0.5*(1/paramShrinkGlobal[[j]][h])*sum(shrinkGlobal[[j]][h:nLatent[j]]*colSums(as.matrix(paramLatentw[,h:nLatent[j]]))) ### *** should "(1/paramShrinkGlobal[[j]][h])" be left in the model
					paramShrinkGlobal[[j]][h]<-rgamma(1,shape=shapeParamShrinkGlobal,rate=rateParamShrinkGlobal) 
					shrinkGlobal[[j]]<-cumprod(paramShrinkGlobal[[j]])
				}
				
				#_+_+_+_+_+_+_+_+
				### Update shrink
				#_+_+_+_+_+_+_+_+
				shrink[[j]]<-t(t(shrinkLocal[[j]])*shrinkGlobal[[j]])
			}
			#___________________________________________________________
			### Adaption for the number of latent variables in the model
			#___________________________________________________________
			### Probability to adapt
			probAdapt<-1/exp(adapt[1]+adapt[2]*i)
			
			for(j in 1:nRandom){
				### Measurements to evaluate whether shrinkage should ba carried out or not
				toshrinkCheck<-abs(paramLatent[[j]]) < redund[2] # tmp
				if(nsp > 1){
					toshrinkCheck<-colSums(toshrinkCheck)
				}else{
					toshrinkCheck<-as.numeric(toshrinkCheck)
				}
				toshrinkProp<-toshrinkCheck/nsp # lind
				toshrink<-toshrinkProp >= redund[1] # vec
				toshrinkSums<-sum(toshrink) # num
				
				if(runif(1) <= probAdapt){
					### Add an additional latent variable
					if(i > 20 && toshrinkSums == 0 && all(toshrinkProp < 0.995)){
						nLatent[j]<-nLatent[j]+1
						latent[[j]] <- cbind(latent[[j]],rnorm(nRandomLev[j]))
						paramLatent[[j]] <- cbind(paramLatent[[j]],0)
						shrinkLocal[[j]] <- cbind(shrinkLocal[[j]],rgamma(nsp,shape=priors$shrinkLocal/2,rate=priors$shrinkLocal/2))
						paramShrinkGlobal[[j]] <- c(paramShrinkGlobal[[j]],rgamma(1,shape=priors$shrinkSpeed[1],rate=priors$shrinkSpeed[2]))
						shrinkGlobal[[j]] <- cumprod(paramShrinkGlobal[[j]])
						shrink[[j]] <- t(t(shrinkLocal[[j]])*shrinkGlobal[[j]])
					### Remove one latent variable that is redundant with one or more of the others
					}else{
						if(toshrinkSums > 0 && nLatent[j] > 2){
							if(toshrinkSums > (nLatent[j] - 2)){
								nonRedundLatent<-1:2
							}else{
								nonRedundLatent<-setdiff(1:nLatent[j],which(toshrink))
							}
							nLatent[j]<-length(nonRedundLatent)
							latent[[j]] <- latent[[j]][,nonRedundLatent]
							paramLatent[[j]] <- paramLatent[[j]][,nonRedundLatent]
							shrinkLocal[[j]] <- shrinkLocal[[j]][,nonRedundLatent]
							paramShrinkGlobal[[j]] <- paramShrinkGlobal[[j]][nonRedundLatent]
							shrinkGlobal[[j]] <- cumprod(paramShrinkGlobal[[j]])
							shrink[[j]] <- t(t(shrinkLocal[[j]])*shrinkGlobal[[j]])
						}
					}
				}
			}
		}
		### Save sampled results
		if(any(i == saveIter)){
			if(verbose){
				print(paste(i,"iterations performed"))
			}
			if(i <= nburn){
				### Parameters of the model
				paramX.burn[,,burnCounter]<-param$paramX
				paramTr.burn[,,burnCounter]<-param$paramTr
				varX.burn[,,burnCounter]<-param$varX
				
				### Latent variables and their parameters (associated to the random effect)
				if(!is.null(data$Random)){
					latent.burn[[burnCounter]]<-latent
					paramLatent.burn[[burnCounter]]<-paramLatent
				}
				
				### Weight to estimate outlier species
				if(!is.logical(outlierSp)){
					outlierSp.burn[burnCounter,]<-outlierSp
				}
				
				burnCounter<-burnCounter+1
			}else{
				### Parameters of the model
				paramX.res[,,interCounter]<-param$paramX
				paramTr.res[,,interCounter]<-param$paramTr
				varX.res[,,interCounter]<-param$varX
				
				### Latent variables and their parameters (associated to the random effect)
				if(!is.null(data$Random)){
					latent.res[[interCounter]]<-latent
					paramLatent.res[[interCounter]]<-paramLatent
				}
				
				### Weight to estimate outlier species
				if(!is.logical(outlierSp)){
					outlierSp.res[interCounter,]<-outlierSp
				}
				interCounter<-interCounter+1
			}
		}
	}
	#-----------------------------------------------
	### Give name to each dimensions of each objects
	#-----------------------------------------------
	if(nburn > 1){
		dimnames(paramX.burn)[[3]]<-paste("burn",BurnValues,sep="")
		dimnames(paramTr.burn)[[3]]<-paste("burn",BurnValues,sep="")
		dimnames(varX.burn)[[3]]<-paste("burn",BurnValues,sep="")
	}else{
		dimnames(paramX.burn)[[3]]<-list("burn1")
		dimnames(paramTr.burn)[[3]]<-list("burn1")
		dimnames(varX.burn)[[3]]<-list("burn1")
	}
	
	dimnames(paramX.burn)[[1]]<-colnames(data$Y)
	dimnames(paramX.burn)[[2]]<-colnames(data$X)

	dimnames(paramTr.burn)[[1]]<-colnames(data$X)
	dimnames(paramTr.burn)[[2]]<-rownames(data$Tr)
	
	dimnames(varX.burn)[[1]]<-colnames(data$X)
	dimnames(varX.burn)[[2]]<-colnames(data$X)
	
	if(!is.null(data$Random)){
		names(latent.burn)<-paste("burn",BurnValues,sep="")
		names(paramLatent.burn)<-paste("burn",BurnValues,sep="")
	}
	
	if(!is.logical(outlierSp)){
		rownames(outlierSp.burn)<-paste("burn",BurnValues,sep="")
		colnames(outlierSp.burn)<-colnames(data$Y)
	}
	
	#-----------------------------------------
	### Gather burning object under one entity
	#-----------------------------------------
	if(is.null(data$Random) & is.logical(outlierSp)){
		burn<-list(paramX=paramX.burn,
				   paramTr=paramTr.burn,
				   varX=varX.burn)
	}
	
	if(!is.null(data$Random) & is.logical(outlierSp)){
		burn<-list(paramX=paramX.burn,
				   paramTr=paramTr.burn,
				   varX=varX.burn,
				   latent=latent.burn,
				   paramLatent=paramLatent.burn)
	}
	
	if(is.null(data$Random) & !is.logical(outlierSp)){
		burn<-list(paramX=paramX.burn,
				   paramTr=paramTr.burn,
				   varX=varX.burn,
				   outlierSp=outlierSp.burn)
	}
	
	if(!is.null(data$Random) & !is.logical(outlierSp)){
		burn<-list(paramX=paramX.burn,
				   paramTr=paramTr.burn,
				   varX=varX.burn,
				   latent=latent.burn,
				   paramLatent=paramLatent.burn,
				   outlierSp=outlierSp.burn)
	}
	
	#-----------------------------------------------
	### Give name to each dimensions of each objects
	#-----------------------------------------------
	dimnames(paramX.res)[[3]]<-paste("estim",IterValues,sep="")
	dimnames(paramX.res)[[1]]<-colnames(data$Y)
	dimnames(paramX.res)[[2]]<-colnames(data$X)

	dimnames(paramTr.res)[[3]]<-paste("estim",IterValues,sep="")
	dimnames(paramTr.res)[[1]]<-colnames(data$X)
	dimnames(paramTr.res)[[2]]<-rownames(data$Tr)

	dimnames(varX.res)[[3]]<-paste("estim",IterValues,sep="")
	dimnames(varX.res)[[1]]<-colnames(data$X)
	dimnames(varX.res)[[2]]<-colnames(data$X)
	
	if(!is.null(data$Random)){
		names(latent.res)<-paste("estim",IterValues,sep="")
		names(paramLatent.res)<-paste("estim",IterValues,sep="")
	}
	
	if(!is.logical(outlierSp)){
		rownames(outlierSp.res)<-paste("estim",IterValues,sep="")
		colnames(outlierSp.res)<-colnames(data$Y)
	}
	
	#-----------------------------------------
	### Gather results object under one entity
	#-----------------------------------------
	if(is.null(data$Random) & is.logical(outlierSp)){
		res<-list(paramX=paramX.res,
				   paramTr=paramTr.res,
				   varX=varX.res)
	}
	
	if(!is.null(data$Random) & is.logical(outlierSp)){
		res<-list(paramX=paramX.res,
				   paramTr=paramTr.res,
				   varX=varX.res,
				   latent=latent.res,
				   paramLatent=paramLatent.res)
	}
	
	if(is.null(data$Random) & !is.logical(outlierSp)){
		res<-list(paramX=paramX.res,
				   paramTr=paramTr.res,
				   varX=varX.res,
				   outlierSp=outlierSp.res)
	}
	
	if(!is.null(data$Random) & !is.logical(outlierSp)){
		res<-list(paramX=paramX.res,
				   paramTr=paramTr.res,
				   varX=varX.res,
				   latent=latent.res,
				   paramLatent=paramLatent.res,
				   outlierSp=outlierSp.res)
	}
	
	
	#=================
	### Output results
	#=================
	finalRes<-list(burning=burn,estimation=res)
	class(finalRes)<-"hmsc"
	return(finalRes)
}
