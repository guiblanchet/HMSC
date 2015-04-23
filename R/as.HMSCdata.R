as.HMSCdata <-
function(Y,X,Tr=NULL,Random=NULL,scaleX=TRUE,scaleTr=TRUE,interceptX=TRUE,interceptTr=TRUE){
#### F. Guillaume Blanchet - July 2014
##########################################################################################
	Ynames<-colnames(Y)
	
	### Check for NAs
	if(any(is.na(Y))){
		stop("There is at least one NA in 'Y'")
	}
	if(any(is.na(X))){
		stop("There is at least one NA in 'X'")
	}
	if(!is.null(Tr)){
		if(any(is.na(Tr))){
			stop("There is at least one NA in 'Tr'")
		}
	}
	if(!is.null(Random)){
		if(any(is.na(Random))){
			stop("There is at least one NA in 'Random'")
		}
	}
	
	#### Check format
	if(length(dim(Y))!=2){
		stop("'Y' shoulds be a table")
	}
	if(length(dim(X))!=2){
		stop("'X' shoulds be a table")
	}
	if(!is.null(Tr)){
		if(length(dim(Tr))!=2){
			stop("'Tr' shoulds be a table")
		}
	}
	
	if(!is.null(Random)){
		if(is.factor(Random)){
			Random<-data.frame(random=Random)
			print("'Random' was included into a data.frame")
		}else{
			if(is.data.frame(Random)){
				if(!all(mapply(is.factor,Random))){
					stop("If 'Random' is a data.frame, it should only include factors")
				}
			}else{
				stop("'Random' should be a factor or a data.frame")
			}
		}
	}
	
	#### Check if dimensions of all tables match
	if(nrow(Y)!=nrow(X)){
		stop("'X' and 'Y' should have the same number of rows")
	}
	
	if(!is.null(Tr)){
		if(ncol(Y)!=ncol(Tr)){
			stop("'Y' and 'Tr' should have the same number of columns")
		}
	}
	if(!is.null(Random)){
		if(nrow(Random)!=nrow(Y)){
			stop("'Random' and 'Y' should have a the same number of rows")
		}
	}
	
	#### Check column names
	if(is.null(Ynames)){
		colnames(Y)<-paste("y",1:ncol(Y),sep="")
		print(paste("column names were added to 'Y'"))
	}
	
	if(is.null(colnames(X))){
		colnames(X)<-paste("x",1:ncol(X),sep="")
		print("column names were added to 'X'")
	}
	
	if(!is.null(Tr)){
		if(is.null(colnames(Tr))){
			colnames(Tr)<-paste("y",1:ncol(Tr),sep="")
			print("column names were added to 'Tr'")
		}
	}
	
	if(!is.null(colnames(Random))){
		colnames(Random)<-paste("random",1:ncol(Random),sep="")
		print("column names were added to 'Random'")
	}

	#### Check row names
	if(is.null(rownames(Y))){
		rownames(Y)<-paste("site",1:nrow(Y),sep="")
		print(paste("row names were added to 'Y'"))
	}
	
	if(is.null(rownames(X))){
		rownames(X)<-paste("site",1:nrow(X),sep="")
		print("row names were added to 'X'")
	}

	if(!is.null(Tr)){
		if(is.null(rownames(Tr))){
			rownames(Tr)<-paste("t",1:nrow(Tr),sep="")
			print("row names were added to 'Tr'")
		}
	}
	if(!is.null(Random)){
		rownames(Random)<-paste("site",1:nrow(Random),sep="")
	}
	
	### Add an intercept to X and scale
	if(interceptX){
		if(scaleX){
			### Check if any columns have a null variance
			zeroVar<-which(apply(X,2,sd)==0)
			if(length(zeroVar)!=0){
				X[,-zeroVar]<-scale(X[,-zeroVar])
				warning(paste(colnames(X)[zeroVar],"are explanatory variable(s) with a variance of 0, for this reason no intercept were added, check to make sure this is OK"))
			}else{
				X<-cbind(1,scale(X))
				colnames(X)[1]<-"Intercept"
			}
		}else{
			X<-cbind(1,X)
			colnames(X)[1]<-"Intercept"
		}
	}else{
		if(scaleX){
			### Check if any columns have a null variance
			zeroVar<-which(apply(X,2,sd)==0)
			if(length(zeroVar)!=0){
				X[,-zeroVar]<-scale(X[,-zeroVar])
				warning(paste(colnames(X)[zeroVar],"are explanatory variable(s) with a variance of 0, make sure this is OK"))
			}else{
				X<-scale(X)
			}
		}
	}
	
	### Add an intercept to Tr
	if(!is.null(Tr)){
		if(interceptTr){
			if(scaleTr){
				### Check if any columns have a null variance
				zeroVar<-which(apply(Tr,1,sd)==0)
				if(length(zeroVar)!=0){
					Tr[-zeroVar,]<-scale(Tr[-zeroVar,])
					warning(paste(rownames(Tr)[zeroVar],"are trait(s) with a variance of 0, for this reason no intercept were added, check to make sure this is OK"))
				}else{
					Tr<-rbind(1,t(scale(t(Tr))))
					rownames(Tr)[1]<-"Intercept"
				}
			}else{
				Tr<-rbind(1,Tr)
				rownames(Tr)[1]<-"Intercept"
			}
		}else{
			if(scaleTr){
				### Check if any columns have a null variance
				zeroVar<-which(apply(Tr,1,sd)==0)
				if(length(zeroVar)!=0){
					Tr[-zeroVar,]<-scale(Tr[-zeroVar,])
					warning(paste(rownames(Tr)[zeroVar],"are traits with a variance of 0, make sure this is OK"))
				}
			}
		}
	}
	
	#### Check classes
	if(!is.matrix(Y)){
		Y<-as.matrix(Y)
		print("'Y' was converted to a matrix")
	}
	
	if(!is.matrix(X)){
		X<-as.matrix(X)
		print("'X' was converted to a matrix")
	}

	if(!is.null(Tr)){
		if(!is.matrix(Tr)){
			Tr<-as.matrix(Tr)
			print("'Tr' was converted to a matrix")
		}
	}
	
	#### Return results
	if(is.null(Tr) & is.null(Random)){
		Tr<-matrix(1,ncol=ncol(Y),nrow=1)
		rownames(Tr)<-paste("t",1:nrow(Tr),sep="")
		colnames(Tr)<-paste("y",1:ncol(Y),sep="")
		res<-list(Y=Y,X=X,Tr=Tr)
		attributes(res)<-list(names=c("Y","X","Tr"))
	}
	
	if(is.null(Random) & !is.null(Tr)){
		res<-list(Y=Y,X=X,Tr=Tr)
		attributes(res)<-list(names=c("Y","X","Tr"))
	}
	
	if(!is.null(Random) & is.null(Tr)){
		Tr<-matrix(1,ncol=ncol(Y),nrow=1)
		rownames(Tr)<-paste("t",1:nrow(Tr),sep="")
		colnames(Tr)<-paste("y",1:ncol(Y),sep="")
		res<-list(Y=Y,X=X,Tr=Tr,Random=Random)
		attributes(res)<-list(names=c("Y","X","Tr","Random"))
	}

	if(!is.null(Tr) & !is.null(Random)){
		res<-list(Y=Y,X=X,Tr=Tr,Random=Random)
		attributes(res)<-list(names=c("Y","X","Tr","Random"))
	}
	
	class(res)<-"HMSCdata"
	return(res)
}
