mixing<-function(x,toplot="paramX",param=c(1,1),add=FALSE,...){
#### F. Guillaume Blanchet - Avril 2013
##########################################################################################
	type<-c("paramX","paramTr","varX","outlierSp")
	toplot<-match.arg(toplot,type)
	
	if(toplot!="outlierSp"){
		if(length(param)!=2){
			stop("'param' should have a length of 2")
		}
		mix<-x$estimation[[toplot]][param[1],param[2],]
	}else{
		mix<-x$estimation[[toplot]][,param]
	}
	
	if(add){
		lines(mix,...)
	}else{
		plot(mix,type="l",...)
	}
}