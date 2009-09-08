## S4 object definitions and assigment/accessor functions for the slots.
##
## created  10.04.09 nicolas turenne
## updated   

## findSvcModel Class 
setClass(  
Class = "findSvcModel", 
representation = representation(		lagrangeCoeff	= "list",
						Matrice		= "list",
						MatriceK	= "vector",
						Data		= "matrix",
						MinMaxXY	= "vector",
						MisClass	= "vector",
						dataFrame	= "list",
						fileIn		= "character",
						ClassPoints	= "vector",
						Cx		= "numeric",
						Cy		= "numeric",
						Nu		= "numeric",
						KNN		= "numeric",
						SizeGrid	= "numeric",
						AroundNullVA	= "numeric",
						NumPoints	= "vector"),

prototype=prototype (lagrangeCoeff=list(), Matrice=list(), MatriceK=matrix(), Data=matrix(), MinMaxXY=vector(), dataFrame=list(),
                      fileIn="a", ClassPoints=vector(), Cx=1, Cy=1, Nu=1, KNN=1, SizeGrid=2, AroundNullVA=1, NumPoints=vector() ),
validity=function(object){return(TRUE)}
)

## kernelMatrix Class 
setClass(
Class= "kernelMatrix", 
representation = representation(		matrixKernel	= "vector", 
						matrixK		= "matrix"),

prototype=prototype ( matrixKernel=vector(), matrixK=matrix() ),
validity=function(object){return(TRUE)}
)

## modelSV Class 
setClass(
Class= "ModelSV", 
representation = representation(		lagrangeCoeff	= "list", 
						RadiusC		= "numeric",
						SmallR		= "numeric"),

prototype=prototype ( lagrangeCoeff=list(), RadiusC=0, SmallR=0 ),
validity=function(object){return(TRUE)}
)

## Labelling Class 
setClass(
Class= "Labelling", 
representation = representation(		ClassPoints	= "vector", 
						NumPoints	= "vector"), 

prototype=prototype ( ClassPoints=vector(), NumPoints=vector() ),
validity=function(object){return(TRUE)}
)

## ClusterEval Class 
setClass(
Class= "ClusterEval", 
representation = representation(		Precision	= "numeric"),

prototype=prototype ( Precision=0 ),
validity=function(object){return(TRUE)}
)

########################################################################
# Accessor functions
########################################################################

# Usage:
#   getNumPoints(ret); getlagrangeCoeff(ret); getMatriceK"(ret); getClassPoints(ret);
#

setGeneric(name="getNumPoints",def=function(object){standardGeneric("getNumPoints")})
setMethod(f="getNumPoints",signature="findSvcModel",
definition=function(object=new("findSvcModel")){return(object@NumPoints)}
)

setGeneric(name="getlagrangeCoeff",def=function(object){standardGeneric("getlagrangeCoeff")})
setMethod(f="getlagrangeCoeff",signature="findSvcModel",
definition=function(object=new("findSvcModel")){return(object@lagrangeCoeff)}
)

setGeneric(name="getMatriceK",def=function(object){standardGeneric("getMatriceK")})
setMethod(f="getMatriceK",signature="findSvcModel",
definition=function(object=new("findSvcModel")){return(object@MatriceK)}
)

setGeneric(name="getClassPoints",def=function(object){standardGeneric("getClassPoints")})
setMethod(f="getClassPoints",signature="findSvcModel",
definition=function(object=new("findSvcModel")){return(object@ClassPoints)}
)

########################################################################
##End of Class##
########################################################################
