## Model functions.
##
## created  10.04.09 nicolas turenne
## updated   


setGeneric("ModelSV",function(x, ...) standardGeneric("ModelSV"))

ModelSV.compute <- function (  x ,
	  MatriceKernel = NULL,
	  MatriceK      = NULL,
	  Nu		= 1,
          nlin		= 1,
          MaxIter	= 2,
          MaxValA	= 2,
          AroundNull	= 0.01,
	  AroundNullVA  = 0.01  )
{ 

WYA <- list(W="",Y="",A="")

if( x == 1 ){ 
	VectorWA =   .C("CalcWcluster_C",
			as.vector(MatriceKernel), 
			as.integer(MaxIter), 
			as.integer(nlin), 
			as.double(MaxValA), 
			as.double(Nu), 
			as.integer(-100000), 
			as.integer(100000), 
			as.double(AroundNull),
			iVectorsYA = numeric(2*nlin+1) )$iVectorsYA ;
	WYA$A = VectorWA[1:nlin]; 
	#print(WVectorsYA$A);
	#print(MatriceKernel);
}
else
if( x == 2)
	WYA <- OptimQuadProgWcluster(MatriceK,Nu,MaxValA,0.0001);
Alert("", "ok", "\n");

Alert("", "radius computation...", "\t\t"); 
RadiusC =   .C("RadiusCluster",     
			as.vector(VectorWA[1:nlin]), 
			as.integer(nlin), 
			as.double(MaxValA), 
			as.double(AroundNullVA), 
			as.vector(MatriceKernel) , 
			iR = numeric(1) )$iR ;
SmallR =   .C("SmallR",     
		as.integer(nlin), 
		as.double(RadiusC), 
		as.double(Nu), 
		as.vector(VectorWA) , 
		as.vector(MatriceKernel) , 
		iResu = numeric(1) )$iResu ;

    return( new( "ModelSV", VectorWA = WYA, RadiusC = RadiusC, SmallR = SmallR ) )
}
setMethod("ModelSV",signature(x="numeric"), ModelSV.compute)



########################################################################
# Compute Criterium for Clustering
# 
########################################################################

# Usage:
#   OptimQuadProgWcluster(MatriceKern=M,nu=1,MaxValA=1,MinW=0.01)
#

setGeneric("OptimQuadProgWcluster",function(MatriceKern, ...) standardGeneric("OptimQuadProgWcluster"))

setMethod("OptimQuadProgWcluster",signature(MatriceKern="matrix"), 
function (
    MatriceKern="",		# matrix of kernel product
    nu="",			# nu
    MaxValA="",			# MaxValA
    MinW=""			# MinW
    ) {

WYA	<- list(W="",Y="",A="");
WYA$W	<- MinW;
N	<- ncol(MatriceKern)

# computation of WYA$A bound sup
BoundSup <- 1E+10;
if( nu && N )
	#BoundSup = 1000 / ( nu * N );
	BoundSup = 500*MaxValA ; # 1 / ( nu * N );
#cat("BoundSup", BoundSup, "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));

Dmat <- 2*MatriceKern;
dvec <- diag(MatriceKern);
Amat <- matrix(0, 2*N+1, N);
for( i in 1:N){ 
		Amat[1,i]      = 1;
		Amat[i+1,i]    = 1;
		Amat[i+N+1,i]  = -1;
}
bvec = as.vector( c(1, matrix(0,1,N),matrix(-BoundSup,1,N) ) );
#  min(-t(dvec).b + 1/2 t(b).Dmat.b) with the constraints t(Amat).b >= b_0.
#  min(-Diag(K).b + 1/2.(b.2K.b)) ou max(+Diag(K).b - 1/2.(b.2K.b)) with the constraints A.b >= b_0.
S <- solve.QP(Dmat,t(dvec),t(Amat),t(bvec), meq=0, factorized=TRUE);
WYA$A  <- S$solution;
WYA$W <- S$value;

#cat( "WYA$A", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(WYA$A) , file=get("FileOutput", env=SvcEnv, inherits = FALSE));
#cat( "WYA$W", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(WYA$W) , file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(WYA);
})


########################################################################
##End of Class##
########################################################################
