########################################################################
# kern.R: interface for kernel computing 
#   kernel is computed in C file to speed up
#
########################################################################

## file calling this library is svcR.R
##  svcR # 2005-2008
##
##   written  by Nicolas Turenne  
##                 # v1.5     beta  release      # Apr-13-09   
##

setGeneric("kernelMatrix",function(x, ...) standardGeneric("kernelMatrix"))

kernelMatrix.compute <- function (x,
          SymMat	= 1,
          q		= 1,
          ncol		= 2,
          nlin		= 2,
          KernChoice	= 1,
          ...)
{ 

Matvec  = .C("calcKernelMatrix_C",
                as.vector(x),
                as.integer(SymMat),
                as.integer(q),
                as.integer(ncol),
		as.integer(nlin),
		as.integer(KernChoice),
                iMatriceKernel = numeric(nlin*nlin))$iMatriceKernel ;

MatTab	= matrix(data = 0, nrow = nlin, ncol = nlin, byrow = FALSE, dimnames = NULL)
for(i in 1:nlin ){	# we fill the full matrix
	MatTab[i,] = Matvec[(i*nlin-nlin+1):(i*nlin)];
}#finfori

    return( new( "kernelMatrix", matrixKernel = Matvec, matrixK = MatTab ) )
}
setMethod("kernelMatrix",signature(x="vector"), kernelMatrix.compute)


########################################################################
# Define Kernel
# THIS FUNCTION IS OBSOLETE AND HAS BEEN REPLACED BY calcKernelMatrix_C
########################################################################

# Usage:
#   Kernel(Vec1=v1, Vec2=v2, Choice=0)
#

Kernel<- function (
    Vec1="",                        # first  vector
    Vec2="",                        # second vector
    Choice=0			    # choice of the kernel
    ) {

if( Choice == 0 )
	Res = KernelLinear(V1=Vec1,V2=Vec2)
else
if( Choice == 1)
	Res = KernelGaussian(V1=Vec1,V2=Vec2)
else
if( Choice == 2)
	Res = KernelGaussianDist(V1=Vec1)

return (Res);
}

########################################################################
# Define KernelLinear
# THIS FUNCTION IS OBSOLETE AND HAS BEEN REPLACED BY MinMaxMat_C
########################################################################

# Usage:
#   KernelLinear(V1=v1, V2=v2)
#

KernelLinear<- function (
    V1,                          # first  vector
    V2                           # second vector
    ) {

Res <- 0;

Res <- as.numeric( as.vector(V1, mode="numeric")%*%as.vector(V2, mode="numeric") );

return(Res);
}

########################################################################
# Define KernelGaussian
# THIS FUNCTION IS OBSOLETE AND HAS BEEN REPLACED BY calcKernelMatrix_C
########################################################################

# Usage:
#   KernelGaussian(V1=v1, V2=v2)
#

KernelGaussian<- function (
    V1,                          # first  vector
    V2                           # second vector
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
Res	<- 0;
q	<- get("q", env=SvcEnv, inherits = FALSE);

Res <- as.vector(V1, mode="numeric") - as.vector(V2, mode="numeric") ;
Res = exp( - q * as.numeric( Res %*% Res ) );

return(Res);
}

########################################################################
# Define KernelGaussian with element computing a distance
# THIS FUNCTION IS OBSOLETE AND HAS BEEN REPLACED BY calcKernelMatrix_C
########################################################################

# Usage:
#   KernelGaussianDist(V1=v1)
#

KernelGaussianDist<- function (
    V1                          # distance element
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
Res	<- 0;
q	<- get("q", env=SvcEnv, inherits = FALSE);

V	<- as.vector(V1, mode="numeric");

Res = exp( - q * as.numeric( sqrt(V%*%V) ) );

return(Res);
}

########################################################################
# Define Radius of a vector
#
########################################################################

# Usage:
#   RadiusPoint(Vec=v, VA=wva, Mat=Mat, MatK=MatK)
#

RadiusPoint<- function (
    Vec="" ,               # vector
    VA=""  ,               # vector of Lagrange multipliers
    Mat="" ,               # matrix of data
    MatK=""              # kernel matrix
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
Dim = ncol(Mat) - 1;

Rad <- 0 ;
KC <- get("KChoice", env=SvcEnv, inherits = FALSE);

for( i in 1:nrow(Mat) ) {

	for( j in 1:nrow(Mat) ) {
		 Rad = Rad + VA[i]*VA[j]*MatK[i,j];
	} #fin for j

	Rad = Rad - 2*VA[i]*Kernel( Vec1=Mat[i,1:Dim], Vec2=Vec[1:Dim] , Choice=KC );

} #fin for i

Rad = Rad + Kernel( Vec1=Vec, Vec2=Vec , Choice=KC) ;

if( Rad > 0 ){
	Rad <- sqrt( Rad );
}
else {
	Rad <- 0;
} #finif

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(Rad);
}

########################################################################
##End of Class##
########################################################################
