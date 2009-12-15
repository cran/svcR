########################################################################
# svcR: an open source library in R 
#   for SVC computation
#
########################################################################

## History of this library 
##  svcR # 2005-2009
##   written  by Nicolas Turenne  
##                 # v1.62     beta  release      # Dec-15-09   
##                 # v1.6      beta  release      # Sep-10-09   
##                 # v1.5      beta  release      # Apr-20-09   
##                 # v1.4      beta  release      # May-20-08   
##                 # v1.3      beta  release      # Jun-20-07   
##                 # v1.2      beta  release      # Apr-24-07   
##                 # v1.0      beta  release      # Jun-26-06   
##                 # v0.9      alpha release      # Sep-26-05   
##                 # v0.0      alpha release      # Apr-27-05   
## source("D:\\rbuild\\svcR\\R\\svcR2.txt")
## load("d:\\r\\library\\svc\\svc.RData")
## save.image("d:\\r\\library\\svc\\svc.RData")

########################################################################
# Load library
########################################################################

library(quadprog)
library(ade4)
library(spdep)

########################################################################
# Main function (SVC)
########################################################################

# Usage:
#   ret = findSvcModel( iris, MetOpt="optimStoch", MetLab="gridLabeling", KernChoice="KernGaussian", Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2)
#



setGeneric("findSvcModel", function(x,...) standardGeneric("findSvcModel"))

setMethod("findSvcModel", signature(x="list"),
function (
  x		= iris,		# data frame
  MetOpt	="optimStoch",			
  MetLab	="gridLabeling",			
  KernChoice	="KernGaussian",			 
  Nu		=0.8,			
  q		=20,			# q parameter
  K		=1,			# k parameter, k nearest neigbours for grid
  G		=10,			# g parameter, grid size
  Cx		=1.,			# choice of x component to display
  Cy		=2. 			# choice of y component to display
 ) {

dataFrame=x;

mOpt=kChoice=mLab=0;
if( is.character(MetOpt) ) 
	MetOpt <- match.arg(MetOpt,c("optimStoch","optimQuad"))

if(       MetOpt == "optimStoch"	){ mOpt = 1 }
else if ( MetOpt == "optimQuad"		){ mOpt = 2 }

if( is.character(MetLab) ) 
	MetLab <- match.arg(MetLab,c("gridLabeling","mstLabeling", "knnLabeling"))

if(       MetLab == "gridLabeling"	){ mLab = 1 }
else if ( MetLab == "mstLabeling"	){ mLab = 2 }
else if ( MetLab == "knnLabeling"	){ mLab = 3 }

if( is.character(KernChoice) ) 
	KernChoice <- match.arg(KernChoice,c("KernLinear","KernGaussian", "KernGaussianDist", "KernDist"))

if(       KernChoice == "KernLinear"       ){ kChoice = 0 }
else if ( KernChoice == "KernGaussian"     ){ kChoice = 1 }
else if ( KernChoice == "KernGaussianDist" ){ kChoice = 2 }
else if ( KernChoice == "KernDist"         ){ kChoice = 3 }

ret <- new("findSvcModel");

ret@Cx		= Cx;
ret@Cy		= Cy;
ret@Nu		= Nu;
ret@dataFrame	= dataFrame;

# Global Mem Heap
GMHmax = 0;

# Removing current objects
#rm(list = ls());
#SvcEnv <- new.env(parent = globalenv()); #Local Environnment
#envir=globalenv();
#assign("SvcEnv", SvcEnv, env=globalenv(), inherits = FALSE);
#assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
#nu=get("nu", env=SvcEnv, inherits = FALSE)

# Constants
MaxIter = 20 ;		# Max number of loops
MinW = -1E+10;		# Value min of w
MaxW = +1E+10;		# Value max of w 

# Value of AroundNull 
AroundNull		= 0.05;
AroundNullVA		= 0.000005;
ret@AroundNullVA	= AroundNullVA;

# Data Structure of Criterium
lagrangeCoeff <- list(W="",Y="",A="")

#Data Matrice structure
Matrice             <- list(Sparse="",Mat="",Var="",Att="");

Knn		= 1;		# K neighbours 
ret@KNN		= Knn;		
ret@SizeGrid	= G;		# Grid Size

# Cluster radius
RadiusC <- 0;

# Delta radius
r <- 0;

# min max of matrix
MaxX = MaxY = 0;
MinX = MinY = 10000000;					

#Control file
fileName   = 0;
FileOutput = 0;						  # output connection

Precision  = 0;
Error = findSvcModel.IsError(mOpt, mLab, kChoice, Nu, q, K, G, Cx, Cy, dataFrame);

##parameters init
TimeNow    <- proc.time();			# catch time reference

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- findSvcModel.loadMat( dataFrame );		# data matrix structure loading
Alert("", "ok", "\n");

ret@Data = Matrice$Mat;

Alert("", "two-feature selection...", "\t");
NumberPCA = 2;
if( cx != 0 ) { 
		Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
}
else { 
		if( sign(min(Matrice$Mat)) < 0 ) {
			MatAdjCoa    = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]),  scannf = FALSE, nf = NumberPCA);
		}
		else
			MatAdjCoa    = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]),  scannf = FALSE, nf = NumberPCA);
		Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
		Matrice$Mat = as.matrix( Matrice$Mat );
		NumberPCA = MatAdjCoa$nf;
} #EndIf
ret@Matrice = Matrice;
Alert("", "ok", "\n");

Alert("", "min max calculation...", "\t");
nlin = nrow(Matrice$Mat);
ncol = ncol(Matrice$Mat);
if( nlin < 2 ) {
	print( " not enough data - less than 2 lines ");
	return (ret);
} #endif
pp = c();
for(i in 1:nlin ){	# we fill the full matrix
	pp = c(pp, as.vector(Matrice$Mat[i,]) );
}#endfori

MinMaxXY  = .C("MinMaxMat_C",
                as.vector(pp),
                as.integer(nlin),
                as.integer(ncol),
                as.integer(1),
		as.integer(2),
                iMinMaxXY = numeric(4))$iMinMaxXY ;
MxX = MinMaxXY[1];
MnX = MinMaxXY[2];
MxY = MinMaxXY[3];
MnY = MinMaxXY[4];

ret@MinMaxXY = MinMaxXY;

#Symmetrical case
SymMat = 0;

NbClassInData = max( Matrice$Mat[,(NumberPCA+1)] );
Alert("", "ok", "\n");

Alert("", "kernel matrix...", "\t\t");
if( kChoice == 2 || kChoice == 3 ) SymMat = 1;

MK            = kernelMatrix.compute(pp, SymMat, q, ncol, nlin, kChoice);
MatriceKernel = MK@matrixKernel;
MatriceK      = MK@matrixK;
ret@MatriceK  = MK@matrixK;

Alert("", "ok", "\n"); 

# lagrange multiplier computation
Alert("", "lagrange coefficients...", "\t");

MaxValA		=  1 / ( Nu * nrow(ret@Matrice$Mat) ); # MaxVal A[i]  valeur quasi optimale pour 2*1/N 

Model			= ModelSV.compute(mOpt, MatriceKernel, MatriceK, Nu, nlin, MaxIter, MaxValA, AroundNull,AroundNullVA);
ret@lagrangeCoeff	= Model@lagrangeCoeff;
RadiusC			= Model@RadiusC;
r			= Model@SmallR;

Alert("", "ok", "\n");
cat("\t\t\t", "radiusC ", RadiusC, "\t", "smallr ", r, "\n");					

Lab = Labelling.compute ( ret, mLab, MatriceKernel, MatriceK, pp,  Nu, G, q, ncol, nlin, RadiusC, r , kChoice, NbClassInData) ;
ret@ClassPoints	= Lab@ClassPoints;
ret@NumPoints	= Lab@NumPoints;

print("time consuming");print(proc.time() - TimeNow);	# output time consuming

invisible(gc());				# freeing memory

return( ret );
} )

########################################################################
# Main function (SVC) - Ligth usage variant
########################################################################

# Usage:
#   ret = findSvcModel.Test()
#
   
findSvcModel.Test <- function (
  MetOpt	=1,			
  MetLab	=1,			
  KernChoice	=1,			 
  Nu		=0.8,			
  q		=40,			# q parameter
  K		=1,			# k parameter, k nearest neigbours for grid
  G		=10,			# g parameter, grid size
  Cx		=0,			# choice of x component to display
  Cy		=0,			# choice of y component to display
  dataFrame	=iris			# data name
 ) {
ret <- new("findSvcModel")

ret@Cx		= Cx;
ret@Cy		= Cy;
ret@Nu		= Nu;
ret@dataFrame	= dataFrame;

# Global Mem Heap
GMHmax = 0;

# Constants
MaxIter = 20 ;		# Max number of loops
MinW = -1E+10;		# Value min of w
MaxW = +1E+10;		# Value max of w 

# Value of AroundNull 
AroundNull		= 0.05;
AroundNullVA		= 0.000005;
ret@AroundNullVA	= AroundNullVA;

# Data Structure of Criterium
lagrangeCoeff <- list(W="",Y="",A="")

#Data Matrice structure
Matrice             <- list(Sparse="",Mat="",Var="",Att="");

Knn		= 1;		# K neighbours 
ret@KNN		= Knn;		
ret@SizeGrid	= G;		# Grid Size

# Cluster radius
RadiusC <- 0;

# Delta radius
r <- 0;

# min max of matrix
MaxX = MaxY = 0;
MinX = MinY = 10000000;					

#Control file
fileName   = 0;
FileOutput = 0;						  # output connection

Precision  = 0;

Error = findSvcModel.IsError(MetOpt, MetLab, KernChoice, Nu, q, K, G, Cx, Cy, dataFrame);

##parameters init
TimeNow    <- proc.time();			# catch time reference
fileName = file.path(tempdir(), "sortie.txt");

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- findSvcModel.loadMat( dataFrame );		# data matrix structure loading
Alert("", "ok", "\n");

ret@Data =  Matrice$Mat ;
 
Alert("", "two-feature selection...", "\t");
NumberPCA = 2;
if( cx != 0 ) { 
		Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
}
else { 
		if( sign(min(Matrice$Mat)) < 0 ) {
			MatAdjCoa    = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]),  scannf = FALSE, nf = NumberPCA);
		}
		else
			MatAdjCoa    = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]),  scannf = FALSE, nf = NumberPCA);
		Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
		Matrice$Mat = as.matrix( Matrice$Mat );
		NumberPCA = MatAdjCoa$nf;
} #EndIf
ret@Matrice = Matrice;
Alert("", "ok", "\n");

Alert("", "min max calculation...", "\t");
nlin = nrow(Matrice$Mat);
ncol = ncol(Matrice$Mat);
if( nlin < 2 ) {
	print( " not enough data - less than 2 lines ");
	return (ret);
} #endif
pp = c();
for(i in 1:nlin ){	# we fill the full matrix
	pp = c(pp, as.vector(Matrice$Mat[i,]) );
}#endfori

MinMaxXY  = .C("MinMaxMat_C",
                as.vector(pp),
                as.integer(nlin),
                as.integer(ncol),
                as.integer(1),
		as.integer(2),
                iMinMaxXY = numeric(4))$iMinMaxXY ;
MxX = MinMaxXY[1];
MnX = MinMaxXY[2];
MxY = MinMaxXY[3];
MnY = MinMaxXY[4];

ret@MinMaxXY = MinMaxXY;

#Symmetrical case
SymMat = 0;

NbClassInData = max( Matrice$Mat[,(NumberPCA+1)] );
Alert("", "ok", "\n");

Alert("", "kernel matrix...", "\t\t");
if( KernChoice == 2 || KernChoice == 3 ) SymMat = 1;

MK            = kernelMatrix.compute(pp, SymMat, q, ncol, nlin, KernChoice);
MatriceKernel = MK@matrixKernel;
MatriceK      = MK@matrixK;
ret@MatriceK  = MK@matrixK;

Alert("", "ok", "\n"); 

# lagrange multiplier computation
Alert("", "lagrange coefficients...", "\t");
#assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ; 

MaxValA			=  1 / ( Nu * nrow(ret@Matrice$Mat) ); # MaxVal A[i]  valeur quasi optimale pour 2*1/N 

Model			= ModelSV.compute(MetOpt, MatriceKernel, MatriceK, Nu, nlin, MaxIter, MaxValA, AroundNull,AroundNullVA);
ret@lagrangeCoeff	= Model@lagrangeCoeff;
RadiusC			= Model@RadiusC;
r			= Model@SmallR;

Alert("", "ok", "\n");
cat("\t\t\t", "radiusC ", RadiusC, "\t", "smallr ", r, "\n");					

Lab = Labelling.compute ( ret, MetLab, MatriceKernel, MatriceK, pp,  Nu, G, q, ncol, nlin, RadiusC, r , KernChoice, NbClassInData) ;
ret@ClassPoints	= Lab@ClassPoints;
ret@NumPoints	= Lab@NumPoints;

print("time consuming");print(proc.time() - TimeNow);	# output time consuming

invisible(gc());				# freeing memory

return( ret );
} 
setMethod("findSvcModel",signature( x="missing" ), findSvcModel.Test );

########################################################################
# Main function (SVC) - Evaluation usage variant
########################################################################

# Usage:
#   findSvcModel.Eval(x)
#
   
findSvcModel.Eval <- function (
  x,
  MetOpt	=1,			
  MetLab	=1,			
  KernChoice	=1,			 
  Nu		=0.8,			
  q		=40,			# q parameter
  K		=1,			# k parameter, k nearest neigbours for grid
  G		=10,			# g parameter, grid size
  Cx		=0,			# choice of x component to display
  Cy		=0,			# choice of y component to display
  dataFrame	="iris"			# data name
 ) {
DatMat = x;
ret <- new("findSvcModel")

ret@Cx		= Cx;
ret@Cy		= Cy;
ret@Nu		= Nu;
ret@dataFrame	= dataFrame;

# Constants
MaxIter = 20 ;		# Max number of loops
MinW = -1E+10;		# Value min of w
MaxW = +1E+10;		# Value max of w 

# Value of AroundNull 
AroundNull		= 0.05;
AroundNullVA		= 0.000005;
ret@AroundNullVA	= AroundNullVA;

# Data Structure of Criterium
lagrangeCoeff <- list(W="",Y="",A="")

#Data Matrice structure
Matrice             <- list(Sparse="",Mat="",Var="",Att="");

Knn		= 1;		# K neighbours 
ret@KNN		= Knn;		
ret@SizeGrid	= G;		# Grid Size

# Cluster radius
RadiusC <- 0;

# Delta radius
r <- 0;

# min max of matrix
MaxX = MaxY = 0;
MinX = MinY = 10000000;					

Precision  = 0;

##parameters init
TimeNow    <- proc.time();			# catch time reference

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");

Matrice$Mat = DatMat
ret@Data = Matrice$Mat;
Alert("", "ok", "\n");

Alert("", "two-feature selection...", "\t");
NumberPCA = 2;
if( cx != 0 ) { 
		Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
}
else { 
		if( sign(min(Matrice$Mat)) < 0 ) {
			MatAdjCoa    = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]),  scannf = FALSE, nf = NumberPCA);
		}
		else
			MatAdjCoa    = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]),  scannf = FALSE, nf = NumberPCA);
		Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
		Matrice$Mat = as.matrix( Matrice$Mat );
		NumberPCA = MatAdjCoa$nf;
} #EndIf
ret@Matrice = Matrice;
Alert("", "ok", "\n");

Alert("", "min max calculation...", "\t");
nlin = nrow(Matrice$Mat);
ncol = ncol(Matrice$Mat);
if( nlin < 2 ) {
	print( " not enough data - less than 2 lines ");
	return (ret);
} #endif
pp = c();
for(i in 1:nlin ){	# we fill the full matrix
	pp = c(pp, as.vector(Matrice$Mat[i,]) );
}#endfori

MinMaxXY  = .C("MinMaxMat_C",
                as.vector(pp),
                as.integer(nlin),
                as.integer(ncol),
                as.integer(1),
		as.integer(2),
                iMinMaxXY = numeric(4))$iMinMaxXY ;
MxX = MinMaxXY[1];
MnX = MinMaxXY[2];
MxY = MinMaxXY[3];
MnY = MinMaxXY[4];

ret@MinMaxXY = MinMaxXY;

#Symmetrical case
SymMat = 0;

NbClassInData = max( Matrice$Mat[,(NumberPCA+1)] );
Alert("", "ok", "\n");

Alert("", "kernel matrix...", "\t\t");
if( KernChoice == 2 || KernChoice == 3 ) SymMat = 1;

MK            = kernelMatrix.compute(pp, SymMat, q, ncol, nlin, KernChoice);
MatriceKernel = MK@matrixKernel;
MatriceK      = MK@matrixK;
ret@MatriceK  = MatriceKernel;

Alert("", "ok", "\n"); 

# lagrange multiplier computation
Alert("", "lagrange coefficients...", "\t");
#assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ; 

MaxValA			=  1 / ( Nu * nrow(ret@Matrice$Mat) ); # MaxVal A[i]  valeur quasi optimale pour 2*1/N 

Model			= ModelSV.compute(MetOpt, MatriceKernel, MatriceK, Nu, nlin, MaxIter, MaxValA, AroundNull,AroundNullVA);
ret@lagrangeCoeff	= Model@lagrangeCoeff;
RadiusC			= Model@RadiusC;
r			= Model@SmallR;

Alert("", "ok", "\n");
cat("\t\t\t", "radiusC ", RadiusC, "\t", "smallr ", r, "\n");					

Lab = Labelling.compute ( ret, MetLab, MatriceKernel, MatriceK, pp,  Nu, G, q, ncol, nlin, RadiusC, r , KernChoice, NbClassInData) ;
ret@ClassPoints	= Lab@ClassPoints;
ret@NumPoints	= Lab@NumPoints;

print("time consuming");print(proc.time() - TimeNow);	# output time consuming

invisible(gc());				# freeing memory

return( ret );
} 
setMethod("findSvcModel",signature(  x="matrix" ), findSvcModel.Eval );

########################################################################
# Testing function for parameters values
#
########################################################################

# Usage:
#   findSvcModel.IsError(MetOpt=1, MetLab=1, KernChoice=1, Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2, dataFrame="iris")
#

findSvcModel.IsError <- function (
  x=1,				# method stoch (1) or quadprog (2)
  MetLab=1,			# method grid  (1) or mst      (2) or knn (3)
  KernChoice=1,			# kernel choice 0,1 or 2
  Nu=1,			     	# nu parameter
  q=1,				# q parameter
  K=1,				# k parameter, k nearest neigbours for grid
  G=1,				# g parameter, grid size
  Cx=1, 			# choice of x component to display
  Cy=1, 			# choice of y component to display
  dataFrame=NULL		# data name
 ) {
MetOpt=x;

if( MetOpt > 2 || MetOpt < 1 ) {
	print( " parameter MetOpt between 1 and 2 - try again ");
	return (1);
} #endif
if( MetLab > 3 || MetLab < 1 ) {
	print( " parameter MetLab between 1 and 3 - try again ");
	return (1);
} #endif
if( KernChoice > 3 || KernChoice < 0 ) {
	print( " parameter KernChoice between 0 and 3 - try again ");
	return (1);
} #endif
if( Nu > 1000 || Nu < 0 ) {
	print( " parameter Nu between 0 and 1000 - try again ");
	return (1);
} #endif
if( K > 10 || K < 0 ) {
	print( " parameter K between 0 and 10 - try again ");
	return (1);
} #endif
if( q > 1000000 || q < 0 ) {
	print( " parameter q between 0 and 1000000 - try again ");
	return (1);
} #endif
if( G > 150 || G < 0 ) {
	print( " parameter G between 0 and 500 - try again ");
	return (1);
} #endif
if( Cx > 100 || Cx < 0 ) {
	print( " parameter Cx between 0 and 100 - try again ");
	return (1);
} #endif
if( Cy > 100 || Cy < 0 ) {
	print( " parameter Cy between 0 and 100 - try again ");
	return (1);
} #endif
if( typeof(dataFrame) != "list" ) {
	print( " Data object is not a frame (i.e. a list) - try again ");
	return (1);
} #endif

return (0);
}
setMethod("findSvcModel",signature(x="numeric"), findSvcModel.IsError);


########################################################################
# Load Data Matrix
#
########################################################################

# Usage:
#   findSvcModel.loadMat(iris)
#

findSvcModel.loadMat <- function (
  x
    ) {
Mat	<- list(Sparse="",Mat="",Var="",Att="");

Mat$Var        = dimnames(x)[[1]]  ;
Mat$Att        = dimnames(x)[[2]]  ;
Mat$Mat        = data.matrix(x)  ;


return (Mat);
}
setMethod("findSvcModel",signature(x="character"), findSvcModel.loadMat);

########################################################################
# chargeMatrix 
# THIS FUNCTION IS OBSOLETE AND HAS BEEN REPLACED BY findSvcModel / signature(x="list")
########################################################################


########################################################################
# Load Data Matrix From HardDisk
#
########################################################################

# Usage:
#   findSvcModel.loadMatD("term", fileIn="D:\\rbuild\\test\\")
#

findSvcModel.loadMatD <- function (
  x,
  fileIn=""				# a file path of data
    ) {
dataFrame=x;
Mat	<- list(Sparse="",Mat="",Var="",Att="");

if( nchar(fileIn) < 2 ){
	Mat$Sparse     = read.table( system.file("data", "term_mat.txt", package = "svcR") , sep=" ");
	Mat$Att        = read.table( system.file("data", "term_att.txt", package = "svcR") , quote="", sep="\n" );
	Mat$Var        = read.table( system.file("data", "term_var.txt", package = "svcR") , quote="", sep="\n" );

	NRows   = max( Mat$Sparse[[1]] ); # IndiceColMax; # we initialize the full matrix
	NCols   = max( Mat$Sparse[[2]] ); # IndiceLineMax; 
	Mat$Mat = matrix(data = 0, nrow = NRows, ncol = NCols, byrow = FALSE, dimnames = NULL)

	for(i in 1:length(Mat$Sparse[,1]) ){	# we fill the full matrix
		IndiceLine = Mat$Sparse[[1]][i];
		IndiceCol  = Mat$Sparse[[2]][i];
		Val        = as.numeric( Mat$Sparse[[3]][i] );
		Mat$Mat[IndiceLine,IndiceCol] = Val;
	}#endfori

} #endif
else { 
	NomFile  = paste(fileIn, x, "_mat", ".txt", sep=""); #"D:\\R\\library\\svcR\\data\\iris_mat.txt";
	List_max = .C( "SizeMat_C",
		as.character(NomFile),
		ListMax = numeric(2) )$ListMax ; 

	MaxLin   = List_max[1];
	MaxCol   = List_max[2];

	ReadMat  = .C("LireMat_C",
                as.character(NomFile),
                as.integer(MaxLin),
                as.integer(MaxCol),
                iRetVec = numeric(MaxLin*MaxCol))$iRetVec

	Mat$Att	= read.table( paste(fileIn, dataFrame, "_att", ".txt", sep=""), quote="", sep="\n" );
	Mat$Var	= read.table( paste(fileIn, dataFrame, "_var", ".txt", sep=""), quote="", sep="\n" );

	Mat$Mat	= matrix(data = 0, nrow = MaxLin, ncol = MaxCol, byrow = FALSE, dimnames = NULL)

	for(i in 1:MaxLin ){	# we fill the full matrix
		Mat$Mat[i,] = ReadMat[(i*MaxCol-MaxCol+1):(i*MaxCol)];
	}#endfori
} # endif 

#write( "Matrice$Mat \n", file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(Mat$Mat), sep="\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));

zz = data.frame( Mat$Mat );
dimnames(zz)[[1]] <-   unlist(Mat$Var)  ;
#dimnames(zz)[[2]] <-   unlist( c(Mat$Att, "class") ) ;

return (zz);
}


########################################################################
##End of Class##
########################################################################
