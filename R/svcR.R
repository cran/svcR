########################################################################
# svcR: an open source library in R 
#   for SVC computation
#
########################################################################

## History of this library 
##  svcR # 2005-2008
##   written  by Nicolas Turenne  
##                 # v1.4     beta  release      # May-20-08   
##                 # v1.3     beta  release      # Jun-20-07   
##                 # v1.2     beta  release      # Apr-24-07   
##                 # v1.0     beta  release      # Jun-26-06   
##                 # v0.9     alpha release      # Sep-26-05   
##                 # v0.0     alpha release      # Apr-27-05   
## source("D:\\rbuild\\svcR\\R\\svcR2.txt")
## load("d:\\r\\library\\svc\\svc.RData")
## save.image("d:\\r\\library\\svc\\svc.RData")

########################################################################
# Load library
########################################################################

rm(list = ls());

library(quadprog)
library(ade4)
library(spdep)

#dyn.load("D:\\RBuild\\svcR\\src\\svcR.dll");

#envir=globalenv();

########################################################################
# Constants
########################################################################

# Global Mem Heap
GMHmax = 0;

# Max number of loops
MaxIter = 20 ;

# Value min of w 
MinW = -1E+10;

# Value max of w 
MaxW = +1E+10;

# MaxVal A[i] 
MaxValA = 0.028; # 1 / ( nu * N ); # valeur quasi optimale pour 2*1/N 

# Value of nu 
nu = 0.8;

# Value of q 
q = 5;

# Value of Cx & Cy 
cx = 1.; cy = 2.;

# Value of AroundNull 
AroundNull   = 0.05;
AroundNullVA = 0.000005;

# Data Structure of Criterium
WVectorsYA <- list(W="",Y="",A="")

#Choice of kernel
KChoice <- 0; # 0: eucli 1: radial 2: radial+dist 

#Symmetrical case
SymMat <- 0;

########################################################################
# Global Variable
########################################################################

#Data Matrice structure
Matrice             <- list(Sparse="",Mat="",Var="",Att="");
Matrice             <- 0;
ListMatrixLearnTest <- list(ListMatLearn="", ListMatTest="", ListMatLearnTest="", ListMatEval="");
ListMatrixLearnTest <- 0;

# Grid Size 
Ngrid = 5;
# K neighbours 
Knn   = 1;

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

#Local Environnment
SvcEnv = 0;

########################################################################
# Main function (SVC)
########################################################################

# Usage:
#   findModelCluster(MetOpt=1, MetLab=1, KernChoice=1, Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2, DName="iris", fileIn="D:\\R\\library\\svcR\\")
#

findModelCluster<- function (
  MetOpt="",				# method stoch (1) or quadprog (2)
  MetLab="",				# method grid  (1) or mst      (2) or knn (3)
  KernChoice="",			# kernel choice 1,2 or 3
  Nu="",			     	# nu parameter
  q="",					# q parameter
  K="",					# k parameter, k nearest neigbours for grid
  G="",					# g parameter, grid size
  Cx="", 				# choice of x component to display
  Cy="", 				# choice of y component to display
  DName="",				# data name
  fileIn=""				# a file path of data
 ) {

Error = IsError(MetOpt, MetLab, KernChoice, Nu, q, K, G, Cx, Cy, DName, fileIn);
if( Error ) return("end-of-routine");

SvcEnv <- new.env(parent = globalenv());
assign("SvcEnv", SvcEnv, env=globalenv(), inherits = FALSE);

#parameters init
assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
assign("q",     q,  env=SvcEnv, inherits = FALSE); 
assign("Ngrid", G,  env=SvcEnv, inherits = FALSE);
assign("Knn",   K,  env=SvcEnv, inherits = FALSE);
assign("KChoice",  KernChoice,  env=SvcEnv, inherits = FALSE);
TimeNow    <- proc.time();			# catch time reference
#MemBeg     <- memory.size(max = FALSE);		# catch memory reference
fileName = file.path(tempdir(), "sortie.txt");
assign( "FileOutput", file(fileName, "w"), env=SvcEnv, inherits = FALSE );		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), ncol = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), byrow = FALSE, dimnames = NULL), env=SvcEnv, inherits = FALSE);
assign("NumPoints", array( list(), (get("Ngrid", env=SvcEnv, inherits = FALSE)+1)*(get("Ngrid", env=SvcEnv, inherits = FALSE)+1)), env=SvcEnv, inherits = FALSE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=SvcEnv, inherits = FALSE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- chargeMatrix(DataName=DName, PathIn=fileIn);		# data matrix structure loading
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
Alert("", "ok", "\n");

Alert("", "min max calculation...", "\t");
nlin = nrow(Matrice$Mat);
ncol = ncol(Matrice$Mat);
pp = c();
for(i in 1:nlin ){	# we fill the full matrix
	pp = c(pp, as.vector(Matrice$Mat[i,]) );
}#finfori
SymMat = 0;

#MinMaxMat(Mat=Matrice$Mat, Cx=cx, Cy=cy);
#print(pp);

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
assign("MaxX", MxX, env=SvcEnv, inherits = FALSE ) ; 
assign("MinX", MnX, env=SvcEnv, inherits = FALSE ) ; 
assign("MaxY", MxY, env=SvcEnv, inherits = FALSE ) ; 
assign("MinY", MnY, env=SvcEnv, inherits = FALSE ) ; 

#print(MxX);print(MnX);print(MxY);print(MnY);

NbClassInData = max( Matrice$Mat[,(NumberPCA+1)] );
assign("NbClassInData", NbClassInData, env=SvcEnv, inherits = FALSE ) ;
Alert("", "ok", "\n");

Alert("", "kernel matrix...", "\t\t");
if( KernChoice == 2) SymMat = 1;
MatriceKernel  = .C("calcKernelMatrix_C",
                as.vector(pp),
                as.integer(SymMat),
                as.integer(q),
                as.integer(ncol),
		as.integer(nlin),
		as.integer(KernChoice),
                iMatriceKernel = numeric(nlin*nlin))$iMatriceKernel ;

MatriceK	= matrix(data = 0, nrow = nlin, ncol = nlin, byrow = FALSE, dimnames = NULL)
for(i in 1:nlin ){	# we fill the full matrix
	MatriceK[i,] = MatriceKernel[(i*nlin-nlin+1):(i*nlin)];
}#finfori
Alert("", "ok", "\n"); 

# lagrange multiplier computation
Alert("", "lagrange coefficients...", "\t");
assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ; 

if( MetOpt == 1 ){ 
	VectorWA =   .C("CalcWcluster_C",
			as.vector(MatriceKernel), 
			as.integer(MaxIter), 
			as.integer(nlin), 
			as.double(MaxValA), 
			as.double(nu), 
			as.integer(-100000), 
			as.integer(100000), 
			as.double(AroundNull),
			iVectorsYA = numeric(2*nlin+1) )$iVectorsYA ;
	WVectorsYA$A = VectorWA[1:nlin]; 
	#print(WVectorsYA$A);
	#print(MatriceKernel);
}
else
if( MetOpt == 2)
	WVectorsYA <- OptimQuadProgWcluster(MatriceKern=MatriceK);
Alert("", "ok", "\n");

Alert("", "radius computation...", "\t\t"); 
RadiusC =   .C("RadiusCluster",     
			as.vector(VectorWA[1:nlin]), 
			as.integer(nlin), 
			as.double(MaxValA), 
			as.double(AroundNullVA), 
			as.vector(MatriceKernel) , 
			iR = numeric(1) )$iR ;
r =   .C("SmallR",     
		as.integer(nlin), 
		as.double(RadiusC), 
		as.double(nu), 
		as.vector(VectorWA) , 
		as.vector(MatriceKernel) , 
		iResu = numeric(1) )$iResu ;
assign("RadiusC", RadiusC, env=SvcEnv, inherits = FALSE ) ; 
assign("r", r, env=SvcEnv, inherits = FALSE ) ; 
Alert("", "ok", "\n");
cat("\t\t\t", "radiusC ", RadiusC, "\t", "smallr ", r, "\n");					

if(MetLab == 1 ) {
	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", ""); 

	NumPoints =   .C("ClusterLabeling_C",     
			as.vector(pp) ,  
			as.vector(MatriceKernel) ,  
			as.integer(G), 
			as.integer(ncol), 
			as.integer(nlin), 
			as.integer(q), 
			as.double(MxX), 
			as.double(MnX), 
			as.double(MxY), 
			as.double(MnY), 
			as.double(RadiusC), 
			as.double(r), 
			as.integer(KernChoice), 
			as.vector(VectorWA) , 
			iNumPoints = numeric( G * G * 7 ) )$iNumPoints;

	NbCluster =   .C("NbCluster_C",     
			as.vector(NumPoints), 
			as.integer(G), 
			iNCluster = integer(1) )$iNCluster;

	Alert("\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 

	ClassPoints =   .C("MatchGridPoint_C",     
			as.vector(pp) , 
			as.vector(MatriceKernel) , 
			as.integer(G), 
			as.integer(ncol), 
			as.integer(nlin), 
			as.double(MxX), 
			as.double(MnX), 
			as.double(MxY), 
			as.double(MnY), 
			as.integer(NbCluster), 
			as.integer(Knn), 
			as.integer(0), 
			as.integer(1), 
			as.vector(NumPoints) , 
			iClassPoints = numeric(nlin) )$iClassPoints;
	
	Alert("\t\t", "ok", "\n");

	cat("\t\t\t", "NbCluster ", NbCluster, "\n");					

	if( NbClassInData > 0 ) {
		Alert("\t\t", "evaluation...", "");

		MisClass =   .C("Evaluation_C",     
			as.vector(pp) ,  
			as.integer(nlin), 
			as.integer(ncol), 
			as.integer(NbCluster), 
			as.vector(ClassPoints) , 
			iMisClass = numeric(nlin) )$iMisClass;

		Alert("\t\t", "ok", "\n");
		cat("\t\t\t", "Precision ", (100*(nlin-sum(MisClass[]>0))/nlin), "\n");					
	}
	Alert("\t\t\t\t" ,"ok", "\n");
}
else if( MetLab == 2 ){
	Alert("", "mst labeling/eval...", "\n"); 
	ClassPoints = MST_labelling(Mat=Matrice$Mat, MatK= MatriceK, WYA=WVectorsYA, Cx=cx, Cy=cy);
	if( NbClassInData > 0 ) {
		Alert("\t\t", "evaluation...", "");
		MisClass    = Evaluation(Mat=Matrice$Mat, NBClass=max(ClassPoints), Cx=cx, Cy=cy, ClassPoints=ClassPoints);
		Alert("\t\t", "ok", "\n");
	}
	Alert("\t\t\t\t", "ok", "\n");
}
else if( MetLab == 3 ){
	Alert("", "knn labeling/eval...", "\n");
	ClassPoints    = KNN_labelling(Mat=Matrice$Mat, MatK= MatriceK, WYA=WVectorsYA, Cx=cx, Cy=cy);
	if( NbClassInData > 0 ) {
		Alert("\t\t", "evaluation...", "");
		MisClass = Evaluation(Mat=Matrice$Mat, NBClass=max(ClassPoints), Cx=cx, Cy=cy, ClassPoints=ClassPoints);
		Alert("\t\t", "...ok", "\n");
	}
	Alert("\t\t\t\t", "ok", "\n");
}

Alert("", "export...", "\t\t\t");
ExportClusters(MatriceVar=Matrice$Var, CPoints=ClassPoints, DName=DName, pathOut=tempdir());
Alert("", "ok", "\n");

Alert("", "display...", "\t\t\t");
DisplayData(Matrice$Mat, MatriceK, WVectorsYA, 0, 1, MisClass, NumPoints, ClassPoints);	# output results
Alert("", "ok", "\n");

print("time consuming");print(proc.time() - TimeNow);	# output time consuming
#print("Max Memory");print(GMHmax);
#print("Memory At Beginning");print(MemBeg);		# output memory consuming
#print("Memory Consuming");print(GMHmax-MemBeg);		# output memory consuming
close( get("FileOutput", env=SvcEnv, inherits = FALSE) );

invisible(gc());				# freeing memory
rm( list=ls( get("SvcEnv", env=globalenv(), inherits = TRUE) ), envir=SvcEnv );
#rm(list=ls(all=TRUE))					# delete all in workspace

return("end-of-routine");
}

########################################################################
# Main function (SVC)
########################################################################

# Usage:
#   findModelCluster(MetOpt=1, MetLab=1, KernChoice=1, Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2, DName="iris", fileIn="D:\\R\\library\\svcR\\")
#

IsError<- function (
  MetOpt="",				# method stoch (1) or quadprog (2)
  MetLab="",				# method grid  (1) or mst      (2) or knn (3)
  KernChoice="",			# kernel choice 0,1 or 2
  Nu="",			     	# nu parameter
  q="",					# q parameter
  K="",					# k parameter, k nearest neigbours for grid
  G="",					# g parameter, grid size
  Cx="", 				# choice of x component to display
  Cy="", 				# choice of y component to display
  DName="",				# data name
  fileIn=""				# a file path of data
 ) {

if( MetOpt > 2 || MetOpt < 1 ) {
	print( " parameter MetOpt between 1 and 2 - try again ");
	return (1);
} #endif
if( MetLab > 3 || MetLab < 1 ) {
	print( " parameter MetLab between 1 and 3 - try again ");
	return (1);
} #endif
if( KernChoice > 2 || KernChoice < 0 ) {
	print( " parameter KernChoice between 0 and 2 - try again ");
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
if( G > 500 || G < 0 ) {
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

return (0);
}

########################################################################
# Afffiche Data
#
########################################################################

# Usage:
#   DisplayData(Mat=matrice,MatK=matriceK,WYA=wya,Cx=1,Cy=2,ListMis=MisClas)
#

DisplayData<- function (
    Mat="" ,			# data matrix
    MatK="" ,			# kernel matrix
    WYA="" ,			# vectors classs/ coefficients
    Cx="", 			# choice of x component to display
    Cy="", 			# choice of y component to display
    ListMis="",			# list of misclassified data
    NumPoints="",		# list of grid membership to clusters
    ClassPoints=""		# list of data points to clusters
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
# framing display
par(mfrow = c(2, 2));					

# plot the data matrix
#cat("matrix", "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(as.matrix(Mat), sep="\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));					

# plot Y class in the case of an svm
if( WYA$Y[1] != "" )
	plot( 1:length(WYA$Y), WYA$Y , xlab="point", ylab="Classe");

# plot W history
#TAB = get("TabW", env=SvcEnv, inherits = FALSE);
#for(Ind in 1:length(TAB)){
#	Val = get( paste("TabW[",Ind,"]", sep="") , env=SvcEnv, inherits = FALSE);
#	TAB[Ind] = Val;
#}
#plot(TAB, xlab="#iterations", ylab="W");

# plot Lagrange parameters
if( WYA$A[1] != "" )
	plot( 1:length(WYA$A), WYA$A , xlab="point", ylab="Coefficient de Lagrange");

# plots data numbers being SV, make SV in red
plot(Mat[,1],Mat[,ncol(Mat)], xlab="data points", ylab="classes in data matrix")

#SV = ListSVPoints(VA=WYA$A);			# list of support vectors
MxValA = get("MaxValA", env=SvcEnv, inherits = FALSE ) ; 

SV =   .C("ListSVPoints_C",     
		as.vector(WYA$A) ,  
		as.integer( nrow(Mat) ), 
		as.numeric( MxValA ), 
		as.numeric( AroundNullVA ), 
		iListPoints = numeric( nrow(Mat) ) )$iListPoints;

#cat("list des SV", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(t(SV), file=get("FileOutput", env=SvcEnv, inherits = FALSE)); 
for(i in 1:length(SV) ) {				
		ISV = SV[i];
		if( !is.na(ISV) )
			points(Mat[ISV,1], Mat[ISV,ncol(Mat)], pch = 24, col = "red", bg = "yellow", cex = 1)
} #fin for i

Grid(Mat=Mat, MatK= MatK, WYA=WYA$A, ListSV=SV , 1, 2, ListMis, NumPoints, ClassPoints);

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

}

########################################################################
# Load Data Matrix
#
########################################################################

# Usage:
#   chargeMatrix(DataName="iris", PathIn="D:/R/library/svcR/")
#

chargeMatrix<- function (
    DataName="",                # name of data
    PathIn=""                   # path of the data matrix
    ) {
SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
Mat              <- list(Sparse="",Mat="",Var="",Att="");

if( nchar(PathIn) < 2 ){

	Mat$Sparse     = read.table( system.file("data", "iris_mat.txt", package = "svcR") , sep=" ");
	Mat$Att        = read.table( system.file("data", "iris_att.txt", package = "svcR") , quote="", sep=" " );
	Mat$Var        = read.table( system.file("data", "iris_var.txt", package = "svcR") , quote="", sep=" " );
	
	NRows   = max( Mat$Sparse[[1]] ); # IndiceColMax; # we initialize the full matrix
	NCols   = max( Mat$Sparse[[2]] ); # IndiceLineMax; 
	Mat$Mat = matrix(data = 0, nrow = NRows, ncol = NCols, byrow = FALSE, dimnames = NULL)

	for(i in 1:length(Mat$Sparse[,1]) ){	# we fill the full matrix
		IndiceLine = Mat$Sparse[[1]][i];
		IndiceCol  = Mat$Sparse[[2]][i];
		Val        = as.numeric( Mat$Sparse[[3]][i] );
		Mat$Mat[IndiceLine,IndiceCol] = Val;
	}#finfori

	return (Mat);
} #endif

NomFile  = paste(PathIn, DataName, "_mat", ".txt", sep=""); #"D:\\R\\library\\svcR\\data\\iris_mat.txt";
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

Mat$Att	= read.table( paste(PathIn, DataName, "_att", ".txt", sep=""), quote="", sep="\n" );
Mat$Var	= read.table( paste(PathIn, DataName, "_var", ".txt", sep=""), quote="", sep="\n" );

Mat$Mat	= matrix(data = 0, nrow = MaxLin, ncol = MaxCol, byrow = FALSE, dimnames = NULL)

for(i in 1:MaxLin ){	# we fill the full matrix
	Mat$Mat[i,] = ReadMat[(i*MaxCol-MaxCol+1):(i*MaxCol)];
}#finfori

#Samp = 50 ;
#Mat$Mat = Mat$Mat[ c( c(1:Samp,(1+50):(Samp+50),(1+100):(Samp+100)) ) , ];
#Mat$Var = as.matrix(Mat$Var[c( c(1:Samp,(1+50):(Samp+50),(1+100):(Samp+100)) ),1]);

#write( "Matrice$Mat \n", file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(Mat$Mat), sep="\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
		
#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return (Mat);
}

########################################################################
# Calc Kernel Matrix
#
########################################################################

# Usage:
#   calcKernelMatrix(matrix=M)
#

calcKernelMatrix<- function (
    matrix="",                           # matrix
    KernChoice=""			# kernel choice
   ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
NCols = ncol(matrix)-1;
NRows = nrow(matrix);

MatriceKernel = matrix(data = 0, nrow = NRows, ncol = NRows, byrow = FALSE, dimnames = NULL)

#cat( "NCols", '\t', NCols, '\t', "NRows", NRows, file=get("FileOutput", env=SvcEnv, inherits = FALSE)); 

i <- 1; j <- 1;
	
for(i in 1:NRows) {
	while( j >= i && j <= NCols ) {
		MatriceKernel[i,j] <- 0;
		if( SymMat == 0 ) {
			MatriceKernel[i,j] <- as.numeric( Kernel(Vec1=matrix[i,1:NCols], Vec2=matrix[j,1:NCols], Choice=KernChoice) );
		} else {	# case of symmetrical data matrix
			MatriceKernel[i,j] <- as.numeric( Kernel(Vec1=c(matrix[i,j]), Vec2=c(matrix[i,j]), Choice=KernChoice) );
		} #endif
		MatriceKernel[j,i] <- MatriceKernel[i,j];
		j=j+1;
	} #endwhile
}#endfori

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(MatriceKernel);
}

########################################################################
# Define Kernel
#
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
#
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
#
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
#
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
# Define Criterium for Clustering
#
########################################################################

# Usage:
#   CritereWcluster(VecteurA=Va, MatrixK=M)
#

CritereWcluster<- function (
    VecteurA="" ,                   # vector of Lagrange multipliers
    MatrixK=""                       # matrix of kernel product
    ) {

# W est le critere a maximiser
W <- 0 ; SomIJ <- 0 ;
N <- ncol(MatrixK);

for( i in 1:N ) {
        SomIJ = SomIJ + VecteurA[i]*MatrixK[i,i];
	for( j in 1:ncol(MatrixK) ) {
		SomIJ = SomIJ - VecteurA[i]*VecteurA[j]*MatrixK[i,j];
	} #fin for j
} #fin for i

W <- SomIJ ;

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(W);
}

########################################################################
# Define Constraint Cluster 1
#
########################################################################

# Usage:
#   ConstraintCluster1(BoundS=bs, VecteurA=Va)
#

ConstraintCluster1<- function (
   BoundS="",                        # bound sup of Lagrange multipliers
   VecteurA=""                       # vector of Lagrange multipliers
    ) {

N <- length(VecteurA);
for( i in 1:N ) {
	if( VecteurA[i] < 0 || VecteurA[i]  > BoundS ){
		#write("beta_i \n", file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(VecteurA[i], file=get("FileOutput", env=SvcEnv, inherits = FALSE));  
		#cat("\n BoundS \t", BoundS, '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
		return( 0 );
	} #endif
} #fin for i

return( 1 );

}

########################################################################
# Define Constraint Cluster 2
#
########################################################################

# Usage:
#   ConstraintCluster2(VecteurA=Va)
#

ConstraintCluster2<- function (
    VecteurA=""                      # vector of Lagrange multipliers
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
Sum   <- 0;
Sum = sum( VecteurA );
#cat("Sum", '\n', Sum, file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#if( Sum != 1 )  #en theorie
if( (Sum >  (1 + AroundNull)) || (Sum < (1 - AroundNull)) )  #en pratique
	return( 0 );

return( 1 );

}

########################################################################
# Define Make Vector A
#
########################################################################

# Usage:
#   MakeA(Dim=N)
#

MakeA<- function (
	Dim=""			#vector dimension
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
#WVectorsYA$A <- runif(Dim);
#WVectorsYA$A = MaxValA * WVectorsYA$A ;
#WVectorsYA$A = 2*WVectorsYA$A / Dim ;

nu           = get("nu", env=SvcEnv, inherits = FALSE);
N            = Dim;
COEFFS       = runif(N);
Z            = runif(1);
for( i in 1:length(COEFFS) ){
	K = 1;
	if( COEFFS[i] <= Z ) 
		COEFFS[i] = 0;
	if( COEFFS[i] > 0.9 && K <= nu*N) {
		COEFFS[i]	=  1/(N*nu);
		K		<- K+1;
	} 
	else if( COEFFS[i] > Z ) 
		COEFFS[i] = COEFFS[i]/(nu*N);
} #finfor

i = 1;
while( sum(COEFFS) > 1 && i <= Dim ) {
	COEFFS[i]	= 0;
	i		= i+1;
} #finwhile
while( sum(COEFFS) < (1 - AroundNull) && i <= Dim ) {
	COEFFS[i]	= get("MaxValA", env=SvcEnv, inherits = FALSE);
	i		= i+1;
} #finwhile

#cat( "MaxValA", get("MaxValA", env=SvcEnv, inherits = FALSE), '\t', "Sum MaxValA", length(COEFFS)*get("MaxValA", env=SvcEnv, inherits = FALSE), '\n');
#cat( "sum(COEFFS)", sum(COEFFS), '\t', "1 + AroundNull", (1 + AroundNull), '\t', "1 - AroundNull", (1 - AroundNull), '\n' );

return(COEFFS);
}

########################################################################
# Compute Criterium for Clustering
#
########################################################################

# Usage:
#   CalcWcluster(MatriceKern=M)
#

CalcWcluster<- function (
    MatriceKern=""                          # matrix of kernel product
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
WYA	<- list(W="",Y="",A="")
W	<- MaxW; PreviousW <- MaxW ; Iter <- 1;
WYA$W	<- W;
# creation of a random vector (randomization uniforme runif) / (randomization gaussienne rnorm)
N	<- ncol(MatriceKern)
WYA$A	<- MakeA(Dim=N);

# computation of WYA$A bound sup
BoundSup	<- 1E+10;
nu		<- get("nu", env=SvcEnv, inherits = FALSE);
MaxValA		<- get("MaxValA", env=SvcEnv, inherits = FALSE);

if( nu && N )
	BoundSup = MaxValA ; # 1 / ( nu * N );
#cat("BoundSup ", BoundSup, "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));

TabWPrec = MinW;

while( (W > TabWPrec &&  Iter <= MaxIter) ||  (Iter <= MaxIter) ) {
	
	PreviousW = W;
	W         = CritereWcluster(VecteurA=WYA$A, MatrixK=MatriceKern);
	#cat("Iter=", Iter, '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));

	if( Iter > 1 ){
		ValW   = paste("TabW[",Iter-1, "]", sep=""); 
		TabWPrec <- get(ValW, env=SvcEnv, inherits = FALSE);
		}
	
	if( W > TabWPrec ){
		ValW   = paste("TabW[",Iter,"]", sep="");
		assign(ValW, W, env=SvcEnv, inherits = FALSE);
	} 
	else {  
		ValW   = paste("TabW[",Iter,"]", sep="");
		assign(ValW, TabWPrec, env=SvcEnv, inherits = FALSE);
		WYA$A      <- MakeA(Dim=N); 
		while( 
			(ConstraintCluster1(BoundS=BoundSup, VecteurA=WYA$A) == 0) 
			||
			(ConstraintCluster2(VecteurA=WYA$A) == 0) 
		) {
				WYA$A <- MakeA(Dim=N);
		} #fin while
	} #finif

	PreviousW = W;
	W         = CritereWcluster(VecteurA=WYA$A, MatrixK=MatriceKern);
	Iter      = Iter + 1;
} #finwhile

WYA$W	<- W;

#cat("A", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(t(WYA$A), file=get("FileOutput", env=SvcEnv, inherits = FALSE));
#cat("W", '\n', t(TabW[Iter-1]) , file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(t(TabW[Iter-1]) , file=get("FileOutput", env=SvcEnv, inherits = FALSE));
#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(WYA);
}

########################################################################
# Compute Criterium for Clustering
#
########################################################################

# Usage:
#   OptimQuadProgWcluster(MatriceKern=M)
#

OptimQuadProgWcluster<- function (
    MatriceKern=""                          # matrix of kernel product
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
WYA	<- list(W="",Y="",A="");
WYA$W	<- MinW;
N	<- ncol(MatriceKern)
nu	<- get("nu", env=SvcEnv, inherits = FALSE);
MaxValA	<- get("MaxValA", env=SvcEnv, inherits = FALSE);

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

for( Iter in 1:MaxIter){ 
	ValW   = paste("TabW[",Iter,"]", sep="");
	assign(ValW, S$value, env=SvcEnv, inherits = FALSE);
}

#cat( "WYA$A", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(WYA$A) , file=get("FileOutput", env=SvcEnv, inherits = FALSE));
#cat( "WYA$W", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(WYA$W) , file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(WYA);
}

########################################################################
# Define Support Vector List
#
########################################################################

# Usage:
#   ListSVPoints(VA=va)
#

ListSVPoints<- function (
    VA=""                        # vector of Lagrange multipliers
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
N		<- length(VA);
ListP		<- vector(); length(ListP) = N;
MaxValA		<- get("MaxValA", env=SvcEnv, inherits = FALSE);

Signe = 0;
for( i in 1:N ) {
	if( VA[i] > AroundNullVA && VA[i] < MaxValA ){
		ListP[i] <- i;	Signe = 1;
	}
} #fin for i

if( Signe == 0){
	ValVA = 1; IndVal = 1;
	for( i in 1:N ) {
		if( VA[i] < ValVA ){
			ValVA  <- VA[i];
			IndVal <- i
		} #endif
	} #fin for i
	ListP[IndVal] <- IndVal;
} #endif

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(ListP);
}

########################################################################
# Calcule Vecteur W for SVC
#
########################################################################

# Usage:
#   vectorWcluster(Mat=Matrice,WYA=wya);
#

vectorWcluster<- function (
    Mat="" ,                         # data matrix
    WYA=""                           # Lagrange coefficients
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
VecW <- 0;
for(i in 1:nrow(Mat) ) {
	VecW = VecW + WYA$A[i]*Mat[i,1:(ncol(Mat)-1)];
} #fin for i
#cat("VecW", '\t', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(t(VecW), file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(VecW);
}

########################################################################
# Calcule Scalar RO for SVC
#
########################################################################

# Usage:
#   scalarRO(Mat=Matrice,VecW=w,WYA=wya);
#

scalarRO<- function (
    Mat="" ,                         # data matrix
    VecW="",				# W vector
    WYA=""                           # Lagrange coefficients
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
ro <- 0;
KC <- get("KChoice", env=SvcEnv, inherits = FALSE);

for(i in 1:nrow(Mat) ) {
	ro = ro + WYA$A[i]*Kernel(Vec1=VecW, Vec2=Mat[i,1:(ncol(Mat)-1)], KC);
} #fin for i
#cat("ro", '\n', ro, file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(ro);
}

########################################################################
# Define Radius of the cluster
#
########################################################################

# Usage:
#   RadiusCluster(VA=va, MatrixK=Mat)
#

RadiusCluster<- function (
    VA="" ,                        # vector of Lagrange multipliers
    MatrixK=""                  # matrix of kernel product
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
R        <- 0 ;
IndiceSV <- 0;
N        <- length(VA);
i        <- 1;
MaxValA  <- get("MaxValA", env=SvcEnv, inherits = FALSE);
 
while( IndiceSV == 0 && i <= N ) {
        if( VA[i] > AroundNullVA && VA[i] < MaxValA )
		IndiceSV <- i ;
	i = i +1;
} #fin while i

Min = 1000;
if( IndiceSV == 0) for(i in 1:N ) {
        if( VA[i] < Min ){
		IndiceSV <- i ;
		Min = VA[i];
	} #endif
} #finfori

#cat(IndiceSV, '\n');
if( IndiceSV != 0 ) {
	R = RadiusData(IndicePoint=IndiceSV, VA=VA, MatK=MatrixK);
	#cat("indiceSV=", IndiceSV, '\t', "R=", R, '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
}
else {	
	R = 0;
	#cat("IndiceSV = 0", "\t", "AroundNullVA", AroundNullVA, "\t", "N=", N, '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
} #finif

#if( (M = memory.size()) > GMHmax ) GMHmax <- M;

return(R);
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
# Define Radius of a point
#
########################################################################

# Usage:
#   RadiusData(IndicePoint=1, VA=wva, MatK=Mat)
#

RadiusData<- function (
    IndicePoint="" ,               # indice of a vector
    VA="" ,                        # vector of Lagrange multipliers
    MatK=""                     # matrix of kernel product
    ) {

Rad      <- 0 ;
N        <- length(VA);

for( i in 1:N ) {

	for( j in 1:N ) {
		Rad = Rad + VA[i]*VA[j]*MatK[i,j];
	} #fin for j

        Rad = Rad - 2*VA[i]*MatK[i,IndicePoint];

} #fin for i

Rad = Rad + MatK[IndicePoint,IndicePoint] ;

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
# Define small variation of r around R
#
########################################################################

# Usage:
#   SmallR(VA=wva, MatK=Mat)
#

SmallR<- function (
    VA="" ,                    # vector of Lagrange multipliers
    MatK=""                    # matrix of kernel product
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
N  <- length(VA);
r  <- vector(); length(r) = N;
nu <- get("nu", env=SvcEnv, inherits = FALSE)
C  <- 1 / ( nu * N );

for( i in 1:N ) {
        
	R    = 0; RC = 0;
	r[i] = 0;
	if( VA[i] > C * 0.98 ) {
		R  = RadiusData(IndicePoint=i,VA=VA,MatK=MatK);
		RC = get("RadiusC", env=SvcEnv, inherits = FALSE);
		if( R > RC )
			r[i] =  R*R - RC*RC;
	} #finif
	#cat("RadiusC", RC, "\t", "R", R, "\t", "ri", r[i], '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
} #finfor

# calculation of the mean and the max
rmean = 0; rmax = 0; compt = 0; rmin = 0;
for( i in 1:N ) {

	if( r[i] != 0 ){
		rmean = as.numeric( rmean + r[i] );
		compt = compt + 1;
	} #finif
	if( abs(r[i]) > rmax )	
		rmax = abs(r[i]);
	if( abs(r[i]) < rmin )	
		rmin = abs(r[i]);

} #finfor
if( compt != 0) 
	rmean = rmean / compt;

resu = as.numeric( abs(rmean) );
resu = as.numeric( rmax );

return ( resu ) ;
}

########################################################################
# Compute adjacency between two points 
#
########################################################################

# Usage:
#   AdjacencyPP(Mat=Matrice, MatK=MatriceK, Vec1=v1, Vec2=v2,VA=va);
#

AdjacencyPP<- function (
    Mat="" ,                    # data matrix
    MatK="",                    # matrix of kernel product
    Vec1="" ,                   # first vector
    Vec2="" ,                   # second vector
    WYA=""                      # Lagrange coefficients
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
adj_flag = 1; # unless a point on the path exits the sphere - pair is adjacent
		
interval = 0.0;
while( interval < 1 && adj_flag ){

	z = Vec1 + interval * (Vec2 - Vec1);	    
	interval = interval + 0.3;
	R  = RadiusPoint(Vec=z, VA=WYA, Mat=Mat, MatK=MatK);	
	RC = get("RadiusC", env=SvcEnv, inherits = FALSE);
	r  = get("r", env=SvcEnv, inherits = FALSE);
	if(  (RC*RC + r) <  R*R ){
		adj_flag = 0;
		interval = 1;
	} #finif
		
} #finwhile interval

return ( adj_flag );
}

########################################################################
# Compute adjacency matrix 
#
########################################################################

# Usage:
#   Adjacency(Mat=matrice, MatK= matriceK, WYA=va);
#

Adjacency<- function (
    Mat="" ,                   # data matrix
    MatK="" ,                  # kernel matrix
    WYA=""                     # Lagrange coefficients
    ) {

N          = nrow(MatK) * 0.1 ;
AdjacencyM = matrix(data = 0, nrow = N, ncol = N); #initialize adjacency matrix

ListIndice <- list()	# we make the list of data indices
for(i in 1:N )
	ListIndice[i] = i;	

j <- 0; 
while( j < length(ListIndice) ) {

	IndList <- 1;
	j = j + 1;
	while( IndList <= length(ListIndice) ) {
	    
	  IndTarget = as.numeric(ListIndice[IndList]);

          if( IndTarget >= j )
	  if ( AdjacencyM[j,IndTarget] != 1 ) { # if adjacency already found - no point in checking again
		
		adj_flag = AdjacencyPP(Mat, MatK, Vec1=Mat[IndTarget,1:(ncol(Mat)-1)], Vec2=Mat[j,1:(ncol(Mat)-1)], WYA$A);

		if( adj_flag == 1 ) {
			AdjacencyM[IndTarget,j] = 1;
			AdjacencyM[j,IndTarget] = 1;
			#ListIndice[[IndList]] <- NULL;
			#IndList = IndList - 1;
		} #finif

	   } #finif

         IndList = IndList +1;

	} # while IndList
	
} # finwhile ListIndice

#cat( "AdjacencyM", '\n', AdjacencyM, file=get("FileOutput", env=SvcEnv, inherits = FALSE));

print("fin adjacency");

return ( AdjacencyM );
}

########################################################################
# Compute MST proximity clustering 
#
########################################################################

# Usage:
#   MST_labelling(Mat=matrice, MatK= matriceK, WYA=wya, Cx=cx, Cy=cy);
#

MST_labelling<- function (
    Mat="" ,            # data matrix
    MatK="" ,           # kernel matrix
    WYA="" ,            # Lagrange coefficients
    Cx=""  , 		# choice of x component to display
    Cy="" 		# choice of y component to display
    ) {

Alert("\t\t", "begin mst...", "");

if( sign(min(Mat)) < 0 ) {
	MatAdjCoa    = dudi.pca(as.data.frame(Mat), scan = FALSE);
}
else
	MatAdjCoa    = dudi.coa(as.data.frame(Mat), scan = FALSE);
	
MatAdjMst    = mstree(dist.dudi(MatAdjCoa), 1);
MatAdjMST01  = neig2mat(MatAdjMst);

MatAdj01     = matrix(0, ncol=nrow(Mat), nrow=nrow(Mat));
IndRowMatMST = IndColMatMST = 0;
for(IndRowMatMST in 1:nrow(MatAdjMST01) ){
	MatAdj01[IndRowMatMST,IndRowMatMST] = 1;
	for(IndColMatMST in 1:ncol(MatAdjMST01) ){ 
		i = as.numeric(IndRowMatMST);
		j = as.numeric(MatAdjMST01[i,IndColMatMST]);
		if( j != 0 ) {
			adj_flag = AdjacencyPP(Mat, MatK, Vec1=Mat[i,1:(ncol(Mat)-1)], Vec2=Mat[j,1:(ncol(Mat)-1)], WYA$A);
			if( adj_flag == 1 ){
				MatAdj01[i,j] = 1; MatAdj01[j,i] = 1;
			} #EndIf
		} #EndIf
	} #endFor
} #endFor

IndCol = IndRow = 0;
ListItemCluster = vector(); length(ListItemCluster) = nrow(MatAdj01);
for( i in 1:nrow(MatAdj01) )	ListItemCluster[i] <- 0;	

NumCluster      = 0;
while( IndRow < nrow(MatAdj01) ){

	IndRow = IndRow + 1;

	if( ListItemCluster[IndRow] == 0 ) {
		NumCluster              = NumCluster + 1;
		ListItemCluster[IndRow] = NumCluster;
	}#EndIf
	IndCol = 1;
	while( IndCol <= ncol(MatAdj01) ) {
			
			if( MatAdj01[IndRow,IndCol] == 1 && ListItemCluster[IndCol] == 0 ){
				ListItemCluster = MineLineMat(ListItemCluster, IndCol, MatAdj01);
			}
			IndCol = IndCol + 1;
	} #EndWhilel

} #EndWhilen

Alert("\t\t", "...ok", "\n");

#cat("ListItemCluster", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(ListItemCluster), file=get("FileOutput", env=SvcEnv, inherits = FALSE));

return ( ListItemCluster );
}

########################################################################
# Compute KNN proximity clustering 
#
########################################################################

# Usage:
#   KNN_labelling(Mat=matrice, MatK= matriceK, WYA=wya, Cx=cx, Cy=cy);
#

KNN_labelling<- function (
    Mat="" ,            # data matrix
    MatK="" ,           # kernel matrix
    WYA="" ,            # Lagrange coefficients
    Cx=""  , 		# choice of x component to display
    Cy="" 		# choice of y component to display
    ) {

Alert("\t\t", "begin knn...", "");

MatKNN <- knearneigh(as.matrix(Mat[,c(Cx,Cy)]), k=4);

MatAdj01  = matrix(0, ncol=nrow(Mat), nrow=nrow(Mat));
IndRowMatKNN = IndColMatKNN = 0;
for(IndRowMatKNN in 1:nrow(MatKNN$nn) ){
	MatAdj01[IndRowMatKNN,IndRowMatKNN] = 1;
	for(IndColMatKNN in 1:ncol(MatKNN$nn) ){ 
		i = as.numeric(IndRowMatKNN);
		j = as.numeric(MatKNN$nn[i,IndColMatKNN]);
		adj_flag = AdjacencyPP(Mat, MatK, Vec1=Mat[i,1:(ncol(Mat)-1)], Vec2=Mat[j,1:(ncol(Mat)-1)], WYA$A);
		if( adj_flag == 1 ){
			MatAdj01[i,j] = 1; MatAdj01[j,i] = 1;
		} #EndIf
	} #endFor
} #endFor

IndCol = IndRow = 0;
ListItemCluster = vector(); length(ListItemCluster) = nrow(MatAdj01);
for( i in 1:nrow(MatAdj01) )	ListItemCluster[i] <- 0;	

NumCluster      = 0;
while( IndRow < nrow(MatAdj01) ){

	IndRow = IndRow + 1;

	if( ListItemCluster[IndRow] == 0 ) {
		NumCluster              = NumCluster + 1;
		ListItemCluster[IndRow] = NumCluster;
	}#EndIf
	IndCol = 1;
	while( IndCol <= ncol(MatAdj01) ) {
			
			if( MatAdj01[IndRow,IndCol] == 1 && ListItemCluster[IndCol] == 0 ){
				ListItemCluster = MineLineMat(ListItemCluster, IndCol, MatAdj01);
			}
			IndCol = IndCol + 1;
	} #EndWhilel

} #EndWhilen

Alert("\t\t", "...end knn", "\n");

#cat("ListItemCluster", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(ListItemCluster), file=get("FileOutput", env=SvcEnv, inherits = FALSE));

return ( ListItemCluster );
}

MineLineMat<- function (
    ListCluster="",                 # matrix adjacency
    NumRow="",
    Mat=""
    ) {
	NumCluster = ListCluster[NumRow];
	IndCol = 1;
	while( IndCol <= ncol(Mat) ) {
			
			if( Mat[NumRow,IndCol] == 1 && ListCluster[IndCol] == 0 ){
				ListCluster = MineLineMat(ListCluster, IndCol, Mat);
			}
			else if( ListCluster[IndCol] != 0 )
				ListCluster[NumRow] = NumCluster = ListCluster[IndCol];

			IndCol = IndCol + 1;
	} #EndWhilel

return(ListCluster);
}

########################################################################
# Compute Evaluation 
#
########################################################################

# Usage:
#   Evaluation(Mat=matrice, NBClass=nb, Cx=1, Cy=2, ClassPoints=cl);
#

Evaluation<- function (
    Mat="" ,            # data matrix
    NBClass="" ,        # number of  clusters
    Cx="", 		# choice of x component to display
    Cy="", 		# choice of y component to display
    ClassPoints=""      # Value of paires point/class
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
#output points class
NumColClass = ncol(Mat);
for(i in 1:nrow(Mat) ){
 cat("N=", '\t', i, "\t", "i=", Mat[i, Cx], "\t", "j=", Mat[i, Cy], "\t", "C=", ClassPoints[i], "\t", "Cdata=", Mat[i, NumColClass], '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
} #finfori

#return(0);

ListSortedItemsByClass  = sort( Mat[, NumColClass] );

LabelClass              = vector(); length(LabelClass)     = nrow(Mat); 
CardLabelClass          = vector(); length(CardLabelClass) = nrow(Mat); 
NBClass                 = 1; 
LabelClass[NBClass]     = ListSortedItemsByClass[1];
CardLabelClass[NBClass] = 1;
if( nrow(Mat) > 1 )
for(i in 2:nrow(Mat) ){
	if( ListSortedItemsByClass[i] == ListSortedItemsByClass[i-1] ){
		CardLabelClass[NBClass] = CardLabelClass[NBClass] + 1;
	}
	else {
		NBClass		        = NBClass + 1;
		LabelClass[NBClass]     = ListSortedItemsByClass[i];
		CardLabelClass[NBClass] = 1;
	} #endif

} #finfori

#cat("labelclass", '\n', LabelClass, '\n', "Cardlabelclass", '\n' , CardLabelClass, file=get("FileOutput", env=SvcEnv, inherits = FALSE));

BestClass <- matrix(data = NA, nrow = NBClass, ncol = 2, byrow = FALSE, dimnames = NULL); 
ListVecClass = list();

for( IndClass in 1:nrow(Mat) )
	ListVecClass[[ IndClass ]]  = vector(); 
for(i in 1:nrow(Mat) ){
	L = length( ListVecClass[[ Mat[i, NumColClass] ]] );
	ListVecClass[[ Mat[i, NumColClass] ]][L+1] = ClassPoints[i];
} #finFor

#cat( "ok ok", '\n' , file=get("FileOutput", env=SvcEnv, inherits = FALSE));

for( IndClass in 1:NBClass ) {
	VecClassSorted = sort(  ListVecClass[[ as.numeric(LabelClass[IndClass]) ]]  );
	Counter = 1; BestClass[IndClass, 2] = 1; BestClass[IndClass, 1] = VecClassSorted[ 1 ];
	#cat( "VecClassSorted", '\n' , file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(VecClassSorted), file=get("FileOutput", env=SvcEnv, inherits = FALSE));
	if( length(VecClassSorted) >= 2 )
	for( IndVecSorted in 2:length(VecClassSorted) ) {
		if( VecClassSorted[ IndVecSorted ] == VecClassSorted[ IndVecSorted - 1] ) {
			Counter = Counter + 1;
			if( Counter > BestClass[IndClass, 2] ){
				BestClass[IndClass, 2] = Counter;
				BestClass[IndClass, 1] = VecClassSorted[ IndVecSorted - 1];
			} #finif
			#cat("counter", Counter, "\t", "BestClass[IndClass, 2]", BestClass[IndClass, 2], "\t", "BestClass[IndClass, 1]", BestClass[IndClass, 1], '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
		} 
		else {
			Counter = 1;
		} #finif

	} #finIndVecSorted
}#finIndClass

#cat("BestClass", '\n', BestClass, file=get("FileOutput", env=SvcEnv, inherits = FALSE));

Precision  <- 0;
for(IndClass in 1:NBClass )  if( BestClass[IndClass, 1] != 0 )
	Precision <- Precision + BestClass[IndClass, 2];
Precision <- ( Precision/nrow(Mat) ) * 100;
assign("Precision", Precision, env=SvcEnv, inherits = FALSE);
options(digits=3); cat("Precision=", Precision, "% \n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#extraction of misclassified items
ListMis = vector(); IndListMis = 1;
for( i in 1:nrow(Mat) ) {
	for( IndClass in 1:NBClass ) {
		if( Mat[i, NumColClass] == as.numeric(LabelClass[IndClass]) ){
			#cat("best", '\t', BestClass[IndClass, 1], "\t class", ClassPoints[i], '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));
			if( BestClass[IndClass, 1] != ClassPoints[i] || BestClass[IndClass, 1] == 0 ) {
				ListMis[IndListMis] = i ;
				IndListMis          = IndListMis + 1;
			} #endIf
		} #endIf
	}#endforIndVec
} #finFor

write("ListMis \n", file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(ListMis), file=get("FileOutput", env=SvcEnv, inherits = FALSE));

return(ListMis);

}

########################################################################
# Compute Grid 
#
########################################################################

# Usage:
#   Grid(Mat=matrice, MatK= matriceK, WYA=va, ListSV=SV, Cx=1,Cy=2);
#

Grid<- function (
    Mat="" ,            # data matrix
    MatK="" ,           # kernel matrix
    WYA="" ,            # Lagrange coefficients
    ListSV="",		# list of support vectors
    Cx="", 		# choice of x component to display
    Cy="", 		# choice of y component to display
    ListMis="",		# list of misclassified data
    NumPoints="",	# list of grid membership to clusters
    ClassPoints=""	# list of data points to clusters
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
N       =  get("Ngrid", env=SvcEnv, inherits = FALSE);					
ListX   <- array(0, N); # we make the list of X values	
ListY   <- array(0, N); # we make the list of Y values
ListVec <- array(0, 3); # we make random vector
NRow = nrow(Mat); NCol = ncol(Mat) ;
							 
MaxX = get("MaxX", env=SvcEnv, inherits = FALSE);
MinX = get("MinX", env=SvcEnv, inherits = FALSE);
MaxY = get("MaxY", env=SvcEnv, inherits = FALSE);
MinY = get("MinY", env=SvcEnv, inherits = FALSE);

#zoom+
plot(NA,NA,xlim=c(MinX,MaxX), ylim=c(MinY,MaxY), xlab="Xgrid", ylab="Ygrid") ; # setting up co
							
KC = get("KChoice", env=SvcEnv, inherits = FALSE);
q  = get("q", env=SvcEnv, inherits = FALSE);


RC = get("RadiusC", env=SvcEnv, inherits = FALSE);
r  = get("r", env=SvcEnv, inherits = FALSE);

# plots grid
for(i in 1:(max(NumPoints)+1) ){
	x		= NumPoints[ (i-1)*6 +3]; 
	y		= NumPoints[ (i-1)*6 +4];
	c		= NumPoints[ (i-1)*6 +6];
	ListX[i]	= ( MaxX - MinX )*( (x) / (N-1) ) + MinX ;
	ListY[i]	= ( MaxY - MinY )*( (y) / (N-1) ) + MinY ;

	if(   c != 0  ){
		points(ListX[i], ListY[i], pch = 24, col = "yellow", bg = "yellow", cex = 1);
	}
	else
		points(ListX[i], ListY[i], pch = 21, col = "blue", bg = "black", cex = 1);
} #endfori

# plots data points
for(i in 1:NRow ){

	if( ClassPoints[i] != 0 ){
		points(Mat[i,1], Mat[i,2], pch = 21, col = "red", bg = "red", cex = 1);
	}
	#else
	#	points(Mat[i,1], Mat[i,2], pch = 24, col = "green", bg = "yellow", cex = 1)
}#finfori

# plots data numbers being SV, make SV in red
for(i in 1:length(ListSV) ) {				
		ISV = ListSV[i];
		if( !is.na(ISV) )
			points(Mat[ISV,1], Mat[ISV,2], pch = 24, col = "red", bg = "yellow", cex = 1)
} #fin for i

# plots data misclassified
if( sum( ListMis[]>0 ) != 0 )
for(i in 1: sum( ListMis[]>0 ) ){
	#cat("ListMis[i]", ListMis[i], "\t", "x", Mat[ ListMis[i]+1 ,1], "\t", "y", Mat[ ListMis[i]+1 ,2],  "\n");
	if( ListMis[i] != 0 ){
		points(Mat[ ListMis[i]+1 ,1], Mat[ ListMis[i]+1 ,2], pch = 24, col = "green", bg = "green", cex = 1);
	} #endif

}#finfori

}

########################################################################
# Compute min max of matrix
#
########################################################################

# Usage:
#   MinMaxMat(Mat=matrice, Cx=1, Cy=2);
#

MinMaxMat<- function (
    Mat="", 			# data matrix
    Cx="", 			# choice of x component to display
    Cy="" 			# choice of y component to display
    ) {
				
SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);

assign("MaxX", max(Mat[, Cx]), env=SvcEnv, inherits = FALSE) ; 
assign("MaxY", max(Mat[, Cy]), env=SvcEnv, inherits = FALSE ) ; 
assign("MinX", min(Mat[, Cx]), env=SvcEnv, inherits = FALSE ) ; 
assign("MinY", min(Mat[, Cy]), env=SvcEnv, inherits = FALSE ) ; 

MaxX = get("MaxX", env=SvcEnv, inherits = FALSE);
MinX = get("MinX", env=SvcEnv, inherits = FALSE);
MaxY = get("MaxY", env=SvcEnv, inherits = FALSE);
MinY = get("MinY", env=SvcEnv, inherits = FALSE);

#cat("MaxX", '\t', MaxX, "MinX", '\t', MinX, '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); 

if( MinX == MaxX || MinY == MaxY ){
	print("impossible to compute the grid");
	return(0);
} # EndIf

}

########################################################################
# Compute labeling clusters 
#
########################################################################

# Usage:
#   ClusterLabeling(Mat=matrice, MatK= matriceK, Cx=1,Cy=2, WYA=va);
#

ClusterLabeling<- function (
    Mat="" ,			# data matrix
    MatK="" ,                   # kernel matrix
    Cx="", 			# choice of x component to display
    Cy="", 			# choice of y component to display
    WYA=""			# Lagrange coefficients
    ) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
N	= get("Ngrid", env=SvcEnv, inherits = FALSE);					
ListVec <- array(0, ncol(Mat)-1); # we make random vector
MaxX	= get("MaxX", env=SvcEnv, inherits = FALSE);
MinX	= get("MinX", env=SvcEnv, inherits = FALSE);
MaxY	= get("MaxY", env=SvcEnv, inherits = FALSE);
MinY	= get("MinY", env=SvcEnv, inherits = FALSE);

InP = 1;
for(i in 1:N ){ 
	for(j in 1:N ){ 
		x  = ( MaxX - MinX )*( (i-1) / N ) + MinX ; ListVec[Cx] = x;
		y  = ( MaxY - MinY )*( (j-1) / N ) + MinY ; ListVec[Cy] = y; 
		R  = RadiusPoint(Vec=ListVec, VA=WYA, Mat=Mat, MatK=MatK);
		RC = get("RadiusC", env=SvcEnv, inherits = FALSE);
		r  = get("r", env=SvcEnv, inherits = FALSE);
		NM_InP = paste("NumPoints[[",InP,"]]", sep="");
		assign(NM_InP, list(IndP="", IndX="", IndY="", InBall="", NumClus=""), env=SvcEnv, inherits = FALSE);
		PG  = paste("PointGrid[",i,",",j,"]", sep="");
		assign( PG , InP, env=SvcEnv, inherits = FALSE); 
		NM1 = paste("NumPoints[",InP,"][[1]][1]", sep="");
		NM2 = paste("NumPoints[",InP,"][[1]][2]", sep="");
		NM3 = paste("NumPoints[",InP,"][[1]][3]", sep="");
		NM4 = paste("NumPoints[",InP,"][[1]][4]", sep="");
		NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
		assign( NM1 , InP, env=SvcEnv, inherits = FALSE); 
		assign( NM2 , i, env=SvcEnv, inherits = FALSE); 
		assign( NM3 , j, env=SvcEnv, inherits = FALSE); 
		assign( NM5 , 0, env=SvcEnv, inherits = FALSE); 
		if(  (RC*RC + r) >=  R*R ) {
	#	if(  (RC*RC + 0) >=  R*R ) {
			assign(NM4, 1, env=SvcEnv, inherits = FALSE);
		}
		else {
			assign(NM4, 0, env=SvcEnv, inherits = FALSE);
		}
		InP = InP +1;
	}#finforj
} #finfori

ListClusters	= list(); 
IndLC		= 0;
Signal		= 0;
for(i in 1:N ){
	for(j in 1:N ){
		PG  = paste("PointGrid[",i,",",j,"]", sep="");
		InP = get( PG , env=SvcEnv, inherits = FALSE);
		NM4 = paste("NumPoints[",InP,"][[1]][4]", sep="");
		NM2 = paste("NumPoints[",InP,"][[1]][2]", sep="");
		NM3 = paste("NumPoints[",InP,"][[1]][3]", sep="");
		NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
		if( get(NM4, env=SvcEnv, inherits = FALSE) == 1) {
			Signal	= 0;
			if( IndLC == 0 ) {	# create a new cluster
				ListClusters[[1]]        = list();
				ListClusters[[1]][1]     = InP; 
				assign( NM5 , 1, env=SvcEnv, inherits = FALSE);
				IndLC = IndLC + 1;
				Signal	= 1;
			}
			else {					# look for assign to an existing cluster
				for(ii in (i-1):(i+1)){
					for(jj in (j-1):(j+1)){

						if( ii < 1 ) ii = 1;if( jj < 1 ) jj = 1;
						if( ii > N ) ii = N;if( jj > N ) jj = N ;

						if( ii == i && jj == j ){ }
						else {
							if( Signal == 0) for( IndListCluster in 1:length(ListClusters) ) {

								ActualCluster = ListClusters[[ IndListCluster ]];
								NbElemCluster = length(ActualCluster);
							
								if( Signal == 0) for( IndMembListClust in 1:NbElemCluster ) {

									ActualPoint = ActualCluster[[ IndMembListClust ]];
									NM2_cur = paste("NumPoints[",ActualPoint,"][[1]][2]", sep="");
									NM3_cur = paste("NumPoints[",ActualPoint,"][[1]][3]", sep="");

									if( get( NM2_cur , env=SvcEnv, inherits = FALSE) == ii && get( NM3_cur , env=SvcEnv, inherits = FALSE) == jj ){
										ListClusters[[ IndListCluster ]][ NbElemCluster+1 ] = InP;
										assign( NM5 , IndListCluster, env=SvcEnv, inherits = FALSE);
										Signal	= 1;
									}#endif

								} #finforIndMembListClust

							} #finforIndListCluster
						}#finifiijj
					} #endforjj
				} #endforii
			} #endif

			if( Signal == 0 ){ #on cree un nouveau cluster
				num = length(ListClusters) + 1;
				ListClusters[[ num ]]	 = list();
				ListClusters[[num]][1]   = InP;
				assign( NM5 , num, env=SvcEnv, inherits = FALSE);
			} #endif

		} #endif
		
	}#finforj
} #finfori

#sortie test
for(i in 1:N ){
	for(j in 1:N ){
		PG  = paste("PointGrid[",i,",",j,"]", sep="");
		InP = get( PG , env=SvcEnv, inherits = FALSE);
		NM3 = paste("NumPoints[",InP,"][[1]][3]", sep="");
		NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
		#cat("N=", InP, "i=", i, "\t j=", j, "\t", "NumPoints=", get(NM3, env=SvcEnv, inherits = FALSE), "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
		#cat("NumCluster=", get(NM5, env=SvcEnv, inherits = FALSE), "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
	}
}

#calcul de la matrice de connexite entre classes
NumberCluster = length(ListClusters);

ClassConnex <- matrix(data = NA, nrow = NumberCluster, ncol = NumberCluster, byrow = FALSE, dimnames = NULL); # we init the matrix of class connexity

if( NumberCluster > 1 ) {

for(i in 1:NumberCluster ) for(j in 1:NumberCluster )	ClassConnex[i,j] = 0;

for(i in 1:N ){
	for(j in 1:N ){
		PG  = paste("PointGrid[",i,",",j,"]", sep="");
		InP = get( PG , env=SvcEnv, inherits = FALSE);
		for(ii in (i-1):(i+1)){
			for(jj in (j-1):(j+1)) { 

				if( ii < 1 ) ii = 1;if( jj < 1 ) jj = 1;
				if( ii > N ) ii = N;if( jj > N ) jj = N ;
				PG  = paste("PointGrid[",ii,",",jj,"]", sep="");
				ActualPoint = get( PG , env=SvcEnv, inherits = FALSE);
				NM5_InP = paste("NumPoints[",InP,"][[1]][5]", sep="");
				NM5_Cur = paste("NumPoints[",ActualPoint,"][[1]][5]", sep="");
				if( get(NM5_InP, env=SvcEnv, inherits = FALSE) != 0 && get(NM5_Cur, env=SvcEnv, inherits = FALSE) != 0 ) { 
					#cat("ii", ii, "\t", "jj", jj, "\t", "C1=", get(NM5_InP, env=SvcEnv, inherits = FALSE), "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
					#cat("C2=", get(NM5_Cur, env=SvcEnv, inherits = FALSE), "\t", "i=", InP, "\t j=", ActualPoint, "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
					k = as.integer(get(NM5_InP, env=SvcEnv, inherits = FALSE));
					l = as.integer(get(NM5_Cur, env=SvcEnv, inherits = FALSE));
					#cat("k=", k, "\t", "l=", l, "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
					ClassConnex[ k, l ] = 1;
				}#endif
			} #endforjj
		} #endforii
	}#finforj
} #finfori

} #endifNmat

#cat( "ClassConnex", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write( t(ClassConnex), sep='\t', file=get("FileOutput", env=SvcEnv, inherits = FALSE));

#deletion of bad classes
IndClassFusIn  = 0;
IndClassFusOut = NumberCluster;
while( IndClassFusOut > 1) {
	
	#cat("IndClassFusOut", IndClassFusOut, "length(ListClusters)", length(ListClusters), '\n');

	for(indColMatConnex in 1:(IndClassFusOut-1) ) if( ClassConnex[ indColMatConnex, IndClassFusOut ] == 1 ) {
	
		ListClusters[[ indColMatConnex ]] = c( ListClusters[[ indColMatConnex ]], ListClusters[[ IndClassFusOut ]] );
		ListClusters[[ IndClassFusOut ]]  = NULL;

	} #finforindColMatConnex
	
	IndClassFusOut = IndClassFusOut - 1;
	
} #finwhile

if( length(ListClusters) > 0 )
for( indListClusters in 1:length(ListClusters) ) {	# deletion of copied clusters (list)
	ActualCluster	= ListClusters[[ indListClusters ]];
	NbElemCluster   = length(ActualCluster);
	for( IndMembListClust in 1:NbElemCluster ) {
		ActualPoint     = ActualCluster[[ IndMembListClust ]];
		NM5		= paste("NumPoints[",ActualPoint,"][[1]][5]", sep="");
		assign( NM5 , indListClusters, env=SvcEnv, inherits = FALSE); 
	} #finforIndMembListClust
} #finforindListClusters

#output clusters
for(i in 1:(N*N) ){
	NM4 = paste("NumPoints[",i,"][[1]][4]", sep="");
	NM2 = paste("NumPoints[",i,"][[1]][2]", sep="");
	NM3 = paste("NumPoints[",i,"][[1]][3]", sep="");
	NM5 = paste("NumPoints[",i,"][[1]][5]", sep="");
	#cat("N=", i, "\t", "i=", get(NM2, env=SvcEnv, inherits = FALSE), "\t", "j=", get(NM3, env=SvcEnv, inherits = FALSE), "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
	#cat("InBal=", get(NM4, env=SvcEnv, inherits = FALSE), "\t", "C=", get(NM5, env=SvcEnv, inherits = FALSE), "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
} #finfori

##### end of computation

return (NumberCluster);
}

########################################################################
# Compute Matching Points to Grid 
#
########################################################################

# Usage:
#   MatchGridPoint(Mat=Mat, NumCluster=length(ListClusters), Cx= cx, Cy=cy, Knn=1);
#

MatchGridPoint<- function (
	Mat="",				# data matrix to match
	NumCluster="",			# number of clusters
	Cx="",				# 1st coordinate
	Cy="",				# 2nd coordinate
	Knn				# number of neighbour points
) {

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
NPoint = nrow(Mat);
NG     = get("Ngrid", env=SvcEnv, inherits = FALSE);
MaxX   = get("MaxX", env=SvcEnv, inherits = FALSE);
MinX   = get("MinX", env=SvcEnv, inherits = FALSE);
MaxY   = get("MaxY", env=SvcEnv, inherits = FALSE);
MinY   = get("MinY", env=SvcEnv, inherits = FALSE);

ClassPoints      <- array(0, NPoint ); 
#cat("N ", NPoint, '\t', "NumCluster ", NumCluster, '\t', "knn ", Knn, '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE));

if( NumCluster > 0)
for(i in 1:NPoint ){
	Xi  =  1 + round(  abs( Mat[i, Cx] - MinX )*NG / ( MaxX - MinX )  ) ;
	Yi  =  1 + round(  abs( Mat[i, Cy] - MinY )*NG / ( MaxY - MinY )  ) ; 
	if( Xi > NG) Xi = NG; if(Yi > NG) Yi = NG;

	PG  = paste("PointGrid[",Xi,",",Yi,"]", sep="");
	InP = get( PG , env=SvcEnv, inherits = FALSE);
	ScoreNeighbours  <- array(0, NumCluster );

        NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
	#cat("N=", i, "\t", "Xi=", Xi, "\t", "Yi=", Yi, "\t", "c1=", Mat[i, Cx], "\t", "c2=", Mat[i, Cy], "\t", "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
	#cat("InP=", InP, "\t", "Class=", get(NM5, env=SvcEnv, inherits = FALSE), "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));

	for(ii in (Xi-1):(Xi+1)){
			for(jj in (Yi-1):(Yi+1)) {
				if( jj == Yi && ii == Xi){
					#we leave
				}
				else {
					if( ii > NG) ii = NG; if(jj > NG) jj = NG;
					if( ii <= 0) ii = 1; if(jj <= 0) jj = 1;
					PG          = paste("PointGrid[",ii,",",jj,"]", sep="");
					ActualPoint = get( PG , env=SvcEnv, inherits = FALSE);
					#cat("PG=", PG, "\t", "ActualPoint=",ActualPoint, "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
					NM_Cur      = paste("NumPoints[",ActualPoint,"][[1]][5]", sep="");
					if( get(NM_Cur, env=SvcEnv, inherits = FALSE) != 0 ){
						NumClass                  = as.integer(get(NM_Cur, env=SvcEnv, inherits = FALSE));
						ScoreNeighbours[NumClass] <- ScoreNeighbours[NumClass] + 1 ;
						#cat("N=", i, "\t", "Class=", NumClass, "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
						#cat("ScoreNeighbours=", ScoreNeighbours[ NumClass ], "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
					}#endif
				}#endif
			} #endforjj
	} #endforii

	for(IndNumClass in 1:NumCluster ){
		#cat("N=", IndNumClass, "\t", "length(ListClusters)=", NumCluster, "\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
		#cat("i=", i, "\t", "ScoreNeighbours=", ScoreNeighbours[ IndNumClass ], "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
		if( ScoreNeighbours[ IndNumClass ] >= Knn )
			ClassPoints[i] = IndNumClass;	
	} #endforIndNumClass
	
	rm(ScoreNeighbours);
} #finfori

return( ClassPoints );
}

########################################################################
# Compute Data Sampling 
#
########################################################################

# Usage:
#   DataSampling(Mat=matrice);
#

DataSampling<- function (
    Mat="" 			# data matrix
    ) {
ListMatrixLearnTest                  <- list();
ListMatrixLearnTest$ListMatLearn     <- list(); 
ListMatrixLearnTest$ListMatTest      <- list();
ListMatrixLearnTest$ListMatLearnTest <- list();
ListMatrixLearnTest$ListMatEval      <- list();

ListMatrixLearnTest$ListMatLearn[[1]]     <- Mat[ c( c(5:44) , c((5+50):(44+50)) , c((5+100):(44+100)) ),];
ListMatrixLearnTest$ListMatTest[[1]]      <- Mat[ c( c(1:4,(1+50):(4+50),(1+100):(4+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[1]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[1]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[2]]     <- Mat[ c( c(1:4,9:44) , c((1+50):(4+50),(9+50):(44+50)), c((1+100):(4+100),(9+100):(44+100))) , ];
ListMatrixLearnTest$ListMatTest[[2]]      <- Mat[ c( c(5:8,(5+50):(8+50),(5+100):(8+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[2]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[2]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[3]]     <- Mat[ c( c(1:8,13:44) , c((1+50):(8+50),(13+50):(44+50)) , c((1+100):(8+100),(13+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[3]]      <- Mat[ c( c(9:12,(9+50):(12+50),(9+100):(12+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[3]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[3]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[4]]     <- Mat[ c( c(1:12,17:44) , c((1+50):(12+50),(17+50):(44+50)) , c((1+100):(12+100),(17+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[4]]      <- Mat[ c( c(13:16,(13+50):(16+50),(13+100):(16+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[4]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[4]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[5]]     <- Mat[ c( c(1:17,21:44) , c((1+50):(17+50),(21+50):(44+50)) , c((1+100):(17+100),(21+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[5]]      <- Mat[ c( c(17:20,(17+50):(20+50),(17+100):(20+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[5]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[5]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[6]]     <- Mat[ c( c(1:20,25:44) , c((1+50):(20+50),(25+50):(44+50)) , c((1+100):(20+100),(25+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[6]]      <- Mat[ c( c(21:24,(21+50):(24+50),(21+100):(24+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[6]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[6]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[7]]     <- Mat[ c( c(1:24,30:44) , c((1+50):(24+50),(30+50):(44+50)) , c((1+100):(24+100),(30+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[7]]      <- Mat[ c( c(25:29,(25+50):(29+50),(25+100):(29+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[7]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[7]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[8]]     <- Mat[ c( c(1:29,34:44) , c((1+50):(29+50),(34+50):(44+50)) , c((1+100):(29+100),(34+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[8]]      <- Mat[ c( c(30:33,(30+50):(33+50),(30+100):(33+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[8]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[8]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[9]]     <- Mat[ c( c(1:33,38:44) , c((1+50):(33+50),(38+50):(44+50)) , c((1+100):(33+100),(38+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[9]]      <- Mat[ c( c(34:37,(34+50):(37+50),(34+100):(37+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[9]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[9]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[10]]     <- Mat[ c( c(1:37,42:44) , c((1+50):(37+50),(42+50):(44+50)) , c((1+100):(37+100),(42+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[10]]      <- Mat[ c( c(38:41,(38+50):(41+50),(38+100):(41+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[10]] <- Mat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[10]]      <- Mat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

return(ListMatrixLearnTest);
}

########################################################################
# Compute Test Evaluation 
#
########################################################################

# Usage:
#   ClusterEvalTest(MetOpt=1, MetLab=1, Nu=0.5, q=40, K=1, G=15, Cx=0, Cy=0, DName="iris", fileIn="D:\\R\\library\\svcR\\");
#

ClusterEvalTest<- function (
  MetOpt="",				# method stoch (1) or quadprog (2)
  MetLab="",				# method grid  (1) or mst      (2) or knn (3)
  Nu="",			     	# nu parameter
  q="",					# q parameter
  K="",					# k parameter, k nearest neigbours for grid
  G="",					# g parameter, grid size
  Cx="", 				# choice of x component to display
  Cy="", 				# choice of y component to display
  DName="",				# data name
  fileIn=""				# a file path of data
    ) {

#parameters init
assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
assign("q",     q,  env=SvcEnv, inherits = FALSE); 
assign("Ngrid", G,  env=SvcEnv, inherits = FALSE);
assign("Knn",   K,  env=SvcEnv, inherits = FALSE);
GlobalTime = 0;
fileName = paste(tempdir(),"\\","sortie.txt", sep='');
assign( "FileOutput", file(fileName, "w") , env=SvcEnv, inherits = FALSE);		# open an output file connection
AvPrecision = 0;
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), ncol = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), byrow = FALSE, dimnames = NULL), env=SvcEnv, inherits = FALSE);
assign("NumPoints", array( list(), (get("Ngrid", env=SvcEnv, inherits = FALSE)+1)*(get("Ngrid", env=SvcEnv, inherits = FALSE)+1)), env=SvcEnv, inherits = FALSE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=SvcEnv, inherits = FALSE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- chargeMatrix(DataName=DName, PathIn=fileIn);						# data matrix structure loading
Alert("", "ok", "\n");

TimeNow <- proc.time();					# catch time reference
#MemBeg  <- memory.size(max = FALSE);			# catch memory reference
	
Alert("", "two-feature selection...", "\t");
if( cx != 0 ) {
	Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
}
else { 
	if( sign(min(Matrice$Mat)) < 0 ) {
		MatAdjCoa   = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
	}
	else
		MatAdjCoa   = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
	Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
} #EndIf
Alert("", "ok", "\n");
	
Alert("", "sampling matrix...", "\t\t");
ListMatrixLearnTest = DataSampling(Mat=Matrice$Mat);
cx = 1; cy = 2;
MinMaxMat(Mat=Matrice$Mat, Cx=cx, Cy=cy);
Alert("", "ok", "\n");

for( IndMat in 1:10 ) {

	Alert("\n", paste("test number...",IndMat) , "\n\n");

	Matrice$Mat = as.data.frame( ListMatrixLearnTest$ListMatLearn[IndMat] );
	MatriceTest = as.data.frame( ListMatrixLearnTest$ListMatTest[IndMat]  );				

	Alert("", "kernel matrix...", "\t\t");
	assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ;  
	MatriceK   <- calcKernelMatrix(Matrice$Mat);						# kernel matrix computation
	Alert("", "ok", "\n");

	# langrange multiplier computation
	Alert("", "lagrange coefficients...", "\t");
	if( MetOpt == 1 )
			WVectorsYA <- CalcWcluster(MatriceKern=MatriceK)
	else
		if( MetOpt == 2)
			WVectorsYA <- OptimQuadProgWcluster(MatriceKern=MatriceK)
	Alert("", "ok", "\n");

	Alert("", "radius computation...", "\t\t"); 
	assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=SvcEnv, inherits = FALSE);  # computation a cluster radius
	assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=SvcEnv, inherits = FALSE);  	          # computation of delta R
	Alert("", "ok", "\n");

	if(MetLab == 1 ) {
		Alert("", "grid labeling...", "\n"); 
		Alert("\t\t", "grid clustering...", "");
		NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, cx, cy, WYA=WVectorsYA$A);	# clusters assignment
		Alert("\t", "ok", "\n");
		Alert("\t\t", "match grid...", ""); 
		ClassPoints   = MatchGridPoint(Mat=MatriceTest, NumCluster=NumberCluster, Cx=cx, Cy=cy, Knn=K);  # cluster assignment
		Alert("\t\t", "ok", "\n");
		Alert("\t\t", "evaluation...", "");
		MisClass      = Evaluation(Mat=MatriceTest, NBClass=NumberCluster, Cx=cx, Cy=cy, ClassPoints=ClassPoints);					# evaluation
		AvPrecision   = AvPrecision + Precision;
	}
	else if( MetLab == 2 ){
		Alert("", "mst labeling/eval...", "\n"); 
		ClassPoints = MST_labelling(Mat=Matrice$Mat, MatK= MatriceK, WYA=WVectorsYA, Cx=cx, Cy=cy);
		Alert("\t\t", "evaluation...", "");
		MisClass    = Evaluation(Mat=Matrice$Mat, NBClass=max(ClassPoints), Cx=cx, Cy=cy, ClassPoints=ClassPoints);
		Alert("\t\t", "ok", "\n");
		Alert("\t\t\t\t", "ok", "\n");
	}
	else if( MetLab == 3 ){
		Alert("", "knn labeling/eval...", "\n");
		ClassPoints    = KNN_labelling(Mat=Matrice$Mat, MatK= MatriceK, WYA=WVectorsYA, Cx=cx, Cy=cy);
		Alert("\t\t", "evaluation...", "");
		MisClass = Evaluation(Mat=Matrice$Mat, NBClass=max(ClassPoints), Cx=cx, Cy=cy, ClassPoints=ClassPoints);
		Alert("\t\t", "...ok", "\n");
		Alert("\t\t\t\t", "ok", "\n");
	}

	invisible(gc());					# freeing memory

	GlobalTime = GlobalTime + ( proc.time() - TimeNow ) ;

} #EndforIndMat

options(digits=3); cat ("PrecisionGlobale=", AvPrecision/10, "% \t", "\n", sep="", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
print("time consuming per process"); print(GlobalTime/10);	# output time consuming
close(get("FileOutput", env=SvcEnv, inherits = FALSE));

return("end-of-routine");

}

########################################################################
# Compute Final Evaluation 
#
########################################################################

# Usage:
#   ClusterEvalFinal(MetOpt=1, MetLab=1, Nu=0.5, q=40, K=1, G=15, Cx=0, Cy=0, DName="iris", fileIn="D:\\R\\library\\svcR\\");
#

ClusterEvalFinal<- function (
  MetOpt="",				# method stoch (1) or quadprog (2)
  MetLab="",				# method grid  (1) or mst      (2) or knn (3)
  Nu="",			     	# nu parameter
  q="",					# q parameter
  K="",					# k parameter, k nearest neigbours for grid
  G="",					# g parameter, grid size
  Cx="", 				# choice of x component to display
  Cy="", 				# choice of y component to display
  DName="",				# data name
  fileIn=""				# a file path of data
    ) {

#parameters init
assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
assign("q",     q,  env=SvcEnv, inherits = FALSE); 
assign("Ngrid", G,  env=SvcEnv, inherits = FALSE);
assign("Knn",   K,  env=SvcEnv, inherits = FALSE);
GlobalTime = 0;
fileName = paste(tempdir(),"\\","sortie.txt", sep='');
assign( "FileOutput", file(fileName, "w") , env=SvcEnv, inherits = FALSE);		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), ncol = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), byrow = FALSE, dimnames = NULL), env=SvcEnv, inherits = FALSE);
assign("NumPoints", array( list(), (get("Ngrid", env=SvcEnv, inherits = FALSE)+1)*(get("Ngrid", env=SvcEnv, inherits = FALSE)+1)), env=SvcEnv, inherits = FALSE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=SvcEnv, inherits = FALSE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- chargeMatrix(DataName=DName, PathIn=fileIn);		# data matrix structure loading
Alert("", "ok", "\n");

TimeNow <- proc.time();					# catch time reference
#MemBeg  <- memory.size(max = FALSE);			# catch memory reference
	
Alert("", "two-feature selection...", "\t");
if( cx != 0 ) { 
		Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
}
else { 
		if( sign(min(Matrice$Mat)) < 0 ) {
			MatAdjCoa    = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
		}
		else
			MatAdjCoa    = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
		Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
} #EndIf
Alert("", "ok", "\n");

Alert("", "sampling matrix...", "\t\t");
ListMatrixLearnTest = DataSampling(Mat=Matrice$Mat);
Matrice$Mat = as.data.frame( ListMatrixLearnTest$ListMatLearnTest[1] );
MatriceEval = as.data.frame( ListMatrixLearnTest$ListMatEval[1] );
cx = 1; cy = 2;
MinMaxMat(Mat=Matrice$Mat, Cx=cx, Cy=cy);				
Alert("", "ok", "\n");

Alert("", "kernel matrix...", "\t\t");
assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ; 
MatriceK   <- calcKernelMatrix(Matrice$Mat);						# kernel matrix computation
Alert("", "ok", "\n");

# lagrange multiplier computation
Alert("", "lagrange coefficients...", "\t");
if( MetOpt == 1 )
		WVectorsYA <- CalcWcluster(MatriceKern=MatriceK)
else
	if( MetOpt == 2)
		WVectorsYA <- OptimQuadProgWcluster(MatriceKern=MatriceK)
Alert("", "ok", "\n");

Alert("", "radius computation...", "\t\t"); 
assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=SvcEnv, inherits = FALSE);  # computation a cluster radius
assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=SvcEnv, inherits = FALSE);  	          # computation of delta R
Alert("", "ok", "\n");

if(MetLab == 1 ) {
	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, cx, cy, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=MatriceEval, NumCluster=NumberCluster, Cx=cx, Cy=cy, Knn=K);  # cluster assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "evaluation...", "");
	MisClass      = Evaluation(Mat=MatriceEval, NBClass=NumberCluster, Cx=cx, Cy=cy, ClassPoints=ClassPoints);					# evaluation
	Alert("\t\t", "ok", "\n");
	Alert("\t\t\t\t", "ok", "\n");
}
else if( MetLab == 2 ){
	Alert("", "mst labeling/eval...", "\n"); 
	ClassPoints = MST_labelling(Mat=Matrice$Mat, MatK= MatriceK, WYA=WVectorsYA, Cx=cx, Cy=cy);
	Alert("\t\t", "evaluation...", "");
	MisClass    = Evaluation(Mat=Matrice$Mat, NBClass=max(ClassPoints), Cx=cx, Cy=cy, ClassPoints=ClassPoints);
	Alert("\t\t", "ok", "\n");
	Alert("\t\t\t\t", "ok", "\n");
}
else if( MetLab == 3 ){
	Alert("", "knn labeling/eval...", "\n");
	ClassPoints    = KNN_labelling(Mat=Matrice$Mat, MatK= MatriceK, WYA=WVectorsYA, Cx=cx, Cy=cy);
	Alert("\t\t", "evaluation...", "");
	MisClass = Evaluation(Mat=Matrice$Mat, NBClass=max(ClassPoints), Cx=cx, Cy=cy, ClassPoints=ClassPoints);
	Alert("\t\t", "...ok", "\n");
	Alert("\t\t\t\t", "ok", "\n");
}

invisible(gc());					# freeing memory
GlobalTime = ( proc.time() - TimeNow ) ;
options(digits=3); cat ("PrecisionGlobale=", Precision, "% \t", "\n", sep="", file=get("FileOutput", env=SvcEnv, inherits = FALSE));
print("time consuming"); print(GlobalTime);	# output time consuming
close(get("FileOutput", env=SvcEnv, inherits = FALSE));

return("end-of-routine");

}

########################################################################
# Test scalability 
#
########################################################################

# Usage:
#   ModelScalable(MetOpt=1, MetLab=1, Nu=0.5, q=40, K=1, G=15, Cx=0, Cy=0, DName="iris", fileIn="D:\\R\\library\\svcR\\");
#

ModelScalable<- function (
  MetOpt="",				# method stoch (1) or quadprog (2)
  MetLab="",				# method grid  (1) or mst      (2) or knn (3)
  Nu="",			     	# nu parameter
  q="",					# q parameter
  K="",					# k parameter, k nearest neigbours for grid
  G="",					# g parameter, grid size
  Cx="", 				# choice of x component to display
  Cy="", 				# choice of y component to display
  DName="",				# data name
  fileIn=""				# a file path of data
    ) {

#parameters init
assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
assign("q",     q,  env=SvcEnv, inherits = FALSE); 
assign("Ngrid", G,  env=SvcEnv, inherits = FALSE);
assign("Knn",   K,  env=SvcEnv, inherits = FALSE);
GlobalTime = 0;
fileName = paste(tempdir(),"\\","sortie.txt", sep='');
assign( "FileOutput", file(fileName, "w"), env=SvcEnv, inherits = FALSE );		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), ncol = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), byrow = FALSE, dimnames = NULL), env=SvcEnv, inherits = FALSE);
assign("NumPoints", array( list(), (get("Ngrid", env=SvcEnv, inherits = FALSE)+1)*(get("Ngrid", env=SvcEnv, inherits = FALSE)+1)), env=SvcEnv, inherits = FALSE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=SvcEnv, inherits = FALSE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- chargeMatrix(DataName=DName, PathIn=fileIn);		# data matrix structure loading
MM = Matrice$Mat
Alert("", "ok", "\n");

TimeNow        <- proc.time();			# catch time reference
FileOutputTime <- file(paste(tempdir(),"\\","timesvc.txt", sep='') , "w");

for(NLine in 3:300 ) {
	
	TimeNow <- proc.time();

	Alert("", "sampling matrix...", "\t\t");
	Matrice$Mat = MM[1:NLine,];
	Alert("\t\t", "ok", "\n");

	Alert("", "two-feature selection...", "\t");
	if( cx != 0 ) { 
		Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
	}
	else { 
		if( sign(min(Matrice$Mat)) < 0 ) {
			MatAdjCoa    = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
		}
		else
			MatAdjCoa    = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
		Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
	} #EndIf
	Alert("", "ok", "\n");
	MinMaxMat(Mat=Matrice$Mat, Cx=1, Cy=2);				
	Alert("", "ok", "\n");

	Alert("", "kernel matrix...", "\t\t");
	assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ; 
	MatriceK   <- calcKernelMatrix(Matrice$Mat);						# kernel matrix computation
	Alert("", "ok", "\n");

	# lagrange multiplier computation
	Alert("", "lagrange coefficients...", "\t");
	WVectorsYA <- CalcWcluster(MatriceKern=MatriceK)
	Alert("", "ok", "\n");

	Alert("", "radius computation...", "\t\t"); 
	assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=SvcEnv, inherits = FALSE);  # computation a cluster radius
	assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=SvcEnv, inherits = FALSE);  	          # computation of delta R
	Alert("", "ok", "\n");

	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, 1, 2, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=Matrice$Mat, NumCluster=NumberCluster, Cx=1, Cy=2, Knn=K);  # cluster assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "evaluation...", "");
	MisClass      = Evaluation(Mat=Matrice$Mat, NBClass=NumberCluster, Cx=1, Cy=2, ClassPoints=ClassPoints);					# evaluation	
	Alert("\t\t\t\t", "ok", "\n");

	GlobalTime = ( proc.time() - TimeNow ) ;
	options(digits=3); 
	cat("NLine=", NLine, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=SvcEnv, inherits = FALSE), '\n', file=FileOutputTime);	# output time consuming
	cat("NLine=", NLine, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=SvcEnv, inherits = FALSE), '\n');	# output time consuming

} #Endfori

invisible(gc());					# freeing memory
close(get("FileOutput", env=SvcEnv, inherits = FALSE));
close(FileOutputTime);

return("end-of-routine");

}

########################################################################
# Test clusterability 
#
########################################################################

# Usage:
#   ModelClusterable(MetOpt=1, MetLab=1, Nu=0.5, q=40, K=1, G=15, Cx=0, Cy=0, DName="iris", fileIn="D:\\Rbuild\\test\\");
#

ModelClusterable<- function (
  MetOpt="",				# method stoch (1) or quadprog (2)
  MetLab="",				# method grid  (1) or mst      (2) or knn (3)
  Nu="",			     	# nu parameter
  q="",					# q parameter
  K="",					# k parameter, k nearest neigbours for grid
  G="",					# g parameter, grid size
  Cx="", 				# choice of x component to display
  Cy="", 				# choice of y component to display
  DName="",				# data name
  fileIn=""				# a file path of data
    ) {

#parameters init
assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
assign("q",     q,  env=SvcEnv, inherits = FALSE); 
assign("Ngrid", G,  env=SvcEnv, inherits = FALSE);
assign("Knn",   K,  env=SvcEnv, inherits = FALSE);
GlobalTime = 0;
fileName = file.path(tempdir(), "sortie.txt");
assign( "FileOutput", file(fileName, "w"), env=SvcEnv, inherits = FALSE );		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), ncol = (get("Ngrid", env=SvcEnv, inherits = FALSE)+1), byrow = FALSE, dimnames = NULL), env=SvcEnv, inherits = FALSE);
assign("NumPoints", array( list(), (get("Ngrid", env=SvcEnv, inherits = FALSE)+1)*(get("Ngrid", env=SvcEnv, inherits = FALSE)+1)), env=SvcEnv, inherits = FALSE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=SvcEnv, inherits = FALSE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

FileOutputTime <- file( file.path(fileIn, "timesvc.txt"), "w"); #file( paste(fileIn, "timesvc.txt", sep="") , "w");

ListNu = c(0.1, 0.5, 1);
ListQ  = c(0.5, 1, 10, 100);

for(IndNu in 1:length(ListNu) ) {
	
   Nu = ListNu[IndNu];
   assign("nu",    Nu, env=SvcEnv, inherits = FALSE);
	
   for(IndQ in 1:length(ListQ)) {
	
	q = ListQ[IndQ];
	assign("q",     q,  env=SvcEnv, inherits = FALSE); 
	
	TimeNow <- proc.time();

	Alert("", "loading matrix...", "\t\t");
	Matrice    <- chargeMatrix(DataName=DName, PathIn=fileIn);		# data matrix structure loading
	Alert("", "ok", "\n");
 
	Alert("", "two-feature selection...", "\t");
	if( cx != 0 ) { 
		Matrice$Mat = Matrice$Mat[,c(cx,cy,ncol(Matrice$Mat))];
	}
	else { 
		if( sign(min(Matrice$Mat)) < 0 ) {
			MatAdjCoa    = dudi.pca(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
		}
		else
			MatAdjCoa    = dudi.coa(as.data.frame(Matrice$Mat[,1:(ncol(Matrice$Mat)-1)]), scan = FALSE);
		Matrice$Mat = as.data.frame( c( MatAdjCoa$li , as.data.frame(Matrice$Mat[,ncol(Matrice$Mat)]) ) );
	} #EndIf
	MinMaxMat(Mat=Matrice$Mat, Cx=1, Cy=2);				
	Alert("", "ok", "\n");

	Alert("", "kernel matrix...", "\t\t");
	assign("MaxValA", 1/(get("nu", env=SvcEnv, inherits = FALSE)*nrow(Matrice$Mat)), env=SvcEnv, inherits = FALSE ) ; 
	MatriceK   <- calcKernelMatrix(Matrice$Mat);						# kernel matrix computation
	Alert("", "ok", "\n");

	# lagrange multiplier computation
	Alert("", "lagrange coefficients...", "\t");
	WVectorsYA <- CalcWcluster(MatriceKern=MatriceK)
	Alert("", "ok", "\n");

	Alert("", "radius computation...", "\t\t"); 
	assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=SvcEnv, inherits = FALSE);  # computation a cluster radius
	assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=SvcEnv, inherits = FALSE);  	          # computation of delta R
	Alert("", "ok", "\n");

	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, 1, 2, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=Matrice$Mat, NumCluster=NumberCluster, Cx=1, Cy=2, Knn=K);  # cluster assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "evaluation...", "");
	MisClass      = Evaluation(Mat=Matrice$Mat, NBClass=NumberCluster, Cx=1, Cy=2, ClassPoints=ClassPoints);					# evaluation	
	Alert("\t\t\t\t", "ok", "\n");

	GlobalTime = ( proc.time() - TimeNow ) ;
	options(digits=3); 
	cat("q=",  get("q", env=SvcEnv, inherits = FALSE), '\t', "nu=", get("nu", env=SvcEnv, inherits = FALSE), '\t', "#cluster", NumberCluster, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=SvcEnv, inherits = FALSE), '\n', file=FileOutputTime);	# output time consuming
	cat("q=",  get("q", env=SvcEnv, inherits = FALSE), '\t' , "nu=", get("nu", env=SvcEnv, inherits = FALSE), '\t', "#cluster", NumberCluster, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=SvcEnv, inherits = FALSE), '\n');	# output time consuming
	
	invisible(gc());					# freeing memory

   } #EndforIndQ

} #EndforIndNu

close(get("FileOutput", env=SvcEnv, inherits = FALSE));
close(FileOutputTime);

return("end-of-routine");

}

########################################################################
# Exporting clusters in a text file 
#
########################################################################

# Usage:
#   ExportClusters(MatriceVar=Matrice$Var, CPoints=ClassPoints, DName="iris", pathOut=fileIn);
#

ExportClusters<- function (
  MatriceVar="",			# var name
  CPoints="",				# vector points / cluster
  DName="",				# name of output
  pathOut=""				# name of output path
    ) {

#opening output stream
CluOutput <- file( file.path(pathOut, paste(DName, "_clu.txt", sep="")) , "w");

#sorting by cluster index
SortedClassPoints = sort(CPoints, method = "sh", index.return = TRUE);

#output points class
for(i in 1:nrow(MatriceVar) ){
 if( i == 1 || SortedClassPoints$x[i] != SortedClassPoints$x[i-1] )
	cat('\n', "cluster", '\t', SortedClassPoints$x[i], '\n', file=CluOutput);
 cat("item", '\t', as.character(MatriceVar[ SortedClassPoints$ix[i] , 1 ]), '\n', file=CluOutput);

} #finfori

close(CluOutput)
}

########################################################################
# message tracing 
#
########################################################################

# Usage:
#   Alert("here", "\n");
#

Alert<- function (
  bTab="",				# pre tabulations
  Mes="",				# message
  aTab=""				# post tabulations
    ) {
cat(bTab, Mes, aTab);flush.console();
}


########################################################################
## test function
########################################################################

#
TestWrapper <- function()
{

print("en route...");

#findModelCluster(MetOpt=1, MetLab=1, KernChoice=1, Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2, DName="iris", fileIn="D:\\R\\library\\svcR\\data\\");
findModelCluster(MetOpt=1, MetLab=1, KernChoice=1, Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2, DName="iris2", fileIn="D:\\R\\library\\svcR\\data\\");

print("fin de svcR");

# done. 
#
}

########################################################################
##End of program##
########################################################################

########################################################################
## call test function
########################################################################

#TestWrapper();
