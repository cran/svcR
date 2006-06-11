########################################################################
# svcR: an open source library in R 
#   for SVC computation
#
########################################################################

## History of this library 
##  svcR # Apr-27-05
##   written  by Nicolas Turenne  
##                 # v0.0     alpha release      # Apr-27-05   
##                 # v0.9     alpha release      # Sep-26-05   
##                 # v1.0     beta  release      # Jun-26-06   
## source("D:\\R\\library\\svc\\svcR.txt")
## load("d:\\r\\library\\svc\\svc.RData")
## save.image("d:\\r\\library\\svc\\svc.RData")

########################################################################
# Load library
########################################################################
library(quadprog)
library(ade4)
library(spdep)

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
AroundNull   = 0.1;
AroundNullVA = 0.00005;

# Data Structure of Criterium
WVectorsYA <- list(W="",Y="",A="")

# framing display
par(mfrow = c(2, 2));					

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
fileName   = "sortie.txt";
FileOutput = 0;
FileOutput <- file( paste(tempdir(),"\\",fileName, sep='') , "w");		# open an output file connection

Precision  = 0;

########################################################################
# Main function (SVC)
########################################################################

# Usage:
#   findModelCluster(MetOpt=1, MetLab=1, Nu=0.5, q=40, K=1, G=15, Cx=1, Cy=2, DName="iris", fileIn="D:\\R\\library\\svcR\\")
#
envir=globalenv();

findModelCluster<- function (
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
assign("nu",    Nu, env=envir, inherits = TRUE);
assign("q",     q,  env=envir, inherits = TRUE); 
assign("Ngrid", G,  env=envir, inherits = TRUE);
assign("Knn",   K,  env=envir, inherits = TRUE);
TimeNow    <- proc.time();			# catch time reference
#MemBeg     <- memory.size(max = FALSE);		# catch memory reference
FileOutput <- file(fileName, "w");		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=envir, inherits=TRUE)+1), ncol = (get("Ngrid", env=envir, inherits=TRUE)+1), byrow = FALSE, dimnames = NULL), env=envir, inherits= TRUE);
assign("NumPoints", array( list(), (get("Ngrid", env=envir, inherits=TRUE)+1)*(get("Ngrid", env=envir, inherits=TRUE)+1)), env=envir, inherits= TRUE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=envir, inherits= TRUE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

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
Alert("", "ok", "\n");

Alert("", "min max calculation...", "\t");
cx = 1; cy = 2;
MinMaxMat(Mat=Matrice$Mat, Cx=cx, Cy=cy);				
Alert("", "ok", "\n");

Alert("", "kernel matrix...", "\t\t");
assign("MaxValA", 1/(get("nu", env=envir, inherits=TRUE)*nrow(Matrice$Mat)), env=envir, inherits= TRUE ) ; 
#cat("MaxValA", '\n', get("MaxValA", env=envir, inherits=TRUE), file=FileOutput);
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
assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=envir, inherits= TRUE);  # computation a cluster radius
assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=envir, inherits= TRUE);  	          # computation of delta R
Alert("", "ok", "\n");

if(MetLab == 1 ) {
	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, cx, cy, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=Matrice$Mat, NumCluster=NumberCluster, Cx=cx, Cy=cy, Knn=K);  # cluster assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "evaluation...", "");
	MisClass      = Evaluation(Mat=Matrice$Mat, NBClass=NumberCluster, Cx=cx, Cy=cy, ClassPoints=ClassPoints);					# evaluation
	Alert("\t\t", "ok", "\n");
	Alert("\t\t\t\t" ,"ok", "\n");
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
else {
	#MatAdj = Adjacency(Mat=Matrice$Mat, MatK=MatriceK, WYA=WVectorsYA);
}

Alert("", "export...", "\t\t\t");
ExportClusters(MatriceVar=Matrice$Var, CPoints=ClassPoints, DName=DName, pathOut=fileIn);
Alert("", "ok", "\n");

Alert("", "display...", "\t\t\t");
DisplayData(Matrice$Mat, MatriceK, WVectorsYA, cx, cy, MisClass);			# output results
Alert("", "ok", "\n");

#rm( MatriceK, VW, RO ); 
invisible(gc());					# freeing memory
#rm(list=ls(all=TRUE))					# delete all in workspace

print("time consuming");print(proc.time() - TimeNow);	# output time consuming
#print("Max Memory");print(GMHmax);
#print("Memory At Beginning");print(MemBeg);		# output memory consuming
#print("Memory Consuming");print(GMHmax-MemBeg);		# output memory consuming
close(FileOutput);

return("end-of-routine");
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
    ListMis=""			# list of misclassified data
    ) {

# plot the data matrix
#cat("matrix", "\n", file=FileOutput); write(as.matrix(Mat), sep="\t", file=FileOutput);					

# plot Y class in the case of an svm
if( WYA$Y[1] != "" )
	plot( 1:length(WYA$Y), WYA$Y , xlab="point", ylab="Classe");

# plot W history
TAB = get("TabW", env=envir, inherits=FALSE);
for(Ind in 1:length(TAB)){
	Val = get( paste("TabW[",Ind,"]", sep="") , env=envir, inherits=FALSE);
	TAB[Ind] = Val;
}
plot(TAB, xlab="#iterations", ylab="W");

# plot Lagrange parameters
if( WYA$A[1] != "" )
	plot( 1:length(WYA$A), WYA$A , xlab="point", ylab="Coefficient de Lagrange");

# plots data numbers being SV, make SV in red
plot(Mat[,Cx],Mat[,ncol(Mat)], xlab="data points", ylab="classes in data matrix")

SV = ListSVPoints(VA=WYA$A);			# list of support vectors
#cat("list des SV", '\n', file=FileOutput); write(t(SV), file=FileOutput); 
for(i in 1:length(SV) ) {				
		ISV = SV[i];
		if( !is.na(ISV) )
			points(Mat[ISV,Cx], Mat[ISV,ncol(Mat)], pch = 24, col = "red", bg = "yellow", cex = 1)
} #fin for i

Grid(Mat=Mat, MatK= MatK, WYA=WYA$A, ListSV=SV , Cx, Cy, ListMis);

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

Mat              <- list(Sparse="",Mat="",Var="",Att="");

if( nchar(PathIn) < 2 ){
	data(iris_mat); data(iris_att); data(iris_var);
	Mat$Sparse     = iris_mat; 
	Mat$Att        = iris_att; 
	Mat$Var        = iris_var;
}
else {
	Mat$Sparse     = read.table( paste(PathIn, DataName, "_mat", ".txt", sep="") , sep=" "); 
	Mat$Att        = read.table( paste(PathIn, DataName, "_att", ".txt", sep="") ); 
	Mat$Var        = read.table( paste(PathIn, DataName, "_var", ".txt", sep="") );
}
IndiceLineMax = IndiceColMax = 0;	# we look for max line and column indices
for(i in 1:length(Mat$Sparse[,1]) ){
		IndiceLine = Mat$Sparse[[1]][i];
		IndiceCol  = Mat$Sparse[[2]][i];
		if( IndiceLineMax < IndiceLine ) IndiceLineMax =  IndiceLine;
		if( IndiceColMax  < IndiceCol )  IndiceColMax  =  IndiceCol;			
} #finfori

NCols   = IndiceColMax; # we initialize the full matrix
NRows   = IndiceLineMax;
Mat$Mat = matrix(data = 0, nrow = NRows, ncol = NCols, byrow = FALSE, dimnames = NULL)

for(i in 1:length(Mat$Sparse[,1]) ){	# we fill the full matrix
		IndiceLine = Mat$Sparse[[1]][i];
		IndiceCol  = Mat$Sparse[[2]][i];
		Val        = as.numeric( Mat$Sparse[[3]][i] );
		Mat$Mat[IndiceLine,IndiceCol] = Val;
}#finfori

#write( "Matrice$Mat \n", file=FileOutput); write( t(Mat$Mat), sep="\t", file=FileOutput);
		
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
    matrix=""                           # matrix
    ) {

NCols = ncol(matrix)-1;
NRows = nrow(matrix);

MatriceKernel = matrix(data = 0, nrow = NRows, ncol = NRows, byrow = FALSE, dimnames = NULL)

#cat( "NCols", '\t', NCols, '\t', "NRows", NRows, file=FileOutput); 
i <- 1; j <- 1;
for(i in 1:NRows) {
  for( j in i:NRows) {
	MatriceKernel[i,j] <- 0;
	MatriceKernel[i,j] <- as.numeric( Kernel(Vec1=matrix[i,1:NCols], Vec2=matrix[j,1:NCols], Choice=1) );
	MatriceKernel[j,i] <- MatriceKernel[i,j];
  }

}

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

Res <- as.numeric( V1%*%V2 );

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

Res	<- 0;
q	<- get("q", env=envir, inherits=TRUE);

Res <- as.vector(V1, mode="numeric") - as.vector(V2, mode="numeric") ;
Res = exp( - q * as.numeric( Res %*% Res ) );

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
		#write("beta_i \n", file=FileOutput); write(VecteurA[i], file=FileOutput);  
		#cat("\n BoundS \t", BoundS, '\n', file=FileOutput);
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

Sum   <- 0;
Sum = sum( VecteurA );
#cat("Sum", '\n', Sum, file=FileOutput);

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

#WVectorsYA$A <- runif(Dim);
#WVectorsYA$A = MaxValA * WVectorsYA$A ;
#WVectorsYA$A = 2*WVectorsYA$A / Dim ;

nu           = get("nu", env=envir, inherits=TRUE);
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
	COEFFS[i]	= get("MaxValA", env=envir, inherits=TRUE);
	i		= i+1;
} #finwhile

cat( "MaxValA", get("MaxValA", env=envir, inherits=TRUE), '\t', "Sum MaxValA", length(COEFFS)*get("MaxValA", env=envir, inherits=TRUE), '\n');
cat( "sum(COEFFS)", sum(COEFFS), '\t', "1 + AroundNull", (1 + AroundNull), '\t', "1 - AroundNull", (1 - AroundNull), '\n' );

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

WYA	<- list(W="",Y="",A="")
W	<- MaxW; PreviousW <- MaxW ; Iter <- 1;
WYA$W	<- W;
# creation of a random vector (randomization uniforme runif) / (randomization gaussienne rnorm)
N	<- ncol(MatriceKern)
WYA$A	<- MakeA(Dim=N);

# computation of WYA$A bound sup
BoundSup	<- 1E+10;
nu		<- get("nu", env=envir, inherits=TRUE);
MaxValA		<- get("MaxValA", env=envir, inherits=TRUE);

if( nu && N )
	BoundSup = MaxValA ; # 1 / ( nu * N );
#cat("BoundSup ", BoundSup, "\n", file=FileOutput);

TabWPrec = MinW;

while( (W > TabWPrec &&  Iter <= MaxIter) ||  (Iter <= MaxIter) ) {
	
	PreviousW = W;
	W         = CritereWcluster(VecteurA=WYA$A, MatrixK=MatriceKern);
	#cat("Iter=", Iter, '\n', file=FileOutput);

	if( Iter > 1 ){
		ValW   = paste("TabW[",Iter-1, "]", sep=""); print("just");
		TabWPrec <- get(ValW, env=envir, inherits=TRUE);
		}
	
	if( W > TabWPrec ){
		ValW   = paste("TabW[",Iter,"]", sep="");
		assign(ValW, W, env=envir, inherits=TRUE);
	} 
	else {  
		ValW   = paste("TabW[",Iter,"]", sep="");
		assign(ValW, TabWPrec, env=envir, inherits=TRUE);
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

#cat("A", '\n', file=FileOutput); write(t(WYA$A), file=FileOutput);
#cat("W", '\n', t(TabW[Iter-1]) , file=FileOutput); write(t(TabW[Iter-1]) , file=FileOutput);

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

WYA	<- list(W="",Y="",A="");
WYA$W	<- MinW;
N	<- ncol(MatriceKern)
nu	<- get("nu", env=envir, inherits=TRUE);
MaxValA	<- get("MaxValA", env=envir, inherits=TRUE);

# computation of WYA$A bound sup
BoundSup <- 1E+10;
if( nu && N )
	BoundSup = MaxValA ; # 1 / ( nu * N );
#cat("BoundSup", BoundSup, "\n", file=FileOutput);

Dmat <- 2*MatriceKern;
dvec <- diag(MatriceKern);
Amat <- matrix(0, 2*N+1, N);
for( i in 1:N){ 
		Amat[1,i]      = 1;
		Amat[i+1,i]       = 1;
		Amat[i+N+1,i]  = -1;
}
bvec = as.vector( c(1, matrix(0,1,N),matrix(-BoundSup,1,N) ) );
S <- solve.QP(Dmat,-dvec,t(Amat),t(bvec), meq=1, factorized=TRUE);
WYA$A  <- S$solution;
WYA$W <- S$value;

#cat( "WYA$A", '\n', file=FileOutput); write( t(WYA$A) , file=FileOutput);
#cat( "WYA$W", '\n', file=FileOutput); write( t(WYA$W) , file=FileOutput);

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

N		<- length(VA);
ListP		<- vector(); length(ListP) = N;
MaxValA		<- get("MaxValA", env=envir, inherits=TRUE);

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

VecW <- 0;
for(i in 1:nrow(Mat) ) {
	VecW = VecW + WYA$A[i]*Mat[i,1:(ncol(Mat)-1)];
} #fin for i
#cat("VecW", '\t', file=FileOutput); write(t(VecW), file=FileOutput);

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

ro <- 0;
for(i in 1:nrow(Mat) ) {
	ro = ro + WYA$A[i]*Kernel(Vec1=VecW, Vec2=Mat[i,1:(ncol(Mat)-1)], Choice=1);
} #fin for i
#cat("ro", '\n', ro, file=FileOutput);

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

R        <- 0 ;
IndiceSV <- 0;
N        <- length(VA);
i        <- 1;
MaxValA  <- get("MaxValA", env=envir, inherits=TRUE);
 
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

cat(IndiceSV, '\n');
if( IndiceSV != 0 ) {
	R = RadiusData(IndicePoint=IndiceSV, VA=VA, MatK=MatrixK);
	#cat("indiceSV=", IndiceSV, '\t', "R=", R, '\n', file=FileOutput);
}
else {	
	R = 0;
	#cat("IndiceSV = 0", "\t", "AroundNullVA", AroundNullVA, "\t", "N=", N, '\n', file=FileOutput);
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

Dim = ncol(Mat) - 1;

Rad <- 0 ;
Rad = Kernel( Vec1=Vec, Vec2=Vec , Choice=1) ;

for( i in 1:nrow(Mat) ) {

        Rad = Rad - 2*VA[i]*Kernel( Vec1=Mat[i,1:Dim], Vec2=Vec[1:Dim] , Choice=1 );

	for( j in 1:nrow(Mat) ) {
		 Rad = Rad + VA[i]*VA[j]*MatK[i,j];
	} #fin for j

} #fin for i

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

Rad = MatK[IndicePoint,IndicePoint] ;

for( i in 1:N ) {
        Rad = Rad - 2*VA[i]*MatK[i,IndicePoint];
	for( j in 1:N ) {
		Rad = Rad + VA[i]*VA[j]*MatK[i,j];
	} #fin for j
} #fin for i

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

N  <- length(VA);
r  <- vector(); length(r) = N;
nu <- get("nu", env=envir, inherits=TRUE)
C  <- 1 / ( nu * N );

for( i in 1:N ) {
        
	R    = 0; RC = 0;
	r[i] = 0;
	if( VA[i] > C * 0.98 ) {
		R  = RadiusData(IndicePoint=i,VA=VA,MatK=MatK);
		RC = get("RadiusC", env=envir, inherits=TRUE);
		if( R > RC )
			r[i] =  R*R - RC*RC;
	} #finif
	#cat("RadiusC", RC, "\t", "R", R, "\t", "ri", r[i], '\n', file=FileOutput);
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

adj_flag = 1; # unless a point on the path exits the sphere - pair is adjacent
		
interval = 0.0;
while( interval < 1 && adj_flag ){

	z = Vec1 + interval * (Vec2 - Vec1);	    
	interval = interval + 0.3;
	R  = RadiusPoint(Vec=z, VA=WYA, Mat=Mat, MatK=MatK);	
	RC = get("RadiusC", env=envir, inherits=TRUE);
	r  = get("r", env=envir, inherits=TRUE);
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

#cat( "AdjacencyM", '\n', AdjacencyM, file=FileOutput);

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

#cat("ListItemCluster", '\n', file=FileOutput); write( t(ListItemCluster), file=FileOutput);

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

#cat("ListItemCluster", '\n', file=FileOutput); write( t(ListItemCluster), file=FileOutput);

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

#output points class
NumColClass = ncol(Mat);
#for(i in 1:nrow(Mat) ){
# cat("N=", '\t', i, "\t", "i=", Mat[i, Cx], "\t", "j=", Mat[i, Cy], "\t", "C=", ClassPoints[i], "\t", "Cdata=", Mat[i, NumColClass], '\n', file=FileOutput);
#} #finfori

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

#cat("labelclass", '\n', LabelClass, '\n', "Cardlabelclass", '\n' , CardLabelClass, file=FileOutput);

BestClass <- matrix(data = NA, nrow = NBClass, ncol = 2, byrow = FALSE, dimnames = NULL); 
ListVecClass = list();

for( IndClass in 1:nrow(Mat) )
	ListVecClass[[ IndClass ]]  = vector(); 
for(i in 1:nrow(Mat) ){
	L = length( ListVecClass[[ Mat[i, NumColClass] ]] );
	ListVecClass[[ Mat[i, NumColClass] ]][L+1] = ClassPoints[i];
} #finFor

for( IndClass in 1:NBClass ) {
	VecClassSorted = sort(  ListVecClass[[ as.numeric(LabelClass[IndClass]) ]]  );
	Counter = 1; BestClass[IndClass, 2] = 1; BestClass[IndClass, 1] = VecClassSorted[ 1 ];
	#cat( "VecClassSorted", '\n' , file=FileOutput); write( t(VecClassSorted), file=FileOutput);
	if( length(VecClassSorted) >= 2 )
	for( IndVecSorted in 2:length(VecClassSorted) ) {
		if( VecClassSorted[ IndVecSorted ] == VecClassSorted[ IndVecSorted - 1] ) {
			Counter = Counter + 1;
			if( Counter > BestClass[IndClass, 2] ){
				BestClass[IndClass, 2] = Counter;
				BestClass[IndClass, 1] = VecClassSorted[ IndVecSorted - 1];
			} #finif
			#cat("counter", Counter, "\t", "BestClass[IndClass, 2]", BestClass[IndClass, 2], "\t", "BestClass[IndClass, 1]", BestClass[IndClass, 1], '\n', file=FileOutput);
		} 
		else {
			Counter = 1;
		} #finif

	} #finIndVecSorted
}#finIndClass

#cat("BestClass", '\n', BestClass, file=FileOutput);

Precision  <- 0;
for(IndClass in 1:NBClass )  if( BestClass[IndClass, 1] != 0 )
	Precision <- Precision + BestClass[IndClass, 2];
Precision <- ( Precision/nrow(Mat) ) * 100;
assign("Precision", Precision, env=envir, inherits = TRUE);
#options(digits=3); cat("Precision=", Precision, "% \n", file=FileOutput);

#extraction of misclassified items
ListMis = vector(); IndListMis = 1;
for( i in 1:nrow(Mat) ) {
	for( IndClass in 1:NBClass ) {
		if( Mat[i, NumColClass] == as.numeric(LabelClass[IndClass]) ){
			#cat("best", '\t', BestClass[IndClass, 1], "\t class", ClassPoints[i], '\n', file=FileOutput);
			if( BestClass[IndClass, 1] != ClassPoints[i] || BestClass[IndClass, 1] == 0 ) {
				ListMis[IndListMis] = i ;
				IndListMis          = IndListMis + 1;
			} #endIf
		} #endIf
	}#endforIndVec
} #finFor

#write("ListMis \n", file=FileOutput); write( t(ListMis), file=FileOutput);

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
    ListMis=""		# list of misclassified data
    ) {

N       =  get("Ngrid", env=envir, inherits=TRUE);					
ListX   <- array(0, N); # we make the list of X values	
ListY   <- array(0, N); # we make the list of Y values
ListVec <- array(0, ncol(Mat)-1); # we make random vector
							 
MaxX = get("MaxX", env=envir, inherits=TRUE);
MinX = get("MinX", env=envir, inherits=TRUE);
MaxY = get("MaxY", env=envir, inherits=TRUE);
MinY = get("MinY", env=envir, inherits=TRUE);
for(i in 1:N )	{
	ListX[i] = ( MaxX - MinX )*( (i-1) / N ) + MinX;
	ListY[i] = ( MaxY - MinY )*( (i-1) / N ) + MinY;
}
							 
#zoom+
plot(NA,NA,xlim=c(MinX,MaxX), ylim=c(MinY,MaxY), xlab="Xgrid", ylab="Ygrid") ; # setting up co
							
for(i in 1:N ){
	for(j in 1:N ){
		x  = ListX[i]; ListVec[Cx] = x;
		y  = ListY[j]; ListVec[Cy] = y;
		R  = RadiusPoint(Vec=ListVec, VA=WYA, Mat=Mat, MatK=MatK);
		RC = get("RadiusC", env=envir, inherits=TRUE);
		r  = get("r", env=envir, inherits=TRUE);

		#cat("i", i,"\t", "j", j,"  ", "RadiusC", RC, "\t", file=FileOutput);
		#cat("r", r, "\t", "x", x, "\t", "y", y, "\t", "R", R, "\n", file=FileOutput);

		if(   (RC*RC + r) >=  R*R ){
			points(x, y, pch = 24, col = "yellow", bg = "yellow", cex = 1);
		}
		else
			points(x, y, pch = 21, col = "blue", bg = "black", cex = 1)
	}#finforj
		
}#finfori
							 
for(i in 1:length(WYA) ){
	R = RadiusPoint(Vec=Mat[i,1:(ncol(Mat)-1)], VA=WYA, Mat=Mat, MatK=MatK);
	if(   (RC*RC + r) >=  R*R ){
		points(Mat[i,Cx], Mat[i,Cy], pch = 21, col = "red", bg = "red", cex = 1);
	}
	else
		points(Mat[i,Cx], Mat[i,Cy], pch = 21, col = "blue", bg = "white", cex = 1)
}#finfori

# plots data numbers being SV, make SV in red
for(i in 1:length(ListSV) ) {				
		ISV = ListSV[i];
		if( !is.na(ISV) )
			points(Mat[ISV,Cx], Mat[ISV,Cy], pch = 24, col = "red", bg = "yellow", cex = 1)
} #fin for i

# plots misclassified data numbers being IMC, make in green
for(i in 1:length(ListMis)){
	IMC = ListMis[i];
	points(Mat[IMC,Cx], Mat[IMC,Cy], pch = 24, col = "green", bg = "yellow", cex = 1)
} #fin for i

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
				
assign("MaxX", max(Mat[, Cx]), env=envir, inherits= TRUE ) ; 
assign("MaxY", max(Mat[, Cy]), env=envir, inherits= TRUE ) ; 
assign("MinX", min(Mat[, Cx]), env=envir, inherits= TRUE ) ; 
assign("MinY", min(Mat[, Cy]), env=envir, inherits= TRUE ) ; 

MaxX = get("MaxX", env=envir, inherits=TRUE);
MinX = get("MinX", env=envir, inherits=TRUE);
MaxY = get("MaxY", env=envir, inherits=TRUE);
MinY = get("MinY", env=envir, inherits=TRUE);

#cat("MaxX", '\t', MaxX, "MinX", '\t', MinX, '\n', file=FileOutput); 

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

N	= get("Ngrid", env=envir, inherits=TRUE);					
ListVec <- array(0, ncol(Mat)-1); # we make random vector
MaxX	= get("MaxX", env=envir, inherits=TRUE);
MinX	= get("MinX", env=envir, inherits=TRUE);
MaxY	= get("MaxY", env=envir, inherits=TRUE);
MinY	= get("MinY", env=envir, inherits=TRUE);

InP = 1;
for(i in 1:N ){ 
	for(j in 1:N ){ 
		x  = ( MaxX - MinX )*( (i-1) / N ) + MinX ; ListVec[Cx] = x;
		y  = ( MaxY - MinY )*( (j-1) / N ) + MinY ; ListVec[Cy] = y;
		R  = RadiusPoint(Vec=ListVec, VA=WYA, Mat=Mat, MatK=MatK);
		RC = get("RadiusC", env=envir, inherits=TRUE);
		r  = get("r", env=envir, inherits=TRUE);
		NM_InP = paste("NumPoints[[",InP,"]]", sep="");
		assign(NM_InP, list(IndP="", IndX="", IndY="", InBall="", NumClus=""), env=envir, inherits= TRUE);
		PG  = paste("PointGrid[",i,",",j,"]", sep="");
		assign( PG , InP, env=envir, inherits= TRUE); 
		NM1 = paste("NumPoints[",InP,"][[1]][1]", sep="");
		NM2 = paste("NumPoints[",InP,"][[1]][2]", sep="");
		NM3 = paste("NumPoints[",InP,"][[1]][3]", sep="");
		NM4 = paste("NumPoints[",InP,"][[1]][4]", sep="");
		NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
		assign( NM1 , InP, env=envir, inherits= TRUE); 
		assign( NM2 , i, env=envir, inherits= TRUE); 
		assign( NM3 , j, env=envir, inherits= TRUE); 
		assign( NM5 , 0, env=envir, inherits= TRUE); 
		if(  (RC*RC + r) >=  R*R ) {
			assign(NM4, 1, env=envir, inherits= TRUE);
		}
		else {
			assign(NM4, 0, env=envir, inherits= TRUE);
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
		InP = get( PG , env=envir, inherits=TRUE);
		NM4 = paste("NumPoints[",InP,"][[1]][4]", sep="");
		NM2 = paste("NumPoints[",InP,"][[1]][2]", sep="");
		NM3 = paste("NumPoints[",InP,"][[1]][3]", sep="");
		NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
		if( get(NM4, env=envir, inherits=TRUE) == 1) {
			Signal	= 0;
			if( IndLC == 0 ) {	# create a new cluster
				ListClusters[[1]]        = list();
				ListClusters[[1]][1]     = InP; 
				assign( NM5 , 1, env=envir, inherits= TRUE);
				IndLC = IndLC + 1;
				Signal	= 1;
			}
			else {					# look for assign to an existing cluster
				for(ii in (i-1):(i+1)){
					for(jj in (j-1):(j+1)){

						if( ii == i && jj == j ){ }
						else {
							if( Signal == 0) for( IndListCluster in 1:length(ListClusters) ) {

								ActualCluster = ListClusters[[ IndListCluster ]];
								NbElemCluster = length(ActualCluster);

								if( Signal == 0) for( IndMembListClust in 1:NbElemCluster ) {

									ActualPoint = ActualCluster[[ IndMembListClust ]];
									NM2_cur = paste("NumPoints[",ActualPoint,"][[1]][2]", sep="");
									NM3_cur = paste("NumPoints[",ActualPoint,"][[1]][3]", sep="");
									if( get( NM2_cur , env=envir, inherits=TRUE) == ii && get( NM3_cur , env=envir, inherits=TRUE) == jj ){
										ListClusters[[ IndListCluster ]][ NbElemCluster+1 ] = InP;
										assign( NM5 , IndListCluster, env=envir, inherits= TRUE);
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
				assign( NM5 , num, env=envir, inherits= TRUE);
			} #endif

		} #endif
		
	}#finforj
} #finfori

#sortie test
for(i in 1:N ){
	for(j in 1:N ){
		PG  = paste("PointGrid[",i,",",j,"]", sep="");
		InP = get( PG , env=envir, inherits=TRUE);
		NM3 = paste("NumPoints[",InP,"][[1]][3]", sep="");
		NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
		#cat("N=", InP, "i=", i, "\t j=", j, "\t", "NumPoints=", get(NM3, env=envir, inherits=TRUE), "\t", file=FileOutput);
		#cat("NumCluster=", get(NM5, env=envir, inherits=TRUE), "\n", file=FileOutput);
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
		InP = get( PG , env=envir, inherits=TRUE);
		for(ii in (i-1):(i+1)){
			for(jj in (j-1):(j+1)) if( ii != 0 && jj != 0)if( ii <= N && jj <= N ) {
				PG  = paste("PointGrid[",ii,",",jj,"]", sep="");
				ActualPoint = get( PG , env=envir, inherits=TRUE);
				NM5_InP = paste("NumPoints[",InP,"][[1]][5]", sep="");
				NM5_Cur = paste("NumPoints[",ActualPoint,"][[1]][5]", sep="");
				if( get(NM5_InP, env=envir, inherits=TRUE) != 0 && get(NM5_Cur, env=envir, inherits=TRUE) != 0 ) { 
					#cat("ii", ii, "\t", "jj", jj, "\t", "C1=", get(NM5_InP, env=envir, inherits=TRUE), "\t", file=FileOutput);
					#cat("C2=", get(NM5_Cur, env=envir, inherits=TRUE), "\t", "i=", InP, "\t j=", ActualPoint, "\n", file=FileOutput);
					k = as.integer(get(NM5_InP, env=envir, inherits=TRUE));
					l = as.integer(get(NM5_Cur, env=envir, inherits=TRUE));
					#cat("k=", k, "\t", "l=", l, "\n", file=FileOutput);
					ClassConnex[ k, l ] = 1;
				}#endif
			} #endforjj
		} #endforii
	}#finforj
} #finfori

} #endifNmat

#cat( "ClassConnex", '\n', file=FileOutput); write( t(ClassConnex), sep='\t', file=FileOutput);

#deletion of bad classes
IndClassFusIn  = 0;
IndClassFusOut = NumberCluster;
while( IndClassFusOut > 1) {
	
	cat("IndClassFusOut", IndClassFusOut, "length(ListClusters)", length(ListClusters), '\n');

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
		assign( NM5 , indListClusters, env=envir, inherits= TRUE); 
	} #finforIndMembListClust
} #finforindListClusters

#output clusters
for(i in 1:(N*N) ){
	NM4 = paste("NumPoints[",i,"][[1]][4]", sep="");
	NM2 = paste("NumPoints[",i,"][[1]][2]", sep="");
	NM3 = paste("NumPoints[",i,"][[1]][3]", sep="");
	NM5 = paste("NumPoints[",i,"][[1]][5]", sep="");
	#cat("N=", i, "\t", "i=", get(NM2, env=envir, inherits=TRUE), "\t", "j=", get(NM3, env=envir, inherits=TRUE), "\t", file=FileOutput);
	#cat("InBal=", get(NM4, env=envir, inherits=TRUE), "\t", "C=", get(NM5, env=envir, inherits=TRUE), "\n", file=FileOutput);
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

NPoint = nrow(Mat);
NG     = get("Ngrid", env=envir, inherits=TRUE);
MaxX   = get("MaxX", env=envir, inherits=TRUE);
MinX   = get("MinX", env=envir, inherits=TRUE);
MaxY   = get("MaxY", env=envir, inherits=TRUE);
MinY   = get("MinY", env=envir, inherits=TRUE);

ClassPoints      <- array(0, NPoint ); 
#cat("N ", NPoint, '\t', "NumCluster ", NumCluster, '\t', "knn ", Knn, '\n', file=FileOutput);

if( NumCluster > 0)
for(i in 1:NPoint ){
	Xi  =  1 + round(  abs( Mat[i, Cx] - MinX )*NG / ( MaxX - MinX )  ) ;
	Yi  =  1 + round(  abs( Mat[i, Cy] - MinY )*NG / ( MaxY - MinY )  ) ; 
	if( Xi > NG) Xi = NG; if(Yi > NG) Yi = NG;

	PG  = paste("PointGrid[",Xi,",",Yi,"]", sep="");
	InP = get( PG , env=envir, inherits=TRUE);
	ScoreNeighbours  <- array(0, NumCluster );

        NM5 = paste("NumPoints[",InP,"][[1]][5]", sep="");
	#cat("N=", i, "\t", "Xi=", Xi, "\t", "Yi=", Yi, "\t", file=FileOutput);
	#cat("InP=", InP, "\t", "Class=", get(NM5, env=envir, inherits=TRUE), "\n", file=FileOutput);

	for(ii in (Xi-1):(Xi+1)){
			for(jj in (Yi-1):(Yi+1)) {
				if( jj == Yi && ii == Xi){
					#we leave
				}
				else {
					if( ii > NG) ii = NG; if(jj > NG) jj = NG;
					if( ii <= 0) ii = 1; if(jj <= 0) jj = 1;
					PG          = paste("PointGrid[",ii,",",jj,"]", sep="");
					ActualPoint = get( PG , env=envir, inherits=TRUE);
					NM_Cur      = paste("NumPoints[",ActualPoint,"][[1]][5]", sep="");
					if( get(NM_Cur, env=envir, inherits=TRUE) != 0 ){
						NumClass                  = as.integer(get(NM_Cur, env=envir, inherits=TRUE));
						ScoreNeighbours[NumClass] <- ScoreNeighbours[NumClass] + 1 ;
						#cat("N=", i, "\t", "Class=", NumClass, "\t", file=FileOutput);
						#cat("ScoreNeighbours=", ScoreNeighbours[ NumClass ], "\n", file=FileOutput);
					}#endif
				}#endif
			} #endforjj
	} #endforii

	for(IndNumClass in 1:NumCluster ){
		#cat("N=", IndNumClass, "\t", "length(ListClusters)=", NumCluster, "\t", file=FileOutput);
		#cat("i=", i, "\t", "ScoreNeighbours=", ScoreNeighbours[ IndNumClass ], "\n", file=FileOutput);
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
assign("nu",    Nu, env=envir, inherits = TRUE);
assign("q",     q,  env=envir, inherits = TRUE); 
assign("Ngrid", G,  env=envir, inherits = TRUE);
assign("Knn",   K,  env=envir, inherits = TRUE);
GlobalTime = 0;
FileOutput <- file(fileName, "w");		# open an output file connection
AvPrecision = 0;
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=envir, inherits=TRUE)+1), ncol = (get("Ngrid", env=envir, inherits=TRUE)+1), byrow = FALSE, dimnames = NULL), env=envir, inherits= TRUE);
assign("NumPoints", array( list(), (get("Ngrid", env=envir, inherits=TRUE)+1)*(get("Ngrid", env=envir, inherits=TRUE)+1)), env=envir, inherits= TRUE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=envir, inherits= TRUE);

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
	assign("MaxValA", 1/(get("nu", env=envir, inherits=TRUE)*nrow(Matrice$Mat)), env=envir, inherits= TRUE ) ;  
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
	assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=envir, inherits= TRUE);  # computation a cluster radius
	assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=envir, inherits= TRUE);  	          # computation of delta R
	Alert("", "ok", "\n");

	if(MetLab == 1 ) {
		Alert("", "grid labeling...", ""); 
		Alert("\t\t", "grid clustering...", "\n");
		NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, cx, cy, WYA=WVectorsYA$A);	# clusters assignment
		Alert("\t\t", "ok", "\n");
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

options(digits=3); cat ("PrecisionGlobale=", AvPrecision/10, "% \t", "\n", sep="", file=FileOutput);
print("time consuming per process"); print(GlobalTime/10);	# output time consuming
close(FileOutput);

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
assign("nu",    Nu, env=envir, inherits = TRUE);
assign("q",     q,  env=envir, inherits = TRUE); 
assign("Ngrid", G,  env=envir, inherits = TRUE);
assign("Knn",   K,  env=envir, inherits = TRUE);
GlobalTime = 0;
FileOutput <- file(fileName, "w");		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=envir, inherits=TRUE)+1), ncol = (get("Ngrid", env=envir, inherits=TRUE)+1), byrow = FALSE, dimnames = NULL), env=envir, inherits= TRUE);
assign("NumPoints", array( list(), (get("Ngrid", env=envir, inherits=TRUE)+1)*(get("Ngrid", env=envir, inherits=TRUE)+1)), env=envir, inherits= TRUE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=envir, inherits= TRUE);

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
assign("MaxValA", 1/(get("nu", env=envir, inherits=TRUE)*nrow(Matrice$Mat)), env=envir, inherits= TRUE ) ; 
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
assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=envir, inherits= TRUE);  # computation a cluster radius
assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=envir, inherits= TRUE);  	          # computation of delta R
Alert("", "ok", "\n");

if(MetLab == 1 ) {
	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, cx, cy, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=MatriceEval, Points=NumPoints, Grid=PointGrid, NumCluster=NumberCluster, Cx=cx, Cy=cy, Knn=K);  # cluster assignment
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
options(digits=3); cat ("PrecisionGlobale=", Precision, "% \t", "\n", sep="", file=FileOutput);
print("time consuming"); print(GlobalTime);	# output time consuming
close(FileOutput);

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
assign("nu",    Nu, env=envir, inherits = TRUE);
assign("q",     q,  env=envir, inherits = TRUE); 
assign("Ngrid", G,  env=envir, inherits = TRUE);
assign("Knn",   K,  env=envir, inherits = TRUE);
GlobalTime = 0;
FileOutput <- file(fileName, "w");		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=envir, inherits=TRUE)+1), ncol = (get("Ngrid", env=envir, inherits=TRUE)+1), byrow = FALSE, dimnames = NULL), env=envir, inherits= TRUE);
assign("NumPoints", array( list(), (get("Ngrid", env=envir, inherits=TRUE)+1)*(get("Ngrid", env=envir, inherits=TRUE)+1)), env=envir, inherits= TRUE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=envir, inherits= TRUE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

Alert("", "loading matrix...", "\t\t");
Matrice    <- chargeMatrix(DataName=DName, PathIn=fileIn);		# data matrix structure loading
MM = Matrice$Mat
Alert("", "ok", "\n");

TimeNow        <- proc.time();			# catch time reference
FileOutputTime <- file("timesvc.txt", "w");

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
	assign("MaxValA", 1/(get("nu", env=envir, inherits=TRUE)*nrow(Matrice$Mat)), env=envir, inherits= TRUE ) ; 
	MatriceK   <- calcKernelMatrix(Matrice$Mat);						# kernel matrix computation
	Alert("", "ok", "\n");

	# lagrange multiplier computation
	Alert("", "lagrange coefficients...", "\t");
	WVectorsYA <- CalcWcluster(MatriceKern=MatriceK)
	Alert("", "ok", "\n");

	Alert("", "radius computation...", "\t\t"); 
	assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=envir, inherits= TRUE);  # computation a cluster radius
	assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=envir, inherits= TRUE);  	          # computation of delta R
	Alert("", "ok", "\n");

	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, 1, 2, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=Matrice$Mat, NumCluster=NumberCluster, Cx=1, Cy=2, Knn=K);  # cluster assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "evaluation...", "");
	MisClass      = Evaluation(Mat=Matrice$Mat, NBClass=NumberCluster, Cx=1, Cy=2, ClassPoints=ClassPoints);					# evaluation	
	Alert("\t\t\t\t", "ok", "\n");

	GlobalTime = ( proc.time() - TimeNow ) ;
	options(digits=3); 
	cat("NLine=", NLine, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=envir, inherits=TRUE), '\n', file=FileOutputTime);	# output time consuming
	cat("NLine=", NLine, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=envir, inherits=TRUE), '\n');	# output time consuming

} #Endfori

invisible(gc());					# freeing memory
close(FileOutput);
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
assign("nu",    Nu, env=envir, inherits = TRUE);
assign("q",     q,  env=envir, inherits = TRUE); 
assign("Ngrid", G,  env=envir, inherits = TRUE);
assign("Knn",   K,  env=envir, inherits = TRUE);
GlobalTime = 0;
FileOutput <- file( paste(fileIn, fileName, sep="") , "w");		# open an output file connection
#Data Grid structure
assign("PointGrid", matrix(data = NA, nrow = (get("Ngrid", env=envir, inherits=TRUE)+1), ncol = (get("Ngrid", env=envir, inherits=TRUE)+1), byrow = FALSE, dimnames = NULL), env=envir, inherits= TRUE);
assign("NumPoints", array( list(), (get("Ngrid", env=envir, inherits=TRUE)+1)*(get("Ngrid", env=envir, inherits=TRUE)+1)), env=envir, inherits= TRUE);
assign("TabW",seq(length=MaxIter,from=MinW,by=0), env=envir, inherits= TRUE);

if( Cx != "" && Cy != "" ) { cx <- Cx ; cy <- Cy; }

FileOutputTime <- file( paste(fileIn, "timesvc.txt", sep="") , "w");

ListNu = c(0.1, 0.5, 1);
ListQ  = c(0.5, 1, 10, 100);

for(IndNu in 1:length(ListNu) ) {
	
   Nu = ListNu[IndNu];
   assign("nu",    Nu, env=envir, inherits = TRUE);
	
   for(IndQ in 1:length(ListQ)) {
	
	q = ListQ[IndQ];
	assign("q",     q,  env=envir, inherits = TRUE); 
	
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
	assign("MaxValA", 1/(get("nu", env=envir, inherits=TRUE)*nrow(Matrice$Mat)), env=envir, inherits= TRUE ) ; 
	MatriceK   <- calcKernelMatrix(Matrice$Mat);						# kernel matrix computation
	Alert("", "ok", "\n");

	# lagrange multiplier computation
	Alert("", "lagrange coefficients...", "\t");
	WVectorsYA <- CalcWcluster(MatriceKern=MatriceK)
	Alert("", "ok", "\n");

	Alert("", "radius computation...", "\t\t"); 
	assign("RadiusC", RadiusCluster(VA=WVectorsYA$A, MatrixK=MatriceK) , env=envir, inherits= TRUE);  # computation a cluster radius
	assign("r", SmallR(VA=WVectorsYA$A, MatK=MatriceK) , env=envir, inherits= TRUE);  	          # computation of delta R
	Alert("", "ok", "\n");

	Alert("", "grid labeling...", "\n"); 
	Alert("\t\t", "grid clustering...", "");
	NumberCluster = ClusterLabeling(Mat=Matrice$Mat, MatK= MatriceK, 1, 2, WYA=WVectorsYA$A);	# clusters assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "match grid...", ""); 
	ClassPoints   = MatchGridPoint(Mat=Matrice$Mat, NumCluster=NumberCluster, Cx=1, Cy=2, Knn=K);  # cluster assignment
	Alert("\t\t", "ok", "\n");
	Alert("\t\t", "evaluation...", "");
	MisClass      = Evaluation(Mat=Matrice$Mat, NBClass=NumberCluster, Cx=1, Cy=2, ClassPoints=ClassPoints);					# evaluation	
	Alert("\t\t\t\t", "ok", "\n");

	GlobalTime = ( proc.time() - TimeNow ) ;
	options(digits=3); 
	cat("q=",  get("q", env=envir, inherits=TRUE), '\t', "nu=", get("nu", env=envir, inherits=TRUE), '\t', "#cluster", NumberCluster, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=envir, inherits=TRUE), '\n', file=FileOutputTime);	# output time consuming
	cat("q=",  get("q", env=envir, inherits=TRUE), '\t' , "nu=", get("nu", env=envir, inherits=TRUE), '\t', "#cluster", NumberCluster, '\t', "T=", GlobalTime[1], '\t', "P=", get("Precision", env=envir, inherits=TRUE), '\n');	# output time consuming
	
	invisible(gc());					# freeing memory

   } #EndforIndQ

} #EndforIndNu

close(FileOutput);
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
CluOutput <- file( paste(pathOut, DName, "_clu.txt", sep="") , "w");

#sorting by cluster index
SortedClassPoints = sort(CPoints, index = TRUE);

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
##End of program##
########################################################################