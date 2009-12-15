## Display functions.
##
## created  10.04.09 nicolas turenne
## updated   

########################################################################
# Afffiche Data
#
########################################################################

# Usage:
#   plot(fmc)
#

setMethod("plot", signature(x = "findSvcModel", y = "missing"),
function(x, data = NULL, Multi = FALSE , slice = list(), ...) {

Mat		= x@Matrice$Mat		# data matrix
WYA		= x@lagrangeCoeff	# vectors classs/ coefficients
Cx		= x@Cx 			# choice of x component to display
Cy		= x@Cy	 		# choice of y component to display
ListMis		= x@MisClass		# list of misclassified data
NumPoints	= x@NumPoints		# list of grid membership to clusters
ClassPoints	= x@ClassPoints		# list of data points to clusters
MinMaxXY	= x@MinMaxXY		# minmax  vector
nu		= x@Nu			# minmax  vector
AroundNullVA	= x@AroundNullVA	
sizegrid	= x@SizeGrid

# framing display
if( Multi == TRUE ) par(mfrow = c(2, 2));					

# plot the data matrix
#cat("matrix", "\n", file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(as.matrix(Mat), sep="\t", file=get("FileOutput", env=SvcEnv, inherits = FALSE));					

# plot Y class in the case of an svm
if( WYA$Y[1] != "" && Multi)
	plot( 1:length(WYA$Y), WYA$Y , xlab="point", ylab="Classe");

# plot Lagrange parameters
if( WYA$A[1] != "" && Multi)
	plot( 1:length(WYA$A), WYA$A , xlab="point", ylab="Coefficient de Lagrange");

# plots data numbers being SV, make SV in red
if(  Multi )
	plot(Mat[,1],Mat[,ncol(Mat)], xlab="data points", ylab="classes in data matrix")

MxValA = 1 / ( nu * length(WYA$A) ); # valeur quasi optimale pour 2*1/N 

SV =   .C("ListSVPoints_C",     
		as.vector(WYA$A) ,  
		as.integer( nrow(Mat) ), 
		as.numeric( MxValA ), 
		as.numeric( AroundNullVA ), 
		iListPoints = numeric( nrow(Mat) ) )$iListPoints;

#cat("list des SV", '\n', file=get("FileOutput", env=SvcEnv, inherits = FALSE)); write(t(SV), file=get("FileOutput", env=SvcEnv, inherits = FALSE)); 
if(  Multi )
for(i in 1:length(SV) ) {				
		ISV = SV[i];
		if( !is.na(ISV) )
			points(Mat[ISV,1], Mat[ISV,ncol(Mat)], pch = 24, col = "red", bg = "yellow", cex = 1)
} #fin for i

Grid(fmc=x, ListSV=SV);

})

########################################################################
# Compute Grid 
#
########################################################################

# Usage:
#   Grid( fmc=ret, ListSV=SV );
#

setGeneric("Grid", function(fmc, ListSV) standardGeneric("Grid"))

setMethod("Grid", signature(fmc = "findSvcModel"),
function(fmc=new("findSvcModel"), ListSV = NULL) {

ListMis		= fmc@MisClass		# list of misclassified data
Mat		= fmc@Matrice$Mat	# data matrix
NumPoints	= fmc@NumPoints		# list of grid membership to clusters
MinMaxXY	= fmc@MinMaxXY		# minmax  vector
ClassPoints	= fmc@ClassPoints
Ngrid		= fmc@SizeGrid
Cx		= fmc@Cx
Cy		= fmc@Cy
Att		= fmc@Matrice$Att	# data matrix

NRow	= nrow(Mat); NCol = ncol(Mat) ;
N	= Ngrid;

ListX		<- c(); # we make the list of X values	
ListY		<- c(); # we make the list of Y values
ListX_clu	<- c(); # we make the list of X values	
ListY_clu	<- c(); # we make the list of Y values
ListX_bg	<- c(); # we make the list of X values	
ListY_bg	<- c(); # we make the list of Y values
							 
MaxX = MinMaxXY[1];
MinX = MinMaxXY[2];
MaxY = MinMaxXY[3];
MinY = MinMaxXY[4];

if( Cx == 0 ) xname = "1st Principal Component" else xname = as.character( Att[[1]][Cx] )
if( Cy == 0 ) yname = "2nd Principal Component" else yname = as.character( Att[[1]][Cy] )

#zoom+
plot(NA,NA,xlim=c(MinX,MaxX), ylim=c(MinY,MaxY), xlab=xname, ylab=yname, main = "SVM clustering plot") ; # setting up co
							
# plots grid
for(i in 1:(max(NumPoints)+1) ){
	x		= NumPoints[ (i-1)*6 +3]; 
	y		= NumPoints[ (i-1)*6 +4];
	c		= NumPoints[ (i-1)*6 +6];

	if(   c != 0  ){
		ListX_clu	= c(ListX_clu, ( MaxX - MinX )*( (x) / (N-1) ) + MinX );
		ListY_clu	= c(ListY_clu, ( MaxY - MinY )*( (y) / (N-1) ) + MinY );
	}
	else {
		ListX_bg	= c(ListX_bg,  ( MaxX - MinX )*( (x) / (N-1) ) + MinX ) ;
		ListY_bg	= c(ListY_bg,  ( MaxY - MinY )*( (y) / (N-1) ) + MinY ) ;
	}
} #endfori

points(ListX_bg, ListY_bg, pch = 21, col = "blue", bg = "black", cex = 1);
points(ListX_clu, ListY_clu, pch = 24, col = "yellow", bg = "yellow", cex = 1);

# plots data points
for(i in 1:NRow ){

	if( ClassPoints[i] != 0 ){
		points(Mat[i,1], Mat[i,2], pch = 21, col = "red", bg = "red", cex = 1);
	}
	#else
	#	points(Mat[i,1], Mat[i,2], pch = 24, col = "green", bg = "yellow", cex = 1)
}#endfori

# plots data numbers being SV, make SV in red
#for(i in 1:length(ListSV) ) {				
#		ISV = ListSV[i];
#		if( !is.na(ISV) )
#			points(Mat[ISV,1], Mat[ISV,2], pch = 24, col = "red", bg = "yellow", cex = 1)
#} #fin for i

# plots misclassified data 
if( sum( ListMis[]>0 ) != 0 )
for(i in 1: sum( ListMis[]>0 ) ){
	#cat("ListMis[i]", ListMis[i], "\t", "x", Mat[ ListMis[i]+1 ,1], "\t", "y", Mat[ ListMis[i]+1 ,2],  "\n");
	if( ListMis[i] != 0 ){
		ListX	= c( ListX, Mat[ ListMis[i]+1 ,1] );
		ListY	= c( ListY, Mat[ ListMis[i]+1 ,2] );
	} #endif
}#endfori
points(ListX, ListY, pch = 24, col = "green", bg = "green", cex = 1);

})


########################################################################
# Exporting clusters in a text file 
#
########################################################################

# Usage:
#   ExportClusters(fmc = findMC, NameFile="nf" );
#

setGeneric("ExportClusters", function(fmc,NameFile) standardGeneric("ExportClusters"))

setMethod("ExportClusters", signature(fmc = "findSvcModel"),
function(fmc=new("findSvcModel"), NameFile="nf") {

MatriceVar	= fmc@Matrice$Var	# variables data matrix name 
dataFrame	= NameFile		# prefix name of data
pathOut		= tempdir()		# output directory
CPoints		= fmc@ClassPoints

#opening output stream
path = file.path(pathOut, paste(dataFrame, "_clu.txt", sep="") );
CluOutput <- file( path , "w");

#sorting by cluster index
SortedClassPoints = sort(CPoints, method = "sh", index.return = TRUE);

#output points class
for(i in 1:length(MatriceVar) ){ 
 if( i == 1 || SortedClassPoints$x[i] != SortedClassPoints$x[i-1] )
	cat('\n', "cluster", '\t', SortedClassPoints$x[i], '\n', file=CluOutput);
 cat("item", '\t', as.character(MatriceVar[ SortedClassPoints$ix[i] ]), '\n', file=CluOutput);

} #endfori

cat("fichier ", '\t', path, " created -- use read.table(\"path//filename\") to load ");

close(CluOutput)
})


########################################################################
# Summary clusters content 
#
########################################################################

# Usage:
#   findSvcModel.summary(x=ret);
#

setGeneric("findSvcModel.summary", function(x) standardGeneric("findSvcModel.summary"))

setMethod("findSvcModel.summary", signature(x = "findSvcModel"),
function(x=new("findSvcModel")) {

Matrice		= x@Data		# original data matrix
CPoints		= x@ClassPoints
ListAtt		= x@Matrice$Att	# original data matrix

NbClusters = max(CPoints);
NbAtt	   = length( t(ListAtt) );

cat("CLUSTER ID ", '\t', "SIZE", '\n');
for(i in 0:NbClusters ){
	SizeCluster = length( x@ClassPoints[ x@ClassPoints[] ==i ] ) ;
	if( i == 0 ) {
		cat("isolated", '\t', SizeCluster, '\n');
	} 
	else {
		cat( "  ", i , '\t\t', SizeCluster, '\n');
	} #endif
} #endfor

cat("average attributes per cluster", '\n');
print( c(t(x@Matrice$Att)) );
for(i in 1:NbClusters ){
	MatriceCluster = x@Data[  x@ClassPoints[] ==i  , ] ;
	if( !is.null(nrow(MatriceCluster)) ) { cat( "cluster ", i, "\t", mean( as.data.frame(MatriceCluster[,1:NbAtt]) ) , '\n\n'); }
	else                           { cat( "cluster ", i, "\t", MatriceCluster  , '\n\n'); }
} #endfor


})

########################################################################
# Accessing to a cluster by its id 
#
########################################################################

# Usage:
#   GetClusterID(fmc=ret, Id=1);
#

setGeneric("GetClusterID", function(fmc, Id=1) standardGeneric("GetClusterID"))

setMethod("GetClusterID", signature(fmc = "findSvcModel"),
function(fmc=new("findSvcModel"), Id=1) {

if( Id < 0 || Id > max(fmc@ClassPoints) || (Id %% 1) != 0 ){
	print("Id is not valid"); 
	return (1);
} #endif

MatriceVar	= fmc@Matrice$Var	# variables data matrix name 

#output points class
MatriceCluster = MatriceVar[  fmc@ClassPoints[] == Id   ] ;
cat( "cluster ", Id, '\n'); 
cat(  MatriceCluster , sep='\n'); 
cat(  '\n'); 
 
})

########################################################################
# Show all clusters
#
########################################################################

# Usage:
#   ShowClusters(fmc=ret);
#

setGeneric("ShowClusters", function(fmc) standardGeneric("ShowClusters"))

setMethod("ShowClusters", signature(fmc = "findSvcModel"),
function(fmc=new("findSvcModel")) {

MatriceVar	= fmc@Matrice$Var	# variables data matrix name 
CPoints		= fmc@ClassPoints

#sorting by cluster index
SortedClassPoints = sort(CPoints, method = "sh", index.return = TRUE);

#output points class
for(i in 1:length(MatriceVar) ){ 
cat("i", i, "\n");
 if( i == 1 || SortedClassPoints$x[i] != SortedClassPoints$x[i-1] )
	cat('\n', "cluster", '\t', SortedClassPoints$x[i], '\n' );
 cat("item", '\t', as.character(MatriceVar[ SortedClassPoints$ix[i] ]), '\n' );

} #endfori


})

########################################################################
# Get the clusters to which an items belongs 
#
########################################################################

# Usage:
#   GetClustersTerm(fmc=ret, term="home");
#

setGeneric("GetClustersTerm", function(fmc,term) standardGeneric("GetClustersTerm"))

setMethod("GetClustersTerm", signature(fmc = "findSvcModel"),
function(fmc=new("findSvcModel"), term="home") {

CPoints		= fmc@ClassPoints
ListVar		= fmc@Matrice$Var	# variables data matrix name 

print(ListVar)
IndiceVectorTerm <- grep(term, ListVar)  # indices
print(IndiceVectorTerm)
ListCluster = unique( CPoints[ IndiceVectorTerm ] )
print(ListCluster)

for(k in 1:length(ListCluster) ) GetClusterID(fmc, ListCluster[k] );

})



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
##End of Class##
########################################################################
