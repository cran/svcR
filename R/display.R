## Display functions.
##
## created  10.04.09 nicolas turenne
## updated   

########################################################################
# Afffiche Data
#
########################################################################

# Usage:
#   DisplayData(Mat=matrice,MatK=matriceK,WYA=wya,Cx=1,Cy=2,ListMis=MisClas)
#

setMethod("plot", signature(x = "findModelCluster", y = "missing"),
function(x, data = NULL, Multi = FALSE , slice = list(), ...) {

Mat		= x@Matrice$Mat		# data matrix
WYA		= x@WVectorsYA		# vectors classs/ coefficients
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
#   Grid(Mat=matrice, MatK= matriceK, WYA=va, ListSV=SV, Cx=1,Cy=2);
#

setGeneric("Grid", function(fmc, ListSV) standardGeneric("Grid"))

setMethod("Grid", signature(fmc = "findModelCluster"),
function(fmc = NULL, ListSV = NULL) {

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
}#finfori

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
}#finfori
points(ListX, ListY, pch = 24, col = "green", bg = "green", cex = 1);

})


########################################################################
# Exporting clusters in a text file 
#
########################################################################

# Usage:
#   ExportClusters(fmc = findMC);
#

setGeneric("ExportClusters", function(fmc) standardGeneric("ExportClusters"))

setMethod("ExportClusters", signature(fmc = "findModelCluster"),
function(fmc = NULL) {

MatriceVar	= fmc@Matrice$Var	# variables data matrix name 
DName		= fmc@DName		# prefix name of data
pathOut		= tempdir()		# output directory
CPoints		= fmc@ClassPoints
Ngrid		= fmc@SizeGrid		# size grid

#opening output stream
path = file.path(pathOut, paste(DName, "_clu.txt", sep="") );
CluOutput <- file( path , "w");

#sorting by cluster index
SortedClassPoints = sort(CPoints, method = "sh", index.return = TRUE);

#output points class
for(i in 1:nrow(MatriceVar) ){
 if( i == 1 || SortedClassPoints$x[i] != SortedClassPoints$x[i-1] )
	cat('\n', "cluster", '\t', SortedClassPoints$x[i], '\n', file=CluOutput);
 cat("item", '\t', as.character(MatriceVar[ SortedClassPoints$ix[i] , 1 ]), '\n', file=CluOutput);

} #finfori

cat("fichier ", '\t', path, " created -- use read.table(\"path//filename\") to load ");

close(CluOutput)
})


########################################################################
# Summary clusters content 
#
########################################################################

# Usage:
#   Summary(fmc=ret);
#

setGeneric("Summary", function(fmc) standardGeneric("Summary"))

setMethod("Summary", signature(fmc = "findModelCluster"),
function(fmc = NULL) {

Matrice		= fmc@Data		# original data matrix
CPoints		= fmc@ClassPoints
ListAtt		= fmc@Matrice$Att	# original data matrix

NbClusters = max(CPoints);
NbAtt	   = length( t(ListAtt) );

cat("CLUSTER ID ", '\t', "SIZE", '\n');
for(i in 0:NbClusters ){
	SizeCluster = length( fmc@ClassPoints[ fmc@ClassPoints[] ==i ] ) ;
	if( i == 0 ) {
		cat("isolated", '\t', SizeCluster, '\n');
	} 
	else {
		cat( "  ", i , '\t\t', SizeCluster, '\n');
	} #endif
} #endfor

cat("average attributes per cluster", '\n');
print( c(t(fmc@Matrice$Att)) );
for(i in 1:NbClusters ){
	MatriceCluster = fmc@Data[ fmc@ClassPoints[ fmc@ClassPoints[] ==i ], ] ;
	cat( "  ", mean( as.data.frame(MatriceCluster[,1:NbAtt]) ) , '\n');
} #endfor

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
