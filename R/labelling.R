## Cluster Labelling functions.
##
## created  10.04.09 nicolas turenne
## updated   

setGeneric("Labelling",function(x, ...) standardGeneric("Labelling"))

Labelling.compute <- function (  x ,
	  MetLab	= 1,
	  MatriceKernel = NULL,
	  MatriceK      = NULL,
	  pp            = NULL,
	  Nu		= 1,
          G		= 1,
          q		= 1,
          ncol		= 1,
          nlin		= 1,
          RadiusC	= 2,
          r		= 2,
          KernChoice	= 0.01,
	  NbClassInData	= 0.01  )
{ 

MxX = x@MinMaxXY[1];
MnX = x@MinMaxXY[2];
MxY = x@MinMaxXY[3];
MnY = x@MinMaxXY[4];

NumPoints = c();
NbCluster = 0;
ClassPoints = vector();

if( MetLab == 1 ) {
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
			as.vector(x@WVectorsYA$A) , 
			iNumPoints = numeric( G * G * 7 ) )$iNumPoints;

	NbCluster =   .C("NbCluster_C",     
			as.vector(NumPoints), 
			as.integer(G), 
			iNCluster = integer(1) )$iNCluster;

	Alert("\t", "ok", "\n");
	cat("\t\t\t\t", "NbCluster ", NbCluster, "\n");					
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
			as.integer(x@KNN), 
			as.integer(0), 
			as.integer(1), 
			as.vector(NumPoints) , 
			iClassPoints = numeric(nlin) )$iClassPoints;
	
	Alert("\t\t", "ok", "\n");

}
else if( MetLab == 2 ){
	Alert("", "mst labeling/eval...", "\n"); 
	ClassPoints = MST_labelling(x=x, MatriceKernel= MatriceK);
	Alert("\t\t\t\t", "ok", "\n");
}
else if( MetLab == 3 ){
	Alert("", "knn labeling/eval...", "\n");
	ClassPoints    = KNN_labelling(x=x, MatriceKernel= MatriceK);
	Alert("\t\t\t\t", "ok", "\n");
}

    return( new( "Labelling", ClassPoints = ClassPoints, NumPoints = NumPoints ) )
}
setMethod("Labelling",signature(x="findModelCluster"), Labelling.compute)


########################################################################
# Compute adjacency between two points 
#
########################################################################

# Usage:
#   AdjacencyPP(x=x, MatriceKernel=MatriceKernel, Vec1=v1, Vec2=v2);
#

setGeneric("AdjacencyPP", function(x, MatriceKernel, Vec1, Vec2) standardGeneric("AdjacencyPP"))

setMethod("AdjacencyPP", signature(x = "findModelCluster"),
function(x,
    MatriceKernel=matrix(),            # matrix of kernel product
    Vec1=vector() ,                   # first vector
    Vec2=vector()                     # second vector
    ) {

Mat  = x@Data ;
WYA  = x@WVectorsYA$A ;                      # Lagrange coefficients

SvcEnv <- get("SvcEnv", env=globalenv(), inherits = FALSE);
adj_flag = 1; # unless a point on the path exits the sphere - pair is adjacent
		
interval = 0.0;
while( interval < 1 && adj_flag ){

	z = Vec1 + interval * (Vec2 - Vec1);	    
	interval = interval + 0.3;
	R  = RadiusPoint(Vec=z, VA=WYA, Mat=Mat, MatK=MatriceKernel);	
	RC = get("RadiusC", env=SvcEnv, inherits = FALSE);
	r  = get("r", env=SvcEnv, inherits = FALSE);
	if(  (RC*RC + r) <  R*R ){
		adj_flag = 0;
		interval = 1;
	} #finif
		
} #finwhile interval

return ( adj_flag );
})

########################################################################
# Compute adjacency matrix 
#
########################################################################

# Usage:
#   Adjacency(x=x, MatriceKernel=MatriceKernel);
#

setGeneric("Adjacency", function(x, MatriceKernel) standardGeneric("Adjacency"))

setMethod("Adjacency", signature(x = "findModelCluster"),
function(x,
    MatriceKernel=matrix()            # matrix of kernel product
    ) {

Mat  = x@Data ;
WYA  = x@WVectorsYA$A ;                      # Lagrange coefficients

N          = nrow(MatriceKernel) * 0.1 ;
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
		
		adj_flag = AdjacencyPP(x, MatriceKernel, Vec1=Mat[IndTarget,1:(ncol(Mat)-1)], Vec2=Mat[j,1:(ncol(Mat)-1)]);

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
})

########################################################################
# Compute MST proximity clustering 
#
########################################################################

# Usage:
#   MST_labelling(Mat=matrice, MatK= matriceK, WYA=wya, Cx=cx, Cy=cy);
#

setGeneric("MST_labelling", function(x, MatriceKernel) standardGeneric("MST_labelling"))

setMethod("MST_labelling", signature(x = "findModelCluster"),
function(x,
    MatriceKernel=matrix()            # matrix of kernel product
    ) {

Mat	= x@Data ;		# data matrix
WYA	= x@WVectorsYA$A ;    # Lagrange coefficients
Cx	= x@Cx;   		# choice of x component to display
Cy	= x@Cy; 		# choice of y component to display

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
			adj_flag = AdjacencyPP(x, MatriceKernel, Vec1=Mat[i,1:(ncol(Mat)-1)], Vec2=Mat[j,1:(ncol(Mat)-1)]);
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
})

########################################################################
# Compute KNN proximity clustering 
#
########################################################################

# Usage:
#   KNN_labelling(x=x, MatriceKernel= matriceK);
#

setGeneric("KNN_labelling", function(x, MatriceKernel) standardGeneric("KNN_labelling"))

setMethod("KNN_labelling", signature(x = "findModelCluster"),
function(x,
    MatriceKernel=matrix()            # matrix of kernel product
    ) {

Mat	= x@Data ;		# data matrix
WYA	= x@WVectorsYA$A ;    # Lagrange coefficients
Cx	= x@Cx;   		# choice of x component to display
Cy	= x@Cy; 		# choice of y component to display

Alert("\t\t", "begin knn...", "");

MatKNN <- knearneigh(as.matrix(Mat[,c(Cx,Cy)]), k=4);

MatAdj01  = matrix(0, ncol=nrow(Mat), nrow=nrow(Mat));
IndRowMatKNN = IndColMatKNN = 0;
for(IndRowMatKNN in 1:nrow(MatKNN$nn) ){
	MatAdj01[IndRowMatKNN,IndRowMatKNN] = 1;
	for(IndColMatKNN in 1:ncol(MatKNN$nn) ){ 
		i = as.numeric(IndRowMatKNN);
		j = as.numeric(MatKNN$nn[i,IndColMatKNN]);
		adj_flag = AdjacencyPP(x, MatriceKernel, Vec1=Mat[i,1:(ncol(Mat)-1)], Vec2=Mat[j,1:(ncol(Mat)-1)]);
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
})

setGeneric("MineLineMat", function(ListCluster, NumRow, Mat ) standardGeneric("MineLineMat"))

setMethod("MineLineMat", signature(ListCluster = "vector"),
function(ListCluster,   # matrix adjacency
    NumRow=1,            # num row
    Mat=NULL            # matrix of kernel product
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
} )


########################################################################
##End of Class##
########################################################################
