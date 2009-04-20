## Sampling & Evaluation functions.
##
## created  10.04.09 nicolas turenne
## updated   


setGeneric("ClusterEval",function(x, ...) standardGeneric("ClusterEval"))

ClusterEval.Pvalue <- function (  x  )
{ 

nlin = nrow(x@Matrice$Mat);
ncol = ncol(x@Matrice$Mat);
pp   = c();
for(i in 1:nlin ){	# we fill the full matrix
	pp = c(pp, as.vector(x@Matrice$Mat[i,]) );
}#finfori

NBClass		= max( x@ClassPoints );
NbClassInData	= max( x@Matrice$Mat[,ncol] );
P               = 0;
MisClass	= c();

#cat("NBClass", NBClass, '\t', "NbClassInData", NbClassInData, "\n");

if( NbClassInData > 0 ) {
		Alert("evaluation...", "");

		MisClass =   .C("Evaluation_C",     
			as.vector(pp) ,  
			as.integer(nlin), 
			as.integer(ncol), 
			as.integer(NBClass), 
			as.vector(x@ClassPoints) , 
			iMisClass = numeric( (NBClass+1)*nlin ) )$iMisClass;

		Alert("\t\t", "ok", "\n");
		P= (100*(nlin-sum(MisClass[]>0))/nlin) ;
		cat("\t\t\t", "Precision ", P, "\n");	
}

rm(pp,MisClass);
invisible(gc());					# freeing memory

return ( new("ClusterEval",  Precision= c(P) ) )
}
setMethod("ClusterEval",signature(x="findModelCluster"), ClusterEval.Pvalue)


########################################################################
# Compute Data Sampling . Sample Iris Data in Ten training/test parts
#
########################################################################

# Usage:
#   DataSampling(DatMat=matrice);
#

setGeneric("DataSampling", function(DatMat) standardGeneric("DataSampling"))

setMethod("DataSampling", signature(DatMat = "matrix"),
function (
    DatMat=matrix() 			# data matrix
    ) {
ListMatrixLearnTest                  <- list();
ListMatrixLearnTest$ListMatLearn     <- list(); 
ListMatrixLearnTest$ListMatTest      <- list();
ListMatrixLearnTest$ListMatLearnTest <- list();
ListMatrixLearnTest$ListMatEval      <- list();

ListMatrixLearnTest$ListMatLearn[[1]]     <- DatMat[ c( c(5:44) , c((5+50):(44+50)) , c((5+100):(44+100)) ),];
ListMatrixLearnTest$ListMatTest[[1]]      <- DatMat[ c( c(1:4,(1+50):(4+50),(1+100):(4+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[1]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[1]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[2]]     <- DatMat[ c( c(1:4,9:44) , c((1+50):(4+50),(9+50):(44+50)), c((1+100):(4+100),(9+100):(44+100))) , ];
ListMatrixLearnTest$ListMatTest[[2]]      <- DatMat[ c( c(5:8,(5+50):(8+50),(5+100):(8+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[2]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[2]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[3]]     <- DatMat[ c( c(1:8,13:44) , c((1+50):(8+50),(13+50):(44+50)) , c((1+100):(8+100),(13+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[3]]      <- DatMat[ c( c(9:12,(9+50):(12+50),(9+100):(12+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[3]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[3]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[4]]     <- DatMat[ c( c(1:12,17:44) , c((1+50):(12+50),(17+50):(44+50)) , c((1+100):(12+100),(17+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[4]]      <- DatMat[ c( c(13:16,(13+50):(16+50),(13+100):(16+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[4]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[4]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[5]]     <- DatMat[ c( c(1:17,21:44) , c((1+50):(17+50),(21+50):(44+50)) , c((1+100):(17+100),(21+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[5]]      <- DatMat[ c( c(17:20,(17+50):(20+50),(17+100):(20+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[5]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[5]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[6]]     <- DatMat[ c( c(1:20,25:44) , c((1+50):(20+50),(25+50):(44+50)) , c((1+100):(20+100),(25+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[6]]      <- DatMat[ c( c(21:24,(21+50):(24+50),(21+100):(24+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[6]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[6]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[7]]     <- DatMat[ c( c(1:24,30:44) , c((1+50):(24+50),(30+50):(44+50)) , c((1+100):(24+100),(30+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[7]]      <- DatMat[ c( c(25:29,(25+50):(29+50),(25+100):(29+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[7]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[7]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[8]]     <- DatMat[ c( c(1:29,34:44) , c((1+50):(29+50),(34+50):(44+50)) , c((1+100):(29+100),(34+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[8]]      <- DatMat[ c( c(30:33,(30+50):(33+50),(30+100):(33+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[8]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[8]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[9]]     <- DatMat[ c( c(1:33,38:44) , c((1+50):(33+50),(38+50):(44+50)) , c((1+100):(33+100),(38+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[9]]      <- DatMat[ c( c(34:37,(34+50):(37+50),(34+100):(37+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[9]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[9]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

ListMatrixLearnTest$ListMatLearn[[10]]     <- DatMat[ c( c(1:37,42:44) , c((1+50):(37+50),(42+50):(44+50)) , c((1+100):(37+100),(42+100):(44+100)) ) ,];
ListMatrixLearnTest$ListMatTest[[10]]      <- DatMat[ c( c(38:41,(38+50):(41+50),(38+100):(41+100)) ) , ];
ListMatrixLearnTest$ListMatLearnTest[[10]] <- DatMat[ c( c(1:44,(1+50):(44+50),(1+100):(44+100)) ) , ];
ListMatrixLearnTest$ListMatEval[[10]]      <- DatMat[ c( c((45:50),(45+50):(50+50),(45+100):(50+100)) ) , ];

return(ListMatrixLearnTest);
})

########################################################################
# Compute Test Evaluation 
#
########################################################################

# Usage:
#   ListP = ClusterEval.crossval( x );
#

ClusterEval.crossval <- function (  x  )
{ 

#parameters init
GlobalTime = 0;
AvPrecision = 0;
ListP = c();

TimeNow <- proc.time();		# catch time reference
	
Alert("", "sampling matrix...", "\t\t");
ListMatrixLearnTest = DataSampling(DatMat=x@Matrice$Mat);
cx = 0; cy = 0;
Alert("", "ok", "\n");

for( IndMat in 1:10 ) {

	Alert("\n", paste("test number...",IndMat) , "\n\n");

	MatriceMat  = as.matrix( ListMatrixLearnTest$ListMatLearn[[IndMat]] );
	MatriceTest = as.matrix( ListMatrixLearnTest$ListMatTest[[IndMat]]  );				

	nlin = nrow(MatriceTest);
	ncol = ncol(MatriceTest);
	pp = c();
	for(i in 1:nlin ){	# we fill the full matrix
		pp = c(pp, as.vector(MatriceTest[i,]) );
	}#finfori

	if( cx != 0 ) {
		MatriceMat = MatriceMat[,c(cx,cy,ncol(MatriceMat))];
	}
	else { 
		if( sign(min(MatriceMat)) < 0 ) {
			MatAdjCoa   = dudi.pca(as.data.frame(MatriceMat[,1:(ncol(MatriceMat)-1)]), scan = FALSE);
		}
		else
			MatAdjCoa   = dudi.coa(as.data.frame(MatriceMat[,1:(ncol(MatriceMat)-1)]), scan = FALSE);
		MatriceMat = as.data.frame( c( MatAdjCoa$li , as.data.frame(MatriceMat[,ncol(MatriceMat)]) ) );
	} #EndIf

	x = findModelCluster.Eval( as.matrix(MatriceMat) );

	NBClass		= max( x@ClassPoints );
	NbClassInData	= max( x@Matrice$Mat[,ncol] );

	ClassPoints =   .C("MatchGridPoint_C",     
			as.vector(pp) , 
			as.vector(x@MatriceK) , 
			as.integer(x@SizeGrid), 
			as.integer(ncol), 
			as.integer(nlin), 
			as.double(x@MinMaxXY[1]), 
			as.double(x@MinMaxXY[2]), 
			as.double(x@MinMaxXY[3]), 
			as.double(x@MinMaxXY[4]), 
			as.integer(NBClass), 
			as.integer(x@KNN), 
			as.integer(0), 
			as.integer(1), 
			as.vector(x@NumPoints) , 
			iClassPoints = numeric(nlin) )$iClassPoints;

	if( NbClassInData > 0 ) {
		Alert("\t\t", "evaluation...", "");

		MisClass = c();

		MisClass =   .C("Evaluation_C",     
			as.vector(pp) ,  
			as.integer(nlin), 
			as.integer(ncol), 
			as.integer(NBClass), 
			as.vector(ClassPoints ) , 
			iMisClass = numeric(nlin) )$iMisClass;

		Alert("\t\t", "ok", "\n");
		P= (100*(nlin-sum(MisClass[]>0))/nlin) ;
		cat("\t\t\t", "Precision ", P, "\n");					
	}
	Alert("\t\t\t\t" ,"ok", "\n");

	AvPrecision = AvPrecision + P;

	ListP = c(ListP, as.numeric(P) );

} #EndforIndMat

AvPrecision = AvPrecision/10;
ListP = c(ListP, AvPrecision);

GlobalTime <- proc.time() - TimeNow ;
print("time consuming per process"); print(GlobalTime/10);	# output time consuming

return( new("ClusterEval",  Precision=ListP ) );

}
setMethod("ClusterEval",signature(x="findModelCluster"), ClusterEval.crossval)


########################################################################
# Compute Final Evaluation 
#
########################################################################

# Usage:
#  P =  ClusterEval.final ( x )
#

ClusterEval.final <- function (  x  ) 
{

#parameters init
GlobalTime = 0;
AvPrecision = 0;
ListP = c();

TimeNow <- proc.time();		# catch time reference
	
Alert("", "sampling matrix...", "\t\t");
ListMatrixLearnTest = DataSampling(DatMat=x@Matrice$Mat);
cx = 0; cy = 0;
Alert("", "ok", "\n");

	MatriceMat  = as.matrix( ListMatrixLearnTest$ListMatLearn[[1]] );
	MatriceTest = as.matrix( ListMatrixLearnTest$ListMatTest[[1]]  );				

	nlin = nrow(MatriceTest);
	ncol = ncol(MatriceTest);
	pp = c();
	for(i in 1:nlin ){	# we fill the full matrix
		pp = c(pp, as.vector(MatriceTest[i,]) );
	}#finfori

	if( cx != 0 ) {
		MatriceMat = MatriceMat[,c(cx,cy,ncol(MatriceMat))];
	}
	else { 
		if( sign(min(MatriceMat)) < 0 ) {
			MatAdjCoa   = dudi.pca(as.data.frame(MatriceMat[,1:(ncol(MatriceMat)-1)]), scan = FALSE);
		}
		else
			MatAdjCoa   = dudi.coa(as.data.frame(MatriceMat[,1:(ncol(MatriceMat)-1)]), scan = FALSE);
		MatriceMat = as.data.frame( c( MatAdjCoa$li , as.data.frame(MatriceMat[,ncol(MatriceMat)]) ) );
	} #EndIf

	x = findModelCluster.Eval( as.matrix(MatriceMat) );

	NBClass		= max( x@ClassPoints );
	NbClassInData	= max( x@Matrice$Mat[,ncol] );

	ClassPoints =   .C("MatchGridPoint_C",     
			as.vector(pp) , 
			as.vector(x@MatriceK) , 
			as.integer(x@SizeGrid), 
			as.integer(ncol), 
			as.integer(nlin), 
			as.double(x@MinMaxXY[1]), 
			as.double(x@MinMaxXY[2]), 
			as.double(x@MinMaxXY[3]), 
			as.double(x@MinMaxXY[4]), 
			as.integer(NBClass), 
			as.integer(x@KNN), 
			as.integer(0), 
			as.integer(1), 
			as.vector(x@NumPoints) , 
			iClassPoints = numeric(nlin) )$iClassPoints;

	if( NbClassInData > 0 ) {
		Alert("\t\t", "evaluation...", "");

		MisClass = c();

		MisClass =   .C("Evaluation_C",     
			as.vector(pp) ,  
			as.integer(nlin), 
			as.integer(ncol), 
			as.integer(NBClass), 
			as.vector(ClassPoints ) , 
			iMisClass = numeric(nlin) )$iMisClass;

		Alert("\t\t", "ok", "\n");
		P= (100*(nlin-sum(MisClass[]>0))/nlin) ;
		cat("\t\t\t", "Precision ", P, "\n");					
	}
	Alert("\t\t\t\t" ,"ok", "\n");

GlobalTime <- proc.time() - TimeNow ;
print("time consuming per process"); print(GlobalTime/10);	# output time consuming

return( new("ClusterEval",  Precision=c(P) ) );

}
setMethod("ClusterEval",signature(x="findModelCluster"), ClusterEval.final)

########################################################################
# Test of scalability . running from 1% of data to 100%
#
########################################################################

# Usage:
#  ListP = ClusterEval.scalable (x );
#

ClusterEval.scalable <- function (  x  )
{ 

#parameters init
GlobalTime	= 0;
cx		= x@Cx;
cy		= x@Cy;
ListP		= c();
AvPrecision	= numeric(0);

Alert("", "loading matrix...", "\t\t");
MM   = x@Matrice$Mat;
nlin = nrow(MM);
Alert("", "ok", "\n");

TimeNow        <- proc.time();			# catch time reference

for(NLine in 2:nlin ) {
	
	Alert("\n", paste("test number...",NLine) , "\n\n");

	Alert("", "sampling matrix...", "\t\t");
	MatriceMat = MM[1:NLine,];
	Alert("\t\t", "ok", "\n");

	if( cx != 0 ) {
		MatriceMat = MatriceMat[,c(cx,cy,ncol(MatriceMat))];
	}
	else { 
		if( sign(min(MatriceMat)) < 0 ) {
			MatAdjCoa   = dudi.pca(as.data.frame(MatriceMat[,1:(ncol(MatriceMat)-1)]), scan = FALSE);
		}
		else
			MatAdjCoa   = dudi.coa(as.data.frame(MatriceMat[,1:(ncol(MatriceMat)-1)]), scan = FALSE);
		MatriceMat = as.data.frame( c( MatAdjCoa$li , as.data.frame(MatriceMat[,ncol(MatriceMat)]) ) );
	} #EndIf

	x_loop = findModelCluster.Eval( as.matrix(MatriceMat) );

	P = ClusterEval.Pvalue( x_loop );

	AvPrecision = AvPrecision + P@Precision;

	ListP = c(ListP, P@Precision );
	
	rm(P,x,MatAdjCoa);
	invisible(gc());					# freeing memory

} #Endfori

invisible(gc());					# freeing memory
ListP = c(ListP, AvPrecision);

GlobalTime <- proc.time() - TimeNow ;
print("time consuming per process"); print(GlobalTime/10);	# output time consuming

return( new("ClusterEval",  Precision=ListP ) );

}
setMethod("ClusterEval",signature(x="findModelCluster"), ClusterEval.scalable)

########################################################################
# Test clusterability 
#
########################################################################

# Usage:
#   ListP = ClusterEval.clusterable();
#

ClusterEval.clusterable <- function (  x  )
{ 

TimeNow <- proc.time();		# catch time reference

#parameters init
GlobalTime	= 0;
ListP		= c();
AvPrecision	= 0;

ListNu = c(0.1, 0.5, 1);
ListQ  = c(0.5, 1, 10, 100);

Lnu = length(ListNu);
Lq  = length(ListQ);

for(IndNu in 1:Lnu ) {
	
   Nu = ListNu[IndNu];
	
   for(IndQ in 1:Lq) {
	
	q = ListQ[IndQ];
	Alert("\n", paste("test number...",((IndNu-1)*Lq + IndQ) , "\t", "Nu ", Nu, "q ", q, "\n\n"));

	x_loop = findModelCluster( as.integer(1), MetLab=1, KernChoice=1, Nu, q, K=1, G=5, Cx=0, Cy=0, DName="iris", fileIn="D:\\rbuild\\test\\")

	P = ClusterEval.Pvalue( x_loop );

	AvPrecision = AvPrecision + P@Precision;

	ListP = c(ListP, P@Precision );
	
	invisible(gc());					# freeing memory

   } #EndforIndQ

} #EndforIndNu

GlobalTime <- proc.time() - TimeNow ;
print("time consuming per process"); print(GlobalTime/Lnu/Lq);	# output time consuming

return( new("ClusterEval",  Precision=ListP ) );

}
setMethod("ClusterEval",signature(x="missing"), ClusterEval.clusterable)

########################################################################
# Compute Evaluation 
# THIS FUNCTION IS OBSOLETE AND HAS BEEN REPLACED BY Evaluation_C
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
##End of Class##
########################################################################
