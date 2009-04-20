## S4 object definitions and assigment/accessor functions for the slots.
##
## created  10.04.09 nicolas turenne
## updated   


setClass("findModelCluster", representation(	WVectorsYA	= "list",
						Matrice		= "list",
						MatriceK	= "vector",
						Data		= "matrix",
						MinMaxXY	= "vector",
						MisClass	= "vector",
						DName		= "character",
						fileIn		= "character",
						ClassPoints	= "vector",
						Cx		= "numeric",
						Cy		= "numeric",
						Nu		= "numeric",
						KNN		= "numeric",
						SizeGrid	= "numeric",
						AroundNullVA	= "numeric",
						NumPoints	= "vector")  )


setClass("kernelMatrix", representation(	matrixKernel	= "vector", 
						matrixK		= "matrix")  )

setClass("ModelSV", representation(	VectorWA	= "list", 
					RadiusC		= "numeric",
					SmallR		= "numeric") )

setClass("Labelling", representation(	ClassPoints	= "vector", 
					NumPoints	= "vector") )

setClass("ClusterEval", representation(	Precision	= "numeric") )


########################################################################
##End of Class##
########################################################################
