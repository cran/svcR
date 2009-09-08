\name{findSvcModel}
\docType{methods}
\alias{findSvcModel}
\alias{findSvcModel-methods}
\alias{findSvcModel,integer-method}
\alias{findSvcModel.chargeMatrix}
\alias{findSvcModel.chargeMatrixHDD}
\alias{findSvcModel.Eval}
\alias{findSvcModel.Test}
\alias{findSvcModel.IsError}
\alias{findSvcModel,character-method}
\alias{findSvcModel,matrix-method}
\alias{findSvcModel,list-method}
\alias{findSvcModel,numeric-method}
\alias{findSvcModel,missing-method}
\alias{findSvcModel.chargeMatrixHDD,character-method}
\alias{findSvcModel.chargeMatrix,list-method}
\alias{findSvcModel.Eval,matrix-method}
\alias{findSvcModel.IsError,numeric-method}
\alias{findSvcModel.Test,missing-method}
\alias{getlagrangeCoeff}
\alias{getMatriceK}
\alias{getClassPoints}
\alias{getNumPoints}
\alias{getlagrangeCoeff,findSvcModel-method}
\alias{getMatriceK,findSvcModel-method}
\alias{getClassPoints,findSvcModel-method}
\alias{getNumPoints,findSvcModel-method}

\title{ Computation of clustering model by support vector machine }

\description{
SvcR implements a clustering algorithm based on separator search in a feature 
space between points described in a data space. Data format is defined by 
an attribute/value table (matrix). The data are transformed within a kernel 
to a feature space into a unic cluster bounded with a ball radius and support vectors. 
We can used the radius of this ball in the data space 
to reconstruct the boundary shaped now in several clusters. 
}

\usage{
\S4method{findSvcModel}{list}( x=iris, MetOpt="optimStoch", MetLab="gridLabeling", KernChoice="KernGaussian", Nu=0.8, q=20, K=1, G=10, Cx=1, Cy=2 )

\S4method{findSvcModel.chargeMatrixHDD}{character}( x="term", fileIn="c:/files" )

\S4method{findSvcModel.chargeMatrix}{list}( x=iris )

\S4method{findSvcModel.Eval}{matrix}( x=matrix() )

\S4method{findSvcModel.Test}{numeric}()

\S4method{getNumPoints}{findSvcModel}( object=new("findSvcModel") )

\S4method{getClassPoints}{findSvcModel}( object=new("findSvcModel") )

\S4method{getMatriceK}{findSvcModel}( object=new("findSvcModel") )

\S4method{getlagrangeCoeff}{findSvcModel}( object=new("findSvcModel") )
}

\arguments{

  \item{x}{ means dataFrame parameter in standard use 
	    means dataFrame in chargeMatrix use
	    means DatMat in Eval use, a Matrix given as unic argument }

  \item{MetOpt}{ option taking value "optimStoch" (stochastic way of optimization) or "optimQuad" (quadratic way of optimization) }

  \item{MetLab}{ option taking value "gridLabeling" (grid labelling) or "mstLabeling" (mst labelling) or "knnLabeling" (knn labelling) }
 
  \item{KernChoice}{ option taking value "KernLinear" (Euclidian) or "KernGaussian" (RBF) or "KernGaussianDist" (Exponential) or "KernDist" (Matrix data as Kernel value) }

  \item{Nu}{ kernel parameter  }
  
  \item{q}{ kernel parameter }
  
  \item{K}{ number of neigbours on the grid }
  
  \item{G}{ size of the grid }
  
  \item{Cx}{ 1st data coordinate to plot for 2D cluster extraction }
 
  \item{Cy}{ 2nd data coordinate to plot for 2D cluster extraction }
    
  \item{fileIn}{ path where to find files as "D:\\R\\library\\svcR\\" }
 
  \item{object}{ a findSvcModel object }

}

\details{

The main function of the package is called findSvcModel. It takes a data set as first argument.
This dataset can be either (1) a data.frame() structure , or (2) text files stored on hard disk drive.

In case (1) data.frame is classical standard list format (see Iris data).

In case (2)  format of \file{dataName\_mat.txt} (data matrix): 
  1 1 5.1	
  1 2 3.5
  2 3 1.4
it mean mat[1, 1] = 5.1, mat[1, 2] = 3.5, mat[2, 3] = 1.4
  
  format of \file{dataName\_att.txt} : 
X1
X2
it mean X1 is the name of first column of the data matrix, X2 is the name of the second column of the data matrix

  format of \file{dataName\_var.txt} : 
v1
v2
it mean v1 is the name of first line of the data matrix, v2 is the name of the second line of the data matrix

For the labeling parameter (MetLab) three choice are available: "gridLabeling", "mstLabeling" and "knnLabeling".

For the kernel parameter (KernelChoice) four choices are available: "KernLinear", "KernGaussian",  "KernGaussianDist" and "KernDist".

}

\value{
  An S4 object of class \code{findSvcModel} 
  The object is the svc model along with 
  the slots :
  \item{lagrangeCoeff}{lagrange coefficients : getlagrangeCoeff\$A}
  \item{Matrice}{variables names Matrice\$var, attributes names Matrice\$Att and data Matrice\$Mat}
  \item{MatriceK}{kernel matrix}
  \item{Data}{Data Matrix}
  \item{MinMaxXY}{min max values for first and second coordinates}
  \item{MisClass}{missclassfied points}
  \item{dataFrame}{prefix name of data for files decoding}
  \item{fileIn}{data file path}
  \item{ClassPoints}{class of grid points}
  \item{Cx}{x column id of data matrix}
  \item{Cy}{y column id of data matrix}
  \item{Nu}{nu value of the svc model}
  \item{KNN}{knn value for labelling}
  \item{SizeGrid}{size grid for labelling}
  \item{AroundNullVA}{almost null value for lagrange coefficient estimation}
  \item{NumPoints}{value fo grid points}

  slots can be accessed by \code{getlagrangeCoeff(object)}, \code{getMatriceK(object)},
  \code{getClassPoints(object)}, \code{getNumPoints(object)}

}

\references{ N.Turenne , Some Heuristics to speed-up Support Vector Clustering , technical report 2006, INRA, France 
\url{http://migale.jouy.inra.fr/~turenne/svc.pdf} 
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

MetOpt     = "optimStoch";      #  optimisation method with randomization
MetLab     = "gridLabeling";    #  grid labelling
KernChoice = "KernGaussian";    #  radial kernel
Nu         = 1.0; 
q          = 2000;   # lot of clusters
K          = 1;      # only 1  nearest neighbour for clustering
Cx = Cy    = 0;      # we use principal component analysis factors
G          = 20;     # size of the grid for cluster labelling
fileIn     = "";     # fileIn migth be such as "D:/R/library/svc/", if NULL it will work on term data

# usage example with a data frame 
data(iris);
fmc = findSvcModel( iris, MetOpt, MetLab, KernChoice, Nu, q, K, G, Cx, Cy); 
plot(fmc);

# usage example with import from textfiles
dataTerm = findSvcModel.chargeMatrixHDD("term", fileIn="");
fmc = findSvcModel( dataTerm, MetOpt , MetLab , KernChoice , Nu, q, K, G, Cx, Cy)
plot(fmc);

#
}

\keyword{cluster}
