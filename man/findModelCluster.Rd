\name{findModelCluster}
\docType{methods}
\alias{findModelCluster}
\alias{findModelCluster-methods}
\alias{findModelCluster,integer-method}
\alias{findModelCluster.chargeMatrix}
\alias{findModelCluster.Eval}
\alias{findModelCluster.Test}
\alias{findModelCluster.IsError}
\alias{findModelCluster,character-method}
\alias{findModelCluster,matrix-method}
\alias{findModelCluster,numeric-method}
\alias{findModelCluster,missing-method}
\alias{findModelCluster.chargeMatrix,character-method}
\alias{findModelCluster.Eval,matrix-method}
\alias{findModelCluster.IsError,numeric-method}
\alias{findModelCluster.Test,missing-method}

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
\S4method{findModelCluster}{integer}(x=as.integer(1), MetLab=1, KernChoice=1, Nu=0.8, q=20, K=1, G=10, Cx=1, Cy=2, DName="iris", fileIn="")

\S4method{findModelCluster.chargeMatrix}{character}( x="term", fileIn="c:/files")

\S4method{findModelCluster.Eval}{matrix}( x=matrix() )

\S4method{findModelCluster.Test}{numeric}()
}

\arguments{

  \item{x}{ means MetOpt parameter in standard use 
	    means DName in chargeMatrix use
	    means DatMat in Eval use, a Matrix given as unic argument }

  \item{MetLab}{ option taking value 1 (grid labelling) or 2 (mst labelling) or 3 (knn labelling) }
 
  \item{KernChoice}{ option taking value 0 (Euclidian) or 1 (RBF) or 2 (Exponential) }
  
  \item{Nu}{ kernel parameter  }
  
  \item{q}{ kernel parameter }
  
  \item{K}{ number of neigbours on the grid }
  
  \item{G}{ size of the grid }
  
  \item{Cx}{ 1st data coordinate to plot for 2D cluster extraction }
 
  \item{Cy}{ 2nd data coordinate to plot for 2D cluster extraction }
  
  \item{DName}{ Name of data which is the prefix of files :
   \file{DName\_mat.txt}, 
   \file{DName\_att.txt}, 
   \file{DName\_var.txt} }
  
  \item{fileIn}{ path where to find files as "D:\\R\\library\\svcR\\" }
 
}

\details{
  format of \file{DName\_mat.txt} (data matrix): 
  1 1 5.1	
  1 2 3.5
  2 3 1.4
it mean mat[1, 1] = 5.1, mat[1, 2] = 3.5, mat[2, 3] = 1.4
  
  format of \file{DName\_att.txt} : 
X1
X2
it mean X1 is the name of first column of the data matrix, X2 is the name of the second column of the data matrix

  format of \file{DName\_var.txt} : 
v1
v2
it mean v1 is the name of first line of the data matrix, v2 is the name of the second line of the data matrix
}

\value{
  An S4 object of class \code{findModelCluster} 
  The object is the svc model along with 
  the slots :
  \item{WVectorsYA}{lagrange coefficients : WVectorsYA\$A}
  \item{Matrice}{variables names Matrice\$var, attributes names Matrice\$Att and data Matrice\$Mat}
  \item{MatriceK}{kernel matrix}
  \item{Data}{Data Matrix}
  \item{MinMaxXY}{min max values for first and second coordinates}
  \item{MisClass}{missclassfied points}
  \item{DName}{prefix name of data for files decoding}
  \item{fileIn}{data file path}
  \item{ClassPoints}{class of grid points}
  \item{Cx}{x column id of data matrix}
  \item{Cy}{y column id of data matrix}
  \item{Nu}{nu value of the svc model}
  \item{KNN}{knn value for labelling}
  \item{SizeGrid}{size grid for labelling}
  \item{AroundNullVA}{almost null value for lagrange coefficient estimation}
  \item{NumPoints}{value fo grid points}

  slots can be accessed by \code{object@slot}
}

\references{ N.Turenne , Some Heuristics to speed-up Support Vector Clustering , technical report 2006, INRA, France 
\url{http://migale.jouy.inra.fr/~turenne/svc.pdf} 
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

MetOpt  = 1;    # optimisation method with randomization
MetLab  = 1;    # grid labelling
KChoice = 1;    # 0: eucli 1: radial 2: radial+dist 
Nu      = 1.0; 
q       = 2000;   # lot of clusters
K       = 1;    # only 1  nearest neighbour for clustering
Cx = Cy = 0; # we use principal component analysis factors
G       = 13; # size of the grid for cluster labelling
DName   = "iris";
fileIn  = ""; # fileIn migth be such as "D:/R/library/svc/", if NULL it will work on iris data

fmc = findModelCluster( as.integer(MetOpt), MetLab, KChoice, Nu, q, K, G, Cx, Cy, DName, fileIn); 

}

\keyword{cluster}
