\name{findModelCluster}
\alias{findModelCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
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
findModelCluster(MetOpt="", MetLab="", KernChoice="", Nu="", q="", K="", G="", Cx="", Cy="", DName="", fileIn="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MetOpt}{ option taking value 1 (randomization) or 2 (quadratic programming) }
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
  \item{fileIn}{ path where to find files }
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
  no return
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
Nu      = 0.7; 
q       = 1200;   # lot of clusters
K       = 1;    # only 1  nearest neighbour for clustering
Cx = Cy = 0; # we use principal component analysis factors
G       = 13; # size of the grid for cluster labelling
DName   = "iris";
fileIn  = ""; # fileIn migth be such as "D:/R/library/svc/", if NULL it will work on iris data

findModelCluster(MetOpt, MetLab, KChoice, Nu, q, K, G, Cx, Cy, DName, fileIn); 

}
\keyword{ cluster }% at least one, from doc/KEYWORDS
