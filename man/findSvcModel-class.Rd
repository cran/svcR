\name{findSvcModel-class}
\docType{class}
\alias{findSvcModel-class}

\title{ Class "findSvcModel" }

\description{ reading data, computation of the svc model and cluster labelling }

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("findSvcModel", ...)}.
  or by calling the  \code{findSvcModel} function.}

\section{Slots}{
  \describe{

    \item{\code{lagrangeCoeff}:}{Object of class \code{"list"} contains
      the lagrange coefficients : lagrangeCoeff\$A }

    \item{\code{Matrice}:}{Object of class \code{"list"} contains
      variables names Matrice\$var, attributes names Matrice\$Att and data Matrice\$Mat}

    \item{\code{MatriceK}:}{Object of class \code{"vector"} contains
      kernel matrix}

    \item{\code{Data}:}{Object of class \code{"matrix"} contains
      Data Matrix }

    \item{\code{MinMaxXY}:}{Object of class \code{"vector"} contains
    min max values for first and second coordinates}

    \item{\code{MisClass}:}{Object of class \code{"vector"} contains
      missclassfied points}

    \item{\code{dataFrame}:}{Object of class \code{"list"} contains
      the prefix name of data for files decoding}

    \item{\code{fileIn}:}{Object of class \code{"character"} contains
      data file path}

    \item{\code{ClassPoints}:}{Object of class \code{"vector"} contains
      class of grid points}

    \item{\code{Cx}:}{Object of class \code{"numeric"} contains
      the x column id of data matrix}

    \item{\code{Cy}:}{Object of class \code{"numeric"} contains
      the y column id of data matrix}

    \item{\code{Nu}:}{Object of class \code{"numeric"} contains
      nu value of the svc model}

    \item{\code{KNN}:}{Object of class \code{"numeric"} contains
      knn value for labelling}

    \item{\code{SizeGrid}:}{Object of class \code{"numeric"} contains
      the size grid for labelling}

    \item{\code{AroundNullVA}:}{Object of class \code{"numeric"} contains
      almost null value for lagrange coefficient estimation}

    \item{\code{NumPoints}:}{Object of class \code{"vector"} contains
      value fo grid points}
  }
}

%\section{Methods}{
%  \describe{ }
%}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\seealso{ \code{\link{findSvcModel}} }

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

data(iris)
fmc = findSvcModel( iris, MetOpt, MetLab, KernChoice, Nu, q, K, G, Cx, Cy); 

NP = getNumPoints(fmc); # values of the grid

}

\keyword{ classes }
