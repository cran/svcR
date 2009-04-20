\name{ModelSV}
\docType{methods}
\alias{ModelSV}
\alias{ModelSV-methods}
\alias{ModelSV.compute}
\alias{OptimQuadProgWcluster}
\alias{ModelSV.compute,numeric-method}
\alias{OptimQuadProgWcluster,matrix-method}

\title{ Computation of lagrange coefficients and radius of the model }

\description{
Computes kernel matrix, lagrange coefficients, support vectors and radius
}

\usage{
\S4method{ModelSV.compute}{numeric}(x ,
	  MatriceKernel = NULL,
	  MatriceK      = NULL,
	  Nu		= 1,
          nlin		= 1,
          MaxIter	= 2,
          MaxValA	= 2,
          AroundNull	= 0.01,
	  AroundNullVA  = 0.01  )

\S4method{OptimQuadProgWcluster}{matrix}(MatriceKernel ,
	  Nu		= 1,
          MaxValA	= 2,
          MinW		= 0.0001 )

}

\arguments{

  \item{x}{ method of lagrange computation 1 or 2  }

  \item{MatriceKernel}{ kernel matrix with vector format }
  
  \item{MatriceK}{ kernel matrix with vector format }
  
  \item{Nu}{ svc parameter }
  
  \item{nlin}{ number   }of line for data matrix
  
  \item{MaxIter}{ maximum iteration for coefficients computation }
  
  \item{MaxValA}{ number of neigbours on the grid }
  
  \item{AroundNull}{ almost null parameter }
  
  \item{AroundNullVA}{ almost null parameter for coefficients }

  \item{MinW}{ min value for coefficients }

}

%\details{}

\value{
  An S4 object of class \code{ModelSV} 
  The object is the svc model along with 
  the slots :
  \item{VectorWA}{lagrange coefficients : VectorsYA\$A}
  \item{RadiusC}{radius of the hypersphere containing data}
  \item{SmallR}{residu of radius}
  OptimQuadProgWcluster return a list containing lagrange coefficients.

  slots can be accessed by \code{object@slot}
}

\references{ N.Turenne , Some Heuristics to speed-up Support Vector Clustering , technical report 2006, INRA, France 
\url{http://migale.jouy.inra.fr/~turenne/svc.pdf} 
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

# model computing 
fmc = findModelCluster.Test(); 

fmc@WVectorsYA$A; # lagrange coefficients

}

\keyword{optimize}
\keyword{classif}
