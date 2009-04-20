\name{kernelMatrix}
\docType{methods}
\alias{kernelMatrix}
\alias{kernelMatrix-methods}
\alias{kernelMatrix.compute}
\alias{kernelMatrix.compute,vector-method}

\title{ Computation of clustering model by support vector machine }

\description{
This function implements three kinds of kernel to build a kernal matrix:
Euclidian, Radial-Base and Exponential.
}

\usage{
\S4method{kernelMatrix.compute}{vector}(x, 
          SymMat	= 1,
          q		= 1,
          ncol		= 2,
          nlin		= 2,
          KernChoice	= 1)

}

\arguments{

  \item{x}{ data matrix with a vector format }
  
  \item{SymMat}{ indicator to know if matrix is symetrical or not }
  
  \item{q}{ q parameter of the model }
  
  \item{ncol}{ number of columns of the data matrix }

  \item{nlin}{ number of lines of the data matrix }

  \item{KernChoice}{ option taking value 0 (Euclidian) or 1 (RBF) or 2 (Exponential) }
}

%\details{}

\value{
  An S4 object of class \code{kernelMatrix} 
  The object is the svc model along with 
  the slots :
  \item{matrixKernel}{kernel matrix with a vector format}
  \item{matrixK}{kernel matrix with a matrix format}

  slots can be accessed by \code{object@slot}
}

\references{ N.Turenne , Some Heuristics to speed-up Support Vector Clustering , technical report 2006, INRA, France 
\url{http://migale.jouy.inra.fr/~turenne/svc.pdf} 
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

fmc = findModelCluster.Test(); 

fmc@MatriceK;

}

\keyword{arith}
