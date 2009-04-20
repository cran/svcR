\name{kernelMatrix-class}
\docType{class}
\alias{kernelMatrix-class}

\title{ Class "kernelMatrix" }

\description{ reading data, computation of the svc model and cluster labelling }

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("kernelMatrix", ...)}.
  or by calling the  \code{kernelMatrix} function.}

\section{Slots}{
  \describe{

    \item{\code{matrixKernel}:}{Object of class \code{"vector"} contains
      kernel matrix with a vector format }

    \item{\code{matrixK}:}{Object of class \code{"matrix"} contains
       kernel matrix with a matrix format}
  }
}

%\section{Methods}{
%  \describe{ }
%}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\seealso{ \code{\link{kernelMatrix}} }

\examples{

## exemple with iris data

fmc = findModelCluster.Test(); 

fmc@MatriceK; # kernel matrix

}

\keyword{ classes }
