\name{ModelSV-class}
\docType{class}
\alias{ModelSV-class}

\title{ Class "ModelSV" }

\description{ reading data, computation of the svc model and cluster labelling }

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("ModelSV", ...)}.
  or by calling the  \code{ModelSV} function.}

\section{Slots}{
  \describe{

    \item{\code{lagrangeCoeff}:}{Object of class \code{"list"} contains
      the lagrange coefficients : VectorWA\$A }

    \item{\code{RadiusC}:}{Object of class \code{"numeric"} contains
      radius of the hypersphere containing data}

    \item{\code{SmallR}:}{Object of class \code{"numeric"} contains
      residu of radius }
  }
}

%\section{Methods}{
%  \describe{ }
%}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\seealso{ \code{\link{ModelSV}} }

\examples{

## exemple with iris data

fmc = findSvcModel.Test(); 

L = getlagrangeCoeff(fmc); # catch lagrange coefficients
L$A;                       # lagrange coefficients computed in the model

}

\keyword{ classes }
