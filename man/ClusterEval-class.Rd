\name{ClusterEval-class}
\docType{class}
\alias{ClusterEval-class}

\title{ Class "ClusterEval" }

\description{ reading data, computation of the svc model and cluster labelling }

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("ClusterEval", ...)}.
  or by calling the  \code{ClusterEval} function.}

\section{Slots}{
  \describe{

    \item{\code{Precision}:}{Object of class \code{"numeric"} contains
      rate of homogenity }

  }
}

%\section{Methods}{
%  \describe{ }
%}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\seealso{ \code{\link{ClusterEval}} }

\examples{

## exemple with iris data

#model computation
fmc = findModelCluster.Test(); 

#precision assessment
P = ClusterEval.Pvalue(fmc);

}

\keyword{classes}
