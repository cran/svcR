\name{ClusterEval}
\docType{methods}
\alias{ClusterEval}
\alias{ClusterEval-methods}
\alias{ClusterEval.scalable}
\alias{ClusterEval.Pvalue}
\alias{DataSampling}
\alias{ClusterEval.crossval}
\alias{ClusterEval.final}
\alias{ClusterEval.clusterable}
\alias{ClusterEval.scalable,findSvcModel-method}
\alias{ClusterEval.Pvalue,findSvcModel-method}
\alias{DataSampling,matrix-method}
\alias{ClusterEval.crossval,findSvcModel-method}
\alias{ClusterEval.final,findSvcModel-method}
\alias{ClusterEval.clusterable,findSvcModel-method}

\title{ Computation of evaluation }

\description{
Set of functions making evaluation. It requires a class column in data.
}

\usage{
\S4method{ClusterEval.Pvalue}{findSvcModel}(x=NULL)

\S4method{DataSampling}{matrix}(DatMat=matrix())

\S4method{ClusterEval.crossval}{findSvcModel}(x=NULL)

\S4method{ClusterEval.final}{findSvcModel}(x=NULL)

\S4method{ClusterEval.scalable}{findSvcModel}(x=NULL)

\S4method{ClusterEval.clusterable}{findSvcModel}(x=NULL)
}

\arguments{

  \item{x}{ a findSvcModel object }

  \item{DatMat}{ iris matrix obtain from a run with findSvcModel function }
  
}

\details{ DataSampling, ClusterEval.crossval, ClusterEval.final, ClusterEval.scalable and 
ClusterEval.clusterable are based only on iris data. ClusterEval.scalable make evaluation
building a matrix adding line by line to current matrix from iris matrix.ClusterEval.clusterable
makes evaluation with a set of q paramter values and a set of nu parameter values.
}

\value{
  return a numeric or a list of numeric meaning a precision value between 0 and 100 %.
  for the method ClusterEval.scalable and ClusterEval.clusterable it return a list of 
  precision values.
}

\references{ N.Turenne , Some Heuristics to speed-up Support Vector Clustering , technical report 2006, INRA, France 
\url{http://migale.jouy.inra.fr/~turenne/svc.pdf} 
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

# generation of the svc model
fmc = findSvcModel.Test(); 

# catch evaluation value
P = ClusterEval.Pvalue( fmc );

}

\keyword{arith}
\keyword{array}
