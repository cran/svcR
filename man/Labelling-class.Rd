\name{Labelling-class}
\docType{class}
\alias{Labelling-class}
\alias{Labelling.compute}

\title{ Class "Labelling" }

\description{ reading data, computation of the svc model and cluster labelling }

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("Labelling", ...)}.
  or by calling the  \code{Labelling} function.}

\section{Slots}{
  \describe{

    \item{\code{ClassPoints}:}{Object of class \code{"vector"} contains
      class of grid points }

    \item{\code{NumPoints}:}{Object of class \code{"vector"} contains
      value fo grid points }
  }
}

%\section{Methods}{
%  \describe{
%
%    \item{Labelling.compute}{\code{signature(x="findSvcModel")}
%					 : computes the cluster labbelling}
%
%    \item{AdjacencyPP}{\code{signature(x = "findSvcModel")}
%					: computes a  pair adjacency}
%
%    \item{Adjacency}{\code{signature(x = "findSvcModel")}
%					: computes a matrix adjacency}
%
%    \item{MST_labelling}{\code{signature(x = "findSvcModel")}
%					: computes a Minimum Spanning Tree Labelling}
%
%    \item{KNN_labelling}{\code{signature(x = "findSvcModel")}
%					: computes a K-Nearest-Neighbours Tree Labelling}
%  }
%}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\seealso{ \code{\link{Labelling}} }

\examples{

## exemple with iris data

fmc = findSvcModel.Test(); 

getNumPoints(fmc); # values of the grid

}

\keyword{ classes }
