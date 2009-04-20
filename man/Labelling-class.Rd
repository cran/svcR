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
%    \item{Labelling.compute}{\code{signature(x="findModelCluster")}
%					 : computes the cluster labbelling}
%
%    \item{AdjacencyPP}{\code{signature(x = "findModelCluster")}
%					: computes a  pair adjacency}
%
%    \item{Adjacency}{\code{signature(x = "findModelCluster")}
%					: computes a matrix adjacency}
%
%    \item{MST_labelling}{\code{signature(x = "findModelCluster")}
%					: computes a Minimum Spanning Tree Labelling}
%
%    \item{KNN_labelling}{\code{signature(x = "findModelCluster")}
%					: computes a K-Nearest-Neighbours Tree Labelling}
%  }
%}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\seealso{ \code{\link{Labelling}} }

\examples{

## exemple with iris data

fmc = findModelCluster.Test(); 

fmc@NumPoints; # values of the grid

}

\keyword{ classes }
