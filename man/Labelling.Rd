\name{Labelling}
\docType{methods}
\alias{Labelling}
\alias{Labelling-methods}
\alias{AdjacencyPP}
\alias{Adjacency}
\alias{MST_labelling}
\alias{KNN_labelling}
\alias{Labelling.compute,findModelCluster-method}
\alias{AdjacencyPP,findModelCluster-method}
\alias{Adjacency,findModelCluster-method}
\alias{MST_labelling,findModelCluster-method}
\alias{KNN_labelling,findModelCluster-method}

\title{ Computation of clustering model by support vector machine }

\description{
Labelling methods try to identify clusters in a grid. 
}

\usage{
\S4method{Labelling.compute}{findModelCluster}(x, 	
		MetLab		= 1, 	
		MatriceKernel	= MK, 	
		MatriceK	= M_K, 	
		pp		= vpp, 	
		Nu		= 1,         
		G		= 1,         
		q		= 1,         
		ncol		= 1,         
		nlin		= 1,         
		RadiusC		= 2,        
		r		= 2,         
		KernChoice	= 0.01, 	
		NbClassInData	= 0.01  )

\S4method{AdjacencyPP}{findModelCluster}(x, MatriceKernel=matrix(), Vec1=vector(), Vec2=vector())

\S4method{Adjacency}{findModelCluster}(x, MatriceKernel=matrix())

\S4method{MST_labelling}{findModelCluster}(x, MatriceKernel=matrix())

\S4method{KNN_labelling}{findModelCluster}(x, MatriceKernel=matrix())
}

\arguments{

  \item{x}{ a findModelCluster object }

  \item{MetLab}{ option taking value 1 (grid labelling) or 2 (mst labelling) or 3 (knn labelling) }

  \item{MatriceKernel}{ kernel matrix with vector format }

  \item{MatriceK}{ kernel matrix with vector format }
  
  \item{pp}{ option taking value 0 (Euclidian) or 1 (RBF) or 2 (Exponential) }
  
  \item{Nu}{ kernel parameter  }
  
  \item{G}{ kernel parameter }

  \item{q}{ kernel parameter }
  
  \item{ncol}{ number of attributes }
  
  \item{nlin}{ size of variables }
  
  \item{RadiusC}{ model  radius }
  
  \item{r}{ residual radius }
  
  \item{KernChoice}{ kernel id }
  
  \item{NbClassInData}{ number of classes in data }

  \item{Vec1}{ vector for computing adjacency }

  \item{Vec2}{ vector for computing adjacency }
}

%\details{ }

\value{
  An S4 object of class \code{labelling} 
  The object is the svc model along with 
  the slots :
  \item{ClassPoints}{class of grid points}
  \item{NumPoints}{value fo grid points}

  slots can be accessed by \code{object@slot}
}

\references{ N.Turenne , Some Heuristics to speed-up Support Vector Clustering , technical report 2006, INRA, France 
\url{http://migale.jouy.inra.fr/~turenne/svc.pdf} 
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

fmc = findModelCluster.Test(); 

fmc@NumPoints; # values of the grid

}

\keyword{arith}
\keyword{array}
\keyword{cluster}
