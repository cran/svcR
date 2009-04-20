\name{plot}
\docType{methods}
\alias{plot}
\alias{ExportClusters}
\alias{Summary}
\alias{plot,findModelCluster-method}
\alias{ExportClusters,findModelCluster-method}
\alias{Summary,findModelCluster-method}

\title{ Display function of an svc object  }

\description{
Plot permit to display a grid and cluster, ExportClusters can get cluster in text file. Summary gives
statistics about clusters.
}

\usage{
\S4method{plot}{findModelCluster}(x)

\S4method{ExportClusters}{findModelCluster}(fmc)

\S4method{Summary}{findModelCluster}(fmc)
}

\arguments{
  \item{x}{ findModelCluster object }

  \item{fmc}{ findModelCluster object }
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

#generate a cluster model
fmc = findModelCluster.Test(); 

#plot the cluster
plot(fmc);

#export in text file on path defined with parameter catched  by findModelCluster function
ExportClusters(fmc);

#give information about basics statistics of clusters (size, attributes mean)
Summary(fmc);

}

\keyword{aplot}
\keyword{list}