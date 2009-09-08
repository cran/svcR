\name{plot}
\docType{methods}
\alias{plot}
\alias{ExportClusters}
\alias{Summary}
\alias{plot,findSvcModel-method}
\alias{ExportClusters,findSvcModel-method}
\alias{Summary,findSvcModel-method}
\alias{GetClustersTerm}
\alias{ShowClusters}
\alias{GetClusterID}
\alias{GetClustersTerm,findSvcModel-method}
\alias{ShowClusters,findSvcModel-method}
\alias{GetClusterID,findSvcModel-method}

\title{ Display function of an svc object  }

\description{
Plot permit to display a grid and cluster, ExportClusters can get cluster in text file. Summary gives
statistics about clusters.
}

\usage{
\S4method{plot}{findSvcModel}(x)

\S4method{ExportClusters}{findSvcModel}(fmc=new("findSvcModel"), NameFile="nf")

\S4method{Summary}{findSvcModel}(fmc=new("findSvcModel"))

\S4method{GetClustersTerm}{findSvcModel}(fmc=new("findSvcModel"), term="home")

\S4method{ShowClusters}{findSvcModel}(fmc=new("findSvcModel"))

\S4method{GetClusterID}{findSvcModel}(fmc=new("findSvcModel"), Id=1)
}

\arguments{
  \item{x}{ findSvcModel object }

  \item{fmc}{ findSvcModel object }

  \item{Id}{ an integer pointing to a cluster }

  \item{term}{ a string naming a variable }

  \item{NameFile}{ a string naming export files }
}

\author{ Nicolas Turenne - INRA France \email{nicolas.turenne@jouy.inra.fr} }

\examples{

## exemple with iris data

#generate a cluster model
fmc = findSvcModel.Test(); 

#plot the cluster
plot(fmc);

#export in text file on path defined with parameter catched  by findSvcModel function
ExportClusters(fmc, "Exportfilename");

#give information about basics statistics of clusters (size, attributes mean)
Summary(fmc);

#show all clusters
ShowClusters(fmc);

#show clusters where names of variables contain the string X2
GetClustersTerm(fmc, term="121");

#show a  clsuter with a specific id, if id is 0, the variables belongs to none cluster
GetClusterID(fmc, Id=1);

}

\keyword{aplot}
\keyword{list}