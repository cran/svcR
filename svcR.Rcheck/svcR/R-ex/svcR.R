### Name: findModelCluster
### Title: Computation of clustering model by support vector machine
### Aliases: findModelCluster
### Keywords: cluster

### ** Examples


## exemple with iris data

MetOpt  = 1;    # optimisation method with randomization
MetLab  = 1;    # grid labelling
Nu      = 0.5; 
q       = 40;   # lot of clusters
K       = 1;    # only 1  nearest neighbour for clustering
Cx = Cy = 0; # we use principal component analysis factors
G       = 15; # size of the grid for cluster labelling
DName   = "iris";
fileIn  = ""; # fileIn migth be such as "D:/R/library/svc/", if NULL it will work on iris data

findModelCluster(MetOpt, MetLab, Nu, q, K, G, Cx, Cy, DName, fileIn); 




