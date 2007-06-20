### Name: findModelCluster
### Title: Computation of clustering model by support vector machine
### Aliases: findModelCluster
### Keywords: cluster

### ** Examples


## exemple with iris data

MetOpt  = 1;    # optimisation method with randomization
MetLab  = 1;    # grid labelling
KChoice = 1;    # 0: eucli 1: radial 2: radial+dist 
Nu      = 0.7; 
q       = 1200;   # lot of clusters
K       = 1;    # only 1  nearest neighbour for clustering
Cx = Cy = 0; # we use principal component analysis factors
G       = 13; # size of the grid for cluster labelling
DName   = "iris";
fileIn  = ""; # fileIn migth be such as "D:/R/library/svc/", if NULL it will work on iris data

findModelCluster(MetOpt, MetLab, KChoice, Nu, q, K, G, Cx, Cy, DName, fileIn); 




