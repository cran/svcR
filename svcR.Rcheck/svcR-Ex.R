### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign(".CheckExEnv", as.environment(2), pos = length(search())) # base
## add some hooks to label plot pages for base and grid graphics
setHook("plot.new", ".newplot.hook")
setHook("persp", ".newplot.hook")
setHook("grid.newpage", ".gridplot.hook")

assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   delayedAssign("T", stop("T used instead of TRUE"),
		  assign.env = .CheckExEnv)
	   delayedAssign("F", stop("F used instead of FALSE"),
		  assign.env = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       env = .CheckExEnv)
assign("..nameEx", "__{must remake R-ex/*.R}__", env = .CheckExEnv) # for now
assign("ptime", proc.time(), env = .CheckExEnv)
grDevices::postscript("svcR-Ex.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), env = .CheckExEnv)
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"), pager="console")
options(warn = 1)    
library('svcR')

assign(".oldSearch", search(), env = .CheckExEnv)
assign(".oldNS", loadedNamespaces(), env = .CheckExEnv)
cleanEx(); ..nameEx <- "svcR"

### * svcR

flush(stderr()); flush(stdout())

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




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", env = .CheckExEnv),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
