path = "/home/tom/Biostats_Research/CellModelPackageRepo/CellBasedModel"

library(Rcpp)
compileAttributes(paste(path,"/Package/",sep=""))
install.packages(paste(path,"/Package/",sep=""),repos=NULL)


