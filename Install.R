#path = "E:/Repository/CellBasedModel"

#library(Rcpp)
#library(testthat)
# compileAttributes(paste(path,"/Package/",sep=""))
# install.packages(paste(path,"/Package/",sep=""),repos=NULL,type="binary")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
install.packages("CellModel_1.0.tar.gz",repos=NULL)
