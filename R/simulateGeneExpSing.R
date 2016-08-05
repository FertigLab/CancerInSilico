#' \code{simulateGeneExpSing} Simulate Gene Expression Data (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A list of pathways, Format:(GtoM, GtoS, Prox)
#' @param perError User defined error for noise calculations
#' @param opt Option for which noise error calculated
#' @param success number of successes in Negative Binomial error model
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGeneExpSing", function(model,pathway,perError = 0.1,opt = 1,success = 1)
    standardGeneric("simulateGeneExpSing"))

setMethod("simulateGeneExpSing", "CellModel",
          
        function(model,pathway,perError = 0.1,opt = 1,success = 1) {
            
            gtompath = pathway[[1]]
            gtospath = pathway[[2]]
            proxpath = pathway[[3]]
            
            #Create simulation data for each pathway
            gtom = simulateGToMPathSing(model,gtompath)
            gtos = simulateGToSPathSing(model,gtospath)
            prox = simulateProxPathSing(model,proxpath)
            
            name = c(gtompath,gtospath,proxpath)
            name = name[!duplicated(name)]
            
            #Output Variable
            output = list()
            
            for(t in 1:model@parameters[2]){
                #Helper Matrices
                t1 = matrix(0,length(name),length(gtom[[t]][1,]))
                colnames(t1)<-1:length(gtom[[t]][1,])
                rownames(t1)<-name
                t2 = matrix(0,length(name),length(gtom[[t]][1,]))
                colnames(t2)<-1:length(gtom[[t]][1,])
                rownames(t2)<-name
                t3 = matrix(0,length(name),length(gtom[[t]][1,]))
                colnames(t3)<-1:length(gtom[[t]][1,])
                rownames(t3)<-name
                                
                #GTOSPATH Check
                temp = intersect(gtospath,name)
                t1[temp,] = t1[temp,] + gtos[[t]][temp,]
                #GTOMPATH Check
                temp = intersect(gtompath,name)
                t2[temp,] = t2[temp,] + gtom[[t]][temp,]
                #ProxPath Check
                temp = intersect(proxpath,name)
                t3[temp,] = t3[temp,] + prox[[t]][temp,]
                
                #Merge Matrix
                temp = which(t1 < t2)
                t1[temp] = t2[temp]
                temp = which(t1 < t3)
                t1[temp] = t3[temp]
                #Noise Calculations
                #Normal Distribution
                if(opt == 1){
                    sdmatrix = matrix(0,length(name),length(gtom[[t]][1,]))
                    for(j in 1:length(gtom[[t]][1,])){
                        sd = pmax(perError * t1[,j], perError)
                        sdmatrix[,j] = sd
                    }
                    distro = rnorm(length(gtom[[t]][1,]) * length(name),t1,sdmatrix)
                    output[[t]] = t1 + distro
                }
                #Negative Binomial
                if(opt == 2){
                    phi = 1/success
                    bimatrix = matrix(0,length(name),length(gtom[[t]][1,]))
                    for(j in 1:length(name)){
                        bi = rnbinom(length(gtom[[t]][1,]),phi,mu = sum(t1[j,])/length(gtom[[t]][j,]))
                        bimatrix[j,] = bi
                    }
                    output[[t]] = t1 + bimatrix
                }
            }
            return(output)
        }
)