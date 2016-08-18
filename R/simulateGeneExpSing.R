#' \code{simulateGeneExpSing} Simulate Gene Expression Data (Per Cell)
#'
#' @param model A CellModel
#' @param pathways A list of pathways, Format:(GtoM, GtoS, Prox)
#' @param perError User defined error for noise calculations
#' @param opt Option for which noise error calculated
#' @param sampFreq The frequency at which data is collected
#' @return the size of the cell population over time

setGeneric("simulateGeneExpSing", function(model,pathways,perError = 0.1,opt = 1, sampFreq = 1)
    standardGeneric("simulateGeneExpSing"))

setMethod("simulateGeneExpSing", "CellModel",
          
        function(model,pathways,perError = 0.1,opt = 1, sampFreq = 1) {
            
            #Get Individual Pathways
            gtompath = pathways[[1]]
            gtospath = pathways[[2]]
            proxpath = pathways[[3]]
            
            #Combined list of genes and exponential
            name = c(gtompath,gtospath,proxpath)
            name = name[!duplicated(name)]
            genes = rexp(length(name),1/3)
            genes = setNames(genes,name)
            
            #Set Gene Exponentials per Path
            temp = intersect(gtompath,name)
            nummgenes = genes[temp]
            temp = intersect(gtospath,name)
            numsgenes = genes[temp]
            temp = intersect(proxpath,name)
            proxgenes = genes[temp]
            
            #Create simulation data for each pathway
            gtom = simulateGToMPathSing(model,nummgenes, sampFreq)
            gtos = simulateGToSPathSing(model,numsgenes, sampFreq)
            prox = simulateProxPathSing(model,proxgenes, sampFreq)
            
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
                    bimatrix = matrix(0,length(name),length(gtom[[t]][1,]))
                    for(j in 1:length(gtom[[t]][1,])){
                        mu = t(t(rowMeans(t1)))
                        bi = NBsim(mu,t1)
                        bimatrix[,j] = bi
                    }
                    output[[t]] = t1 + bimatrix
                }
            }
            return(output)
        }
)
