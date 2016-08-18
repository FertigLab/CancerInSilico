#' \code{simulateGeneExpGroup} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathways A list of pathways, Format:(GtoM, GtoS, Prox)
#' @param perError User defined error for noise calculations
#' @param opt Option for which noise error calculated
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGeneExpGroup", function(model,pathways,perError = 0.1,opt = 1,frequency =1 )
    standardGeneric("simulateGeneExpGroup"))

setMethod("simulateGeneExpGroup", "CellModel",
          
        function(model,pathways,perError = 0.1,opt = 1,frequency = 1) {
                  
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
            gtom = simulateGToMPathGroup(model,nummgenes)
            gtos = simulateGToSPathGroup(model,numsgenes)
            prox = simulateProxPathGroup(model,proxgenes)
            
            #Initialize some helper matrices
            output = matrix(0,model@parameters[2],length(name))
            rownames(output)<-1:model@parameters[2]
            colnames(output)<-name
            t2 = matrix(0,model@parameters[2],length(name))
            rownames(t2)<-1:model@parameters[2]
            colnames(t2)<-name
            t3 = matrix(0,model@parameters[2],length(name))
            rownames(t3)<-1:model@parameters[2]
            colnames(t3)<-name
            
            #GTOSPATH Check
            temp = intersect(gtospath,name)
            output[,temp] = output[,temp] + gtos[,temp]
            #GTOMPATH Check
            temp = intersect(gtompath,name)
            t2[,temp] = t2[,temp] + gtom[,temp]
            #ProxPath Check
            temp = intersect(proxpath,name)
            t3[,temp] = t3[,temp] + prox[,temp]
            
            #Merge
            temp = which(output < t2)
            output[temp] = t2[temp]
            temp = which(output < t3)
            output[temp] = t3[temp]
            
            #Noise Calculations
            #Normal Distribution
            if(opt == 1){
                sdmatrix = matrix(0,model@parameters[2],length(name))
                for(t in 1:model@parameters[2]){
                    sd = pmax(perError * output[t,1:length(name)], perError)
                    sdmatrix[t,] = sd
                }
                distro = rnorm(model@parameters[2] * length(name),output,sdmatrix)
                output = output + distro
            }
            #Negative Binomial
            else if(opt == 2){
                bimatrix = matrix(0,model@parameters[2],length(name))
                for(t in 1:model@parameters[2]){
                    mu = t(output[t,])
                    bi = NBsim(mu,t(output))
                    bimatrix[t,] = t(bi)
                }
                output = output + bimatrix
            }
            
            return(t(output))
                
        }
)