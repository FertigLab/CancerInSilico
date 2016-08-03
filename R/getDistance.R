setGeneric("getDistance", function(model,time,cella,cellb)
    standardGeneric("getDistance"))

setMethod("getDistance", "CellModel",
          
            function(model,time,cella,cellb) {
                ida = vector()
                ida[1] = data@cells[[time]][cella]
                ida[2] = data@cells[[time]][cella+1]
                idb = vector()
                idb[1] = data@cells[[time]][cellb]
                idb[2] = data@cells[[time]][cellb+1]
                test = (ida[1] - idb[1])**2 + (ida[2] - idb[2])**2
                return(test)
            }
)