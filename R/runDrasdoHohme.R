runDrasdoHohme <- function(initialNum,
                         runTime,
                         density,
                         cycleLengthDist,
                         inheritGrowth,
                         outputIncrement,
                         recordIncrement,
                         randSeed,
                         drugEffect,
                         drugTime,
                         ...)
  
{
  
  nG <- list(...)$nG
  if (is.null(nG)) {nG = 24}
  epsilon <- list(...)$epsilon
  if (is.null(epsilon)) {epsilon = 10}
  
  if (density > 0.1) {
    
    message("density too high to seed efficiently\n")
    stop()
    
  }
  
  delta <- 0.2 ## must be less than 4 or calculations break
  
  #timeIncrement is the time between each timestep
  timeIncrement = delta / (4 * nG * (4 - sqrt(2)))
  if (timeIncrement > delta * (min(cycleLengthDist) - 1) / (8 * nG * (sqrt(2) - 1))) {
    timeIncrement = delta * (min(cycleLengthDist) - 1) / (8 * nG * (sqrt(2) - 1))
  }
  maxDeform <- 2 * timeIncrement * nG * (4 - sqrt(2))
  grRates <- 2 * (sqrt(2) - 1) * timeIncrement * nG / (cycleLengthDist - 1)
  mcSteps <- ceiling(runTime / timeIncrement)
  maxTranslation <- delta / 2
  maxRotate <- acos((16 + delta ^ 2 - 4 * delta) / 16)
  outputIncrement2 <- floor(outputIncrement / timeIncrement)
  recordIncrement2 <- floor(recordIncrement / timeIncrement)  
  if (recordIncrement == 0) {
    recordIncrement2 <- 1
  }
  
  for (i in 1:length(drugEffect)) {
    
      drugEffect[[i]][1] <- 2 * (sqrt(2) - 1) * timeIncrement * nG / (drugEffect[[i]][1] - 1)
  
  }
  
  output <- tryCatch({
    
    CellModel(initialNum, mcSteps, density, maxTranslation,
              maxDeform, maxRotate, epsilon, delta, outputIncrement2,
              randSeed, drugEffect, grRates, inheritGrowth, nG, timeIncrement,recordIncrement2,drugTime)
    
  }, error = function(cond) {
    
    message(cond, '\n')
    stop()
    
  })
  
  cellMat <- new("CellModel",cells = output,parameters = c(initialNum,runTime,density,inheritGrowth,outputIncrement,randSeed,epsilon,nG,timeIncrement,recordIncrement,drugTime), paramCycleLengthDist = cycleLengthDist, paramDrugEffect = drugEffect)

  return(cellMat)
  
}
