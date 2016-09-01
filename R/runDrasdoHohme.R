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
  max_incr = delta * (min(cycleLengthDist) - 1) / (8 * nG * (sqrt(2) - 1))

  if (timeIncrement > max_incr) {

    timeIncrement = max_incr

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
  
  output <- CellModel(initialNum, mcSteps, density, maxTranslation,
              maxDeform, maxRotate, epsilon, delta, outputIncrement2,
              randSeed, drugEffect, grRates, inheritGrowth, nG, timeIncrement,recordIncrement2,drugTime)
    
  cellMat <- new("CellModel",
                    m_cells = output,
                    m_initialNumCells = initialNum,
                    m_runTime = runTime,
                    m_initialDensity = density,
                    m_inheritGrowth = inheritGrowth,
                    m_outputIncrement = outputIncrement,
                    m_randSeed = randSeed,
                    m_epsilon = epsilon,
                    m_nG = nG,
                    m_timeIncrement = timeIncrement,
                    m_cycleLengthDist = cycleLengthDist)
  
  return(cellMat)
  
}
