



extractYTargets <- function(Y, coords, xlocations, DEBUGINFOS=F) {

  # Y : a vector of target node values, each element in the vector corresponds to the target value of one location
  # asegment : coordinates of a segment that is to be extracted, format = c(segid, x1,y1,x2,y2)
  # xlocations : the number of location along the X axis (of the grid)

  targets = c()

  for (y in coords[2]:coords[4]) {

      startindex = (y - 1) * xlocations + coords[1]
      endindex = startindex + (coords[3] - coords[1])

      if(DEBUGINFOS) { cat("startindex ", startindex, ", endindex ", endindex, "\n") }

      targets = c(targets, Y[startindex:endindex])
    }

  return(targets)
}


## AA , put this together with extractYTargets ? with check 
extractXPredictors <- function(X, coords, xlocations, DEBUGINFOS=F) {

  # X : a matrix [(xlocs*ylocs) x nr.predictors] of predictor node values, each column corresponds to a node and each row to a location
  # coords : coordinates of a segment that is to be extracted, format = c(x1,y1,x2,y2)
  # xlocations : the number of location along the X axis (of the grid), needed for the offset

  predictors = c()

  for (y in coords[2]:coords[4]) {

      startindex = (y - 1) * xlocations + coords[1]
      endindex = startindex + (coords[3] - coords[1]) 

      if(DEBUGINFOS) { cat("startindex ", startindex, ", endindex ", endindex, "\n") }
      
      predictors = rbind(predictors, X[startindex:endindex,])
    }

  return(predictors)
}


