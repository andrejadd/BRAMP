##
##
## AA: extract values of specific segments defined by coordinates 'coords'
##
##


extractNodes <- function(nodeData, coords, xlocations, DEBUGINFOS=F) {

  ##  nodeData can be the predictor matrix or the target (response) vector.
  ##    if target:    a vector of target node values, each element in the vector corresponds to the target value of one location
  ##    if predictor: a matrix [(xlocs*ylocs) x nr.predictors] of predictor node values, each column corresponds to a node and each row to a location
  ##
  ##  coords : coordinates of a segment that is to be extracted, format = c(segid, x1,y1,x2,y2)
  ##  xlocations : the number of location along the X axis (of the grid)

  node.values = c()

  for (y in coords[2]:coords[4]) {

      startindex = (y - 1) * xlocations + coords[1]
      endindex = startindex + (coords[3] - coords[1])

      if(DEBUGINFOS) { cat("startindex ", startindex, ", endindex ", endindex, "\n") }

      if(is.matrix(nodeData)) {   ## in case its a matrix (predictors)
        node.values = rbind(node.values, Y[startindex:endindex,])
      } else {
        node.values = c(node.values, Y[startindex:endindex])
      }
    }

  return(node.values)
}






