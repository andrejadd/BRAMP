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

  node.values = c()  ## becomes matrix with a rbind() below (in the case nodeData is a matrix)

  for (y in coords[2]:coords[4]) {

      startindex = (y - 1) * xlocations + coords[1]
      endindex = startindex + (coords[3] - coords[1])

      if(DEBUGINFOS) { cat("startindex ", startindex, ", endindex ", endindex, "\n") }

      if(is.matrix(nodeData)) {   ## in case its a matrix (predictors)

        ## get the rows of interest
        rows.tmp = nodeData[startindex:endindex,]

        ## extract rows without any NaNs (these are rows were no measurement was made, hence, should be ignored)
        ## then append to return matrix
        node.values = rbind(node.values, rows.tmp[apply(rows.tmp, 1, function(x) all(!is.nan(x))),] )
      } else {

        ## extract target elements
        elem.tmp = nodeData[startindex:endindex]

        ## exclude all NaNs (unsampled nodes)
        node.values = c(node.values, elem.tmp[which(!is.nan(elem.tmp))])
      }
    }

  return(node.values)
}






