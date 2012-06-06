



newNode <- function(id=NULL, value=NULL, budget=NULL, parent=NULL){

  if(is.null(id)) id = 1
  
  object=new.env(parent=globalenv())  

  ## member parameters
  object$id = id
  object$value=value
  object$budget=budget

  ## these should be also of type 'NewNode()'
  object$parent=parent ## only the root node has a NULL here, all other can look up their parent node
  object$child1=NULL
  object$child2=NULL
  
  class(object)='pointer' 

  return(object)  
} 



addLeafpair <- function(tree.pointer, parent, child1, child2, budget=-1) {

  rec.tree.nodes = list(tree.pointer)

  ## find the leaf with parent id 
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      ## find the leaf with id 'parent' , so the parent that is to be split
      if(rec.node$id == parent && is.null(rec.node$child1) && is.null(rec.node$child2) ) {

        rec.node$child1 = newNode(child1, getFreeValue(tree.pointer), budget=budget, parent=rec.node)
        rec.node$child2 = newNode(child2, getFreeValue(tree.pointer), budget=budget, parent=rec.node)

        rec.tree.nodes.tmp = list() ## this leads to exit of while loop
        break ## exit of for
             
      } else {
        ## this is not the node, proceed with next ones
        
        if(!is.null(rec.node$child1)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
        }
        
        if(!is.null(rec.node$child2)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
        }
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }
}

getNrLeafPairs <- function(tree.pointer) { ## tree.ok

  pair.counter = 0
  
  rec.tree.nodes = list(tree.pointer)

  ## find the leaf with parent id 
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      found.leaf.parent = F
      
      ## check if this is a inner node with two childs
      if((!is.null(rec.node$child1)) && (!is.null(rec.node$child2)) ) {

        ## check if these childs are both leafs
        if(is.null(rec.node$child1$child1) && is.null(rec.node$child1$child2) && is.null(rec.node$child2$child1) && is.null(rec.node$child2$child2)) {
          pair.counter = pair.counter + 1
          found.leaf.parent = T
        }

      }

      if(!found.leaf.parent) {
        if(!is.null(rec.node$child1)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
        }
        
        if(!is.null(rec.node$child2)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
        }
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(pair.counter)

}

##
## I use this to identify existing pairs, later I select from it by chance and extract it with
## another method
##
getLeafPairIDs <- function(tree.pointer) { ## tree.ok

  pairs = matrix(0, nrow=0, ncol=2)
  pair.counter = 0
  
  rec.tree.nodes = list(tree.pointer)

  ## find the leaf with parent id 
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      found.leaf.parent = F
      
      ## check if this is a inner node with two childs
      if((!is.null(rec.node$child1)) && (!is.null(rec.node$child2)) ) {

        ## check if these childs are both leafs
        if(is.null(rec.node$child1$child1) && is.null(rec.node$child1$child2) && is.null(rec.node$child2$child1) && is.null(rec.node$child2$child2)) {

          ## add leaf pair to pairs matrix
          pairs = rbind(pairs, c(rec.node$child1$id, rec.node$child2$id))
          found.leaf.parent = T
        }

      }

      if(!found.leaf.parent) {
        if(!is.null(rec.node$child1)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
        }
        
        if(!is.null(rec.node$child2)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
        }
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(pairs)

}


getLeafNode <- function(tree.pointer, id) {


  rec.tree.nodes = list(tree.pointer)

  ## find the leaf with parent id 
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      ## check if this is the leaf 
      if(rec.node$id == id && is.null(rec.node$child1) && is.null(rec.node$child2) ) {

        return(rec.node)
             
      } else {
        ## this is not the node, proceed with next ones
        
        if(!is.null(rec.node$child1)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
        }
        
        if(!is.null(rec.node$child2)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
        }
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(NULL)
}



addNode <- function(tree.pointer, parentid, value) {

  rec.tree.nodes = list(tree.pointer)

  ## get a new id
  new.id = getFreeID(tree.pointer)

  
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      cat("id:", rec.node$id, "\n")

      ## check if this is the node
      if(rec.node$id == parentid) {
        
        if(is.null(rec.node$child1)) {
          rec.node$child1 = newNode(new.id, value, parent=rec.node)
          break
        }

        if(is.null(rec.node$child2)) {
          rec.node$child2 = newNode(new.id, value, parent=rec.node)
          break
        }
        
      } else {
        ## this is not the node, proceed with next ones
        
        if(!is.null(rec.node$child1)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
        }
        
        if(!is.null(rec.node$child2)) {
          rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
        }
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(new.id)
}


getFreeID <- function(tree.pointer) { ## return a id that has not been used by any node

  rec.tree.nodes = list(tree.pointer)

  ## collect all node ids
  id.vec = c()
  
  while(length(rec.tree.nodes) > 0) {

    ## add encountered nodes into here 
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      id.vec = c(id.vec, rec.node$id)

      if(!is.null(rec.node$child1)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
      }

      if(!is.null(rec.node$child2)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
      }
    }
    
    rec.tree.nodes = rec.tree.nodes.tmp
  }

#  print(id.vec)

  ## get the lowest not used id
  new.id = min( setdiff( c(min(id.vec):(max(id.vec) + 1)), id.vec))
  
  return(new.id)
}

getFreeValue <- function(tree.pointer) { ## return a id that has not been used by any node

  rec.tree.nodes = list(tree.pointer)

  ## collect all node ids
  id.vec = c()
  
  while(length(rec.tree.nodes) > 0) {

    ## add encountered nodes into here 
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      id.vec = c(id.vec, rec.node$value)

      if(!is.null(rec.node$child1)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
      }

      if(!is.null(rec.node$child2)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
      }
    }
    
    rec.tree.nodes = rec.tree.nodes.tmp
  }

#  print(id.vec)

  ## get the lowest not used id
  new.id = min( setdiff( c(min(id.vec):(max(id.vec) + 1)), id.vec))
  
  return(new.id)
}


isLeaf <- function(tree.pointer, id) {

  
  rec.tree.nodes = list(tree.pointer)
  
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      ## check if this is the node
      if(rec.node$id == id) {

        if(is.null(rec.node$child1) && is.null(rec.node$child2)) {
          return(TRUE)
        }
      
        return(FALSE)
      }

      if(!is.null(rec.node$child1)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
      }
      
      if(!is.null(rec.node$child2)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(FALSE)
}




drawTree <- function(tree.pointer, BY.ID=T) {

  outfile = "tmp.tree.dot"

  outstring = "digraph G{\n"

  ## here go all the node infos, start with root of tree
  labelstr = paste("  ", tree.pointer$value, "[ label=\"", tree.pointer$id, " (", tree.pointer$budget, ")\"] \n", sep="")
 
  rec.tree.nodes = list(tree.pointer)
  
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      ## check if child 1 exists
      if(!is.null(rec.node$child1)) {

        outstring = paste(outstring, " ", rec.node$value, " -> ", rec.node$child1$value, "\n", sep="")

        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1

        if(BY.ID) {
          labelstr = paste(labelstr, "  ", rec.node$child1$value, "[ label=\"", rec.node$child1$id, " (", sprintf("%.4f",rec.node$child1$budget), ")\"] \n", sep="")
        } else {
          labelstr = paste(labelstr, "  ", rec.node$child1$value, "[ label=\"", rec.node$child1$value, " (", rec.node$child1$id, ")\"] \n", sep="")
        }
      }

      ## check if child 2 exists
      if(!is.null(rec.node$child2)) {

        outstring = paste(outstring, " ", rec.node$value, " -> ", rec.node$child2$value, "\n", sep="")

        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2

        if(BY.ID) {
          labelstr = paste(labelstr, "  ", rec.node$child2$value, "[ label=\"", rec.node$child2$id, " (", sprintf("%.4f",rec.node$child2$budget), ")\"] \n", sep="")
        } else {
          labelstr = paste(labelstr, "  ", rec.node$child2$value, "[ label=\"", rec.node$child2$value, " (", rec.node$child2$id, ")\"] \n", sep="")
        }

      }
    }


    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  outstring = c(outstring, labelstr)

  outstring = c(outstring, "}\n")
  
  write(outstring, file = outfile) # ,  ncolumns = if(is.character(x)) 1 else 5, append = FALSE, sep = " ")

  
  system('dot -Tpng tmp.tree.dot -o tmp.tree.png')
}

getLeafIDs <- function(tree.pointer) {

  rec.tree.nodes = list(tree.pointer)

  ## collect the leaf ids
  leaf.ids = c()
  
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      has.childs = F
      
      if(!is.null(rec.node$child1)) {
        has.childs = T
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
      }
      
      if(!is.null(rec.node$child2)) {
        has.childs = T
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
      }

      if(!has.childs) {
        leaf.ids = c(leaf.ids, rec.node$id)
      }
    }


    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(leaf.ids)
}

testTree <- function() {

  tree = newNode(1,1,0)

  
  addLeafpair(tree, 1, 1 , 2)
  
#  addNode(tree, 1, value=runif(1))

#  for( g in 1: 20) {
    
#    leafs = getLeafIDs(tree)

#    append.to = sample(c( leafs, leafs), 1)

#    addNode(tree, append.to, value=runif(1))
#  }
  
  # cat("isLeaf id ", 2, ": ", isLeaf(tree, 2), "\n")

  drawTree(tree)
}


SegmentGrid2Tree <- function(grid.obj) {

  start.budget = 1 ## FIXME: get!
  
  tree = newNode(1,value=1)

  ## check if only root node exists
  if( nrow(grid.obj$inside.node.pairs) == 0 && nrow(grid.obj$leaf.pairs) == 0) {
    return(tree)
  }

  ## check the inside nodes
  if(nrow(grid.obj$inside.node.pairs) > 0) {  ## if exist

    for(r in 1:nrow(grid.obj$inside.node.pairs)) {

      pair = grid.obj$inside.node.pairs[r,]

      parent = min(pair)
      child1 = parent
      child2 = max(pair)

      addLeafpair(tree, parent, child1, child2)
    }
  }


  ## check the inside nodes
  if(nrow(grid.obj$leaf.pairs) > 0) {  ## if exist

    for(r in 1:nrow(grid.obj$leaf.pairs)) {

      pair = grid.obj$leaf.pairs[r,]

      parent = min(pair)
      child1 = parent
      child2 = max(pair)

      addLeafpair(tree, parent, child1, child2)
    }
  }

  return(tree)
}



delLeaf <- function(tree.pointer, id) {


  rec.tree.nodes = list(tree.pointer)
  
  while(length(rec.tree.nodes) > 0) {

    ## add new nodes here
    rec.tree.nodes.tmp = list()
    
    for(rec.node in rec.tree.nodes) {

      ## check child 1 if its a leaf
      if(!is.null(rec.node$child1)) {

        if(rec.node$child1$id == id && is.null(rec.node$child1$child1) && is.null(rec.node$child1$child2)) {
          rec.node$child1 = NULL
          return(T)
        }
      }

      ## check child 1 if its a leaf
      if(!is.null(rec.node$child2)) {
      
        if(rec.node$child2$id == id && is.null(rec.node$child2$child1) && is.null(rec.node$child2$child2)) {
          rec.node$child2 = NULL
          return(T)
        }
      }

      
       
      if(!is.null(rec.node$child1)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child1
      }
      
      if(!is.null(rec.node$child2)) {
        rec.tree.nodes.tmp[[length(rec.tree.nodes.tmp) + 1 ]] = rec.node$child2
      }
    }
  
    ## copy the potential childs that we want to recurse
    rec.tree.nodes = rec.tree.nodes.tmp
  }

  return(FALSE)


}
