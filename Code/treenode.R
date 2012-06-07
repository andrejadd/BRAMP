


createTree <- function() {

  my.tree = list(node=NULL)

  my.tree$node = createNode(1)

  my.tree$node$leftchild = createNode(2)
  my.tree$node$rightchild = createNode(2)

}

createNode <- function(budget) {

  node = list(budget=budget, leftchild=NULL, rightchild=NULL)
    
  return(node)
}
