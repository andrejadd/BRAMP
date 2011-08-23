
##
##
## create 22.02.2011
##
## spatAutoCorrelation calculates the value of 'a' which is defined as
##
## a = (sum(inv_d_ij * x_i))/(sum(inv_d_ij)) , where w_i is the inverse euclidian distance to the neighbour cells which are checked over the sum
## 
##  x_i is the density of the target species at location i (the neighbour cells)
##  

spatAutoCorrelation <- function (Y,xlocs,ylocs,nb=4) {



     neighbours = list(c(-1,0), c(1,0), c(0,-1),c(0,1))	       

     ## this vector holds the node values of avg population at neighbouring cells
     vec_a = c()

     if(nb == 4) {


	   for(yj in 1:ylocs) {
	   	  for(xj in 1:xlocs) {

  	      	   	 suma = 0
		   	 sumb = 0


		  	 ## loop all neighbour cells
		  	 for(nbs in neighbours) {			  	 
			 
			      ## coordinates of neighbour cell				    
			      xi = xj + nbs[1]
			      yi = yj + nbs[2]
		
		#	      cat("\nxi: ", xi, ", yi: ", yi, "\t")

			      ## check if in grid, otherwise ignore
			      ## FIXME: take other (alternative) neighbour
			      if(xi > 0 && xi <= xlocs && yi > 0 && yi <= ylocs) {

			        ## calculate inverse euclidian distance
				inv_d_ij = 1 / (sqrt( (xj-xi)^2 + (yj - yi)^2 ) )

			        ## transform in absolute (1D) coordinates
			        abs_pos = (yi - 1) * xlocs + xi

		#	        cat("in grid! abs_pos is ", abs_pos)


			        ## add up nominator and denominator
			        suma = suma + (inv_d_ij * Y[abs_pos])
				sumb = sumb + inv_d_ij
			     }
			 }
			 
			 vec_a = c(vec_a, (suma/sumb))
		  
		     }
              }
     }

     return( vec_a)

}