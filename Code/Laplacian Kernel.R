#### Compute the distance between x and y using laplacian kernel
LaplacianKernel <- function(x,y,beta){
  return(exp(-beta*sum(abs(x - y))))
}