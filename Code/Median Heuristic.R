### Compute beta of Laplacian kernel using median heuristic

MedianHeuristic <- function(X, Y){
  if(is.vector(X)){
    Z <- c(X,Y)
  }else{
    Z <- rbind(X,Y)
  }
  
  
  return(1/median(matrix(dist(Z, method = 'manhattan'))))
}