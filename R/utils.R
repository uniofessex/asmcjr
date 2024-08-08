

#'@description `gray.palette` generates 


gray.palette <- function(n, lower=.3, upper=.7){
  s <- seq(lower, upper, length=n)
  rgb(matrix(rep(s, each=3), ncol=3, byrow=TRUE))
}
