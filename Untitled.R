getDistance<-function(x,w){
  return(abs(sum(x*w)/norm(w,'2')))
}

getDistance(c(1,1),c(1,-2))

