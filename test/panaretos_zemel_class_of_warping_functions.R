warping <- function(x, k){
  if (k==0) {return(x)}
  else x - sin(pi*k*x) / abs(k) / pi
}

x0 <- seq(0, 1, length=100)
plot(warping(x0, 0)~x0, type="l")
lines(warping(x0, -1)~x0)
lines(warping(x0, 1)~x0)
lines(warping(x0, -2)~x0)
lines(warping(x0, 2)~x0)
lines(warping(x0, -3)~x0)
lines(warping(x0, 3)~x0)
lines(warping(x0, -4)~x0)
lines(warping(x0, 4)~x0)
lines(warping(x0, -5)~x0)
lines(warping(x0, 5)~x0)


