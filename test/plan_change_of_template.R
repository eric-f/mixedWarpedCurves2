layout(matrix(1:4, nrow=2))


f1 <- function(x){
  -4*x*(1-x)
}
f1_v2 <- function(x){
  4*(x-0.5)^2-1
}
h <- function(x){
  2*x*(1-x)*(x-0.5) + x
}
f1_star <- function(x){
  16*x^2*(x-1)^2*(x+0.5)*(x+0.5-2)
}
f1_star_v2 <- function(x){
  # -16*x^2*(x-1)^2*(1+(x-0.5))*(1-(x-0.5))
  -16*x^2*(x-1)^2  *  (1-(x-0.5)^2)

}

curve(f1)
curve(f1_v2, add=T, col="red", lty=2)
curve(h(x))
curve(h(x)-x)
curve(f1(h(x)))
curve(f1_star, add=T, col="red", lty=2)
curve(f1_star_v2, add=T, col="blue", lty=3)



f2 <- function(x){
  -16*x^2*(1-x)^2
}
f2_v2 <- function(x){
  -16*(x-0.5)^4 + 8*(x-0.5)^2 - 1
}
curve(f2)
curve(f_star2, add=T, col="red", lty=1)


f1_star2 <- function(x){
  16*(x-0.5)^4 - 1
}
h2 <- function(x){
  2*sign(x-0.5)*(x-0.5)^2+0.5
}
layout(matrix(1:2, nrow=2))
curve(f1_star2)
curve(f1(h2(x)), add=TRUE, col="red", lty=2)
curve(h2)






h2 <- function(x){
  2*sign(x-0.5)*(x-0.5)^2+0.5
}
h2_inv <- function(y){
  sign(y-0.5)*sqrt(abs(y-0.5)/2)+0.5
}
curve(h2)
curve(h2_inv, add=T, col="blue")






h2 <- function(x, a=2){
  (sign(x-0.5) * abs(2*x-1)^a + 1)/2
}
h2_inv <- function(y, a=2){
  # sign(y-0.5)*(0.5^a*abs(2*y-1))^(1/a) + 0.5
  (sign(y-0.5)*abs(2*y-1)^(1/a)+1)/2
}
curve(h2(x, a=1.8))
curve(h2_inv(x, a=1.8), add=T, col="blue")
f1_star3 <- function(y, a=2){
  abs(2*y-1)^(2*a) - 1
}

curve(f1_star3(x, a=1.2))
curve(f1(h2(x, 1.2)), add=TRUE, col="red", lty=2)
curve(f1)
curve(f1_star3(h2_inv(x, 1.2), 1.2), add=TRUE, col="red", lty=2)
