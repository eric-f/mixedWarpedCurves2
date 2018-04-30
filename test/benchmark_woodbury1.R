library(microbenchmark)
ni <- 100
x <- seq(0, 1, length=ni)
f <- -4*(1-x)*x
F.tilde = cbind(1, f)
y <- rnorm(ni)
s2 = 5
s2_sh = 10
s2_sc = 0.1
Sigma = diag(c(s2_sh, s2_sc))
d = s2^2*s2_sh*s2_sc/((s2+ni*s2_sh)*(sum(f^2)*s2_sc+s2)-s2_sh*s2_sc*sum(f)^2)
K = d * matrix(c(sum(f^2)/s2+1/s2_sc, -sum(f)/s2, -sum(f)/s2, 1/s2_sh + ni/s2), 2, 2)
microbenchmark(S_post = solve(diag(s2, ni) + F.tilde %*% Sigma %*% t(F.tilde)), unit = "us")
microbenchmark(S_post2 = diag(ni)/s2 - F.tilde %*% K %*% t(F.tilde) / s2^2, unit="us")
max(abs(S_post-S_post2))

microbenchmark(y %*% solve(diag(s2, ni) + F.tilde %*% Sigma %*% t(F.tilde), y), unit="us")
microbenchmark(y %*% solve(diag(s2, ni) + F.tilde %*% Sigma %*% t(F.tilde)) %*% y, unit="us")
microbenchmark(
  y %*% (diag(ni)/s2 - F.tilde %*% K %*% t(F.tilde) / s2^2) %*% y
  , unit="us")

y %*% (diag(ni)/s2 - F.tilde %*% K %*% t(F.tilde) / s2^2) %*% y
y %*% diag(ni) %*% y /s2 - y %*% F.tilde %*% K %*% t(F.tilde) %*% y / s2^2

{
  term1 = sum(y^2)
  d = s2^2*s2_sh*s2_sc/((s2+ni*s2_sh)*(sum(f^2)*s2_sc+s2)-s2_sh*s2_sc*sum(f)^2)
  K = d * matrix(c(sum(f^2)/s2+1/s2_sc, -sum(f)/s2, -sum(f)/s2, 1/s2_sh + ni/s2), 2, 2)
  Fy = y %*% F.tilde
  term1/s2 - Fy %*% K %*% t(Fy) / s2^2
}
