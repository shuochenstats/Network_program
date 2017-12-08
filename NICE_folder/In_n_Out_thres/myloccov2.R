myloccov2=function (X, X0, i0, f, ests, N) 
{
  d = ests[1]
  s = ests[2]
  p0 = ests[3]
  theo = I(ncol(X0) == 1)
  Xtil <- X[i0, ]
  X0til <- X0[i0, ]
  G <- t(X) %*% (f * X)
  G0 <- t(X0til) %*% X0til
  B0 <- X0 %*% (solve(G0) %*% t(X0til)) %*% Xtil
  C <- B0 - X
  Ilfdr = C %*% solve(G, t(X))
  Cov <- C %*% solve(G) %*% t(C)
  if (theo) 
    D = matrix(1, 1, 1)
  else D = matrix(c(1, 0, 0, d, s^2, 0, s^2 + d^2, 2 * d * 
                      s^2, s^3), 3)
  gam. = solve(G0, t(X0til)) %*% (Xtil %*% solve(G, t(X)))
  pds. = D %*% gam.
  if (theo) 
    pds. = rbind(pds., matrix(0, 2, nrow(X)))
  pds.[1, ] = pds.[1, ] - 1/N
  m1 = pds. %*% f
  m2 = pds.^2 %*% f
  stdev = as.vector(sqrt(m2 - m1^2/N))
  stdev[1] = p0 * stdev[1]
  pds.[1, ] = p0 * pds.[1, ]
  rownames(pds.) = c("p", "d", "s")
  list(Ilfdr = Ilfdr, pds. = pds., stdev = stdev, Cov = Cov)
}