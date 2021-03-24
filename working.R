RAR_sim = function(nt, theta0, good.out = T, nb = 1, maxN = 500, N = 1000, upper = 0.975, uppfut = 0.95, lower = .01,
                   burn = 10*nt, response.type, conjugate_prior = T, padhere = rep(1,nt), adapt = T,
                   compCon = F, MID = 0) {
  if (good.out == F) {
    if (response.type == 'continuous') theta0 = - theta0
    if (response.type == 'binary') theta0 = 1 - theta0
  }
  ng = nb
  j = 0
  x = array(0, dim = c(nt, 1))
  x0 = array(0, dim = c(nt, 1))
  y = NULL
  theta = array(rnorm(N*nt, 0, 10), dim = c(N, nt, 1))
  check = array(0, dim = c(N, nt))
  fcheck = array(0, dim = c(N, nt))
  if (conjugate_prior != T | response.type == 'continuous') post0 = cbind(rep(0, nt), rep(10, nt))
  if (conjugate_prior == T & response.type == 'binary') {
    post0 = cbind(rep(1, nt), rep(1, nt))
    p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
    theta = array(- log((1 - p_new) / p_new), dim = c(N, nt, 1))
  }
  psup = array(rep(1/nt, nt), dim = c(nt, 1))
  psup_out = array(rep(1/nt, nt), dim = c(nt, 1))
  prand = rep(1/nt, nt)
  fut.stop = NULL
  repeat {
    j = j + 1
    if (j*nb>maxN) nb = maxN - (j-1)*nb
    xb = rmultinom(nb, 1, prob = prand)
    xb0 = xb
    for (k in 1:nb) {
      xb0[which(xb0[,k] == 1), k] = rbinom(1, 1, padhere[which(xb0[,k] == 1)])
    }
    if (response.type == 'binary') {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rbinom(1, 1, prob = z))
    } else {
      yb = apply(t(xb0)%*%theta0, 1, function(z) rnorm(1, z))
    }
    x = abind(x, xb, along = 2)
    x0 = abind(x0, xb0, along = 2)
    y = c(y, yb)
    if (response.type == 'binary' & conjugate_prior == T) {
      post0 = post_beta(yb, xb, post0)
      p_new = apply(post0, 1, function(x) rbeta(N, x[1], x[2]))
      theta_new = - log((1 - p_new) / p_new)
    } else {
      ysd = ifelse(length(y) == 1, 1, sd(y))
      post0 = post_Gauss(yb, xb, post0, ysd)
      theta_new = apply(post0, 1, function(x) rnorm(N, x[1], x[2]))
    }
    theta = abind(theta, theta_new, along = 3)
    if (compCon == F) {
      check[,which(prand>0)] = t(apply(theta[,which(prand>0),j+1], 1, sup_check))
      check[,which(prand==0)] = 0
      psup_out = abind(psup_out, apply(check, 2, mean), along = 2)
      if (length(y) < burn) {
        psup = abind(psup, psup[,j], along = 2)
      } else psup = abind(psup, apply(check, 2, mean), along = 2)
      ll = NULL
      if (length(which(psup[,j+1] < lower))>0) {
        ll = which(psup[,j+1] < lower)
      }
      if (length(ll) != 0) {
        fut.stop = c(fut.stop, ll)
      }
      if (!is.null(fut.stop)) {
        psup[fut.stop,j+1] = 0
      }
      ntj = length(which(psup[,j+1]>0)) 
      prand = rep(0, nt)
      if (adapt == T & length(y)>=burn)  {
        prand[which(psup[,j+1]>0)] = ((ntj - 1)/ntj)*sqrt(psup[which(addarmlater<=j & psup[,j+1]>0),j+1])
      } else prand[which(psup[,j+1]>0)] = 1/ntj
      
    } else {
      mat = apply(theta[,which(prand>0),j+1], 1, con_sup_check)
      if (is.null(dim(mat))) mat = matrix(mat, N, length(which(prand>0)) - 1) else mat = t(mat)
      check[,which(prand>0)] = cbind(rep(0, N), mat)
      check[,which(prand==0)] = 0
      psup_out = abind(psup_out, apply(check, 2, mean), along = 2)
      if (length(y) < burn) {
        psup = abind(psup, psup[,j], along = 2)
      } else psup = abind(psup, apply(check, 2, mean), along = 2)
      if (response.type == 'continuous') fmat = apply(theta[,which(prand>0),j+1], 1, con_fut_check, MID = MID)
      if (response.type == 'binary') fmat = apply(p_new[,which(prand>0)], 1, con_fut_check, MID = MID)
      if (is.null(dim(fmat))) fmat = matrix(fmat, N, length(which(prand>0)) - 1) else fmat = t(fmat)
      fcheck[,which(prand>0)] = cbind(rep(0, N), fmat)
      fcheck[,which(prand==0)] = 0
      pfut = apply(fcheck, 2, mean)
      ff = NULL
      if (length(which(pfut > uppfut))>0) ff = which(pfut > uppfut)
      if (!is.null(ff)) {
        fut.stop = c(fut.stop, ff)
      }
      if (!is.null(fut.stop)) {
        psup[fut.stop,j+1] = 0
        }
      psup_out[,1] = c(0,rep(1/(nt-1), (nt - 1)))
      ntj = length(which(psup[,j+1]>0)) + 1
      prand = rep(0, nt)
      prand[which(psup[,j+1]>0)]  = 1
      prand[1] = 1
      if (length(y)>=burn) prand[which(pfut > uppfut)] = 0
    }
    
    if ((((max(psup[,j+1]) > upper) | sum(prand>0) <= 1  | length(y) >= maxN)) & length(y)>=burn) break
  }
    if (response.type == 'continuous') {
      p.est = post0[,1]
      low = apply(post0, 1, function(z) qnorm(.025, z[1], z[2]))
      up = apply(post0, 1, function(z) qnorm(.975, z[1], z[2]))
    }
    if (response.type == 'binary') {
      p.est = post0[,1]/apply(post0,1,sum)
      low = apply(post0, 1, function(z) qbeta(.025, z[1], z[2]))
      up = apply(post0, 1, function(z) qbeta(.975, z[1], z[2]))
    }
    est = data.frame(p.est = p.est, low = low, up = up)
  if (good.out == F) {
    if (response.type == 'continuous') {
      theta = - theta
      est = -est
      low = est$low
      est$low = est$up
      est$up = low
    }
    if (response.type == 'binary') {
      theta = 1 - theta
      est = 1 - est
      low = est$low
      est$low = est$up
      est$up = low
    }
  }
  out = list(psup0 = psup, psup = psup_out, theta = theta, est = est, y = y, x = x[,-1], x0 = x0[,-1])
  class(out) = 'trial'
  return(out)
}





beta = function(x, m, v) c(x[1]/(x[1]+x[2]) - m, x[1]*x[2]/((x[1]+x[2])^2*(x[1]+x[2]+1)) - v)
library(nleqslv)
nleqslv(c(1,1), beta, m = .03, v = .1)
