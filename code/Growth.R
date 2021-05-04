Growth <- function(age,parameters) {
  model = parameters[1]
  if (model == 1) {
    amin = parameters[2]
    amax = parameters[3]
    offset = parameters[4]
    Lmin = parameters[5]
    Lmax = parameters[6]
    c = parameters[7]
    tmp <- 0.0
    tmp <- (1.0-c^(age+offset-amin))/(1.0-c^(amax-amin))
    length  <- Lmin+tmp*(Lmax-Lmin)
  }
  else if (model == 2) {
    offset = parameters[2]
    Linf = parameters[3]
    K = parameters[4]
    t0 = parameters[5]
    tmp <- 0.0
    tmp <- (1.0-exp(-K*(age+offset-t0)))
    length  <- tmp*Linf
  }
  return(length)
}
