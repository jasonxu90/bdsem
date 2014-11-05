######################################
##### functions for solving ODEs #####
######################################
#this is the ODE for our probability generating function (pgf) used for calculating transition probabilities

# arguments: 
# t is a number or vector of times to be evaluated
# param is a vector containing b/s/d rates lam,v,mu and argument s2
# state contains the state vectors at their initial conditions
# returns the rate of change in list form
# for further info see deSolve documentation/vignette

#' Transition probability ODE
#' 
#' Evaluates the ODE for generating function corresponding to transition probabilities. 
#' This is in a format to be solved using zvode from package deSolve; see deSolve documentation/vignette
#' for further details
#' 
#' @param t A number or vector of numbers: evaluation times of ODE
#' @param state A named vector containing initial value of the state variable G, complex valued
#' @param param A named vector of numbers containing the other arguments lam, v, mu, and complex number s2
#' @return The rate of change of G the generating function in list form
#' @examples
#' t = .5; state = c(G=exp(2*1i)); param = c(lam = .2, v = .05, mu = .1, s2 = exp(3*1i))
#' de.trans(t,state,param)
de.trans <- function(t, state, param){  
  with(as.list(c(state,param)), {    
    C <- 1/(s2-1+.5*.Machine$double.eps) + lam/(lam-mu)
    f <- 1 + 1/(lam/(mu-lam) + C*exp((mu-lam)*t))
    dG <- G*(lam*f - lam - v - mu) + v*f + mu
    #return rates of change
    list(dG)
  }) #end with as.list
}

#the following are the ODEs to be used for mean 
#births, mean deaths, mean shifts, and particle times (restricted moments)

#' Expected shifts ODE
#' 
#' Evaluates the ODE for generating function corresponding to transition probabilities. 
#' This is in a format to be solved using zvode from package deSolve; see deSolve documentation/vignette
#' for further details
#' 
#' @param t A number or vector of numbers: evaluation times of ODE
#' @param state A named vector containing initial value of the state variable G, complex valued
#' @param param A named vector of numbers containing the other arguments lam, v, mu, r, and complex number s2
#' @return The rate of change of G, the generating function for expected shifts, in list form
#' @examples
#' t = .5; state = c(G=exp(2*1i)); param = c(lam = .2, v = .05, mu = .1, r = 3, s2 = exp(3*1i))
#' de.shift(t,state,param)
de.shift <- function(t, state, param){  
  with(as.list(c(state,param)), {    
    C <- 1/(s2-1+.5*.Machine$double.eps) + lam/(lam-mu)
    f <- 1 + 1/(lam/(mu-lam) + C*exp((mu-lam)*t))
    dG <- G*(lam*f - lam - v - mu) + v*r*f + mu
    #return rates of change
    list(dG)
  }) #end with as.list
}

#' Expected births ODE
#' 
#' Evaluates the ODE for generating function corresponding to transition probabilities. 
#' This is in a format to be solved using zvode from package deSolve; see deSolve documentation/vignette
#' for further details
#' 
#' @param t A number or vector of numbers: evaluation times of ODE
#' @param state A named vector containing initial value of the state variable G, complex valued
#' @param param A named vector of numbers containing the other arguments lam, v, mu, r, and complex number s2
#' @return The rate of change of G, the generating function for expected births, in list form
#' @examples
#' t = .5; state = c(G=exp(2*1i)); param = c(lam = .2, v = .05, mu = .1, r = 3, s2 = exp(3*1i))
#' de.birth(t,state,param)
de.birth <- function(t, state, param){
  with(as.list(c(state,param)), {    
    yb <- (lam + mu + sqrt(lam^2 + 2*lam*mu + mu^2 - 4*lam*mu*r + .5*.Machine$double.eps)) / (2*lam*r)
    C <- 1/(s2-yb) + lam*r/(2*lam*r*yb - lam - mu)
    f <- yb + 1/( -lam*r/(2*lam*r*yb - lam - mu) + C*exp(-(2*yb*lam*r - lam - mu)*t) )
    dG <- G*(lam*r*f - lam - v - mu) + v*f + mu
    #return rates of change
    list(dG)
  }) #end with as.list
}

#' Expected deaths ODE
#' 
#' Evaluates the ODE for generating function corresponding to transition probabilities. 
#' This is in a format to be solved using zvode from package deSolve; see deSolve documentation/vignette
#' for further details
#' 
#' @param t A number or vector of numbers: evaluation times of ODE
#' @param state A named vector containing initial value of the state variable G, complex valued
#' @param param A named vector of numbers containing the other arguments lam, v, mu, r, and complex number s2
#' @return The rate of change of G, the generating function for expected deaths, in list form
#' @examples
#' t = .5; state = c(G=exp(2*1i)); param = c(lam = .2, v = .05, mu = .1, r = 3, s2 = exp(3*1i))
#' de.death(t,state,param)
de.death <- function(t, state, param){
  with(as.list(c(state,param)), {    
    yd <- (lam + mu + sqrt(lam^2 + 2*lam*mu + mu^2 - 4*lam*mu*r + .5*.Machine$double.eps)) / (2*lam)
    C <- 1/(s2-yd) + lam/(2*lam*yd - lam - mu)
    f <- yd + 1/( -lam/(2*lam*yd - lam - mu) + C*exp(-(2*yd*lam - lam - mu)*t) )
    dG <- G*(lam*f - lam - v - mu) + v*f + mu*r
    #return rates of change
    list(dG)
  }) #end with as.list
}

#' Expected particle time ODE
#' 
#' Evaluates the ODE for generating function corresponding to transition probabilities. 
#' This is in a format to be solved using zvode from package deSolve; see deSolve documentation/vignette
#' for further details
#' 
#' @param t A number or vector of numbers: evaluation times of ODE
#' @param state A named vector containing initial value of the state variable G, complex valued
#' @param param A named vector of numbers containing the other arguments lam, v, mu, r, and complex number s2
#' @return The rate of change of G, the generating function for expected particle time, in list form
#' @examples
#' t = .5; state = c(G=exp(2*1i)); param = c(lam = .2, v = .05, mu = .1, r = 3, s2 = exp(3*1i))
#' de.particleT(t,state,param)
de.particleT <- function(t, state, param){
  with(as.list(c(state,param)), {    
    y <- (lam + mu + r +sqrt( (lam + mu + r)^2 - 4*lam*mu + .5*.Machine$double.eps)) / (2*lam)
    C <- 1/(s2-y+.Machine$double.eps) + lam/(2*lam*y - lam - mu - r)
    f <- y + 1/( -lam/(2*lam*y - lam - mu - r) + C*exp(-(2*y*lam - lam - mu - r)*t) )
    dG <- G*(lam*f - lam - v - mu - r) + v*f + mu
    #return rates of change
    list(dG)
  }) #end with as.list
}

####################################
###The following solve these pgfs###
####################################

#solves the pgf for transitiions: param needs to be a vector of (lam,v,mu), time is the interval length
#numsteps is steps taken by deSolve, s1 s2 are arguments

#' Transition ODE solver
#' 
#' Numerically solves the ODE defined in \code{\link{de.trans}} using \code{zvode} from package 
#' \code{deSolve} given evaluation time, initial state values s1 and s2, and birth/shift/death rates
#' 
#' @param time A number corresponding to the desired evaluation time of ODE
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1 A complex number giving the initial value of the ODE G
#' @param s2 A complex number
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return The function value of the transition probability generating function 
#' @examples 
#' time = 5;  dt = 5; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4
#' solve.trans(time,dt,s1,s2,lam,v,mu)
solve.trans <- function(time, dt, s1, s2, lam,v,mu){
  t <- seq(0, time, by=dt)
  state <- c(G = s1)
  param <- c(lam=lam, v=v, mu=mu, s2=s2)
  out <- zvode(y = state, times = t, func = de.trans, parms = param, atol = 1e-10, rtol = 1e-10)
  tail(out, n=1)[2]
}


#' Expected shifts ODE solver
#' 
#' Numerically solves the ODE defined in \code{\link{de.shift}} using \code{zvode} from package 
#' \code{deSolve} given evaluation time, initial state values s1 and s2, and birth/shift/death rates
#' 
#' @param time A number corresponding to the desired evaluation time of ODE
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1 A complex number giving the initial value of the ODE G
#' @param s2 A complex number
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @param r A real number
#' @return The function value of the shifts generating function 
#' @examples 
#' time = 5;  dt = 1; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4; r = 2
#' solve.shift(time,dt,s1,s2, r, lam,v,mu)
solve.shift <- function(time, dt, s1, s2, r, lam, v, mu){
  t <- seq(0, time, by=dt)
  state <- c(G = s1)
  param <- c(lam=lam, v=v, mu=mu, s2=s2, r = r)
  out <- zvode(y = state, times = t, func = de.shift, parms = param, atol = 1e-10, rtol = 1e-10)
  tail(out, n=1)[2]
}


#' Expected births ODE solver
#' 
#' Numerically solves the ODE defined in \code{\link{de.birth}} using \code{zvode} from package 
#' \code{deSolve} given evaluation time, initial state values s1 and s2, and birth/shift/death rates
#' 
#' @param time A number corresponding to the desired evaluation time of ODE
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1 A complex number giving the initial value of the ODE G
#' @param s2 A complex number
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @param r a real number
#' @return The function value of the births generating function 
#' @examples 
#' time = 5;  dt = 1; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4; r = 3
#' solve.birth(time,dt,s1,s2,r,lam,v,mu)
solve.birth <- function(time, dt, s1, s2, r, lam, v, mu){
  t <- seq(0, time, by=dt)
  state <- c(G = s1)
  param <- c(lam=lam, v=v, mu=mu, s2=s2, r = r)
  out <- zvode(y = state, times = t, func = de.birth, parms = param, atol = 1e-10, rtol = 1e-10)
  tail(out, n=1)[2]
}

#' Expected deaths ODE solver
#' 
#' Numerically solves the ODE defined in \code{\link{de.death}} using \code{zvode} from package 
#' \code{deSolve} given evaluation time, initial state values s1 and s2, and birth/shift/death rates
#' 
#' @param time A number corresponding to the desired evaluation time of ODE
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1 A complex number giving the initial value of the ODE G
#' @param s2 A complex number
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @param r a real number
#' @return The function value of the deaths generating function 
#' @examples 
#' time = 5;  dt = 1; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4; r = 3
#' solve.death(time,dt,s1,s2,r,lam,v,mu)
solve.death <- function(time, dt, s1, s2, r, lam, v, mu){
  t <- seq(0, time, by=dt)
  state <- c(G = s1)
  param <- c(lam=lam, v=v, mu=mu, s2=s2, r = r)
  out <- zvode(y = state, times = t, func = de.death, parms = param, atol = 1e-10, rtol = 1e-10)
  tail(out, n=1)[2]
}

#' Expected births ODE solver
#' 
#' Numerically solves the ODE defined in \code{\link{de.particleT}} using \code{zvode} from package 
#' \code{deSolve} given evaluation time, initial state values s1 and s2, and birth/shift/death rates
#' 
#' @param time A number corresponding to the desired evaluation time of ODE
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1 A complex number giving the initial value of the ODE G
#' @param s2 A complex number
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @param r a real number
#' @return The function value of the particle time generating function 
#' @examples 
#' time = 5;  dt = 1; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4; r = 3
#' solve.particleT(time,dt,s1,s2,r,lam,v,mu)
solve.particleT <- function(time, dt, s1, s2, r, lam, v, mu){
  t <- seq(0, time, by=dt)
  state <- c(G = s1)
  param <- c(lam=lam, v=v, mu=mu, s2=s2, r = r)
  out <- zvode(y = state, times = t, func = de.particleT, parms = param, atol = 1e-10, rtol = 1e-10)
  tail(out, n=1)[2]
}


#' Evaluates transition probability ODE over 2D grid of arguments
#' 
#' Applies the function \code{\link{solve.trans}} to a grid of inputs s1, s2 for a fixed time.
#' 
#' @param time A number corresponding to the desired evaluation time of ODE
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.trans(time,dt,s1.seq,s2.seq,lam,v,mu)

makeGrid.trans <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- solve.trans(time, dt, s1.seq[i], s2.seq[j], lam,v,mu) #put the result in the matrix
    }
  }
  return(result)
}

#' Evaluates expected shifts ODE over 2D grid of arguments
#' 
#' Applies the function \code{\link{solve.shift}} to a grid of inputs s1, s2 at one fixed time and r=1
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.shift.r1(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.shift.r1 <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- solve.shift(time, dt, s1.seq[i], s2.seq[j], 1.0, lam,v,mu) #put the result in the matrix
    }
  }
  return(result)
}

#' Evaluates expected births ODE over 2D grid of arguments
#' 
#' Applies the function \code{\link{solve.birth}} to a grid of inputs s1, s2 at one fixed time and r=1
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.birth.r1(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.birth.r1 <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- solve.birth(time, dt, s1.seq[i], s2.seq[j], 1.0, lam,v,mu) #put the result in the matrix
    }
  }
  return(result)
}

#' Evaluates expected deaths ODE over 2D grid of arguments
#' 
#' Applies the function \code{\link{solve.death}} to a grid of inputs s1, s2 at one fixed time and r=1
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.death.r1(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.death.r1 <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- solve.death(time, dt, s1.seq[i], s2.seq[j], 1.0, lam,v,mu) #put the result in the matrix
    }
  }
  return(result)
}

#' Evaluates expected particle time ODE over 2D grid of arguments
#' 
#' Applies \code{\link{solve.particleT}} to a grid of inputs s1, s2 at one fixed time and r=0
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.particleT.r0(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.particleT.r0 <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- solve.particleT(time, dt, s1.seq[i], s2.seq[j], 0.0, lam,v,mu) #put the result in the matrix
    }
  }
  return(result)
}

#' Evaluates partial derivative of expected shifts ODE over 2D grid
#' 
#' Numerically partially differentiates the solution
#' given in \code{\link{solve.shift}} over a grid of input values s1, s2 at one fixed time, and r=1
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.shift.partial(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.shift.partial <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- grad(solve.shift, 1.0, time = time, dt = dt, s1 = s1.seq[i], s2 = s2.seq[j], lam = lam, v = v, mu=mu, method = "simple")
    }
  }
  return(result)
}

#' Evaluates partial derivative of expected births ODE over 2D grid
#' 
#' Numerically partially differentiates the solution
#' given in \code{\link{solve.birth}} over a grid of input values s1, s2 at one fixed time, and r=1
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.birth.partial(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.birth.partial <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- grad(solve.birth, 1.0, time = time, dt = dt, s1 = s1.seq[i], s2 = s2.seq[j], lam = lam, v = v, mu=mu, method = "simple")
    }
  }
  return(result)
}

#' Evaluates partial derivative of expected deaths ODE over 2D grid
#' 
#' Numerically partially differentiates the solution
#' given in \code{\link{solve.death}} over a grid of input values s1, s2 at one fixed time, and r=1
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.death.partial(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.death.partial <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- grad(solve.death, 1.0, time = time, dt = dt, s1 = s1.seq[i], s2 = s2.seq[j], lam = lam, v = v, mu=mu, method = "simple")
    }
  }
  return(result)
}

#' Evaluates partial derivative of expected particle time ODE over 2D grid
#' 
#' Numerically partially differentiates the solution
#' given in \code{\link{solve.particleT}} over a grid of input values s1, s2 at one fixed time, and r=0
#' 
#' @param time A number corresponding to the desired evaluation time of ODEs
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param lam Birth rate
#' @param v Shift rate
#' @param mu Death rate
#' @return A matrix of dimension length(s1.seq) by length(s2.seq) of the function values
#' 
#' @examples 
#' time = 5;  dt = 5; lam = .5; v = .2; mu = .4
#' gridLength = 32
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' makeGrid.particleT.partial(time,dt,s1.seq,s2.seq,lam,v,mu)
makeGrid.particleT.partial <- function(time, dt, s1.seq, s2.seq, lam, v, mu){
  result <- matrix(nrow = length(s1.seq), ncol = length(s2.seq)) #this will store the grid
  for(i in 1:length(s1.seq)){
    for(j in 1:length(s2.seq)){    
      result[i,j] <- grad(solve.particleT, 0.0, time = time, dt = dt, s1 = s1.seq[i], s2 = s2.seq[j], lam = lam, v = v, mu=mu, method = "simple")
    }
  }
  return(result)
}

############################################################################
### The following functions create transition probability approximations ###
### by transforming solutions from makeGrid functions via FFT            ###
############################################################################

#' Compute transition probabilities over a grid at a list of evaluation times
#' 
#' Transforms the output matrix from \code{\link{makeGrid.trans}} to interpretable 
#' transition probabilities using \code{fft}. Does this at several input times given in tList.
#' Returns a list of matrices, where each list entry corresponds to a time in tList. 
#' The i,j entry of each matrix corresponds to the probability of transitioning from initNum type 1 particles 
#' and 0 type 2 particles to i type 1 particles, j type 2 particles, over the corresponding time length.
#' 
#' @param tList A vector of numbers corresponding to the desired evaluation times
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param initNum An integer giving the number of initial particles
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq), each list entry corresponds to an evaluation time from tList, each matrix is a matrix of transition probabilities
#' 
#' @examples 
#' tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getTrans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
getTrans.timeList <- function(tList, lam, v, mu, initNum, s1.seq,s2.seq, dt){
  tpm.list <- vector("list", length(tList))      #store the trans prob matrices
  for(i in 1:length(tList)){
    grid.de <- makeGrid.trans(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid.de <- grid.de^initNum
    fourier.de <- fft(grid.de)/length(grid.de)
    tpm.list[[i]] <- Re(fourier.de)
  }
  return(tpm.list)
}

#' Compute expected particle time over a grid at a list of evaluation times
#' 
#' Transforms the output matrices from \code{\link{makeGrid.particleT.r0}} and \code{\link{makeGrid.particleT.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several input times given in tList.
#' Returns a list of matrices, where each list entry corresponds to a time in tList. 
#' The i,j entry of each matrix corresponds to the expected particle time spent with i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, over the corresponding time length.
#' 
#' @param tList A vector of numbers corresponding to the desired evaluation times
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param initNum An integer giving the number of initial particles
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an evaluation time from tList
#' @examples 
#' tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getParticleT.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
getParticleT.timeList <- function(tList, lam, v, mu, initNum, s1.seq,s2.seq, dt){
  tpm.list <- vector("list", length(tList))      #store the trans prob matrices
  for(i in 1:length(tList)){
    grid.partial <- makeGrid.particleT.partial(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid <- makeGrid.particleT.r0(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid.de <- initNum*grid.partial*grid^(initNum-1)
    fourier.de <- fft(grid.de)/length(grid.de)
    tpm.list[[i]] <- Re(fourier.de)
  }
  return(tpm.list)
}

#' Compute expected deaths over a grid at a list of evaluation times
#' 
#' Transforms the output matrices from \code{\link{makeGrid.death.r1}} and \code{\link{makeGrid.death.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several input times given in tList.
#' Returns a list of matrices, where each list entry corresponds to a time in tList. 
#' The i,j entry of each matrix corresponds to the number of expected deaths when process has i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, at the end of corresponding time interval.
#' 
#' @param tList A vector of numbers corresponding to the desired evaluation times
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param initNum An integer giving the number of initial particles
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an evaluation time from tList
#' 
#' @examples 
#' tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getDeathMeans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
getDeathMeans.timeList <- function(tList, lam, v, mu, initNum, s1.seq,s2.seq, dt){
  tpm.list <- vector("list", length(tList))      #store the trans prob matrices
  for(i in 1:length(tList)){
    grid.partial <- makeGrid.death.partial(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid <- makeGrid.death.r1(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid.de <- initNum*grid.partial*grid^(initNum-1)
    fourier.de <- fft(grid.de)/length(grid.de)
    tpm.list[[i]] <- Re(fourier.de)
  }
  return(tpm.list)
}

#' Compute expected shifts over a grid at a list of evaluation times
#' 
#' Transforms the output matrices from \code{\link{makeGrid.shift.r1}} and \code{\link{makeGrid.shift.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several input times given in tList.
#' Returns a list of matrices, where each list entry corresponds to a time in tList. 
#' The i,j entry of each matrix corresponds to the number of expected shifts when process has i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, at the end of corresponding time interval.
#' 
#' @param tList A vector of numbers corresponding to the desired evaluation times
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param initNum An integer giving the number of initial particles
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an evaluation time from tList
#' 
#' @examples 
#' tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getShiftMeans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
getShiftMeans.timeList <- function(tList, lam, v, mu, initNum, s1.seq,s2.seq, dt){
  tpm.list <- vector("list", length(tList))      #store the trans prob matrices
  for(i in 1:length(tList)){
    grid.partial <- makeGrid.shift.partial(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid <- makeGrid.shift.r1(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid.de <- initNum*grid.partial*grid^(initNum-1)
    fourier.de <- fft(grid.de)/length(grid.de)
    tpm.list[[i]] <- Re(fourier.de)
  }
  return(tpm.list)
}

#' Compute expected birth over a grid at a list of evaluation times
#' 
#' Transforms the output matrices from \code{\link{makeGrid.birth.r1}} and \code{\link{makeGrid.birth.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several input times given in tList.
#' Returns a list of matrices, where each list entry corresponds to a time in tList. 
#' The i,j entry of each matrix corresponds to the number of expected births when process has i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, at the end of corresponding time interval.
#' 
#' @param tList A vector of numbers corresponding to the desired evaluation times
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param initNum An integer giving the number of initial particles
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an evaluation time from tList
#' 
#' @examples 
#' tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getBirthMeans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
getBirthMeans.timeList <- function(tList, lam, v, mu, initNum, s1.seq,s2.seq, dt){
  tpm.list <- vector("list", length(tList))      #store the trans prob matrices
  for(i in 1:length(tList)){
    grid.partial <- makeGrid.birth.partial(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid <- makeGrid.birth.r1(time = tList[i], dt, s1.seq, s2.seq, lam,v,mu)
    grid.de <- initNum*grid.partial*grid^(initNum-1)
    fourier.de <- fft(grid.de)/length(grid.de)
    tpm.list[[i]] <- Re(fourier.de)
  }
  return(tpm.list)
}

#' Compute transition probabilites over a grid at a list of initial particle counts
#' 
#' Transforms the output matrices from \code{\link{makeGrid.trans}} 
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several initial particle counts from initList
#' Returns a list of matrices, where each list entry corresponds to a number of initial particles
#' The i,j entry of each matrix corresponds to the probability of transitioning to i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, by the end of corresponding time interval.
#' 
#' @param initList A vector of integers corresponding to the desired initial particle counts
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param u A number giving the observation interval length, equivalently the time to evaluate the ODEs
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an initial number of particles from initList
#' 
#' @examples 
#' initList = c(10,11) 
#' #gives matrices of transition probabilities corresponding to 10 initial particles and 11 initial particles
#' u = 1; dt = 1; lam = .5; v = .2; mu = .4
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getTrans.initList(u, initList, lam, v, mu, s1.seq, s2.seq, dt)

getTrans.initList <- function(u, initList, lam, v, mu, s1.seq, s2.seq, dt){  #initList list of different initial sizes, and u = dt
  tpm.list <- vector("list", length(initList))      #store the trans prob matrices
  grid.de <- makeGrid.trans(u, dt, s1.seq, s2.seq, lam,v,mu) 
  for(i in 1:length(initList)){
    grid.temp <- grid.de^initList[i]
    tpm.list[[i]] <- Re(fft(grid.temp)/length(grid.temp))
  }
  return(tpm.list)
}

#below, we use the fact that d/dr f^k evaluated at r=1 is equal to k* d/dr(f at 1)*f(1)^(k-1)

#' Compute expected births over a grid at a list of initial particle counts
#' 
#' Transforms the output matrices from \code{\link{makeGrid.birth.r1}} and \code{\link{makeGrid.birth.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several initial particle counts from initList
#' Returns a list of matrices, where each list entry corresponds to a number of initial particles
#' The i,j entry of each matrix corresponds to the expected births when the process has i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, by the end of corresponding time interval.
#' 
#' @param initList A vector of integers corresponding to the desired initial particle counts
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param u A number giving the observation interval length, equivalently the time to evaluate the ODEs
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an initial number of particles from initList
#' 
#' @examples 
#' initList = c(10,11) 
#' #gives matrices of expected sufficient statistics corresponding to 10 initial particles and 11 initial particles
#' u = 1; dt = 1; lam = .5; v = .2; mu = .4
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getBirthMeans.initList(u, initList, lam, v, mu, s1.seq, s2.seq, dt)
getBirthMeans.initList <- function(u, initList, lam, v, mu, s1.seq, s2.seq, dt){  #initList list of different initial sizes, and u = dt
  tpm.list <- vector("list", length(initList))      #store the trans prob matrices
  grid.partial <- makeGrid.birth.partial(u, dt, s1.seq, s2.seq, lam,v,mu)
  grid <- makeGrid.birth.r1(u, dt, s1.seq, s2.seq, lam,v,mu)
  for(i in 1:length(initList)){
    grid.temp <- initList[i]*grid.partial*grid^(initList[i]-1)
    tpm.list[[i]] <- Re(fft(grid.temp)/length(grid.temp))
  }
  return(tpm.list)
}

#' Compute expected shifts over a grid at a list of initial particle counts
#' 
#' Transforms the output matrices from \code{\link{makeGrid.shift.r1}} and \code{\link{makeGrid.shift.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several initial particle counts from initList
#' Returns a list of matrices, where each list entry corresponds to a number of initial particles
#' The i,j entry of each matrix corresponds to the expected shifts when the process has i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, by the end of corresponding time interval.
#' 
#' @param initList A vector of integers corresponding to the desired initial particle counts
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param u A number giving the observation interval length, equivalently the time to evaluate the ODEs
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an initial number of particles from initList
#' 
#' @examples 
#' initList = c(10,11) 
#' #gives matrices of expected sufficient statistics corresponding to 10 initial particles and 11 initial particles
#' u = 1; dt = 1; lam = .5; v = .2; mu = .4
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getShiftMeans.initList(u, initList, lam, v, mu, s1.seq, s2.seq, dt)
getShiftMeans.initList <- function(u, initList, lam, v, mu, s1.seq, s2.seq, dt){  #initList list of different initial sizes, and u = dt
  tpm.list <- vector("list", length(initList))      #store the trans prob matrices
  grid.partial <- makeGrid.shift.partial(u, dt, s1.seq, s2.seq, lam,v,mu)
  grid <- makeGrid.shift.r1(u, dt, s1.seq, s2.seq, lam,v,mu)
  for(i in 1:length(initList)){
    grid.temp <- initList[i]*grid.partial*grid^(initList[i]-1)
    tpm.list[[i]] <- Re(fft(grid.temp)/length(grid.temp))
  }
  return(tpm.list)
}

#' Compute expected deaths over a grid at a list of initial particle counts
#' 
#' Transforms the output matrices from \code{\link{makeGrid.death.r1}} and \code{\link{makeGrid.death.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several initial particle counts from initList
#' Returns a list of matrices, where each list entry corresponds to a number of initial particles
#' The i,j entry of each matrix corresponds to the expected deaths when the process has i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, by the end of corresponding time interval.
#' 
#' @param initList A vector of integers corresponding to the desired initial particle counts
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param u A number giving the observation interval length, equivalently the time to evaluate the ODEs
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an initial number of particles from initList
#' 
#' @examples 
#' initList = c(10,11) 
#' #gives matrices of expected sufficient statistics corresponding to 10 initial particles and 11 initial particles
#' u = 1; dt = 1; lam = .5; v = .2; mu = .4
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getDeathMeans.initList(u, initList, lam, v, mu, s1.seq, s2.seq, dt)
getDeathMeans.initList <- function(u, initList, lam, v, mu, s1.seq, s2.seq, dt){  #initList list of different initial sizes, and u = dt
  tpm.list <- vector("list", length(initList))      #store the trans prob matrices
  grid.partial <- makeGrid.death.partial(u, dt, s1.seq, s2.seq, lam,v,mu)
  grid <- makeGrid.death.r1(u, dt, s1.seq, s2.seq, lam,v,mu)
  for(i in 1:length(initList)){
    grid.temp <- initList[i]*grid.partial*grid^(initList[i]-1)
    tpm.list[[i]] <- Re(fft(grid.temp)/length(grid.temp))
  }
  return(tpm.list)
}

#' Compute expected particle time over a grid at a list of initial particle counts
#' 
#' Transforms the output matrices from \code{\link{makeGrid.particleT.r0}} and \code{\link{makeGrid.particleT.partial}}
#' to a list of expected sufficient statistics matrices using \code{fft}. Does this at several initial particle counts from initList
#' Returns a list of matrices, where each list entry corresponds to a number of initial particles
#' The i,j entry of each matrix corresponds to the expected particle time spent with i type 1 particles and j type 2 particles,
#' beginning with initNum type 1 particles, by the end of corresponding time interval.
#' 
#' @param initList A vector of integers corresponding to the desired initial particle counts
#' @param dt A number giving the increment length used in solving the ODE
#' @param s1.seq A vector of complex numbers; initial values of the ODE G
#' @param s2.seq A vector of complex numbers as inputs of s2.seq
#' @param u A number giving the observation interval length, equivalently the time to evaluate the ODEs
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @return A list of matrices of dimension length(s1.seq) by length(s2.seq); each list entry corresponds to an initial number of particles from initList
#' 
#' @examples 
#' initList = c(10,11) 
#' #gives matrices of expected sufficient statistics corresponding to 10 initial particles and 11 initial particles
#' u = 1; dt = 1; lam = .5; v = .2; mu = .4
#' gridLength = 16
#' s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
#' getParticleT.initList(u, initList, lam, v, mu, s1.seq, s2.seq, dt)
getParticleT.initList <- function(u, initList, lam, v, mu, s1.seq, s2.seq, dt){  #initList list of different initial sizes, and u = dt
  tpm.list <- vector("list", length(initList))      #store the trans prob matrices
  grid.partial <- makeGrid.particleT.partial(u, dt, s1.seq, s2.seq, lam,v,mu)
  grid <- makeGrid.particleT.r0(u, dt, s1.seq, s2.seq, lam,v,mu)
  for(i in 1:length(initList)){
    grid.temp <- initList[i]*grid.partial*grid^(initList[i]-1)
    tpm.list[[i]] <- -Re(fft(grid.temp)/length(grid.temp))
    #notice the negative sign
  }
  return(tpm.list)
}


#############################
### Simulation Functions  ###
#############################
#the following function simulates from the birth-shift-death process:
#the "indices" correspond to locations; we just let S=100,000
#returns information at end of time interval as ordered pair (new,old)

#' Simulate a birth-shift-death-process
#' 
#' This function simulates one realization of a birth-shift-death CTMC with per-particle birth, death, and shift rates
#' lambda, nu, and mu. 
#' 
#' The process begins with initNum particles, each assigned a random location indexed by a number between 1 and 100,000.
#' Simulation occurs until time t.end. The process returns an ordered pair giving the number of initial locations (indices)
#' still occupied, followed by the number of new indices present, by t.end.
#' 
#' @param t.end A number giving the length of time for the simulation
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @param initNum An integer, initial total number of particles
#' @return A pair of integers giving the number of initial indices still present followed by the number of
#' new indices present in the population. Otherwise, returns '999' as an error code
#' 
#' @examples
#' simulate.true(2,.2,.12,.15,10)
simulate.true <- function(t.end, lam, v, mu, initNum){
  initInd <- sample(1:100000,initNum)            #locations of initial copy generation
  K <- length(initInd)                          #number of particles at current time
  ind <- initInd                                #indices/locations of current collection
  t.cur <- copies <- shifts <- deaths <- 0
  for(i in 1:99999999){                          
    if(K == 0){
      return(c(0,0))  #entire process dies
    }
    
    copy <- shift <- death <- 0   #rates
    copy <-lam*K; shift <- v*K; death <- mu*K; rates <- c(copy, shift, death)
    t.next <- rexp(1, sum(rates)) #time of next event
    t.cur <- t.cur+t.next
    
    if(t.cur > t.end){ #end and return the population
      overlap <- length(intersect(ind,initInd))
      return( c(overlap, length(ind)-overlap) )  #common elements are number of old/type 1 particles; rest are new type 2
    }
    
    decision <- sample( c(1,2,3), 1, prob = rates/sum(rates))  #ratios determine type of event that occurs
    if(decision == 1){ ind <- c(ind, sample(1:100000,1))       #add a new location or index
                       copies <- copies + 1}
    else if(decision == 2){ ind[sample(1:length(ind),1)] <- sample(1:100000,1) #replace an index randomly
                            shifts <- shifts + 1}
    else if(decision == 3){ ind <- ind[-sample(1:length(ind),1)]  #delete an index randomly
                            deaths <- deaths + 1}
    K <- length(ind)  #update the overall rate
    
  }
  return(999) #for catching errors
}


#' Simulate a BSD process with error catch
#'
#' Wrapper for \code{\link{simulate.true}} that catches any possible errors
#' @inheritParams simulate.true
#' @return A pair of integers giving the number of initial indices still present followed by the number of
#' new indices present in the population
#' @examples
#' sim.one.true(2,.2,.12,.15,10)
sim.one.true <- function(t.end, lam, v, mu, initNum){
  res = 999 #this is a check so that nothing weird happens
  while(res[1] == 999){
    res <- simulate.true(t.end, lam, v, mu, initNum)  }
  return(res)
}

#' Simulate N realizations from the BSD process
#'
#' Repeats the function \code{\link{sim.one.true}} N times using \code{replicate}
#' 
#' @inheritParams simulate.true
#' @param N An integer, the number of realizations to simulate
#' @return A 2 by N matrix; each \emph{i}th column is a pair of integers giving the number of initial indices still present followed by the number of
#' new indices present in the population in the \emph{i}th realization
#' @examples
#' sim.N.true(10,2,.2,.12,.15,10)
sim.N.true<- function(N, t.end,lam,v,mu,initNum){
  replicate(N, sim.one.true(t.end, lam, v, mu, initNum))
}

#this simulation is for use with random initial starting sizes across replications;
#same as sim.one.true but keeps track of initial number as well in return statement

#' Simulate a BSD process, returning initial number as well as end state
#'
#' Wrapper for \code{\link{simulate.true}} that catches any possible errors and returns the initial number.
#' 
#' @inheritParams simulate.true
#' @return A vector containing the initial number of particles, followed by a pair of integers giving the number of initial indices still present followed by the number of
#' new indices present in the population
#' @examples
#' sim.one.ran(2,.2,.12,.15,10)
sim.one.ran <- function(t.end, lam, v, mu, initNum){  #also stores initial number
  res = 999 #this is a check so that nothing weird happens
  while(res[1] == 999){
    res <- simulate.true(t.end, lam, v, mu, initNum)  }
  return(c(initNum,res))
}

#' Simulate N realizations from the BSD process, each with random initial number each time
#'
#' Repeats \code{\link{sim.one.ran}} N times using \code{replicate}. Each replication 
#' begins with a random number of initial particles uniformly generated from a specified range.
#' 
#' @inheritParams simulate.true
#' @param N An integer, the number of realizations to simulate
#' @param range A vector containing possible initial populations, typically a sequence of integers
#' @return A 3 by N matrix; each \emph{i}th column corresponds to the \emph{i}th realization. The first row contains
#' initial numbers, the second row contains the number of original indices still present in the population by t.end,
#' and the third row contains the number of new indices present.
#' @examples
#' ransim.N.true(5,2,.2,.12,.15, seq(7,15))
ransim.N.true<- function(N, t.end,lam,v,mu,range){
  replicate(N, sim.one.ran(t.end, lam, v, mu, sample(range,1)))
}

##############################################################################
### The following functions simulate from the same true process, but return the number of shifts, 
### births, deaths as well as the ending state, necessary for verifying restricted moment calculations
### Not necessary for any inference in itself  
##############################################################################

#' Simulate a BSD process and return sufficient statistics
#' 
#' Simulates a birth-shift-death process and records the number of births, deaths, and shifts, as well as total
#' particle time, over the observation interval until time t.end. Used toward simulation functions that check
#' accuracy of expected restricted moment calculations
#' 
#' @inheritParams simulate.true
#' @return A vector of four numbers giving total copies, shifts, deaths, and particle time, respectively. Returns 
#' the integer '999' as an error code if error occurs.
#' @examples
#' sim.eventcount(2,.2,.12,.15,10)
sim.eventcount <- function(t.end, lam, v, mu, initNum){
  initInd <- sample(1:100000,initNum)      	    #locations of initial copy generation, S =10,000 for now
  K <- length(initInd)                          #number of particles at current time
  ind <- initInd                                #indices/locations of current collection
  t.cur <- copies <- shifts <- deaths <- particleT <- 0 # number of copies/shifts/deaths, and total particle time
  for(i in 1:9999999){
    if(K == 0){
      return(c(copies,shifts,deaths, particleT))	#entire process dies
    }
    
    copy <- shift <- death <- 0   #rates
    copy <-lam*K; shift <- v*K; death <- mu*K; rates <- c(copy, shift, death)
    t.next <- rexp(1, sum(rates))
    t.cur <- t.cur+t.next
    
    if(t.cur > t.end){ #end and return the population
      particleT <- particleT + K*(t.end - t.cur + t.next) #add the remaining chunk of time
      return(c(copies,shifts,deaths,particleT))
    }
    
    decision <- sample( c(1,2,3), 1, prob = rates/sum(rates))	#ratios determine type of event that occurs
    if(decision == 1){ ind <- c(ind, sample(1:100000,1))       #add a new location or index
                       copies <- copies + 1}
    else if(decision == 2){ ind[sample(1:length(ind),1)] <- sample(1:100000,1) #replace an index randomly
                            shifts <- shifts + 1}
    else if(decision == 3){ ind <- ind[-sample(1:length(ind),1)]  #delete an index randomly
                            deaths <- deaths + 1}
    particleT <- particleT + K*t.next
    K <- length(ind)  #update the overall rate
  }
  return(999)
}

#' Simulate a BSD process and return sufficient statistics with error catch
#' 
#' Simulates a birth-shift-death process and records the number of births, deaths, and shifts, as well as total
#' particle time, over the observation interval until time t.end. Used toward simulation functions that check
#' accuracy of expected restricted moment calculations
#' 
#' @inheritParams simulate.true
#' @return A vector of four numbers giving total copies, shifts, deaths, and particle time, respectively.
#' @examples 
#' sim.one.eventcount(2,.2,.12,.15,10)
sim.one.eventcount <- function(t.end, lam, v, mu, initNum){
  res = 999 #this is a check so that nothing weird happens
  while(res[1] == 999){
    res <- sim.eventcount(t.end, lam, v, mu, initNum)  }
  return(res)
}

#' Simulate N BSD processes and return sufficient statistics
#' 
#' Replicates \code{\link{sim.one.eventcount}} N times. Each replication 
#' begins with initNum particles.
#' 
#' @inheritParams simulate.true
#' @param N The number of replications, an integer
#' @return A 4 by N matrix, where rows correspond to total copies, shifts, deaths, and particle time per realization, respectively.
#' Each column corresponds to one realization of the process.
#' @examples
#' sim.N.eventcount(3,2,.2,.12,.15,10)
sim.N.eventcount<- function(N, t.end,lam,v,mu,initNum){
  replicate(N, sim.one.eventcount(t.end, lam, v, mu, initNum))
}


#generates a list of simulated data matrices, each entry of the list corresponds to a time in tList
#each column of each matrix stores one "interval", for N total intervals
#row one is initial number of [old] particles, row 2 is number of type 1 [old] at end of interval, 
#row 3 is type 2 [new] particles

#' Generate synthetic dataset for BSD process with simple rates
#' 
#' \code{makedata.simple} creates a list of synthetic observed datasets, each with N observation intervals
#' from the BSD process. Observation interval lengths are entries in tList. This is used in simulation studies checking
#' inferential procedures. Simulates observations using \code{\link{ransim.N.true}}.
#' 
#' @param N An integer specifying the number of observation intervals/realization from the simple BSD process.
#' @param tList A list of observation interval lengths. The number of datasets returned is equal to the length of tList
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @param initList A vector containing possible initial population sizes
#' @return A list of matrices. Each entry in the list is in the format returned by \code{\link{ransim.N.true}}, and 
#' corresponds to equal observation times specified in the corresponding entry in tList
#' @examples
#' makedata.simple(10,c(.2,.4,.6),.2,.12,.15,seq(8,12))
makedata.simple <- function(N, tList, lam, v, mu, initList){
  simDataList <- vector("list", length(tList))
  for(i in 1:length(tList)){
    simDataList[[i]] <- ransim.N.true(N, tList[i],lam,v,mu,initList)
  }
  return(simDataList)
}

#' Create example design matrix in the format required for \code{\link{MakePatientRates}}
#' 
#' This function creates one possible design matrix as an input for the function \code{\link{MakePatientRates}}
#' This particular design matrix features one covariate uniformly sampled between 6 and 10, as well as a row 
#' of 1's corresponding to an intercept term. The design matrix is thus 2 by number of patients
#' 
#' @param num.patients An integer giving the number of patients in the design matrix
#' @return A 2 by num.patients matrix: the first row is constant 1's corresponding to intercepts, and 
#' second row is a patient-specific covariate sampled uniformly between 6 and 10
#' 
PatientDesignExample <- function(num.patients){
  x1 <- runif(num.patients,6,10)
  patients.design <- rbind(1, x1) 
  return(patients.design)
}

#' Makes matrix of patient-specific birth, shift, and death rates
#' 
#' This function returns a matrix containing the birth, shift, and death rates of each patient. The rates are calculated
#' based on the log-linear relationship between regression coefficients 
#' \eqn{\mathbf{\beta} = (\mathbf{\beta^\lambda}, \mathbf{\beta}^\nu, \mathbf{\beta}^\mu)} 
#' and covariates \eqn{\mathbf{z_p}} in the \emph{p}th column of the design matrix given by
#' \eqn{\log(\lambda_p) = \mathbf\beta^\lambda \cdot \mathbf{z_p},  \log(\nu_p) = \mathbf\beta^\nu \cdot \mathbf{z_p}, \log(\mu_p) = \mathbf\beta^\mu \cdot \mathbf{z_p} }
#' 
#' @param patients.design A \eqn{n \times m} matrix, where n is number of covariates (including intercept) in the model and m is number of patients
#' @param betaVec The vector \eqn{\mathbf\beta} of regression coefficients. Note this function assumes 
#' that \eqn{ \mathbf\beta^\lambda, \mathbf\beta^\nu, \mathbf\beta^\mu} are equal in length. That is,
#' each rate depends on the same number of covariates
#' 
#' @return A \eqn{3 \times m} matrix where m is the number of patients, and rows correspond to birth, shift, and death rates respectively. 
#' 
#' @examples
#' num.patients = 10
#' patients.design <- PatientDesignExample(num.patients)
#' beta.lam <- c(log(8), log(.6)); beta.v <- c( log(.5), log(.7)); beta.mu <- c(log(.8), log(.8))
#' betas <- c(beta.lam,beta.v,beta.mu)
#' MakePatientRates(patients.design, betas)
MakePatientRates <- function(patients.design, betaVec){ 
  beta.lam <- betaVec[1:(length(betaVec)/3)]; beta.v = betaVec[(length(betaVec)/3+1):(2*length(betaVec)/3)]
  beta.mu = betaVec[(2*length(betaVec)/3 + 1) : length(betaVec)]
  return(exp( rbind( colSums(beta.lam*patients.design), colSums(beta.v*patients.design), colSums(beta.mu*patients.design)) ))
}

#' Generates synthetic patient data for inference in discretely observed BSD process with covariates
#' 
#' \code{MakePatientData} is the main function for generating a synthetic dataset with covariates. It simulates observations from
#' a birth-shift-death process for a number of synthetic "patients", using \code{\link{ransim.N.true}}. This provides data
#' for simulation studies toward assessing inference in a panel data setting, with rates depending on multiple covariates.
#' 
#' Each observation interval has a common fixed length given by the argument t. For each patient, between 2 and 6 observation intervals
#' are generated, and each begins with between 2 and 12 initial particles. These numbers are initialized uniformly at random, and
#' passed in as arguments to \code{\link{ransim.N.true}}, with the true rates set to the corresponding column of patients.rates.
#' 
#' @param t A number giving the observation interval length
#' @param num.patients An integer, the number of synthetic patients
#' @param patients.rates A matrix in the format returned by \code{\link{MakePatientRates}}
#' @return A \eqn{5 \times m} matrix where m is the number of total observation intervals generated. The first row
#' corresponds to the patient ID, the second row gives the initial number of particles for that observation interval,
#' the third through fifth column are in the format returned by \code{\link{ransim.N.true}}
#' 
#' @examples
#' num.patients = 10; t = .4
#' patients.design <- PatientDesignExample(num.patients)
#' beta.lam <- c(log(8), log(.6)); beta.v <- c( log(.5), log(.7)); beta.mu <- c(log(.8), log(.8))
#' betas <- c(beta.lam,beta.v,beta.mu)
#' patients.rates <- MakePatientRates(patients.design, betas)
#' MakePatientData(t, num.patients, patients.rates)
MakePatientData <- function(t, num.patients, patients.rates){ 
  res <- c()
  #first create design matrix  
  for(i in 1:num.patients){
    obs <- sample(seq(2,6),1) #how many observations we have for this patient, uniform between 2 and 6
    init.patient <- sample(seq(2,12),1)
    #let's say the possible values per patient are only + or - 1 from the init.patient, which is
    #kind of the "average" starting value for a given patient
    #this is more realistic and easier computationally
    res <- cbind(res, rbind(i, init.patient, ransim.N.true(obs, t, patients.rates[1,i], patients.rates[2,i], patients.rates[3,i], seq(init.patient-1, init.patient+1)) ) )
  }
  return(res)
}

#' Calculate empirical transition probabilities
#' 
#' Calculates transition probabilities empirically using Monte Carlo simulation
#' 
#' Note: function can be modified to initialize matrices to be larger sized in the case where rates are large
#' 
#' @param N An integer, the number of MC simulations per element of tList
#' @param tList A list of observation interval lengths to simulate
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @param initNum Integer giving the number of initial particles
#' @return A list of matrices, where each entry of the list corresponds to an element of tList. The i,j entry of
#' each matrix in the list give the probability of the process ending with i type 1 particles and j type 2 particles, beginning
#' with initNum type 1 particles, by the end of the corresponding observation length.
#' 
#' @examples
#' N = 20; tList = c(.5,1); lam = .2; v = .1; mu = .15; initNum = 10
#' getTrans.MC(N,tList,lam,v,mu,initNum)
getTrans.MC <- function(N, tList, lam,v,mu, initNum){
  tpm.list <- vector("list", length(tList))     #this will store 2 by 2 transitions
  for(t in 1:length(tList)){
    result <- sim.N.true(N, tList[t], lam, v, mu, initNum)
    trans.count <- matrix(0, (initNum+1),(initNum+20))  #make big enough to account for everything that may happen: count end states
    for(i in 1:N){
      id <- result[,i]+1 #indices in the resulting transition count: ie if you end at (1,1), you add a count to the (2,2) entry of the count matrix, etc
      trans.count[id[1], id[2]] = trans.count[id[1], id[2]] + 1
    }
    tpm.list[[t]] <- trans.count/sum(trans.count)
  }
  return(tpm.list)
}


#########################################################################################
### The following functions are relevant to EM algorithm and other simulation studies ###
#########################################################################################

##NOTE may include the real data version of functions where times can vary, etc

#' Log likelihood for BSD process with covariates
#' 
#' \code{logFFT.patients} evaluates the log likelihood of a dataset with observations corresponding to "patients" in the setting
#' where rates of the process depend on patient-specific covariates. The transition probabilities given the states of the process at 
#' endpoints of each observation interval are computed using the FFT/generating function method, relying on \code{\link{getTrans.initList}}. 
#' The log likelihood is then the sum of these log transition probabilities.
#' 
#' Note: this function is used so that MLE estimation of the coefficient vector can be accomplished, i.e. using \code{optim}, and 
#' is also used in numerically computing standard errors at the MLE. 
#' 
#' Vectors s1.seq and s2.seq should be of length greater than the total number of particles of either type at any observation interval
#' 
#' @param betas A vector of numbers \eqn{\mathbf\beta = (\mathbf\beta^\lambda, \mathbf\beta^\nu, \mathbf\beta^\mu)}
#' @param t.pat A number, the observation interval length
#' @param num.patients An integer, number of unique patients
#' @param PATIENTDATA A matrix in the form returned by \code{\link{MakePatientData}} containing the set of observation intervals
#' @param patients.design A design matrix in the same form as returned by \code{\link{PatientDesignExample}}
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @return The negative log likelihood of the observations in PATIENTDATA, given rates determined by coefficient vector betas and 
#' covariate values in patients.design
logFFT.patients <- function(betas, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq){
  b.lam <- betas[1:(length(betas)/3)]; b.v = betas[(length(betas)/3+1):(2*length(betas)/3)]; b.mu = betas[(2*length(betas)/3 + 1) : length(betas)]
  logresult <- 0
  for(i in 1:num.patients){
    lam <- exp( sum(b.lam*patients.design[,i])); v <- exp( sum(b.v*patients.design[,i])); mu <- exp(sum(b.mu*patients.design[,i]))
    #this takes the subset of the fake data corresponding to patient i
    dat.i <- PATIENTDATA[,which(PATIENTDATA[1,]==i)]
    tpm.list <- getTrans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq,t.pat) #stores trans matrices for all possible init size for this patient
    for(j in 1:dim(dat.i)[2]){
      initial <- dat.i[3,j] - dat.i[2,1] + 2  #subtract the init.patient to get the proper index in the tpm list
      n.old <- dat.i[4,j]
      n.new <- dat.i[5,j]
      logresult <- logresult + log(tpm.list[[initial]][n.old+1, n.new+1]) #initial is a list index, and then pick the entries in the matrix chosen
    }
  }
  return(-logresult)
}

#' Optimizes the function \code{\link{logFFT.patients}} using \code{optim} package
#' 
#' Function uses Nelder-Mead optimization as implemented in \code{optim} to maximize the log likelihood function
#' \code{\link{logFFT.patients}}. 
#' 
#' Example code is not included here because of runtime; see vignette for tutorial on using this function.
#' 
#' @param betaInit A vector, the initial guess for the algorithm
#' @param t.pat A number, the observation interval length
#' @param num.patients An integer, number of unique patients
#' @param PATIENTDATA A matrix in the form returned by \code{\link{MakePatientData}} containing the set of observation intervals
#' @param patients.design A design matrix in the same form as returned by \code{\link{PatientDesignExample}}
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @param tol A number for setting the relative tolerance for the algorithm (the reltol argument in \code{optim})
#' @param max An integer, the max number of iterations before termination
#' @return An \code{optim} type object
#optimizes the function logFFT.patients beginning with betaInit as initial guess
#uses Nelder-Mead optimization in optim package
FFT.pat.optim <- function(betaInit, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq, tol, max=2000, hess=TRUE){
  optim(betaInit, logFFT.patients, t.pat= t.pat, num.patients = num.patients, PATIENTDATA = PATIENTDATA, patients.design = patients.design, 
        s1.seq = s1.seq, s2.seq = s2.seq, hessian=hess, control = list(maxit = max, reltol = tol))
}


#' Perform one accelerated E-step of the EM algorithm
#' 
#' \code{ESTEP} performs one E-step of the EM algorithm, computing expected sufficient statistics given current settings
#' of the parameters. This function is the "accelerated" version, meaning that intervals with no observed changes are computed
#' more efficiently using closed form expressions, bypassing generating function and FFT calculations on these intervals.
#' 
#' @param betaVec A vector, the setting of beta coefficients
#' @param t.pat A number, the observation interval length
#' @param num.patients An integer, number of unique patients
#' @param PATIENTDATA A matrix in the form returned by \code{\link{MakePatientData}} containing the set of observation intervals
#' @param patients.design A design matrix in the same form as returned by \code{\link{PatientDesignExample}}
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @return A list containing a matrix of the expected sufficient statistics as well as the observed log likelihood value
ESTEP <- function(betaVec, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq){
  bLam <- betaVec[1:(length(betaVec)/3)]; bNu = betaVec[(length(betaVec)/3+1):(2*length(betaVec)/3)]; bMu = betaVec[(2*length(betaVec)/3 + 1) : length(betaVec)]
  U <- D <- V <- P <- rep(0, num.patients) #initialize vectors of expected sufficient statistics
  observedLikelihood <- 0
  for(i in 1:num.patients){
    #get lam, v, mu for current patients
    lam <- exp( sum(bLam*patients.design[,i]))
    v <- exp( sum(bNu*patients.design[,i])) 
    mu <- exp(sum(bMu*patients.design[,i]))
    #this takes the subset of the fake data corresponding to patient i
    dat.i <- PATIENTDATA[,which(PATIENTDATA[1,]==i)]
    
    #generate the transition probabilities and restricted moments for all possible initial sizes of patients
    #all information for the E step is calculated in these functions: this is the expensive part/meat of update
    tpm.list <- getTrans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat) #stores trans matrices for all possible init size for this patient
    restrictedBirth.list <- getBirthMeans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    restrictedShift.list <- getShiftMeans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    restrictedDeath.list <- getDeathMeans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    restrictedParticle.list <- getParticleT.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    
    #loop over number of observations associated with patient i
    for(j in 1:dim(dat.i)[2]){
      #these three items are indices to look up terms in our transition lists
      initial <- dat.i[3,j] - dat.i[2,1] + 2  #subtract the init.patient to get the proper index in the tpm list
      n.old <- dat.i[4,j]
      n.new <- dat.i[5,j]
      #initial is a list index, and then pick the entries in the matrix chosen
      #we divide all the restricted moments by the corresponding transition probability
      #we add to the i'th entry for each interval associated with patient i in this loop
      
      #if no changes occur, use quicker 'FM' methods
      if( dat.i[3,j] == dat.i[4,j] & n.new == 0){
        #print('no event')
        #zero expected counts to add, but particle time needs to be accounted for
        P[i] <- P[i] + dat.i[3,j]*t.pat   #particle time is just initial particles multiplied by interval length
        observedLikelihood <- observedLikelihood - dat.i[3,j]*(lam+v+mu)*t.pat
      } else {
        #NOTE we should check for any divide by zero etc that may happen here
        observedLikelihood <- observedLikelihood + log(tpm.list[[initial]][n.old+1, n.new+1])
        U[i] <- U[i] + restrictedBirth.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1] 
        D[i] <- D[i] + restrictedDeath.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1] 
        V[i] <- V[i] + restrictedShift.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1] 
        P[i] <- P[i] + restrictedParticle.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1]
        #note the P has a minus
      }
    }
  }
  result <- rbind(U,D,V,P)
  result[result < 0] = .000001 #numerical stability for negatives
  resultList <- list( 'matrix' = result, 'loglike' = observedLikelihood)
  return(resultList)  #return the vectors of expected sufficient statistics
}


#' Perform one E-step of the EM algorithm
#' 
#' \code{ESTEP.slow} performs one E-step of the EM algorithm, computing expected sufficient statistics given current settings
#' of the parameters. This function is the un-accelerated version, simply using the generating function approach to compute
#' necessary quantities for all observation intervals.
#' 
#' @inheritParams ESTEP
#' @return A list containing a matrix of the expected sufficient statistics as well as the observed log likelihood value
ESTEP.slow <- function(betaVec, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq){
  bLam <- betaVec[1:(length(betaVec)/3)]; bNu = betaVec[(length(betaVec)/3+1):(2*length(betaVec)/3)]; bMu = betaVec[(2*length(betaVec)/3 + 1) : length(betaVec)]
  U <- D <- V <- P <- rep(0, num.patients) #initialize vectors of expected sufficient statistics
  observedLikelihood <- 0
  for(i in 1:num.patients){
    #get lam, v, mu for current patients
    lam <- exp( sum(bLam*patients.design[,i]))
    v <- exp( sum(bNu*patients.design[,i])) 
    mu <- exp(sum(bMu*patients.design[,i]))
    #this takes the subset of the fake data corresponding to patient i
    dat.i <- PATIENTDATA[,which(PATIENTDATA[1,]==i)]
    #generate the transition probabilities and restricted moments for all possible initial sizes of patients
    #all information for the E step is calculated in these functions: this is the expensive part/meat of update
    tpm.list <- getTrans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat) #stores trans matrices for all possible init size for this patient
    restrictedBirth.list <- getBirthMeans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    restrictedShift.list <- getShiftMeans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    restrictedDeath.list <- getDeathMeans.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    restrictedParticle.list <- getParticleT.initList(t.pat, seq(dat.i[2,1]-1, dat.i[2,1]+1), lam, v, mu, s1.seq, s2.seq, t.pat)
    
    #loop over number of observations associated with patient i
    for(j in 1:dim(dat.i)[2]){
      #these three items are indices to look up terms in our transition lists
      initial <- dat.i[3,j] - dat.i[2,1] + 2  #subtract the init.patient to get the proper index in the tpm list
      n.old <- dat.i[4,j]
      n.new <- dat.i[5,j]
      #initial is a list index, and then pick the entries in the matrix chosen
      #we divide all the restricted moments by the corresponding transition probability
      #we add to the i'th entry for each interval associated with patient i in this loop
      
      #NOTE we should check for any divide by zero etc that may happen here
      observedLikelihood <- observedLikelihood + log(tpm.list[[initial]][n.old+1, n.new+1])
      U[i] <- U[i] + restrictedBirth.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1] 
      D[i] <- D[i] + restrictedDeath.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1] 
      V[i] <- V[i] + restrictedShift.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1] 
      P[i] <- P[i] + restrictedParticle.list[[initial]][n.old+1, n.new+1]/tpm.list[[initial]][n.old+1, n.new+1]
      #note the P has a minus
    }
  }
  result <- rbind(U,D,V,P)
  result[result < 0] = .000001 #numerical stability for negatives
  resultList <- list( 'matrix' = result, 'loglike' = observedLikelihood)
  return(resultList)  #return the vectors of expected sufficient statistics
}

#performs one newton raphson step: takes in matrix returned by E step and current estimate of parameters as arguments
#returns the new betas after a newton raphson update

#' Execute a Newton-Raphson step within M-step of the EM algorithm
#' 
#' \code{MSTEP} executes one iteration of a Newton-Raphson algorithm as part of the maximization (M-step) of the EM algorithm. 
#' Given the matrix of expected sufficient statistics returned by \code{\link{ESTEP}}, this function uses closed form gradient and
#' hessian expressions to efficiently optimize the current settings of the coefficients beta. This is called up to 10 times per M-step 
#' within \code{\link{EM.run}}
#' 
#' @param matrix A matrix in the format returned by \code{\link{ESTEP}} or \code{\link{ESTEP.slow}}
#' @param betaVec A vector of regression coefficients
#' @param num.patients An integer, the number of unique patients
#' @param patients.design A design matrix in the format generated by \code{\link{PatientDesignExample}}
#' @return An updated coefficient vector after one Newton-Raphson step
MSTEP <- function(matrix, betaVec, num.patients, patients.design){
  U <- matrix[1,]; D <- matrix[2,]; V <- matrix[3,]; P <- matrix[4,]
  bLam <- betaVec[1:(length(betaVec)/3)]; bNu = betaVec[(length(betaVec)/3+1):(2*length(betaVec)/3)]; bMu = betaVec[(2*length(betaVec)/3 + 1) : length(betaVec)]
  lamList <- vList <- muList <- rep(0, num.patients) #initialize vectors of expected sufficient statistics
  for(i in 1:num.patients){
    #get lam, v, mu for current patients
    lamList[i] <- exp( sum(bLam*patients.design[,i]))
    vList[i] <- exp( sum(bNu*patients.design[,i])) 
    muList[i] <- exp(sum(bMu*patients.design[,i]))
  }
  
  gradient <- c( patients.design %*% ( -diag(P) %*% lamList + U), patients.design %*% ( -diag(P) %*% vList + V), patients.design %*% ( -diag(P) %*% muList + D) )
  lamBlock <- -patients.design %*% diag(P) %*% diag(lamList) %*% t(patients.design)
  nuBlock <- -patients.design %*% diag(P) %*% diag(vList) %*% t(patients.design)
  muBlock <- -patients.design %*% diag(P) %*% diag(muList) %*% t(patients.design)
  hess <- bdiag(lamBlock,nuBlock,muBlock)
  hessInverse <- solve(hess)
  betaNew <- betaVec - as.vector( hessInverse %*% gradient )
  return(betaNew)
}

###runs the EM algorithm, with initial guess betaInit
#returns log likelihood and beta estimate at convergence, and number of iterations

#' Inference via EM algorithm
#' 
#' This function runs the EM algorithm from an initial guess. Infers the coefficient vector in setting where rates
#' depend on patient-specific covariates. The EM algorithm alternates between calling \code{\link{ESTEP}} and \code{\link{MSTEP}}
#' until the change in observed log-likelihood changes less than a specified relative tolerance between iterations
#' 
#' Examples are not included here due to runtime, but see vignette for usage.
#' 
#' @param betaInit A vector, the initial guess for coefficients beta
#' @param t.pat A number, the observation interval length
#' @param num.patients An integer, number of unique patients
#' @param PATIENTDATA A matrix in the form returned by \code{\link{MakePatientData}} containing the set of observation intervals
#' @param patients.design A design matrix in the same form as returned by \code{\link{PatientDesignExample}}
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @param relTol A number, the relative convergence criterion
#' 
#' @return A list containing the log-likelihood value at convergence, the final beta estimate, and the number of iterations
EM.run <- function(betaInit, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq, relTol, verbose=FALSE){
  llist <- c(999,99) #records log likelihoods, ignore these first two entries
  betacur <- bMat <- betaInit #betacur will keep track of initial estimate, bMat records all estimates
  mat <- mat.or.vec(4,num.patients) 
  iter = 2
  
  while(abs(llist[iter]-llist[iter-1]) > relTol*abs(llist[iter])){
    ###E STEP
    mat <- ESTEP(betacur, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq)
    if(verbose){print(betacur)}
    ###M step
    #initialize betaold, a vector to store previous beta setting after each Newton-Raphson step
    betaold <- rep(99,length(betacur)); numMsteps = 0
    #convergence criterion based on relTol for number of N-R steps, also max set at 10      
    while( sum(abs(betaold - betacur)) > relTol*sum(abs(betacur)) & numMsteps < 10){
      betaold <- betacur
      betacur <- MSTEP(mat$matrix,betacur, num.patients, patients.design)
      numMsteps = numMsteps + 1
    }
    #print(paste("M steps:", numMsteps))
    llist <- c(llist, mat$loglike )
    bMat <- rbind(bMat,as.vector(betacur))
    iter <- iter +1
  }
  
  #list of results: converged log like, final beta estimate, and number of iterations
  resList <- list('loglike' = llist[length(llist)], 'betaEstimate' =  bMat[dim(bMat)[1],] , 'iters' = iter-2)
  return(resList)
}

#calculates fisher info at the value of betas, and returns a vector of standard deviations

#' Calculate numerical standard errors of estimated coefficients
#' 
#' \code{getStandardErrors} uses \code{optimHess} to numerically compute the hessian at the value of the estimated MLE, 
#' which can be obtained using \code{\link{EM.run}}. The function inverts the hessian and takes the square root of diagonal entries
#' to yield standard errors for each coefficient estimate.
#' 
#' @inheritParams logFFT.patients
#' @return A vector the length of the coefficient vector, giving standard errors for each estimated coefficient
getStandardErrors <- function(betas, t.pat, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq){
  hessian <- optimHess(betas, logFFT.patients, t.pat = t.pat, num.patients = num.patients, PATIENTDATA = PATIENTDATA,
                       patients.design = patients.design, s1.seq = s1.seq, s2.seq = s2.seq)
  fisherInfo <- solve(hessian)
  return(sqrt(diag(fisherInfo)))
}


######################################
### Frequent Monitoring functions ####
######################################

#This calculates the four transition probabilities available under frequent monitoring, for comparison

#' Calculates the four transition probabilities defined under Frequent Monitoring
#' 
#' \code{getTrans.FreqMon} uses simple closed form expressions to compute transition probabilities available under the
#' Frequent Monitoring assumption, which allows at most one event to occur per observation interval. Computes the four
#' possible transition possibilities for intervals with lengths corresponding to entries in tList
#' 
#' @param tList A vector of observation interval lengths
#' @param lam Per-particle birth rate
#' @param v Per-particle shift rate
#' @param mu Per-particle death rate
#' @param initNum Integer, the number of initial particles
#' @return A list of \eqn{2 \times 2} matrices, each containing the probabilitiy of one birth, one shift, one death, or no event occurring
#' over a corresponding observation interval length from tList
#' 
#' @examples
#' getTrans.FreqMon(c(.5,1,10),.3,.1,.2,10)
getTrans.FreqMon <- function(tList, lam, v, mu, initNum){
  tpm.list <- vector("list", length(tList))     #this will store 2 by 2 transitions
  theta <- v + mu + lam
  for(u in 1:length(tList)){
    a <-  (1 - exp(-initNum*theta*tList[u]))*mu/theta  #top left (N,0) -> (N-1, 0)
    b <-  exp(-initNum*theta*tList[u])        #bottom left (N,0) -> (N,0)
    c <-  (1 - exp(-initNum*theta*tList[u]))*v/theta   #top right (N,0) -> (N-1,1)
    d <-  (1 - exp(-initNum*theta*tList[u]))*lam/theta   #bototm right (N,0) -> (N,1)
    tpm.list[[u]] <- matrix(c(a,b,c,d), 2, 2)      #append the matrix
  }
  return(tpm.list)
}

#the frequent monitoring log likelihood, u is length of interval, param is (lam, v, mu)
#this can be directly maximized via optim()
#takes in FM object

#' Calculate approximate log-likelihood based on Frequent Monitoring computations
#' 
#' \code{logFM} calculates the log-likelihood where transition probabilities are computed under the frequent monitoring 
#' assumption. Used for simulation studies comparing results using frequent monitoring with our methods.
#' 
#' @param param A vector of three numbers containing the birth, shift, and death rate respectively
#' @param u A number, the observation interval length
#' @param FM An object returned by \code{\link{FM.data}} containing the number of each event under FM
#' @return The negative frequent monitoring log-likelihood
logFM <- function(param, u, FM){
  FM.births <- FM$births; FM.deaths <- FM$deaths; FM.shifts <- FM$shifts; FM.nothing <- FM$nothing
  Bi <- length(FM.births)
  De <- length(FM.deaths)
  Sh <- length(FM.shifts)
  No <- length(FM.nothing)
  lam <- param[1]; v = param[2]; mu = param[3]
  theta <- lam+v+mu
  birthsum <- Bi*log(lam/theta) + sum( log(1 - exp(-FM.births*theta*u)) )
  deathsum <- De*log(mu/theta) + sum( log(1 - exp(-FM.deaths*theta*u)) )
  shiftsum <- Sh*log(v/theta) + sum( log(1 - exp(-FM.shifts*theta*u)) )
  nosum <- sum( log(exp(-FM.nothing*theta*u)))
  return( -(birthsum+deathsum+shiftsum+nosum ) )
}


##This function now returns the quantities necessary for logFM computation
#returns the number of births, deaths. shifts, and no event under FM in a vector called FM

#' Finds all intervals compatible under frequent monitoring assumption in a synthetic dataset
#' 
#' \code{FM.data} takes a simulated dataset in the format generated by \code{\link{makedata.simple}} and finds
#' the subset of observation intervals where at most one event occurred. This is necessary for computation of the
#' log likelihood in \code{\link{logFM}}
#' 
#' @param simDataList A list of synthetic observed datasets, returned by \code{\link{makedata.simple}}
#' @param u The index of the desired entry of simDataList
#' @return A list containing information about each type of FM event
FM.data <- function(simDataList, u){
  data <- simDataList[[u]]
  #takes only the intervals consistent under FM assumption:
  FM <- data[, intersect( which(data[1,] - data[2,] < 2), which(data[3,] < 2) )]  #only keep those that differ by at most 1 old particle
  #split up into births, deaths, shifts, or nothing happening, first row (initial number at that interval)
  FM.births <- FM[, intersect(which(FM[3,]== 1), which(FM[1,] == FM[2,])), drop = FALSE ][1,]
  FM.deaths <- FM[, intersect(which(FM[2,] < FM[1,]), which(FM[3,] == 0)), drop = FALSE ][1,]
  FM.shifts <- FM[, intersect(which(FM[2,] < FM[1,]), which(FM[3,] == 1)), drop = FALSE ][1,]
  FM.nothing <- FM[, intersect(which(FM[3,] == 0), which(FM[1,] == FM[2,])), drop = FALSE ][1,]
  resultList <- list('births' = FM.births, 'deaths' = FM.deaths, 'shifts' = FM.shifts, 'nothing' = FM.nothing)
  return(resultList)
}

#' Optimizes the frequent monitoring log-likelihood
#' 
#' This function converts a dataset generated in the format returned by \code{\link{makedata.simple}} 
#' and finds the intervals compatible under FM using \code{\link{FM.data}}, then infers the MLE birth, 
#' shift, and death rates.
#' 
#' @param simDataList A list of synthetic observed datasets, returned by \code{\link{makedata.simple}}
#' @param u The index of the desired entry of simDataList
#' @param initGuess A vector, initial guess for beta
#' @return An optim object corresponding to maximized frequent monitoring log-likelihood
FM.optim <- function(simDataList, u, initGuess){
  optim(initGuess, logFM, u=tList[u], FM = FM.data(simDataList, u), hessian = TRUE, control = list(trace = 6))
}

#takes in optim type object, returns vector of length 3 of coverage counts for each parameter
#for evaluating coverage post-hoc
#' Calculates coverage indicators for estimated parameters
#' 
#' This function takes in an object returned by \code{optim} and returns logicals on whether the corresponding 95%
#' confidence interval of each estimate contains the true value.
#' 
#' @param opt An optim object
#' @param trueParam A vector containing the true parameters
#' @return A vector of indicators, 1 indicating that the true parameter was included in the corresponding confidence interval
getCover <- function(opt, trueParam){
  sd <- sqrt(diag(solve(opt$hessian)))
  est <- opt$par
  as.numeric( est - 1.96*sd < trueParam)*as.numeric(est + 1.96*sd > trueParam )  #multiply numerics so that its in both interval endpoints
}

#' Get coefficient estimate from optim object
#' 
#' @param opt An optim object
#' @return The estimates contained in the optim object
getEst <- function(opt){
  opt$par
}

#main FM function, generates new fake data, returns a list with each entry an optim object
#each list entry corresponds to an element u from tList

#' Main function for generating observations from birth-shift-death process and inferring parameters under frequent monitoring
#' 
#' This function generates observation intervals from the simple birth-shift-death process without covariates, and generates
#' a dataset for each observation time length in tList. 
#' Next, the frequent monitoring log-likelihood is maximized using \code{\link{FM.optim}} for each dataset.
#' 
#' See accompanying vignette for example usage.
#' 
#' @inheritParams makedata.simple
#' @param initGuess Vector of numbers, initial guess for \code{optim}
#' @return A list of optim objects
FM.run <- function(N, tList, lam, v, mu,initList, initGuess){
  simDataList <- makedata.simple(N,tList,lam,v,mu,initList)  #simulate data
  result <- vector("list",length(tList))
  for(u in 1:length(tList)){
    result[[u]] <- FM.optim(simDataList, u, initGuess)
  }
  return(result)
}

#' Replicate \code{\link{FM.run}}
#' 
#' @param numReps The number of desired replications
#' @inheritParams FM.run
#' @return An array with entries of type returned by \code{\link{FM.run}}. Rows correspond to a dt value in tList, and
#' columns correspond to replications
#' 
#' @examples 
#' tList <- c(.2,.4,.6); initList <- c(1:15)
#' lam = .06; v = .02; mu = .11
#' trueParam <- c(lam,v,mu) 
#' N = 50; numReps = 2
#' example <- FM.replicate(numReps,N,tList,lam,v,mu,initList, trueParam)
FM.replicate <- function(numReps,N,tList,lam,v,mu,initList, initGuess){
  replicate(numReps, FM.run(N,tList,lam,v,mu,initList, initGuess))
}

### Repeat these for FFT method ###

#' Calculate log-likelihood using generating function technique 
#' 
#' \code{logFFT} calculates the log-likelihood where transition probabilities are computed using our FFT generating
#' function method
#' 
#' @param param A vector of three numbers containing the birth, shift, and death rate respectively
#' @param u A number, the observation interval length
#' @param simData An matrix in the format of a list entry returned by \code{\link{makedata.simple}}
#' @param initList A vector of possible initial populations
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @return The negative log-likelihood value
logFFT <- function(param, u, simData, initList, s1.seq, s2.seq){
  lam <- param[1]; v = param[2]; mu = param[3]
  tpm.list <- getTrans.initList(u, initList, lam, v, mu, s1.seq, s2.seq,.01) #this stores transition matrices for all initial sizes
  logresult <- 0
  for(i in 1:dim(simData)[2]){ #fake is a global holding the fake data
    initial <- simData[1,i]
    n.old <- simData[2,i]
    n.new <- simData[3,i]
    logresult <- logresult + log(tpm.list[[initial]][n.old+1,n.new+1])
  }
  return(-logresult)
}

#' Maximize the log-likelihood in \code{\link{logFFT}}
#' 
#' \code{FFT.optim} maximizes the log-likelihood using \code{optim}, with hessian = TRUE.
#' 
#' @param initGuess A vector containing the initial guess of birth, shift, and death rates
#' @param param A vector of three numbers containing the birth, shift, and death rate respectively
#' @param u A number, the observation interval length
#' @param simDataList A list of matrices in the format returned by \code{\link{makedata.simple}}
#' @param initList A vector of possible initial populations
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @return An optim type object 
FFT.optim <- function(simDataList, u, initGuess, initList, s1.seq, s2.seq){
  optim(initGuess, logFFT, u=tList[u], simData = simDataList[[u]], initList = initList, 
        s1.seq = s1.seq, s2.seq = s2.seq, hessian=T, control = list(trace = 6))
}

#' Generate synthetic data from simple BSD process and infer rates using generating function method
#' 
#' The main function for simulation studies assessing our generating function approach in the the discretely observed
#' simple birth-shift-death-process without covariates. Generates synthetic datasets using \code{\link{makedata.simple}}, and
#' infers the MLE rates for each using \code{\link{FFT.optim}}.
#' 
#' @inheritParams makedata.simple
#' @inheritParams FFT.optim
#' @return A list of optim objects
FFT.run <- function(N, tList, lam, v, mu, initList, initGuess, s1.seq, s2.seq ){
  simDataList <- makedata.simple(N,tList,lam,v,mu,initList)  #simulate data
  result <- vector("list",length(tList))
  for(u in 1:length(tList)){
    result[[u]] <- FFT.optim(simDataList, u, initGuess, initList, s1.seq, s2.seq)
  }
  return(result)
}

#' Replicate the function \code{\link{FFT.run}}
#' 
#' This function does the same as \code{\link{FM.replicate}}, but optimizes the likelihood based on transition
#' probabilities computed by the generating function method instead of using frequent monitoring. An example
#' of usage is included in the vignette.
#' 
#' @param numReps The number of replications
#' @inheritParams FFT.run
#' @return An array of optim objects in the same layout as \code{\link{FM.replicate}}
FFT.replicate <- function(numReps, N, tList, lam, v, mu, initList, initGuess, s1.seq, s2.seq){
  replicate(numReps, FFT.run(N,tList,lam,v,mu,initList,initGuess, s1.seq, s2.seq))
}


#' Perform one E-step of the EM algorithm for unevenly spaced observations
#' 
#' \code{ESTEP.realdata} performs one E-step of the EM algorithm similarly to \code{\link{ESTEP}},
#' but for datasets with unevenly spaced time intervals between observations. In particular, the second
#' row of the PATIENTDATA argument, a matrix returned by \code{\link{MakePatientData}} for simulated data
#' contains no relevant information for the algorithm. When working with real data, we may create a matrix
#' of the same format, except replace this row to instead contain information about the time between observations.
#' 
#' 
#' @param betaVec A vector, the setting of beta coefficients
#' @param num.patients An integer, number of unique patients
#' @param PATIENTDATA A matrix in the form returned by \code{\link{MakePatientData}} containing the set of observation intervals,
#' but the second row now contains the observation interval lengths 
#' @param patients.design A design matrix in the same form as returned by \code{\link{PatientDesignExample}}
#' @param s1.seq A vector of complex arguments evenly spaced along the unit circle
#' @param s2.seq A vector of complex arguments evenly spaced along the unit circle
#' @return A list containing a matrix of the expected sufficient statistics as well as the observed log likelihood value
ESTEP.realdata <- function(betaVec, num.patients, PATIENTDATA, patients.design, s1.seq, s2.seq){  
  bLam <- betaVec[1:(length(betaVec)/3)]; bNu = betaVec[(length(betaVec)/3+1):(2*length(betaVec)/3)]; bMu = betaVec[(2*length(betaVec)/3 + 1) : length(betaVec)]
  U <- D <- V <- P <- rep(0, num.patients) #initialize vectors of expected sufficient statistics
  observedLikelihood <- 0
  for(i in 1:num.patients){
    #get lam, v, mu for current patient
    lam <- exp( sum(bLam*patients.design[,i]))
    v <- exp( sum(bNu*patients.design[,i])) 
    mu <- exp(sum(bMu*patients.design[,i]))
    #this takes the subset of the fake data corresponding to patient i
    dat.i <- as.matrix(PATIENTDATA[,which(PATIENTDATA[1,]==UIDs[i])])
    #loop over number of observations associated with current patient:
    for(j in 1:dim(dat.i)[2]){
      #these three items are indices to look up terms in our transition lists
      n.initial <- dat.i[3,j]  #now this is actually initial number of particles, not an index
      n.old <- dat.i[4,j]
      n.new <- dat.i[5,j]
      timeInterval <- dat.i[2,j] #interval length
      #set the dt for numerical solver
      dt = timeInterval
      if( n.initial == n.old & n.new == 0){
        #print('no event in this interval')
        #zero expected counts to add, but particle time needs to be accounted for
        P[i] <- P[i] + n.initial*timeInterval  #particle time is just initial particles multiplied by interval length in row 2
        observedLikelihood <- observedLikelihood -n.initial*(lam+v+mu)*timeInterval
      } else {
        tpmMat <- getTrans(timeInterval, n.initial, lam, v, mu, s1.seq, s2.seq, dt)
        restrictedBirthMat <- getBirthMeans(timeInterval, n.initial, lam, v, mu, s1.seq, s2.seq, dt)
        restrictedDeathMat <- getDeathMeans(timeInterval, n.initial, lam, v, mu, s1.seq, s2.seq, dt)
        restrictedShiftMat <- getShiftMeans(timeInterval, n.initial, lam, v, mu, s1.seq, s2.seq, dt)
        restrictedParticleMat <- getParticleT(timeInterval, n.initial, lam, v, mu, s1.seq, s2.seq, dt)       
        #we divide all the restricted moments by the corresponding transition probability
        #we add to the i'th entry for each interval associated with patient i in this loop
        observedLikelihood <- observedLikelihood + log(tpmMat[n.old+1, n.new+1])
        U[i] <- U[i] + restrictedBirthMat[n.old+1, n.new+1]/tpmMat[n.old+1, n.new+1]
        D[i] <- D[i] + restrictedDeathMat[n.old+1, n.new+1]/tpmMat[n.old+1, n.new+1]
        V[i] <- V[i] + restrictedShiftMat[n.old+1, n.new+1]/tpmMat[n.old+1, n.new+1]
        P[i] <- P[i] + restrictedParticleMat[n.old+1, n.new+1]/tpmMat[n.old+1, n.new+1]
      }
    }
  }
  result <- rbind(U,D,V,P)
  result[result < 0] = .000001 #numerical stability for negatives
  resultList <- list( 'matrix' = result, 'loglike' = observedLikelihood)
  return(resultList)  #return the vectors of expected sufficient statistics
}