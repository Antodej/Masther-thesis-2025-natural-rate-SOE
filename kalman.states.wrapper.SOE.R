#------------------------------------------------------------------------------#
# File:        kalman.states.wrapper.R
#
# Description: This is a wrapper function for kalman.states.R that specifies
#              inputs based on the estimation stage.
#------------------------------------------------------------------------------#
kalman.states.wrapper.SOE <- function(parameters, y.data, x.data, stage = NA,
                                  lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA,
                                  use.kappa=FALSE, kappa.inputs=NA, param.num){

  out <- unpack.parameters.SOE(parameters=parameters, y.data=y.data, x.data=x.data,
                                         xi.00=xi.00, P.00=P.00,
                                        use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)

  for (n in names(out)) {
      eval(parse(text=paste0(n, "<-out$", n)))
  }
  t.end <- dim(y.data)[1]
  # Run Kalman filter and smoother 
  states <- kalman.states(xi.tm1tm1=xi.00, P.tm1tm1=P.00, F=F, Q=Q, A=A, H=H, R=R, kappa=kappa.vec, y=y.data, x=x.data)
  if (stage == 1) {
      states$filtered$xi.tt <- states$filtered$xi.tt + cbind(1:t.end,0:(t.end-1),-1:(t.end-2)) * parameters[param.num["g"]]
      states$smoothed$xi.tT <- states$smoothed$xi.tT + cbind(1:t.end,0:(t.end-1),-1:(t.end-2)) * parameters[param.num["g"]]
  }
  return(states)
}