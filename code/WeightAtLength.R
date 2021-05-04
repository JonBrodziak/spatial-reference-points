################################################################################################
# WeightAtLength.R
################################################################################################
WeightAtLength <- function(length,parameters) {
  model = parameters[1]
  A = parameters[2]
  B = parameters[3]
  weight <- A*(length^B)
  return(weight)
}