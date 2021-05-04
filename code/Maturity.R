################################################################################################
# Maturity
################################################################################################
Maturity <- function(age,parameters) {
  model = parameters[1]
  if (model == 1) {
    a50 = parameters[2]
    slope = parameters[3]
    fraction.mature <- 1.0/(1.0+exp(-(age-a50)/slope))
    return(fraction.mature)
  }
}