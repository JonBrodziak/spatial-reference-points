################################################################################################
# StockRecruitment
################################################################################################
StockRecruitment <- function(SB,parameters) {
  model = parameters[1]
  SB0 = parameters[2]
  R0 = parameters[3]
  h = parameters[4]
  R <- 4*h*R0*SB/(SB0*(1-h)+SB*(5*h-1))
  return(R)
}