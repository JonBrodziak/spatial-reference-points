################################################################################################
# Calculate_Maturity_Probability_At_Age.R
################################################################################################

for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
    {
      for (a in 1:NAge) 
      {
        if (RecAge == 0)
          age <- a-1
        else if (RecAge == 1)
          age <- a
        if (Maturity.model[p,d] == 1)
          parameters <- c(Maturity.model[p,d],Maturity.a50[p,d,g],Maturity.slope[p,d,g])
        MaturityProbabilityAtAge[p,d,g,a]  <- Maturity(age,parameters)
      }
    }
