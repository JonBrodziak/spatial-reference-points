################################################################################################
# Calculate_Survey_Selectivity_At_Age.R
################################################################################################

# Calculate survey selectivity at age by gender, area, and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (I in 1:NSurvey)
    if (SurveySelectivity.model[p,SurveyArea[I]] == 1) {
      parameters <- c(SurveySelectivity.model[p,SurveyArea[I]],SurveySelectivity.a50[p,SurveyArea[I],g],SurveySelectivity.slope[p,SurveyArea[I],g])
      for (g in 1:NGender) {
        tmp <- 0.0
        for (a in 1:NAge) {
          if (RecAge == 0)
            age <- a-1
          else if (RecAge == 1)
            age <- a
          SurveySelectivityAtAge[p,SurveyArea[I],g,a]  <- Selectivity(age,parameters)
        }
        # Rescale to set maximum selectivity at age to be 1      
        tmp <- max(SurveySelectivityAtAge[p,SurveyArea[I],g,])
        for (a in 1:NAge) {
          SurveySelectivityAtAge[p,SurveyArea[I],g,a] <- SurveySelectivityAtAge[p,SurveyArea[I],g,a]/tmp
          }
      }
    }

