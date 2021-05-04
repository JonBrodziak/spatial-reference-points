################################################################################################
# Calculate_Fished_Equilibrium_Mean_Size_At_NAge.R
################################################################################################

# Calculate fished equilibrium lengths and weights at age 
# for the plus group, which is stored in array position NAge
# using the MaxPlusGroupAge adjustment
#-----------------------------------------------------------------
# MeanLengthStartOfYear
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender) 
    {
      tmp <- 0.0
      if (RecAge == 0) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,d] == 1)
            parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],StartZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
          else if (Growth.model[p,d] == 2)
            parameters <- c(Growth.model[p,d],StartZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,d,g,NAge])*Growth((NAge+k-1),parameters)
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,d] == 1)
            parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],StartZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
          else if (Growth.model[p,d] == 2)
            parameters <- c(Growth.model[p,d],StartZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,d,g,NAge])*Growth((NAge+k),parameters)
        }		  
      }
      tmp <- tmp*(1.0-exp(-EquilibriumTotalMortalityAtAge[p,d,g,NAge]))
      if (RecAge == 0)
        MeanLengthStartOfYear[p,d,g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,d,g,NAge]*(MaxPlusGroupAge-NAge+1)))
      else if (RecAge == 1)
        MeanLengthStartOfYear[p,d,g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,d,g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# MeanLengthSpawning
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender) 
    {
      tmp <- 0.0
      if (RecAge == 0) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,d] == 1)
            parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],SpawningZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
          else if (Growth.model[p,d] == 2)
            parameters <- c(Growth.model[p,d],SpawningZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,d,g,NAge])*Growth((NAge+k-1),parameters)
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,d] == 1)
            parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],SpawningZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
          else if (Growth.model[p,d] == 2)
            parameters <- c(Growth.model[p,d],SpawningZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,d,g,NAge])*Growth((NAge+k),parameters)
        }		  
      }
      tmp <- tmp*(1.0-exp(-EquilibriumTotalMortalityAtAge[p,d,g,NAge]))
      if (RecAge == 0)
        MeanLengthSpawning[p,d,g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,d,g,NAge]*(MaxPlusGroupAge-NAge+1)))
      else if (RecAge == 1)
        MeanLengthSpawning[p,d,g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,d,g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# MeanLengthCatch
for (p in 1:NPopulation)
  for (v in 1:NFleet)
    for (g in 1:NGender) 
    {
      tmp <- 0.0
      if (RecAge == 0) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,FleetArea[v]] == 1)
            parameters <- c(Growth.model[p,FleetArea[v]],Length.Amin[p,FleetArea[v],g],Length.Amax[p,FleetArea[v],g],CatchZFraction[p,FleetArea[v]],Length.Lmin[p,FleetArea[v],g],Length.Lmax[p,FleetArea[v],g],Length.c[p,FleetArea[v],g])
          else if (Growth.model[p,FleetArea[v]] == 2)
            parameters <- c(Growth.model[p,FleetArea[v]],CatchZFraction[p,FleetArea[v]],Length.Linf[p,FleetArea[v],g],Length.K[p,FleetArea[v],g],Length.t0[p,FleetArea[v],g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,FleetArea[v],g,NAge])*Growth((NAge+k-1),parameters) 
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,FleetArea[v]] == 1)
            parameters <- c(Growth.model[p,FleetArea[v]],Length.Amin[p,FleetArea[v],g],Length.Amax[p,FleetArea[v],g],CatchZFraction[p,FleetArea[v]],Length.Lmin[p,FleetArea[v],g],Length.Lmax[p,FleetArea[v],g],Length.c[p,FleetArea[v],g])
          else if (Growth.model[p,FleetArea[v]] == 2)
            parameters <- c(Growth.model[p,FleetArea[v]],CatchZFraction[p,FleetArea[v]],Length.Linf[p,FleetArea[v],g],Length.K[p,FleetArea[v],g],Length.t0[p,FleetArea[v],g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,FleetArea[v],g,NAge])*Growth((NAge+k),parameters) 
        }		  
      }
      tmp <- tmp*(1.0-exp(-EquilibriumTotalMortalityAtAge[p,FleetArea[v],g,NAge]))
      if (RecAge == 0)
        MeanLengthCatch[p,FleetArea[v],g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,FleetArea[v],g,NAge]*(MaxPlusGroupAge-NAge+1)))
      else if (RecAge == 1)
        MeanLengthCatch[p,FleetArea[v],g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,FleetArea[v],g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# MeanLengthSurvey
for (p in 1:NPopulation)
  for (I in 1:NSurvey)
    for (g in 1:NGender) 
    {
      tmp <- 0.0
      if (RecAge == 0) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,SurveyArea[I]] == 1)
            parameters <- c(Growth.model[p,SurveyArea[I]],Length.Amin[p,SurveyArea[I],g],Length.Amax[p,SurveyArea[I],g],SurveyZFraction[p,SurveyArea[I]],Length.Lmin[p,SurveyArea[I],g],Length.Lmax[p,SurveyArea[I],g],Length.c[p,SurveyArea[I],g])
          else if (Growth.model[p,SurveyArea[I]] == 2)
            parameters <- c(Growth.model[p,SurveyArea[I]],SurveyZFraction[p,SurveyArea[I]],Length.Linf[p,SurveyArea[I],g],Length.K[p,SurveyArea[I],g],Length.t0[p,SurveyArea[I],g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,SurveyArea[I],g,NAge])*Growth((NAge+k-1),parameters)
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,SurveyArea[I]] == 1)
            parameters <- c(Growth.model[p,SurveyArea[I]],Length.Amin[p,SurveyArea[I],g],Length.Amax[p,SurveyArea[I],g],SurveyZFraction[p,SurveyArea[I]],Length.Lmin[p,SurveyArea[I],g],Length.Lmax[p,SurveyArea[I],g],Length.c[p,SurveyArea[I],g])
          else if (Growth.model[p,SurveyArea[I]] == 2)
            parameters <- c(Growth.model[p,SurveyArea[I]],SurveyZFraction[p,SurveyArea[I]],Length.Linf[p,SurveyArea[I],g],Length.K[p,SurveyArea[I],g],Length.t0[p,SurveyArea[I],g])
          tmp <- tmp + exp(-k*EquilibriumTotalMortalityAtAge[p,SurveyArea[I],g,NAge])*Growth((NAge+k),parameters)
        }		  
      }
      tmp <- tmp*(1.0-exp(-EquilibriumTotalMortalityAtAge[p,SurveyArea[I],g,NAge]))
      if (RecAge == 0)
        MeanLengthSurvey[p,SurveyArea[I],g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,SurveyArea[I],g,NAge]*(MaxPlusGroupAge-NAge+1)))
      else if (RecAge == 1)
        MeanLengthSurvey[p,SurveyArea[I],g,NAge] <- tmp/(1.0-exp(-EquilibriumTotalMortalityAtAge[p,SurveyArea[I],g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# Compute fished equilibrium mean weights of the plus group for the start of year, 
# spawning, catch and survey by gender, area and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender) {
      if (WeightAtLength.model[p,d] == 1)
        parameters <- c(WeightAtLength.model[p,d],WeightAtLength.A[p,d,g],WeightAtLength.B[p,d,g])
      MeanWeightStartOfYear[p,d,g,NAge]  <-  WeightAtLength(MeanLengthStartOfYear[p,d,g,NAge],parameters)
      MeanWeightSpawning[p,d,g,NAge] <- WeightAtLength(MeanLengthSpawning[p,d,g,NAge],parameters)
    }
for (p in 1:NPopulation)
  for (v in 1:NFleet)
    for (g in 1:NGender) {
      if (WeightAtLength.model[p,d] == 1)
        parameters <- c(WeightAtLength.model[p,d],WeightAtLength.A[p,FleetArea[v],g],WeightAtLength.B[p,FleetArea[v],g])
      MeanWeightCatch[p,FleetArea[v],g,NAge] <-  WeightAtLength(MeanLengthCatch[p,FleetArea[v],g,NAge],parameters)
      }
for (p in 1:NPopulation)
  for (I in 1:NSurvey)
    for (g in 1:NGender) {
      if (WeightAtLength.model[p,d] == 1)
        parameters <-c(WeightAtLength.model[p,d],WeightAtLength.A[p,SurveyArea[I],g],WeightAtLength.B[p,SurveyArea[I],g])
      MeanWeightSurvey[p,SurveyArea[I],g,NAge] <-  WeightAtLength(MeanLengthSurvey[p,SurveyArea[I],g,NAge],parameters)
    }
