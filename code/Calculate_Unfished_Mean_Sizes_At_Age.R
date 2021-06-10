################################################################################################
# Calculate_Unfished_Mean_Sizes_At_Age.R
################################################################################################

# Compute mean length at age at start of year by gender, area, and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
      for (a in 2:(NAge-1)) {
        if (RecAge == 0)
		      age <- a-1
		    else if (RecAge == 1)
		      age <- a
        if (Growth.model[p,d] == 1)
          parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],StartZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
        else if (Growth.model[p,d] == 2)
          parameters <- c(Growth.model[p,d],StartZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
        MeanLengthStartOfYear[p,d,g,a]  <- Growth(age,parameters) 
      }

# If RecAge=0, use linear interpolation for age-0 fish (index a=1)
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
      {
      a <- 1
	  if (RecAge == 0)
        MeanLengthStartOfYear[p,d,g,a]  <-  StartZFraction[p,d]*MeanLengthStartOfYear[p,d,g,(a+1)]
	  else if (RecAge == 1)
	    MeanLengthStartOfYear[p,d,g,a]  <- Growth(a,parameters) 
    }

# Use MaxPlusGroup age adjustment for the plus group by RecAge
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
          tmp <- tmp + exp(-k*NaturalMortality[p,d,g,NAge])*Growth((NAge+k-1),parameters)
        }
	  }
    else if (RecAge == 1) {
	    for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,d] == 1)
            parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],StartZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
          else if (Growth.model[p,d] == 2)
            parameters <- c(Growth.model[p,d],StartZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
		  tmp <- tmp + exp(-k*NaturalMortality[p,d,g,NAge])*Growth((NAge+k),parameters)
      }		  
	  }
      tmp <- tmp*(1.0-exp(-NaturalMortality[p,d,g,NAge]))
      MeanLengthStartOfYear[p,d,g,NAge] <- tmp/(1.0-exp(-NaturalMortality[p,d,g,NAge]*(MaxPlusGroupAge-NAge)))
    }
      
# Compute mean length during spawning season age by gender, area, and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
      for (a in 2:(NAge-1)) {
        if (RecAge == 0)
          age <- a-1
        else if (RecAge == 1)
          age <- a
        if (Growth.model[p,d] == 1)
          parameters <-c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],SpawningZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
        else if (Growth.model[p,d] == 2)
          parameters <-c(Growth.model[p,d],SpawningZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
        MeanLengthSpawning[p,d,g,a]  <-  Growth(age,parameters) 
      }

# If RecAge=0, use linear interpolation for age-0 fish (index a=1)
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
    {
      a <- 1
      if (RecAge == 0)
        MeanLengthSpawning[p,d,g,a]  <-  SpawningZFraction[p,d]*MeanLengthStartOfYear[p,d,g,(a+1)]
      else if (RecAge == 1)
        MeanLengthSpawning[p,d,g,a]  <- Growth(a,parameters) 
    }

# Use MaxPlusGroup age adjustment for the plus group by RecAge
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
          tmp <- tmp + exp(-k*NaturalMortality[p,d,g,NAge])*Growth((NAge+k-1),parameters)
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,d] == 1)
            parameters <- c(Growth.model[p,d],Length.Amin[p,d,g],Length.Amax[p,d,g],SpawningZFraction[p,d],Length.Lmin[p,d,g],Length.Lmax[p,d,g],Length.c[p,d,g])
          else if (Growth.model[p,d] == 2)
            parameters <- c(Growth.model[p,d],SpawningZFraction[p,d],Length.Linf[p,d,g],Length.K[p,d,g],Length.t0[p,d,g])
          tmp <- tmp + exp(-k*NaturalMortality[p,d,g,NAge])*Growth((NAge+k),parameters)
        }		  
      }
      tmp <- tmp*(1.0-exp(-NaturalMortality[p,d,g,NAge]))
      MeanLengthSpawning[p,d,g,NAge] <- tmp/(1.0-exp(-NaturalMortality[p,d,g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# Compute mean length of fishery catch at age by gender, fleet, and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (v in 1:NFleet)
    for (g in 1:NGender)
      for (a in 2:(NAge-1)) {
        if (RecAge == 0)
          age <- a-1
        else if (RecAge == 1)
          age <- a
        if (Growth.model[p,d] == 1)
          parameters <- c(Growth.model[p,FleetArea[v]],Length.Amin[p,FleetArea[v],g],Length.Amax[p,FleetArea[v],g],CatchZFraction[p,FleetArea[v]],Length.Lmin[p,FleetArea[v],g],Length.Lmax[p,FleetArea[v],g],Length.c[p,FleetArea[v],g])
        else if (Growth.model[p,FleetArea[v]] == 2)
          parameters <- c(Growth.model[p,d],CatchZFraction[p,FleetArea[v]],Length.Linf[p,FleetArea[v],g],Length.K[p,FleetArea[v],g],Length.t0[p,FleetArea[v],g])
        MeanLengthCatch[p,FleetArea[v],g,a]  <-  Growth(age,parameters) 
      }

# If RecAge=0, use linear interpolation for age-0 fish (index a=1)
for (p in 1:NPopulation)
  for (v in 1:NFleet)
    for (g in 1:NGender)
    {
      a <- 1
      if (RecAge == 0)
        MeanLengthCatch[p,FleetArea[v],g,a]  <-  CatchZFraction[p,FleetArea[v]]*MeanLengthStartOfYear[p,FleetArea[v],g,(a+1)]
      else if (RecAge == 1)
        MeanLengthCatch[p,FleetArea[v],g,a]  <- Growth(a,parameters) 
    }

# Use MaxPlusGroup age adjustment for the plus group by RecAge
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
          tmp <- tmp + exp(-k*NaturalMortality[p,FleetArea[v],g,NAge])*Growth((NAge+k-1),parameters) 
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,FleetArea[v]] == 1)
            parameters <- c(Growth.model[p,FleetArea[v]],Length.Amin[p,FleetArea[v],g],Length.Amax[p,FleetArea[v],g],CatchZFraction[p,FleetArea[v]],Length.Lmin[p,FleetArea[v],g],Length.Lmax[p,FleetArea[v],g],Length.c[p,FleetArea[v],g])
          else if (Growth.model[p,FleetArea[v]] == 2)
            parameters <- c(Growth.model[p,FleetArea[v]],CatchZFraction[p,FleetArea[v]],Length.Linf[p,FleetArea[v],g],Length.K[p,FleetArea[v],g],Length.t0[p,FleetArea[v],g])
          tmp <- tmp + exp(-k*NaturalMortality[p,FleetArea[v],g,NAge])*Growth((NAge+k),parameters) 
        }		  
      }
      tmp <- tmp*(1.0-exp(-NaturalMortality[p,FleetArea[v],g,NAge]))
      MeanLengthCatch[p,FleetArea[v],g,NAge] <- tmp/(1.0-exp(-NaturalMortality[p,FleetArea[v],g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# Compute mean length of survey catch at age by gender, survey, and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (I in 1:NSurvey)
    for (g in 1:NGender)
      for (a in 2:(NAge-1)) {
        if (RecAge == 0)
          age <- a-1
        else if (RecAge == 1)
          age <- a
        if (Growth.model[p,SurveyArea[I]] == 1)
          parameters <- c(Growth.model[p,SurveyArea[I]],Length.Amin[p,SurveyArea[I],g],Length.Amax[p,SurveyArea[I],g],SurveyZFraction[p,SurveyArea[I]],Length.Lmin[p,SurveyArea[I],g],Length.Lmax[p,SurveyArea[I],g],Length.c[p,SurveyArea[I],g])
        else if (Growth.model[p,SurveyArea[I]] == 2)
          parameters <- c(Growth.model[p,SurveyArea[I]],SurveyZFraction[p,SurveyArea[I]],Length.Linf[p,SurveyArea[I],g],Length.K[p,SurveyArea[I],g],Length.t0[p,SurveyArea[I],g])
        MeanLengthSurvey[p,SurveyArea[I],g,a]  <-  Growth(age,parameters) 
      }

# If RecAge=0, use linear interpolation for age-0 fish (index a=1)
for (p in 1:NPopulation)
  for (I in 1:NSurvey)
    for (g in 1:NGender)
    {
      a <- 1
      if (RecAge == 0)
        MeanLengthSurvey[p,SurveyArea[I],g,a] <-  SurveyZFraction[p,SurveyArea[I]]*MeanLengthStartOfYear[p,SurveyArea[I],g,(a+1)]
      else if (RecAge == 1)
        MeanLengthSurvey[p,SurveyArea[I],g,a] <- Growth(a,parameters)
    }

# Use MaxPlusGroup age adjustment for the plus group by RecAge
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
          tmp <- tmp + exp(-k*NaturalMortality[p,SurveyArea[I],g,NAge])*Growth((NAge+k-1),parameters)
        }
      }
      else if (RecAge == 1) {
        for (k in 0:(MaxPlusGroupAge-NAge-1))  { 
          if (Growth.model[p,SurveyArea[I]] == 1)
            parameters <- c(Growth.model[p,SurveyArea[I]],Length.Amin[p,SurveyArea[I],g],Length.Amax[p,SurveyArea[I],g],SurveyZFraction[p,SurveyArea[I]],Length.Lmin[p,SurveyArea[I],g],Length.Lmax[p,SurveyArea[I],g],Length.c[p,SurveyArea[I],g])
          else if (Growth.model[p,SurveyArea[I]] == 2)
            parameters <- c(Growth.model[p,SurveyArea[I]],SurveyZFraction[p,SurveyArea[I]],Length.Linf[p,SurveyArea[I],g],Length.K[p,SurveyArea[I],g],Length.t0[p,SurveyArea[I],g])
          tmp <- tmp + exp(-k*NaturalMortality[p,SurveyArea[I],g,NAge])*Growth((NAge+k),parameters)
        }		  
      }
      tmp <- tmp*(1.0-exp(-NaturalMortality[p,SurveyArea[I],g,NAge]))
      MeanLengthSurvey[p,SurveyArea[I],g,NAge] <- tmp/(1.0-exp(-NaturalMortality[p,SurveyArea[I],g,NAge]*(MaxPlusGroupAge-NAge)))
    }

# Compute unfished equilibrium mean weights at age at start of year, spawning, catch, and survey
# by gender, area, and population
#-----------------------------------------------------------------------------------------------

for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender) {
      if (WeightAtLength.model[p,d] == 1)
        parameters <- c(WeightAtLength.model[p,d],WeightAtLength.A[p,d,g],WeightAtLength.B[p,d,g])
      for (a in 1:NAge) {
        MeanWeightStartOfYear[p,d,g,a]  <-  WeightAtLength(MeanLengthStartOfYear[p,d,g,a],parameters)
        MeanWeightSpawning[p,d,g,a] <- WeightAtLength(MeanLengthSpawning[p,d,g,a],parameters)
      }
    }
for (p in 1:NPopulation)
  for (v in 1:NFleet)
    for (g in 1:NGender) {
      if (WeightAtLength.model[p,FleetArea[v]] == 1)
        parameters <- c(WeightAtLength.model[p,FleetArea[v]],WeightAtLength.A[p,FleetArea[v],g],WeightAtLength.B[p,FleetArea[v],g])
      for (a in 1:NAge) {
        MeanWeightCatch[p,FleetArea[v],g,a] <-  WeightAtLength(MeanLengthCatch[p,FleetArea[v],g,a],parameters)
      }
    }
for (p in 1:NPopulation)
  for (I in 1:NSurvey)
    for (g in 1:NGender) {
      if (WeightAtLength.model[p,SurveyArea[I]] == 1)
        parameters <- c(WeightAtLength.model[p,SurveyArea[I]],WeightAtLength.A[p,SurveyArea[I],g],WeightAtLength.B[p,SurveyArea[I],g])
      for (a in 1:NAge)
        MeanWeightSurvey[p,SurveyArea[I],g,a] <-  WeightAtLength(MeanLengthSurvey[p,SurveyArea[I],g,a],parameters)
    }
      
