################################################################################################
# Calculate_Assessment_Period.R
################################################################################################
# The R program file Calculate_Assessment_Period.R calculates
# the fishery system conditions during the assessment time horizon 
# for MAS_Data_Simulator_v1 of the MAS project.
# Jon Brodziak, PIFSC, jon.brodziak@noaa.gov 18-Jun-2020
# Updated to include RecAge 21-Apr-2021
################################################################################################

# (2.3.1) SET INITIAL POPULATION NUMBERS AT AGE BY POPULATION, AREA, AND GENDER ON JANUARY 1ST
#-----------------------------------------------------------------------------------------------
NumbersAtAge[1,,,,] <- FishedEquilibriumNumbersAtAge

################################################################################################
# LOOP OVER ASSESSMENT TIME HORIZON TO COMPUTE NUMBERS AT AGE
# BY POPULATION, AREA, AND GENDER, AND OTHER QOI
#-----------------------------------------------------------------------------------------------
print('_____________________________________________________________________________________________________')
print('Calculating population numbers at age and spawning biomasses by area through assessment time horizon ...')
print('_____________________________________________________________________________________________________')
for (y in 1:NYear)
{
  
  # Store mean lengths and weights for each year in the assessment time horizon
  # If year y=1 use the fished equilibrium sizes at age
  # Else y>1 use the previous year (y-1) values
  #-----------------------------------------------------------------------------------------------
  if (y == 1)
  {
    AssessmentMeanLengthStartOfYear[y,,,,] <- FishedMeanLengthStartOfYear
    AssessmentMeanLengthSpawning[y,,,,] <- FishedMeanLengthSpawning
    AssessmentMeanLengthCatch[y,,,,] <- FishedMeanLengthCatch
    AssessmentMeanLengthSurvey[y,,,,] <- FishedMeanLengthSurvey
    AssessmentMeanWeightStartOfYear[y,,,,] <- FishedMeanWeightStartOfYear 
    AssessmentMeanWeightSpawning[y,,,,] <- FishedMeanWeightSpawning
    AssessmentMeanWeightCatch[y,,,,] <- FishedMeanWeightCatch
    AssessmentMeanWeightSurvey[y,,,,] <- FishedMeanWeightSurvey
  }
  else
  {
    AssessmentMeanLengthStartOfYear[y,,,,] <- AssessmentMeanLengthStartOfYear[(y-1),,,,]
    AssessmentMeanLengthSpawning[y,,,,] <- AssessmentMeanLengthSpawning[(y-1),,,,]
    AssessmentMeanLengthCatch[y,,,,] <- AssessmentMeanLengthCatch[(y-1),,,,]
    AssessmentMeanLengthSurvey[y,,,,] <- AssessmentMeanLengthSurvey[(y-1),,,,]
    AssessmentMeanWeightStartOfYear[y,,,,] <- AssessmentMeanWeightStartOfYear[(y-1),,,,] 
    AssessmentMeanWeightSpawning[y,,,,] <- AssessmentMeanWeightSpawning[(y-1),,,,]
    AssessmentMeanWeightCatch[y,,,,] <- AssessmentMeanWeightCatch[(y-1),,,,]
    AssessmentMeanWeightSurvey[y,,,,] <- AssessmentMeanWeightSurvey[(y-1),,,,]
  }
  
  # (2.3.2) CALCULATE FISHERY SELECTIVITY FOR YEAR Y
  #-----------------------------------------------------------------------------------------------
  FisherySelectivityAtAge <- EquilibriumFisherySelectivityAtAge
  
  # (2.3.3) CALCULATE SURVEY SELECTIVITY FOR YEAR Y
  #-----------------------------------------------------------------------------------------------
  SurveySelectivityAtAge <- SurveySelectivityAtAge
  
  # (2.3.4) CALCULATE TOTAL MORTALITY FOR YEAR y
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (d in 1:NArea)
      for (g in 1:NGender)
        for (a in 1:NAge)
        {          
          TotalMortality[y,p,d,g,a] <- NaturalMortality[p,d,g,a]
          for (v in 1:NFleet)
            TotalMortality[y,p,FleetArea[v],g,a] <- TotalMortality[y,p,FleetArea[v],g,a]+AssessmentFishingMortality[y,FleetArea[v],d]*FisherySelectivityAtAge[p,FleetArea[v],g,a]
        }
  
  # (2.3.5) CALCULATE SPAWNING NUMBERS AT AGE AND SPAWNING BIOMASS FOR YEAR y
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (d in 1:NArea)
      for (g in 1:NGender)
        for (a in 1:NAge)
          SpawningNumbersAtAge[p,d,g,a] <- exp(-SpawningZFraction[p,d]*TotalMortality[y,p,d,g,a])*NumbersAtAge[y,p,d,g,a]
  
  for (p in 1:NPopulation)
    for (d in 1:NArea)
      for (g in 1:NGender)
      {
        tmp <- 0.0
        for (a in 1:NAge)
        {
          tmp <- tmp + SpawningNumbersAtAge[p,d,g,a]*AssessmentMeanWeightSpawning[y,p,d,g,a]*MaturityProbabilityAtAge[p,d,g,a]
        }
        SpawningBiomass[y,p,d,g] <- tmp
        
        # Rescale to Spawning Biomass Units
        SpawningBiomass[y,p,d,g] <- SpawningBiomass[y,p,d,g]/Recruitment.SpawningBiomassUnits[p,d] 
      }
  
  # (2.3.6) CALCULATE RECRUITMENT BY POPULATION, AREA, AND GENDER FOR YEAR y
  # ----------------------------------------------------------------------------------------------
  # Compute recruitment production by population and area for year=y
  #-----------------------------------------------------------------
  if (y == 1)
    for (p in 1:NPopulation)
      for (d in 1:NArea)
      {
        parameters <- c(Recruitment.model[p,d],UnfishedSpawningBiomass[p,d,1], Recruitment.UnfishedR[p,d], Recruitment.steepness[p,d])
        RecruitmentProduction[y,p,d] <- StockRecruitment(FishedEquilibriumSpawningBiomass[p,d,1],parameters)
      }
  else if (y > 1)
    for (p in 1:NPopulation)
      for (d in 1:NArea)
      {
        parameters <- c(Recruitment.model[p,d],UnfishedSpawningBiomass[p,d,1], Recruitment.UnfishedR[p,d], Recruitment.steepness[p,d])
		if (RecAge == 0)
          RecruitmentProduction[y,p,d] <- StockRecruitment(SpawningBiomass[(y),p,d,1],parameters)
		else if (RecAge == 1)
		  RecruitmentProduction[y,p,d] <- StockRecruitment(SpawningBiomass[(y-1),p,d,1],parameters)
      }
  
  # Compute annual recruitment strength by population, 
  # area, and gender for year=y and age-0 index a=1
  #-----------------------------------------------------------------
  a <- 1
  for (p in 1:NPopulation) {
    for (d in 1:NArea)
      for (g in 1:NGender)
      {   
        tmp <- 0.0
        for (k in 1:NArea)
        {
          tmp <- tmp + Recruitment.GenderFraction[p,k,g]*RecruitmentDistribution[p,k,d]*RecruitmentProduction[y,p,k]
        }
        Recruitment[y,p,d,g] <- tmp
      }
  }
  
  # Update NumbersAtAge at age index (a=1) with predicted recruitment
  # by population, area, and gender for year=y
  #-----------------------------------------------------------------
  for (p in 1:NPopulation)
    for (d in 1:NArea)
      for (g in 1:NGender)
        for (a in 1:NAge)
          NumbersAtAge[y,p,d,g,1] <- Recruitment[y,p,d,g]
  
  # (2.3.7) CALCULATE FISHERY OBSERVATIONS
  # CALCULATE FISHERY CATCH NUMBERS AT AGE FOR YEAR y
  # BY POPULATION, FLEET, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (v in 1:NFleet)
      for (d in 1:NArea)
        if (FleetArea[v]==d)
          for (g in 1:NGender)
            for (a in 1:NAge)
            {
              FisheryCatchNumbersAtAgeByPopulation[y,p,FleetArea[v],d,g,a] <- NumbersAtAge[y,p,d,g,a]*AssessmentFishingMortality[y,FleetArea[v],d]*FisherySelectivityAtAge[p,FleetArea[v],g,a]
              FisheryCatchNumbersAtAgeByPopulation[y,p,v,d,g,a] <-  FisheryCatchNumbersAtAgeByPopulation[y,p,FleetArea[v],d,g,a]*(1.0-exp(-TotalMortality[y,p,d,g,a]))/TotalMortality[y,p,d,g,a]
            }
  
  # CALCULATE FISHERY CATCH PROPORTIONS AT AGE FOR YEAR y
  # BY POPULATION, FLEET, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (v in 1:NFleet)
      for (d in 1:NArea)
        if (FleetArea[v]==d)
          for (g in 1:NGender)
          {
            tmp <- sum(FisheryCatchNumbersAtAgeByPopulation[y,p,FleetArea[v],d,g,])
            for (a in 1:NAge)
            {
              FisheryCatchProportionAtAgeByPopulation[y,p,FleetArea[v],d,g,a] <- FisheryCatchNumbersAtAgeByPopulation[y,p,FleetArea[v],d,g,a]/tmp
            }
          }
  
  # CALCULATE FISHERY CATCH BIOMASS AT AGE FOR YEAR y
  # BY POPULATION, FLEET, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (v in 1:NFleet)
      for (d in 1:NArea)
        if (FleetArea[v]==d)
          for (g in 1:NGender)
            for (a in 1:NAge)
            {
              tmp <- 0.0
              FisheryCatchBiomassAtAgeByPopulation[y,p,FleetArea[v],d,g,a] <- FisheryCatchNumbersAtAgeByPopulation[y,p,FleetArea[v],d,g,a]*AssessmentMeanWeightCatch[y,p,FleetArea[v],g,a]
            }
  
  # (2.3.8) CALCULATE SURVEY OBSERVATIONS
  # CALCULATE SURVEY CATCH NUMBERS AT AGE
  # BY YEAR, POPULATION, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (I in 1:NSurvey)
      for (g in 1:NGender)
        for (a in 1:NAge)
          SurveyCatchNumbersAtAgeByPopulation[y,p,SurveyArea[I],g,a] <- SurveyCatchability[p,SurveyArea[I]]*SurveySelectivityAtAge[p,SurveyArea[I],g,a]*exp(-SurveyZFraction[p,SurveyArea[I]]*TotalMortality[y,p,SurveyArea[I],g,a])*NumbersAtAge[y,p,SurveyArea[I],g,a]
  
  # CALCULATE SURVEY CATCH PROPORTIONS AT AGE FOR YEAR y
  # BY AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (I in 1:NSurvey)
      for (g in 1:NGender)
      {
        tmp <- sum(SurveyCatchNumbersAtAgeByPopulation[y,p,SurveyArea[I],g,])
        for (a in 1:NAge)
        {
          SurveyCatchProportionAtAgeByPopulation[y,p,SurveyArea[I],g,a] <- SurveyCatchNumbersAtAgeByPopulation[y,p,SurveyArea[I],g,a]/tmp
        }
      }
  
  # CALCULATE SURVEY CATCH BIOMASS AT AGE FOR YEAR y
  # BY POPULATION, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (p in 1:NPopulation)
    for (I in 1:NSurvey)
      for (g in 1:NGender)
        for (a in 1:NAge)
        {
          tmp <- 0.0
          SurveyCatchBiomassAtAgeByPopulation[y,p,SurveyArea[I],g,a] <- SurveyCatchNumbersAtAgeByPopulation[y,p,SurveyArea[I],g,a]*AssessmentMeanWeightSurvey[y,p,SurveyArea[I],g,a]
        }
  
  # (2.3.9) CALCULATE QUANTITIES OF INTEREST
  # CALCULATE FISHERY CATCH NUMBERS AT AGE FOR YEAR y
  # BY FLEET, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (v in 1:NFleet)
    for (d in 1:NArea)
      if (FleetArea[v]==d)
        for (g in 1:NGender)
          for (a in 1:NAge)
          {
            FisheryCatchNumbersAtAge[y,FleetArea[v],d,g,a] <-  sum(FisheryCatchNumbersAtAgeByPopulation[y,,FleetArea[v],d,g,a])
          }
  
  # CALCULATE FISHERY CATCH PROPORTIONS AT AGE FOR YEAR y
  # BY FLEET, AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (v in 1:NFleet)
    for (d in 1:NArea)
      if (FleetArea[v]==d)
        for (g in 1:NGender)
        {
          tmp <- sum(FisheryCatchNumbersAtAge[y,FleetArea[v],d,g,])
          for (a in 1:NAge)
          {
            FisheryCatchProportionAtAge[y,FleetArea[v],d,g,a] <- FisheryCatchNumbersAtAge[y,FleetArea[v],d,g,a]/tmp
          }
        }
  
  # CALCULATE FISHERY CATCH BIOMASS FOR YEAR y
  # BY FLEET AND AREA
  # ----------------------------------------------------------------------------------------------
  for (v in 1:NFleet)
    for (d in 1:NArea)
      if (FleetArea[v]==d)
        FisheryCatchBiomass[y,FleetArea[v],d] <- sum(FisheryCatchBiomassAtAgeByPopulation[y,,FleetArea[v],d,,])
  
  # CALCULATE SURVEY CATCH NUMBERS AT AGE FOR YEAR y
  # BY AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (I in 1:NSurvey)
    for (g in 1:NGender)
      for (a in 1:NAge)
      {
        SurveyCatchNumbersAtAge[y,SurveyArea[I],g,a] <-  sum(SurveyCatchNumbersAtAgeByPopulation[y,,SurveyArea[I],g,a])
      }
  
  # CALCULATE SURVEY CATCH PROPORTIONS AT AGE FOR YEAR y
  # BY AREA, AND GENDER
  # ----------------------------------------------------------------------------------------------
  for (I in 1:NSurvey)
    for (g in 1:NGender)
    {
      tmp <- sum(SurveyCatchNumbersAtAge[y,SurveyArea[I],g,])
      for (a in 1:NAge)
      {
        SurveyCatchProportionAtAge[y,SurveyArea[I],g,a] <- SurveyCatchNumbersAtAge[y,SurveyArea[I],g,a]/tmp
      }
    }
  
  # CALCULATE SURVEY CATCH BIOMASS FOR YEAR y
  # BY AREA
  # ----------------------------------------------------------------------------------------------
  for (I in 1:NSurvey)
    SurveyCatchBiomass[y,SurveyArea[I]] <- sum(SurveyCatchBiomassAtAgeByPopulation[y,,SurveyArea[I],,])
  
  # (2.3.10) CALCULATE NUMBERS AT AGE IN YEAR Y+1
  # Compute numbers at age by population, area, and gender
  # for true ages [a=2:(NAge-1)]
  #-----------------------------------------------------------------
  for (p in 1:NPopulation) {
    for (d in 1:NArea)
      for (g in 1:NGender)
        for (a in 2:(NAge-1))
        {
          tmp <- 0.0
          for (k in 1:NArea)
          {
            tmp <- tmp + MovementProbability[p,k,d]*NumbersAtAge[y,p,k,g,(a-1)]*exp(-TotalMortality[y,p,k,g,(a-1)])
          }
          NumbersAtAge[(y+1),p,d,g,a] <- tmp
        }
  }
  
  # Compute fished numbers at age by population, area, and gender
  # For the plus-group age  [a=NAge]
  #-----------------------------------------------------------------
  a <- NAge
  
  for (p in 1:NPopulation) 
  {
    for (d in 1:NArea)
      for (g in 1:NGender)
      { 
        tmp <- 0.0
        for (k in 1:NArea)
        {
          # Age-(A-1) survivors that remained or immigrated to area d
          tmp <- tmp + MovementProbability[p,k,d]*NumbersAtAge[y,p,k,g,(a-1)]*exp(-TotalMortality[y,p,k,g,(a-1)])
          # Age-A survivors that remained or immigrated to area d
          tmp <- tmp + MovementProbability[p,k,d]*NumbersAtAge[y,p,k,g,a]*exp(-TotalMortality[y,p,k,g,a])
        }
        NumbersAtAge[(y+1),p,d,g,a] <- tmp
      }
  }
 
  # Update mean lengths and weights of the plus group, age index=NAge
  #-------------------------------------------------------------------
  if (y > 1) 
  {
    # Calculate time-varying mean lengths at age ONLY
    # for the plus group using the MaxPlusGroupAge adjustment
	# tmp1 is weighting for mean lengths of age (NAge-1) survivors
	# tmp2 is weighting for mean lengths of age NAge+ survivors
	# tmp3 is the weighted mean length of the true age class (NAge-1)
	# using weight tmp1/(tmp1+tmp2) plus the weighted mean length
	# of the plus group NAge+ using weight tmp2/(tmp1+tmp2)
    # This calculation assumes that growth is contant in time
    # and that the plus group mean length changes with time-varying F
    # As a result one can use the unfished mean lengths for true age (NAge-1)	
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender) 
        {
          tmp1 <- NumbersAtAge[(y-1),p,d,g,(NAge-1)]*exp(-TotalMortality[(y-1),p,d,g,(NAge-1)])
          tmp2 <- NumbersAtAge[(y-1),p,d,g,NAge]*exp(-TotalMortality[(y-1),p,d,g,NAge])
		  
          tmp3 <- tmp1*FishedMeanLengthStartOfYear[p,d,g,NAge]+tmp2*AssessmentMeanLengthStartOfYear[(y-1),p,d,g,NAge]
          AssessmentMeanLengthStartOfYear[y,p,d,g,NAge] <- tmp3/(tmp1+tmp2)
		  
          tmp3 <- tmp1*FishedMeanLengthSpawning[p,d,g,NAge]+tmp2*AssessmentMeanLengthSpawning[(y-1),p,d,g,NAge]
          AssessmentMeanLengthSpawning[y,p,d,g,NAge] <- tmp3/(tmp1+tmp2)
		  
          tmp3 <- tmp1*FishedMeanLengthCatch[p,d,g,NAge]+tmp2*AssessmentMeanLengthCatch[(y-1),p,d,g,NAge]
          AssessmentMeanLengthCatch[y,p,d,g,NAge] <- tmp3/(tmp1+tmp2)
		  
          tmp3 <- tmp1*FishedMeanLengthSurvey[p,d,g,NAge]+tmp2*AssessmentMeanLengthSurvey[(y-1),p,d,g,NAge]
          AssessmentMeanLengthSurvey[y,p,d,g,NAge] <- tmp3/(tmp1+tmp2) 
       }

    # Calculate time-varying mean weights at age 
    # for the plus group using the MaxPlusGroupAge adjustment
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender) {
		  if (WeightAtLength.model[p,d] == 1)
		    parameters <- c(WeightAtLength.model[p,d],WeightAtLength.A[p,d,g],WeightAtLength.B[p,d,g])
          AssessmentMeanWeightStartOfYear[y,p,d,g,NAge]  <-  WeightAtLength(AssessmentMeanLengthStartOfYear[y,p,d,g,NAge],parameters)
          AssessmentMeanWeightSpawning[y,p,d,g,NAge] <- WeightAtLength(AssessmentMeanLengthSpawning[y,p,d,g,NAge],parameters)
        }
    for (p in 1:NPopulation)
      for (v in 1:NFleet)
        for (g in 1:NGender) {
		  if (WeightAtLength.model[p,d] == 1)
		    parameters <- c(WeightAtLength.model[p,d],WeightAtLength.A[p,FleetArea[v],g],WeightAtLength.B[p,FleetArea[v],g])
          AssessmentMeanWeightCatch[y,p,FleetArea[v],g,NAge] <-  WeightAtLength(AssessmentMeanLengthCatch[y,p,FleetArea[v],g,NAge],parameters)
        }
    for (p in 1:NPopulation)
      for (I in 1:NSurvey)
        for (g in 1:NGender) {
		  if (WeightAtLength.model[p,d] == 1)
		    parameters <- c(WeightAtLength.model[p,d],WeightAtLength.A[p,SurveyArea[I],g],WeightAtLength.B[p,SurveyArea[I],g])
          AssessmentMeanWeightSurvey[y,p,SurveyArea[I],g,NAge] <-  WeightAtLength(AssessmentMeanLengthSurvey[y,p,SurveyArea[I],g,NAge],parameters)
		}
  }

}
# END OF FOR LOOP OVER TIME

