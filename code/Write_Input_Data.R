################################################################################################
# Write_Input_Data.R
################################################################################################
# The R program file Write_Input_Data.R writes the input data and
# results of initialization calculations for case study v1 
# of the MAS project to user-named OutputFile
# Jon Brodziak, PIFSC, jon.brodziak@noaa.gov 3-DEC-2020
################################################################################################
Write_Input_Data<-function(filename)
{
  print('_____________________________________________________________________________________________________')
  print('Writing input data to output file ...')
  print('_____________________________________________________________________________________________________')

  sink(filename)

  print('#################################################################################################')
  print('# INITIIALIZING ANALYSIS CLASS')
  print('#################################################################################################')
  print('#')
  
  print('#################################################################################################')
  print('# SETTING MODEL DOMAIN PARAMETERS')
  print('#################################################################################################')
  
  print('# Number of populations (1:NPopulation)')
  print(NPopulation)
  
  print('# Number of age classes (1:NAge)')
  print(NAge)
  
  print('# Recruitment age or first age in the population numbers at age vector')
  print(RecAge)
  
  print('# Maximum plus group age')
  print(MaxPlusGroupAge)
  
  print('# Number of genders (1:NGender)')
  print(NGender)
  
  print('# Number of areas (1:NArea)')
  print(NArea)
  
  print('# Number of fleets (1:NFleet)')
  print(NFleet)
  
  for (I in 1:NFleet)
  print(c("# Fleet =",I,"operates in Area=",FleetArea[I]))
  
  print('# Number of surveys (1:NSurvey)')
  print(NSurvey)
  
  for (I in 1:NSurvey)
  print(c("# Survey=",I,"operates in Area=",SurveyArea[I]))
  
  print('# Number of years in assessment time horizon (1:NYear)')
  print(NYear)
  
  print('# Number of seasons in a year (1:NSeason)')
  print(NSeason)
  
  print('# Maximum number of iterations for computing initial equilibrium (1:MaxIteration)')
  print(MaxIteration)
  
  print('# Convergence criterion for computing initial equilibrium: EquilibriumConvergenceCriterion')
  print(EquilibriumConvergenceCriterion)
  
  print('# Population numbers units for output (numbers of fish)')
  print(OutputNumbersUnits)
  
  print('# Population biomass units for output (kilograms of fish)')
  print(OutputBiomassUnits)
  
  print('# Maximum fishing mortality for reference point calculations')
  print(MaxF)
  
  print('# Fishing mortality mesh for reference point calculations')
  print(FMesh)
  
  print('# Number of grid points for reference point calculations')
  print(NGrid)
  
  print('# Fishing mortality grid for reference point calculations')
  print(FGrid)
  
  print('#################################################################################################')
  print('# INITIIALIZING POPULATION CLASS')
  print('#################################################################################################')
  
  print('# Timing of start of year within year: StartZFraction[NPopulation,NArea]')
  print(StartZFraction)
  
  print('# Timing of spawning within year: SpawningZFraction[NPopulation,NArea]')
  print(SpawningZFraction)
  
  print('# Timing of catch within year: CatchZFraction[NPopulation,NFleet]')
  print(CatchZFraction)
  
  print('# Timing of survey within year: SurveyZFraction[NPopulation,NSurvey]')
  print(SurveyZFraction)
  
  print('# Population movement probability arrays: MovementProbability[NPopulation,NArea,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Movement by area array for population:',i))
    print(MovementProbability[i,,])
  }
  
  print('# Population recruitment distribution arrays: RecruitmentDistribution[NPopulation,NArea,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Recruitment distribution by area array for population:',i))
    print(RecruitmentDistribution[i,,])
  }
  
  print('# Population natural mortality arrays: NaturalMortality[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Natural mortality at age array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(NaturalMortality[i,j,k,])
  }
  
  print('# Population maturity models: Maturity.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Maturity model array for population:',i))
    print(Maturity.model[i,])
  }
  
  if (Maturity.model[1,1] == 1)
	{
	print('# Population maturity parameters: Maturity.a50[NPopulation,NArea,NGender]')
	for (i in 1:NPopulation) {
	  print(c('# Maturity a50 parameter array for population:',i))
	  for (j in 1:NArea)
		print(Maturity.a50[i,j,])
	  }
  
	print('# Population maturity parameters: Maturity.slope[NPopulation,NArea,NGender]')
	for (i in 1:NPopulation) {
	  print(c('# Maturity slope parameter array for population:',i))
	  for (j in 1:NArea)
		print(Maturity.slope[i,j,])
	  }
	}
  
  print('# Population maturity probability at age arrays: MaturityProbabilityAtAge[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Maturity proabability at age array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) {
        print(c("Area:",j,"Gender:",k))
		print(MaturityProbabilityAtAge[i,j,k,]) 
      }		
  }
    
  print('# Population growth models: Growth.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Growth model array for population:',i))
    print(Growth.model[i,])
  }

  if (Growth.model[1,1] == 1)
    {
# Growth.model == 1 is the modified von Bertalanffy growth curve
	  print('# Population growth parameters: Minimum reference age Length.Amin[NPopulation,NArea,NGender]:')
	  for (i in 1:NPopulation) {
		print(c('# Age at amin (yr) array for population:',i))
		for (j in 1:NArea)
		  print(Length.Amin[i,j,])
	  }
	  
	  print('# Population growth parameters: Maximum reference age Length.Amax[NPopulation,NArea,NGender]:')
	  for (i in 1:NPopulation) {
		print(c('# Age at amax (yr) array for population:',i))
		for (j in 1:NArea)
		  print(Length.Amax[i,j,])
	  }
	  
	  print('# Population growth parameters: Length.Lmin[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length at age amin (cm) array for population:',i))
		for (j in 1:NArea)
		  print(Length.Lmin[i,j,])
	  }
	  
	  print('# Population growth parameters: Length.Lmax[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length at age amax (cm) array for population:',i))
		for (j in 1:NArea)
		  print(Length.Lmax[i,j,])
	  }
	  
	  print('# Population growth parameters: Length.c[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length at age curvature parameter array for population:',i))
		for (j in 1:NArea)
		  print(Length.c[i,j,])
	  }
  }
  else if (Growth.model[1,1] == 2) {
# Growth.model == 2 is the standard von Bertalanffy growth curve
	  print('# Population growth parameters: Length.Linf[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length at age Linf (cm) array for population:',i))
		for (j in 1:NArea)
		  print(Length.Linf[i,j,])
	  }
	  
	  print('# Population growth parameters: Length.K[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length at age K array for population:',i))
		for (j in 1:NArea)
		  print(Length.K[i,j,])
	  }

	  print('# Population growth parameters: Length.t0[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length at age t0 array for population:',i))
		for (j in 1:NArea)
		  print(Length.t0[i,j,])
	  }	  
  } 
  
  print('# Unfished population mean length at age on January 1st arrays: UnfishedMeanLengthStartOfYear[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean length at age (cm) on January 1st array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(UnfishedMeanLengthStartOfYear[i,j,k,]) 
  }
  
  print('# Unfished population mean length at age during spawning season arrays: UnfishedMeanLengthSpawning[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean length at age (cm) during spawning season array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(UnfishedMeanLengthSpawning[i,j,k,]) 
  }
  
  print('# Unfished population mean length at age of fishery catch arrays: UnfishedMeanLengthCatch[NPopulation,Nfleet,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean length at age (cm) of fishery catch array for population:',i))
    for (j in 1:NFleet)
      for (k in 1:NGender) 
        print(UnfishedMeanLengthCatch[i,j,k,]) 
  }
  
  print('# Unfished population mean length at age during survey arrays: UnfishedMeanLengthSurvey[NPopulation,NSurvey,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean length at age (cm) during survey array for population:',i))
    for (j in 1:NSurvey)
      for (k in 1:NGender) 
        print(UnfishedMeanLengthSurvey[i,j,k,]) 
  }

  print('# Population weight at length models: WeightAtLength.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Weight at length model array for population:',i))
    print(WeightAtLength.model[i,])
  }

  if (WeightAtLength.model[1,1] == 1) {
	  print('# Population length-weight parameters: WeightAtLength.A[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length-weight scalar parameter A array for population:',i))
		for (j in 1:NArea)
		  print(WeightAtLength.A[i,j,])
	  }
	  
	  print('# Population length-weight parameters: WeightAtLength.B[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length-weight exponent parameter B array for population:',i))
		for (j in 1:NArea)
		  print(WeightAtLength.B[i,j,])
	  }
  }
  
  print('# Unfished population mean weight at age on January 1st arrays: UnfishedMeanWeightStartOfYear[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean weight at age (kg) on January 1st array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(UnfishedMeanWeightStartOfYear[i,j,k,]) 
  }
  
  print('# Unfished population mean weight at age during spawning season arrays: UnfishedMeanWeightSpawning[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean weight at age (kg) during spawning season array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(UnfishedMeanWeightSpawning[i,j,k,]) 
  }
  
  print('# Unfished population mean weight at age of fishery catch arrays: UnfishedMeanWeightCatch[NPopulation,NFleet,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean weight at age (kg) of fishery catch array for population:',i))
    for (j in 1:NFleet)
      for (k in 1:NGender) 
        print(UnfishedMeanWeightCatch[i,j,k,]) 
  }
  
  print('# Unfished population mean weight at age during survey arrays: UnfishedMeanWeightSurvey[NPopulation,NSurvey,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Unfished mean weight at age (kg) during survey array for population:',i))
    for (j in 1:NSurvey)
      for (k in 1:NGender) 
        print(UnfishedMeanWeightSurvey[i,j,k,]) 
  }

  print('# Fished equilibrium population mean length at age on January 1st arrays: FishedMeanLengthStartOfYear[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean length at age (cm) on January 1st array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(FishedMeanLengthStartOfYear[i,j,k,]) 
  }
  
  print('# Fished equilibrium population mean length at age during spawning season arrays: FishedMeanLengthSpawning[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean length at age (cm) during spawning season array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(FishedMeanLengthSpawning[i,j,k,]) 
  }
  
  print('# Fished equilibrium population mean length at age of fishery catch arrays: FishedMeanLengthCatch[NPopulation,Nfleet,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean length at age (cm) of fishery catch array for population:',i))
    for (j in 1:NFleet)
      for (k in 1:NGender) 
        print(FishedMeanLengthCatch[i,j,k,]) 
  }
  
  print('# Fished equilibrium population mean length at age during survey arrays: FishedMeanLengthSurvey[NPopulation,NSurvey,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean length at age (cm) during survey array for population:',i))
    for (j in 1:NSurvey)
      for (k in 1:NGender) 
        print(FishedMeanLengthSurvey[i,j,k,]) 
  }

  print('# Population weight at length models: WeightAtLength.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Weight at length model array for population:',i))
    print(WeightAtLength.model[i,])
  }

  if (WeightAtLength.model[1,1] == 1) {
	  print('# Population length-weight parameters: WeightAtLength.A[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length-weight scalar parameter A array for population:',i))
		for (j in 1:NArea)
		  print(WeightAtLength.A[i,j,])
	  }
	  
	  print('# Population length-weight parameters: WeightAtLength.B[NPopulation,NArea,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Length-weight exponent parameter B array for population:',i))
		for (j in 1:NArea)
		  print(WeightAtLength.B[i,j,])
	  }
  }
  
  print('# Fished equilibrium population mean weight at age on January 1st arrays: FishedMeanWeightStartOfYear[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean weight at age (kg) on January 1st array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(FishedMeanWeightStartOfYear[i,j,k,]) 
  }
  
  print('# Fished equilibrium population mean weight at age during spawning season arrays: FishedMeanWeightSpawning[NPopulation,NArea,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean weight at age (kg) during spawning season array for population:',i))
    for (j in 1:NArea)
      for (k in 1:NGender) 
        print(FishedMeanWeightSpawning[i,j,k,]) 
  }
  
  print('# Fished equilibrium population mean weight at age of fishery catch arrays: FishedMeanWeightCatch[NPopulation,NFleet,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean weight at age (kg) of fishery catch array for population:',i))
    for (j in 1:NFleet)
      for (k in 1:NGender) 
        print(FishedMeanWeightCatch[i,j,k,]) 
  }
  
  print('# Fished equilibrium population mean weight at age during survey arrays: FishedMeanWeightSurvey[NPopulation,NSurvey,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Fished mean weight at age (kg) during survey array for population:',i))
    for (j in 1:NSurvey)
      for (k in 1:NGender) 
        print(FishedMeanWeightSurvey[i,j,k,]) 
  }

  
  print('# Units of female spawning biomass (kg) for recruitment process models: Recruitment.SpawningBiomassUnits[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Units of female spawning biomass (kg) by area array for population:',i))
    print(Recruitment.SpawningBiomassUnits[i,])
  }

  print('# Fraction of recruits by gender parameters for recruitment process models: Recruitment.GenderFraction[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Fraction of recruits by gender array for population:',i))
    for (j in 1:NArea)
      print(Recruitment.GenderFraction[i,j,])
  }

  print('# Population reruitment models: Recruitment.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Recruitment model array for population:',i))
    print(Recruitment.model[i,])
  }

  if (Recruitment.model[1,1] == 1) {  
	  print('# Steepness parameter for recruitment process models: Recruitment.steepness[NPopulation,NArea]')
	  for (i in 1:NPopulation) {
		print(c('# Steepness parameter by area array for population:',i))
		print(Recruitment.steepness[i,])
	  }
	  
	  print('# Unfished recruitment parameter for recruitment process models: Recruitment.UnfishedR[NPopulation,NArea]')
	  for (i in 1:NPopulation) {
	  print(c('# Unfished recruitment parameter by area array for population:',i))
	  print(Recruitment.UnfishedR[i,])
	  }
	  
	  print('# Logarithm of unfished recruitment parameter for recruitment process models: Recruitment.LogUnfishedR[NPopulation,NArea]')
	  for (i in 1:NPopulation) {
		print(c('# Logarithm of unfished recruitment parameter by area array for population:',i))
		print(Recruitment.LogUnfishedR[i,])
	  }
	  
	  print('# Standard deviation of log-scale R for recruitment process models: Recruitment.sigma.R[NPopulation,NArea]')
	  for (i in 1:NPopulation) {
		print(c('# Standard deviation of log-scale R parameter by area array for population:',i))
		print(Recruitment.sigmaR[i,])
	  }
  }
  print('#')
  
  print('#################################################################################################')
  print('# INITIIALIZING OBSERVATION CLASS')
  print('#################################################################################################')
  
  print('# Equilibrium fishing mortality by fleet: SimEquilibriumFishingMortality[NFleet,NArea]')
  for (i in 1:NFleet) {
    print(c('# Equilibrium fishing mortality by area array for fleet:',i))
    print(SimEquilibriumFishingMortality[i,])
  }  
  
  print('# Fishing mortality by year and fleet: SimFishingMortality[NYear,NFleet,NArea]')
  for (i in 1:NFleet) {
    print(c('# Fishing mortality by year and area arrays for fleet:',i))
    for (j in 1:NArea)
      print(SimFishingMortality[,i,j])
  }   
    
  print('# Population fishery selectivity models: FisherySelectivity.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Fishery selectivity model array for population:',i))
    print(FisherySelectivity.model[i,])
  }

  if (FisherySelectivity.model[p,d] == 1) {  
	  print('# Fishery selectivity parameters for populations by fleet/area: FisherySelectivity.a50[NPopulation,NFleet,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Fishery age at 50% selectivity by fleet/area and gender arrays for population:',i))
		for (j in 1:NFleet)
		  print(FisherySelectivity.a50[i,j,])  
	  }
	  
	  print('# Fishery selectivity parameters for populations by fleet/area: FisherySelectivity.slope[NPopulation,Nfleet,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Fishery selectivity slope parameter by area and gender arrays for population:',i))
		print(FisherySelectivity.slope[i,,])  
	  }	
  }
  
  print('# Population fishery selectivity at age arrays: FisherySelectivityAtAge[NPopulation,NFleet,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Initial fishery selectivity at age by area and gender arrays for population:',i))
    for (j in 1:NFleet)
      for (k in 1:NGender) 
        print(FisherySelectivityAtAge[i,j,k,])  
  }
     
  print('# Population survey selectivity models: SurveySelectivity.model[NPopulation,NArea]')
  for (i in 1:NPopulation) {
    print(c('# Survey selectivity model array for population:',i))
    print(SurveySelectivity.model[i,])
  } 
  
  if (SurveySelectivity.model[1,1] == 1) {
	  print('# Survey selectivity parameters for populations by fleet/area: SurveySelectivity.a50[NPopulation,NSurvey,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Survey age at 50% selectivity by fleet/area and gender arrays for population:',i))
		for (j in 1:NSurvey)
		  print(SurveySelectivity.a50[i,j,]) 
	  }
	  
	  print('# Survey selectivity parameters for populations by fleet/area: SurveySelectivity.slope[NPopulation,NSurvey,NGender]')
	  for (i in 1:NPopulation) {
		print(c('# Survey selectivity slope parameter by area and gender arrays for population:',i))
		print(SurveySelectivity.slope[i,,])  
	  }
  }
  
  print('# Population survey selectivity at age arrays: SurveySelectivityAtAge[NPopulation,NSurvey,NGender,NAge]')
  for (i in 1:NPopulation) {
    print(c('# Initial survey selectivity at age by area and gender arrays for population:',i))
    for (j in 1:NSurvey)
      for (k in 1:NGender) 
        print(SurveySelectivityAtAge[i,j,k,])  
  }
  print('#')
  
  print('#################################################################################################')
  print('# INITIIALIZING ENVIRONMENT CLASS')
  print('#################################################################################################')
  print('#')
  
  sink()

  return()
}






