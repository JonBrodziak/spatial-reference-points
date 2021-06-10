################################################################################################
# SRP.R
################################################################################################
# The R program file SRP.R extends case_study_v1 to include:
# (1) Standard functions with parameters
# (2) Spatial Reference Points (SRP)
# Jon Brodziak, PIFSC, jon.brodziak@noaa.gov
################################################################################################
# 27-April-2021
# Included RecAge parameter for first age in population NAA vector and recruitment lag
# RecAge can be set to age-0 or age-1
# If RecaAge=0 then array positions 1:NAge represent ages 0 to (NAge-1) & NAge-1 is a plusgroup
# If RecaAge=1 then array positions 1:NAge represent ages 1 to NAge & NAge is a plusgroup
# 
################################################################################################
#
# (1) MODEL INITIALIZATION
#
################################################################################################
# (1.1) INITIALIZE ANALYSIS CLASS
################################################################################################
# NOTE the age groups in each population consist of age=0 to age=9+ IF RecAge=0
# and the age group information is stored in array positions a=1 to a=10.
# NOTE the age groups in each population consist of age=1 to age=10+ IF RecAge=1
# and the age group information is stored in array positions a=1 to a=10.
# NOTE that females are indexed by g=1 and males by g=2 in arrays.
# NOTE that fleet 1 fishes in area 1 and fleet 2 fishes in area 2.
# NOTE that survey 1 operates in area 1 and survey 2 operates in area 2.
#-----------------------------------------------------------------------------------------------

# Remove all variable objects from the session memory
rm(list=ls(all.names=F))

# Close any open graphics windows
graphics.off()

# NOTE: NEED TO UPDATE WORKSTATION SETUP
# workstation setup for source.folder
# source.folder <- 'D:/path'

# laptop setup for source.folder
source.folder <- 'C:/Users/Jon.Brodziak/Desktop/MAS/2021_MAS_Spatial_Reference_Points/SRP/code/'

# Source the R functions for submodels

source(paste(source.folder,'Maturity.R',sep=""))

source(paste(source.folder,'Growth.R',sep=""))

source(paste(source.folder,'WeightAtLength.R',sep=""))

source(paste(source.folder,'StockRecruitment.R',sep=""))

source(paste(source.folder,'Selectivity.R',sep=""))

source(paste(source.folder,'Write_Input_Data.R',sep=""))

################################################################################################
# (1.2) SET MODEL DOMAIN PARAMETERS
################################################################################################

source(paste(source.folder,'Read_Input_Data.R',sep=""))

################################################################################################
# (1.3) INITIALIZE UNFISHED POPULATION EQUILIBRIUM CLASS
################################################################################################

# Compute maturity probabilities at age by gender, area, and population
#-----------------------------------------------------------------------------------------------
MaturityProbabilityAtAge <- array(rep(1.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))

source(paste(source.folder,'Calculate_Maturity_Probability_At_Age.R',sep=""))

# Dimension unfished mean lengths at age arrays
#-----------------------------------------------------------------------------------------------
MeanLengthStartOfYear <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
MeanLengthSpawning <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
MeanLengthCatch <- array(rep(0.0,DimPopulationFleetGenderAge),c(NPopulation,NFleet,NGender,NAge))
MeanLengthSurvey <- array(rep(0.0,DimPopulationSurveyGenderAge),c(NPopulation,NSurvey,NGender,NAge))

# Dimension unfished mean weights at age arrays
#-----------------------------------------------------------------------------------------------
MeanWeightStartOfYear <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
MeanWeightSpawning <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
MeanWeightCatch <- array(rep(0.0,DimPopulationFleetGenderAge),c(NPopulation,NFleet,NGender,NAge))
MeanWeightSurvey <- array(rep(0.0,DimPopulationSurveyGenderAge),c(NPopulation,NSurvey,NGender,NAge))

# Calculate unfished mean lengths and weights at age
#-----------------------------------------------------------------------------------------------
source(paste(source.folder,'Calculate_Unfished_Mean_Sizes_At_Age.R',sep=""))

# Store unfished equilibrium mean lengths and mean weights at age
#-----------------------------------------------------------------------------------------------
UnfishedMeanLengthStartOfYear <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
UnfishedMeanLengthSpawning <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
UnfishedMeanLengthCatch <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
UnfishedMeanLengthSurvey <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
UnfishedMeanWeightStartOfYear <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
UnfishedMeanWeightSpawning <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
UnfishedMeanWeightCatch <- array(rep(0.0,DimPopulationFleetGenderAge),c(NPopulation,NFleet,NGender,NAge))
UnfishedMeanWeightSurvey <- array(rep(0.0,DimPopulationSurveyGenderAge),c(NPopulation,NSurvey,NGender,NAge))

UnfishedMeanLengthStartOfYear <- MeanLengthStartOfYear
UnfishedMeanLengthSpawning <- MeanLengthSpawning
UnfishedMeanLengthCatch <- MeanLengthCatch
UnfishedMeanLengthSurvey <- MeanLengthSurvey
UnfishedMeanWeightStartOfYear <- MeanWeightStartOfYear
UnfishedMeanWeightSpawning <- MeanWeightSpawning
UnfishedMeanWeightCatch <- MeanWeightCatch
UnfishedMeanWeightSurvey <- MeanWeightSurvey

################################################################################################
# (1.4) INITIALIZE FISHED EQUILIBRIUM POPULATION CLASS
################################################################################################

# Calculate equilibrium fishery selectivity at age by gender, area (fleet), and population
#-----------------------------------------------------------------------------------------------
EquilibriumFisherySelectivityAtAge <- array(rep(1.0,DimPopulationFleetGenderAge),c(NPopulation,NFleet,NGender,NAge))

source(paste(source.folder,'Calculate_Equilibrium_Fishery_Selectivity_At_Age.R',sep=""))

# Set fishery selectivity during assessment time horizon to equal equilibrium fishery selectivity
#-----------------------------------------------------------------------------------------------
FisherySelectivityAtAge <- array(rep(1.0,DimPopulationFleetGenderAge),c(NPopulation,NFleet,NGender,NAge))

FisherySelectivityAtAge <- EquilibriumFisherySelectivityAtAge

# Calculate survey selectivity at age by gender, area, and population
#-----------------------------------------------------------------------------------------------
SurveySelectivityAtAge <- array(rep(1.0,DimPopulationSurveyGenderAge),c(NPopulation,NSurvey,NGender,NAge))

source(paste(source.folder,'Calculate_Survey_Selectivity_At_Age.R',sep=""))

# Dimension arrays
#-----------------------------------------------------------------------------------------------
EquilibriumFishingMortalityAtAge <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
EquilibriumTotalMortalityAtAge <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))

# Calculate fished equilibrium fishing mortality at age
#-----------------------------------------------------------------
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
      for (a in 1:NAge) {
        EquilibriumFishingMortalityAtAge[p,d,g,a] <- EquilibriumFisherySelectivityAtAge[p,d,g,a]*FishedEquilibriumFishingMortality[d,d]
        EquilibriumTotalMortalityAtAge[p,d,g,a] <- NaturalMortality[p,d,g,a]+EquilibriumFishingMortalityAtAge[p,d,g,a]
      }

source(paste(source.folder,'Calculate_Fished_Equilibrium_Mean_Size_At_NAge.R',sep=""))

# Store fished equilibrium mean lengths and mean weights for the plus group, all true ages are the same
#-----------------------------------------------------------------------------------------------
FishedMeanLengthStartOfYear <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
FishedMeanLengthSpawning <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
FishedMeanLengthCatch <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
FishedMeanLengthSurvey <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
FishedMeanWeightStartOfYear <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
FishedMeanWeightSpawning <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
FishedMeanWeightCatch <- array(rep(0.0,DimPopulationFleetGenderAge),c(NPopulation,NFleet,NGender,NAge))
FishedMeanWeightSurvey <- array(rep(0.0,DimPopulationSurveyGenderAge),c(NPopulation,NSurvey,NGender,NAge))

FishedMeanLengthStartOfYear <- MeanLengthStartOfYear
FishedMeanLengthSpawning <- MeanLengthSpawning
FishedMeanLengthCatch <- MeanLengthCatch
FishedMeanLengthSurvey <- MeanLengthSurvey
FishedMeanWeightStartOfYear <- MeanWeightStartOfYear
FishedMeanWeightSpawning <- MeanWeightSpawning
FishedMeanWeightCatch <- MeanWeightCatch
FishedMeanWeightSurvey <- MeanWeightSurvey

# Set stock-recruitment submodel parameters by gender, area, and population
#-----------------------------------------------------------------------------------------------
Recruitment.LogUnfishedR <- log(Recruitment.UnfishedR)

################################################################################################
# (1.5) INITIALIZE ENVIRONMENT CLASS
################################################################################################
# Placeholder
#-----------------------------------------------------------------------------------------------

################################################################################################
# WRITE INPUT DATA TO FILE
################################################################################################

Write_Input_Data(OutputFile)

################################################################################################
#
# (2) POPULATION LOOP OVER THREE MODEL TIME PERIODS:
# 1. UNFISHED EQUILIBRIUM
# 2. FISHED EQUILIBRIUM
# 3. ASSESSMENT TIME HORIZON
#
################################################################################################
#
# (2.1) CALCULATE UNFISHED EQUILIBRIUM
#       POPULATION NUMBERS AND SPAWNING BIOMASS BY POPULATION, AREA, AND GENDER
#
################################################################################################
#
# Dimension arrays for unfished equilibrium calculation
#-----------------------------------------------------------------------------------------------
EquilibriumNumbersAtAge <- array(rep(0.0,DimIterationPopulationAreaGenderAge),c(MaxIteration,NPopulation,NArea,NGender,NAge)) 
EquilibriumSpawningBiomass <- array(rep(0.0,DimIterationPopulationAreaGender),c(MaxIteration,NPopulation,NArea,NGender))
R.Iteration <- array(rep(0.0,DimIterationPopulationArea),c(MaxIteration,NPopulation,NArea))
 
# Calculate unfished equilibrium spawning biomasses
#-----------------------------------------------------------------------------------------------

source(paste(source.folder,'Calculate_Unfished_Equilibrium.R',sep=""))
 
################################################################################################
# (2.1) END
################################################################################################

################################################################################################
#
# (2.2) CALCULATE FISHED EQUILIBRIUM
#       POPULATION NUMBERS AND SPAWNING BIOMASS BY POPULATION, AREA, AND GENDER
#
################################################################################################

# Dimension arrays for fished equilibrium calculation
#-----------------------------------------------------------------
EquilibriumNumbersAtAge <- array(rep(0.0,DimIterationPopulationAreaGenderAge),c(MaxIteration,NPopulation,NArea,NGender,NAge)) 
EquilibriumSpawningBiomass <- array(rep(0.0,DimIterationPopulationAreaGender),c(MaxIteration,NPopulation,NArea,NGender))
R.Iteration <- array(rep(0.0,DimIterationPopulationArea),c(MaxIteration,NPopulation,NArea))
EquilibriumRecruitmentProduction <- array(rep(0.0,DimPopulationArea),c(NPopulation,NArea))

# Calculate fished equilibrium numbers at age 
#-----------------------------------------------------------------
source(paste(source.folder,'Calculate_Fished_Equilibrium.R',sep=""))

################################################################################################
#
# (2.3) CALCULATE ASSESSMENT PERIOD
#       POPULATION NUMBERS, SPAWNING BIOMASS, AND RECRUITMENT BY POPULATION, AREA, AND GENDER
#
################################################################################################
# Dimension arrays for assessment period calculation
#-----------------------------------------------------------------------------------------------
NumbersAtAge <- array(rep(0.0,DimYearPlusOnePopulationAreaGenderAge),c(NYearPlusOne,NPopulation,NArea,NGender,NAge))
AssessmentMeanLengthStartOfYear <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanLengthSpawning <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanLengthCatch <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanLengthSurvey <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanWeightStartOfYear <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanWeightSpawning <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanWeightCatch <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))
AssessmentMeanWeightSurvey <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))

FisheryCatchNumbersAtAgeByPopulation <- array(rep(0.0,DimYearPopulationFleetAreaGenderAge),c(NYear,NPopulation,NFleet,NArea,NGender,NAge))
FisheryCatchProportionAtAgeByPopulation <- array(rep(0.0,DimYearPopulationFleetAreaGenderAge),c(NYear,NPopulation,NFleet,NArea,NGender,NAge))
FisheryCatchNumbersAtAge <- array(rep(0.0,DimYearFleetAreaGenderAge),c(NYear,NFleet,NArea,NGender,NAge))
FisheryCatchProportionAtAge <- array(rep(0.0,DimYearFleetAreaGenderAge),c(NYear,NFleet,NArea,NGender,NAge))
FisheryCatchBiomassAtAgeByPopulation <- array(rep(0.0,DimYearPopulationFleetAreaGenderAge),c(NYear,NPopulation,NFleet,NArea,NGender,NAge))
FisheryCatchBiomass <- array(rep(0.0,DimYearFleetArea),c(NYear,NFleet,NArea))

TotalMortality  <- array(rep(0.0,DimYearPopulationAreaGenderAge),c(NYear,NPopulation,NArea,NGender,NAge))

SpawningNumbersAtAge <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
SpawningBiomass <- array(rep(0.0,DimYearPopulationAreaGender),c(NYear,NPopulation,NArea,NGender))

RecruitmentProduction <- array(rep(0.0,DimYearPopulationArea),c(NYear,NPopulation,NArea))
Recruitment <- array(rep(0.0,DimYearPopulationAreaGender),c(NYear,NPopulation,NArea,NGender))

SurveyCatchability <- array(rep(0.00001,DimPopulationSurvey),c(NPopulation,NSurvey))
SurveyCatchNumbersAtAgeByPopulation <- array(rep(0.0,DimYearPopulationSurveyGenderAge),c(NYear,NPopulation,NSurvey,NGender,NAge))
SurveyCatchProportionAtAgeByPopulation <- array(rep(0.0,DimYearPopulationSurveyGenderAge),c(NYear,NPopulation,NSurvey,NGender,NAge))
SurveyCatchNumbersAtAge <- array(rep(0.0,DimYearSurveyGenderAge),c(NYear,NSurvey,NGender,NAge))
SurveyCatchProportionAtAge <- array(rep(0.0,DimYearSurveyGenderAge),c(NYear,NSurvey,NGender,NAge))
SurveyCatchBiomassAtAgeByPopulation <- array(rep(0.0,DimYearPopulationSurveyGenderAge),c(NYear,NPopulation,NSurvey,NGender,NAge))
SurveyCatchBiomass <- array(rep(0.0,DimYearSurvey),c(NYear,NSurvey))

# Calculate numbers at age and QOI during assessment period 
#-----------------------------------------------------------------
source(paste(source.folder,'Calculate_Assessment_Period.R',sep=""))

################################################################################################
#
# (2.4) CALCULATE EQUILIBRIUM REFERENCE POINTS
#       MSY-BASED REFERENCE POINTS BY POPULATION, AREA, AND GLOBAL
#
################################################################################################
# Dimension arrays for reference point calculations
#-----------------------------------------------------------------
EquilibriumNumbersAtAge <- array(rep(0.0,DimIterationPopulationAreaGenderAge),c(MaxIteration,NPopulation,NArea,NGender,NAge)) 
EquilibriumSpawningBiomass <- array(rep(0.0,DimIterationPopulationAreaGender),c(MaxIteration,NPopulation,NArea,NGender))
R.Iteration <- array(rep(0.0,DimIterationPopulationArea),c(MaxIteration,NPopulation,NArea))
EquilibriumRecruitmentProduction <- array(rep(0.0,DimPopulationArea),c(NPopulation,NArea))

DimGridGrid <- NGrid*NGrid
DimGridGridPopulationAreaGender <- NGrid*NGrid*NPopulation*NArea*NGender
DimGridGridPopulationAreaGenderAge <- NGrid*NGrid*NPopulation*NArea*NGender*NAge

FGridFinalIteration <- array(rep(0,DimGridGrid),c(NGrid,NGrid))
FGridSpawningBiomass <- array(rep(0.0,DimGridGridPopulationAreaGender),c(NGrid,NGrid,NPopulation,NArea,NGender))
FGridNumbersAtAge <- array(rep(0.0,DimGridGridPopulationAreaGenderAge),c(NGrid,NGrid,NPopulation,NArea,NGender,NAge)) 

sink(file=OutputFile,append=TRUE,type="output")
# Check calculation of FGrid
#-----------------------------------------------------------------------------------------------
print(c("FGrid is ", FGrid))
print('_____________________________________________________________________________________________________')

FbyArea <- rep(0,length=NArea)

# Calculate equilibrium fishing mortality at age by area for 2 areas
#-----------------------------------------------------------------------------------------------
for (g1 in 1:NGrid)
  for (g2 in 1:NGrid)
    {
    FbyArea[1] <- FGrid[g1]
    FbyArea[2] <- FGrid[g2]
#    print(c("FbyArea is ", FbyArea))

# Calculate fished equilibrium fishing mortality by area over grid
#-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        {
        for (g in 1:NGender)
          for (a in 1:NAge) {
            EquilibriumFishingMortalityAtAge[p,d,g,a] <- EquilibriumFisherySelectivityAtAge[p,d,g,a]*FbyArea[d]
            }  
#     print(c("FbyArea for population p=",p , " in area=",d, " for gender=",g ," is ", EquilibriumFishingMortalityAtAge[p,d,g,])) 

#        print(c("FbyArea for population p=",p , " in area=",d, " for gender=",g ," is ", FbyArea[d]))
        }
    
    # Calculate equilibrium numbers at age for FbyArea
    #################################################################################################
    print('Calculating FGrid population numbers at age and spawning biomasses ...')
    print('_____________________________________________________________________________________________________')
    
    # ITERATION 1
    
    # Compute initial estimates of FGrid numbers at age 
    # and spawning biomasses by population, area, and 
    # gender where iteration i=1
    #-----------------------------------------------------------------
    
    i <- 1
    a <- 1
    
    R.Iteration[i,,] <- Recruitment.UnfishedR
    
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender)
        {        
          tmp <- 0.0
          for (k in 1:NArea)
          {
            tmp <- tmp + Recruitment.GenderFraction[p,k,g]*RecruitmentDistribution[p,k,d]*R.Iteration[i,p,k]
          }
          EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
        }
    
    # Compute initial fished numbers at age by population, area, and gender
    # For iteration i=1 and true age indexes [a=2:(NAge-1)]
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender)
          for (a in 2:(NAge-1))
          {
            EquilibriumNumbersAtAge[i,p,d,g,a] <- EquilibriumNumbersAtAge[i,p,d,g,(a-1)]*exp(-NaturalMortality[p,d,g,(a-1)]-EquilibriumFishingMortalityAtAge[p,d,g,(a-1)])
          }
    
    # Compute initial fished numbers at age by population, area, and gender
    # For iteration i=1 and plus-group age index [a=NAge]
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender)
          for (a in NAge:NAge)
          {
            EquilibriumNumbersAtAge[i,p,d,g,a] <- EquilibriumNumbersAtAge[i,p,d,g,(a-1)]*exp(-NaturalMortality[p,d,g,(a-1)]-EquilibriumFishingMortalityAtAge[p,d,g,(a-1)])
            EquilibriumNumbersAtAge[i,p,d,g,a] <- EquilibriumNumbersAtAge[i,p,d,g,a]/(1-exp(-NaturalMortality[p,d,g,a]-EquilibriumFishingMortalityAtAge[p,d,g,a]))
          }
    
    # Compute initial fished spawning biomass by population, area, and gender
    # For iteration i=1 and age indexes [a=1:NAge]
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender)
        {
          for (a in 1:NAge)
          {
            tmp <- EquilibriumNumbersAtAge[i,p,d,g,a]*FishedMeanWeightSpawning[p,d,g,a]*MaturityProbabilityAtAge[p,d,g,a]*exp(-SpawningZFraction[p,d]*(NaturalMortality[p,d,g,a]+EquilibriumFishingMortalityAtAge[p,d,g,a]))
            EquilibriumSpawningBiomass[i,p,d,g] <- EquilibriumSpawningBiomass[i,p,d,g] + tmp
          }
          # Rescale to Spawning Biomass Units
          EquilibriumSpawningBiomass[i,p,d,g] <- EquilibriumSpawningBiomass[i,p,d,g]/Recruitment.SpawningBiomassUnits[p,d]      
        }
    
    if (VerboseEquilibriumOutput == 1) {
      print('_____________________________________________________________________________________________________')
      print("Iteration:")
      print(i)
      print("Population 1 Female Numbers at Age and Spawning Biomass by Area")
      print(EquilibriumNumbersAtAge[i,1,1,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,1,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      print(EquilibriumNumbersAtAge[i,1,2,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,1,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      
      print("Population 2 Female Numbers at Age and Spawning Biomass by Area")
      print(EquilibriumNumbersAtAge[i,2,1,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,2,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      print(EquilibriumNumbersAtAge[i,2,2,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,2,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      print('_____________________________________________________________________________________________________')
    }
    # ITERATION 2
    
    i <- i+1
    
    # Compute recruitment production by population and area for i=2
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
      {
        if (Recruitment.model[p,d] == 1)
          parameters <- c(Recruitment.model[p,d],UnfishedSpawningBiomass[p,d,1], Recruitment.UnfishedR[p,d], Recruitment.steepness[p,d])
        R.Iteration[i,p,d] <- StockRecruitment(EquilibriumSpawningBiomass[(i-1),p,d,1],parameters)
      }
    
    # Compute fished recruitment strength by population, 
    # area, and gender for i=2
    #-----------------------------------------------------------------
    a <- 1
    for (p in 1:NPopulation) {
      for (d in 1:NArea)
        for (g in 1:NGender)
        {   
          tmp <- 0.0
          for (k in 1:NArea)
          {
            tmp <- tmp + Recruitment.GenderFraction[p,k,g]*RecruitmentDistribution[p,k,d]*R.Iteration[i,p,k]
          }
          EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
        }
    }
    
    # Compute fished numbers at age by population, area, and gender
    # for i=2 and true age indexes [a=2:(NAge-1)]
    #-----------------------------------------------------------------
    for (p in 1:NPopulation) {
      for (d in 1:NArea)
        for (g in 1:NGender)
          for (a in 2:(NAge-1))
          {
            tmp <- 0.0
            for (k in 1:NArea)
            {
              tmp <- tmp + MovementProbability[p,k,d]*EquilibriumNumbersAtAge[(i-1),p,k,g,(a-1)]*exp(-NaturalMortality[p,k,g,(a-1)]-EquilibriumFishingMortalityAtAge[p,k,g,(a-1)])
            }
            EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
          }
    }
    
    # Compute fished numbers at age by population, area, and gender
    # For iteration i=2 and plus-group age index [a=NAge]
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
            # Age-(A-1) survivors that remained or immigrated to area s
            tmp <- tmp + MovementProbability[p,k,d]*EquilibriumNumbersAtAge[(i-1),p,k,g,(a-1)]*exp(-NaturalMortality[p,k,g,(a-1)]-EquilibriumFishingMortalityAtAge[p,k,g,(a-1)])
            # Age-A survivors that remained or immigrated to area s
            tmp <- tmp + MovementProbability[p,k,d]*EquilibriumNumbersAtAge[(i-1),p,k,g,a]*exp(-NaturalMortality[p,k,g,a]-EquilibriumFishingMortalityAtAge[p,k,g,a])
          }
          EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
        }
    }
    
    # Compute fished spawning biomass by population, area, and gender
    # For iteration i=2 and age indexes [a=1:NAge]
    #-----------------------------------------------------------------
    for (p in 1:NPopulation)
      for (d in 1:NArea)
        for (g in 1:NGender)
        {
          tmp <- 0.0
          for (a in 1:NAge)
          {
            tmp <- tmp + EquilibriumNumbersAtAge[i,p,d,g,a]*FishedMeanWeightSpawning[p,d,g,a]*MaturityProbabilityAtAge[p,d,g,a]*exp(-SpawningZFraction[p,d]*(NaturalMortality[p,d,g,a]+EquilibriumFishingMortalityAtAge[p,d,g,a]))
          }   
          EquilibriumSpawningBiomass[i,p,d,g] <- tmp
          # Rescale to Spawning Biomass Units
          EquilibriumSpawningBiomass[i,p,d,g] <- EquilibriumSpawningBiomass[i,p,d,g]/Recruitment.SpawningBiomassUnits[p,d] 
        }
    
    if (VerboseEquilibriumOutput == 1) {
      print('_____________________________________________________________________________________________________')
      print("Iteration:")
      print(i)
      print("Population 1 Female Numbers at Age and Spawning Biomass by Area")
      print(EquilibriumNumbersAtAge[i,1,1,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,1,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      print(EquilibriumNumbersAtAge[i,1,2,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,1,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      
      print("Population 2 Female Numbers at Age and Spawning Biomass by Area")
      print(EquilibriumNumbersAtAge[i,2,1,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,2,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      print(EquilibriumNumbersAtAge[i,2,2,1,]/OutputNumbersUnits)
      print(EquilibriumSpawningBiomass[i,2,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      print('_____________________________________________________________________________________________________')
    }
    
    # BEGIN WHILE LOOP
    
    HasConverged <- 0
    
    while ((i<MaxIteration) && (HasConverged==0))
    {
      i <- i+1
      
      # Compute recruitment production by population and area for i
      #-----------------------------------------------------------------
      for (p in 1:NPopulation)
        for (d in 1:NArea)
        {
          if (Recruitment.model[p,d] == 1)
            parameters <- c(Recruitment.model[p,d],UnfishedSpawningBiomass[p,d,1], Recruitment.UnfishedR[p,d], Recruitment.steepness[p,d])
          R.Iteration[i,p,d] <- StockRecruitment(EquilibriumSpawningBiomass[(i-1),p,d,1],parameters)
        }
      
      # Compute fished recruitment strength by population, 
      # area, and gender for i
      #-----------------------------------------------------------------
      a <- 1
      for (p in 1:NPopulation) {
        for (d in 1:NArea)
          for (g in 1:NGender)
          {   
            tmp <- 0.0
            for (k in 1:NArea)
            {
              tmp <- tmp + Recruitment.GenderFraction[p,k,g]*RecruitmentDistribution[p,k,d]*R.Iteration[i,p,k]
            }
            EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
          }
      }
      
      # Compute fished numbers at age by population, area, and gender
      # for i and true age indexes [a=2:(NAge-1)]
      #-----------------------------------------------------------------
      for (p in 1:NPopulation) {
        for (d in 1:NArea)
          for (g in 1:NGender)
            for (a in 2:(NAge-1))
            {
              tmp <- 0.0
              for (k in 1:NArea)
              {
                tmp <- tmp + MovementProbability[p,k,d]*EquilibriumNumbersAtAge[(i-1),p,k,g,(a-1)]*exp(-NaturalMortality[p,k,g,(a-1)]-EquilibriumFishingMortalityAtAge[p,k,g,(a-1)])
              }
              EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
            }
      }
      
      # Compute fished numbers at age by population, area, and gender
      # For iteration i and plus-group age index [a=NAge]
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
              # Age-(A-1) survivors that remained or immigrated to area s
              tmp <- tmp + MovementProbability[p,k,d]*EquilibriumNumbersAtAge[(i-1),p,k,g,(a-1)]*exp(-NaturalMortality[p,k,g,(a-1)]-EquilibriumFishingMortalityAtAge[p,k,g,(a-1)])
              # Age-A survivors that remained or immigrated to area s
              tmp <- tmp + MovementProbability[p,k,d]*EquilibriumNumbersAtAge[(i-1),p,k,g,a]*exp(-NaturalMortality[p,k,g,a]-EquilibriumFishingMortalityAtAge[p,k,g,a])
            }
            EquilibriumNumbersAtAge[i,p,d,g,a] <- tmp
          }
      }
      
      # Compute fished spawning biomass by population, area, and gender
      # For iteration i and age indexes [a=1:NAge]
      #-----------------------------------------------------------------
      for (p in 1:NPopulation)
        for (d in 1:NArea)
          for (g in 1:NGender)
          {
            tmp <- 0.0
            for (a in 1:NAge)
            {
              tmp <- tmp + EquilibriumNumbersAtAge[i,p,d,g,a]*FishedMeanWeightSpawning[p,d,g,a]*MaturityProbabilityAtAge[p,d,g,a]*exp(-SpawningZFraction[p,d]*(NaturalMortality[p,d,g,a]+EquilibriumFishingMortalityAtAge[p,d,g,a]))
            }   
            EquilibriumSpawningBiomass[i,p,d,g] <- tmp
            # Rescale to Spawning Biomass Units
            EquilibriumSpawningBiomass[i,p,d,g] <- EquilibriumSpawningBiomass[i,p,d,g]/Recruitment.SpawningBiomassUnits[p,d] 
          }
      
      # Compute L1 distance between iterates of unfished spawning
      # biomass by population, area, and gender and check for convergence
      #-----------------------------------------------------------------
      Distance <- 0.0
      
      for (p in 1:NPopulation)
        for (d in 1:NArea)
          for (g in 1:NGender)
          {
            Distance <- Distance + abs(EquilibriumSpawningBiomass[i,p,d,g]-EquilibriumSpawningBiomass[(i-1),p,d,g])
          }
      
      if ((i %% 100 == 0) && (VerboseEquilibriumOutput != 1))
      {
        print('_____________________________________________________________________________________________________')
        print("Iteration:")
        print(i)
        print("Distance:")
        print(Distance)
        print("Population 1 FGrid Female Numbers at Age and Spawning Biomass by Area")
        print(EquilibriumNumbersAtAge[i,1,1,1,]/OutputNumbersUnits)
        print(EquilibriumSpawningBiomass[i,1,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
        print(EquilibriumNumbersAtAge[i,1,2,1,]/OutputNumbersUnits)
        print(EquilibriumSpawningBiomass[i,1,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
        
        print("Population 2 FGrid Female Numbers at Age and Spawning Biomass by Area")
        print(EquilibriumNumbersAtAge[i,2,1,1,]/OutputNumbersUnits)
        print(EquilibriumSpawningBiomass[i,2,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
        print(EquilibriumNumbersAtAge[i,2,2,1,]/OutputNumbersUnits)
        print(EquilibriumSpawningBiomass[i,2,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
      }  
      
      if (Distance < EquilibriumConvergenceCriterion)
      {
        HasConverged <- 1
        FinalIteration <- i
        print('_____________________________________________________________________________________________________')
        print('Calculation of FGrid population numbers at age and spawning biomasses has converged at iteration:')
        print(FinalIteration)
      }
      
      if (i == MaxIteration)
      {
        print('_____________________________________________________________________________________________________')
        print('Calculation of FGrid population numbers at age and spawning biomasses did not converge with distance:')
        print(Distance)
      }
      
    }
    
    # END WHILE LOOP
    
    # Set FGrid population numbers and spawning
    # biomass by population, area, and gender
    #-----------------------------------------------------------------
    FGridFinalIteration[g1,g2] <- FinalIteration
    FGridSpawningBiomass[g1,g2,,,] <- EquilibriumSpawningBiomass[FinalIteration,,,]
    FGridNumbersAtAge[g1,g2,,,,] <- EquilibriumNumbersAtAge[FinalIteration,,,,]

    print('_____________________________________________________________________________________________________')
    print(c("FbyArea is ", FbyArea))
    print('FGrid Results')
    print('_____________________________________________________________________________________________________')
    print('Population 1 FGrid Female Numbers at Age in Area 1')
    print(FGridNumbersAtAge[g1,g2,1,1,1,]/OutputNumbersUnits)
    print('Population 1 FGrid Female Spawning Biomass in Area 1')
    print(FGridSpawningBiomass[g1,g2,1,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
    print('_____________________________________________________________________________________________________')
    print('Population 1 FGrid Female Numbers at Age in Area 2')
    print(FGridNumbersAtAge[g1,g2,1,2,1,]/OutputNumbersUnits)
    print('Population 1 FGrid Female Spawning Biomass in Area 2')
    print(FGridSpawningBiomass[g1,g2,1,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
    print('_____________________________________________________________________________________________________')
    print('Population 2 FGrid Female Numbers at Age in Area 1')
    print(FGridNumbersAtAge[g1,g2,2,1,1,]/OutputNumbersUnits)
    print('Population 2 FGrid Female Spawning Biomass in Area 1')
    print(FGridSpawningBiomass[g1,g2,2,1,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
    print('_____________________________________________________________________________________________________')
    print('Population 2 FGrid Female Numbers at Age in Area 2')
    print(FGridNumbersAtAge[g1,g2,2,2,1,]/OutputNumbersUnits)
    print('Population 2 FGrid Female Spawning Biomass in Area 2')
    print(FGridSpawningBiomass[g2,g2,2,2,1]*Recruitment.SpawningBiomassUnits/OutputBiomassUnits)
    print('_____________________________________________________________________________________________________')
     
  }
# END of F GRID LOOP
####################################################################################################

sink()

################################################################################################
#
# (3)	WRITE MODEL RESULTS TO OUTPUT FILE
#  
################################################################################################

# Output numbers at age and QOI during assessment period 
#-----------------------------------------------------------------
source(paste(source.folder,'Write_Output_Data.R',sep=""))

################################################################################################
# END 
################################################################################################
