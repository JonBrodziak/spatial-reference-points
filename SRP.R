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
source.folder <- 'C:/Users/Jon.Brodziak/Desktop/MAS/2021_MAS_Reference_Points/SRP/code/'

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
# (1.3) INITIALIZE POPULATION CLASS
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

# Set stock-recruitment submodel parameters by gender, area, and population
#-----------------------------------------------------------------------------------------------
Recruitment.LogUnfishedR <- log(Recruitment.UnfishedR)

################################################################################################
# (1.4) INITIALIZE OBSERVATION CLASS
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

################################################################################################
# (1.5) INITIALIZE ENVIRONMENT CLASS
################################################################################################
# Placeholder
#-----------------------------------------------------------------------------------------------

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
# Dimension arrays
#-----------------------------------------------------------------------------------------------
 EquilibriumNumbersAtAge <- array(rep(0.0,DimIterationPopulationAreaGenderAge),c(MaxIteration,NPopulation,NArea,NGender,NAge)) 
 EquilibriumSpawningBiomass <- array(rep(0.0,DimIterationPopulationAreaGender),c(MaxIteration,NPopulation,NArea,NGender))

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

# Dimension arrays
#-----------------------------------------------------------------
EquilibriumNumbersAtAge <- array(rep(0.0,DimIterationPopulationAreaGenderAge),c(MaxIteration,NPopulation,NArea,NGender,NAge)) 
EquilibriumSpawningBiomass <- array(rep(0.0,DimIterationPopulationAreaGender),c(MaxIteration,NPopulation,NArea,NGender))
EquilibriumFishingMortalityAtAge <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
EquilibriumTotalMortalityAtAge <- array(rep(0.0,DimPopulationAreaGenderAge),c(NPopulation,NArea,NGender,NAge))
EquilibriumRecruitmentProduction <- array(rep(0.0,DimPopulationArea),c(NPopulation,NArea))

# Calculate equilibrium fishing mortality at age
#-----------------------------------------------------------------
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for (g in 1:NGender)
      for (a in 1:(NAge)) {
        EquilibriumFishingMortalityAtAge[p,d,g,a] <- EquilibriumFisherySelectivityAtAge[p,d,g,a]*SimEquilibriumFishingMortality[d,d]
        EquilibriumTotalMortalityAtAge[p,d,g,a] <- NaturalMortality[p,k,g,a]+EquilibriumFishingMortalityAtAge[p,k,g,a]
      }

source(paste(source.folder,'Calculate_Fished_Equilibrium_Mean_Size_At_NAge.R',sep=""))

# Store fished equilibrium mean lengths and mean weights at age
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

# Calculate fished equilibrium numbers at age 
#-----------------------------------------------------------------
source(paste(source.folder,'Calculate_Fished_Equilibrium.R',sep=""))

################################################################################################
# WRITE INPUT DATA TO FILE
################################################################################################

Write_Input_Data(OutputFile)

################################################################################################
#
# (2.3) CALCULATE ASSESSMENT TIME PERIOD
#       POPULATION NUMBERS, SPAWNING BIOMASS, AND RECRUITMENT BY POPULATION, AREA, AND GENDER
#
################################################################################################
# Dimension Arrays
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
# (3)	WRITE MODEL RESULTS TO OUTPUT FILE
#  
################################################################################################

# Output numbers at age and QOI during assessment period 
#-----------------------------------------------------------------
source(paste(source.folder,'Write_Output_Data.R',sep=""))

################################################################################################
# END 
################################################################################################
