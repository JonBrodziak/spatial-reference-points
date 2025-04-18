################################################################################################
# Read_Input_Data.R
################################################################################################
# The R program file Read_Input_Data.R reads the input data and output file name and
# prints the input data to screen
# Jon Brodziak, PIFSC, jon.brodziak@noaa.gov 20-APR-2021
################################################################################################

# Read input and output file names from console
#-----------------------------------------------------------------------------------------------
print("Enter the input file name (with no spaces):")
InputFile <- scan(what="",quiet=T)
print("Input file name is:")
print(InputFile)
print("Enter the output file name (with no spaces):")
OutputFile <- scan(what="",quiet=T)
print("Output file name is:")
print(OutputFile)
SkipLines <- 1

print('_____________________________________________________________________________________________________')
print('Reading input data from input file ...')
print('_____________________________________________________________________________________________________')

# Read the verbose output to console flag 
#-----------------------------------------------------------------------------------------------
VerboseEquilibriumOutput <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Verbose console output (1=True)")
print(VerboseEquilibriumOutput) 
SkipLines <- SkipLines+2

# Read the number of populations
#-----------------------------------------------------------------------------------------------
NPopulation <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of populations")
print(NPopulation)
SkipLines <- SkipLines+2

# Read the number of age classes including the plus group (age-0 to age-(NAge-1))
#-----------------------------------------------------------------------------------------------
NAge <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of ages")
print(NAge)
SkipLines <- SkipLines+2

# Read the recruitment age, or first age in the NAA vector (Either RecAge=0 or RecAge=1)
#-----------------------------------------------------------------------------------------------
RecAge <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Recruitment age or first age in population numbers at age vector")
print(RecAge)
SkipLines <- SkipLines+2

# Read the maximum age in the plus group
#-----------------------------------------------------------------------------------------------
MaxPlusGroupAge <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Maximum plus group age")
print(MaxPlusGroupAge)
SkipLines <- SkipLines+2

# Read the number of genders
#-----------------------------------------------------------------------------------------------
NGender <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of genders")
print(NGender)
SkipLines <- SkipLines+2

# Read the number of areas
#-----------------------------------------------------------------------------------------------
NArea <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of areas")
print(NArea)
SkipLines <- SkipLines+2

# Read the number of fleets
#-----------------------------------------------------------------------------------------------
NFleet <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of fleets")
print(NFleet)
SkipLines <- SkipLines+2

# Read the fishing area for each fleet (1 area per fleet)
#-----------------------------------------------------------------------------------------------
FleetArea <- rep(0,length=NFleet)
FleetArea <- scan(InputFile, skip=SkipLines, n=2, quiet=T)
for (v in 1:NFleet) {
  print(c("Area for fleet",v))
  print(FleetArea[v]) 
}
SkipLines <- SkipLines+2

# Read the number of surveys
#-----------------------------------------------------------------------------------------------
NSurvey <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of surveys")
print(NSurvey)
SkipLines <- SkipLines+2

# Read the sampling area for each survey (1 area per survey)
#-----------------------------------------------------------------------------------------------
SurveyArea <- rep(0,length=NSurvey)
SurveyArea <- scan(InputFile, skip=SkipLines, n=2, quiet=T)
for (v in 1:NSurvey) {
  print(c("Area for survey",v))
  print(SurveyArea[v]) 
}
SkipLines <- SkipLines+2

# Read the number of years in the assessment time horizon
#-----------------------------------------------------------------------------------------------
NYear <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of years")
print(NYear)
SkipLines <- SkipLines+2

# Read the number of seasons in a year
#-----------------------------------------------------------------------------------------------
NSeason <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Number of seasons")
print(NSeason)
SkipLines <- SkipLines+2

# Read the maximum number of iterations to compute unfished and fished equilibrium
#-----------------------------------------------------------------------------------------------
MaxIteration <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Maximum number of iterations for equilibrium calculations")
print(MaxIteration)
SkipLines <- SkipLines+2

# Set convergence criterion for calculating equilibrium spawning biomasses
#-----------------------------------------------------------------------------------------------
EquilibriumConvergenceCriterion <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Convergence criterion for equilibrium calculations")
print(EquilibriumConvergenceCriterion)
SkipLines <- SkipLines+2

# Read the output units for population numbers (numbers of fish)
#-----------------------------------------------------------------------------------------------
OutputNumbersUnits <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Population numbers units for output")
print(OutputNumbersUnits)
SkipLines <- SkipLines+2

# Read the output units for population biomass (kilograms of fish)
#-----------------------------------------------------------------------------------------------
OutputBiomassUnits <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Population biomass units for output")
print(OutputBiomassUnits)
SkipLines <- SkipLines+2

# Read the maximum value of fishing mortality for reference point calculations
#-----------------------------------------------------------------------------------------------
MaxF <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Maximum fishing mortality for reference points")
print(MaxF)
SkipLines <- SkipLines+2

# Read the mesh for fishing mortality (e.g., FMesh=0.01) used in reference point calculations
#-----------------------------------------------------------------------------------------------
FMesh <- scan(InputFile, skip=SkipLines, n=1, quiet=T)
print("Fishing mortality mesh used for reference point calculations")
print(FMesh)
SkipLines <- SkipLines+2

# Set number of grid points for calculating reference points
#-----------------------------------------------------------------------------------------------
NGrid <- floor(MaxF/FMesh) + 1

################################################################################################
# Set array dimensions
#-----------------------------------------------------------------------------------------------
NYearPlusOne <- NYear+1
DimGrid <- NGrid                                                                                                         
DimAreaGender <- NArea*NGender
DimIterationPopulationAreaGenderAge <- MaxIteration*NPopulation*NArea*NGender*NAge
DimIterationPopulationAreaGender <- MaxIteration*NPopulation*NArea*NGender
DimIterationPopulationArea <- MaxIteration*NPopulation*NArea
DimFleetArea <- NFleet*NArea
DimPopulationAreaGenderAge <- NPopulation*NArea*NGender*NAge
DimPopulationFleetGenderAge <- NPopulation*NFleet*NGender*NAge
DimPopulationFleetGender <- NPopulation*NFleet*NGender
DimPopulationAreaGender <- NPopulation*NArea*NGender
DimPopulationAreaAge <- NPopulation*NArea*NAge
DimPopulationAreaArea <- NPopulation*NArea*NArea
DimPopulationArea <- NPopulation*NArea
DimPopulationFleet <- NPopulation*NFleet
DimPopulationSurveyGenderAge <- NPopulation*NSurvey*NGender*NAge
DimPopulationSurveyGender <- NPopulation*NSurvey*NGender
DimPopulationSurvey <- NPopulation*NSurvey
DimSurveyGender <- NSurvey*NGender
DimYearSurvey <- NYear*NSurvey
DimYearPopulationSurveyGenderAge <- NYear*NPopulation*NSurvey*NGender*NAge
DimYearPopulationFleetAreaGenderAge <- NYear*NPopulation*NFleet*NArea*NGender*NAge
DimYearPopulationAreaGenderAge <- NYear*NPopulation*NArea*NGender*NAge
DimYearPlusOnePopulationAreaGenderAge <- NYearPlusOne*NPopulation*NArea*NGender*NAge
DimYearPopulationAreaGender <- NYear*NPopulation*NArea*NGender
DimYearSurveyGenderAge <- NYear*NSurvey*NGender*NAge
DimYearPopulationArea <- NYear*NPopulation*NArea
DimYearFleetAreaGenderAge <- NYear*NFleet*NArea*NGender*NAge
DimYearAreaGenderAge <- NYear*NArea*NGender*NAge
DimYearFleetArea <- NYear*NFleet*NArea
DimYearArea <- NYear*NArea
################################################################################################

# Set F grid by area for calculating reference points
#-----------------------------------------------------------------------------------------------
FGrid <- array(rep(0.0,DimGrid),c(NGrid))
FGrid <- seq(0.0,MaxF,FMesh)

# Read 2D array for timing of start of year by area and population
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
StartZFraction <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation) {
  print(c("Z fraction by area for population",p))
  print(StartZFraction[p,]) 
}
SkipLines <- SkipLines+NPopulation+1

# Read 2D array for timing of spawning by area and population
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
SpawningZFraction <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation) {
  print(c("Start of year Z fraction by area for population",p))
  print(SpawningZFraction[p,]) 
}
SkipLines <- SkipLines+NPopulation+1

# Read 2D array for timing of catch by area and population
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
CatchZFraction <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation) {
  print(c("Catch Z fraction by area for population",p))
  print(CatchZFraction[p,]) 
}
SkipLines <- SkipLines+NPopulation+1

# Read 2D array for timing of surveys by area and population
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
SurveyZFraction <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation) {
  print(c("Survey Z fraction by area for population",p))
  print(SurveyZFraction[p,]) 
}
SkipLines <- SkipLines+NPopulation+1

# Read 3D array for population movement probability parameters
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationAreaArea, skip=SkipLines, quiet=T)
MovementProbability <- array(tmp,c(NPopulation,NArea,NArea))
for (p in 1:NPopulation) {
  print(c("Movement probabilities by area for population",p))
  print(MovementProbability[p,,]) 
}
SkipLines <- SkipLines+DimPopulationArea+1

# Read 3D array for recruitment distribution probability parameters
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationAreaArea, skip=SkipLines, quiet=T)
RecruitmentDistribution <- array(tmp,c(NPopulation,NArea,NArea))
for (p in 1:NPopulation) {
  print(c("Recruitment Distribution by area for population",p))
  print(RecruitmentDistribution[p,,]) 
}
SkipLines <- SkipLines+DimPopulationArea+1

# Read 3D array for natural mortality at age parameters
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationAreaGenderAge, skip=SkipLines, quiet=T)
NaturalMortality <- array(tmp,c(NPopulation,NArea,NGender,NAge))
for (p in 1:NPopulation)
  for (d in 1:NArea)
    for(g in 1:NGender) {
      print(c("Natural mortality at age for population",p,"in area",d,"and gender",g))
      print(NaturalMortality[p,d,g,])
    }
SkipLines <- SkipLines+NAge+1

# Read maturity model type
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
Maturity.model <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation){
    print(c("Maturity model by area for population",p))
    print(Maturity.model[p,])
  }
SkipLines <- SkipLines+NArea+1

if (Maturity.model[1,1] == 1) {
# Read 3D array for maturity a50 parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Maturity.a50 <- array(tmp,c(NPopulation,NArea,NGender))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Maturity a50 parameter by gender for population",p,"in area",d))
      print(Maturity.a50[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for maturity slope parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Maturity.slope <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Maturity slope parameter by gender for population",p,"in area",d))
      print(Maturity.slope[p,d,])
    }
  SkipLines <- SkipLines+NGender+1
}

# Read 2D array for growth models by population and area
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
Growth.model <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation){
    print(c("Growth model by area for population",p))
    print(Growth.model[p,])
  }
SkipLines <- SkipLines+NArea+1

if (Growth.model[1,1] == 1) {
# Growth.model == 1 is the modified von Bertalanffy growth curve
# Read 3D array for length at age Amin parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.Amin <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age Amin parameter by gender for population",p,"in area",d))
      print(Length.Amin[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for length at age Amax parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.Amax <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age Amax parameter by gender for population",p,"in area",d))
      print(Length.Amax[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for length at age Lmin parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.Lmin <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age Lmin parameter by gender for population",p,"in area",d))
      print(Length.Lmin[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for length at age Lmax parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.Lmax <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age Lmax parameter by gender for population",p,"in area",d))
      print(Length.Lmax[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for length at age c parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.c <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age c parameter by gender for population",p,"in area",d))
      print(Length.c[p,d,])
    }
  SkipLines <- SkipLines+NGender+1
}

if (Growth.model[1,1] == 2) {
# Growth.model == 2 is the standard von Bertalanffy growth curve
# Read 3D array for length at age Linf parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.Linf <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age Linf parameter by gender for population",p,"in area",d))
      print(Length.Linf[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for length at age K parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.K <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age K parameter by gender for population",p,"in area",d))
      print(Length.K[p,d,])
    }
  SkipLines <- SkipLines+NGender+1

# Read 3D array for length at age t0 parameters
#-----------------------------------------------------------------------------------------------
  tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
  Length.t0 <- array(tmp,c(NPopulation,NArea,NArea))
  for (p in 1:NPopulation)
    for (d in 1:NArea) {
      print(c("Length at age t0 parameter by gender for population",p,"in area",d))
      print(Length.t0[p,d,])
    }
  SkipLines <- SkipLines+NGender+1
}

# Read weight at length model type
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
WeightAtLength.model <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation){
    print(c("Weight at length model by area for population",p))
    print(WeightAtLength.model[p,])
  }
SkipLines <- SkipLines+NArea+1

if (WeightAtLength.model[1,1] == 1) {
	# Read 3D array for weight at length A parameters
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
	WeightAtLength.A <- array(tmp,c(NPopulation,NArea,NArea))
	for (p in 1:NPopulation)
	  for (d in 1:NArea) {
		print(c("Weight at length at age A parameter by gender for population",p,"in area",d))
		print(WeightAtLength.A[p,d,])
	  }
	SkipLines <- SkipLines+NGender+1

	# Read 3D array for weight at length B parameters
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
	WeightAtLength.B <- array(tmp,c(NPopulation,NArea,NArea))
	for (p in 1:NPopulation)
	  for (d in 1:NArea) {
		print(c("Weight at length at age B parameter by gender for population",p,"in area",d))
		print(WeightAtLength.B[p,d,])
	  }
	SkipLines <- SkipLines+NGender+1
}

# Read 2D array for units of spawning used in recruitment models
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
Recruitment.SpawningBiomassUnits <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation) {
  print(c("Spawning biomass units for recruitment by area for population",p))
  print(Recruitment.SpawningBiomassUnits[p,]) 
}
SkipLines <- SkipLines+NArea+1

# Read 3D array for gender fractions of recruitment output
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationAreaGender, skip=SkipLines, quiet=T)
Recruitment.GenderFraction <- array(tmp,c(NPopulation,NArea,NArea))
for (p in 1:NPopulation)
  for (d in 1:NArea) {
    print(c("Gender fraction for recruitment by gender for population",p,"in area",d))
    print(Recruitment.GenderFraction[p,d,])
  }
SkipLines <- SkipLines+NGender+1

# Read recruitment model type
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
Recruitment.model <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation){
    print(c("Weight at length model by area for population",p))
    print(Recruitment.model[p,])
  }
SkipLines <- SkipLines+NArea+1

if (Recruitment.model[1,1] == 1) {
	# Read 2D array for steepness parameter used in recruitment models
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
	Recruitment.steepness <- array(tmp,c(NPopulation,NArea))
	for (p in 1:NPopulation) {
	  print(c("Steepness parameters for recruitment by area for population",p))
	  print(Recruitment.steepness[p,]) 
	}
	SkipLines <- SkipLines+NArea+1

	# Read 2D array for unfished recruitment parameter used in recruitment models
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
	Recruitment.UnfishedR <- array(tmp,c(NPopulation,NArea))
	for (p in 1:NPopulation) {
	  print(c("Unfished recruitment parameters for recruitment by area for population",p))
	  print(Recruitment.UnfishedR[p,]) 
	}
	SkipLines <- SkipLines+NArea+1

	# Read 2D array for sigmaR parameters used in recruitment models
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
	Recruitment.sigmaR <- array(tmp,c(NPopulation,NArea))
	for (p in 1:NPopulation) {
	  print(c("Recruitment sigmaR parameters for recruitment by area for population",p))
	  print(Recruitment.sigmaR[p,]) 
	}
	SkipLines <- SkipLines+NArea+1
}

# Read 2D array for fished equilibrium fishing mortality parameters to initialize model
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimFleetArea, skip=SkipLines, quiet=T)
FishedEquilibriumFishingMortality <- array(tmp,c(NFleet,NArea))
for (p in 1:NPopulation) {
  print(c("Fished equilibrium fishing mortality parameters by area for population",p))
  print(FishedEquilibriumFishingMortality[p,]) 
}
SkipLines <- SkipLines+NArea+1
#
# Read 3D array for fishing mortality parameters in the assessment time horizon
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimYearFleetArea, skip=SkipLines, quiet=T)
AssessmentFishingMortality <- array(tmp,c(NYear,NFleet,NArea))
for (v in 1:NFleet)
  for (d in 1:NArea) {
    print(c("Assessment fishing mortality parameters by year for fleet",v,"in area",d))
    print(AssessmentFishingMortality[,v,d]) 
}
SkipLines <- SkipLines+NArea+1

# Read fishery selectivity model type
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
FisherySelectivity.model <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation){
    print(c("Fishery selectivity model by area for population",p))
    print(FisherySelectivity.model[p,])
  }
SkipLines <- SkipLines+NArea+1

if (FisherySelectivity.model[1,1] == 1) {
	# Read 3D array for fishery selectivity a50 parameters
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationFleetGender, skip=SkipLines, quiet=T)
	FisherySelectivity.a50 <- array(tmp,c(NPopulation,NFleet,NGender))
	for (p in 1:NPopulation)
	  for (v in 1:NFleet) {
		print(c("Fishery selectivity a50 parameters for population",p,"and fleet",v))
		print(FisherySelectivity.a50[p,v,])
	  }
	SkipLines <- SkipLines+NGender+1

	# Read 3D array for fishery selectivity slope parameters
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationFleetGender, skip=SkipLines, quiet=T)
	FisherySelectivity.slope <- array(tmp,c(NPopulation,NFleet,NGender))
	for (p in 1:NPopulation)
	  for (v in 1:NFleet) {
		print(c("Fishery selectivity slope parameters for population",p,"and fleet",v))
		print(FisherySelectivity.slope[p,v,])
	  }
	SkipLines <- SkipLines+NGender+1
}

# Read survey selectivity model type
#-----------------------------------------------------------------------------------------------
tmp <- scan(InputFile, nmax=DimPopulationArea, skip=SkipLines, quiet=T)
SurveySelectivity.model <- array(tmp,c(NPopulation,NArea))
for (p in 1:NPopulation){
    print(c("Survey selectivity model by area for population",p))
    print(SurveySelectivity.model[p,])
  }
SkipLines <- SkipLines+NArea+1

if (SurveySelectivity.model[1,1] == 1) {
	# Read 3D array for survey selectivity a50 parameters
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationSurveyGender, skip=SkipLines, quiet=T)
	SurveySelectivity.a50 <- array(tmp,c(NPopulation,NSurvey,NGender))
	for (p in 1:NPopulation)
	  for (I in 1:NSurvey) {
		print(c("Survey selectivity a50 parameters for population",p,"and survey",I))
		print(SurveySelectivity.a50[p,I,])
	  }
	SkipLines <- SkipLines+NGender+1

	# Read 3D array for survey selectivity slope parameters
	#-----------------------------------------------------------------------------------------------
	tmp <- scan(InputFile, nmax=DimPopulationSurveyGender, skip=SkipLines, quiet=T)
	SurveySelectivity.slope <- array(tmp,c(NPopulation,NSurvey,NGender))
	for (p in 1:NPopulation)
	  for (I in 1:NSurvey) {
		print(c("Survey selectivity slope parameters for population",p,"and survey",I))
		print(SurveySelectivity.slope[p,I,])
	  }
	SkipLines <- SkipLines+NGender+1
}



