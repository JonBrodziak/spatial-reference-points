################################################################################################
# Calculate_Equilibrium_Fishery_Selectivity_At_Age.R
################################################################################################

# Compute equilibrium fishery selectivity at age by gender, area (fleet), and population
#-----------------------------------------------------------------------------------------------
for (p in 1:NPopulation)
  for (v in 1:NFleet)
    if (FisherySelectivity.model[p,FleetArea[v]] == 1) {
      parameters <- c(FisherySelectivity.model[p,FleetArea[v]],FisherySelectivity.a50[p,FleetArea[v],g],FisherySelectivity.slope[p,FleetArea[v],g])
      for (g in 1:NGender) {
        tmp <- 0.0
        for (a in 1:NAge) {
          if (RecAge == 0)
            age <- a-1
          else if (RecAge == 1)
            age <- a
          EquilibriumFisherySelectivityAtAge[p,FleetArea[v],g,a]  <- Selectivity(age,parameters)
        }
# Rescale to set maximum selectivity at age to be 1      
        tmp <- max(EquilibriumFisherySelectivityAtAge[p,FleetArea[v],g,])
        for (a in 1:NAge) {
          EquilibriumFisherySelectivityAtAge[p,FleetArea[v],g,a] <- EquilibriumFisherySelectivityAtAge[p,FleetArea[v],g,a]/tmp
          }
      }
    }

