# Surge Capacity Model with Sensitivity

# Martina Curran 2016
# Information Technology Department, National University of Ireland Galway

# Libraries needed
library(deSolve)    # Needed for ode function
library(ggplot2)    # Plotting
library(reshape2)   # Manipulating data
library(FME)        # Needed for sensitivity function

# Time sequence for model
times = seq(0,150,by = 0.0625)

# Compartments
Susceptible=99999; Exposed=0;
Infected_Asymptomatic=0; Infected_Mild=1;
Infected_Critical=0; Critical_in_Hospital=0; Critical_in_Community=0;
Death=0; Recovered=0; 

#  Parameters
R0=2; Recovery_Delay=2; Total_Population=100000;
Total_Capacity=500;PropM=.6;

# Vector of paramaters
par <- c(R0=R0,Recovery_Delay=Recovery_Delay, Total_Population=Total_Population,
         PropM=PropM,Total_Capacity=Total_Capacity)

# Model
SurgeCapacityModel <- function(par) {
  # This is the callback function for ode
  derivs <- function(t, state, par) {
    with(as.list(c(state, par)),{
      # Proportion Asymptomatic
      PropA=.3;
      # Proportion Critical
      PropC= 1-(PropA+PropM);  
      
      # Effective contact rate
      Beta <- R0/(Recovery_Delay*Total_Population)
      # Total number of people infected in the population
      Total_Infected <- Infected_Mild+Infected_Asymptomatic+
        Infected_Critical+Critical_in_Hospital+
        Critical_in_Community
      
      # Force of infection 
      Lambda <- Beta*Total_Infected
      # Available Capacity of Hospital
      Available_Capacity <- Total_Capacity-Critical_in_Hospital
      # Pressure put on the system from Critical Infections based on the Total Capacity of Hospital
      System_Pressure <- (Infected_Critical+Critical_in_Community+
                            Critical_in_Hospital)/Total_Capacity
      
      IR <- Lambda*Susceptible # Infection Rate
      PC <- PropC*(Exposed/2)  # Proportion Critical Rate
      PM <- PropM*(Exposed/2)  # Proportion Mild Rate
      PA <- PropA*(Exposed/2)  # Proportion Asymptomatic Rate
      CIH <- min(Available_Capacity, Infected_Critical) # Critical in Hospital Rate
      CIC <- Infected_Critical-(min(Available_Capacity, Infected_Critical)) # Critical in Community Rate
      RRM <- Infected_Mild/2 # Recovery Rate of Mild Infection
      RRA <- Infected_Asymptomatic/2 # Recovery Rate of Asymptomatic Infection
      RRC <- Critical_in_Community/2 # Recovery Rate of Critical Infection 
      DRC <- Critical_in_Community *.07 # Death Rate of Critical Infection in Community
      DRH <- Critical_in_Hospital * .01 # Death Rate of Critical Infection in Hospital
      RRH <- Critical_in_Hospital/2 # Recovery Rate of Critical Infection in Hospital 
      
      dS_dt <- -IR #Rate people leave Susceptible
      dE_dt <- IR-PA-PC-PM #Rate people enter and leave Exposed 
      dIA_dt <- PA-RRA # Rate people enter and leave Infected Asymptomatic
      dIM_dt <- PM-RRM # Rate people enter and leave Infected Mild
      dIC_dt <- PC-CIC-CIH # Rate people enter and leave Infected Critical
      dCIH_dt <- CIH-DRH-RRH # Rate people enter and leave Critical in Hospital
      dCIC_dt <- CIC-DRC-RRC # Rate people enter and leave Critical in Community
      dD_dt <- DRC+DRH # Rate of Deaths 
      dR_dt <- RRA+RRC+RRH+RRM # Rate of Recoveries
      
      # Details to return:
      return (list(
        # (Compartment values for each iteration), extra paramaters wanted/needed
        c(dS_dt,dE_dt,dIA_dt,dIM_dt,dIC_dt,dCIH_dt,dCIC_dt,dD_dt,dR_dt),Total_Capacity=Total_Capacity,Total_Infected=Total_Infected)
        )
    })
  }
  
  # Input of compartments for model. Returning list from ode needs to be in the same order
  state <-c(Susceptible=Susceptible, Exposed=Exposed, Infected_Asymptomatic=Infected_Asymptomatic,
            Infected_Mild=Infected_Mild, Infected_Critical=Infected_Critical,
            Critical_in_Hospital=Critical_in_Hospital, Critical_in_Community=Critical_in_Community, Death=Death, Recovered=Recovered)
  
  # Call the ode function, passing in all paramaters, and return values
  return(as.data.frame(ode(y = state, times, func = derivs, parms = par, method = "euler")))
  
}

# Call base model
BaseModel <- SurgeCapacityModel(par)

# Ranges for first sensitivity test
capacityRanges <- data.frame(min=c(20), max=c(1000))
rownames(capacityRanges) <- c("Total_Capacity")
# Sensitivity function for Total Capacity on Deaths
TotalCapacityImpact <- sensRange(SurgeCapacityModel, # Name of model to call
                                 parms=par, # Parameters used in the model
                                 dist="latin", # Distribution of values for sensitivity
                                 c("Death"), # Compartment we want results for 
                                 parRange=capacityRanges, # Ranges of parameter we want to test
                                 num=200) # Number of sensitivity runs

# Get only values for Death compartment
results <- TotalCapacityImpact[grepl("Death", names(TotalCapacityImpact))]
# Transpose and use as data frame for plotting
results <- as.data.frame(t(results))
# Add time values for plotting
results$time <- times
# Reshape into tidy data format on times
resultsM <- melt(results, id.vars="time")
# Add run number
resultsM$run <- rep(1:200, each=length(times))
# Plot results for sensitivity, (x-axis, y-axis, color)
ggplot(resultsM, aes(time, value, col=run))+geom_line()+
  # Add base model result for comparison
  geom_line(aes(BaseModel$time, BaseModel$Death), col="black")

# Ranges for second sensitivity test
proportionRanges <- data.frame(min=c(0.1), max=c(0.7))
rownames(proportionRanges) <- c("PropM")
# Sensitivity function for Proportion Mild (and Critical by proxy) on Total number of people Infected
ProportionImpact <- sensRange(SurgeCapacityModel, # Name of model to call
                              parms=par, # Parameters used in the model
                              dist="latin", # Distribution of values for sensitivity
                              c("Total_Infected"), # Compartment we want results for
                              parRange=proportionRanges, # Ranges of parameter we want to test
                              num=200) # Number of sensitivity runs

# Get only values for Total Infected parameter
results2 <- ProportionImpact[grepl("Total_Infected", names(ProportionImpact))]
# Transpose and use as data frame for plotting
results2 <- as.data.frame(t(results2))
# Add time values for plotting
results2$time <- times
# Reshape into tidy data format on times
resultsM2 <- melt(results2, id.vars="time")
# Add run number
resultsM2$run <- rep(1:200, each=length(times))
# Plot results for sensitivity (x-axis, y-axis, color)
ggplot(resultsM2, aes(time, value, col=run))+geom_line()+
  # Add base model result for comparison
  geom_line(aes(BaseModel$time, BaseModel$Total_Infected), col="black")

