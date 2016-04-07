library(shiny)
library(deSolve)    # Needed for ode function
library(ggplot2)    # Plotting
library(reshape2)   # Manipulating data
library(FME)
times = seq(0,150,by = 0.0625)

# Compartments
Susceptible=99999; Exposed=0;
Infected_Asymptomatic=0; Infected_Mild=1;
Infected_Critical=0; Critical_in_Hospital=0; Critical_in_Community=0;
Death=0; Recovered=0; 

#  Parameters
R0=2; Recovery_Delay=2; Total_Population=100000;
Total_Capacity=500;PropA=.3;

# Vector of paramaters
par <- c(R0=R0,Recovery_Delay=Recovery_Delay, Total_Population=Total_Population,
         PropM=PropM,Total_Capacity=Total_Capacity)

# Model
SurgeCapacityModel <- function(par) {
  # This is the callback function for ode
  derivs <- function(t, state, par) {
    with(as.list(c(state, par)),{
      # Proportion Asymptomatic
      PropM=.6;
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
BaseModel <<- SurgeCapacityModel(par)
ui <- fluidPage(
  tabsetPanel(
    tabPanel("Base Model", 
             sidebarPanel(
              # conditionalPanel('input.dataset == "BaseModel"',
                 checkboxGroupInput('display', 'Curves in Base Model to show:',
                                    names(BaseModel[2:length(names(BaseModel))]), selected=names(BaseModel[2:length(names(BaseModel))]))
               ),
             mainPanel(fluidRow(column(10, offset=0, plotOutput(outputId="BaseModel"))))),
    tabPanel("Deaths",
             fluidRow(column(2, offset=2,sliderInput(inputId = "totalCapacitySingle",
                                                     label="Total Capacity Single Value",
                                                     value=500, min = 20, max = 1000)),
                      column(2, offset=4,sliderInput(inputId = "totalCapacitySens",
                                                     label="Total Capacity Sensitivity",
                                                     dragRange = c(20,1000), min=20, max=1000, value=c(100,110)))),
             
             fluidRow(column(6, plotOutput(outputId = "CapacityPlotSingle")),
                      column(6, plotOutput(outputId = "DeathSensitivity"))),
             
             fluidRow(column(4, offset=8, numericInput(inputId="numRun", label="Number of Runs", value=1, min=1),
                             actionButton(inputId="runCapacitySens", label="Run Sensitivity")))
    ),
    # have variable to show current propC too!!!
    tabPanel("Total Infected",
             fluidRow(column(2, offset=2, sliderInput(inputId = "propM",
                                                      label="Proportion Mild Single Value",
                                                      value = .6, min = 0, max = .7)),
                      column(2, offset=4,sliderInput(inputId = "totalProportionSens",
                                                     label="Proportion Mild Sensitivity",
                                                     dragRange = c(0,.7), min=0, max=.7, value=c(.2,.6)))),
             fluidRow(column(2, offset=2, p("Proportion Critical")),
                      column(1, textOutput("PropC")),
                      column(3,  offset=2,textOutput("PropTitle")),
                      column(2, textOutput("rangeProp"))),
             
             fluidRow(column(6, plotOutput(outputId = "ProportionPlotSingle")),
                      column(6, plotOutput(outputId = "ProportionSensitivity"))),
             
             fluidRow(column(4, offset=8, numericInput(inputId="numRun2", label="Number of Runs", value=1, min=1),
                             actionButton(inputId="runProportionSens", label="Run Sensitivity")))
             
             #  plotOutput()
             
    )
  )
)

server <- function(input, output)
{
  
  
  output$BaseModel <- renderPlot({plotData <- BaseModel[c(input$display)]
  plotData$time <- BaseModel$time
    BaseMelt <- melt(plotData, id.vars="time")
  ggplot(BaseMelt, aes(time, value, col=variable))+geom_line()
  })
  
  output$CapacityPlotSingle <- renderPlot({
    par <- c(R0=R0,Recovery_Delay=Recovery_Delay, Total_Population=Total_Population,
             PropM=PropM,Total_Capacity=input$totalCapacitySingle)
    NewModel1 <- SurgeCapacityModel(par)
    ggplot(NewModel1, aes(time, Death))+geom_line()
  })
  
  output$DeathSensitivity <- renderPlot({sensPlot()})
  
  
  sensPlot <- eventReactive(input$runCapacitySens, { 
    par <- c(R0=R0,Recovery_Delay=Recovery_Delay, Total_Population=Total_Population,
             PropM=PropM,Total_Capacity=Total_Capacity)
    capacityRanges <- data.frame(min=c(input$totalCapacitySens[1]), max=c(input$totalCapacitySens[2]))
    #capacityRanges <- data.frame(min=c(40), max=c(50))
    
    rownames(capacityRanges) <- c("Total_Capacity")
    #Sensitivity function for Total Capacity on Deaths
    TotalCapacityImpact <- sensRange(SurgeCapacityModel, # Name of model to call
                                     parms=par, # Parameters used in the model
                                     dist="latin", # Distribution of values for sensitivity
                                     c("Death"), # Compartment we want results for 
                                     parRange=capacityRanges, # Ranges of parameter we want to test
                                     num=input$numRun) # Number of sensitivity runs
    
    # Get only values for Death compartment
    results <- TotalCapacityImpact[grepl("Death", names(TotalCapacityImpact))]
    # Transpose and use as data frame for plotting
    results <- as.data.frame(t(results))
    # Add time values for plotting
    results$time <- times
    # Reshape into tidy data format on times
    resultsM <- melt(results, id.vars="time")
    # Add run number
    resultsM$run <- rep(1:input$numRun, each=length(times))
    # Plot results for sensitivity, (x-axis, y-axis, color)
    BaseModel <- SurgeCapacityModel(par)
    ggplot(resultsM, aes(time, value, col=run))+geom_line(aes(group=run))+
      # Add base model result for comparison
      geom_line(aes(BaseModel$time, BaseModel$Death), col="red") })
  
  output$PropTitle <- renderText("Proportion Critical Range")
  
  output$ProportionPlotSingle <- renderPlot({
    par <- c(R0=R0,Recovery_Delay=Recovery_Delay, Total_Population=Total_Population,
             PropM=input$propM,Total_Capacity=Total_Capacity)
    NewModel2 <- SurgeCapacityModel(par)
    ggplot(NewModel2, aes(time, Total_Infected))+geom_line()
  })
  
  #totalProportionSens
  output$PropC <- renderText(1-(.3+input$propM))
  output$rangeProp  <- renderText({paste((1-(.3+input$totalProportionSens[2])), (1-(.3+input$totalProportionSens[1])), sep="-")})
  
  
  output$ProportionSensitivity<- renderPlot({sensPlot2()})
  
  sensPlot2 <- eventReactive(input$runProportionSens, {
    
    par <- c(R0=R0,Recovery_Delay=Recovery_Delay, Total_Population=Total_Population,
             PropM=PropM,Total_Capacity=Total_Capacity)
    
    # Ranges for second sensitivity test
    proportionRanges <- data.frame(min=c(input$totalProportionSens[1]), max=c(input$totalProportionSens[2]))
    rownames(proportionRanges) <- c("PropM")
    # Sensitivity function for Proportion Mild (and Critical by proxy) on Total number of people Infected
    ProportionImpact <- sensRange(SurgeCapacityModel, # Name of model to call
                                  parms=par, # Parameters used in the model
                                  dist="latin", # Distribution of values for sensitivity
                                  c("Total_Infected"), # Compartment we want results for
                                  parRange=proportionRanges, # Ranges of parameter we want to test
                                  num=input$numRun2) # Number of sensitivity runs
    
    # Get only values for Total Infected parameter
    results2 <- ProportionImpact[grepl("Total_Infected", names(ProportionImpact))]
    # Transpose and use as data frame for plotting
    results2 <- as.data.frame(t(results2))
    # Add time values for plotting
    results2$time <- times
    # Reshape into tidy data format on times
    resultsM2 <- melt(results2, id.vars="time")
    # Add run number
    resultsM2$run <- rep(1:input$numRun2, each=length(times))
    # Plot results for sensitivity (x-axis, y-axis, color)
    ggplot(resultsM2, aes(time, value, col=run))+geom_line(aes(group=run))+
      # Add base model result for comparison
      geom_line(aes(BaseModel$time, BaseModel$Total_Infected), col="red")
  })
  
  
  
  
  
}

shinyApp(ui=ui, server=server)
