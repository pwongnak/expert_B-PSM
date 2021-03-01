##########################################################################################
#Bayesian parametric survival model
#Shiny Web-application for exert opinion interview
############################################################################################
library(shiny)
library(ggplot2)
library(shinythemes)
library(ggpubr)
#############################################################################################
#A function for opimize a log-normal distribution given mean, and upper and lower values within an interval
f_lognormal <- function(data, params){
  mean <- data[1] #T_bar_m
  q1 <- data[2] #T_bar_l
  q3 <- data[3] #T_bar_h
  ci <- data[4] #C
  
  sdlog <- params[1]
  meanlog <- log(mean) - ((sdlog^2)/2)
  
  cumeanlog1 <- plnorm(q1, meanlog = meanlog, sdlog  = sdlog)
  cumeanlog3 <- plnorm(q3, meanlog = meanlog, sdlog = sdlog)
  area <- cumeanlog3 - cumeanlog1
  abs(area - ci)
}

#User input####################################################################################
ui <- fluidPage(theme = shinytheme("cerulean"),
                # Application title
                titlePanel( h1(div(HTML("Expert opinions :<br> A questionaire on an average survival time of female adult <em> Ixodes ricinus </em> exposed to constant temperature and relative humidity</br>" )))),
                textInput(inputId = "ExpertName", label = "Name and Last name : ", value = "Annonymous"),
                #Condition1---------------------------------------------------------------------------
                titlePanel(h3("Condition 1: Low humidity (dry) and low temperature (cold)")),
                sidebarLayout(
                  sidebarPanel(
                    numericInput(inputId = "Tp1", label = "Temperature (°C)", value = 5),
                    numericInput(inputId = "U1", label = "Relative humidity [0-1]", value = 0.1, min = 0, max = 1, step = 0.1),
                    numericInput(inputId = "mean1", label = "Average survival time", value = ""),
                    numericInput(inputId = "min1", label = "Lower value", value = ""),
                    numericInput(inputId = "max1", label = "Higher value", value = ""),
                    numericInput(inputId = "CI1", label = "Confidence level [0-1]", value = 0.95, min = 0, max = 1, step = 0.1),
                    width = 3
                  ),
                  mainPanel(
                    plotOutput("plot1"),
                    h4("Parameters"),
                    tableOutput("dat1"),
                    width = 6
                  )
                ),
                #Condition2---------------------------------------------------------------------------
                titlePanel(h3("Condition 2: Low humidity (dry) and High temperature (warm)")),
                sidebarLayout(
                  sidebarPanel(
                    numericInput(inputId = "Tp2", label = "Temperature (°C)", value = 25),
                    numericInput(inputId = "U2", label = "Relative humidity [0-1]", value = 0.1, min = 0, max = 1, step = 0.1),
                    numericInput(inputId = "mean2", label = "Average survival time", value = ""),
                    numericInput(inputId = "min2", label = "Lower value", value = ""),
                    numericInput(inputId = "max2", label = "Higher value", value = ""),
                    numericInput(inputId = "CI2", label = "Confidence level [0-1]", value = 0.95, min = 0, max = 1, step = 0.1),
                    width = 3
                  ),
                  mainPanel(
                    plotOutput("plot2"),
                    h4("Parameters"),
                    tableOutput("dat2"),
                    width = 6
                  )
                ),
                
                #Condition 3-----------------------------------------------------------------------
                titlePanel(h3("Condition 3: Humidité haute (=humide) and température basse (=froid)")),
                sidebarLayout(
                  sidebarPanel(
                    numericInput(inputId = "Tp3", label = "Temperature (°C)", value = 5),
                    numericInput(inputId = "U3", label = "Humidité Relative [0-1]", value = 0.95, min = 0, max = 1, step = 0.1),
                    numericInput(inputId = "mean3", label = "Moyenne", value = ""),
                    numericInput(inputId = "min3", label = "Fourchette basse", value = ""),
                    numericInput(inputId = "max3", label = "Fourchette haute", value = ""),
                    numericInput(inputId = "CI3", label = "Niveau de Confiance [0-1]", value = 0.95, min = 0, max = 1, step = 0.1),
                    width = 3
                  ),
                  mainPanel(
                    plotOutput("plot3"),
                    h4("Paramètres"),
                    tableOutput("dat3"),
                    width = 6
                  )
                ),
                
                #Condition 4-----------------------------------------------------------------------
                titlePanel(h3("Condition 4: Humidité haute (=humide) and température haute (=chaud)")),
                sidebarLayout(
                  sidebarPanel(
                    numericInput(inputId = "Tp4", label = "Temperature (°C)", value = 25),
                    numericInput(inputId = "U4", label = "Humidité Relative [0-1]", value = 0.95, min = 0, max = 1, step = 0.1),
                    numericInput(inputId = "mean4", label = "Moyenne", value = ""),
                    numericInput(inputId = "min4", label = "Fourchette basse", value = ""),
                    numericInput(inputId = "max4", label = "Fourchette haute", value = ""),
                    numericInput(inputId = "CI4", label = "Niveau de Confiance [0-1]", value = 0.95, min = 0, max = 1, step = 0.1),
                    width = 3
                  ),
                  mainPanel(
                    plotOutput("plot4"),
                    h4("Paramètres"),
                    tableOutput("dat4"),
                    width = 6
                  ),
                ),
                
                titlePanel(h3("Summary")),
                plotOutput("AllPlot", width = "800px", height = "600px"),
                
                titlePanel("Enregistrer vos réponses"),
                actionButton(inputId = "Submit", label = "Soumettre")
)

#Server#######################################################################################################################
server <- function(input, output) {
  #Condition 1---------------------------------------------------------------------------------------------------------------------------
  #Data
  dataT1 <- reactive({
    mean <- input$mean1
    min <- input$min1
    max <- input$max1
    ci <- input$CI1
    Tp <- input$Tp1
    U <- input$U1
    
    data.frame(Tp, U, mean, min, max, ci)
  })
  #response to update button #use for 1. plot, 2. save
  paramT1 <- reactive({
    inputData <- as.data.frame(dataT1())
    min <- inputData$min
    max <- inputData$max
    mean <- inputData$mean
    prob <- inputData$ci
    
    data <- c(mean, min, max, prob)
    fit <- optim(par = 1, fn = f_lognormal, data = data, method = "Brent", lower = 0.0001, upper = 1)
    sdlog  <- fit$par
    meanlog <- log(mean - ((sdlog^2)/2))
    
    data.frame(meanlog, sdlog)
  })
  #output1
  output$dat1 <- renderTable({
    paramT1()
  })
  
  G1 <- reactive({
    param <- as.data.frame(paramT1())
    meanlog <- param$meanlog
    sdlog <- param$sdlog
    
    data <- as.data.frame(dataT1())
    mean <- data$mean
    min <- data$min
    max <- data$max
    ci <- data$ci
    
    x <- seq(0, 2*max, 0.1)
    density <- as.vector(dlnorm(x, meanlog, sdlog))
    df <- as.data.frame(cbind(x, density))
    plot <- ggplot(df, aes(x = x, y = density)) +
      geom_area(data = subset(df, x >= min & x <= max), aes(x=x,y=density), fill = 'lightblue') +
      geom_line(size = 1) +
      geom_segment(x = mean, y = 0, xend = mean , yend = df[df$x == mean,2], col = "red", lwd = 0.5) +
      labs(title = paste0("Temperature = ", input$Tp1, " °C, Relative humidity = ", input$U1),  x =   expression(paste("Average survival time (days), ", italic(bar("T")))), y  = "Density") +
      theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    plot
  })
  
  #plot
  output$plot1 <- renderPlot({
    G1()
  })
  
  #Condition 2---------------------------------------------------------------------------------------------------------------------------
  #Data
  dataT2 <- reactive({
    mean <- input$mean2
    min <- input$min2
    max <- input$max2
    ci <- input$CI2
    Tp <- input$Tp2
    U <- input$U2
    
    data.frame(Tp, U, mean, min, max, ci)
  })
  
  #response to update button #use for 1. plot, 2. save
  paramT2 <- reactive({
    inputData <- as.data.frame(dataT2())
    min <- inputData$min
    max <- inputData$max
    mean <- inputData$mean
    prob <- inputData$ci
    
    data <- c(mean, min, max, prob)
    fit <- optim(par = 1, fn = f_lognormal, data = data, method = "Brent", lower = 0.0001, upper = 1)
    sdlog  <- fit$par
    meanlog <- log(mean - ((sdlog^2)/2))
    
    data.frame(meanlog, sdlog)
  })
  
  #output2
  output$dat2 <- renderTable({
    paramT2()
  })
  
  G2 <- reactive({
    param <- as.data.frame(paramT2())
    meanlog <- param$meanlog
    sdlog <- param$sdlog
    
    data <- as.data.frame(dataT2())
    mean <- data$mean
    min <- data$min
    max <- data$max
    ci <- data$ci
    
    x <- seq(0, 2*max, 0.1)
    density <- as.vector(dlnorm(x, meanlog, sdlog))
    df <- as.data.frame(cbind(x, density))
    plot <- ggplot(df, aes(x = x, y = density)) +
      geom_area(data = subset(df, x >= min & x <= max), aes(x=x,y=density), fill = 'lightblue') +
      geom_line(size = 1) +
      geom_segment(x = mean, y = 0, xend = mean , yend = df[df$x == mean,2], col = "red", lwd = 0.5) +
      labs(title = paste0("Temperature = ", input$Tp2, " °C, Relative humidity = ", input$U2),  x = "Mean survival time (days)", y  = "Density") +
      theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    plot
  })
  
  
  #plot
  output$plot2 <- renderPlot({
    G2()
  })
  
  #Condition 3---------------------------------------------------------------------------------------------------------------------------
  #Data
  dataT3 <- reactive({
    mean <- input$mean3
    min <- input$min3
    max <- input$max3
    ci <- input$CI3
    Tp <- input$Tp3
    U <- input$U3
    
    data.frame(Tp, U, mean, min, max, ci)
  })
  
  #response to update button #use for 1. plot, 2. save
  paramT3 <- reactive({
    inputData <- as.data.frame(dataT3())
    min <- inputData$min
    max <- inputData$max
    mean <- inputData$mean
    prob <- inputData$ci
    
    data <- c(mean, min, max, prob)
    fit <- optim(par = 1, fn = f_lognormal, data = data, method = "Brent", lower = 0.0001, upper = 1)
    sdlog  <- fit$par
    meanlog <- log(mean - ((sdlog^2)/2))
    
    data.frame(meanlog, sdlog)
  })
  
  #output3
  output$dat3 <- renderTable({
    paramT3()
  })
  
  G3 <- reactive({
    param <- as.data.frame(paramT3())
    meanlog <- param$meanlog
    sdlog <- param$sdlog
    
    data <- as.data.frame(dataT3())
    mean <- data$mean
    min <- data$min
    max <- data$max
    ci <- data$ci
    
    x <- seq(0, 2*max, 0.1)
    density <- as.vector(dlnorm(x, meanlog, sdlog))
    df <- as.data.frame(cbind(x, density))
    plot <- ggplot(df, aes(x = x, y = density)) +
      geom_area(data = subset(df, x >= min & x <= max), aes(x=x,y=density), fill = 'lightblue') +
      geom_line(size = 1) +
      geom_segment(x = mean, y = 0, xend = mean , yend = df[df$x == mean,2], col = "red", lwd = 0.5) +
      labs(title = paste0("Temperature = ", input$Tp3, " °C, Relative humidity = ", input$U3),  x = "Mean survival time (days)", y  = "Density") +
      theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    plot
  })
  
  
  #plot
  output$plot3 <- renderPlot({
    G3()
  }) 
  
  #Condition 4---------------------------------------------------------------------------------------------------------------------------
  #Data
  dataT4 <- reactive({
    mean <- input$mean4
    min <- input$min4
    max <- input$max4
    ci <- input$CI4
    Tp <- input$Tp4
    U <- input$U4
    
    data.frame(Tp, U, mean, min, max, ci)
  })
  
  #response to update button #use for 1. plot, 2. save
  paramT4 <- reactive({
    inputData <- as.data.frame(dataT4())
    min <- inputData$min
    max <- inputData$max
    mean <- inputData$mean
    prob <- inputData$ci
    
    data <- c(mean, min, max, prob)
    fit <- optim(par = 1, fn = f_lognormal, data = data, method = "Brent", lower = 0.0001, upper = 1)
    sdlog  <- fit$par
    meanlog <- log(mean - ((sdlog^2)/2))
    
    data.frame(meanlog, sdlog)
  })
  
  #output4
  output$dat4 <- renderTable({
    paramT4()
  })
  
  G4 <- reactive({
    param <- as.data.frame(paramT4())
    meanlog <- param$meanlog
    sdlog <- param$sdlog
    
    data <- as.data.frame(dataT4())
    mean <- data$mean
    min <- data$min
    max <- data$max
    ci <- data$ci
    
    x <- seq(0, 2*max, 0.1)
    density <- as.vector(dlnorm(x, meanlog, sdlog))
    df <- as.data.frame(cbind(x, density))
    plot <- ggplot(df, aes(x = x, y = density)) +
      geom_area(data = subset(df, x >= min & x <= max), aes(x=x,y=density), fill = 'lightblue') +
      geom_line(size = 1) +
      geom_segment(x = mean, y = 0, xend = mean , yend = df[df$x == mean,2], col = "red", lwd = 0.5) +
      labs(title = paste0("Temperature = ", input$Tp4, " °C, Relative humidity = ", input$U4),  x = "Mean survival time (days)", y  = "Density") +
      theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    plot
  })
  
  
  
  #plot
  output$plot4 <- renderPlot({
    G4()
  })
  
  
  output$AllPlot <- renderPlot({
    P1 <- G1() 
    P2 <- G2()
    P3 <- G3()
    P4 <- G4()
    
    ggarrange(P1, P2, P3, P4, ncol = 2, nrow = 2)
  })
  
  
  
  #Save the answer
  observeEvent(input$Submit, {
    dataT1 <- as.data.frame(dataT1())
    dataT2 <- as.data.frame(dataT2())
    dataT3 <- as.data.frame(dataT3())
    dataT4 <- as.data.frame(dataT4())
    
    paramT1 <- as.data.frame(paramT1())
    paramT2 <- as.data.frame(paramT2())
    paramT3 <- as.data.frame(paramT3())
    paramT4 <- as.data.frame(paramT4())
    
    data <- rbind(dataT1, dataT2, dataT3, dataT4)
    param <- rbind(paramT1, paramT2, paramT3, paramT4)
    Time <- rep(Sys.time(), 4)
    Expert <- rep(input$ExpertName,4)
    allData <- data.frame(Expert,data, param, Time)
    colnames(allData) <- c("Expert", "Tp", "U", "T_mean", "T_min", "T_max", "Ci", "meanLog", "sdLog", "Time")
    write.csv(allData, file = paste0("ExpertData_", input$ExpertName,".csv"), row.names = F)
  })
  
}
#Application execution#########################################################################################
shinyApp(ui = ui, server = server)

