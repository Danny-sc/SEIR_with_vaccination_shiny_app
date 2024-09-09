# 
# Shiny app for SIR model with vaccinations
#
# Created by Claus Ekstrøm 2019
# @ClausEkstrom
#

library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")
library("markdown")

## Create an SIR function
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # dS <- -beta * S * I
    # dI <-  beta * S * I - gamma * I
    # dR <-                 gamma * I
    # dV <- 0
    #return(list(c(dS, dI, dR, dV)))
    dS <- -beta * S * I
    dE <-  beta * S * I - alpha *E
    dI <-                 alpha * E - gamma * I
    dR <-                             gamma * I
    #dV <- 0
    #return(list(c(dS, dE, dI, dR, dV)))
    return(list(c(dS, dE, dI, dR)))
  })
}


#
# Define UI 
#

ui <- sidebarLayout(
#  dashboardHeader(disable = TRUE),
#  dashboardSidebar(
    
#    sidebarLayout(
      
      sidebarPanel(
        
        fluidRow(
          column(width=10,  
    
    sliderInput("popsize",
                "Population size (millions):",
                min = 1, max = 300, value = 6
    ),
    sliderInput("connum",
                "Basic reproductive number (R0, # persons):",
                min = .5, max = 20, value = 5
    ),
    sliderInput("pinf",
                "Number of infectious individuals at the start of the outbreak:",
                min = 1, max = 50, value = 2
    ),
    sliderInput("pvac",
                "Proportion vaccinated (%):",
                min = 0, max = 100, value = 0
    ),
    sliderInput("vaceff",
                "Vaccine efficacy (%):",
                min = 0, max = 100, value = 0
    ),
    sliderInput("infper",
                "Infectious period (days):",
                min = 1, max = 30, value = 7
    ),
    sliderInput("expper",
                "Latent period (days):",
                min = 1, max = 15, value = 5
    ),
    sliderInput("timeframe",
                "Time frame (days):",
                min = 1, max = 400, value = 200
    )
          ))),
    
  # ),
  # dashboardBody(
  #   tags$head(tags$style(HTML('
  #                             /* body */
  #                             .content-wrapper, .right-side {
  #                             background-color: #fffff8;
  #                             }                              
  #                             '))),
    
    mainPanel(
      navbarPage("Output:",
    tabPanel("Spread", br(),
               tags$head(tags$style(HTML('
                                         /* body */
                                         .content-wrapper, .right-side {
                                         background-color: #fffff8;
                                         }
                                         '))),
             fluidRow(column(12,
                             h2("Simulate the course of an epidemic (parameters can be varied with sliders, left panel)"),
                             #plotOutput("plot4", height=200),
                             #includeMarkdown("SEIR.Rmd"),
                             #h3("Equations"),
                             #br(),
                             # h2("Output"),
                             # h3("Rate parameters of dynamic model"),
                             # p(HTML("These parameters can be changed using the sliders in the other tabs. The values in this table represent the current values chosen via the sliders. Note that the transmission rates chosen by the sliders are always scaled by \\(N\\), so that \\(\\beta*N\\) is constant as \\(N\\) changes.")),
                             # tableOutput("ParameterTable"),br(),
                             # h3("Ratios of cases during early growth phase"),
                             # p(HTML("These values are calculated based on the current model parameters")),
                             # tableOutput("RatioTable"),br(),
             )),         
    fluidRow(plotOutput("distPlot")),
    br(),
    fluidRow(
      # Dynamic valueBoxes
      valueBoxOutput("progressBox", width = 6),
      valueBoxOutput("approvalBox", width = 6),
      valueBoxOutput("BRRBox", width = 6),
      valueBoxOutput("HIBox", width = 6)
    ),
    br(),
    br()
    ),
    tabPanel("Vaccine impact", br(),
             tags$head(tags$style(HTML('
                                         /* body */
                                         .content-wrapper, .right-side {
                                         background-color: #fffff8;
                                         }
                                         '))),
             fluidRow(column(12,
                             h2("Explore the impact of a vaccination campaign on the number of infectious individuals"),
                             #plotOutput("plot4", height=200),
                             #includeMarkdown("SEIR.Rmd"),
                             #h3("Equations"),
                             #br(),
                             # h2("Output"),
                             # h3("Rate parameters of dynamic model"),
                             # p(HTML("These parameters can be changed using the sliders in the other tabs. The values in this table represent the current values chosen via the sliders. Note that the transmission rates chosen by the sliders are always scaled by \\(N\\), so that \\(\\beta*N\\) is constant as \\(N\\) changes.")),
                             # tableOutput("ParameterTable"),br(),
                             # h3("Ratios of cases during early growth phase"),
                             # p(HTML("These values are calculated based on the current model parameters")),
                             # tableOutput("RatioTable"),br(),
             )),  
             fluidRow(plotOutput("distPlotNoInt")),
              br(),
              fluidRow(
                # Dynamic valueBoxes
                valueBoxOutput("progressBox", width = 6),
                valueBoxOutput("progressBoxNoInt", width = 6),
             #   valueBoxOutput("BRRBox", width = 6),
             #   valueBoxOutput("HIBox", width = 6)
              ),
             br(),
             br()
    ),
    tabPanel("Model", br(),
             fluidRow(column(12,
                             withMathJax(),
                             h2("Model Description"),
                             plotOutput("plot4", height=200),
                             includeMarkdown("SEIR.Rmd"),
                             #h3("Equations"),
                             br(),
                             # h2("Output"),
                             # h3("Rate parameters of dynamic model"),
                             # p(HTML("These parameters can be changed using the sliders in the other tabs. The values in this table represent the current values chosen via the sliders. Note that the transmission rates chosen by the sliders are always scaled by \\(N\\), so that \\(\\beta*N\\) is constant as \\(N\\) changes.")),
                             # tableOutput("ParameterTable"),br(),
                             # h3("Ratios of cases during early growth phase"),
                             # p(HTML("These values are calculated based on the current model parameters")),
                             # tableOutput("RatioTable"),br(),
             )))
  )
)

)

#
# Define server 
#
server <- function(input, output) {
  # Create reactive input
  dataInput <- reactive({
    init       <-
      c(
        S = 1 - input$pinf / (input$popsize*1000000) - input$pvac / 100 * input$vaceff / 100,
        E = 0,
        I = input$pinf /  (input$popsize*1000000),
        #R = 0,
        #V = input$pvac / 100 * input$vaceff / 100
        R = 0
        #V = input$pvac / 100 * input$vaceff / 100
      )
    ## beta: infection parameter; gamma: recovery parameter
    parameters <-
      c(beta = input$connum * 1 / input$infper,
        alpha = 1 / input$expper,
        # * (1 - input$pvac/100*input$vaceff/100),
        gamma = 1 / input$infper)
    ## Time frame
    times      <- seq(0, input$timeframe, by = .2)
    
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    out <- ode(
      y = init,
      times = times,
      func = sir,
      parms = parameters
    )   
    #    out
    as.data.frame(out)
  })
  
  dataInputNoInt <- reactive({
    init       <-
      c(
        S = 1 - input$pinf / (input$popsize*1000000),
        E = 0,
        I = input$pinf /  (input$popsize*1000000),
        #R = 0,
        #V = 0
        R = 0
      )
    ## beta: infection parameter; gamma: recovery parameter
    parameters <-
      c(beta = input$connum * 1 / input$infper,
        alpha = 1 / input$expper,
        # * (1 - input$pvac/100*input$vaceff/100),
        gamma = 1 / input$infper)
    ## Time frame
    times      <- seq(0, input$timeframe, by = .2)
    
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    out <- ode(
      y = init,
      times = times,
      func = sir,
      parms = parameters
    )   
    #    out
    as.data.frame(out)
  })
  
  
  output$distPlot <- renderPlot({
    out <-
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        key2 = recode(
          key,
          S = "Susceptible (S)",
          #E = "Exposed (E)",
          I = "Infectious (I)",
          #R = "Recovered (R)",
          #V = "Vaccinated (V)"
          R = "Recovered (R)"
          #V = "Vaccinated (V)"
        ),
        keyleft = recode(
          key,
          S = "Susceptible (S)",
          E = "Exposed (E)",
          I = "",
          #R = "",
          #V = "Vaccinated (V)"
          R = ""
          #V = "Vaccinated (V)"
        ),
        keyright = recode(
          key,
          S = "",
          E = "",
          I = "Infectious (I)",
          #R = "Recovered (R)",
          #V = ""
          R = "Recovered (R)"
          #V = ""
        )
      )
    
    ggplot(data = out,
           aes(
             x = time,
             y = value,
             group = key2,
             col = key2,
             label = key2,
             data_id = id
           )) + # ylim(0, 1) +
      ylab("Proportion of full population") + xlab("Time (days)") +
      geom_line(size = 2) +
      geom_text_repel(
        data = subset(out, time == max(time)),
        aes(label = keyright),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      geom_text_repel(
        data = subset(out, time == min(time)),
        aes(label = keyleft),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 0,
        direction = "y"
      ) +
      theme(legend.position = "none") +
      scale_colour_manual(values = c("orange", "green4", "black", "blue")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      theme(
        axis.text.x=element_text(size=rel(2), angle=0),
        axis.text.y=element_text(size=rel(2), angle=0),
        axis.title.x=element_text(size=rel(2), angle=0),
        axis.title.y=element_text(size=rel(2), angle=90),
        rect=element_rect(size=0),
        legend.position="none",
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent")
      )
    
  })
  
  output$distPlotNoInt <- renderPlot({
    outNoInt <-
      dataInputNoInt() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        # key2 = recode(
        #   key,
        #   S = "Susceptible (S)",
        #   #E = "Exposed (E)",
        #   I = "Infectious (I)",
        #   R = "Recovered (R)",
        #   V = "Vaccinated (V)"
        # ),
        # keyleft = recode(
        #   key,
        #   S = "Susceptible (S)",
        #   E = "Exposed (E)",
        #   I = "",
        #   R = "",
        #   V = "Vaccinated (V)"
        # ),
        keyright = recode(
          key,
          S = "",
          E = "",
          I = "Infectious (I), without vaccination",
          #R = "",
          #V = ""
          R = ""
          #V = ""
        )
      )
    ciao <- dataInput() %>% 
      gather(key, value, -time) %>%
      filter(key == "I") %>%
      mutate(
        id = row_number(),
        # key2 = recode(
        #   key,
        #   S = "Susceptible (S)",
        #   #E = "Exposed (E)",
        #   I = "Infectious (I)",
        #   R = "Recovered (R)",
        #   V = "Vaccinated (V)"
        # ),
        keyleft = recode(
          key,
          # S = "",
          # E = "",
          I = "Infectious (I), with vaccination",
          #R = "",
          #V = ""
          #R = ""
          #V = ""
        ),
        keyright = recode(
          key,
          #S = "",
          #E = "",
          I = "Infectious (I), with vaccination",
          #R = "",
          #V = ""
          #R = ""
          #V = ""
        )
      )
    out3 <-
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        # key2 = recode(
        #   key,
        #   S = "Susceptible (S)",
        #   #E = "Exposed (E)",
        #   I = "Infectious (I)",
        #   R = "Recovered (R)",
        #   V = "Vaccinated (V)"
        # ),
        keyleft = recode(
          key,
          S = "",
          E = "",
          I = "Infectious (I), with vaccination",
          #R = "",
          #V = ""
          R = ""
          #V = ""
        ),
        keyright = recode(
          key,
          S = "",
          E = "",
          I = "Infectious (I), with vaccination",
          #R = "",
          #V = ""
          R = ""
          #V = ""
        )
      )
    
    ggplot(data = outNoInt,
           aes(
             x = time,
             y = value,
             group = keyright,
             col = keyright,
             label = keyright,
             data_id = id
           )) + # ylim(0, 1) +
      ylab("Proportion of full population") + xlab("Time (days)") +
      geom_line(size = 2) + geom_line(data = ciao, aes(x=time, y=value), alpha = 0.5, size =2) +
      geom_text_repel(
        data = subset(outNoInt, time == max(time)),
        aes(label = keyright),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      #geom_text(x=max(out3$time)/2, y=0.97, label="Infectious (I), with vaccination", fontface =7, colour = "green4", size =6 , alpha =0.5) +
      #geom_text(x=max(out3$time)/2, y=0.89, label="Infectious (I), without vaccination", fontface =7, colour = "purple", size =6, alpha =0.5) +
      geom_text_repel(
        data = subset(ciao, time == min(time)),
        aes(label = keyleft),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      # geom_text_repel(
      #   data = subset(out, time == min(time)),
      #   aes(label = keyleft),
      #   size = 6,
      #   segment.size  = 0.2,
      #   segment.color = "grey50",
      #   nudge_x = 0,
      #   hjust = 0,
      #   direction = "y"
      # ) +
      theme(legend.position = "none") +
      #scale_colour_manual(values = c("orange", "green4", "black", "blue", "orange")) +
      scale_colour_manual(values = c("white","purple", "green4")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      theme(
        axis.text.x=element_text(size=rel(2), angle=0),
        axis.text.y=element_text(size=rel(2), angle=0),
        axis.title.x=element_text(size=rel(2), angle=0),
        axis.title.y=element_text(size=rel(2), angle=90),
        rect=element_rect(size=0),
        legend.position="none",
        panel.background=element_rect(fill="transparent", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent")
      )
    
  })
  
  
  output$plot4 <- renderImage({
    #filename <- normalizePath(file.path('./images',"model_diagram_1.png"))
    #filename <- file.path('./images',"model_diagram_1.png") non funziona
    filename <- "images/model_diagram_1.png"
    
    list(src = filename, height=80, width=500)
    
  }, deleteFile = FALSE)
  
  output$progressBox <- renderValueBox({
    valueBox(
      dataInput() %>% filter(time == max(time)) %>% select(R) %>% mutate(R = round(100 * R, 2)) %>% paste0("%"),
      subtitle = tags$p("Proportion of full population that got the disease by end of time frame in scenario with vaccination", style = "font-size: 140%;"),
      #"Proportion of full population that got the disease by end of time frame in scenario with vaccination",
      #icon = icon("thumbs-up", lib = "glyphicon"),
      color = "black"
    )
  })
  
  output$progressBoxNoInt <- renderValueBox({
    valueBox(
      dataInputNoInt() %>% filter(time == max(time)) %>% select(R) %>% mutate(R = round(100 * R, 2)) %>% paste0("%"),
      subtitle = tags$p("Proportion of full population that got the disease by end of time frame in scenario without vaccination", style = "font-size: 140%;"),
      #"Proportion of full population that got the disease by end of time frame in scenario without vaccination",
      #icon = icon("thumbs-up", lib = "glyphicon"),
      color = "black"
    )
  })
  
  # output$approvalBox <- renderValueBox({
  #   valueBox(
  #     paste0(round(
  #       100 * (dataInput() %>% filter(row_number() == n()) %>% mutate(res = (R + I) / (S + I + R)) %>% pull("res")), 2), "%"),
  #     "Proportion of susceptibles that will get the disease by end of time frame",
  #     #icon = icon("thermometer-full"),
  #     color = "black"
  #   )
  # })
  
  output$BRRBox <- renderValueBox({
    valueBox(
      paste0(round(input$connum *
                     (1 - input$pvac / 100 * input$vaceff / 100), 2), ""),
      subtitle = tags$p("Effective R0 (for population at outbreak, when immunity is taken into account)", style = "font-size: 140%;"),
      #"Effective R0 (for population at outbreak, when immunity is taken into account)",
      #icon = icon("arrows-alt"),
      color = "orange"
    )
  })
  
  output$HIBox <- renderValueBox({
    valueBox(
      paste0(round(100 * (1 - 1 / (input$connum)), 2), "%"),
      subtitle = tags$p("Proportion of population that needs to be immune for herd immunity", style = "font-size: 140%;"),
      #"Proportion of population that needs to be immune for herd immunity",
      #icon = icon("medkit"),
      color = "blue"
    )
  })  
}

# Run the application
shinyApp(ui = ui, server = server)
