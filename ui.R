library(shiny)
shinyUI(fluidPage(
  titlePanel(HTML("<h3>Initiation Strategies for Renal-Replacement Therapy in the Intensive Care Unit:<h3/><h5>a Precision Medicine Approach<h5/>")),
  sidebarLayout(
    sidebarPanel(
      p(HTML('Plug values from the time point when severe acute kidney injury occurs<br/>(KDIGO III or RIFLE failure stage)')),
      
      numericInput("sofa", "SOFA", 
                   value = 10, min = 0, max = 24, step = 1),
      numericInput("ph", "pH", 
                   value = 7.3, min = 6.85, max = 7.6, step = .01),
      numericInput("potassium", "Potassium (mmol/L)", 
                   value = 4.5, min = 2, max = 7.1, step = .1),
      numericInput("urea", "Urea (mmol/L)", 
                   value = 20, min = 1, max = 60, step = 1),
      numericInput("weight", "Weight (kg)", 
                   value = 80, min = 30, max = 210, step = 1),
      radioButtons("drug", HTML("Immunosuppressive Drug<br/>(non-corticosteroid)"),
                   c("Yes" = 1, "No" = 0), 
                   selected = 0),
      submitButton("Predict")
      ),
    mainPanel(
      #h4("Predicted Probablity of RRT Initiation Within 48 hours:"),
      br(), 
      h5(textOutput("prediction")),
      #p('Effects of early vs. delayed strategy on death at day 60 at different risks for RRT initiation'),
      imageOutput("myImage"), #, height = "100%", width = "100%")
      textOutput("dr_late")
    )
  )
))
