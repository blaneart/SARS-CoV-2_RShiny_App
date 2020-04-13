library(shiny)

ui <- fluidPage(tabsetPanel(
  tabPanel(title = "GC content",
           plotOutput("GC")), 
  tabPanel(title = "Amino Acid Table",
           DT::dataTableOutput("AA")), 
  tabPanel(title = "Protein Weights",
           plotOutput("PW"))))
