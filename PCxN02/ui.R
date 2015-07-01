
shinyUI(
  pageWithSidebar(
    headerPanel(title="Pathway Coexpression Network (Canonical Pathways)",windowTitle="PCxN (CP)"),
    
    sidebarPanel(  
      # Select pathway
      includeHTML("www/js/tools.js"),
      selectInput("tool1", label = "Select pathway:", choices = path.names, selected = NULL, multiple = FALSE),
      
      br(),
      
      # Select the number top n connected pathways
      sliderInput("top.n", "Pathways to display:", 
                  min=0, max=50, value=15),

      br(),
      
      # select the BIC cut-off
      sliderInput("p.cut","p-value cut-off:",
                  min=0, max=0.1, value=0.05,step=0.001),

      br(),
      
      # select the correlation cut-off
      sliderInput("cor.cut","PathCor cut-off:",
                  min=0, max=1, value=0.05,step=0.025),
      
      
      br(),
      
      selectInput("method", "Top edges by:", 
                  choices = c("Absolute Value", "Decreasing", "Increasing")),

      br(),
      
      submitButton("Submit")#,
      

    ),
    
    mainPanel(
      h3(textOutput("caption")),
      tabsetPanel(
        tabPanel("Table", tableOutput("view"))
      )
      
    )
  )
)