####################################
# Coagulation-flocculation model   #
# Jose A. D. Lopez-Coronado        #
# Agriculture Victoria - 2026      #
####################################

# Load R packages
library(shiny)
library(shinythemes)
library(markdown)

# Read models

get_rds <- function(file) {
  
  base <- "https://raw.githubusercontent.com/jlopezco-wm/tanflocpfsmodel/main/"
  dest <- tempfile(fileext = ".rds")
  
  download.file(
    url      = paste0(base, file),
    destfile = dest,
    mode     = "wb",
    quiet    = TRUE
  )
  
  readRDS(dest)
}

dose_mod    <- get_rds("dose_mod.rds") #Function to estimate the dose
tanfloc_rem <- get_rds("tanfloc_rem.rds") #Model with tanfloc data
pfs_rem <- get_rds("pfs_rem.rds") #Model with pfs data

environment(tanfloc_rem)    <- .GlobalEnv
environment(pfs_rem)    <- .GlobalEnv
environment(dose_mod)    <- .GlobalEnv

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "Tanfloc-PFS",
                  tabPanel("Calculator",
                           sidebarPanel(
                             selectInput(
                               inputId = "slc",
                               label = "Coagulant:",
                               choices = c("Tanfloc","Polyferric sulphate")
                             ),
                             numericInput("txt1", "Effluent volume (L):", value = NA, min = 0),
                             numericInput("txt2", "Effluent turbidity (NTU):", value = NA, min = 1000, max = 30000),
                             numericInput("txt3", "Coagulant concentration (g/L)", value = NA, min = 1),
                             helpText("Leave this field blank if unknown."),
                             actionButton("actionButton", "Calculate", class = "btn-primary")
                           ), # sidebarPanel
                           mainPanel(
                             verbatimTextOutput("contents"),
                             verbatimTextOutput("contents_model"),
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel(
                    "About",
                    uiOutput("about_md")
                  )
                ) # navbarPage
) # fluidPage

# server
# Define server function  
server <- function(input, output, session) {
  
  # Input Data
  model_output <- eventReactive(input$actionButton,{  
    
    # Validation
    validate(
      need(input$txt1 >= 0,
           "Effluent volume cannot be negative.")
    )
    
    validate(
      need(
        !is.na(input$txt2) && (input$txt2 >= 1000 && input$txt2 <= 30000),
        "Turbidity should be between 1000 and 30000 NTU."
      )
    )
    
    validate(
      need(is.na(input$txt3) || input$txt3 > 0,
           "Concentration should be higher than zero. Leave blank if unknown.")
    )
    
    withProgress(message = "Estimating dose…", value = 0, {
      
      incProgress(0.1)

    df <- data.frame(
      Name = c("effluent",
               "turbidity",
               "concentration"),
      
      Value = c(
        as.numeric(input$txt1),
        as.numeric(input$txt2),
        if (is.na(input$txt3)) {
          if (input$slc == "Tanfloc") {
            207
          } else {
            475
          }
        } else {
          as.numeric(input$txt3)
        }
      ),
      stringsAsFactors = FALSE)
    incProgress(0.3)
    
    
    model_selected <- if (input$slc == "Tanfloc") {
      tanfloc_rem
    } else {
      pfs_rem
    }
    
    #Call for model
    Output2 <- dose_mod(model = model_selected,
                           effluent = df[1,2],
                           turbidity = df[2,2],
                           concentration = df[3,2])
    incProgress(1)
    Output2
    })
  })
    
    # Status/Output Text Box
  
  output$contents <- renderPrint({
    
    # BEFORE first click
    if (input$actionButton == 0) {
      return("Server is ready for calculation.")
    }
    
    # AFTER click
    result <- try(model_output(), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      "Error!"
    } else {
      "Calculation complete."
    }
  })

    # Model text output
    output$contents_model <- renderText({
      req(input$actionButton)
      model_output()
    })
    
    output$about_md <- renderUI({
      
      about_url <- "https://raw.githubusercontent.com/jlopezco-wm/tanflocpfsmodel/main/About.md"
      
      tmp <- file.path(tempdir(), "About.md")
      
      if (!file.exists(tmp)) {
        download.file(
          about_url,
          tmp,
          mode = "wb",
          quiet = TRUE
        )
      }
      
      div(
        includeMarkdown(tmp),
        align = "justify"
      )
    })
} 

# Create Shiny object
shinyApp(ui = ui, server = server)
