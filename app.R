library(shiny)
library(VariantAnnotation)
library(shinyWidgets)

# Define UI
ui <- fluidPage(
  titlePanel("VCF File Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("vcf_file", "Choose a VCF file"),
      uiOutput("allele_freq_selector"),
      uiOutput("variant_type_selector"),
      uiOutput("pathogenicity_key_selector")
    ),
    mainPanel(
      h4("Analysis Results:"),
      fluidRow(
        column(6, tableOutput("variant_counts")),
        column(6, plotOutput("variant_types_plot", width = "100%")),
        column(12, tableOutput("pathogenic_variants"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Initialize reactive values to store INFO keys
  allele_freq_keys <- reactiveVal(character(0))
  variant_type_keys <- reactiveVal(character(0))
  pathogenic_keys <- reactiveVal(character(0))
  
  # Read VCF file
  vcf_data <- reactive({
    req(input$vcf_file)
    vcf <- readVcf(input$vcf_file$datapath)
    
    # Extract INFO field
    info_fields <- info(vcf)
    info_keys <- names(info_fields)
    
    # Update reactive values with INFO keys
    allele_freq_keys(info_keys)
    variant_type_keys(info_keys)
    pathogenic_keys(info_keys)
    
    return(vcf)
  })
  
  # Update UI elements for selecting keys
  observe({
    if (!is.null(vcf_data())) {
      output$allele_freq_selector <- renderUI({
        shinyWidgets::pickerInput("allele_freq_key", "Select Allele Frequency Key:",
                                  choices = allele_freq_keys(), selected = NULL,
                                  multiple = TRUE)
      })
      
      output$variant_type_selector <- renderUI({
        shiny::selectInput("variant_type_key", "Select Variant Type Key:",
                           choices = c("", variant_type_keys()), selected = "")
      })
      
      output$pathogenicity_key_selector <- renderUI({
        shiny::selectInput("pathogenic_key", "Select Pathogenicity Key:",
                           choices = c("", pathogenic_keys()), selected = "")
      })
    }
  })
  
  # Perform analysis when VCF file is loaded and keys are selected
  observeEvent({
    input$allele_freq_key
    input$variant_type_key
    input$pathogenic_key
  }, {
    # Perform analysis if all keys are selected
    if (!is.null(input$allele_freq_key) && input$variant_type_key != "" && input$pathogenic_key != ""){
      analyze_vcf(vcf_data(), input$allele_freq_key, input$variant_type_key, input$pathogenic_key)
    } else {
      # Clear plot if no option is selected
      output$variant_types_plot <- renderPlot(NULL)
    }
  })
  
  # Function to perform analysis
  analyze_vcf <- function(vcf, allele_freq_keys, variant_type_key, pathogenic_key) {
    # Variant counts
    num_variants <- nrow(vcf)
    
    # Perform analysis for each allele frequency key
    num_variants_high_freq <- sapply(allele_freq_keys, function(allele_freq_key) {
      tryCatch({
        sum(info(vcf)[[allele_freq_key]] > 0.5, na.rm = TRUE)
      }, error = function(e) {
        return(0)
      })
    })
    
    # Filter pathogenic variants based on pathogenic key
    pathogenic_variants <- tryCatch({
      clnsig_values <- info(vcf)[[pathogenic_key]]
      print(clnsig_values)
      subset(vcf, clnsig_values %in% c("Pathogenic", "Likely_pathogenic"))
    }, error = function(e) {
      return(data.frame())
    })
    
    # Extract all pathogenic key values into a dataframe
    df_pathogenic <- as.data.frame(info(vcf)[[pathogenic_key]])
 
    # Count Pathogenic and Likely_pathogenic variants separately
    num_pathogenic <- nrow(subset(df_pathogenic, df_pathogenic$value == "Pathogenic"))
    num_likely_pathogenic <- nrow(subset(df_pathogenic, df_pathogenic$value == "Likely_pathogenic"))
    
    # Create separate table for pathogenic variants
    pathogenic_table <- data.frame(
      "Pathogenic" = num_pathogenic,
      "Likely_pathogenic" = num_likely_pathogenic
    )
    
    # Perform analysis for variant type key
    variant_types <- tryCatch({
      if (variant_type_key %in% names(info(vcf))) {
        as.data.frame(table(info(vcf)[[variant_type_key]]))
      } else {
        data.frame()
      }
    }, error = function(e) {
      return(data.frame())
    })
    
    # Update outputs
    output$variant_counts <- renderTable({
      data.frame(
        "Number of Variants" = num_variants,
        "Variants with Frequency > 0.5" = sum(num_variants_high_freq)
      )
    }, rownames = FALSE)
    
    output$pathogenic_variants <- renderTable({
      pathogenic_table
    }, rownames = FALSE)
    
    output$variant_types_plot <- renderPlot({
      if (nrow(variant_types) > 0) {
        barplot(variant_types$Freq, names.arg = variant_types$Var1,
                xlab = "Variant Type", ylab = "Frequency", main = "Variant Types")
      } else {
        plot(1, type = "n", xlab = "", ylab = "", main = "Variant Types")
        text(1, 1, "No variant types found", cex = 1.2, col = "red", font = 2)
      }
    })
  }
}

# Run the application
shinyApp(ui = ui, server = server)