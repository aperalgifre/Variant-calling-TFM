library(shiny)
library(VariantAnnotation)
library(data.table)
library(shinyWidgets)
library(RTCGA)
library(DT)
library(ggplot2)

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        body {
          background-color: #f9f9f9;
        }
        .title {
          font-family: 'Arial', sans-serif;
          color: #333333;
          font-size: 24px;
          font-weight: bold;
          margin-bottom: 20px;
        }
        .subtitle {
          font-family: 'Arial', sans-serif;
          color: #666666;
          font-size: 18px;
          margin-bottom: 10px;
        }
        .panel-heading {
          background-color: #f0f0f0;
          border: 1px solid #cccccc;
          padding: 10px;
        }
        .sidebar-panel {
          background-color: #ffffff;
          border: 1px solid #cccccc;
          padding: 20px;
          margin-bottom: 20px;
        }
        .main-panel {
          background-color: #ffffff;
          border: 1px solid #cccccc;
          padding: 20px;
          margin-bottom: 20px;
        }
        .btn-primary {
          background-color: #007bff;
          border-color: #007bff;
        }
        .btn-primary:hover {
          background-color: #0056b3;
          border-color: #0056b3;
        }
      ")
    )
  ),
  titlePanel("VCF File Analysis"),
  sidebarLayout(
    sidebarPanel(
      class = "sidebar-panel",
      fileInput("vcf_file", "Choose a VCF file"),
      h3("Select Analysis Parameters", class = "subtitle"),
      uiOutput("allele_freq_selector"),
      uiOutput("variant_type_selector"),
      uiOutput("clinical_significance_selector"),
      uiOutput("clinical_significance_categories"),
      uiOutput("cancer_type_selector"),  
      actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
      downloadButton("download_csv", "Download Filtered Variants CSV", class = "btn-primary")
    ),
    mainPanel(
      class = "main-panel",
      h4("Analysis Results:", class = "title"),
      fluidRow(
        column(6, 
               conditionalPanel(condition = "input.run_analysis",
                                DTOutput("variant_counts"))),
        column(6, 
               conditionalPanel(condition = "input.run_analysis",
                                plotOutput("variant_types_plot", width = "100%"))),
        column(12, 
               conditionalPanel(condition = "input.run_analysis",
                                DTOutput("clinical_significance_variants")))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Define a reactive value for tcga_data
  tcga_data <- reactiveVal(NULL)
  
  # Initialize reactive values to store INFO keys
  allele_freq_keys <- reactiveVal(character(0))
  variant_type_keys <- reactiveVal(character(0))
  clinical_significance_keys <- reactiveVal(character(0))
  clinical_significance_categories <- reactiveVal(character(0))
  
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
    clinical_significance_keys(info_keys)
    
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
      
      output$clinical_significance_selector <- renderUI({
        shiny::selectInput("clinical_significance_key", "Select Clinical Significance Key:",
                           choices = c("", clinical_significance_keys()), selected = "")
      })
      
      output$cancer_type_selector <- renderUI({
        selectInput("cancer_type", "Select TCGA Cancer Type:",
                    choices = c("No Cancer Selected", "ACC","BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBMLGG", "GBM", "HNSC", 
                                "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC",
                                "SKCM", "STAD", "STES", "TGCT", "THCA", "UCEC", "UCS", "UVM"), selected = "No Cancer Selected")
      })
    }
  })
  
  # Update UI elements for selecting clinical significance categories
  observe({
    if (!is.null(vcf_data()) && !is.null(input$clinical_significance_key)) {
      # Extract clinical significance categories of interest
      clinical_significance_categories(unique(info(vcf_data())[[input$clinical_significance_key]]))
      
      output$clinical_significance_categories <- renderUI({
        selectizeInput("clinical_significance_categories", "Select Clinical Significance Categories:",
                       choices = clinical_significance_categories(), multiple = TRUE)
      })
    }
  })
  
  # Perform analysis when "Run analysis" button is clicked
  observeEvent(input$run_analysis, {
    # Perform analysis if all keys are selected
    if (!is.null(input$allele_freq_key) && input$variant_type_key != "" && input$clinical_significance_key != "" && !is.null(input$clinical_significance_categories)){
      analyze_vcf(vcf_data(), input$allele_freq_key, input$variant_type_key, input$clinical_significance_key, input$clinical_significance_categories)
    }
  })
  
  # Function to perform analysis
  analyze_vcf <- function(vcf, allele_freq_keys, variant_type_key, clinical_significance_key, clinical_significance_categories) {
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
    
    # Filter clinical significance variants based on selected categories
    clinical_significance_variants <- tryCatch({
      interest_values <- info(vcf)[[clinical_significance_key]]
      subset(vcf, interest_values %in% input$clinical_significance_categories)
    }, error = function(e) {
      return(data.frame())
    })
    
    # Create a data frame for clinical signficance variants
    pathogenic_table <- data.frame(Category = clinical_significance_categories, Count = numeric(length(clinical_significance_categories)), 
                                   stringsAsFactors = FALSE)
    
    # Count occurrences for each selected category
    category_counts <- table(factor(info(vcf)[[clinical_significance_key]], levels = clinical_significance_categories))
    pathogenic_table$Count <- as.numeric(category_counts)
    
    # Add a row for the total count of each category
    total_count <- colSums(pathogenic_table[, "Count", drop = FALSE])
    total_row <- c("Total", total_count)
    pathogenic_table <- rbind(pathogenic_table, total_row)
    
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
    output$variant_counts <- renderDT({
      data.frame(
        "Number of Variants" = num_variants,
        "Variants with Frequency > 0.5" = sum(num_variants_high_freq)
      )
    }, rownames = FALSE,
    options = list(
      dom ='t',
      searching = FALSE
    ))
    
    output$clinical_significance_variants <- renderDT({
      pathogenic_table
    }, rownames = FALSE,
    options = list(
      dom ='t',
      searching = FALSE,
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().container()).css({'padding-top': '50px'});",
        "}")
    ))
    
    output$variant_types_plot <- renderPlot({
      if (nrow(variant_types) > 0) {
        ggplot(variant_types, aes(x = Var1, y = Freq)) +
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
          labs(x = "Variant Type", y = "Frequency", title = "Variant Types")
      } else {
        plot(1, type = "n", xlab = "", ylab = "", main = "Variant Types")
        text(1, 1, "No variant types found", cex = 1.2, col = "red", font = 2)
      }
    })
  }
  
  # Define a reactive expression to filter variants and return the filtered dataframe
  filtered_variants <- reactive({
    req(input$vcf_file, input$clinical_significance_key, input$clinical_significance_categories)
    
    # Filter clinical significance variants based on selected key and categories
    clinical_significance_variants <- tryCatch({
      interest_values <- info(vcf)[[input$clinical_significance_key]]
      subset(vcf, interest_values %in% input$clinical_significance_categories)
    }, error = function(e) {
      return(data.frame())
    })
    
    # Process the filtered variants
    chrom <- as.data.frame(seqnames(rowRanges(clinical_significance_variants)))
    ranges <- (as.data.frame(ranges(rowRanges(clinical_significance_variants)))[, -c(2,3)])
    ref <- as.data.frame(ref(clinical_significance_variants))
    alt <- as.data.frame(alt(clinical_significance_variants))[, -c(1,2)]
    qual <- as.data.frame(qual(clinical_significance_variants))
    filter <- as.data.frame(filt(clinical_significance_variants))
    info <- as.data.frame(info(clinical_significance_variants))
    
    partial_df <- cbind(chrom, ranges, ref, alt, qual, filter)
    colnames(partial_df) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
    completed_df <- cbind(partial_df, info)
    rownames(completed_df) <- NULL
    
    return(completed_df)
  })
  
  # Function to download TCGA data and generate CSV
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("filtered_variants", ".csv")
    },
    content = function(file) {
      req(input$cancer_type)
      
      # Check if "No Cancer Selected" option is chosen
      if (input$cancer_type == "No Cancer Selected") {
        # Return filtered dataframe without adding the extra column
        fwrite(filtered_variants(), file, row.names = FALSE)
        return()
      }
      
      # Download TCGA data
      temp_dir <- tempdir()
      downloadTCGA(input$cancer_type, dataSet = 'Mutation_Packager_Calls.Level', destDir = temp_dir)
      
      # Construct file path for downloaded data
      data_folder <- file.path(temp_dir, paste0("gdac.broadinstitute.org_", input$cancer_type, ".Mutation_Packager_Calls.Level_3.2016012800.0.0"))
      
      # Read the data 
      data <- readTCGA(data_folder, 'mutations')  
      
      # Get filtered variants dataframe
      filtered_df <- filtered_variants()
      
      # Compare with TCGA data and add extra column
      if (!is.null(data)) {
        filtered_df$Present_in_TCGA <- ifelse(apply(filtered_df, 1, function(row) {
          any(apply(data, 1, function(tcga_row) {
            paste(row["CHROM"], row["POS"], row["ALT"]) == paste(tcga_row["Chromosome"], tcga_row["Start_Position"], tcga_row["Tumor_Seq_Allele2"])
          }))
        }), "Yes", "No")
      } else {
        filtered_df$Present_in_TCGA <- NA
      }
      
      # Write filtered dataframe to CSV file
      fwrite(filtered_df, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui, server)