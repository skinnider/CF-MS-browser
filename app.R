options(stringsAsFactors = FALSE)
library(shiny)
library(tidyverse)
library(magrittr)
library(shinythemes)

# read all chromatograms
tidy = readRDS('/Users/michaelskinnider/git/CF-MS-browser/data/tidy_chromatograms.rds')
# read identifier choices
identifiers = readRDS('~/git/CF-MS-browser/data/identifier_choices.rds')
# read identifier map
idmap = readRDS('~/git/CF-MS-browser/data/identifier_map.rds')

# initially, populate with gene names
gene_names = idmap %>%
  filter(identifier == 'Gene name') %>%
  pull(id) %>%
  unique() %>% sort()
## A first
a_idx = dplyr::first(which(startsWith(gene_names, 'A')))
special = gene_names[1:(a_idx - 1)]
gene_names = c(gene_names[a_idx:length(gene_names)], special)

# load plotting functions
source("plots.R")

################################################################################
### Server
################################################################################
server = function(input, output, session) {
  # update values choices when identifier is selected
  # choices = reactive({
  #   choices = filter(idmap, identifier == input$identifier) %>%
  #         pull(id) %>%
  #         unique() %>%
  #         sort()
  #       ## special characters at the end
  #       a_idx = dplyr::first(which(startsWith(choices, 'A')))
  #       special = choices[1:(a_idx - 1)]
  #       choices = c(choices[a_idx:length(choices)], special)
  #   return(choices)
  # })
  output$choose_values = renderUI({
    choices = filter(idmap, identifier == input$identifier) %>%
      pull(id) %>%
      unique() %>%
      sort()
    ## special characters at the end
    a_idx = dplyr::first(which(startsWith(choices, 'A')))
    if (!is.na(a_idx)) {
      special = choices[1:(a_idx - 1)]
      choices = c(choices[a_idx:length(choices)], special)
    }
    
    selectizeInput(
      inputId = 'values', 
      label = 'Protein(s)', 
      choices = choices,
      multiple = TRUE,
      selected = intersect(c('PSMA1', 'PSMA2', 'PSMA3'),
                           choices),
      options = list(maxOptions = 10))
  })
  
  # observeEvent(input$identifier, {
  #   choices = filter(idmap, identifier == input$identifier) %>%
  #     pull(id) %>%
  #     unique() %>%
  #     sort()
  #   ## special characters at the end
  #   a_idx = dplyr::first(which(startsWith(choices, 'A')))
  #   special = choices[1:(a_idx - 1)]
  #   choices = c(choices[a_idx:length(choices)], special)
  #   
  #   updateSelectizeInput(session, "values", choices = choices)
  # })
  
  output$chromatogram_plot = renderPlot({
    withProgress(
      plot_chromatograms(input$identifier,
                         input$values, 
                         tidy,
                         map = idmap,
                         mode = input$plot_mode,
                         log_transform = input$log_transform),
      message = 'Updating plot...'
    )
  })
  
  output$download_data = downloadHandler(
    filename = function() {
      names = input$values[order(input$values)] %>% unique()
      paste0(paste(names, collapse = "-"), ".csv")
    },
    content = function(file) {
      map0 = filter(idmap, identifier == input$identifier, id %in% input$values)
      chroms = filter(tidy, uniprot %in% map0$uniprot) %>%
        left_join(map0, by = 'uniprot') %>%
        distinct(accession, replicate, species, protein_group, identifier, id,
                 fraction, intensity)
      write.csv(chroms, file, row.names = F)
    }
  )
}

################################################################################
### UI
################################################################################
ui = navbarPage("CF-MS explorer",
                theme = shinytheme("yeti"),
                tabPanel(
                  "Browse",
                  sidebarLayout(
                    sidebarPanel(
                      helpText("Visualize protein chromatograms from up to 206",
                               " CF-MS experiments.",
                               "Specify an identifier (e.g., gene name or",
                               " UniProt ID), then select the proteins in the",
                               " combined dataset you would like to plot."),
                      # Select identifier
                      selectInput(
                        inputId = 'identifier', 
                        label = 'Identifier(s)', 
                        choices = identifiers$identifier,
                        selected = 'Gene name'),
                      # Select gene(s)
                      uiOutput('choose_values'),
                      # selectizeInput(
                      #   inputId = 'values', 
                      #   label = 'Protein(s)', 
                      #   choices = gene_names,
                      #   multiple = TRUE,
                      #   selected = c('PSMA1', 'PSMA2', 'PSMA3'),
                      #   options = list(maxOptions = 10)),
                      # Plot mode
                      radioButtons("plot_mode", "Plot mode", 
                                   list("Scatterplot", "Heatmap"),
                                   selected = "Scatterplot"),
                      # Log-transform?
                      checkboxInput(
                        inputId = 'log_transform',
                        label = 'Log-transform?'
                      ),
                      # Download
                      downloadButton("download_data", "Download data")
                    ),
                    mainPanel(
                      plotOutput("chromatogram_plot", height = "1600px")
                    )
                  )
                ),
                tabPanel(
                  "About",
                  fluidRow(
                    shinydashboard::box(width = 6, 
                                        includeMarkdown("About.md"))
                  )
                )
)

shinyApp(ui, server)
