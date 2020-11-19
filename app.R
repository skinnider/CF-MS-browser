options(stringsAsFactors = FALSE)
library(shiny)
library(tidyverse)
library(magrittr)
library(shinythemes)

# read all chromatograms
# tidy = readRDS('~/git/CF-MS-browser/data/tidy_chromatograms.rds')
# read identifier choices
identifiers = readRDS('~/git/CF-MS-browser/data/identifier_choices.rds')
# read experiment choices
expts = readRDS('~/git/CF-MS-browser/data/experiment_choices.rds')
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
    suppressWarnings({
      # read the idmap
      accession = gsub(": .*$", "", input$experiment)
      replicate = gsub("^.*: ", "", input$experiment)
      idmap = readRDS(paste0("data/idmaps/", accession, "-", replicate, ".rds"))
      
      if (is.null(input$identifier))
        return(character(0))
      
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
  })
  
  output$choose_identifier = renderUI({
    # read the idmap
    accession = gsub(": .*$", "", input$experiment)
    replicate = gsub("^.*: ", "", input$experiment)
    idmap = readRDS(paste0("data/idmaps/", accession, "-", replicate, ".rds"))
    choices = unique(idmap$identifier)
    # gene name, then alphabetically, then eggNOGs
    choices = c('Gene name',
                sort(choices[!grepl("eggNOG", choices)]),
                sort(choices[grepl("eggNOG", choices)])) %>%
      intersect(choices)
    
    selectInput(
      inputId = 'identifier',
      label = 'Identifiers',
      choices = choices,
      selected = 'Gene name')
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
      plot_chromatograms(input$experiment,
                         input$identifier,
                         input$values, 
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
      # read the id map
      accession = gsub(": .*$", "", input$experiment)
      replicate = gsub("^.*: ", "", input$experiment)
      idmap = readRDS(paste0("data/idmaps/", accession, "-", replicate, ".rds"))
      # filter to the identifier in question
      idmap %<>% filter(identifier == !!input$identifier)
      
      # read the chromatograms
      chroms = readRDS(paste0("data/chromatograms/", accession, "-", replicate, 
                              ".rds")) %>%
        left_join(idmap, by = 'uniprot') %>%
        drop_na(id)
      # filter to the values in question
      chroms %<>% filter(id %in% input$values)
      # flag protein group AND query identifier value
      chroms %<>%
        mutate(group = paste0(protein_group, ' (', id, ')'))
      # return 
      chroms %<>% 
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
                               "Specify a CF-MS experiment and an identifier ",
                               "(e.g., gene name or UniProt ID), then select ",
                               "the proteins in the dataset you would",
                               " like to plot."),
                      
                      # Select experiment
                      selectInput(
                        inputId = 'experiment', 
                        label = 'Experiment', 
                        choices = expts,
                        selected = 'PXD001220: rep1'),
                      
                      # Select identifier
                      uiOutput('choose_identifier'),
                      # selectInput(
                      #   inputId = 'identifier', 
                      #   label = 'Identifiers', 
                      #   choices = identifiers$identifier,
                      #   selected = 'Gene name'),
                      
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
                      plotOutput("chromatogram_plot", height = "400px")
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
