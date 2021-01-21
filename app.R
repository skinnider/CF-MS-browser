options(stringsAsFactors = FALSE)
library(shiny)
library(tidyverse)
library(magrittr)
library(shinythemes)
library(DT)
library(RColorBrewer)

# read identifier choices
identifiers = readRDS('data/identifier_choices.rds')
# read experiment metadata
meta = readRDS("data/experiment_metadata.rds") %>%
  dplyr::rename(`Sample description` = Sample.description)
# read network species
network_species = readRDS("data/network_species.rds")

# load plotting functions
source("plots.R")
source("plot_network.R")

################################################################################
### Server
################################################################################
server = function(input, output, session) {
  output$choose_species = renderUI({
    selectInput(
      inputId = 'species',
      label = 'Species',
      choices = unique(meta$Species) %>% sort(),
      selected = 'Homo sapiens'
    )
  })
  
  output$choose_identifier = renderUI({
    # catch non-set species
    if (is.null(input$species))
      return(character(0))
    # read a list of choices from a file
    choices_file = paste0('data/identifiers/choices_per_species/',
                          chartr(' ', '_', input$species), '.rds')
    choices = readRDS(choices_file) %>%
      ## manually filter out EMBL and PDB
      setdiff(c('EMBL', 'PDB'))
    selectInput(
      inputId = 'identifier',
      label = 'Identifiers',
      choices = choices,
      selected = 'Gene name'
    )
  })
  
  output$choose_values = renderUI({
    suppressWarnings({
      # catch non-set species or identifier 
      if (is.null(input$species) | is.null(input$identifier))
        return(character(0))
      
      # read possible identifier values
      identifier_clean = identifiers %>%
        filter(identifier == input$identifier) %>%
        pull(clean)
      species_clean = chartr(' ', '_', input$species)
      value_file = paste0("data/identifiers/values_per_species/", species_clean,
                          '/', identifier_clean, ".rds")
      idmap = readRDS(value_file)
      
      # pull out unique IDs
      choices = unique(idmap$id) %>% sort()
      
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
        options = list(maxOptions = 12, maxItems = 5))
    })
  })
  
  output$choose_experiments = renderUI({
    # catch non-set species or values 
    if (is.null(input$species) || is.null(input$values) ||
        length(input$values) == 0)
      return(character(0))
    
    # get experiments
    expts0 = filter(meta, Species == input$species)
    
    # read a list of potential identifiers
    species_clean = chartr(' ', '_', input$species)
    identifier_clean = identifiers %>%
      filter(identifier == input$identifier) %>%
      pull(clean)
    ids_per_expt_file = paste0("data/identifiers/values_per_experiment/",
                               species_clean, "/", identifier_clean, ".rds")
    ids_per_expt = readRDS(ids_per_expt_file)
    choices = ids_per_expt %>%
      filter(id %in% input$values) %>%
      # filter on _all_, not _any_
      group_by(Accession, Replicate) %>%
      filter(n_distinct(id) == n_distinct(input$values)) %>%
      ungroup() %>%
      distinct(Accession, Replicate) %>% 
      left_join(meta, by = c('Accession', 'Replicate')) %$%
      paste0(Accession, ': ', Replicate, " (", Fractionation, ")") %>%
      sort()
    if (length(choices) == 1 && choices == ":  ()")
      choices = character(0)
    
    selectizeInput(
      inputId = 'experiments', 
      label = 'Experiment(s)', 
      choices = choices,
      multiple = TRUE,
      selected = character(0),
      options = list(maxItems = 4))
  })
  
  output$chromatogram_plot = renderPlot({
    withProgress(
      plot_chromatograms(input$experiments,
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
      experiments = gsub(" \\(.*$", "", input$experiments)
      accessions = gsub(": .*$", "", experiments)
      replicates = gsub("^.*: ", "", experiments)
      idmap_files = paste0("data/idmaps/", accessions, "-", replicates, ".rds")
      idmap = map(idmap_files, readRDS) %>% bind_rows()
      # filter to the identifier in question
      idmap %<>% filter(identifier == !!input$identifier)
      
      # read the chromatograms
      chrom_files = paste0("data/chromatograms/", accessions, "-", replicates,
                           ".rds")
      chroms = map(chrom_files, ~ readRDS(.) %>%
                     left_join(idmap, by = 'uniprot') %>%
                     drop_na(id) %>%
                     filter(id %in% input$values)) %>%
        bind_rows()
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
  
  output$metadata = renderDT(
    meta #, options = list(lengthChange = FALSE)
  )
  
  output$choose_network_species = renderUI({
    selectInput(
      inputId = 'network_species',
      label = 'Species',
      choices = network_species,
      selected = 'Homo sapiens'
    )
  })
  
  output$choose_network_identifier = renderUI({
    # catch non-set species
    if (is.null(input$network_species))
      return(character(0))
    # read a list of choices from a file
    choices_file = paste0('data/networks/identifier_choices/',
                          chartr(' ', '_', input$network_species), '.rds')
    choices = readRDS(choices_file)
    selectInput(
      inputId = 'network_identifier',
      label = 'Identifiers',
      choices = choices,
      selected = ifelse(input$network_species %in% c('BOP clade',
                                                     'Deuterostomia',
                                                     'Euarchontoglires',
                                                     'Eukaryota',
                                                     'Mesangiospermae',
                                                     'Opiskothonta',
                                                     'Plasmodium',
                                                     'Tetrapoda'),
                        'eggNOG (euk)', 'Gene name')
    )
  })
  
  output$choose_network_values = renderUI({
    suppressWarnings({
      # catch non-set species or identifier 
      if (is.null(input$network_species) | is.null(input$network_identifier))
        return(character(0))
      
      # read map
      map_file = paste0('data/networks/identifier_maps/',
                        chartr(' ', '_', input$network_species), '.rds')
      maps = readRDS(map_file)
      
      # extract the relevant one
      map = maps[[input$network_identifier]]
      
      # pull out unique IDs
      choices = unique(map$id) %>% sort()
      
      selectizeInput(
        inputId = 'network_values', 
        label = 'Protein(s)', 
        choices = choices,
        multiple = TRUE,
        # selected = intersect(c('STAT1'), choices),
        options = list(maxOptions = 12, maxItems = 5))
    })
  })
  
  output$network_plot <- renderVisNetwork({
    plot_network(species = input$network_species,
                 identifier = input$network_identifier,
                 values = input$network_values, 
                 plot = input$network_plot_opt)
  })
  
  output$download_network = downloadHandler(
    filename = function() {
      names = input$network_values[order(input$network_values)] %>% unique()
      paste0(paste(names, collapse = "-"), ".csv")
    },
    content = function(file) {
      # first, read the species network
      species_clean = chartr(' ', '_', input$network_species)
      network_file = paste0('data/networks/', species_clean, '.rds')
      network = readRDS(network_file)
      # strip eggNOG scopes
      network %<>% mutate_at(vars(protein_A, protein_B), ~ gsub("@.*$", "", .))
      
      # convert to graph
      g = network %>%
        dplyr::select(protein_A, protein_B) %>%
        graph_from_data_frame(directed = FALSE)
      
      # second, read the ID map
      map_file = paste0('data/networks/identifier_maps/', species_clean, '.rds')
      maps = readRDS(map_file)
      map = maps[[input$network_identifier]] %>% dplyr::select(-identifier)
      
      # find the proteins of interest from the ID map
      proteins = filter(map, id %in% input$network_values) %$%
        unique(original_id)
      
      # optionally, filter the entire network down to proteins/neighbors
      if (input$network_plot_opt == "Only query proteins") {
        g %<>% induced_subgraph(proteins)
      } else if (input$network_plot_opt == 
                 "Query proteins and their neighbors") {
        neighbors = map(proteins, ~ names(neighbors(g, .))) %>%
          unlist() %>%
          unique()
        g %<>% induced_subgraph(union(proteins, neighbors))
      }
      
      # set up edge list
      map_reduced = map %>%
        group_by(original_id) %>%
        summarise(all_ids = paste0(sort(id), collapse = ';')) %>%
        ungroup() %>%
        dplyr::rename(id = original_id) %$%
        setNames(all_ids, id)
      edges = g %>%
        as_data_frame()
      if (!(input$network_species == 'Homo sapiens' &
            input$network_identifier == 'Gene name') &
          !(input$network_species != 'Homo sapiens' &
            input$network_identifier == 'eggNOG (euk)')) {
        edges %<>%
          mutate(alt_A = map_reduced[from], alt_B = map_reduced[to]) %>%
          set_colnames(c("Protein A", "Protein B", 
                         "Alternate identifier(s) A", 
                         "Alternate identifier(s) B"))
      } else {
        edges %<>% set_colnames(c("Protein A", "Protein B"))
      }
      
      # return 
      write.csv(edges, file, row.names = FALSE)
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
                               p(),
                               "First, select a species.", 
                               p(),
                               "Next, specify an identifier (e.g., gene name,",
                               " UniProt ID, or eggNOG orthogroup) and",
                               " choose up to five proteins to plot.", 
                               p(), 
                               "This will generate a list of CF-MS datasets", 
                               " containing all of the selected proteins.", 
                               " Choose up to four datasets at once in which",
                               " to plot chromatograms for the selected",
                               " proteins.",
                               p(),
                               "To see a full list of datasets and their", 
                               " associated metadata, such as species or", 
                               " fractionation approach, see the \"Metadata\"",
                               " tab."),
                      
                      # Select species
                      uiOutput('choose_species'),
                      
                      # Select identifier
                      uiOutput('choose_identifier'),
                      
                      # Select gene(s)
                      uiOutput('choose_values'),
                      
                      # Select experiment
                      uiOutput('choose_experiments'),
                      
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
                  "Metadata",
                  fluidPage(DTOutput('metadata'))
                ),
                tabPanel(
                  "Networks",
                  sidebarLayout(
                    sidebarPanel(
                      helpText(
                        "Explore consensus CF-MS interactomes for 27 species",
                        " or clades.", 
                        p(),
                        "Specify a species or clade, and up to ten proteins",
                        " for which to plot protein interaction networks."),
                      
                      # Select species
                      uiOutput('choose_network_species'),
                      
                      # Select identifier
                      uiOutput('choose_network_identifier'),
                      
                      # Select gene(s)
                      uiOutput('choose_network_values'),
                      
                      selectInput("network_plot_opt", 
                                  label = "What to plot?", 
                                  width = "100%",
                                  choices = c(
                                    "Only query proteins", 
                                    "Query proteins and their neighbors"),
                                  selected = "Query proteins and their neighbors"),
                      
                      downloadButton("download_network", 
                                     "Download data")
                    ),
                    mainPanel(
                      visNetworkOutput("network_plot",
                                       width = "700px",
                                       height = "700px")
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
