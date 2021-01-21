library(tidyverse)
library(magrittr)
library(igraph)
library(visNetwork)

plot_network = function(species,
                        identifier,
                        values, 
                        plot = c("Only query proteins",
                                 "Query proteins and their neighbors",
                                 "Entire network"),
                        remove = NULL,
                        interactive = T) {
  if (is.null(species) ||
      is.null(identifier) ||
      is.null(values) ||
      length(values) == 0)
    return(ggplot())
  
  # first, read the species network
  species_clean = chartr(' ', '_', species)
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
  map = maps[[identifier]] %>% dplyr::select(-identifier)
  
  # find the proteins of interest from the ID map
  proteins = filter(map, id %in% values) %$%
    unique(original_id)
  
  # optionally, filter the entire network down to proteins/neighbors
  if (plot == "Only query proteins") {
    g %<>% induced_subgraph(proteins)
  } else if (plot == "Query proteins and their neighbors") {
    neighbors = map(proteins, ~ names(neighbors(g, .))) %>%
      unlist() %>%
      unique()
    g %<>% induced_subgraph(union(proteins, neighbors))
  }
  
  # set gene attributes
  g %<>% set_vertex_attr('Highlight', value = F)
  g %<>% set_vertex_attr('Highlight', index = proteins, value = T)
  
  # set up color palette
  node_cols = c('TRUE' = brewer.pal(8, 'Set2')[3], 'FALSE' = '#cccccc')
  edge_cols = c('Between neighbors' = '#cccccc', 
                'Query-neighbor' = brewer.pal(8, 'Set2')[3],
                'Between query genes' = '#000000')
  
  # set up objects to plot
  ## nodes
  map_reduced = map %>%
    mutate(id = ifelse(id %in% values, paste0('<strong>', id, '</strong>'),
                       id)) %>%
    group_by(original_id) %>%
    summarise(all_ids = paste0(sort(id), collapse = '; ')) %>%
    ungroup() %>%
    dplyr::rename(id = original_id)
  nodes = data.frame(id = names(V(g)),
                     highlight = vertex_attr(g, 'Highlight')) %>%
    left_join(map_reduced, by = 'id') %>%
    mutate(color = node_cols[as.character(highlight)],
           label = id)
  if ((species == 'Homo sapiens' & identifier == 'Gene name') | 
      (species != 'Homo sapiens' & identifier == 'eggNOG (euk)')) {
    nodes %<>% mutate(
      title = paste0("<span style='font-family: Helvetica, Arial, ",
                     "sans-serif; font-size: 10px;'><strong>", id,
                     "</strong></span>"))
  } else {
    nodes %<>% mutate(
      title = paste0("<span style='font-family: Helvetica, Arial, ",
                     "sans-serif; font-size: 10px;'><strong>", id,
                     "</strong><br />(", identifier, ": ", all_ids,
                     ")</span>"))
  }
  
  edges = g %>%
    as_data_frame() %>%
    mutate(class = ifelse(from %in% proteins,
                          ifelse(to %in% proteins, 'Between query genes',
                                 'Query-neighbor'),
                          ifelse(to %in% proteins, 'Query-neighbor',
                                 'Between neighbors')),
           color = edge_cols[class],
           width = ifelse(class == 'Between neighbors', 1.5,
                          ifelse(class == "Query-neighbor", 3, 5)),
           # color = edge_cols[class],
           # width = ifelse(class == 'Stimulated', 5, 2),
           title = paste0("<span style='font-family: Helvetica, Arial, ",
                          "sans-serif; font-size: 10px;'>",
                          from, "&mdash;", to, 
                          "</span>")
    )
  
  # plot
  visNetwork(nodes, edges, width = "100%") %>%
    # disable stabilization
    visPhysics(stabilization = FALSE) %>%
    # disable smooth edges
    visEdges(smooth = FALSE) %>%
    # use igraph layout
    visIgraphLayout(layout = "layout_with_graphopt", randomSeed = 0)
}
