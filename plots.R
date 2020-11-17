plot_chromatograms = function(identifier, values, tidy, map,
                              mode = c('Scatterplot', 'Heatmap'),
                              log_transform = FALSE) {
  # subset the map
  map0 = filter(map, identifier == !!identifier, id %in% values)
  # link to chromatograms
  chroms = filter(tidy, uniprot %in% map0$uniprot) %>%
    left_join(map0, by = 'uniprot')
  # flag protein group AND query identifier value
  chroms %<>%
    mutate(group = paste0(protein_group, ' (', id, ')'))
  
  # optionally, log-transform
  if (log_transform)
    chroms %<>% mutate(intensity = log(intensity))
  
  # create facets
  chroms %<>%
    mutate(facet = paste0(accession, ": ", replicate, "\n", species))
  
  # plot 
  size_sm = 6
  size_lg = 7
  if (tolower(mode) == 'scatterplot') {
    pal = pals::polychrome(n = 36) %>% extract(-2) %>% unname()
    pal = pals::alphabet() %>% unname()
    p = chroms %>%
      ggplot(aes(x = fraction, y = intensity, color = group)) +
      facet_wrap(~ facet, ncol = 4, scales = 'free') +
      geom_point(size = 0.5) + 
      geom_line() +
      scale_color_manual('Protein group', values = pal) +
      scale_x_continuous('Fraction') +
      scale_y_continuous('Intensity (iBAQ)') +
      guides(color = guide_legend(override.aes = list(size = 0.9),
                                  ncol = 1)) +
      theme_bw() +
      theme(axis.title.x = element_text(size = size_lg), 
            axis.title.y = element_text(size = size_lg),
            axis.text.x = element_text(size = size_sm), 
            axis.text.y = element_text(size = size_sm),
            panel.grid = element_blank(), 
            strip.text = element_text(size = size_lg),
            strip.background = element_blank(), 
            axis.line.y = element_blank(),
            axis.line.x = element_blank(), 
            legend.position = "top", 
            legend.text = element_text(size = size_sm), 
            legend.title = element_text(size = size_sm), 
            legend.key.size = unit(0.25, "lines"), 
            legend.margin = margin(rep(0, 4)), 
            legend.background = element_blank(),
            plot.title = element_text(size = size_lg, hjust = 0.5),
            aspect.ratio = 0.75)
    p
  } else if (tolower(mode) == 'heatmap') {
    chroms0 = chroms %>%
      group_by(accession, replicate, protein_group) %>%
      mutate(
        intensity_w = ifelse(intensity > quantile(intensity, probs = 0.95, na.rm = TRUE),
                             quantile(intensity, probs = 0.95, na.rm = TRUE),
                             ifelse(intensity < quantile(intensity, probs = 0.05, na.rm = TRUE),
                                    quantile(intensity, probs = 0.05, na.rm = TRUE),
                                    intensity)),
        norm = scale(intensity_w)) %>%
      ungroup()
    range = range(chroms0$norm, na.rm = TRUE) %>% round(digits = 1)
    p = chroms0 %>%
      ggplot(aes(x = fraction, y = group, fill = norm)) +
      facet_wrap(~ facet, ncol = 2, scales = 'free') +
      geom_tile() +
      scale_fill_distiller(palette = 'RdBu', na.value = 'grey94',
                           name = 'Scaled abundance  ',
                           breaks = range,
                           limits = range) + 
      scale_x_continuous('Fraction', expand = c(0, 0)) +
      scale_y_discrete('Protein', expand = c(0, 0)) +
      guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
      theme_bw() +
      theme(axis.title.x = element_text(size = size_lg), 
            axis.title.y = element_text(size = size_lg),
            axis.text.x = element_text(size = size_sm), 
            axis.text.y = element_text(size = size_sm),
            panel.grid = element_blank(), 
            strip.text = element_text(size = size_lg),
            strip.background = element_blank(), 
            axis.line.y = element_blank(),
            axis.line.x = element_blank(), 
            legend.position = "top", 
            legend.text = element_text(size = size_sm), 
            legend.title = element_text(size = size_sm), 
            legend.key.size = unit(0.25, "lines"), 
            legend.margin = margin(rep(0, 4)), 
            legend.background = element_blank(),
            plot.title = element_text(size = size_lg, hjust = 0.5))
    p
  }
  
  print(p)
}