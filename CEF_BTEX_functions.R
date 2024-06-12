# all custom functions used with CEF_BTEX

hello_world = function(){
  print('Hello_world!')
}


#' Generate Summary Statistics
#'
#' @param data, tibble/dataframe with columns c(mean, sd, N, n_BDL, %_BDL, min, p25, p50, p75, max)
#' @return flextable object
generate_summary_statistics = function(data, species_levels){
  require(flextable)
  require(dplyr)
  
  # set season and species to factor
  # note min and max are some kind of special keywords... you will get a bug by naming your variables min or max
  data = data %>% dplyr::mutate(stat = 'Statistic',
                                species = factor(species, levels = species_levels),
                                season = factor(season, levels = c('fall', 'winter')),
                                min_val = as.numeric(min),
                                max_val = as.numeric(max)) %>%
    dplyr::rename(Season =season, Species = species)
  
  # use tabulate to set up aggregation
  tab = data %>% 
    flextable::tabulator(x = .,
                         rows = c('Species', 'Season'),
                         columns = 'stat',
                         n = flextable::as_paragraph(N),
                         `n BDL` = flextable::as_paragraph(n_BDL),
                         `% BDL` = flextable::as_paragraph(`%_BDL`),
                         `μ` = flextable::as_paragraph(mean),
                         `σ` = flextable::as_paragraph(sd),
                         min = flextable::as_paragraph(as.numeric(min_val)),
                         Q1 = flextable::as_paragraph(p25),
                         median = flextable::as_paragraph(p50),
                         Q3 = flextable::as_paragraph(p75),
                         max = flextable::as_paragraph(as.numeric(max_val)))
  
  # return as flextable object, ready to render
  return(flextable::as_flextable(tab))
}



seasonal_comparison_plot = function(data, stat_res, sigma = 0.01, logbase = 10) {
  require(ggpubr)
  require(ggplot2)
  require(rstatix)
  
  y_name = bquote(paste("Concentration", " (", n * g, "/", m ^ 3, ")"))
  y_range = c(0,3)
  y_breaks = c(0.1, 0.5, 1, 2, 3)
  y_breaks_minor = c(seq(0.1, 5, 0.1))
  
  stat_res$y.position = asinh(stat_res$y.position / (2 * sigma)) / log(logbase)
  
  plot = data %>%
    dplyr::filter(species %in% stat_res$species) %>%
    dplyr::mutate(species = factor(species, levels = stat_res$species)) %>%
    ggpubr::ggboxplot(data = .,
                      x = 'species',
                      y = 'conc',
                      fill = 'season') +
    ggplot2::scale_y_continuous(name = y_name,
                                trans = scales::pseudo_log_trans(sigma = 0.01, base = 10),
                                breaks = y_breaks,
                                minor_breaks = y_breaks_minor, labels = y_breaks) +
    ggplot2::scale_x_discrete(name = 'Species') +
    ggplot2::scale_fill_manual(name = 'Season', values = c('skyblue3', 'tomato3'), labels = c('Fall', 'Winter')) +
    ggpubr::stat_pvalue_manual(data = stat_res, tip.length = 0) +
    ggplot2::theme(legend.position = 'right',
                   panel.grid.major.y = element_line(colour = 'grey75'),
                   panel.grid.minor.y = element_line(colour = 'grey95'),
                   axis.text.x = element_text(
                     angle = 45,
                     vjust = 1,
                     hjust = 1,
                     size = 8
                   ), 
                   axis.text.x.bottom = element_text(size = 8),
                   axis.text.y.left = element_text(size = 8),
                   text = element_text(size = 8))
  
  return(plot)
}
