# Phill Pham
# 2024-06-06


# Load libraries + source functions
library(tidyverse)
library(stringr)
library(flextable)
library(moments)
library(rstatix)
library(ggpubr)
library(readxl)
library(rstatix)

source('CEF_BTEX_functions.R')


# Loading Data ------------------------------------------------------------

# load data & remove data where voc data is missing (can easily do this by filtering out any VOC species by na values)
# Note - units assumed to be ug/m3
voc = readxl::read_xlsx(path = 'data/EGROW_passive_for_analysis.xlsx', sheet = 2, col_names = T) %>%
  dplyr::filter(!is.na(Dichloromethane))

# list of btex species - apply name transformations
voc_species = name_transform(c('Dichloromethane', 'Hexane', 'Chloroform', '__2_Dichloroethane', 'Benzene', 'Trichloroethylene',
                'Toluene', 'Tetrachloroethylene', 'Ethylbenzene', '_m_p__Xylene', 'o_Xylene', 'Styrene', 'Cumene', 
                'a_Pinene', '__1_2_2_Tetrchloroethane', 'n_Decane', '__3_5_Trimethylbenzene', '__2_4_Trimethylbenzene', 
                'Pentachloroethane', 'd_Limonene', 'p_Cymene', '__3_Dichlorobenzene', '__4_Dichlorobenzene', 'Hexachloroethane',
                '__2_4_Trichlorobenzene', 'Naphthalene'))

names(voc_species) = c('Dichloromethane', 'Hexane', 'Chloroform', '__2_Dichloroethane', 'Benzene', 'Trichloroethylene',
                       'Toluene', 'Tetrachloroethylene', 'Ethylbenzene', '_m_p__Xylene', 'o_Xylene', 'Styrene', 'Cumene', 
                       'a_Pinene', '__1_2_2_Tetrchloroethane', 'n_Decane', '__3_5_Trimethylbenzene', '__2_4_Trimethylbenzene', 
                       'Pentachloroethane', 'd_Limonene', 'p_Cymene', '__3_Dichlorobenzene', '__4_Dichlorobenzene', 'Hexachloroethane',
                       '__2_4_Trichlorobenzene', 'Naphthalene')


# generate list of variables (VOCs) to be renamed
rename_dict = c(voc_species, paste0(voc_species,'_BDL'))
names(rename_dict) = c(names(voc_species), paste0(names(voc_species), '_BDL'))

# identify which names to replace
voc_names_to_replace = names(voc) %in% names(rename_dict)

# replace names by matching key-value pairs (name : value)
names(voc)[voc_names_to_replace] = rename_dict[names(voc)[voc_names_to_replace]]


# simplify voc_species since dictionary matching no longer needed
voc_species = name_transform(c('Dichloromethane', 'Hexane', 'Chloroform', '__2_Dichloroethane', 'Benzene', 'Trichloroethylene',
                               'Toluene', 'Tetrachloroethylene', 'Ethylbenzene', '_m_p__Xylene', 'o_Xylene', 'Styrene', 'Cumene', 
                               'a_Pinene', '__1_2_2_Tetrchloroethane', 'n_Decane', '__3_5_Trimethylbenzene', '__2_4_Trimethylbenzene', 
                               'Pentachloroethane', 'd_Limonene', 'p_Cymene', '__3_Dichlorobenzene', '__4_Dichlorobenzene', 'Hexachloroethane',
                               '__2_4_Trichlorobenzene', 'Naphthalene'))

# convert data to long format - for calculations and plotting later
voc_long = voc %>%
  dplyr::select(season, tidyselect::all_of(voc_species)) %>%
  tidyr::pivot_longer(cols = tidyselect::all_of(voc_species), 
                      names_to = 'species', 
                      values_to = 'conc')

# examine data distribution - does not look normally distributed. Should use Dunn Test or other non parametric approaches
# data is also BDL-skewed; dataset is set such that < BDL is set to BDL
voc_long %>%
  ggplot(aes(x = conc, y = after_stat(density))) +
  geom_density() +
  facet_wrap(~species) +
  coord_cartesian(xlim = c(0, 0.5)) +
  theme_bw()


# Data Wrangling ----------------------------------------------------------

# compute number of below detection limits
voc_BDL = voc %>% dplyr::select(season, contains("BDL")) %>%
  tidyr::pivot_longer(cols = all_of(contains("BDL")), names_to = 'species', values_to = 'flag') %>%
  dplyr::mutate(species = stringr::str_remove(string = species, pattern = '_BDL')) %>%
  dplyr::group_by(season, species) %>%
  dplyr::summarise(N = n(),
                   n_BDL = sum(flag %in% c('1', '2')),
                   `%_BDL` = (n_BDL/N) * 100) %>%
  dplyr::select(season, species, `%_BDL`, n_BDL)

# compute summary statstics: N, mean, SD, min, p25, p50/median, p75, max from voc_long
voc_summary_stats = voc_long %>%
  dplyr::group_by(season, species) %>%
  dplyr::summarise(N = n(),
                   mean = mean(conc),
                   sd = sd(conc),
                   min = min(conc),
                   p25 = quantile(conc, probs = 0.25),
                   p50 = quantile(conc, probs = 0.5),
                   p75 = quantile(conc, probs = 0.75),
                   max = max(conc)) %>%
  dplyr::left_join(., voc_BDL, by = c('species', 'season')) %>%
  dplyr::arrange(species, season) %>%
  dplyr::select(species, season, N, n_BDL, `%_BDL`, min, max, mean, sd, p25, p50, p75) %>%
  dplyr::ungroup()

# generate and modify table
summary_stats_flextable = generate_summary_statistics(voc_summary_stats, voc_species)
summary_stats_flextable = summary_stats_flextable %>%
  flextable::width(x = ., j = 3, unit = 'mm', width = 1) %>%
  flextable::width(x = ., j = c(1,2), unit = 'mm', width = c(35,15)) %>%
  flextable::width(x = ., j = c(4,5,6,7,8,9,10, 11, 12, 13),
                   width = c(10, 13, 13, 10, 10, 10, 10, 15, 10, 10), unit = 'mm') %>%
  flextable::fontsize(x = ., size = 8, part = c('all'))

# save table (if doesn't exist)
if(!file.exists('tables/summary_stats.docx')){
  summary_stats_flextable %>% flextable::save_as_docx(., path = 'tables/summary_stats.docx')
}


#------------------------------------------------------------#
#### Generate separate data frame with >80% BDL set to NA ####
#------------------------------------------------------------#

# General Rationale - the results where (nearly) everything in a season is below detection limits are not interesting; 
# this can be removed from visualizations


# run conditional mutate to screen for >80% data bdl. Data table already stratified species/season
voc_species_bdl = voc_summary_stats %>% 
  dplyr::filter(`%_BDL` > 80) %>% 
  dplyr::select(species, season) %>% 
  unique()

# set values to -999 to not render
voc_long_hide_BDL = voc_long
matches = merge(voc_long, voc_species_bdl, by = c('season', 'species'))
voc_long_hide_BDL$conc[with(voc_long_hide_BDL, paste(season, species)) %in% with(matches, paste(season, species))] <- -999

# set values to NA to 'remove' in calculation or lm plotting
voc_long_NA_BDL = voc_long
matches = merge(voc_long, voc_species_bdl, by = c('season', 'species'))
voc_long_NA_BDL$conc[with(voc_long_NA_BDL, paste(season, species)) %in% with(matches, paste(season, species))] <- -NA

####



# Effect of season on BTEX species concentration --------------------------

# note that only 13 species had values above detection limits
# species where all measurements were below detection limit were automatically excluded from analysis (due to 0 difference in means)
voc_season_dunn = voc_long %>%
  dplyr::group_by(species) %>%
  rstatix::dunn_test(data = ., formula = conc ~ season, detailed = TRUE) %>%
  rstatix::add_xy_position() %>%
  dplyr::select(-groups)

# some adjustments for plotting stats comparisons
voc_season_dunn = voc_season_dunn %>%
  dplyr::mutate(x = seq(1:nrow(.)),
                xmin = x - 0.4,
                xmax = x + 0.4,
                group1 = 'fall',
                group2 = 'winter')

# we will only look at BTEX species where concentrations were above detection limit (list from voc_season_dunn)
seasonal_comparison_boxplot = seasonal_comparison_plot(data = voc_long, stat_res = voc_season_dunn)
if(!file.exists('figures/seasonal_comparison_dunn.png')){
  png(filename = 'figures/seasonal_comparison_dunn.png', width = 8, height = 6, units = 'in', res = 600)
  seasonal_comparison_boxplot
  dev.off()
}



# some adjustments for plotting stats comparisons
voc_season_dunn_hide_BDL = voc_season_dunn %>%
  dplyr::filter(!species %in% c("Chloroform", "Trichloroethylene", "d_Limonene")) %>% 
  dplyr::mutate(x = seq(1:nrow(.)),
                xmin = x - 0.4,
                xmax = x + 0.4,
                group1 = 'fall',
                group2 = 'winter')

season_comparison_boxplot_na_rm = seasonal_comparison_plot(data = voc_long_hide_BDL, stat_res = voc_season_dunn_hide_BDL)
if(!file.exists('figures/seasonal_comparison_dunn_hide80BDL.png')){
  png(filename = 'figures/seasonal_comparison_dunn_hide80BDL.png', width = 8, height = 6, units = 'in', res = 600)
  season_comparison_boxplot_na_rm
  dev.off()
}


# Effect of distance to gas station ---------------------------------------

# select only compounds where % BDL is not 100%
voc_detectable = voc_summary_stats %>%
  dplyr::filter(`%_BDL` < 100) %>%
  dplyr::select(species) %>%
  unique() %>%
  unlist(use.names = FALSE)
  

# examine the effect of gas stations on BTEX species concentrations
# concentration ~ dist_to_closest_gas_m

seasons = c('fall', 'winter')
voc_seasons = expand.grid(voc_detectable, seasons)
voc_seasons = paste(voc_seasons$Var1, voc_seasons$Var2, sep = '_')

##** NOTE - We may want to consider taking all regression results and running a P-value correction since this can be viewed as multiple testing?**##
lm_results = vector(mode = 'list', length = (length(voc_detectable) * length(seasons)))
lm_plots = vector(mode = 'list', length = (length(voc_detectable) * length(seasons)))

for(i in 1:length(voc_detectable)){
  for(j in 1:length(seasons)){
    
    index = (i-1)*length(seasons) + j
    names(lm_results)[index] = paste(voc_detectable[i], seasons[j], sep = '_')
    data_subset = voc %>% dplyr::filter(season == seasons[j])
    form = as.formula(paste(paste0('`', voc_detectable[i], '`'), '~ dist_to_closest_gas_m'))
    
    # run lm and save results to list
    lm_results[[index]]$season = seasons[j]
    lm_results[[index]]$voc = voc_detectable[i]
    lm_results[[index]]$lm_result = stats::lm(formula = form, data = data_subset)
    lm_results[[index]]$lm_summary = summary(lm_results[[index]]$lm_result)
    lm_results[[index]]$slope = lm_results[[index]]$lm_summary$coefficients[2,1]
    lm_results[[index]]$slope_signif = lm_results[[index]]$lm_summary$coefficients[2,4]
    
    ## Graphical Representation of Linear Regression
    
    # plot vars
    y_pos = max(data_subset[,voc_detectable[i]]) * 1.05
    y_name = bquote(paste("Concentration", " (", mu * g, "/", m ^ 3, ")"))
    plot_title = paste(stringr::str_to_sentence(seasons[j]), voc_detectable[i])
    
    eqn = lm_equation(lm_results[[index]]$lm_result)
    
    lm_plots[[index]] = ggplot2::ggplot(data = data_subset, 
                                        aes(x = dist_to_closest_gas_m,
                                            y = !!rlang::sym(voc_detectable[i])))+
      ggplot2::stat_smooth(method = 'lm',
                           colour = 'black', alpha = 0.2) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(name = "Distance to Closest Gas Station (m)") +
      ggplot2::scale_y_continuous(name = y_name) +
      ggplot2::annotate(geom = 'text', x = 0, y = y_pos, label = eqn, parse = TRUE, hjust = 0) +
      ggplot2::ggtitle(label = plot_title) +
      ggpubr::theme_pubr()
    lm_plots[[index]]
  }
}

lm_plots_combined = vector(mode = 'list', length = length(voc_detectable))
for(i in 1:length(voc_detectable)){
  y_pos_fall = max(voc[,voc_detectable[i]]) * 1.1
  y_pos_winter = max(voc[,voc_detectable[i]]) * 1.05
  y_name = bquote(paste("Concentration", " (", mu * g, "/", m ^ 3, ")"))
  plot_title = voc_detectable[i]
  
  form = as.formula(paste(paste0('`', voc_detectable[i], '`'), '~ dist_to_closest_gas_m'))
  fall_regression = lm(formula = form, data = voc, subset = season == 'fall')
  winter_regression = lm(formula = form, data = voc, subset = season == 'winter')
  eqn_fall = lm_equation(fall_regression)
  eqn_winter = lm_equation(winter_regression)
  
  lm_plots_combined[[i]] = ggplot2::ggplot(data = voc, 
                                           aes(x = dist_to_closest_gas_m,
                                               y = !!rlang::sym(voc_detectable[i]),
                                               color = season,
                                               fill = season)) +
    ggplot2::geom_point() +
    # add regression lines
    ggplot2::stat_smooth(method = 'lm', 
                         alpha = 0.1,
                         inherit.aes = TRUE, linewidth = 1) +
    # adjust colour/fill scales
    ggplot2::scale_fill_manual(name = 'Season', values = c('tomato3', 'skyblue3'), labels = c('Fall', 'Winter')) +
    ggplot2::scale_color_manual(name = 'Season', values = c('tomato3', 'skyblue3'), labels = c('Fall', 'Winter')) +
    ggplot2::scale_x_continuous(name = "Distance to Closest Gas Station (m)") +
    ggplot2::scale_y_continuous(name = y_name) +
    ggplot2::annotate(geom = 'text', x = 0, y = y_pos_fall, label = eqn_fall, parse = TRUE, hjust = 0, color = 'tomato4') +
    ggplot2::annotate(geom = 'text', x = 0, y = y_pos_winter, label = eqn_winter, parse = TRUE, hjust = 0, color = 'skyblue4') +
    ggplot2::ggtitle(label = plot_title) +
    ggpubr::theme_pubr()
  lm_plots_combined[[i]]
}



# Second set of linear regressions - only render > 80% detectable ---------

# in other words, drop the species-season combinations where > 80% of samples are below detection limit

# tabulate linear regression results
lm_results_df = data.frame(
  species = character(length(lm_results)),
  season = character(length(lm_results)),
  y_intercept = numeric(length(lm_results)),
  y_intercept_sd = numeric(length(lm_results)),
  y_intercept_p_value = numeric(length(lm_results)),
  slope = numeric(length(lm_results)),
  slope_sd = numeric(length(lm_results)),
  slope_p_value = numeric(length(lm_results))
)
for(i in 1:length(lm_results)){
  lm_results_df[i,1] = lm_results[[i]]$voc
  lm_results_df[i,2] = lm_results[[i]]$season
  lm_results_df[i,3] = as.numeric(lm_results[[i]]$lm_summary$coefficients[1,1])
  lm_results_df[i,4] = as.numeric(lm_results[[i]]$lm_summary$coefficients[1,2])
  lm_results_df[i,5] = as.numeric(lm_results[[i]]$lm_summary$coefficients[1,4])
  lm_results_df[i,6] = as.numeric(lm_results[[i]]$lm_summary$coefficients[2,1])
  lm_results_df[i,7] = as.numeric(lm_results[[i]]$lm_summary$coefficients[2,2])
  lm_results_df[i,8] = as.numeric(lm_results[[i]]$lm_summary$coefficients[2,4])
}

# all significant slopes - relationships primarily occur in winter months
lm_results_df %>% 
  dplyr::filter(slope_p_value <= 0.05)

selected_species = c("(m,p)_Xylene", "Benzene", "Ethylbenzene", "Hexane", "o_Xylene")
names(selected_species) = c("(m,p)-Xylene", "Benzene", "Ethylbenzene", "Hexane", "o-Xylene")

lm_plots_combined_detectable80 = vector(mode = 'list', length = length(selected_species))
for(i in 1:length(selected_species)){
  y_pos_fall = max(voc[,selected_species[i]]) * 1.1
  y_pos_winter = max(voc[,selected_species[i]]) * 1.05
  y_name = bquote(paste("Concentration", " (", mu * g, "/", m ^ 3, ")"))
  plot_title = selected_species[i]
  
  form = as.formula(paste(paste0('`', selected_species[i], '`'), '~ dist_to_closest_gas_m'))
  fall_regression = lm(formula = form, data = voc, subset = season == 'fall')
  winter_regression = lm(formula = form, data = voc, subset = season == 'winter')
  eqn_fall = lm_equation(fall_regression)
  eqn_winter = lm_equation(winter_regression)
  lm_plots_combined_detectable80[[i]] = ggplot2::ggplot(data = voc, 
                                           aes(x = dist_to_closest_gas_m,
                                               y = !!rlang::sym(selected_species[i]),
                                               color = season,
                                               fill = season)) +
    ggplot2::geom_point() +
    # add regression lines
    ggplot2::stat_smooth(method = 'lm', 
                         alpha = 0.1,
                         inherit.aes = TRUE, linewidth = 1) +
    # adjust colour/fill scales
    ggplot2::scale_fill_manual(name = 'Season', values = c('tomato3', 'skyblue3'), labels = c('Fall', 'Winter')) +
    ggplot2::scale_color_manual(name = 'Season', values = c('tomato3', 'skyblue3'), labels = c('Fall', 'Winter')) +
    ggplot2::scale_x_continuous(name = "Distance to Closest Gas Station (m)") +
    ggplot2::scale_y_continuous(name = y_name) +
    ggplot2::annotate(geom = 'text', x = 0, y = y_pos_fall, label = eqn_fall, parse = TRUE, hjust = 0, color = 'tomato4') +
    ggplot2::annotate(geom = 'text', x = 0, y = y_pos_winter, label = eqn_winter, parse = TRUE, hjust = 0, color = 'skyblue4') +
    ggplot2::ggtitle(label = plot_title) +
    ggpubr::theme_pubr()
  
  
}

# save results - use the first plot as a check if plots have already been generated
if(!file.exists(stringr::str_glue("figures/{selected_species[1]}_linear_regression.png"))){
  for(i in 1:length(selected_species)){
    file_name = stringr::str_glue("figures/{selected_species[i]}_linear_regression.png")
    plot = lm_plots_combined_detectable80[[i]] +
      ggtitle(label = names(selected_species[i])) +
      theme(legend.position = 'right')
    
    png(filename = file_name, width = 8, height = 6, res = 600, units = 'in')
    print(plot)
    dev.off()
  }
}

# Table of significant linear regressions --------------------------------

# table configuration for flextable
options(scipen=999) # want to display all digits
lm_signif_table = lm_results_df %>% 
  dplyr::filter(species %in% selected_species) %>% 
  dplyr::select(species, season, slope, slope_sd, slope_p_value) %>% 
  dplyr::mutate(species = factor(species, levels = selected_species),
                season = factor(season, levels = c('fall', 'winter')),
                slope_round = ifelse(round(slope, 4) == 0, '<0.0001', as.character(round(slope,4))),
                slope_sd_round = ifelse(round(slope_sd, 4) == 0, '<0.0001', as.character(round(slope,4))),
                slope_and_sd = as.character(stringr::str_glue("{slope_round} ± {slope_sd_round}")),
                p_value = ifelse(round(slope_p_value, 4) == 0, '<0.0001', as.character(round(slope_p_value, 4)))) %>% 
  dplyr::select(species, season, slope_and_sd, p_value) %>% 
  as.data.frame()

# generate the flextable
lm_signif_ft = lm_signif_table %>% 
  flextable::flextable() %>% 
  flextable::merge_v(., j='Species', ) %>% 
  flextable::fix_border_issues(., part = 'all') %>% 
  flextable::width(., width = c(1,1,2,1)) %>% 
  flextable::valign(., j = 1, valign = 'top') %>% 
  flextable::set_header_labels(
    species = "Species",
    season = "Season",
    slope_and_sd = "Slope ± SD",
    p_value = "P-value"
  ) %>% 
  flextable::fontsize(x = ., size = 10) %>% 
  # building the superscripts for column headers Slope +- SD with units
  flextable::compose(
    part = 'header',
    j = "slope_and_sd",
    value = as_paragraph(
      "Slope ± SD (µg/m", 
      as_sup("3"),       # Superscript for 3
      " m", 
      as_sup("-1"),      # Superscript for -1
      ")"
    )
  )

if(!file.exists('tables/lm_results_signif.docx')){
  lm_signif_ft %>% flextable::save_as_docx(., path = 'tables/lm_results_signif.docx')
}
