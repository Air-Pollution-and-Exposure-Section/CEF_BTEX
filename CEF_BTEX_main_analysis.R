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

# list of btex species
voc_species = c('Dichloromethane', 'Hexane', 'Chloroform', '__2_Dichloroethane', 'Benzene', 'Trichloroethylene',
                'Toluene', 'Tetrachloroethylene', 'Ethylbenzene', '_m_p__Xylene', 'o_Xylene', 'Styrene', 'Cumene', 
                'a_Pinene', '__1_2_2_Tetrchloroethane', 'n_Decane', '__3_5_Trimethylbenzene', '__2_4_Trimethylbenzene', 
                'Pentachloroethane', 'd_Limonene', 'p_Cymene', '__3_Dichlorobenzene', '__4_Dichlorobenzene', 'Hexachloroethane',
                '__2_4_Trichlorobenzene', 'Naphthalene')


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
  dplyr::select(species, season, N, n_BDL, `%_BDL`, min, max, mean, sd, p25, p50, p75)

# generate and modify table
summary_stats_flextable = generate_summary_statistics(voc_summary_stats, voc_species)
summary_stats_flextable = summary_stats_flextable %>%
  flextable::width(x = ., j = 3, unit = 'mm', width = 1) %>%
  flextable::width(x = ., j = c(1,2), unit = 'mm', width = c(40,17)) %>%
  flextable::width(x = ., j = c(4,5,6,7,8,9,10, 11, 12, 13),
                   width = c(10, 13, 13, 10, 10, 10, 10, 15, 10, 10), unit = 'mm') %>%
  flextable::fontsize(x = ., size = 8, part = c('all'))

# save table
summary_stats_flextable %>% flextable::save_as_docx(., path = 'tables/summary_stats.docx')



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
seasonal_comparison_plot(data = voc_long, stat_res = voc_season_dunn)
