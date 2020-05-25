library(tidyverse)

wide_data <- read_csv(here::here('cleaned_numeric_data_wide.csv'), col_types = cols(Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor(),
                                                                                    keyword_batch_comp = col_factor()))

dep_var <- names(wide_data[,68:86])

supplements <- tibble(
                  output_file = stringr::str_c(dep_var, "_supplement.html"),
                  params = map(dep_var, ~list(dep_var = .))
)

supplements %>%
  pwalk(rmarkdown::render, 'dv_supplements.Rmd')

  