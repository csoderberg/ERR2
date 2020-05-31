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


# supplements by variable for between data
long_data <- read_csv(here::here('cleaned_numeric_data_long.csv'), col_types = cols(article_type = col_factor(),
                                                                                    Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor(),
                                                                                    keyword_batch_comp = col_factor())) %>%
  mutate(question = case_when(article_type == 'RR' ~ str_sub(question, 3),
                              article_type == 'nonRR' ~ str_sub(question, 4)))

btw_dep_var <- long_data %>% select(question) %>% distinct(question) %>% pull(question)

btw_supplements <- tibble(
  output_file = stringr::str_c(btw_dep_var, "btwsubj_supplement.html"),
  params = map(btw_dep_var, ~list(btw_dep_var = .))
)

btw_supplements %>%
  pwalk(rmarkdown::render, 'dv_supplements_between.Rmd')
