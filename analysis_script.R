# load libraries
library(brms)
library(tidyverse)
library(ggplot2)

# load data
long_data <- read_csv(here::here('cleaned_numeric_data_long.csv'), col_types = cols(article_type = col_factor(),
                                                                                    Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor()))
wide_data <- read_csv(here::here('cleaned_numeric_data_wide.csv'), col_types = cols(Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor()))

options(contrasts = rep("contr.sum", 2))

## check distribution of question quality variable
ggplot(long_data %>% filter(grepl('QuestionQuality', question)), aes(x = response)) +
  geom_histogram(breaks=seq(-5, 4, by = 1), col = 'blue', fill = 'grey') +
  scale_x_continuous(breaks=seq(-5, 4, by = 1))


# bayesian models for single DV
priors <- c(set_prior("normal(0,2)", "b"),
            set_prior("normal(0, 2.5)", "sd"))

### within-subjects models

# difference score model
within_model_diffs <- brm(diff_question_quality ~ Field + Order + Match +
                                                    Order*Match +
                                                    (1|RR),
                           data = wide_data,
                           prior = priors, 
                           family = 'gaussian',
                           chains = 4)

summary(within_model_diffs)
plot(within_model)
pp_check(within_model_diffs)
pp_check(within_model_diffs, type = "stat", stat = 'median', nsamples = 100)
WAIC(within_model_diffs)
loo(within_model_diffs)

# no different score model
within_model <- brm(response ~ Field + article_type + Order + Match +
                                  article_type*Order + article_type*Match + Order*Match +
                                  article_type*Order*Match + 
                                          (1|participant_id) + (article_type|RR),
                     data = long_data %>% filter(grepl('QuestionQuality', question)),
                     prior = priors,
                     family = 'gaussian',
                     chains = 4)

summary(within_model)
plot(within_model)
pp_check(within_model)
pp_check(within_model, type = "stat", stat = 'median', nsamples = 100)
WAIC(within_model)
loo(within_model)

### between-subjects models

between_model <- brm(response ~ Field + article_type + Match + article_type*Match +
                                            (article_type|RR),
                      data = long_data %>% filter(grepl('QuestionQuality', question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                      prior = priors,
                      family = 'gaussian',
                      chains = 4)

summary(between_model)
pp_check(between_model)
WAIC(between_model)
loo(between_model)




