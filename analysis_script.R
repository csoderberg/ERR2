# load libraries
library(brms)
library(tidyverse)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(broom)
library(dotwhisker)

# load data
long_data <- read_csv(here::here('cleaned_numeric_data_long.csv'), col_types = cols(article_type = col_factor(),
                                                                                    Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor(),
                                                                                    keyword_batch_comp = col_factor()))
wide_data <- read_csv(here::here('cleaned_numeric_data_wide.csv'), col_types = cols(Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor(),
                                                                                    keyword_batch_comp = col_factor()))

# set up contrasts codes
contrasts(wide_data$Field) <- contr.sum(3)
contrasts(wide_data$keyword_batch_comp) <- contr.sum(2)
contrasts(wide_data$Order) <- contr.sum(2)
contrasts(wide_data$Match) <- contr.sum(2)
contrasts(long_data$article_type) <- contr.sum(2)
contrasts(long_data$Field) <- contr.sum(3)
contrasts(long_data$keyword_batch_comp) <- contr.sum(2)
contrasts(long_data$Order) <- contr.sum(2)
contrasts(long_data$Match) <- contr.sum(2)

## check distribution of question quality variable
ggplot(long_data %>% filter(grepl('QuestionQuality', question)), aes(x = response)) +
  geom_histogram(breaks=seq(-5, 4, by = 1), col = 'blue', fill = 'grey') +
  scale_x_continuous(breaks=seq(-5, 4, by = 1))


# bayesian models for single DV
priors <- c(set_prior("normal(0,2)", "b"),
            set_prior("normal(0, 2.5)", "sd"))




### functions for within models

# set up function to run within difference pooled models on all dvs
within_diff_pooled_model <- function(dv, set_priors) {
  within_model_diffs <- brm(as.formula(paste(dv, "~ Field + keyword_batch_comp + Order + Match +
                                                    Order*Match +
                                                    (1|RR)")),
                            data = wide_data,
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_diffs)
}

# set up function to run within difference model on all dvs for those in first batch
within_diff_keywords1_model <- function(dv, set_priors) {
  within_model_keywords1 <- brm(as.formula(paste(dv, "~ Field + Order + Match +
                                                    Order*Match +
                                                    (1|RR)")),
                            data = wide_data %>% filter(keyword_batch_comp == 1),
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_keywords1)
}

# set up function to run within difference model on all dvs for those in 2+3 batch
within_diff_keywords2_model <- function(dv, set_priors) {
  within_model_keywords2 <- brm(as.formula(paste(dv, "~ Field + Order + Match +
                                                    Order*Match +
                                                    (1|RR)")),
                                data = wide_data %>% filter(keyword_batch_comp == 2),
                                prior = set_priors, 
                                family = 'gaussian',
                                chains = 4)
  return(within_model_keywords2)
}

# function to create graph of posterior samples for intercept
create_posteriors <- function(results, variable){
  posteriors <- list(posterior_samples(results) %>%
                       select('b_Intercept') %>%
                       rename(setNames('b_Intercept', variable)))
  
  return(posteriors)
}

# Set up which model/prior/dv combinations to run
within_models <- crossing(dv = names(wide_data[,68:86]),
                    set_priors = c(list(priors))) %>%
  mutate(within_pooled_model_results = pmap(list(dv, set_priors), within_diff_pooled_model)) %>%
  mutate(posteriors = pmap(list(within_pooled_model_results, variable = dv), create_posteriors))


# get all intercepts into wide format for graphing
intercepts <- c()

for (i in 1:nrow(within_models)) {
  intercepts <- bind_cols(intercepts, within_models$posteriors[[i]])
}

mcmc_areas(intercepts,
           prob=.95)

mcmc_intervals(intercepts, prob = .95)



### functions for between subjects models

# between-subjects models comparing pooled, batch 1, batch 2 + 3
between_pooled_model <- function(dv, set_priors) {
  between_model <- brm(response ~ Field + keyword_batch_comp + article_type + Match + article_type*Match +
                         (article_type|RR),
                       data = long_data %>% filter(grepl(as.string(dv), question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                       prior = priors,
                       family = 'gaussian',
                       chains = 4)
  return(between_model)
  }


# between-subjects model for batch 1
between_keywords1_model <- function(dv, set_priors) {
  between_keywords1 <- brm(response ~ Field + article_type + Match + article_type*Match +
                                   (article_type|RR),
                                 data = long_data %>% 
                                   filter(keyword_batch_comp == 1) %>%
                                   filter(grepl(as.string(dv), question)) %>% 
                                   filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                                 prior = priors,
                                 family = 'gaussian',
                                 chains = 4)
  return(between_keywords1)
}

# differenc score model for batched 2+3
between_keywords2_model <- function(dv, set_priors) {
  between_keywords2 <- brm(response ~ Field + article_type + Match + article_type*Match +
                                   (article_type|RR),
                                 data = long_data %>% 
                                   filter(keyword_batch_comp == 2) %>%
                                   filter(grepl(as.string(dv), question)) %>% 
                                   filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                                 prior = priors,
                                 family = 'gaussian',
                                 chains = 4)
}

# Set up which model/prior/dv combinations to run for between models
between_models <- crossing(dv = c(long_data %>% select(question) %>% distinct(question)),
                          set_priors = c(list(priors))) %>%
  mutate(between_pooled_model_results = pmap(list(dv, set_priors), within_diff_pooled_model)) %>%
  mutate(posteriors = pmap(list(between_pooled_model_results, variable = dv), create_posteriors))



# graphical comparisons
between_posteriors_keywords2 <- suppressMessages( 
  mcmc_areas(posterior_samples(between_model_keywords2),
             regex_pars = "b_",
             prob=.9) +
    xlim(-2, 2) +
    labs(title = 'Batch 2')
)

between_posteriors_keywords1 <- suppressMessages( 
  mcmc_areas(posterior_samples(between_model_keywords1),
             regex_pars = "b_",
             prob=.9) +
    xlim(-2, 2) +
    labs(title = 'Batch 1')
)

posteriors_pooled <- suppressMessages( 
  mcmc_areas(posterior_samples(between_model),
             regex_pars = "b_",
             prob=.9) +
    xlim(-2, 2) +
    labs(title = 'Pooled Batched')
)

bind_rows(tidy(between_model) %>% mutate(model = 'between_model'),
          tidy(between_model_keywords1) %>% mutate(model = 'batch1'),
          tidy(between_model_keywords2) %>% mutate(model = 'batch2')) %>%
  filter(grepl('b_', term)) %>%
  dotwhisker::dwplot() 

between_posteriors_keywords1 / between_posteriors_keywords2
