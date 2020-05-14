# load libraries
library(brms)
library(tidyverse)
library(ggplot2)
library(bayesplot)
library(tidybayes)

# load data
long_data <- read_csv(here::here('cleaned_numeric_data_long.csv'), col_types = cols(article_type = col_factor(),
                                                                                    Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor()))
wide_data <- read_csv(here::here('cleaned_numeric_data_wide.csv'), col_types = cols(Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor()))

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


### Deciding between 2 possible within subjects model specifications:

# no different score model
within_model <- brm(response ~ Field + article_type + Order + Match +
                      article_type*Order + article_type*Match + Order*Match +
                      article_type*Order*Match + 
                      (1|participant_id) + (article_type|RR),
                    data = long_data %>% filter(grepl('QuestionQuality', question)),
                    prior = priors,
                    family = 'gaussian',
                    chains = 4,
                    seed = 1)

summary(within_model)
plot(within_model)
pp_check(within_model)
pp_check(within_model, type = "stat", stat = 'median', nsamples = 100)
WAIC(within_model)
loo(within_model)

# difference score model
within_model_diffs <- brm(diff_question_quality ~ Field + Order + Match +
                            Order*Match +
                            (1|RR),
                          data = wide_data,
                          prior = priors, 
                          family = 'gaussian',
                          chains = 4,
                          seed = 2)

summary(within_model_diffs)
plot(within_model_diffs)
pp_check(within_model_diffs)
pp_check(within_model_diffs, type = "stat", stat = 'median', nsamples = 100)
WAIC(within_model_diffs)
loo(within_model_diffs)

### within-subjects models

# difference score model, pooled across batches
within_model_diffs <- brm(diff_question_quality ~ Field + keyword_batch_comp + Order + Match +
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


# differenc score model for batch 1
within_model_diffs_keywords1 <- brm(diff_question_quality ~ Field + Order + Match +
                            Order*Match +
                            (1|RR),
                          data = wide_data %>% filter(keyword_batch_comp == 1),
                          prior = priors, 
                          family = 'gaussian',
                          chains = 4)

summary(within_model_diffs_keywords1)
pp_check(within_model_diffs_keywords1)
WAIC(within_model_diffs_keywords1)
loo(within_model_diffs_keywords1)

# differenc score model for batched 2+3
within_model_diffs_keywords2 <- brm(diff_question_quality ~ Field + Order + Match +
                                      Order*Match +
                                      (1|RR),
                                    data = wide_data %>% filter(keyword_batch_comp == 2),
                                    prior = priors, 
                                    family = 'gaussian',
                                    chains = 4)

summary(within_model_diffs_keywords2)
pp_check(within_model_diffs_keywords2)
WAIC(within_model_diffs_keywords2)
loo(within_model_diffs_keywords2)


posteriors_keywords2 <- suppressMessages( 
  mcmc_areas(posterior_samples(within_model_diffs_keywords2),
             regex_pars = "b_",
             prob=.9) +
             xlim(-2, 2) +
             labs(title = 'Batch 2')
)

posteriors_keywords1 <- suppressMessages( 
  mcmc_areas(posterior_samples(within_model_diffs_keywords1),
             regex_pars = "b_",
             prob=.9) +
             xlim(-2, 2) +
             labs(title = 'Batch 1')
)

posteriors_pooled <- suppressMessages( 
  mcmc_areas(posterior_samples(within_model_diffs),
             regex_pars = "b_",
             prob=.9) +
             xlim(-2, 2) +
             labs(title = 'Pooled Batched')
)

bind_rows(tidy(within_model_diffs) %>% mutate(model = 'within_diff_pooled'),
          tidy(within_model_diffs_keywords1) %>% mutate(model = 'batch1'),
          tidy(within_model_diffs_keywords2) %>% mutate(model = 'batch2')) %>%
  filter(grepl('b_', term)) %>%
  dotwhisker::dwplot() 

posteriors_keywords1 / posteriors_keywords2

# no difference score model, pooled across batches
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



