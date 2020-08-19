# load libraries
library(brms)
library(tidyverse)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(broom)
library(dotwhisker)
library(skimr)
library(broom.mixed)

# load data
long_data <- read_csv(here::here('cleaned_numeric_data_long.csv'), col_types = cols(article_type = col_factor(),
                                                                                    Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor(),
                                                                                    keyword_batch_comp = col_factor())) %>%
             mutate(question = case_when(article_type == 'RR' ~ str_sub(question, 3),
                                         article_type == 'nonRR' ~ str_sub(question, 4))) %>%
            mutate(article_type = fct_relevel(article_type, c('nonRR', 'RR'))) %>%
            mutate(guessed_right_first = case_when(Order == 'RRFirst' & BelieveFirstRR == 1 ~ 1,
                                             Order == 'RRFirst' & (BelieveFirstRR == 2 | BelieveFirstRR == 3) ~ 0,
                                             Order == 'RRSecond' & (BelieveFirstRR == 1 | BelieveFirstRR == 2) ~ 0,
                                             Order == 'RRSecond' & BelieveFirstRR == 3 ~ 1),
                   guessed_right_second = case_when(Order == 'RRSecond' & BelieveSecondRR == 1 ~ 1,
                                                   Order == 'RRSecond' & (BelieveSecondRR == 2 | BelieveSecondRR == 3) ~ 0,
                                                   Order == 'RRFirst' & (BelieveSecondRR == 1 | BelieveSecondRR == 2) ~ 0,
                                                   Order == 'RRFirst' & BelieveSecondRR == 3 ~ 1)) %>%
            mutate(behavior_familiar = rowMeans(across(c(RRFamiliar,PreregFamiliar)), na.rm = T),
                   believe_improve = rowMeans(across(c(BelieveRigor,BelieveQuality)), na.rm = T))

wide_data <- read_csv(here::here('cleaned_numeric_data_wide.csv'), col_types = cols(Field = col_factor(),
                                                                                    Match = col_factor(),
                                                                                    Order = col_factor(),
                                                                                    keyword_batch_comp = col_factor())) %>%
                mutate(behavior_familiar = rowMeans(across(c(RRFamiliar,PreregFamiliar)), na.rm = T),
                        believe_improve = rowMeans(across(c(BelieveRigor,BelieveQuality)), na.rm = T)) %>%
                mutate(behavior_familiar_c = behavior_familiar - mean(behavior_familiar),
                       believe_improve_c = believe_improve - mean(believe_improve),
                       firstqualified_c = FirstQualified - mean(FirstQualified),
                       secondqualified_c = SecondQualified - mean(SecondQualified)) %>%
                mutate(guessed_RR_right = case_when(Order == 'RRFirst' & BelieveFirstRR == 1 ~ 1,
                                         Order == 'RRFirst' & (BelieveFirstRR == 2 | BelieveFirstRR == 3) ~ 0,
                                         Order == 'RRSecond' & BelieveSecondRR == 1 ~ 0,
                                         Order == 'RRSecond' & (BelieveSecondRR == 3 | BelieveSecondRR == 2) ~ 1),
                      guessed_nonRR_right = case_when(Order == 'RRSecond' & BelieveFirstRR == 3 ~ 1,
                                                Order == 'RRSecond' & (BelieveFirstRR == 2 | BelieveFirstRR == 1) ~ 0,
                                                Order == 'RRFirst' & (BelieveSecondRR == 1 | BelieveSecondRR == 2) ~ 0,
                                                Order == 'RRFirst' & (BelieveSecondRR == 3) ~ 1),
                      guessed_right = case_when(guessed_RR_right == 1 & guessed_nonRR_right == 1 ~ 'both',
                                              guessed_RR_right == 0 & guessed_nonRR_right == 0 ~ 'neither',
                                              TRUE ~ 'half'))
            
                       

# set up contrasts codes
contrasts(wide_data$Field) <- contr.sum(3)
contrasts(wide_data$keyword_batch_comp) <- contr.sum(2)
contrasts(wide_data$Order) <- contr.sum(2)
contrasts(wide_data$Match) <- contr.sum(2)
contrasts(long_data$article_type) <- contr.treatment(2)
contrasts(long_data$Field) <- contr.sum(3)
contrasts(long_data$keyword_batch_comp) <- contr.sum(2)
contrasts(long_data$Order) <- contr.sum(2)
contrasts(long_data$Match) <- contr.sum(2)

# check means & SDs of each question by article type
long_data %>%
  group_by(question, article_type) %>%
  summarize(mean = mean(response, na.rm = T), sd = sd(response, na.rm = T))


# bayesian models for single DV
priors <- c(set_prior("normal(0,2)", "b"),
            set_prior("normal(0, 2.5)", "sd"))




#### functions for within models ####

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



#### functions for between subjects models ####

# between-subjects models comparing pooled, batch 1, batch 2 + 3
between_pooled_model <- function(dv, set_priors) {
  between_model <- brm(response ~ Field + keyword_batch_comp + article_type + Match + article_type*Match +
                         (article_type|RR),
                       data = long_data %>% filter(grepl(as.character(dv), question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
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
                                   filter(grepl(as.character(dv), question)) %>% 
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
                                   filter(grepl(as.character(dv), question)) %>% 
                                   filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                                 prior = priors,
                                 family = 'gaussian',
                                 chains = 4)
}

# function to create graph of posterior samples for intercept
create_posteriors_btw <- function(results, variable){
  posteriors <- list(posterior_samples(results) %>%
                       select('b_article_type2') %>%
                       rename(setNames('b_article_type2', variable)))
  
  return(posteriors)
}

# Set up which model/prior/dv combinations to run for between models
between_models <- crossing(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                          set_priors = c(list(priors))) %>%
  mutate(between_pooled_model_results = pmap(list(dv, set_priors), between_pooled_model)) %>%
  mutate(posteriors = pmap(list(between_pooled_model_results, variable = dv), create_posteriors_btw))



#graphs for between models
# get all intercepts into wide format for graphing
intercepts <- c()

for (i in 1:nrow(between_models)) {
  intercepts <- bind_cols(intercepts, between_models$posteriors[[i]])
}

mcmc_areas(intercepts,
           prob=.95)

mcmc_intervals(intercepts, prob = .95)

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


#### demographics and sample descriptives ####
for_descriptive <- wide_data %>%
  select(-starts_with('diff'), -starts_with('Alt'), -c(StartDate, EndDate)) %>%
  select(Field, SPSubfield, Keyword, Gender, Education, ProfTitle, FirstFamiliar, FirstQualified,
         SecondFamiliar, SecondQualified, RRFamiliar, EverRR, PreregFamiliar, EverPrereg,
         BelieveRigor, BelieveQuality, BelieveFirstRR, BelieveSecondRR, RR, Order, Match)



## Overall descriptives
for_descriptive %>%
  select(-c(Field, SPSubfield, Keyword, Gender, Education, ProfTitle, BelieveFirstRR, BelieveSecondRR, RR, Order, Match, EverRR, EverPrereg), -starts_with('First'), -starts_with('Second')) %>%
  skim_to_wide() %>%
  mutate(missing = as.numeric(missing),
         complete = as.numeric(complete),
         n = as.numeric(n),
         mean = as.numeric(mean),
         sd = as.numeric(mean),
         p0 = as.numeric(p0),
         p25 = as.numeric(p25),
         p50 = as.numeric(p50),
         p75 = as.numeric(p75),
         p100 = as.numeric(p100)) %>%
  gt()

for_descriptive %>%
  select(RRFamiliar) %>%
  group_by(RRFamiliar) %>%
  tally() %>%
  mutate(RRFamiliar = as.factor(RRFamiliar),
         RRFamiliar = fct_recode(RRFamiliar, `Not at all familiar` = "1", `Slightly familiary` = "2", `Somewhat familiary` = "3",
                                     `Moderately familiary` = "4", `Substantially familiar` = "5"))

for_descriptive %>%
  select(PreregFamiliar) %>%
  group_by(PreregFamiliar) %>%
  tally() %>%
  mutate(PreregFamiliar = as.factor(PreregFamiliar),
         PreregFamiliar = fct_recode(PreregFamiliar, `Not at all familiar` = "1", `Slightly familiary` = "2", `Somewhat familiary` = "3",
                                 `Moderately familiary` = "4", `Substantially familiar` = "5"))

for_descriptive %>%
  select(BelieveRigor) %>%
  group_by(BelieveRigor) %>%
  tally() %>%
  mutate(BelieveRigor = as.factor(BelieveRigor),
         BelieveRigor = fct_recode(BelieveRigor,
                                   `Subtantially less` = "-4", `Much less` = "-3", `Moderately less` = "-2", `Slightly less` = "-1", `Not change` = "0",
                                   `Subtantially more` = "4", `Much more` = "3", `Moderately more` = "2", `Slightly more` = "1" ))

for_descriptive %>%
  select(BelieveQuality) %>%
  group_by(BelieveQuality) %>%
  tally() %>%
  mutate(BelieveQuality = as.factor(BelieveQuality),
         BelieveQuality = fct_recode(BelieveQuality,
                                   `Subtantially less` = "-4", `Much less` = "-3", `Moderately less` = "-2", `Slightly less` = "-1", `Not change` = "0",
                                   `Subtantially more` = "4", `Much more` = "3", `Moderately more` = "2", `Slightly more` = "1" ))


#### RR guesses

## RR guesses by order
for_descriptive %>%
  mutate(RR_isRR = case_when(Order == 'RRFirst' ~ BelieveFirstRR,
                                  Order == 'RRSecond' ~ BelieveSecondRR),
         Alt_isRR = case_when(Order == 'RRFirst' ~ BelieveSecondRR,
                                   Order == 'RRSecond' ~ BelieveFirstRR)) %>%
  summarize(RR_yesRR = sum(RR_isRR == 1), RR_notRR = sum(RR_isRR == 3), RR_unsure = sum(RR_isRR == 2),
            Alt_yesRR = sum(Alt_isRR == 1), Alt_notRR = sum(Alt_isRR == 3), Alt_unsure = sum(Alt_isRR == 2))

# RR guesses by article & order
for_descriptive %>%
  mutate(RR_isRR = case_when(Order == 'RRFirst' ~ BelieveFirstRR,
                             Order == 'RRSecond' ~ BelieveSecondRR),
         Alt_isRR = case_when(Order == 'RRFirst' ~ BelieveSecondRR,
                              Order == 'RRSecond' ~ BelieveFirstRR)) %>%
  group_by(RR) %>%
  summarize(RR_yesRR = sum(RR_isRR == 1), RR_notRR = sum(RR_isRR == 3), RR_unsure = sum(RR_isRR == 2),
            Alt_yesRR = sum(Alt_isRR == 1), Alt_notRR = sum(Alt_isRR == 3), Alt_unsure = sum(Alt_isRR == 2))

#### Article review qualifications

# by order
for_descriptive %>%
  mutate(RR_qualified = case_when(Order == 'RRFirst' ~ FirstQualified,
                                  Order == 'RRSecond' ~ SecondQualified),
         Alt_qualified = case_when(Order == 'RRFirst' ~ SecondQualified,
                                   Order == 'RRSecond' ~ FirstQualified)) %>%
  select(Order, RR_qualified, Alt_qualified) %>%
  summarize(RR_qualified_mean = mean(RR_qualified), RR_qualified_sd = sd(RR_qualified),
            Alt_qualified_mean = mean(Alt_qualified), Alt_qualified_sd = sd(Alt_qualified))

# by order & article
for_descriptive %>%
  mutate(RR_qualified = case_when(Order == 'RRFirst' ~ FirstQualified,
                                  Order == 'RRSecond' ~ SecondQualified),
         Alt_qualified = case_when(Order == 'RRFirst' ~ SecondQualified,
                                   Order == 'RRSecond' ~ FirstQualified)) %>%
  select(RR, RR_qualified, Alt_qualified) %>%
  group_by(RR) %>%
  summarize(RR_qualified_mean = mean(RR_qualified), RR_qualified_sd = sd(RR_qualified),
            Alt_qualified_mean = mean(Alt_qualified), Alt_qualified_sd = sd(Alt_qualified))

#### RR/Prereg behavior

for_descriptive %>%
  group_by(EverPrereg) %>%
  tally() %>%
  mutate(EverPrereg = as.factor(EverPrereg),
         EverPrereg = fct_recode(EverPrereg,
                                 No = '1', `Yes, once` = '2', `Yes, a few times` = '3', `Yes, many times` = '4'))
  
for_descriptive %>%
  group_by(EverRR) %>%
  tally() %>%
  mutate(EverRR = as.factor(EverRR),
         EverRR = fct_recode(EverRR,
                                 No = '2', Yes = '1'))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

for_descriptive %>%
  select(BelieveQuality, BelieveRigor, PreregFamiliar, RRFamiliar, FirstQualified, SecondQualified, Order) %>%
  mutate(RR_qualified = case_when(Order == 'RRFirst' ~ FirstQualified,
                                  Order == 'RRSecond' ~ SecondQualified),
         Alt_qualified = case_when(Order == 'RRFirst' ~ SecondQualified,
                                   Order == 'RRSecond' ~ FirstQualified)) %>%
  select(-c("Order", "FirstQualified", "SecondQualified")) %>%
  cor(use = 'everything') %>%
  corrplot(method="color", col=col(200),  
           type="lower", order="hclust", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation
           # hide correlation coefficient on the principal diagonal
           diag=FALSE 
  )
  

#### Exploratory Analyses ####

### between subjects models 

## add covariates for 'Qualified to review', behavior_familiar (item average score of Pre-reg adn RR familiarity items), believe_improve (item average score of BelieveRigor and BelieveQuality items)
between_pooled_model_covariates <- function(dv, set_priors) {
  between_model <- brm(response ~ Field + keyword_batch_comp + FirstQualified + behavior_familiar + believe_improve
                                      + article_type + Match + article_type*Match +
                         (article_type|RR),
                       data = long_data %>% filter(grepl(as.character(dv), question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                       prior = priors,
                       family = 'gaussian',
                       chains = 4)
  return(between_model)
}


between_models_covariates <- crossing(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                           set_priors = c(list(priors))) %>%
  mutate(between_pooled_covariates_results = pmap(list(dv, set_priors), between_pooled_model_covariates)) %>%
  mutate(posteriors = pmap(list(between_pooled_covariates_results, variable = dv), create_posteriors_btw))

intercepts <- c()

for (i in 1:nrow(between_models_covariates)) {
  intercepts <- bind_cols(intercepts, between_models_covariates$posteriors[[i]])
}

mcmc_areas(intercepts,
           prob=.95)

mcmc_intervals(intercepts, prob = .95)

## split data and look at results (with covariates) for those who guess right and didn't guess right
between_pooled_model_covariates_right <- function(dv, set_priors) {
  between_model <- brm(response ~ Field + keyword_batch_comp + FirstQualified + behavior_familiar + believe_improve
                       + article_type + Match + article_type*Match +
                         (article_type|RR),
                       data = long_data %>% filter(grepl(as.character(dv), question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')) %>% filter(guessed_right_first == 1),
                       prior = priors,
                       family = 'gaussian',
                       chains = 4)
  return(between_model)
}

between_pooled_model_covariates_wrong <- function(dv, set_priors) {
  between_model <- brm(response ~ Field + keyword_batch_comp + FirstQualified + behavior_familiar + believe_improve
                       + article_type + Match + article_type*Match +
                         (article_type|RR),
                       data = long_data %>% filter(grepl(as.character(dv), question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')) %>% filter(guessed_right_first == 0),
                       prior = priors,
                       family = 'gaussian',
                       chains = 4)
  return(between_model)
}

# run both models for all DVs
between_models_covariates_right <- crossing(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                                      set_priors = c(list(priors))) %>%
  mutate(between_pooled_covariates_right_results = pmap(list(dv, set_priors), between_pooled_model_covariates_right)) %>%
  mutate(posteriors = pmap(list(between_pooled_covariates_right_results, variable = dv), create_posteriors_btw))

between_models_covariates_wrong <- crossing(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                                      set_priors = c(list(priors))) %>%
  mutate(between_pooled_covariates_wrong_results = pmap(list(dv, set_priors), between_pooled_model_covariates_wrong)) %>%
  mutate(posteriors = pmap(list(between_pooled_covariates_wrong_results, variable = dv), create_posteriors_btw))

## function to keep only article_type coefficient for graphing
extract_coef <- function(coef, fitted_model){
  tidied_result <- tidy(fitted_model) %>%
    filter(grepl(as.character(coef), term) & effect == 'fixed')
  
  return(tidied_result)
}

## create graph showing intervals side by side
bind_rows(between_models_covariates_right %>%
  mutate(tidied = pmap(list(between_pooled_covariates_right_results), extract_coef)) %>%
  select(dv, tidied) %>%
  mutate(model = 'Yes') %>%
  unnest_wider(tidied) %>%
  rbind() %>%
  rename(coef = term,
         term = dv),
between_models_covariates_wrong %>%
  mutate(tidied = pmap(list(between_pooled_covariates_wrong_results), extract_coef)) %>%
  select(dv, tidied) %>%
  mutate(model = 'No') %>%
  unnest_wider(tidied) %>%
  rbind() %>%
  rename(coef = term,
         term = dv)) %>%
  dotwhisker::dwplot()

## descriptives for guessing

# means by article type and guessed right
View(long_data %>% 
       filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')) %>% 
       group_by(question, article_type, guessed_right_first) %>% 
       summarize(mean = mean(response, na.rm = T), n = n()))

# what is correlated with guessing right?
long_data %>% 
  filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')) %>% 
  distinct(participant_id, .keep_all = T) %>%
  mutate(Field = as.numeric(Field),
         EverPrereg = case_when(is.na(EverPrereg) ~ 0,
                                !is.na(EverPrereg) ~ EverPrereg),
         EverRR = case_when(is.na(EverRR) ~ 0,
                            !is.na(EverRR) ~ EverPrereg)) %>%
  select(guessed_right_first, Field, RRFamiliar, EverRR, PreregFamiliar, EverPrereg, BelieveRigor, BelieveQuality) %>%
  cor(method = 'spearman') %>%
  as_tibble() %>%
  slice(1L)

# what is correct guess rate by article?
bind_rows(long_data %>% 
  filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR') & article_type == 'RR') %>% 
  distinct(participant_id, .keep_all = T) %>%
  group_by(RR, article_type, guessed_right_first) %>%
  tally() %>%
  group_by(RR, article_type) %>%
  mutate(perc_guess = round(100 *n/sum(n),2)) %>%
  select(-n) %>%
  pivot_wider(names_from = 'guessed_right_first', values_from = 'perc_guess') %>%
  rename(guessed_right = `1`,
         guessed_not_right = `0`) %>%
  mutate(guessed_right = case_when(is.na(guessed_right) ~ 0,
                                   TRUE ~ guessed_right),
         guessed_not_right = case_when(is.na(guessed_not_right) ~ 0,
                                       TRUE ~ guessed_not_right)) %>%
  arrange(article_type, guessed_right),
  long_data %>% 
              filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR') & article_type == 'nonRR') %>% 
              distinct(participant_id, .keep_all = T) %>%
              group_by(RR, Match, guessed_right_first) %>%
              tally() %>%
              group_by(RR, Match) %>%
              mutate(perc_guess = round(100 *n/sum(n),2)) %>%
              select(-n) %>%
              pivot_wider(names_from = 'guessed_right_first', values_from = 'perc_guess') %>%
              rename(guessed_right = `1`,
                     guessed_not_right = `0`,
                     article_type = Match) %>%
              mutate(guessed_right = case_when(is.na(guessed_right) ~ 0,
                                               TRUE ~ guessed_right),
                     guessed_not_right = case_when(is.na(guessed_not_right) ~ 0,
                                                   TRUE ~ guessed_not_right)) %>%
              arrange(article_type, guessed_right)) %>%
  write_csv('guessed_first_right_article_perc.csv')

### within subjects models 
within_diff_pooled_covariate_model <- function(dv, set_priors) {
  within_model_diffs <- brm(as.formula(paste(dv, "~ Field + keyword_batch_comp + firstqualified_c + secondqualified_c + behavior_familiar_c + believe_improve_c + 
                                                      Order + Match + Order*Match +
                                                    (1|RR)")),
                            data = wide_data,
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_diffs)
}


within_diff_pooled_covariate_guessed_model <- function(dv, set_priors, guessed) {
  within_model_diffs <- brm(as.formula(paste(dv, "~ Field + keyword_batch_comp + firstqualified_c + secondqualified_c + behavior_familiar_c + believe_improve_c + 
                                                      Order + Match + Order*Match +
                                                    (1|RR)")),
                            data = wide_data %>% filter(guessed_right == as.character(guessed)),
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_diffs)
}


within_covariate_models <- crossing(dv = names(wide_data[,68:86]),
                          set_priors = c(list(priors))) %>%
  mutate(within_pooled_model_covariate_results = pmap(list(dv, set_priors), within_diff_pooled_covariate_model)) %>%
  mutate(posteriors = pmap(list(within_pooled_model_covariate_results, variable = dv), create_posteriors))


within_diff_covariate_guessed_models <- crossing(dv = names(wide_data[,68:86]),
                                                 set_priors = c(list(priors)),
                                                 guessed = c('both', 'half', 'neither')) %>%
  mutate(within_pooled_model_covariate_guessed_results = pmap(list(dv, set_priors,guessed), within_diff_pooled_covariate_guessed_model)) %>%
  mutate(posteriors = pmap(list(within_pooled_model_covariate_guessed_results, variable = dv), create_posteriors))


# get all intercepts into wide format for graphing
intercepts <- c()

for (i in 1:nrow(within_covariate_models)) {
  intercepts <- bind_cols(intercepts, within_covariate_models$posteriors[[i]])
}

mcmc_areas(intercepts,
           prob=.95)

mcmc_intervals(intercepts, prob = .95)

# visualization for guessing models

within_diff_covariate_guessed_models %>%
  mutate(tidied = pmap(list('(Intercept)', within_pooled_model_covariate_guessed_results), extract_coef)) %>%
  select(dv, guessed, tidied) %>%
  unnest_wider(tidied) %>%
  rbind() %>%
  rename(coef = term,
         term = dv,
         model = guessed) %>%
  dotwhisker::dwplot()  


# descriptives for wide_data 
wide_data %>% 
  group_by(Order, guessed_nonRR_right, guessed_RR_right) %>% 
  tally()

# descriptives for wide_data 
wide_data %>% 
  group_by(Order, guessed_right) %>% 
  tally()

