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
library(lavaan)
library(psych)
library(patchwork)
library(ggdist)
library(gt)
library(corrplot)

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
                
                # add categorical predictors for improve and familiar variables
                mutate(improve_6L = as.factor(case_when(believe_improve < 0 ~ 'negative',
                                              believe_improve == 0 ~ 'neutral',
                                              believe_improve > 0 & believe_improve <= 1 ~ 'slightly_more',
                                              believe_improve > 1 & believe_improve <= 2 ~ 'moderately_more',
                                              believe_improve > 2 & believe_improve <= 3 ~ 'much_more',
                                              believe_improve > 3 & believe_improve <= 4 ~ 'substantially_more')),
                       familiar_5L = as.factor(floor(behavior_familiar)),
                       improve_6L = fct_relevel(improve_6L, c('negative', 'neutral', 'slightly_more', 'moderately_more', 'much_more', 'substantially_more')),
                       familiar_5L = fct_relevel(familiar_5L, c('1', '2', '3', '4', '5'))) %>%
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
                                              TRUE ~ 'half'),
                      guessed_right = as.factor(guessed_right),
                      guessed_right = fct_relevel(guessed_right, c('neither', 'half', 'both')))
            
                       

# set up contrasts codes
contrasts(wide_data$Field) <- contr.sum(3)
contrasts(wide_data$keyword_batch_comp) <- contr.sum(2)
contrasts(wide_data$Order) <- contr.sum(2)
contrasts(wide_data$Match) <- contr.sum(2)
contrasts(wide_data$improve_6L) <- contr.treatment(6)
contrasts(wide_data$familiar_5L) <- contr.treatment(5)
contrasts(wide_data$guessed_right) <- contr.treatment(3)
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
priors <- c(set_prior("normal(0,2", "Intercept"),
            set_prior("normal(0,2)", "b"),
            set_prior("normal(0, 2.5)", "sd"),
            set_prior("normal(0, 2)", "sigma"))

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


# create more flexible version of posteriors function
create_posteriors_term <- function(results, variable, term){
  posteriors <- list(posterior_samples(results) %>%
                       select(term) %>%
                       rename(setNames(term, variable)))
  
  return(posteriors)
}

# Set up which model/prior/dv combinations to run
within_models <- crossing(dv = names(wide_data[,68:86]),
                    set_priors = c(list(priors))) %>%
  mutate(within_pooled_model_results = pmap(list(dv, set_priors), within_diff_pooled_model)) %>%
  mutate(posteriors = pmap(list(within_pooled_model_results, variable = dv, term = 'b_Intercept'), create_posteriors_term))


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
between_pooled_model <- function(dv, seed_num) {
  between_model <- brm(response ~ Field + keyword_batch_comp + article_type + Match + article_type*Match +
                         (article_type|RR),
                       data = long_data %>% filter(grepl(as.character(dv), question)) %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                       prior = priors,
                       family = 'gaussian',
                       iter = 6000,
                       seed = seed_num,
                       control = list(adapt_delta = 0.95),
                       chains = 4)
  return(between_model)
  }


# between-subjects model for batch 1
between_keywords1_model <- function(dv, seed_num) {
  between_keywords1 <- brm(response ~ Field + article_type + Match + article_type*Match +
                                   (article_type|RR),
                                 data = long_data %>% 
                                   filter(keyword_batch_comp == 1) %>%
                                   filter(grepl(as.character(dv), question)) %>% 
                                   filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                                 prior = priors,
                                 family = 'gaussian',
                                 iter = 6000,
                                 seed = seed_num,
                                 control = list(adapt_delta = 0.95),
                                 chains = 4)
  return(between_keywords1)
}

between_models_keywords1 <- cbind(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                                  seed_num = c(201:219)) %>%
                              as_tibble() %>%
                              mutate(seed_num = as.numeric(seed_num)) %>%
                              mutate(between_pooled_model_results = pmap(list(dv, seed_num), between_keywords1_model))

# differenc score model for batched 2+3
between_keywords2_model <- function(dv, seed_num) {
  between_keywords2 <- brm(response ~ Field + article_type + Match + article_type*Match +
                                   (article_type|RR),
                                 data = long_data %>% 
                                   filter(keyword_batch_comp == 2) %>%
                                   filter(grepl(as.character(dv), question)) %>% 
                                   filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                                 prior = priors,
                                 family = 'gaussian',
                                 iter = 6000,
                                 seed = seed_num,
                                 control = list(adapt_delta = 0.99),
                                 chains = 4)
}

between_models_keywords2 <- cbind(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                                  seed_num = c(301:319)) %>%
                              as_tibble() %>%
                              mutate(seed_num = as.numeric(seed_num)) %>%
                              as_tibble() %>%
                              mutate(seed_num = as.numeric(seed_num)) %>%
                              mutate(between_pooled_model_results = pmap(list(dv, seed_num), between_keywords2_model))

# Set up which run all btw subj individual dv models
between_models <- cbind(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                          seed_num = c(101:119)) %>%
  as_tibble() %>%
  mutate(seed_num = as.numeric(seed_num)) %>%
  mutate(between_pooled_model_results = pmap(list(dv, seed_num), between_pooled_model))

# create function to get posteriors for individual models
btw_posteriors <- function(results) {
  results %>%
    spread_draws(b_article_type2) %>%
    mean_qi(article_effect = b_article_type2, .width = c(.95, .80))
}

#get and save posteriors for each dv
btw_posteriors_individual_dvs <- between_models %>%
  mutate(posteriors = map(between_pooled_model_results, btw_posteriors)) %>%
  unnest(posteriors) %>%
  select(-c(seed_num, between_pooled_model_results)) %>%
  mutate(model = 'Individual DVs')

btw_posteriors_keyowrd1 <- between_models_keywords1  %>%
  mutate(posteriors = map(between_pooled_model_results, btw_posteriors)) %>%
  unnest(posteriors) %>%
  select(-c(seed_num, between_pooled_model_results)) %>%
  mutate(model = 'Initial Keywords')

btw_posteriors_keyword2 <- between_models_keywords2  %>%
  mutate(posteriors = map(between_pooled_model_results, btw_posteriors)) %>%
  unnest(posteriors) %>%
  select(-c(seed_num, between_pooled_model_results)) %>%
  mutate(model = 'Reduced Keywords')

### create graph comparing individual DV btw subjs keyword models
compare_btw_keyword_models <- rbind(btw_posteriors_individual_dvs,
                                    btw_posteriors_keyowrd1,
                                    btw_posteriors_keyword2) %>%
                                as.data.frame() %>%
                                mutate(dv = fct_rev(dv)) %>%
                                mutate(model = case_when(model == 'Individual DVs' ~ 'Pooled Keywords',
                                                         TRUE ~ model)) %>%
                                rename(Model = 'model') %>%
                                ggplot(aes(y = dv, x = article_effect, xmin = .lower, xmax = .upper, color = Model)) +
                                geom_pointinterval(position=position_dodge(width=0.85)) +
                                geom_vline(xintercept = 0) +
                                scale_x_continuous(breaks=seq(-1, 1.5, .5),
                                                   limits = c(-1, 1.75),
                                                   name = 'Difference between RR and non-RR articles') +
                                theme_minimal() +
                                theme(axis.title.y = element_blank(),
                                      axis.title.x = element_text(size = 14),
                                      axis.text = element_text(size = 14),
                                      panel.grid.minor.y = element_blank(),
                                      strip.text = element_text(size=14))

compare_btw_keyword_models

## create graph comparing posteriors across models for each DV
compare_btw_DV_models_graph <- rbind(btw_posteriors_individual_dvs,
                            posteriors_btw_mlm %>%
                              mutate(model = 'All DVs') %>%
                              rename(dv = "DV") %>%
                              select(-article_type) %>%
                              ungroup()) %>%
  as.data.frame() %>%
  mutate(dv = fct_rev(dv)) %>%
  rename(Model = 'model') %>%
  ggplot(aes(y = dv, x = article_effect, xmin = .lower, xmax = .upper, color = Model)) +
  geom_pointinterval(position=position_dodge(width=0.75)) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks=seq(-.5, 1.5, .5),
                     limits = c(-.5, 1.5),
                     name = 'Difference between RR and non-RR articles') +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size=14))

compare_btw_DV_models_graph



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

#### means & sds for Familiarity and Rigor questions
for_descriptive %>%
  select(RRFamiliar, PreregFamiliar, BelieveRigor, BelieveQuality) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))

#### cor for rigor and quality beliefs
cor(for_descriptive$BelieveRigor, for_descriptive$BelieveQuality)
cor(for_descriptive$RRFamiliar, for_descriptive$PreregFamiliar)

### every RR or Pre-reg tallys
for_descriptive %>%
  group_by(EverRR) %>%
  tally()

for_descriptive %>%
  group_by(EverPrereg) %>%
  tally()

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

### within model with DVs as ML component
mlm_dvs_data <- wide_data %>%
  select(starts_with('diff'), participant_id, RR, Field, keyword_batch_comp, Order, Match, improve_6L, familiar_5L, guessed_right) %>%
  pivot_longer(cols = starts_with('diff'), names_to = 'questions', values_to = 'response') %>%
  mutate(questions = as.factor(questions))

within_alldvs <-  brm(response ~ Field + keyword_batch_comp + 
                            Order + Match + Order*Match +
                            (1|RR) + (1|participant_id) + (1|questions),
                          data = mlm_dvs_data,
                          prior = priors,
                          family = 'gaussian',
                          chains = 4,
                          iter = 4000,
                          seed = 25,
                          control = list(adapt_delta = 0.95, max_treedepth = 15))


# main figure for RR vs. non-RR difference across DVs with partial pooling
with_alldvs_graph_nums <- within_alldvs %>%
  spread_draws(b_Intercept, r_questions[dv,]) %>%
  mutate(dv_estimates = b_Intercept + r_questions) %>%
  left_join(within_alldvs %>%
         spread_draws(b_Intercept, r_questions[dv,]) %>% 
         mean_qi(dv_mean = b_Intercept + r_questions) %>% 
         select(dv, dv_mean), by = 'dv') %>%
  ungroup()

## setup function to create base graph for each section
main_graph_creation <- function(data) {
  data %>%  
    mutate(dv = as.factor(dv),
           dv = fct_reorder(dv, dv_mean)) %>%
    ggplot(aes(y = dv, x = dv_estimates, fill = stat(x <= 0))) +
    stat_halfeye(point_interval = mean_qi, .width = c(.95, .8)) +
    scale_x_continuous(breaks=seq(-.5, 1.5, .5),
                       limits = c(-.75, 1.8),
                       name = 'Difference between RR and non-RR articles') +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 16),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

## graph for intro questions
intro_qs <- with_alldvs_graph_nums %>%
  filter(grepl('question', dv) |
         grepl('method', dv) |
         dv == 'diff_aligned' | 
         dv == 'diff_will_learn' | 
         dv == 'diff_intro_importance') %>%
  mutate(dv = case_when(dv == 'diff_aligned' ~ 'Question & methods aligned',
                      dv == 'diff_will_learn' ~ 'Amount will learn',
                      dv == 'diff_intro_importance' ~ 'Important research',
                      dv == 'diff_method_rigor' ~ 'Methods rigor',
                      dv == 'diff_method_quality' ~ 'Quality of methods',
                      dv == 'diff_question_quality' ~ 'Quality of question',
                      dv == 'diff_method_creative' ~ 'Creativity of methods',
                      dv == 'diff_question_novel' ~ 'Novelty of question')) %>%
  main_graph_creation() +
  scale_fill_manual(values = c("#fbb4ae", "gray80")) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle('Evaluation before knowing study outcomes') +
  theme(plot.title = element_text(face = 'bold', size = 18))

## graph for results questions
results_qs <- with_alldvs_graph_nums %>%
                filter(grepl('result', dv) |
                         dv == "diff_analysis_rigor" |
                         dv == "diff_overall_import" |
                         dv == "diff_did_learn" |
                         dv == "diff_discussion_quality" |
                         dv == "diff_justificed") %>%
  mutate(dv = case_when(dv == 'diff_analysis_rigor' ~ 'Analysis rigor',
                        dv == 'diff_overall_import' ~ 'Important findings',
                        dv == 'diff_justificed' ~ 'Conclusions justified',
                        dv == 'diff_result_quality' ~ 'Quality of results',
                        dv == 'diff_discussion_quality' ~ 'Qualtiy of discussion',
                        dv == 'diff_did_learn' ~ 'Amount learned',
                        dv == 'diff_result_innovative' ~ 'Innovative results')) %>%
  main_graph_creation() +
  scale_fill_manual(values = c("#b3cde3", "gray80")) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle('Evaluation after knowing study outcomes') +
  theme(plot.title = element_text(face = 'bold', size = 18))

## graph fo abstract questions
abstract_qs <- with_alldvs_graph_nums %>%
                  filter(dv == 'diff_inspire' |
                           dv == 'diff_abstract_aligned' |
                           dv == 'diff_field_importance' |
                           dv == "diff_overall_quality") %>%
  mutate(dv = case_when(dv == 'diff_overall_quality' ~ 'Overall quality of paper',
                        dv == 'diff_field_importance' ~ 'Important discoveries',
                        dv == 'diff_abstract_aligned' ~ 'Abstract & findings aligned',
                        dv == 'diff_inspire' ~ 'Inspire new research')) %>%
                main_graph_creation() +
  scale_fill_manual(values = c("#ccebc5", "gray80"))+
  ggtitle('Evaluation after finishing the paper') +
  theme(plot.title = element_text(face = 'bold', size = 18),
        axis.title.x = element_text(hjust = .15, vjust = 0))

## combine graphs
combined_plot <- intro_qs / results_qs / abstract_qs + plot_layout(heights = c(8, 7, 4))
combined_plot


### compare model estimates from hierarchical and non-hierarchical DV models
rbind(within_models %>% 
  mutate(draws = map(within_pooled_model_results, spread_draws, b_Intercept),
         posterior_info = map(draws, median_qi, .width = c(.95, .90))) %>%
  select(dv, posterior_info) %>%
  unnest(posterior_info) %>%
  rename(question_mean = b_Intercept) %>%
  mutate(model = 'no dv pooling'),
  within_alldvs %>%
    spread_draws(b_Intercept, r_questions[dv,]) %>%
    median_qi(question_mean = b_Intercept + r_questions, .width = c(.95, .90)) %>%
    mutate(model = 'partial dv pooling')) %>%
  ggplot(aes(y = dv, x = question_mean, xmin = .lower, xmax = .upper, color = model)) +
  geom_pointinterval(position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0) +
  theme_classic()
  

### Btw subjects model with DV as ML
between_model_mlm <- brm(response ~ Field + keyword_batch_comp + article_type + Match + article_type*Match +
                       (1|RR) + (1|question) + (1|participant_id),
                     data = long_data %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                     prior = priors,
                     family = 'gaussian',
                     seed = 53,
                     chains = 4,
                     iter = 8000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))

between_model_mlm_slopes <- brm(response ~ Field + keyword_batch_comp + article_type + Match + article_type*Match +
                           (article_type|RR) + (article_type|question) + (1|participant_id),
                         data = long_data %>% filter((Order == 'RRFirst' & article_type == 'RR') | (Order == 'RRSecond' & article_type == 'nonRR')),
                         prior = priors,
                         family = 'gaussian',
                         seed = 54,
                         chains = 4,
                         iter = 8000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15))

# get posteriors by DV
posteriors_btw_mlm <- between_model_mlm_slopes %>%
                          spread_draws(b_article_type2, r_question[DV, article_type]) %>%
                          filter(article_type == 'article_type2') %>%
                          mean_qi(article_effect = b_article_type2 + r_question, .width = c(.95, .80))

# get overall, min and max across DVs
posteriors_btw_mlm %>%
  ungroup() %>%
  filter(.width == 0.95) %>%
  summarize(mean = round(mean(article_effect),2), min = round(min(article_effect),2), max = round(max(article_effect),2))

# partial pooling across all DVs and improve

within_alldvs_improve <-  brm(response ~ Field + keyword_batch_comp + 
                        Order + Match + Order*Match +
                        (1|RR) + (1|participant_id) + (1|questions) + (1|improve_6L),
                      data = mlm_dvs_data,
                      prior = priors,
                      family = 'gaussian',
                      chains = 4,
                      iter = 3000,
                      seed = 28,
                      control = list(adapt_delta = 0.99, max_treedepth = 15))

summary(within_alldvs_improve)
pp_check(within_alldvs_improve)
WAIC(within_alldvs_improve)
loo(within_alldvs_improve)

# parameters by improvement level
within_alldvs_improve %>%
  spread_draws(b_Intercept, r_improve_6L[improve_level,]) %>%
  mean_qi(cond_mean = b_Intercept + r_improve_6L, .width = c(.95)) %>%
  mutate_if(is.numeric, round, 2)


### graph collapsed across DVs
within_alldvs_improve %>%
  spread_draws(b_Intercept, r_improve_6L[improve_level,]) %>%
  median_qi(cond_mean = b_Intercept + r_improve_6L, .width = c(.95, .90)) %>%
  ungroup() %>%
  mutate(improve_level = as.factor(improve_level),
         improve_level = fct_relevel(improve_level, c('negative', 'neutral', 'slightly_more', 'moderately_more', 'much_more', 'substantially_more'))) %>%
  ggplot(aes(y = improve_level, x = cond_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  theme_classic()

# partial pooling across all DVs and familiar

within_alldvs_familiar <-  brm(response ~ Field + keyword_batch_comp + 
                                Order + Match + Order*Match +
                                (1|RR) + (1|participant_id) + (1|questions) + (1|familiar_5L),
                              data = mlm_dvs_data,
                              prior = priors,
                              family = 'gaussian',
                              chains = 4,
                              iter = 4000,
                              seed = 27,
                              control = list(adapt_delta = .99, max_treedepth = 15))

# get table of estimates for familiarity levels across all DVs
within_alldvs_familiar %>%
  spread_draws(b_Intercept,r_familiar_5L[familiar_level,]) %>%
  mean_qi(fam_level_mean = b_Intercept + r_familiar_5L, .width = c(.95)) %>%
  mutate_if(is.numeric, round, 2)

### graph collapsed across DVs
within_alldvs_familiar %>%
  spread_draws(b_Intercept, r_familiar_5L[familiar_level,]) %>%
  median_qi(cond_mean = b_Intercept + r_familiar_5L, .width = c(.95, .90)) %>%
  ungroup() %>%
  mutate(familiar_level = as.factor(familiar_level),
         familiar_level = fct_relevel(familiar_level, c('1', '2', '3', '4', '5'))) %>%
  ggplot(aes(y = familiar_level, x = cond_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  theme_classic()

### model for guessing (random slope for guessing by question)
within_diff_pooled_guessed_model_slopes <- brm(response ~ Field + keyword_batch_comp + guessed_right +
                                          Order + Match + Order*Match +
                                          (1|RR) + (1|participant_id) + (guessed_right|questions),
                                        data = mlm_dvs_data,
                                        prior = priors,
                                        family = 'gaussian',
                                        chains = 4,
                                        iter = 4000,
                                        seed = 31,
                                        control = list(adapt_delta = .99, max_treedepth = 15))

# get posteriors for guessing across DVs
posteriors_by_guessing <- within_diff_pooled_guessed_model_slopes %>% 
  spread_draws(b_Intercept, b_guessed_right2, b_guessed_right3, r_questions[DV,guessed]) %>%
  pivot_longer(cols = b_Intercept:b_guessed_right3,
               names_to = "fixed_guessed",
               values_to = 'fixed_value') %>%
  mutate(fixed_guessed = case_when(fixed_guessed == 'b_Intercept' ~ 'Intercept',
                                   fixed_guessed == 'b_guessed_right2' ~ 'guessed_right2',
                                   fixed_guessed == 'b_guessed_right3' ~ 'guessed_right3')) %>%
  filter(guessed == fixed_guessed) %>%
  mean_qi(mean = fixed_value + r_questions, .width = c(.95, .80)) %>%
  mutate(guessed = case_when(guessed == 'Intercept' ~ 'None',
                             guessed == 'guessed_right2' ~ 'Half',
                             guessed == 'guessed_right3' ~ 'Both'))

# calculate number in the same direction as main model [where all are positive]
posteriors_by_guessing %>%
  select(-c(.point, .interval)) %>%
  filter(.width == 0.95) %>%
  group_by(guessed) %>%
  summarize(num_same_direction = sum(mean > 0),
            mean_effect = mean(mean),
            min = min(mean),
            max = max(mean))
  

#### Supplemental Analyses

### Correlations between DVs
png(height = 600, width = 600, file = 'dv_corrs.png')

wide_data %>%
  select(starts_with('diff')) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(type = 'upper', order="hclust", method = 'circle')
dev.off() 

### Model fit and summary results from within-subjects ML pooled DV model

# table of results
within_alldvs_model_results <- summary(within_alldvs) 

within_alldvs_model_table <- rbind(within_alldvs_model_results$fixed %>%
        as.data.frame() %>% 
        rownames_to_column(var = "Predictor"), 
      within_alldvs_model_results$random$participant_id %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "Predictor") %>% 
        mutate(Predictor = 'sd(Intercept): Participant_id'),
      within_alldvs_model_results$random$questions %>%
        as.data.frame() %>% 
        rownames_to_column(var = "Predictor") %>% 
        mutate(Predictor = 'sd(Intercept): Question'),
      within_alldvs_model_results$random$RR %>%
        as.data.frame() %>% 
        rownames_to_column(var = "Predictor") %>% 
        mutate(Predictor = 'sd(Intercept): RR'),
      within_alldvs_model_results$spec_pars %>%
        as.data.frame() %>%
        rownames_to_column(var = "Predictor")) %>%
  mutate(effect_type = case_when(grepl('sd', Predictor) | grepl('sigma', Predictor) ~ 'Random Effects',
                                 TRUE ~ 'Fixed Effects')) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(-c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  rename(SE = 'Est.Error',
         lower_ci = "l-95% CI",
         upper_ci = "u-95% CI") %>%
  gt(groupname_col = 'effect_type') %>%
  tab_stubhead("label") %>% 
  cols_merge(columns = vars(lower_ci,upper_ci),
             hide_columns = vars(upper_ci),
             pattern = "[{1}, {2}]") %>%
  cols_label(lower_ci = '95% CrI',
             Predictor = '') %>%
  cols_align(align = 'center',
             columns = 2:4) %>%
  cols_align(align = 'left',
             columns = 1) %>%
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = list(
      cells_row_groups(),
      cells_column_labels(everything())
    )
  ) %>% 
  tab_options(
    row_group.border.top.width = px(3),
    row_group.border.top.color = "black",
    row_group.border.bottom.color = "black",
    table_body.hlines.color = "white",
    table.border.top.color = "white",
    table.border.top.width = px(3),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )

gtsave(within_alldvs_model_table, 'within_alldvs_model_table.rtf')  

# table of point estimates & CrI by DV
within_alldvs_dv_posteriors <- within_alldvs %>%
  spread_draws(b_Intercept, r_questions[dv,]) %>% 
  mean_qi(dv_mean = b_Intercept + r_questions) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(dv, dv_mean, .lower, .upper) %>%
  gt() %>%
  cols_merge(columns = vars(.lower, .upper),
             hide_columns = vars(.upper),
             pattern = "[{1}, {2}]") %>%
  cols_label(dv = 'DV',
             dv_mean = 'Posterior Mean',
             .lower = '95% CrI') %>%
  cols_align(align = 'center',
             columns = 2:3) %>%
  cols_align(align = 'left',
             columns = 1) %>%
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = list(
      cells_column_labels(everything())
    )
  ) %>% 
  tab_options(
    table_body.hlines.color = "white",
    table.border.top.color = "white",
    table.border.top.width = px(3),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )

gtsave(within_alldvs_dv_posteriors, 'within_alldvs_dv_posteriors.rtf')


#### Model fit and summary statistics for exploratory model with familiarity

within_familiar_model_results <- summary(within_alldvs_familiar) 

within_familiar_model_table <- rbind(within_familiar_model_results$fixed %>%
                                     as.data.frame() %>% 
                                     rownames_to_column(var = "Predictor"), 
                                     within_familiar_model_results$random$familiar_5L %>% 
                                       as.data.frame() %>% 
                                       rownames_to_column(var = "Predictor") %>% 
                                       mutate(Predictor = 'sd(Intercept): Familiarity'),
                                     within_familiar_model_results$random$participant_id %>% 
                                     as.data.frame() %>% 
                                     rownames_to_column(var = "Predictor") %>% 
                                     mutate(Predictor = 'sd(Intercept): Participant_id'),
                                     within_familiar_model_results$random$questions %>%
                                     as.data.frame() %>% 
                                     rownames_to_column(var = "Predictor") %>% 
                                     mutate(Predictor = 'sd(Intercept): Question'),
                                     within_familiar_model_results$random$RR %>%
                                     as.data.frame() %>% 
                                     rownames_to_column(var = "Predictor") %>% 
                                     mutate(Predictor = 'sd(Intercept): RR'),
                                     within_familiar_model_results$spec_pars %>%
                                     as.data.frame() %>%
                                     rownames_to_column(var = "Predictor")) %>%
  mutate(effect_type = case_when(grepl('sd', Predictor) | grepl('sigma', Predictor) ~ 'Random Effects',
                                 TRUE ~ 'Fixed Effects')) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(-c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  rename(SE = 'Est.Error',
         lower_ci = "l-95% CI",
         upper_ci = "u-95% CI") %>%
  gt(groupname_col = 'effect_type') %>%
  tab_stubhead("label") %>% 
  cols_merge(columns = vars(lower_ci,upper_ci),
             hide_columns = vars(upper_ci),
             pattern = "[{1}, {2}]") %>%
  cols_label(lower_ci = '95% CrI',
             Predictor = '') %>%
  cols_align(align = 'center',
             columns = 2:4) %>%
  cols_align(align = 'left',
             columns = 1) %>%
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = list(
      cells_row_groups(),
      cells_column_labels(everything())
    )
  ) %>% 
  tab_options(
    row_group.border.top.width = px(3),
    row_group.border.top.color = "black",
    row_group.border.bottom.color = "black",
    table_body.hlines.color = "white",
    table.border.top.color = "white",
    table.border.top.width = px(3),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )

gtsave(within_familiar_model_table, 'within_familiar_model_table.rtf')

pp_check(within_alldvs_familiar)
loo(within_alldvs_familiar)

# graph of familiarity by DV

familiarity_by_dv_graph <- within_alldvs_familiar %>%
  spread_draws(b_Intercept, r_familiar_5L[familiar_level,], r_questions[DV,]) %>%
  mean_qi(cond_mean = b_Intercept + r_familiar_5L + r_questions, .width = c(.95, .80)) %>%
  ungroup() %>%
  mutate(familiar_level = as.factor(familiar_level),
         familiar_level = fct_relevel(familiar_level, c('1', '2', '3', '4', '5')),
         familiar_level = fct_recode(familiar_level, `Not familiar` = '1',
                                     `Slightly familiar` = '2',
                                     `Moderately familiar` = '3',
                                     `Very familiar` = '4',
                                     `Substantially familiar` = '5')) %>%
  mutate(DV = case_when(DV == 'diff_abstract_aligned' ~ 'Abstract Aligned',
                        DV == 'diff_aligned' ~ 'Methods Aligned',
                        DV == 'diff_analysis_rigor' ~ 'Analysis Rigor',
                        DV == 'diff_did_learn' ~ 'Amt Learned',
                        DV == 'diff_discussion_quality' ~ 'Disc Quality',
                        DV == 'diff_field_importance' ~ 'Impt Discovery',
                        DV == 'diff_inspire' ~ 'Inspire Research',
                        DV == 'diff_intro_importance' ~ 'Impt Research',
                        DV == 'diff_method_rigor' ~ 'Methdos Rigor',
                        DV == 'diff_will_learn' ~ 'Amt will Learn',
                        DV == 'diff_method_quality' ~ 'Methdos Quality',
                        DV == 'diff_question_quality' ~ 'Quest Quality',
                        DV == 'diff_question_novel' ~ 'Quest Novelty',
                        DV == 'diff_method_creative' ~ 'Creative Methods',
                        DV == 'diff_overall_quality' ~ 'Overall Quality',
                        DV == 'diff_overall_import' ~ 'Impt Findings',
                        DV == 'diff_justificed' ~ 'Conclusion Justified',
                        DV == 'diff_result_quality' ~ 'Qualt Results',
                        DV == 'diff_result_innovative' ~ 'Innovative Result')) %>%
  ggplot(aes(y = familiar_level, x = cond_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ DV) +
  scale_x_continuous(breaks=seq(-1, 2, 1),
                     limits = c(-1, 2),
                     name = 'Difference between RR and non-RR articles') +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size=10))

familiarity_by_dv_graph

#### Model fit and summary for exploratory analysis of improve variable
within_improve_model_results <- summary(within_alldvs_improve)

within_improve_model_table <- rbind(within_improve_model_results$fixed %>%
                                       as.data.frame() %>% 
                                       rownames_to_column(var = "Predictor"), 
                                     within_improve_model_results$random$improve_6L %>% 
                                       as.data.frame() %>% 
                                       rownames_to_column(var = "Predictor") %>% 
                                       mutate(Predictor = 'sd(Intercept): Improve'),
                                     within_improve_model_results$random$participant_id %>% 
                                       as.data.frame() %>% 
                                       rownames_to_column(var = "Predictor") %>% 
                                       mutate(Predictor = 'sd(Intercept): Participant_id'),
                                     within_improve_model_results$random$questions %>%
                                       as.data.frame() %>% 
                                       rownames_to_column(var = "Predictor") %>% 
                                       mutate(Predictor = 'sd(Intercept): Question'),
                                     within_improve_model_results$random$RR %>%
                                       as.data.frame() %>% 
                                       rownames_to_column(var = "Predictor") %>% 
                                       mutate(Predictor = 'sd(Intercept): RR'),
                                     within_improve_model_results$spec_pars %>%
                                       as.data.frame() %>%
                                       rownames_to_column(var = "Predictor")) %>%
  mutate(effect_type = case_when(grepl('sd', Predictor) | grepl('sigma', Predictor) ~ 'Random Effects',
                                 TRUE ~ 'Fixed Effects')) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(-c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  rename(SE = 'Est.Error',
         lower_ci = "l-95% CI",
         upper_ci = "u-95% CI") %>%
  gt(groupname_col = 'effect_type') %>%
  tab_stubhead("label") %>% 
  cols_merge(columns = vars(lower_ci,upper_ci),
             hide_columns = vars(upper_ci),
             pattern = "[{1}, {2}]") %>%
  cols_label(lower_ci = '95% CrI',
             Predictor = '') %>%
  cols_align(align = 'center',
             columns = 2:4) %>%
  cols_align(align = 'left',
             columns = 1) %>%
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = list(
      cells_row_groups(),
      cells_column_labels(everything())
    )
  ) %>% 
  tab_options(
    row_group.border.top.width = px(3),
    row_group.border.top.color = "black",
    row_group.border.bottom.color = "black",
    table_body.hlines.color = "white",
    table.border.top.color = "white",
    table.border.top.width = px(3),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )

gtsave(within_improve_model_table, 'within_improve_model_table.rtf')
  
# graph of improvement by DV
improve_by_dv_graph <- within_alldvs_improve %>%
  spread_draws(b_Intercept, r_improve_6L[improve_level,], r_questions[DV,]) %>%
  mean_qi(cond_mean = b_Intercept + r_improve_6L + r_questions, .width = c(.95, .8)) %>%
  ungroup() %>%
  mutate(improve_level = as.factor(improve_level),
         improve_level = fct_relevel(improve_level, c('negative', 'neutral', 'slightly_more', 'moderately_more', 'much_more', 'substantially_more'))) %>%
  mutate(DV = case_when(DV == 'diff_abstract_aligned' ~ 'Abstract Aligned',
                        DV == 'diff_aligned' ~ 'Methods Aligned',
                        DV == 'diff_analysis_rigor' ~ 'Analysis Rigor',
                        DV == 'diff_did_learn' ~ 'Amt Learned',
                        DV == 'diff_discussion_quality' ~ 'Disc Quality',
                        DV == 'diff_field_importance' ~ 'Impt Discovery',
                        DV == 'diff_inspire' ~ 'Inspire Research',
                        DV == 'diff_intro_importance' ~ 'Impt Research',
                        DV == 'diff_method_rigor' ~ 'Methdos Rigor',
                        DV == 'diff_will_learn' ~ 'Amt will Learn',
                        DV == 'diff_method_quality' ~ 'Methdos Quality',
                        DV == 'diff_question_quality' ~ 'Quest Quality',
                        DV == 'diff_question_novel' ~ 'Quest Novelty',
                        DV == 'diff_method_creative' ~ 'Creative Methods',
                        DV == 'diff_overall_quality' ~ 'Overall Quality',
                        DV == 'diff_overall_import' ~ 'Impt Findings',
                        DV == 'diff_justificed' ~ 'Conclusion Justified',
                        DV == 'diff_result_quality' ~ 'Qualt Results',
                        DV == 'diff_result_innovative' ~ 'Innovative Result')) %>%
  ggplot(aes(y = improve_level, x = cond_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ DV) +
  scale_x_continuous(breaks=seq(-1, 2, 1),
                     limits = c(-1, 2),
                     name = 'Difference between RR and non-RR articles') +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size=10))

improve_by_dv_graph 
  
# exploration of guessing %>%
wide_data %>%
  group_by(guessed_right) %>%
  tally()

# full guessing model summary table

guessing_model_results <- summary(within_diff_pooled_guessed_model_slopes)

guessing_model_table <-  rbind(guessing_model_results$fixed %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = 'Predictor'),
                               guessing_model_results$random$RR %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor") %>% 
                                 mutate(Predictor = case_when(Predictor == 'sd(Intercept)' ~ 'sd(Intercept): RR')),
                               guessing_model_results$random$question %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor") %>% 
                                 mutate(Predictor = case_when(Predictor == 'sd(Intercept)' ~ 'sd(Intercept): question',
                                                              Predictor == 'sd(guessed_right2)' ~ 'sd(guessed_right2): question',
                                                              Predictor == 'sd(guessed_right3)' ~ 'sd(guessed_right3): question',
                                                              Predictor == 'cor(Intercept,guessed_right2)' ~ 'cor(Intercept,guessed_right2): question',
                                                              Predictor == 'cor(Intercept,guessed_right3)' ~ 'cor(Intercept,guessed_right3): question',
                                                              Predictor == 'cor(guessed_right2,guessed_right3)' ~ 'cor(guessed_right2,guessed_right3): question')),
                               guessing_model_results$random$participant_id %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor") %>% 
                                 mutate(Predictor = 'sd(Intercept): participant_id'),
                               guessing_model_results$spec_pars %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor")) %>%
  mutate(effect_type = case_when(grepl('sd', Predictor) | grepl('sigma', Predictor) | grepl('cor', Predictor) ~ 'Random Effects',
                                 TRUE ~ 'Fixed Effects')) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(-c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  rename(SE = 'Est.Error',
         lower_ci = "l-95% CI",
         upper_ci = "u-95% CI") %>%
  gt(groupname_col = 'effect_type') %>%
  tab_stubhead("label") %>% 
  cols_merge(columns = vars(lower_ci,upper_ci),
             hide_columns = vars(upper_ci),
             pattern = "[{1}, {2}]") %>%
  cols_label(lower_ci = '95% CrI',
             Predictor = '') %>%
  cols_align(align = 'center',
             columns = 2:4) %>%
  cols_align(align = 'left',
             columns = 1) %>%
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = list(
      cells_row_groups(),
      cells_column_labels(everything())
    )
  ) %>% 
  tab_options(
    row_group.border.top.width = px(3),
    row_group.border.top.color = "black",
    row_group.border.bottom.color = "black",
    table_body.hlines.color = "white",
    table.border.top.color = "white",
    table.border.top.width = px(3),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )

gtsave(guessing_model_table, 'guessing_model_table.rtf')  

# graph for guessing by DV

guessed_by_dv_graph <- posteriors_by_guessing %>%
  ungroup() %>%
  mutate(DV = case_when(DV == 'diff_abstract_aligned' ~ 'Abstract Aligned',
                        DV == 'diff_aligned' ~ 'Methods Aligned',
                        DV == 'diff_analysis_rigor' ~ 'Analysis Rigor',
                        DV == 'diff_did_learn' ~ 'Amt Learned',
                        DV == 'diff_discussion_quality' ~ 'Disc Quality',
                        DV == 'diff_field_importance' ~ 'Impt Discovery',
                        DV == 'diff_inspire' ~ 'Inspire Research',
                        DV == 'diff_intro_importance' ~ 'Impt Research',
                        DV == 'diff_method_rigor' ~ 'Methdos Rigor',
                        DV == 'diff_will_learn' ~ 'Amt will Learn',
                        DV == 'diff_method_quality' ~ 'Methdos Quality',
                        DV == 'diff_question_quality' ~ 'Quest Quality',
                        DV == 'diff_question_novel' ~ 'Quest Novelty',
                        DV == 'diff_method_creative' ~ 'Creative Methods',
                        DV == 'diff_overall_quality' ~ 'Overall Quality',
                        DV == 'diff_overall_import' ~ 'Impt Findings',
                        DV == 'diff_justificed' ~ 'Conclusion Justified',
                        DV == 'diff_result_quality' ~ 'Qualt Results',
                        DV == 'diff_result_innovative' ~ 'Innovative Result')) %>%
  ggplot(aes(y = guessed, x = mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ DV) +
  scale_x_continuous(breaks=seq(-1, 2, 1),
                     limits = c(-1, 2),
                     name = 'Difference between RR and non-RR articles') +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size=16))

guessed_by_dv_graph


##### full model report for between subjects ML model

# table of model results
btw_mlm_slopes_results <- summary(between_model_mlm_slopes)

btw_mlm_slopes_table <-  rbind(btw_mlm_slopes_results$fixed %>%
                                as.data.frame() %>%
                                 rownames_to_column(var = 'Predictor'),
                               btw_mlm_slopes_results$random$RR %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor") %>% 
                                 mutate(Predictor = case_when(Predictor == 'sd(Intercept)' ~ 'sd(Intercept): RR',
                                                               Predictor == 'sd(article_type2)' ~ 'sd(article_type2): RR',
                                                               Predictor == 'cor(Intercept,article_type2)' ~ 'cor(Intercept,article_type2): RR')),
                               btw_mlm_slopes_results$random$question %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor") %>% 
                                 mutate(Predictor = case_when(Predictor == 'sd(Intercept)' ~ 'sd(Intercept): question',
                                                              Predictor == 'sd(article_type2)' ~ 'sd(article_type2): question',
                                                              Predictor == 'cor(Intercept,article_type2)' ~ 'cor(Intercept,article_type2): question')),
                               btw_mlm_slopes_results$random$participant_id %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor") %>% 
                                 mutate(Predictor = 'sd(Intercept): participant_id'),
                               btw_mlm_slopes_results$spec_pars %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = "Predictor")) %>%
                          mutate(effect_type = case_when(grepl('sd', Predictor) | grepl('sigma', Predictor) | grepl('cor', Predictor) ~ 'Random Effects',
                                 TRUE ~ 'Fixed Effects')) %>%
                          mutate_if(is.numeric, round, 2) %>%
                          select(-c(Bulk_ESS, Tail_ESS, Rhat)) %>%
                          rename(SE = 'Est.Error',
                                 lower_ci = "l-95% CI",
                                 upper_ci = "u-95% CI") %>%
                          gt(groupname_col = 'effect_type') %>%
                          tab_stubhead("label") %>% 
                          cols_merge(columns = vars(lower_ci,upper_ci),
                                     hide_columns = vars(upper_ci),
                                     pattern = "[{1}, {2}]") %>%
                          cols_label(lower_ci = '95% CrI',
                                     Predictor = '') %>%
                          cols_align(align = 'center',
                                     columns = 2:4) %>%
                          cols_align(align = 'left',
                                     columns = 1) %>%
                          tab_style(
                            style = cell_text(color = "black", weight = "bold"),
                            locations = list(
                              cells_row_groups(),
                              cells_column_labels(everything())
                            )
                          ) %>% 
                          tab_options(
                            row_group.border.top.width = px(3),
                            row_group.border.top.color = "black",
                            row_group.border.bottom.color = "black",
                            table_body.hlines.color = "white",
                            table.border.top.color = "white",
                            table.border.top.width = px(3),
                            table.border.bottom.color = "white",
                            table.border.bottom.width = px(3),
                            column_labels.border.bottom.color = "black",
                            column_labels.border.bottom.width = px(2)
                          )

gtsave(btw_mlm_slopes_table, 'btw_mlm_slopes_table.rtf')  

# posterior predictive check graph
btw_mlm_slopes_pp_check_graph <- pp_check(between_model_mlm_slopes)
btw_mlm_slopes_pp_check_graph

# table of posteriors for each DV
btw_mlm_slopes_dv_posteriors <- posteriors_btw_mlm %>%
  ungroup() %>%
  filter(.width == 0.95) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(DV, article_effect, .lower, .upper) %>%
  gt() %>%
  cols_merge(columns = vars(.lower, .upper),
             hide_columns = vars(.upper),
             pattern = "[{1}, {2}]") %>%
  cols_label(article_effect = 'Posterior Mean',
             .lower = '95% CrI') %>%
  cols_align(align = 'center',
             columns = 2:3) %>%
  cols_align(align = 'left',
             columns = 1) %>%
  tab_style(
    style = cell_text(color = "black", weight = "bold"),
    locations = list(
      cells_column_labels(everything())
    )
  ) %>% 
  tab_options(
    table_body.hlines.color = "white",
    table.border.top.color = "white",
    table.border.top.width = px(3),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2)
  )

gtsave(btw_mlm_slopes_dv_posteriors, 'btw_mlm_slopes_dv_posteriors.rtf')
  