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
                                              TRUE ~ 'half'))
            
                       

# set up contrasts codes
contrasts(wide_data$Field) <- contr.sum(3)
contrasts(wide_data$keyword_batch_comp) <- contr.sum(2)
contrasts(wide_data$Order) <- contr.sum(2)
contrasts(wide_data$Match) <- contr.sum(2)
contrasts(wide_data$improve_6L) <- contr.treatment(6)
contrasts(wide_data$familiar_5L) <- contr.treatment(5)
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


# Set up which model/prior/dv combinations to run for between models
between_models <- crossing(dv = long_data %>% select(question) %>% distinct(question) %>% pull(question),
                          set_priors = c(list(priors))) %>%
  mutate(between_pooled_model_results = pmap(list(dv, set_priors), between_pooled_model)) %>%
  mutate(posteriors = pmap(list(between_pooled_model_results, variable = dv, term = 'b_article_type2'), create_posteriors_term))



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

### within subjects models ###
within_diff_pooled_improve_model <- function(dv, set_priors) {
  within_model_diffs <- brm(as.formula(paste(dv, "~ Field + keyword_batch_comp + 
                                                      Order + Match + Order*Match +
                                                    (1|RR) + (1|improve_6L)")),
                            data = wide_data,
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_diffs)
}


# example with 1 dv & partial pooling across improve levels
diff_rigor_model <- brm(diff_analysis_rigor ~ Field + keyword_batch_comp + 
                          Order + Match + Order*Match +
                          (1|RR) + (1|improve_6L),
                        data = wide_data,
                        prior = priors,
                        family = 'gaussian',
                        chains = 4,
                        control = list(adapt_delta = 0.95))

diff_rigor_model %>%
  spread_draws(b_Intercept, r_improve_6L[improve_level,]) %>%
  median_qi(cond_mean = b_Intercept + r_improve_6L, .width = c(.95, .90)) %>%
  mutate(improve_level = fct_relevel(improve_level, c('negative', 'neutral', 'slightly_more', 'moderately_more', 'much_more', 'substantially_more'))) %>%
  ggplot(aes(y = improve_level, x = cond_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval()

diff_rigor_model %>% 
  spread_draws(r_improve_6L[improve_level,]) %>%
  compare_levels(r_improve_6L, by = improve_level) %>%
  ggplot(aes(y = improve_level, x = r_improve_6L)) +
  stat_halfeye()

# example with 1 dv & no partial pooling across improve levels
diff_rigor_model_fixed <- brm(diff_analysis_rigor ~ Field + keyword_batch_comp + 
                          Order + Match + Order*Match + improve_6L +
                          (1|RR),
                        data = wide_data,
                        prior = priors,
                        family = 'gaussian',
                        chains = 4,
                        control = list(adapt_delta = 0.95))


diff_rigor_model_fixed %>%
  spread_draws(b_Intercept) %>%
  median_qi(cond_mean = b_Intercept, .width = c(.95))

test <- brm(diff_analysis_rigor ~ Field + keyword_batch_comp + 
              Order + Match + Order*Match + believe_improve +
              (1|RR),
            data = wide_data,
            prior = priors,
            family = 'gaussian',
            chains = 4,
            control = list(adapt_delta = 0.95))

test %>%
  spread_draws(b_Intercept) %>%
  median_qi(cond_mean = b_Intercept, .width = c(.95))

## compare estimates from 2 models:

rbind(diff_rigor_model %>%
  spread_draws(b_Intercept, r_improve_6L[improve_level,]) %>%
  median_qi(cond_mean = b_Intercept + r_improve_6L, .width = c(.95, .90)) %>%
  mutate(model = 'RE variable'),
    rbind(diff_rigor_model_fixed %>% spread_draws(b_Intercept) %>%
            median_qi(cond_mean = b_Intercept, .width = c(.95, .90)) %>%
            mutate(improve_level = 'negative'),
          diff_rigor_model_fixed %>% spread_draws(b_Intercept, b_improve_6L2) %>%
            median_qi(cond_mean = b_Intercept + b_improve_6L2, .width = c(.95, .90)) %>%
            mutate(improve_level = 'neutral'),
          diff_rigor_model_fixed %>% spread_draws(b_Intercept, b_improve_6L3) %>%
            median_qi(cond_mean = b_Intercept + b_improve_6L3, .width = c(.95, .90)) %>%
            mutate(improve_level = 'slightly_more'),
          diff_rigor_model_fixed %>% spread_draws(b_Intercept, b_improve_6L4) %>%
            median_qi(cond_mean = b_Intercept + b_improve_6L4, .width = c(.95, .90)) %>%
            mutate(improve_level = 'moderately_more'),
          diff_rigor_model_fixed %>% spread_draws(b_Intercept, b_improve_6L5) %>%
            median_qi(cond_mean = b_Intercept + b_improve_6L5, .width = c(.95, .90)) %>%
            mutate(improve_level = 'much_more'),
          diff_rigor_model_fixed %>% spread_draws(b_Intercept, b_improve_6L6) %>%
            median_qi(cond_mean = b_Intercept + b_improve_6L6, .width = c(.95, .90)) %>%
            mutate(improve_level = 'substantially_more')) %>% mutate(model = 'FE variable')) %>%
  mutate(improve_level = as.factor(improve_level),
         improve_level = fct_relevel(improve_level, c('negative', 'neutral', 'slightly_more', 'moderately_more', 'much_more', 'substantially_more'))) %>%
  ggplot(aes(y = improve_level, x = cond_mean, xmin = .lower, xmax = .upper, color = model)) +
  geom_pointinterval(position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0) +
  theme_classic()


within_diff_pooled_familiar_model <- function(dv, set_priors) {
  within_model_diffs <- brm(as.formula(paste(dv, "~ Field + keyword_batch_comp + familiar_5L + 
                                                      Order + Match + Order*Match +
                                                    (1|RR)")),
                            data = wide_data,
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_diffs)
}

# model fit checks
summary(diff_rigor_model)
pp_check(diff_rigor_model)
WAIC(diff_rigor_model)
loo(diff_rigor_model)

summary(diff_rigor_model_fixed )
pp_check(diff_rigor_model_fixed )
WAIC(diff_rigor_model_fixed )
loo(diff_rigor_model_fixed )


# example with 1 dv & partial pooling across familiar levels
diff_rigor_model_familiar <- brm(diff_analysis_rigor ~ Field + keyword_batch_comp + 
                          Order + Match + Order*Match +
                          (1|RR) + (1|familiar_5L),
                        data = wide_data,
                        prior = priors,
                        family = 'gaussian',
                        chains = 4,
                        control = list(adapt_delta = 0.99))

diff_rigor_model_familiar %>%
  spread_draws(b_Intercept, r_familiar_5L[familiar_level,]) %>%
  median_qi(cond_mean = b_Intercept + r_familiar_5L, .width = c(.95, .90)) %>%
  ggplot(aes(y = familiar_level, x = cond_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval()

diff_rigor_model_familiar %>% 
  spread_draws(r_familiar_5L[familiar_level,]) %>%
  compare_levels(r_familiar_5L, by = familiar_level) %>%
  ggplot(aes(y = familiar_level, x = r_familiar_5L)) +
  stat_halfeye()


# example with 1 dv & no partial pooling across improve familiar
diff_rigor_model_familiar_fixed <- brm(diff_analysis_rigor ~ Field + keyword_batch_comp + 
                                Order + Match + Order*Match + familiar_5L +
                                (1|RR),
                              data = wide_data,
                              prior = priors,
                              family = 'gaussian',
                              chains = 4,
                              control = list(adapt_delta = 0.99))


diff_rigor_model_familiar_fixed %>%
  spread_draws(b_Intercept) %>%
  median_qi(cond_mean = b_Intercept, .width = c(.95))



within_diff_pooled_guessed_model <- function(dv, set_priors, guessed) {
  within_model_diffs <- brm(as.formula(paste(dv, "~ Field + keyword_batch_comp + guessed_right 
                                                      Order + Match + Order*Match +
                                                    (1|RR)")),
                            data = wide_data,
                            prior = set_priors, 
                            family = 'gaussian',
                            chains = 4)
  return(within_model_diffs)
}


rbind(diff_rigor_model_familiar %>%
        spread_draws(b_Intercept, r_familiar_5L[familiar_level,]) %>%
        median_qi(cond_mean = b_Intercept + r_familiar_5L, .width = c(.95, .90)) %>%
        mutate(model = 'RE variable'),
      rbind(
            diff_rigor_model_familiar_fixed %>% spread_draws(b_Intercept) %>%
              median_qi(cond_mean = b_Intercept, .width = c(.95, .90)) %>%
              mutate(familiar_level = '1'),
            diff_rigor_model_familiar_fixed %>% spread_draws(b_Intercept, b_familiar_5L2) %>%
              median_qi(cond_mean = b_Intercept + b_familiar_5L2, .width = c(.95, .90)) %>%
              mutate(familiar_level = '2'),
            diff_rigor_model_familiar_fixed %>% spread_draws(b_Intercept, b_familiar_5L3) %>%
              median_qi(cond_mean = b_Intercept + b_familiar_5L3, .width = c(.95, .90)) %>%
              mutate(familiar_level = '3'),
            diff_rigor_model_familiar_fixed %>% spread_draws(b_Intercept, b_familiar_5L4) %>%
              median_qi(cond_mean = b_Intercept + b_familiar_5L4, .width = c(.95, .90)) %>%
              mutate(familiar_level = '4'),
            diff_rigor_model_familiar_fixed %>% spread_draws(b_Intercept, b_familiar_5L5) %>%
              median_qi(cond_mean = b_Intercept + b_familiar_5L5, .width = c(.95, .90)) %>%
              mutate(familiar_level = '5')) %>% mutate(model = 'FE variable')) %>%
  ggplot(aes(y = familiar_level, x = cond_mean, xmin = .lower, xmax = .upper, color = model)) +
  geom_pointinterval(position=position_dodge(width=0.5)) +
  geom_vline(xintercept = 0) +
  theme_classic()

### EFA of diff DVs to investigate potential exchangability

#initial correlations
wide_data %>%
  select(starts_with('diff')) %>%
  cor(use = "pairwise.complete.obs") %>%
  corrplot(type = 'upper', order="hclust", method = 'number')

wide_data %>%
  select(participant_id, starts_with('diff')) %>%
  column_to_rownames('participant_id') %>%
  fa.parallel()

wide_data %>%
  select(participant_id, starts_with('diff')) %>%
  select(-c(diff_abstract_aligned, diff_intro_importance, diff_question_novel)) %>%
  column_to_rownames('participant_id') %>%
  fa.parallel()

efa3 <- wide_data %>%
          select(participant_id, starts_with('diff')) %>%
          column_to_rownames('participant_id') %>%
          fa(nfactors = 3, rotate = 'oblimin') 
efa3
fa.diagram(efa3)

efa3_cutdown <- wide_data %>%
  select(participant_id, starts_with('diff')) %>%
  select(-c(diff_abstract_aligned, diff_intro_importance, diff_question_novel)) %>%
  column_to_rownames('participant_id') %>%
  fa(nfactors = 3, rotate = 'oblimin') 

efa2_cutdown <- wide_data %>%
  select(participant_id, starts_with('diff')) %>%
  select(-c(diff_abstract_aligned, diff_intro_importance, diff_question_novel)) %>%
  column_to_rownames('participant_id') %>%
  fa(nfactors = 2, rotate = 'oblimin') 

efa3_cutdown
efa2_cutdown
fa.diagram(efa3_cutdown)
fa.diagram(efa2_cutdown)

### within model with DVs as ML component
mlm_dvs_data <- wide_data %>%
  select(starts_with('diff'), participant_id, RR, Field, keyword_batch_comp, Order, Match, improve_6L, familiar_5L) %>%
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

summary(within_alldvs)
pp_check(within_alldvs)
WAIC(within_alldvs)
loo(within_alldvs)

# main figure
with_alldvs_graph_nums <- within_alldvs %>%
  spread_draws(b_Intercept, r_questions[dv,]) %>%
  mutate(dv_estimates = b_Intercept + r_questions) %>%
  left_join(within_alldvs %>%
         spread_draws(b_Intercept, r_questions[dv,]) %>% 
         mean_qi(dv_mean = b_Intercept + r_questions) %>% 
         select(dv, dv_mean), by = 'dv') %>%
  ungroup()

main_graph_creation <- function(data) {
  data %>%  
    mutate(dv = as.factor(dv),
           dv = fct_reorder(dv, dv_median)) %>%
    ggplot(aes(y = dv, x = dv_estimates, fill = stat(x <= 0))) +
    stat_halfeye(point_interval = mean_qi, .width = c(.95, .8)) +
    scale_x_continuous(breaks=seq(-.5, 1.5, .5),
                       limits = c(-.75, 1.75),
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


combined_plot <- intro_qs / results_qs / abstract_qs + plot_layout(heights = c(8, 7, 4))
combined_plot



scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
                                               width = 24)) +

within_alldvs %>%
  spread_draws(b_Intercept, r_questions[dv,]) %>%
  mutate(dv_estimates = b_Intercept + r_questions) %>%
  left_join(within_alldvs %>%
              spread_draws(b_Intercept, r_questions[dv,]) %>% 
              median_qi(dv_median = b_Intercept + r_questions) %>% 
              select(dv, dv_median), by = 'dv') %>%
  ungroup() %>%
  mutate(article_section = case_when(dv == 'diff_inspire' |
                                       dv == 'diff_abstract_aligned' |
                                       dv == 'diff_field_importance' |
                                       dv == "diff_overall_quality" ~ 'Abstract',
                                     grepl('result', dv) |
                                       dv == "diff_analysis_rigor" |
                                       dv == "diff_overall_import" |
                                       dv == "diff_did_learn" |
                                       dv == "diff_discussion_quality" |
                                       dv == "diff_justificed" ~ 'Results/Discussion',
                                     grepl('question', dv) |
                                       grepl('method', dv) |
                                       dv == 'diff_aligned' | 
                                       dv == 'diff_will_learn' | 
                                       dv == 'diff_intro_importance' ~ 'Intro/Methods')) %>%
  mutate(dv = as.factor(dv),
         dv = fct_reorder(dv, dv_median),
         article_section = as.factor(article_section),
         article_section = fct_relevel(article_section, c('Intro/Methods', 'Results/Discussion', 'Abstract'))) %>%
  ggplot(aes(y = dv, x = dv_estimates, fill = stat(x <= 0))) +
  stat_halfeye(.width = c(.95, .8)) +
  scale_x_continuous(breaks=seq(-.5,2,.5)) +
  facet_wrap(~ article_section, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("skyblue", "gray80")) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title = element_blank())

within_alldvs %>%
  spread_draws(b_Intercept, r_questions[dv,]) %>%
  mutate(dv_estimates = b_Intercept + r_questions) %>%
  left_join(within_alldvs %>%
              spread_draws(b_Intercept, r_questions[dv,]) %>% 
              median_qi(dv_median = b_Intercept + r_questions) %>% 
              select(dv, dv_median), by = 'dv') %>%
  ungroup() %>%
  mutate(article_section = case_when(dv == 'diff_inspire' |
                                       dv == 'diff_abstract_aligned' |
                                       dv == 'diff_field_importance' |
                                       dv == "diff_overall_quality" ~ 'After finishing the paper',
                                     grepl('result', dv) |
                                       dv == "diff_analysis_rigor" |
                                       dv == "diff_overall_import" |
                                       dv == "diff_did_learn" |
                                       dv == "diff_discussion_quality" |
                                       dv == "diff_justificed" ~ 'After knowing study outcomes',
                                     grepl('question', dv) |
                                       grepl('method', dv) |
                                       dv == 'diff_aligned' | 
                                       dv == 'diff_will_learn' | 
                                       dv == 'diff_intro_importance' ~ 'Before knowing study outcomes')) %>%
  mutate(dv = as.factor(dv),
         dv = fct_reorder(dv, dv_median),
         article_section = as.factor(article_section),
         article_section = fct_relevel(article_section, c('Before knowing study outcomes', 'After knowing study outcomes', 'After finishing the paper'))) %>%
  ggplot(aes(y = dv, x = dv_estimates, fill = stat(x <= 0))) +
  stat_halfeye(.width = c(.95, .8)) +
  scale_x_continuous(breaks=seq(-.5,2,.5)) +
  facet_grid(rows = vars(article_section), scales = "free_y", space = 'free_y', switch = "y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("skyblue", "gray80")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        strip.text.y = element_text(angle = 270, face = "bold"),
        strip.placement = "outside",
        axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
  

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
  


# partial pooling across all DVs and improve

within_alldvs_improve <-  brm(response ~ Field + keyword_batch_comp + 
                        Order + Match + Order*Match +
                        (1|RR) + (1|participant_id) + (1|questions) + (1|improve_6L),
                      data = mlm_dvs_data,
                      prior = priors,
                      family = 'gaussian',
                      chains = 4,
                      iter = 3000,
                      control = list(adapt_delta = 0.99, max_treedepth = 15))

summary(within_alldvs_improve)
pp_check(within_alldvs_improve)
WAIC(within_alldvs_improve)
loo(within_alldvs_improve)

within_alldvs_improve %>%
  spread_draws(b_Intercept, r_improve_6L[improve_level,], r_questions[DV,]) %>%
  median_qi(cond_mean = b_Intercept + r_improve_6L + r_questions, .width = c(.95, .90)) %>%
  ungroup() %>%
  mutate(improve_level = as.factor(improve_level),
         improve_level = fct_relevel(improve_level, c('negative', 'neutral', 'slightly_more', 'moderately_more', 'much_more', 'substantially_more'))) %>%
  ggplot(aes(y = DV, x = cond_mean, xmin = .lower, xmax = .upper, color = improve_level)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  theme_classic()

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
                              iter = 3000,
                              control = list(adapt_delta = 0.99, max_treedepth = 15))

summary(within_alldvs_familiar)
pp_check(within_alldvs_familiar)
WAIC(within_alldvs_familiar)
loo(within_alldvs_familiar)

within_alldvs_familiar %>%
  spread_draws(b_Intercept, r_familiar_5L[familiar_level,], r_questions[DV,]) %>%
  median_qi(cond_mean = b_Intercept + r_familiar_5L + r_questions, .width = c(.95, .90)) %>%
  ungroup() %>%
  mutate(familiar_level = as.factor(familiar_level),
         familiar_level = fct_relevel(familiar_level, c('1', '2', '3', '4', '5'))) %>%
  ggplot(aes(y = DV, x = cond_mean, xmin = .lower, xmax = .upper, color = familiar_level)) +
  geom_pointinterval(position=position_dodge(width=1)) +
  geom_vline(xintercept = 0) +
  theme_classic()

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
