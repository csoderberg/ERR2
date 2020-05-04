##loading libraries
library(osfr)
library(tidyverse)
library(here)


#### Creating deidentified raw data set ####

#read in raw data
raw_numeric <- read_csv(here::here('/ERR2_rawdata_numeric.csv'))

#how many will be excluded?
raw_numeric %>%
  slice(-1:-2) %>% #take out first 2 rows added by qualtrics
  filter(FirstFamiliar == '1' | SecondFamiliar == '1') %>%
  nrow()
  
# exclude anyone who got one of their own papers (in line with pre-reg)
raw_deidentified_numeric <- raw_numeric %>%
                                slice(-1:-2) %>% #take out first 2 rows added by qualtrics
                                filter(FirstFamiliar != '1' & SecondFamiliar != '1')

# remove free response variables in line with IRB for public posting on OSF
raw_deidentified_numeric <- raw_deidentified_numeric %>%
                              select(-c(RRintrocomment, RRresultcomment, RRabstractcomment, Altintrocomment, Altresultcomment, Altabstractcomment,
                                        `Altresultcomment - Parent Topics`, `Altresultcomment - Topics`))

# created and upload deidentified raw data
write_csv(raw_deidentified_numeric, 'deidentified_raw_numeric.csv')
osf_upload('https://osf.io/q6pef/', 'deidentified_raw_numeric.csv')


#### Creating Cleaned data file ####

#download & read in deidentifed raw data
osf_retrieve_file("") %>%
  osf_download()

ERR2_numeric_data <- read_csv(here::here('/deidentified_raw_numeric.csv'))

#count how many opened survey but didn't fill out field question
sum(is.na(ERR2_numeric_data))

#count how many did not meet field criterion
sum(ERR2_numeric_data$Field == 4)

#count how many consent by didn't answer any questions after that
ERR2_numeric_data %>%
  filter(Consent == 2) %>%
  mutate(num_NA = rowSums(is.na(ERR2_numeric_data[,15:81])))


#filter out those that never consented & select only pertinent variables
numeric_data <- ERR2_numeric_data %>%
                    filter(Consent == 2) %>%
                    select(-c(Status, Progress, `Duration (in seconds)`, Finished, RecordedDate,
                              ResponseId, DistributionChannel, UserLanguage, Consent, 
                              FromLink, `Altresultcomment - Parent Topics`, `Altresultcomment - Topics`))

#create difference scores
numeric_data_wide <- numeric_data %>%
                    mutate(diff_question_quality = RRQuestionQuality - AltQuestionQuality,
                           diff_method_quality = RRMethodQuality - AltMethodQuality,
                           diff_method_rigor = RRMethodRigor - AltMethodRigor,
                           diff_intro_importance = RRIntroImportance - AltIntroImportance,
                           diff_question_novel = RRQuestionNovelty - AltQuestionNovelty,
                           diff_method_creative = RRMethodCreative - AltMethodCreative,
                           diff_will_learn = RRWillLearn - AltWillLearn,
                           diff_aligned = RRAligned - AltAligned,
                           diff_analysis_rigor = RRAnalysisRigor - AltAnalysisRigor,
                           diff_overall_import = RROverallImport - AltOverallImport,
                           diff_result_innovative = RRResultInnovative - AltResultInnovative,
                           diff_did_learn = RRDidLearn - AltDidLearn,
                           diff_justificed = RRJustified - AltJustified,
                           diff_result_quality = RRResultQuality - AltResultQuality,
                           diff_discussion_quality = RRDiscussionQuality - AltDiscussionQuality,
                           diff_abstract_aligned = RRAbstractAligned - AltAbstractAligned,
                           diff_field_importance = RRFieldImportance - AltFieldImportance,
                           diff_inspire = RRInspire - AltInspire,
                           diff_overall_quality = RROverallQuality - AltOverallQuality)

# write out clean data file in wide format
write_csv(numeric_data_wide, 'cleaned_numeric_data_wide.csv')

# make long_format file
numeric_data_long <- numeric_data %>%
                        pivot_longer(cols = RRQuestionQuality:AltOverallQuality, names_to = 'question', values_to = 'response') %>%
                        mutate(article_type = case_when(grepl('^RR', question) ~ 'RR',
                                                        grepl('^Alt', question) ~ 'nonRR'))

# write out clean data file
write_csv(numeric_data_long, 'cleaned_numeric_data_long.csv')



