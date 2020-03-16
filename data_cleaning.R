##loading libraries
library(osfr)
library(tidyverse)
library(here)


#### Creating deidentified raw data set ####

#download raw data
osf_retrieve_file("https://osf.io/wjkg6/") %>%
  osf_download()

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

# created and upload deidentified raw data
write_csv(raw_deidentified_numeric, 'deidentified_raw_numeric.csv')
osf_upload('https://osf.io/q6pef/', oath = 'deidentified_raw_numeric.csv')


#### Creating Cleaned data file ####

#download & read in deidentifed raw data
osf_retrieve_file("") %>%
  osf_download()

deidentified_raw_numeric <- read_csv(here::here('/deidentified_raw_numeric.csv'))


#create difference scores
numeric_data <- numeric_data %>%
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


