##loading libraries
library(osfr)
library(tidyverse)


#read in data



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


