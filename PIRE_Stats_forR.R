####  User Inputs ####
# rm(list=ls())
# setwd(dir = "/Users/TEST/Texas A&M University - Corpus Christi/Bird, Chris - MartinFrenchMasters/sandbox/")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("../Final Versions/100K_Barcodes/")
# setwd("./")

Aduo_MuseInFile <- "./Aduo_100K_MuseumBarcodes.csv"
Aduo_ContInFile <- "./Aduo_100K_ContemporaryBarcodes.csv"

Sspi_MuseInFile <- "./Sspi_100K_MuseumBarcodes.csv"
Sspi_ContInFile <- "./Sspi_100K_ContemporaryBarcodes.csv"

library(tidyverse)
library(magrittr)
library(splines)
library(lme4)
library(afex)
library(emmeans)
library(expss)
library(car)
library(betareg)
library(gridExtra)
library(RColorBrewer)
library(DiffXTables)
library(rstatix)

#### Make Aduo_sequencing_error_data tibble ####
sink('./Aduo_ErrorCalc_Analysis_Output.txt')
writeLines("> Species = Aduo")
writeLines(c("\n","####################################","\n"))
sink()

Aduo_Museum <- read_csv(Aduo_MuseInFile, col_types = str_c(rep('c', 20), collapse = '')) 
Aduo_Contemporary <- read_csv(Aduo_ContInFile, col_types = str_c(rep('c', 20), collapse = ''))

Aduo_Both <- bind_rows(Aduo_Museum, Aduo_Contemporary, .id = 'Sampling')

a<-1; b<-1

Aduo_sequencing_error_data <- bind_rows(Aduo_Museum, Aduo_Contemporary, .id = 'Sampling') %>%
  mutate(Sampling = case_when(Sampling == '1' ~ 'Museum',
                              Sampling == '2' ~ 'Contemporary',
                              TRUE ~ 'error')) %>%
  mutate(Barcode = str_c('GG', Barcode, "GG", sep = '')) %>% 
  gather(position, base, -Sampling:-Indels)  %>% 
  mutate(groupings = case_when(position == 'BP1' | position == 'BP2' ~ 'GG',
                               position == 'BP3' | position == 'BP4' | position == 'BP5' | position == 'BP6' | position == 'BP7' | position == 'BP8' | position == 'BP9' | position == 'BP10' ~ 'Barcode',
                               position == 'BP11' | position == 'BP12' | position == 'BP13' | position == 'BP14' ~ 'Ligation Site',
                               position == 'BP15' | position == 'BP16' ~ 'End Cut Site',
                               TRUE ~ 'error')) %>%
  mutate(base_position = str_extract(position, "[0-9]+") %>% as.integer,
         true_base = str_sub(Barcode, base_position, base_position),
         true_base = if_else(true_base == "", NA_character_, true_base)) %>%
  group_by(Sampling, groupings, base_position, Barcode, fqBase, Read, Indels) %>%
  summarise(correct = sum(base == true_base),
            error = sum(base != true_base),
            total = n()) %>%
  ungroup %>%
  mutate_if(is.character, as.factor) %>%
  mutate(base_position = as.factor(base_position)) %>%
  mutate(base_position = factor(base_position, levels = 1:16, ordered = TRUE),
         Indels = factor(Indels, levels = c("2del","1del","0Ind","1Ins","2Ins","3Ins","4Ins","5Ins","6Ins","7Ins","8Ins"), ordered = TRUE),
         IndelCategory = factor(Indels, levels = c("Insertions","NoIndels","Deletions"), ordered = TRUE),
         groupings = factor(groupings, levels = c("Barcode","Ligation Site","End Cut Site"), ordered = TRUE)) %>%
  filter(fqBase != "Sub-100000-20190528-PIRE-Adu-A-P1P4-L008", groupings != "GG") 

Aduo_sequencing_error_data %>% mutate(error_rate = error/total) -> Aduo_sequencing_error_data
Aduo_sequencing_error_data %>% mutate(success_rate = correct/total) -> Aduo_sequencing_error_data
Aduo_sequencing_error_data$Indels <- factor(Aduo_sequencing_error_data$Indels, levels = c("2del", "1del", "0Ind", "1Ins", "2Ins", "3Ins", "4Ins", "5Ins", "6Ins", "7Ins", "8Ins"))

Aduo_sequencing_error_data$IndelCategory <- gsub(".Ins", "Insertions", Aduo_sequencing_error_data$Indels)
Aduo_sequencing_error_data$IndelCategory <- gsub(".del", "Deletions", Aduo_sequencing_error_data$IndelCategory)
Aduo_sequencing_error_data$IndelCategory <- gsub("0Ind", "NoIndels", Aduo_sequencing_error_data$IndelCategory)
Aduo_sequencing_error_data$IndelCategory <- factor(Aduo_sequencing_error_data$IndelCategory, levels = c("Insertions","NoIndels","Deletions"))

#### Plot mean error rate using Aduo_sequencing_error_data ####
# summarized data to make data labels
Aduo_SummaryData <- Aduo_sequencing_error_data %>% 
  group_by(Sampling, groupings, fqBase, IndelCategory) %>%
  summarise(correct = sum(correct),
            error = sum(error),
            total = sum(total)) %>% 
  mutate(error_rate = error/total) %>% 
  group_by(Sampling, groupings, IndelCategory) %>%
  summarize(mean_error_rate = mean(error_rate)) 

Aduo_SummaryData$Rounded_mean_error_rate <- sprintf("%0.3f", round(Aduo_SummaryData$mean_error_rate, 3))

ggplot(Aduo_SummaryData, aes(x = groupings, y = mean_error_rate, fill = Sampling)) +
  geom_bar(stat = "identity") +
  ggtitle('Aduo Mean Error Rate') +
  theme_classic() +
  geom_text(label = Aduo_SummaryData$Rounded_mean_error_rate, nudge_y = .15) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0), axis.title.x = element_blank()) +
  facet_grid(IndelCategory ~ Sampling)
ggsave(filename ='./Aduo_MeanErrorRate.png', width = 6, height = 6)

# this version is to make combined plot
Aduo_MeanErrorRate <- ggplot(Aduo_SummaryData, aes(x = groupings, y = mean_error_rate)) +
  geom_bar(stat = "identity", fill = "lightgrey") +
  ggtitle('Aduo') +
  theme_classic() +
  geom_text(label = Aduo_SummaryData$Rounded_mean_error_rate, size = 3, nudge_y = .15) +
  theme(title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12)) +
  facet_grid(IndelCategory ~ Sampling)


#### Make Collapsed Aduo_sequencing_error_data_2 ####
Aduo_sequencing_error_data_2 <- Aduo_sequencing_error_data %>%
  ungroup() %>%
  group_by(Sampling, groupings, fqBase, IndelCategory, Barcode, Read) %>%
  summarize(correct = sum(correct),
            error = sum(error),
            total = sum(total),
            correct_2 = sum(correct) + 1,
            error_2 = sum(error) + 1,
            total_2 = sum(total) + 2) %>%
  mutate(error_rate = error / total,
         success_rate = correct / total,
         error_rate_2 = error_2 / total_2,
         success_rate_2 = correct_2 / total_2)

# Aduo_error_model_2 <- mixed(cbind(correct_2, error_2) ~ Sampling * groupings * IndelCategory
#                        + (1 | fqBase) + (1 | Barcode) + (1 | Read),
#                        data = Aduo_sequencing_error_data_2,
#                        family = binomial,
#                        method = 'LRT',
#                        control = glmerControl(optimizer = "bobyqa",
#                                             optCtrl = list(maxfun = 2e8)))

Aduo_error_model_2 <- mixed(cbind(correct_2, error_2) ~ Sampling * groupings * IndelCategory + (1 | fqBase) + (1 | Barcode) + (1 | Read),
                            data = Aduo_sequencing_error_data_2,
                            family = binomial,
                            method = 'LRT',
                            control = glmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 2e9),
                                                   check.conv.grad = .makeCC("warning", tol = 3e-3, relTol = NULL)))


#### Make Collapsed Aduo_sequencing_error_data_3 ####
# collapsing random factors that we don't care about and should have no effect on the results
# one warning is produced by this model but CEB and JS determined it was ok
Aduo_sequencing_error_data_3 <- Aduo_sequencing_error_data %>%
  separate(fqBase, remove = FALSE, sep = "-", into = c(NA, NA, NA, NA, NA, NA, "plate_pool", NA)) %>%
  mutate(library = str_c(plate_pool, Barcode, sep = "_")) %>%
  ungroup() %>%
  group_by(Sampling, groupings, library, IndelCategory, Read) %>%
  summarize(correct = sum(correct),
            error = sum(error),
            total = sum(total),
            correct_2 = sum(correct) + 1,
            error_2 = sum(error) + 1,
            total_2 = sum(total) + 2) %>%
  mutate(error_rate = error / total,
         success_rate = correct / total,
         error_rate_2 = error_2 / total_2,
         success_rate_2 = correct_2 / total_2) 

#### Make Aduo_error_model_3 ####
Aduo_error_model_3 <- mixed(cbind(correct_2, error_2) ~ Sampling * groupings * IndelCategory 
                       + (1 | library) + (1 | Read),
                       data = Aduo_sequencing_error_data_3, 
                       family = binomial, 
                       method = 'LRT',
                       control = glmerControl(optimizer = "bobyqa", 
                                            optCtrl = list(maxfun = 2e8),
                                            check.conv.grad = .makeCC("warning", tol = 3e-3, relTol = NULL)))
Aduo_error_model_2 <- Aduo_error_model_3

# estimated means from model & test pairwaise comps
Aduo_error_emm <- emmeans(Aduo_error_model_2, ~ Sampling * groupings * IndelCategory)

# using BH becuase it does the best job of estimating false discovery rate.
Aduo_pairwise_postHoc <- contrast(Aduo_error_emm, 
                                  method = 'pairwise', 
                                  simple = 'each', 
                                  combine = FALSE, 
                                  adjust = "bonferroni")

# write statistical result table to output
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n","####################################","\n","\n","> Aduo_error_model"))
Aduo_error_model_2
writeLines(c("\n","####################################","\n","\n","> Aduo_pairwise_postHoc"))
Aduo_pairwise_postHoc
sink()

# format model results
Aduo_tbl_emm <- emmeans(Aduo_error_model_2, 
                        ~ Sampling * groupings * IndelCategory, 
                        type = 'response') %>%
  as_tibble %>%
  mutate(IndelCategory = factor(IndelCategory, levels = c('Insertions', 
                                                          'NoIndels', 
                                                          'Deletions')),
         groupings = factor(groupings, levels = c("GG",
                                                  "Barcode",
                                                  "Ligation Site",
                                                  "End Cut Site"), ordered = TRUE))

# plot model results
Aduo_p_emm <- Aduo_tbl_emm %>%
  ggplot(aes(x = groupings, y = prob, ymin = asymp.LCL, ymax = asymp.UCL,
             fill = Sampling)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.2) +
  geom_col(position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ IndelCategory)
Aduo_p_emm
ggsave(filename = './Aduo_p_emm.png', width = 6, height = 6)

# format means from data
Aduo_tbl_means <- Aduo_sequencing_error_data_2 %>%
  ungroup(Barcode, fqBase) %>%
  summarize(mean_success_rate = mean(success_rate),
            mean_success_rate_2 = mean(success_rate_2)) %>%
  mutate(groupings = factor(groupings, levels = c("GG",
                                                  "Barcode",
                                                  "Ligation Site",
                                                  "End Cut Site"), ordered = TRUE))

# plot model results
Aduo_p_means <- Aduo_tbl_means %>%
  ggplot(aes(x = groupings, y = mean_success_rate, fill = Sampling)) +
  geom_col(position = position_dodge(0.9)) +
  geom_point(aes(x = groupings, y = mean_success_rate_2),
             position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ IndelCategory)
Aduo_p_means
ggsave(filename ='./Aduo_p_means.png', width = 6, height = 6)

# plot estimated marginal means from model and means from data
png(file = './Aduo_p_emm_and_Aduo_p_means.png', width = 6, height = 6, units = "in", res = 500)
grid.arrange(Aduo_p_emm, Aduo_p_means, nrow=2,
             top = 'Aduo_p_emm_and_Aduo_p_means')
dev.off()

# combine and compare means to estimated marginal means
Aduo_tbl_emm_means <- Aduo_tbl_emm %>%
  full_join(Aduo_tbl_means,
            by=c("Sampling", "groupings", "IndelCategory")) %>%
  rename(emm = prob,
         means = mean_success_rate)

# plot comps between means and emm
Aduo_means_VS_emm <- Aduo_tbl_emm_means %>%
  ggplot(aes(x = means, y = emm)) +
  geom_point() +
  geom_abline() +
  geom_smooth(color = "blue")

Aduo_emm_VS_residuals <- Aduo_tbl_emm_means %>%
  mutate(residuals = emm - means) %>%
  ggplot(aes(x = emm, y = residuals)) +
  geom_point() +
  geom_abline() +
  geom_smooth(color = "blue")


png(file = './Aduo_emm_and_means.png', width = 6, height = 6, units = "in", res = 500)
grid.arrange(Aduo_means_VS_emm, Aduo_emm_VS_residuals, nrow = 2,
             top = 'Aduo_emm_and_means')
dev.off()

# this is the plot of emm with the raw mean prop successful and mean prop successful 2 (which has 1 added to the correct and error values)
Aduo_tbl_emm_means %>%
  ggplot(aes(x = groupings, y = emm, ymin = asymp.LCL, ymax = asymp.UCL,
             fill = Sampling)) +
  ggtitle('Aduo emm mean prop successful 1&2') +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.2, color = "grey") +
  geom_point(aes(x = groupings, y = mean_success_rate_2),
             position = position_dodge(0.9),
             shape = 1,
             size = 10) +
  geom_point(aes(x = groupings, y = means),
             position = position_dodge(0.9),
             color = "green",
             shape = 1,
             size = 10) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.title = element_blank()) +
  facet_wrap(. ~ IndelCategory)
ggsave(filename = './Aduo_emm_mean-prop-successful_and_mean-prop-successful-2.png', width = 8, height = 6)

Aduo_error_model <- Aduo_error_model_2

sum_Aduo_error_emm <- Aduo_error_emm %>%
  summary(type = 'response') %>%
  as_tibble %>%
  ggplot(aes(x = groupings, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, colour = Sampling)) +
  geom_linerange(position = position_dodge(0.5), size = 1) +
  geom_point(position = position_dodge(0.5)) +
  labs(title = 'Aduo Cutsite Proportion Correct by Group', x = "DNA Identity", y = "Proportion Correct") +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  theme_classic() +
  facet_grid(IndelCategory ~ .) +
  theme(legend.title = element_blank())

sum_Aduo_error_emm
ggsave(filename= './Aduo_CutsiteGroupProportionCorrect.png', width = 6, height = 6)

#### Indel Proportion Plots ####
Aduo_IndelCounts <- table(Aduo_Both$fqBase, Aduo_Both$Indels)

Aduo_IndelCountsTidy <- rbind(data.frame(table(Aduo_Museum$fqBase, Aduo_Museum$Indels)), data.frame(table(Aduo_Contemporary$fqBase, Aduo_Contemporary$Indels)))
colnames(Aduo_IndelCountsTidy) <- c("Fqbase", "Indel", "Freq")

sink('./Aduo_ErrorCalc_Analysis_Output.txt', append=TRUE)
writeLines(c("\n","####################################","\n","> Aduo_IndelCounts"))
Aduo_IndelCounts
sink()

Indels <- data.frame(rep(c("2del", "1del", "0Ind", "1Ins", "2Ins", "3Ins", "4Ins", "5Ins", "6Ins", "7Ins", "8Ins"), 2))

Aduo_Museum_total <- count(Aduo_Museum)
Aduo_Contemporary_total <- count(Aduo_Contemporary)
Aduo_Both_total <- count(Aduo_Both)

Proportion <- as.numeric(paste(data.frame(c(round((count(Aduo_Museum[Aduo_Museum$Indels=="2del",])/Aduo_Museum_total*100),3),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="1del",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="0Ind",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="1Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="2Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="3Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="4Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="5Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="6Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="7Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Museum[Aduo_Museum$Indels=="8Ins",])/Aduo_Museum_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="2del",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="1del",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="0Ind",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="1Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="2Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="3Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="4Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="5Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="6Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="7Ins",])/Aduo_Contemporary_total*100),3)),
                                            (round((count(Aduo_Contemporary[Aduo_Contemporary$Indels=="8Ins",])/Aduo_Contemporary_total*100),3)))),sep = "\n"))

Group <- data.frame(c(rep("Museum", 11), rep("Contemporary", 11)))

Aduo_Indel_Proportions_Table <- data.frame(Indels, Proportion, Group)
colnames(Aduo_Indel_Proportions_Table) <- c("Indels", "Prop", "Group")

Aduo_Indel_Proportions_Table$Indels <- factor(Aduo_Indel_Proportions_Table$Indels, ordered = TRUE, levels = c("2del", "1del", "0Ind", "1Ins", "2Ins" ,"3Ins" ,"4Ins" ,"5Ins" ,"6Ins" ,"7Ins" ,"8Ins"))

sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n","####################################","\n","> Aduo_Indel_Proportions_Table"))
Aduo_Indel_Proportions_Table
sink()

ggplot(Aduo_Indel_Proportions_Table, aes(x = Indels, y = Proportion, fill = Group)) + 
  geom_bar(position = "dodge", stat = "identity") +
  ggtitle('Aduo Indel Proportion') +
  theme_classic() +
  facet_wrap(. ~ Group) +
  # geom_text(label=Proportion, angle=90, nudge_y = 8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) 
ggsave(filename = './Aduo_Indel_Proportions.png', width = 6, height = 6)


#### Make Aduo_sequencing_error_data_15_16 tibble ####
source('processSeqContentJSON.R')
Aduo_sequencing_error_data_15_16 <- bind_rows(Aduo_Museum, Aduo_Contemporary, .id = 'Sampling') %>%
  mutate(Sampling = case_when(Sampling == '1' ~ 'Museum',
                              Sampling == '2' ~ 'Contemporary',
                              TRUE ~ 'error')) %>%
  mutate(Barcode = str_c('GG', Barcode, "GG", sep = '')) %>% 
  mutate(`End Cut Site` = paste(BP15,BP16, sep = '')) %>%
  select(-BP1:-BP14) %>%
  mutate(only15Err = if_else(BP15 != "G" & BP16 == "G", 1, 0),
         only16Err = if_else(BP16 != "G" & BP15 == "G", 1, 0),
         Err15 = if_else(BP15 != "G", 1, 0),
         Err16 = if_else(BP16 != "G", 1, 0),
         BothErr = if_else(BP15 != "G" & BP16 != "G", 1, 0),
         NeitherErr = if_else(BP15 == "G" & BP16 == "G", 1, 0)) %>%
  group_by(Sampling, Barcode, fqBase, Indels) %>%
  summarise(only15Err = sum(only15Err),
            only16Err = sum(only16Err),
            Err15 = sum(Err15),
            Err16 = sum(Err16),
            BothErr = sum(BothErr),
            NeitherErr = sum(NeitherErr),
            total = n()) %>%
  mutate(only15ErrRate = only15Err / total,
         only16ErrRate = only16Err / total,
         Err15ErrRate = Err15 / total,
         Err16ErrRate = Err16 / total,
         BothErrRate = BothErr / total,
         NeitherErrRate = NeitherErr/ total,
         ProductErr15Err16ErrRate = Err15ErrRate * Err16ErrRate,
         ProductOnly15Only16ErrRate = only15ErrRate * only16ErrRate) %>%
  separate(col=fqBase,
           into = c(NA,NA,NA,NA,"species", "era", "plate_pool", NA),
           remove=FALSE) %>%
  left_join(seq_cont, by=c("species", "era", "plate_pool")) %>%
  filter(fqBase != "Sub-100000-20190528-PIRE-Adu-A-P1P4-L008")

#### bp 15_16 error rate stats ####

sink('./Aduo_ErrorCalc_Analysis_Output.txt', append=TRUE)
writeLines(c("\n","####################################","\n","> mean Err15ErrRate"))
mean(Aduo_sequencing_error_data_15_16$Err15ErrRate)
writeLines(c("\n","####################################","\n","> mean Err16ErrRate"))
mean(Aduo_sequencing_error_data_15_16$Err16ErrRate)
writeLines(c("\n","####################################","\n","> mean only15ErrRate"))
mean(Aduo_sequencing_error_data_15_16$only15Err)
writeLines(c("\n","####################################","\n","> mean only16ErrRate"))
mean(Aduo_sequencing_error_data_15_16$only16Err)
sink()

# ceb: this equation assumes that an error at pos 15 = error at pos 16.  we know that can't be the case because if and error at 15 leads to a randome base at 16, sometimes 16 will be correct
BP15_16ErrLm <- lm(data = Aduo_sequencing_error_data_15_16, formula = Err16ErrRate ~ only16ErrRate + Err15ErrRate)

sink('./Aduo_ErrorCalc_Analysis_Output.txt', append=TRUE)
writeLines(c("\n","####################################","\n","> lm Summary"))
summary(lm(data = Aduo_sequencing_error_data_15_16, formula = Err16ErrRate ~ only16ErrRate + Err15ErrRate))
sink()

fit <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(x = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         y = Err16ErrRate) %>%
  select(y,x) %>%
  lm(formula = y~x) 
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Museum_0Ind"))
linearHypothesis(fit, "x = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# Aduo_sequencing_error_data_15_16 %>%
#   filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
#   mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
#          #rand_exp = (BothErr/total) + only16ErrRate,
#          ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
#          observed = Err16ErrRate) %>%
#   #select(y,x) %>%
#   ggplot(aes(x=observed, y=nonind_pred)) +
#   geom_point(color = 'green4') +
#   geom_smooth(method="lm", fullrange=TRUE, color = "green4") +
#   geom_point(aes(y=ind_exp, x=Err16ErrRate), color="orange3") +
#   geom_smooth(method = "lm",
#               aes(y=ind_exp, x=Err16ErrRate), 
#               color="orange3", 
#               fullrange=TRUE,
#               se=FALSE) +
#   geom_abline(slope=1, intercept=0, linetype="dashed", size=2) +
#   xlim(0,.05) +
#   ylim(0,0.05)
  
# test of similarity between observations and model of non independence of pos 7 and 8 based on eq 1
fit <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(observed, nonind_pred) %>%
  lm(formula = nonind_pred ~ observed) 
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Contemporary_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# test of similarity between observations and model of independence of pos 7 and 8 based on sum and product rules of probability
fit <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(observed, ind_exp) %>%
  lm(formula = ind_exp ~ observed) 
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Contemporary_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()


# test of similarity between observations and model of non independence of pos 7 and 8 based on eq 1
fit <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(observed, nonind_pred) %>%
  lm(formula = nonind_pred ~ observed) 
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Contemporary_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# test of similarity between observations and model of independence of pos 7 and 8 based on sum and product rules of probability
fit <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(observed, ind_exp) %>%
  lm(formula = ind_exp ~observed ) 
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Contemporary_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# There is a difference in the error rates between Museum and Contemporary. We were interested in the behavior of the restriction enzyme so we modeled the effect of having a mistake at the second the last position causing a random base to be insertted in the last position is under-reporting the number of errors which means it 
# The enzyme is reacting differently when interacting with Museum DNA compared with Contemporary DNA at the restriction cutsite. The Error rates are not independent in Museum, they are more independent in Contemporary samples.

fit <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(x = BothErrRate + sqrt(BothErrRate),
         y = Err16ErrRate) %>%
  select(y,x) %>%
  lm(formula = y~x) 
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_0Ind"))
linearHypothesis(fit, "x = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind") %>%
  mutate(Neither = total - Err15 - Err16) %>%
  ungroup() %>%
  group_by(Sampling) %>% 
  summarise(mean_Neither = mean(Neither),
            mean_only15Err = mean(only15Err),
            mean_only16Err = mean(only16Err),
            mean_BothErr = mean(BothErr),
            mean_Err15ErrRate = mean(Err15ErrRate),
            mean_Err16ErrRate = mean(Err16ErrRate)) %>% view()

#### bp 15-16 error rate plots ####
ggplot(Aduo_sequencing_error_data_15_16[Aduo_sequencing_error_data_15_16$Indels == "0Ind",], 
       aes(x = Err15ErrRate * mean_non_g_prop + only16ErrRate, y = Err16ErrRate, color = Sampling)) +
  geom_point() +
  ggtitle('Aduo 0Ind Err15Rate*PropNotG By Only16Err') +
  theme_classic() +
  geom_smooth(method = 'lm',fill = "NA") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme(legend.title = element_blank())


ggsave(filename = './Aduo_0Ind_Err15Rate_PropNotG_by_Only16Err_Linear.png', width = 6, height = 6)

# Museum
# Aduo_Museum15_16Err <- 
#   ggplot(Aduo_sequencing_error_data_15_16[Aduo_sequencing_error_data_15_16$Indels == "0Ind" & Aduo_sequencing_error_data_15_16$Sampling == "Museum",], 
#                             aes(x = Err15ErrRate * mean_non_g_prop + only16ErrRate, y = Err16ErrRate)) +
#   geom_point() +
#   ggtitle('Aduo Museum') +
#   theme_classic() +
#   geom_smooth(method = 'lm', fill = "NA") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   #facet_grid(Sampling ~ ., scales = "free") +
#   coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
#   theme(title = element_blank(), plot.title = element_text(size=12), 
#         legend.title = element_blank(), strip.text = element_blank())

Aduo_Museum15_16Err <- 
Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         #rand_exp = (BothErr/total) + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  #select(y,x) %>%
  ggplot(aes(x=observed, y=nonind_pred)) +
  geom_point(color = 'green4') +
  geom_smooth(method="lm", fullrange=TRUE, color = "green4", se=FALSE) +
  geom_point(aes(y=ind_exp, x=Err16ErrRate), color="orange3") +
  geom_smooth(method = "lm",
              aes(y=ind_exp, x=Err16ErrRate), 
              color="orange3", 
              fullrange=TRUE,
              se=FALSE) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(0,NA) +
  ylim(0,NA) +
  ggtitle('Aduo Museum') +
  theme_classic() +
  theme(title = element_blank(), plot.title = element_text(size=12), 
        legend.title = element_blank(), strip.text = element_blank())


# Contemporary
# Aduo_Contemporary15_16Err <- ggplot(Aduo_sequencing_error_data_15_16[Aduo_sequencing_error_data_15_16$Indels == "0Ind" & Aduo_sequencing_error_data_15_16$Sampling == "Contemporary",], 
#                           aes(x = Err15ErrRate * mean_non_g_prop + only16ErrRate, y = Err16ErrRate)) +
#   geom_point() +
#   ggtitle('Aduo Contemporary') +
#   theme_classic() +
#   geom_smooth(method='lm',fill="NA") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_grid(Sampling ~ ., scales = "free") +
#   coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
#   theme(title = element_blank(), plot.title = element_text(size=12),
#         legend.title = element_blank(), strip.text = element_blank())

Aduo_Contemporary15_16Err <-
Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         #rand_exp = (BothErr/total) + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  #select(y,x) %>%
  ggplot(aes(x=observed, y=nonind_pred)) +
  geom_point(color = 'green4') +
  geom_smooth(method="lm", fullrange=TRUE, color = "green4", se=FALSE) +
  geom_point(aes(y=ind_exp, x=Err16ErrRate), color="orange3") +
  geom_smooth(method = "lm",
              aes(y=ind_exp, x=Err16ErrRate), 
              color="orange3", 
              fullrange=TRUE,
              se=FALSE) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(0,NA) +
  ylim(0,NA) +
  ggtitle('Aduo Contemporary') +
  theme_classic() +
  theme(title = element_blank(), plot.title = element_text(size=12), 
        legend.title = element_blank(), strip.text = element_blank())


png(file = './Aduo_bp15-16_Lm.png', width = 6, height = 6, units = "in", res = 500)
grid.arrange(Aduo_Museum15_16Err, Aduo_Contemporary15_16Err, nrow = 1,
             top = 'Aduo bp 15-16 Lm', left = "Greater Error Rate Than Expected By Random Chance", 
             bottom = "Lesser Error Rate Than Expected By Random Chance")
dev.off()

ggplot(data = Aduo_sequencing_error_data, aes(y = error_rate, fill = Sampling)) + 
  geom_boxplot() +
  ggtitle('Aduo Error Rate') +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  facet_grid(IndelCategory ~ as.factor(base_position))
ggsave(filename = './Aduo_ErrorRateAll.png', width = 6, height = 6)

Aduo_BP_Error_Rate <- Aduo_sequencing_error_data %>% 
  group_by(base_position, IndelCategory) %>%
  summarise(mean(error_rate)) 

sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n", "####################################", "\n", "> BP_Error_Rate"))
Aduo_BP_Error_Rate %>%
  print(n = 1e3)
sink()

Aduo_Museum15MeanError <- aggregate(Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 15 & Aduo_sequencing_error_data$Sampling == "Museum",]$error_rate ~ Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 15 & Aduo_sequencing_error_data$Sampling == "Museum",]$IndelCategory, FUN = mean)
Aduo_Museum16MeanError <- aggregate(Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 16 & Aduo_sequencing_error_data$Sampling == "Museum",]$error_rate ~ Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 16 & Aduo_sequencing_error_data$Sampling == "Museum",]$IndelCategory, FUN = mean)
Aduo_Contemporary15MeanError <- aggregate(Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 15 & Aduo_sequencing_error_data$Sampling == "Contemporary",]$error_rate ~ Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 15 & Aduo_sequencing_error_data$Sampling == "Contemporary",]$IndelCategory, FUN = mean)
Aduo_Contemporary16MeanError <- aggregate(Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 16 & Aduo_sequencing_error_data$Sampling == "Contemporary",]$error_rate ~ Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position == 16 & Aduo_sequencing_error_data$Sampling == "Contemporary",]$IndelCategory, FUN = mean)

Aduo_Error_Rate_15_16 <- ggplot(data = Aduo_sequencing_error_data[Aduo_sequencing_error_data$base_position >= 15 & Aduo_sequencing_error_data$IndelCategory == "NoIndels",],
       aes(y = error_rate, fill = Sampling)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys") +
  ggtitle('Aduo') +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  facet_grid(IndelCategory ~ as.factor(base_position)) +
  theme(legend.position = "none",
  axis.title.y = element_blank(),
  strip.text.y = element_blank())
# ggsave(filename = './Aduo_ErrorRate15_16.png', width = 6, height = 6)

#####################################################################################################################
#### Make Sspi_sequencing_error_data tibble ####
sink('./Sspi_ErrorCalc_Analysis_Output.txt')
writeLines("> Species = Sspi")
writeLines(c("\n","####################################","\n"))
sink()

Sspi_Museum <- read_csv(Sspi_MuseInFile, col_types = str_c(rep('c', 20), collapse = '')) 
Sspi_Contemporary <- read_csv(Sspi_ContInFile, col_types = str_c(rep('c', 20), collapse = ''))

Sspi_Both <- bind_rows(Sspi_Museum, Sspi_Contemporary, .id = 'Sampling')

a<-1; b<-1

Sspi_sequencing_error_data <- bind_rows(Sspi_Museum, Sspi_Contemporary, .id = 'Sampling') %>%
  mutate(Sampling = case_when(Sampling == '1' ~ 'Museum',
                              Sampling == '2' ~ 'Contemporary',
                              TRUE ~ 'error')) %>%
  mutate(Barcode = str_c('GG', Barcode, "GG", sep = '')) %>% 
  gather(position, base, -Sampling:-Indels)  %>% 
  mutate(groupings = case_when(position == 'BP1' | position == 'BP2' ~ 'GG',
                               position == 'BP3' | position == 'BP4' | position == 'BP5' | position == 'BP6' | position == 'BP7' | position == 'BP8' | position == 'BP9' | position == 'BP10' ~ 'Barcode',
                               position == 'BP11' | position == 'BP12' | position == 'BP13' | position == 'BP14' ~ 'Ligation Site',
                               position == 'BP15' | position == 'BP16' ~ 'End Cut Site',
                               TRUE ~ 'error')) %>%
  mutate(base_position = str_extract(position, "[0-9]+") %>% as.integer,
         true_base = str_sub(Barcode, base_position, base_position),
         true_base = if_else(true_base == "", NA_character_, true_base)) %>%
  group_by(Sampling, groupings, base_position, Barcode, fqBase, Read, Indels) %>%
  summarise(correct = sum(base == true_base),
            error = sum(base != true_base),
            total = n()) %>%
  ungroup %>%
  mutate_if(is.character, as.factor) %>%
  mutate(base_position = as.factor(base_position)) %>%
  mutate(base_position = factor(base_position, levels = 1:16, ordered = TRUE),
         Indels = factor(Indels, levels = c("2del","1del","0Ind","1Ins","2Ins","3Ins","4Ins","5Ins","6Ins","7Ins","8Ins"), ordered = TRUE),
         IndelCategory = factor(Indels, levels = c("Insertions","NoIndels","Deletions"), ordered = TRUE),
         groupings = factor(groupings, levels = c("Barcode","Ligation Site","End Cut Site"), ordered = TRUE)) %>%
  filter(groupings != "GG") 

Sspi_sequencing_error_data %>% mutate(error_rate = error/total) -> Sspi_sequencing_error_data
Sspi_sequencing_error_data %>% mutate(success_rate = correct/total) -> Sspi_sequencing_error_data
Sspi_sequencing_error_data$Indels <- factor(Sspi_sequencing_error_data$Indels, levels = c("2del", "1del", "0Ind", "1Ins", "2Ins", "3Ins", "4Ins", "5Ins", "6Ins", "7Ins", "8Ins"))

Sspi_sequencing_error_data$IndelCategory <- gsub(".Ins", "Insertions", Sspi_sequencing_error_data$Indels)
Sspi_sequencing_error_data$IndelCategory <- gsub(".del", "Deletions", Sspi_sequencing_error_data$IndelCategory)
Sspi_sequencing_error_data$IndelCategory <- gsub("0Ind", "NoIndels", Sspi_sequencing_error_data$IndelCategory)
Sspi_sequencing_error_data$IndelCategory <- factor(Sspi_sequencing_error_data$IndelCategory, levels = c("Insertions","NoIndels","Deletions"))

#### Plot mean error rate using Sspi_sequencing_error_data ####
# summarized data to make data labels
Sspi_SummaryData <- Sspi_sequencing_error_data %>% 
  group_by(Sampling, groupings, fqBase, IndelCategory) %>%
  summarise(correct = sum(correct),
            error = sum(error),
            total = sum(total)) %>% 
  mutate(error_rate = error/total) %>% 
  group_by(Sampling, groupings, IndelCategory) %>%
  summarize(mean_error_rate = mean(error_rate)) 

Sspi_SummaryData$Rounded_mean_error_rate <- sprintf("%0.3f", round(Sspi_SummaryData$mean_error_rate, 3))

ggplot(Sspi_SummaryData, aes(x = groupings, y = mean_error_rate, fill = Sampling)) +
  geom_bar(stat = "identity") +
  ggtitle('Sspi Mean Error Rate') +
  theme_classic() +
  geom_text(label = Sspi_SummaryData$Rounded_mean_error_rate, nudge_y = .15) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0), axis.title.x = element_blank()) +
  facet_grid(IndelCategory ~ Sampling)
ggsave(filename ='./Sspi_MeanErrorRate.png', width = 6, height = 6)

# This version is to make combined plot
Sspi_MeanErrorRate <- ggplot(Sspi_SummaryData, aes(x = groupings, y = mean_error_rate)) +
  geom_bar(stat = "identity", fill="lightgrey") +
  ggtitle('Sspi') +
  theme_classic() +
  geom_text(label = Sspi_SummaryData$Rounded_mean_error_rate, size = 3, nudge_y = .15) +
  theme(title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12)) +
  facet_grid(IndelCategory ~ Sampling)


#### Make Collapsed Sspi_sequencing_error_data_2 ####
Sspi_sequencing_error_data_2 <- Sspi_sequencing_error_data %>%
  ungroup() %>%
  group_by(Sampling, groupings, fqBase, IndelCategory, Barcode, Read) %>%
  summarize(correct = sum(correct),
            error = sum(error),
            total = sum(total),
            correct_2 = sum(correct) + 1,
            error_2 = sum(error) + 1,
            total_2 = sum(total) + 2) %>%
  mutate(error_rate = error / total,
         success_rate = correct / total,
         error_rate_2 = error_2 / total_2,
         success_rate_2 = correct_2 / total_2)

# Sspi_error_model_2 <- mixed(cbind(correct_2, error_2) ~ Sampling * groupings * IndelCategory
#                             + (1 | fqBase) + (1 | Barcode) + (1 | Read),
#                             data = Sspi_sequencing_error_data_2,
#                             family = binomial,
#                             method = 'LRT',
#                             control = glmerControl(optimizer = "bobyqa",
#                                                    optCtrl = list(maxfun = 2e8)))

Sspi_error_model_2 <- mixed(cbind(correct_2, error_2) ~ Sampling * groupings * IndelCategory + (1 | fqBase) + (1 | Barcode) + (1 | Read),
                            data = Sspi_sequencing_error_data_2,
                            family = binomial,
                            method = 'LRT',
                            control = glmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 2e9),
                                                   check.conv.grad = .makeCC("warning", tol = 3e-3, relTol = NULL)))


#### Make Collapsed Sspi_sequencing_error_data_3 ####
# collapsing random factors that we don't care about and should have no effect on the results
# one warning is produced by this model but CEB and JS determined it was ok
Sspi_sequencing_error_data_3 <- Sspi_sequencing_error_data %>%
  separate(fqBase, remove = FALSE, sep = "-", into = c(NA, NA, NA, NA, NA, NA, "plate_pool", NA)) %>%
  mutate(library = str_c(plate_pool, Barcode, sep = "_")) %>%
  ungroup() %>%
  group_by(Sampling, groupings, library, IndelCategory, Read) %>%
  summarize(correct = sum(correct),
            error = sum(error),
            total = sum(total),
            correct_2 = sum(correct) + 1,
            error_2 = sum(error) + 1,
            total_2 = sum(total) + 2) %>%
  mutate(error_rate = error / total,
         success_rate = correct / total,
         error_rate_2 = error_2 / total_2,
         success_rate_2 = correct_2 / total_2) 

#### Make Sspi_error_model_3 ####
Sspi_error_model_3 <- mixed(cbind(correct_2, error_2) ~ Sampling * groupings * IndelCategory 
                            + (1 | library) + (1 | Read),
                            data = Sspi_sequencing_error_data_3, 
                            family = binomial, 
                            method = 'LRT',
                            control = glmerControl(optimizer = "bobyqa", 
                                                   optCtrl = list(maxfun = 2e8),
                                                   check.conv.grad = .makeCC("warning", tol = 3e-3, relTol = NULL)))
Sspi_error_model_2 <- Sspi_error_model_3

# estimated means from model & test pairwaise comps
Sspi_error_emm <- emmeans(Sspi_error_model_2, ~ Sampling * groupings * IndelCategory)

# using BH becuase it does the best job of estimating false discovery rate.
Sspi_pairwise_postHoc <- contrast(Sspi_error_emm, 
                                  method = 'pairwise', 
                                  simple = 'each', 
                                  combine = FALSE, 
                                  adjust = "bonferroni")

# write statistical result table to output
sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n","####################################","\n","\n","> Sspi_error_model"))
Sspi_error_model_2
writeLines(c("\n","####################################","\n","\n","> Sspi_pairwise_postHoc"))
Sspi_pairwise_postHoc
sink()

# format model results
Sspi_tbl_emm <- emmeans(Sspi_error_model_2, 
                        ~ Sampling * groupings * IndelCategory, 
                        type = 'response') %>%
  as_tibble %>%
  mutate(IndelCategory = factor(IndelCategory, levels = c('Insertions', 
                                                          'NoIndels', 
                                                          'Deletions')),
         groupings = factor(groupings, levels = c("GG",
                                                  "Barcode",
                                                  "Ligation Site",
                                                  "End Cut Site"), ordered = TRUE))

# plot model results
Sspi_p_emm <- Sspi_tbl_emm %>%
  ggplot(aes(x = groupings, y = prob, ymin = asymp.LCL, ymax = asymp.UCL,
             fill = Sampling)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.2) +
  geom_col(position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ IndelCategory)
Sspi_p_emm
ggsave(filename = './Sspi_p_emm.png', width = 6, height = 6)

# format means from data
Sspi_tbl_means <- Sspi_sequencing_error_data_2 %>%
  ungroup(Barcode, fqBase) %>%
  summarize(mean_success_rate = mean(success_rate),
            mean_success_rate_2 = mean(success_rate_2)) %>%
  mutate(groupings = factor(groupings, levels = c("GG",
                                                  "Barcode",
                                                  "Ligation Site",
                                                  "End Cut Site"), ordered = TRUE))

# plot model results
Sspi_p_means <- Sspi_tbl_means %>%
  ggplot(aes(x = groupings, y = mean_success_rate, fill = Sampling)) +
  geom_col(position = position_dodge(0.9)) +
  geom_point(aes(x = groupings, y = mean_success_rate_2),
             position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ IndelCategory)
Sspi_p_means
ggsave(filename ='./Sspi_p_means.png', width = 6, height = 6)

# plot estimated marginal means from model and means from data
png(file = './Sspi_p_emm_and_Sspi_p_means.png', width = 6, height = 6, units = "in", res = 500)
grid.arrange(Sspi_p_emm, Sspi_p_means, nrow=2,
             top = 'Sspi_p_emm_and_Sspi_p_means')
dev.off()

# combine and compare means to estimated marginal means
Sspi_tbl_emm_means <- Sspi_tbl_emm %>%
  full_join(Sspi_tbl_means) %>%
  rename(emm = prob,
         means = mean_success_rate)

# plot comps between means and emm
Sspi_means_VS_emm <- Sspi_tbl_emm_means %>%
  ggplot(aes(x = means, y = emm)) +
  geom_point() +
  geom_abline() +
  geom_smooth(color = "blue")

Sspi_emm_VS_residuals <- Sspi_tbl_emm_means %>%
  mutate(residuals = emm - means) %>%
  ggplot(aes(x = emm, y = residuals)) +
  geom_point() +
  geom_abline() +
  geom_smooth(color = "blue")


png(file = './Sspi_emm_and_means.png', width = 6, height = 6, units = "in", res = 500)
grid.arrange(Sspi_means_VS_emm, Sspi_emm_VS_residuals, nrow = 2,
             top = 'Sspi_emm_and_means')
dev.off()

# this is the plot of emm with the raw mean prop successful and mean prop successful 2 (which has 1 added to the correct and error values)
Sspi_tbl_emm_means %>%
  ggplot(aes(x = groupings, y = emm, ymin = asymp.LCL, ymax = asymp.UCL,
             fill = Sampling)) +
  ggtitle('Sspi emm mean prop successful 1&2') +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.2, color = "grey") +
  geom_point(aes(x = groupings, y = mean_success_rate_2),
             position = position_dodge(0.9),
             shape = 1,
             size = 10) +
  geom_point(aes(x = groupings, y = means),
             position = position_dodge(0.9),
             color = "green",
             shape = 1,
             size = 10) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.title = element_blank()) +
  facet_wrap(. ~ IndelCategory)
ggsave(filename = './Sspi_emm_mean-prop-successful_and_mean-prop-successful-2.png', width = 8, height = 6)

Sspi_error_model <- Sspi_error_model_2

sum_Sspi_error_emm <- Sspi_error_emm %>%
  summary(type = 'response') %>%
  as_tibble %>%
  ggplot(aes(x = groupings, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, colour = Sampling)) +
  geom_linerange(position = position_dodge(0.5), size = 1) +
  geom_point(position = position_dodge(0.5)) +
  labs(title = 'Sspi Cutsite Proportion Correct by Group', x = "DNA Identity", y = "Proportion Correct") +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  theme_classic() +
  facet_grid(IndelCategory ~ .) +
  theme(legend.title = element_blank())

sum_Sspi_error_emm
ggsave(filename= './Sspi_CutsiteGroupProportionCorrect.png', width = 6, height = 6)

#### Indel Proportion Plots ####
Sspi_IndelCounts <- table(Sspi_Both$fqBase, Sspi_Both$Indels)

Sspi_IndelCountsTidy <- rbind(data.frame(table(Sspi_Museum$fqBase, Sspi_Museum$Indels)), data.frame(table(Sspi_Contemporary$fqBase, Sspi_Contemporary$Indels)))
colnames(Sspi_IndelCountsTidy) <- c("Fqbase", "Indel", "Freq")

sink('./Sspi_ErrorCalc_Analysis_Output.txt', append=TRUE)
writeLines(c("\n","####################################","\n","> Sspi_IndelCounts"))
Sspi_IndelCounts
sink()

Indels <- data.frame(rep(c("2del", "1del", "0Ind", "1Ins", "2Ins", "3Ins", "4Ins", "5Ins", "6Ins", "7Ins", "8Ins"), 2))

Sspi_Museum_total <- count(Sspi_Museum)
Sspi_Contemporary_total <- count(Sspi_Contemporary)
Sspi_Both_total <- count(Sspi_Both)

Proportion <- as.numeric(paste(data.frame(c(round((count(Sspi_Museum[Sspi_Museum$Indels=="2del",])/Sspi_Museum_total*100),3),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="1del",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="0Ind",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="1Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="2Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="3Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="4Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="5Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="6Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="7Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Museum[Sspi_Museum$Indels=="8Ins",])/Sspi_Museum_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="2del",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="1del",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="0Ind",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="1Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="2Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="3Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="4Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="5Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="6Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="7Ins",])/Sspi_Contemporary_total*100),3)),
                                            (round((count(Sspi_Contemporary[Sspi_Contemporary$Indels=="8Ins",])/Sspi_Contemporary_total*100),3)))),sep = "\n"))

Group <- data.frame(c(rep("Museum", 11), rep("Contemporary", 11)))

Sspi_Indel_Proportions_Table <- data.frame(Indels, Proportion, Group)
colnames(Sspi_Indel_Proportions_Table) <- c("Indels", "Prop", "Group")

Sspi_Indel_Proportions_Table$Indels <- factor(Sspi_Indel_Proportions_Table$Indels, ordered = TRUE, levels = c("2del", "1del", "0Ind", "1Ins", "2Ins" ,"3Ins" ,"4Ins" ,"5Ins" ,"6Ins" ,"7Ins" ,"8Ins"))

sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n","####################################","\n","> Sspi_Indel_Proportions_Table"))
Sspi_Indel_Proportions_Table
sink()

ggplot(Sspi_Indel_Proportions_Table, aes(x = Indels, y = Proportion, fill = Group)) + 
  geom_bar(position = "dodge", stat = "identity") +
  ggtitle('Sspi Indel Proportion') +
  theme_classic() +
  facet_wrap(. ~ Group) +
  # geom_text(label=Proportion, angle=90, nudge_y = 8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) 
ggsave(filename = './Sspi_Indel_Proportions.png', width = 6, height = 6)


#### Make Sspi_sequencing_error_data_15_16 tibble ####
source('processSeqContentJSON.R')
Sspi_sequencing_error_data_15_16 <- bind_rows(Sspi_Museum, Sspi_Contemporary, .id = 'Sampling') %>%
  mutate(Sampling = case_when(Sampling == '1' ~ 'Museum',
                              Sampling == '2' ~ 'Contemporary',
                              TRUE ~ 'error')) %>%
  mutate(Barcode = str_c('GG', Barcode, "GG", sep = ''),
         `End Cut Site` = str_c(BP15,BP16, sep = '')) %>%
  select(-BP1:-BP14) %>%
  mutate(only15Err = if_else(BP15 != "G" & BP16 == "G", 1, 0),
         only16Err = if_else(BP16 != "G" & BP15 == "G", 1, 0),
         Err15 = if_else(BP15 != "G", 1, 0),
         Err16 = if_else(BP16 != "G", 1, 0),
         BothErr = if_else(BP15 != "G" & BP16 != "G", 1, 0),
         NeitherErr = if_else(BP15 == "G" & BP16 == "G", 1, 0)) %>%
  group_by(Sampling, Barcode, fqBase, Indels) %>%
  summarise(only15Err = sum(only15Err),
            only16Err = sum(only16Err),
            Err15 = sum(Err15),
            Err16 = sum(Err16),
            BothErr = sum(BothErr),
            NeitherErr = sum(NeitherErr),
            total = n()) %>%
  mutate(only15ErrRate = only15Err / total,
         only16ErrRate = only16Err / total,
         Err15ErrRate = Err15 / total,
         Err16ErrRate = Err16 / total,
         BothErrRate = BothErr / total,
         NeitherErrRate = NeitherErr/ total,
         ProductErr15Err16ErrRate = Err15ErrRate * Err16ErrRate,
         ProductOnly15Only16ErrRate = only15ErrRate * only16ErrRate) %>%
  separate(col=fqBase,
           into = c(NA,NA,NA,NA,"species", "era", "plate_pool", NA),
           remove=FALSE) %>%
  left_join(seq_cont, by=c("species", "era", "plate_pool"))

#### bp 15_16 error rate stats ####

sink('./Sspi_ErrorCalc_Analysis_Output.txt', append=TRUE) 
writeLines(c("\n","####################################","\n","> mean Err15ErrRate"))
mean(Sspi_sequencing_error_data_15_16$Err15ErrRate)
writeLines(c("\n","####################################","\n","> mean Err16ErrRate"))
mean(Sspi_sequencing_error_data_15_16$Err16ErrRate)
writeLines(c("\n","####################################","\n","> mean only15ErrRate"))
mean(Sspi_sequencing_error_data_15_16$only15Err)
writeLines(c("\n","####################################","\n","> mean only16ErrRate"))
mean(Sspi_sequencing_error_data_15_16$only16Err)
sink()

BP15_16ErrLm <- lm(data = Sspi_sequencing_error_data_15_16, formula = Err16ErrRate ~ only16ErrRate + Err15ErrRate)

sink('./Sspi_ErrorCalc_Analysis_Output.txt', append=TRUE)
writeLines(c("\n","####################################","\n","> lm Summary"))
summary(lm(data = Sspi_sequencing_error_data_15_16, formula = Err16ErrRate ~ only16ErrRate + Err15ErrRate))
sink()

# test of eq1 (non independence of pos 7 and 8)
fit <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(nonind_pred,observed) %>%
  lm(formula = nonind_pred ~ observed) 
sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Museum_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# test of null model (independence of pos 7 and 8)
fit <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(ind_exp,observed) %>%
  lm(formula = ind_exp~observed) 
sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Museum_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# test of eq1 (non independence of pos 7 and 8)
fit <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(nonind_pred,observed) %>%
  lm(formula = nonind_pred~observed) 
sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Museum_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

# test of null model (independence of pos 7 and 8)
fit <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  select(ind_exp,observed) %>%
  lm(formula = ind_exp~observed) 
sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
AIC(fit)
BIC(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_Museum_0Ind"))
linearHypothesis(fit, "observed = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()


# There is a difference in the error rates between Museum and Contemporary. We were interested in the behavior of the restriction enzyme so we modeled the effect of having a mistake at the second the last position causing a random base to be insertted in the last position is under-reporting the number of errors which means it 
# The enzyme is reacting differently when interacting with Museum DNA compared with Contemporary DNA at the restriction cutsite. The Error rates are not independent in Museum, they are more independent in Contemporary samples.

fit <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(x = BothErrRate + sqrt(BothErrRate),
         y = Err16ErrRate) %>%
  select(y,x) %>%
  lm(formula = y~x) 
sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
summary(fit)
writeLines(c("\n","####################################","\n","> linearHypothesis_0Ind"))
linearHypothesis(fit, "x = 1")
linearHypothesis(fit, c("(Intercept) = 0"))
sink()

Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind") %>%
  mutate(Neither = total - Err15 - Err16) %>%
  ungroup() %>%
  group_by(Sampling) %>% 
  summarise(mean_Neither = mean(Neither),
            mean_only15Err = mean(only15Err),
            mean_only16Err = mean(only16Err),
            mean_BothErr = mean(BothErr),
            mean_Err15ErrRate = mean(Err15ErrRate),
            mean_Err16ErrRate = mean(Err16ErrRate)) 

#### bp 15-16 error rate plots ####
ggplot(Sspi_sequencing_error_data_15_16[Sspi_sequencing_error_data_15_16$Indels == "0Ind",], 
       aes(x = Err15ErrRate * mean_non_g_prop + only16ErrRate, y = Err16ErrRate, color = Sampling)) +
  geom_point() +
  ggtitle('Sspi 0Ind Err15Rate*PropNotG By Only16Err') +
  theme_classic() +
  geom_smooth(method = 'lm',fill = "NA") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme(legend.title = element_blank())
ggsave(filename = './Sspi_0Ind_Err15Rate_PropNotG_by_Only16Err_Linear.png', width = 6, height = 6)

# Museum
# Sspi_Museum15_16Err <- ggplot(Sspi_sequencing_error_data_15_16[Sspi_sequencing_error_data_15_16$Indels == "0Ind" & Sspi_sequencing_error_data_15_16$Sampling == "Museum",], 
#                                  aes(x = Err15ErrRate * mean_non_g_prop + only16ErrRate, y = Err16ErrRate)) +
#   geom_point() +
#   ggtitle('Sspi Museum') +
#   theme_classic() +
#   geom_smooth(method = 'lm', fill = "NA") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_grid(Sampling ~ ., scales = "free") +
#   coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
#   theme(title = element_blank(), plot.title = element_text(size=12), 
#         legend.title = element_blank(), strip.text = element_blank())

Sspi_Museum15_16Err <- 
  Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Museum") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         #rand_exp = (BothErr/total) + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  #select(y,x) %>%
  ggplot(aes(x=observed, y=nonind_pred)) +
  geom_point(color = 'green4') +
  geom_smooth(method="lm", fullrange=TRUE, color = "green4", se=FALSE) +
  geom_point(aes(y=ind_exp, x=Err16ErrRate), color="orange3") +
  geom_smooth(method = "lm",
              aes(y=ind_exp, x=Err16ErrRate), 
              color="orange3", 
              fullrange=TRUE,
              se=FALSE) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(0,NA) +
  ylim(0,NA) +
  ggtitle('Sspi Museum') +
  theme_classic() +
  theme(title = element_blank(), plot.title = element_text(size=12), 
        legend.title = element_blank(), strip.text = element_blank())

# Contemporary
# Sspi_Contemporary15_16Err <- ggplot(Sspi_sequencing_error_data_15_16[Sspi_sequencing_error_data_15_16$Indels == "0Ind" & Sspi_sequencing_error_data_15_16$Sampling == "Contemporary",], 
#                                aes(x = Err15ErrRate * mean_non_g_prop + only16ErrRate, y = Err16ErrRate)) +
#   geom_point() +
#   ggtitle('Sspi Contemporary') +
#   theme_classic() +
#   geom_smooth(method='lm',fill="NA") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_grid(Sampling ~ ., scales = "free") +
#   coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
#   theme(title = element_blank(), plot.title = element_text(size=12),
#         legend.title = element_blank(), strip.text = element_blank())

Sspi_Contemporary15_16Err <-
  Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind" & Sampling == "Contemporary") %>%
  mutate(nonind_pred = Err15ErrRate * mean_non_g_prop + only16ErrRate,
         #rand_exp = (BothErr/total) + only16ErrRate,
         ind_exp = Err15ErrRate * Err16ErrRate + only16ErrRate,
         observed = Err16ErrRate) %>%
  #select(y,x) %>%
  ggplot(aes(x=observed, y=nonind_pred)) +
  geom_point(color = 'green4') +
  geom_smooth(method="lm", fullrange=TRUE, color = "green4", se=FALSE) +
  geom_point(aes(y=ind_exp, x=Err16ErrRate), color="orange3") +
  geom_smooth(method = "lm",
              aes(y=ind_exp, x=Err16ErrRate), 
              color="orange3", 
              fullrange=TRUE,
              se=FALSE) +
  geom_abline(slope=1, intercept=0, linetype="dashed", size=1) +
  xlim(0,NA) +
  ylim(0,NA) +
  ggtitle('Sspi Contemporary') +
  theme_classic() +
  theme(title = element_blank(), plot.title = element_text(size=12), 
        legend.title = element_blank(), strip.text = element_blank())


png(file = './Sspi_bp15-16_Lm.png', width = 6, height = 6, units = "in", res = 500)
grid.arrange(Sspi_Museum15_16Err, Sspi_Contemporary15_16Err, nrow = 1,
             top = 'Sspi bp 15-16 Lm', left = "Greater Error Rate Than Expected By Random Chance", 
             bottom = "Lesser Error Rate Than Expected By Random Chance")
dev.off()

ggplot(data = Sspi_sequencing_error_data, aes(y = error_rate, fill = Sampling)) + 
  geom_boxplot() +
  ggtitle('Sspi Error Rate') +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  facet_grid(IndelCategory ~ as.factor(base_position))
ggsave(filename = './Sspi_ErrorRateAll.png', width = 6, height = 6)

Sspi_BP_Error_Rate <- Sspi_sequencing_error_data %>% 
  group_by(base_position, IndelCategory) %>%
  summarise(mean(error_rate)) 

sink('./Sspi_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n", "####################################", "\n", "> BP_Error_Rate"))
Sspi_BP_Error_Rate %>%
  print(n = 1e3)
sink()

Sspi_Museum15MeanError <- aggregate(Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 15 & Sspi_sequencing_error_data$Sampling == "Museum",]$error_rate ~ Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 15 & Sspi_sequencing_error_data$Sampling == "Museum",]$IndelCategory, FUN = mean)
Sspi_Museum16MeanError <- aggregate(Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 16 & Sspi_sequencing_error_data$Sampling == "Museum",]$error_rate ~ Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 16 & Sspi_sequencing_error_data$Sampling == "Museum",]$IndelCategory, FUN = mean)
Sspi_Contemporary15MeanError <- aggregate(Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 15 & Sspi_sequencing_error_data$Sampling == "Contemporary",]$error_rate ~ Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 15 & Sspi_sequencing_error_data$Sampling == "Contemporary",]$IndelCategory, FUN = mean)
Sspi_Contemporary16MeanError <- aggregate(Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 16 & Sspi_sequencing_error_data$Sampling == "Contemporary",]$error_rate ~ Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position == 16 & Sspi_sequencing_error_data$Sampling == "Contemporary",]$IndelCategory, FUN = mean)

Sspi_Error_Rate_15_16 <- ggplot(data = Sspi_sequencing_error_data[Sspi_sequencing_error_data$base_position >= 15 & Sspi_sequencing_error_data$IndelCategory == "NoIndels",],
                                aes(y = error_rate, fill = Sampling)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys") +
  ggtitle('Sspi') +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  facet_grid(IndelCategory ~ as.factor(base_position)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text.y = element_blank())
# ggsave(filename = './Sspi_ErrorRate15_16.png', width = 6, height = 6)


######################################################################################################################
#### Combined Species Plots ####

png(file = './Both_MeanErrorRates.png', height = 6.5, width = 10, units = "in", res = 500)
grid.arrange(Aduo_MeanErrorRate, Sspi_MeanErrorRate, nrow = 2)
dev.off()

png(file = './Both_MeanErrorRates15_16.png', height = 6.5, width = 10, units = "in", res = 500)
grid.arrange(Aduo_Error_Rate_15_16, Sspi_Error_Rate_15_16, nrow = 2)
dev.off()

png(file = './Both_bp15-16_Lm.png', height = 6.5, width = 10, units = "in", res = 500)
grid.arrange(Aduo_Museum15_16Err, Aduo_Contemporary15_16Err, Sspi_Museum15_16Err, Sspi_Contemporary15_16Err, nrow = 2,
             left = "Greater Error Rate Than Expected By Random Chance", 
             bottom = "Lesser Error Rate Than Expected By Random Chance")
dev.off()

#### Global Chi-Squared Test ####

Aduo_ChiSQ <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind") %>%
  mutate(Neither = total - Err15 - Err16) %>%
  group_by(Sampling) %>%
  summarise(Observed_NeitherErr = sum(Neither),
            Observed_only15 = sum(only15Err),
            Observed_only16 = sum(only16Err),
            Observed_BothErr = sum(BothErr),
            Expected_NeitherErr = (1 - mean(Err15ErrRate)) * (1 - mean(Err16ErrRate)) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
            Expected_only16 = (1 - mean(Err16ErrRate)) * mean(Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
            Expected_only15 = mean(Err15ErrRate) * (1 - mean(Err16ErrRate)) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
            Expected_BothErr = mean(Err15ErrRate) * mean(Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr))

Sspi_ChiSQ <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind") %>%
  mutate(Neither = total - Err15 - Err16) %>%
  group_by(Sampling) %>%
  summarise(Observed_NeitherErr = sum(Neither),
            Observed_only15 = sum(only15Err),
            Observed_only16 = sum(only16Err),
            Observed_BothErr = sum(BothErr),
            Expected_NeitherErr = (1 - mean(Err15ErrRate)) * (1 - mean(Err16ErrRate)) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
            Expected_only16 = (1 - mean(Err16ErrRate)) * mean(Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
            Expected_only15 = mean(Err15ErrRate) * (1 - mean(Err16ErrRate)) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
            Expected_BothErr = mean(Err15ErrRate) * mean(Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr))

Chi_Table <- bind_rows(Aduo_ChiSQ, Sspi_ChiSQ, .id = 'Species')  %>%
  mutate(Species = case_when(Species == '1' ~ 'Aduo', 
                             Species == '2' ~ 'Sspi')) %>%
  pivot_longer(cols = Observed_NeitherErr:Expected_BothErr, 
               names_to = c("ObsExp", "ErrType"),
               names_pattern = "(.*)_(.*)",
               values_to = "NumReads") %>%
  pivot_wider(names_from = ObsExp, 
              values_from = NumReads) %>%
  mutate(Treatment_Combination = str_c(Species, Sampling, ErrType, sep = "_"),
         Expected_Props = Expected/sum(Expected)) %>%
  group_by(Species)

Pooled_Chi_Table <- Chi_Table %>%
  group_by(ErrType) %>%
  summarise(Observed = sum(Observed), Expected = sum(Expected)) %>%
  mutate(Expected_Props = round(Expected, 0)/sum(round(Expected, 0)))
Pooled_Chi_Test <- chisq_test(x = Pooled_Chi_Table$Observed, p = Pooled_Chi_Table$Expected_Props)
pairwise_chisq_test_against_p(x = Pooled_Chi_Table$Observed, p = Pooled_Chi_Table$Expected_Props, p.adj = "bonferroni")

Aduo_Contemp_Chi_Table <- Chi_Table %>%
  filter(Species == "Aduo", Sampling == "Contemporary") %>%
  mutate(Expected_Props = round(Expected, 0)/sum(round(Expected, 0))) 
Aduo_Contemp_Chi_Test <- chisq_test(x = Aduo_Contemp_Chi_Table$Observed, p = Aduo_Contemp_Chi_Table$Expected_Props)
pairwise_chisq_test_against_p(x = Aduo_Contemp_Chi_Table$Observed, p = Aduo_Contemp_Chi_Table$Expected_Props, p.adj = "bonferroni")

Sspi_Contemp_Chi_Table <- Chi_Table %>%
  filter(Species == "Sspi", Sampling == "Contemporary") %>%
  mutate(Expected_Props = round(Expected, 0)/sum(round(Expected, 0)))
Sspi_Contemp_Chi_Test <- chisq_test(x = Sspi_Contemp_Chi_Table$Observed, p = Sspi_Contemp_Chi_Table$Expected_Props)
pairwise_chisq_test_against_p(x = Sspi_Contemp_Chi_Table$Observed, p = Sspi_Contemp_Chi_Table$Expected_Props, p.adj = "bonferroni")

Aduo_Museum_Chi_Table <- Chi_Table %>%
  filter(Species == "Aduo", Sampling == "Museum") %>%
  mutate(Expected_Props = round(Expected, 0)/sum(round(Expected, 0)))
Aduo_Museum_Chi_Test <- chisq_test(x = Aduo_Museum_Chi_Table$Observed, p = Aduo_Museum_Chi_Table$Expected_Props)
pairwise_chisq_test_against_p(x = Aduo_Museum_Chi_Table$Observed, p = Aduo_Museum_Chi_Table$Expected_Props, p.adj = "bonferroni")

Sspi_Museum_Chi_Table <- Chi_Table %>%
  filter(Species == "Sspi", Sampling == "Museum") %>%
  mutate(Expected_Props = round(Expected, 0)/sum(round(Expected, 0)))
Sspi_Museum_Chi_Test <- chisq_test(x = Sspi_Museum_Chi_Table$Observed, p = Sspi_Museum_Chi_Table$Expected_Props)
pairwise_chisq_test_against_p(x = Sspi_Museum_Chi_Table$Observed, p = Sspi_Museum_Chi_Table$Expected_Props, p.adj = "bonferroni")

Hetero_ChiSQ <- sum(Aduo_Contemp_Chi_Test$statistic,
Sspi_Contemp_Chi_Test$statistic,
Aduo_Museum_Chi_Test$statistic,
Sspi_Museum_Chi_Test$statistic) - Pooled_Chi_Test$statistic 
pchisq(Hetero_ChiSQ, df = 9, lower.tail = FALSE)

Chi_matrix <- as.matrix(bind_rows(Aduo_ChiSQ, Sspi_ChiSQ, .id = 'Species')  %>% 
  mutate(Species = case_when(Species == '1' ~ 'Aduo', 
                             Species == '2' ~ 'Sspi')) %>%
  select(starts_with("Observed")))
chisq_test(Chi_matrix)


#### ChiSq Heterogeneity Test for every combo of Species, Sampling, Barcode, fqBase ####
# We need to make sure that we can combine categories by testing for heterogeneity a result of 0 means we CANNOT combine categories.
# 4 tests rows show be unique by fqbase and barcode 1 test for each combo of species and collection time 
# Chi-Sq 0s replaced with 1s to prevent NaNs

Aduo_ChiSQ <- Aduo_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind") %>%
  rename(Observed_NeitherErr = NeitherErr,
         Observed_only15 = only15Err,
         Observed_only16 = only16Err,
         Observed_BothErr = BothErr) %>%
  mutate(Expected_NeitherErr = (1 - Err15ErrRate) * (1 - Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
         Expected_only16 = (1 - Err16ErrRate) * Err16ErrRate * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
         Expected_only15 = Err15ErrRate * (1 - Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
         Expected_BothErr = Err15ErrRate * Err16ErrRate * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr))

Sspi_ChiSQ <- Sspi_sequencing_error_data_15_16 %>%
  filter(Indels == "0Ind") %>%
  rename(Observed_NeitherErr = NeitherErr,
         Observed_only15 = only15Err,
         Observed_only16 = only16Err,
         Observed_BothErr = BothErr) %>%
  mutate(Expected_NeitherErr = (1 - Err15ErrRate) * (1 - Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
         Expected_only16 = (1 - Err16ErrRate) * Err16ErrRate * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
         Expected_only15 = Err15ErrRate * (1 - Err16ErrRate) * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr),
         Expected_BothErr = Err15ErrRate * Err16ErrRate * (Observed_NeitherErr + Observed_only15 + Observed_only16 + Observed_BothErr))

Chi_Table <- bind_rows(Aduo_ChiSQ, Sspi_ChiSQ, .id = 'Species') %>%
  mutate(Species = case_when(Species == '1' ~ 'Aduo', 
                             Species == '2' ~ 'Sspi'),
         total_observations = sum(Observed_NeitherErr,
                                  Observed_only15,
                                  Observed_only16,
                                  Observed_BothErr)) %>%
  select(Species:Indels,
         total,
         total_observations,
         Err15,
         Err16,
         starts_with("Observed"),
         starts_with("Expected")) %>%
  mutate(ObsRate_only15 = Observed_only15/total_observations,
         ObsRate_only16 = Observed_only16/total_observations,
         ObsRate_BothErr = Observed_BothErr/total_observations,
         ObsRate_NeitherErr = Observed_NeitherErr/total_observations,
         total_proportion = ObsRate_only15 + ObsRate_only16 + ObsRate_BothErr + ObsRate_NeitherErr) #%>%
  # filter(total_observations > 100)

write.csv(Chi_Table, file = "./Chi_Table_ceb.csv", row.names = FALSE)

# view histogram of number of observations 
Chi_Table %>%
  ggplot(aes(x=total_observations)) +
  geom_histogram() +
  facet_wrap(Species ~ Sampling, scales="free")

# Aduo contemp vs hist by error_category visualization
Chi_Table %>%
  ungroup %>%
  filter(Species == "Aduo") %>%
  select(Species:Indels, starts_with("ObsRate")) %>%
  pivot_longer(cols = starts_with("ObsRate"), names_to="error_category", values_to = "proportion") %>%
  ggplot(aes(x=Sampling, y=proportion, fill=Sampling)) +
    geom_boxplot() +
    facet_grid(error_category ~ ., scales = "free")

# Aduo contemp vs hist by error_category visualization
Chi_Table %>%
  ungroup %>%
  # filter(Species == "Aduo") %>%
  select(Species:Indels, starts_with("ObsRate")) %>%
  pivot_longer(cols = starts_with("ObsRate"), names_to="error_category", values_to = "proportion") %>%
  ggplot(aes(x=Barcode, y=proportion, fill=Barcode)) +
  geom_boxplot() +
  facet_grid(error_category ~ Species + Sampling, scales = "free")

Aduo_Contemp_Chi_matrix <- Chi_Table %>%
  ungroup %>%
  filter(Species == "Aduo", Sampling == "Contemporary", Barcode == "GGAAGCCGGTTGCAGG") %>%
  select(starts_with("Observed"))

chisq_test(Aduo_Contemp_Chi_matrix)
chisq.test(Aduo_Contemp_Chi_matrix)
pairwise_chisq_gof_test(Aduo_Contemp_Chi_matrix %>%
                          pivot_longer())
library(rcompanion)
pairwiseNominalIndependence(as.matrix(Aduo_Contemp_Chi_matrix),
                            fisher = FALSE,
                            gtest  = FALSE,
                            chisq  = TRUE,
                            method = "fdr")



Sspi_Contemp_Chi_matrix <- Chi_Table %>%
  ungroup %>%
  filter(Species == "Sspi", Sampling == "Contemporary") %>%
  select(starts_with("Observed"))
chisq_test(Sspi_Contemp_Chi_matrix)

Aduo_Museum_Chi_matrix <- Chi_Table %>%
  ungroup %>%
  filter(Species == "Aduo", Sampling == "Museum") %>%
  select(starts_with("Observed"))
chisq_test(Aduo_Museum_Chi_matrix)

Sspi_Museum_Chi_matrix <- Chi_Table %>%
  ungroup %>%
  filter(Species == "Sspi", Sampling == "Museum") %>%
  select(starts_with("Observed"))
chisq_test(Sspi_Museum_Chi_matrix)

write.csv(Chi_matrix, file = "./Chi_matrix.csv", row.names = FALSE)

write.csv(Aduo_Contemp_Chi_matrix, file = "./Aduo_Contemp_Chi_matrix.csv", row.names = FALSE)
write.csv(Sspi_Contemp_Chi_matrix, file = "./Sspi_Contemp_Chi_matrix.csv", row.names = FALSE)
write.csv(Aduo_Museum_Chi_matrix, file = "./Aduo_Museum_Chi_matrix.csv", row.names = FALSE)
write.csv(Sspi_Museum_Chi_matrix, file = "./Sspi_Museum_Chi_matrix.csv", row.names = FALSE)





