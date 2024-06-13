rm(list = ls())

# Set working directory to location of this script
setwd(dir = "/Users/TEST/Texas A&M University - Corpus Christi/Bird, Chris - MartinFrenchMasters/Final Versions/")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("../Final Versions/")

getwd()

#### Set User Defined Variables ####
DataPath <- "./SequenceCounts/"
OUTDIR <- "./Output"

library(emmeans)
library(tidyverse)
library(RColorBrewer)
library(magrittr)
library(betareg)
library(psych)
library(lmtest)
library(rcompanion)
library(MASS)
library(plyr)
library(car)
library(scales)
library(afex)
library(lme4)
library(optimx)

#### Make tbl_full. Combine all sequence counts files ####
tbl_full <- bind_rows(read_tsv("SequenceCounts/Aduo_SequenceCounts_3Prime.tsv", col_names=FALSE) %>% 
                        dplyr::mutate(Trimming = "Trimmed"), 
                      read_tsv("SequenceCounts/Sspi_SequenceCounts_3PrimeTrim.tsv", col_names=FALSE) %>% 
                        dplyr::mutate(Trimming = "Trimmed"),
                      read_tsv("SequenceCounts/Aduo_SequenceCounts_NotTrimmed.tsv", col_names=FALSE) %>%
                        dplyr::mutate(Trimming = "Untrimmed"),
                      read_tsv("SequenceCounts/Sspi_SequenceCounts_NotTrimmed.tsv", col_names=FALSE) %>%
                        dplyr::mutate(Trimming = "Untrimmed")) %>%
  dplyr::rename(cat = X1, barcode=X2, fqfile=X3, count=X4) %>%
  separate(col=fqfile, 
           into=c("date","project","species", "collectionTime","libID","lane","readDir", "fileType"), 
           remove=FALSE,
           convert=TRUE,
           sep="[-._]") %>%
  dplyr::mutate(collectionTime = case_when(collectionTime == 'A' ~ 'Historical',
                               collectionTime == 'C' ~ 'Contemporary',
                               collectionTime == 'PreFrag' ~ 'Historical',
                               collectionTime == 'NonFrag' ~ 'Contemporary'),
         species = case_when(species == 'Ssp' ~ 'Sspi',
                             species == 'Adu' ~ 'Aduo',
                             species == 'Ola' ~ 'Olat'),
         # barcodeID = case_when(barcode == 'CGATGCTCTGCA' ~ '1st',
         #                       barcode == 'AAGCCGGTTGCA' ~ '2nd'),
         indID = paste(libID, barcode,sep="_"),
         cat = dplyr::recode(cat, 'TotalAfterTrim' = 'Total'),
         cat2 = case_when(cat == '2Del' ~ 'MotifWithIndels',
                          cat == '1Del' ~ 'MotifWithIndels',
                          cat == 'NoInd' ~ 'NoInd',
                          cat == '1Ins' ~ 'MotifWithIndels',
                          cat == '2Ins' ~ 'MotifWithIndels',
                          cat == '3Ins' ~ 'MotifWithIndels',
                          cat == '4Ins' ~ 'MotifWithIndels',
                          cat == '5Ins' ~ 'MotifWithIndels',
                          cat == '6Ins' ~ 'MotifWithIndels',
                          cat == '7Ins' ~ 'MotifWithIndels',
                          cat == '8Ins' ~ 'MotifWithIndels',
                          cat == 'Total' ~ 'Total',
                          cat == 'IndelCat' ~ 'IndelCat')) %>%
  dplyr::select(-fileType,-lane,-project,-date) %>%
  filter(!(cat == "Total" & (readDir == "F" | readDir == "R1"))) %>%
  filter(!(species == 'Aduo' & libID == 'P1P4'))

#### Organize tibble by libID ####
tbl_libID <- 
  tbl_full %>%
  
  dplyr::group_by(cat2,species,collectionTime,libID, Trimming) %>%
  dplyr::summarise(countSum = sum(count)) %>%
  tidyr::pivot_wider(names_from = cat2,
              values_from = countSum,
              values_fill = list(countSum = 0)) %>%
  dplyr::mutate(NoMotif = Total - IndelCat) %>%
  dplyr::select(-Total,-IndelCat) %>%
  tidyr::pivot_longer(cols = MotifWithIndels:NoMotif,
               names_to = 'cat2',
               values_to = 'countSum') %>%
  
  
  tidyr::pivot_wider(names_from = Trimming,
              values_from = countSum) %>%
  ungroup() %>%
  dplyr::group_by(species, collectionTime, libID) %>%
  dplyr::mutate(Trimmed_Filtered = sum(Untrimmed) - sum(Trimmed)) %>%
  tidyr::pivot_wider(names_from = cat2,
              values_from = c(Untrimmed,Trimmed)) %>%
  tidyr::pivot_longer(cols=Trimmed_Filtered:Trimmed_NoMotif,
               names_sep = "_",
               names_to = c("Trimming", "cat2"),
               values_to = "countSum") %>%
  ungroup() %>%
  
  dplyr::group_by(species, collectionTime, libID, Trimming) %>%
  dplyr::mutate(prop = countSum / sum(countSum),
         collectionTime_Trimming = paste(collectionTime, Trimming, sep = "_"),
         cat2 = factor(cat2,c("NoInd", "MotifWithIndels", "NoMotif", "Filtered"))) %>%
  dplyr::arrange(species, collectionTime, libID, Trimming, cat2) %>%
  tidyr::pivot_longer(cols = countSum:prop,
               names_to = "measure",
               values_to = 'value') %>%
  dplyr::mutate(collectionTime_Measure = paste(collectionTime, measure, sep = "_")) %>%
  dplyr::mutate(libID = factor(libID,c("P1P1","P1P2","P1P3","P1P4","P1P5","P1P6","P1P7","P1P8",
                                "P1P9","P1P10","P1P11","P1P12","P1P13","P1P14","P1P15",
                                "P1P16","P1P17","P1P18","P1P19","P1P20","P1P21","P1P22",
                                "P1P23","P1P24","50nMLow","50nMHigh","500nMLow","500nMHigh")),
         collectionTime_Trimming = factor(collectionTime_Trimming,c("Contemporary_Untrimmed","Contemporary_Trimmed","Historical_Untrimmed","Historical_Trimmed")),
         cat = case_when(cat2 == 'MotifWithIndels' ~ 'Motif',
                         cat2 == 'NoInd' ~ 'Motif',
                         cat2 == 'NoMotif' ~ 'NoMotif',
                         cat2 == 'Filtered' ~ 'Filtered'))

#### Create tibble with mean and st.dev ####
tbl_libID_means <- tbl_libID %>%
  dplyr::select(-collectionTime_Measure) %>%
  tidyr::pivot_wider(names_from = measure,
              values_from = c(value)) %>%
  dplyr::group_by(species, collectionTime, Trimming, cat2) %>% 
  dplyr::summarise(tf_Spp_tr = n(),
            mean_countSum = mean(countSum),
            sd_countSum = sd(countSum),
            mean_prop = mean(prop),
            sd_prop = sd(prop)) %>%
  tidyr::pivot_longer(cols = mean_countSum:sd_prop,
               names_pattern = '(.*)_(.*)',
               names_to = c('measure2', 'measure')) %>%
  tidyr::pivot_wider(names_from = measure2,
              values_from = c(value)) %>%
  dplyr::mutate(Trimming = factor(Trimming, c("Untrimmed", "Trimmed")), 
         cat = case_when(cat2 == 'MotifWithIndels' ~ 'Motif',
                         cat2 == 'NoInd' ~ 'Motif',
                         cat2 == 'NoMotif' ~ 'NoMotif',
                         cat2 == 'Filtered' ~ 'Filtered'))


#### Plot Proportions and Counts of Facet grid collectionTime ~ Trimming ####
poolBarPlot_means <- function(data=tbl_libID_means, msr="prop"){
  if(msr == "prop"){
    title = paste("Proportion of Reads", sep = "\n")
  } else if(msr == "countSum"){
    title = paste("Number of Reads", sep = "\n")
  }
  
  data %>%
    filter(measure == msr) %>%
    ggplot(aes(x=species, y=mean, fill=cat2)) +
    geom_bar(stat="identity", width=0.9) +
    scale_fill_brewer(palette="Set1") +
    theme_classic() +
    theme(title = element_text(size = 16), legend.title = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12), strip.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)) +
    labs(title = title) +
    facet_grid(collectionTime ~ Trimming) 
}          


poolBarPlot_means(msr = "prop")
ggsave(filename='Output/Mean Read Proportions - Barplots.png', width = 6, height = 6)

poolBarPlot_means(msr = "countSum")
ggsave(filename='Output/Mean Read Counts - Barplots.png', width = 6, height = 6)

msr = "countSum"
trm = "Trimmed"

sink(sprintf('%s/PIRE_Stat_Results.txt',OUTDIR))
lm_filtered <- lm(data = tbl_libID %>%
                    filter(measure == msr, Trimming == trm, cat == "Filtered"),
                  formula = value ~ collectionTime * species) 
summary(lm_filtered)

lm_NoMotif <- lm(data = tbl_libID %>%
                   filter(measure == msr, Trimming == trm, cat == "NoMotif"),
                 formula = value ~ collectionTime * species) 
summary(lm_NoMotif)

lm_IndelCat <- lm(data = tbl_libID %>%
                    filter(measure == msr, Trimming == trm, cat == "Motif"),
                  formula = value ~ collectionTime * species) 
summary(lm_IndelCat)
sink()

#### Make tbl_libID_cat_means - Means of LibID Tbl by Cat ("Motif"    "NoMotif"  "Filtered") ####
tbl_libID_cat_means <- tbl_libID %>%
  dplyr::select(-collectionTime_Measure) %>%
  tidyr::pivot_wider(names_from = measure,
              values_from = c(value)) %>%
  dplyr::group_by(species, collectionTime, Trimming, cat) %>% 
  dplyr::summarise(tf_Spp_tr = n(),
            mean_countSum = mean(countSum),
            sd_countSum = sd(countSum),
            mean_prop = mean(prop),
            sd_prop = sd(prop)) %>%
  tidyr::pivot_longer(cols = mean_countSum:sd_prop,
               names_pattern = '(.*)_(.*)',
               names_to = c('measure2','measure')) %>%
  tidyr::pivot_wider(names_from = measure2,
              values_from = c(value)) %>%
  dplyr::mutate(Trimming = factor(Trimming, c("Untrimmed","Trimmed"))) %>%
  dplyr::mutate(partially_logged_mean = case_when(measure == "countSum" ~ log10(mean),
                                                     measure == "prop" ~ mean),
                   partially_logged_Lower_StDev = case_when(measure == "countSum" ~ log10(mean - sd),
                                           measure == "prop" ~ mean-sd),
                   partially_logged_Upper_StDev = case_when(measure == "countSum" ~ log10(mean + sd),
                                           measure == "prop" ~ mean + sd))

#### Make tbl_libID_cat2_means - Means of LibID Tbl by Cat2 ("NoInd" "MotifWithIndels" "NoMotif" "Filtered") ####
tbl_libID_cat2_means <- tbl_libID %>%
  dplyr::select(-collectionTime_Measure) %>%
  tidyr::pivot_wider(names_from = measure,
              values_from = c(value)) %>%
  dplyr::group_by(species, collectionTime, Trimming, cat2) %>% 
  dplyr::summarise(tf_Spp_tr = n(),
            mean_countSum = mean(countSum),
            sd_countSum = sd(countSum),
            mean_prop = mean(prop),
            sd_prop = sd(prop)) %>%
  tidyr::pivot_longer(cols = mean_countSum:sd_prop,
               names_pattern = '(.*)_(.*)',
               names_to = c('measure2','measure')) %>%
  tidyr::pivot_wider(names_from = measure2,
              values_from = c(value)) %>%
  dplyr::mutate(Trimming = factor(Trimming, c("Untrimmed","Trimmed"))) 


#### Graph Proportion/Number of Read Pairs from nonDMX files with mean and sd for "cat" level ####
poolBarPlot_cat_means <- function(data = tbl_libID_cat_means, msr = "prop", Category = "Filtered"){
  if(msr == "prop"){
    title = paste("Proportions of Reads", sep = "\n")
  } else if(msr == "countSum"){
    title = paste("Number of Reads", sep = "\n")
  }
  
  data %>%
    filter(measure == msr, cat == Category, species != "Olat") %>%
    ggplot(aes(x=species, y=mean)) +
    geom_bar(stat="identity", width=0.9, fill="lightgrey") +
    geom_errorbar(aes(x = species, ymin = mean-sd, ymax = mean + sd),
                  width = 0.4, colour = "black", alpha = 0.9, size = 1.3) +
    theme_classic() +
    labs(title = title) +
    facet_grid(collectionTime ~ Trimming) 
}          

pdf(file = "Output/Cat Means - Read Proportions And Counts.pdf", width = 14, height = 7)
poolBarPlot_cat_means(Category = "Filtered", msr = "prop")
poolBarPlot_cat_means(Category = "Filtered", msr = "countSum")
poolBarPlot_cat_means(Category = "Motif", msr = "prop")
poolBarPlot_cat_means(Category = "Motif", msr = "countSum")
poolBarPlot_cat_means(Category = "NoMotif", msr ="prop")
poolBarPlot_cat_means(Category = "NoMotif", msr = "countSum")
dev.off()

#### For each species make plots with mean and sd for each "cat2" level. ####
poolBarPlot_cat2_means <- function(data=tbl_libID_cat2_means, msr = "prop", Category = "Filtered"){
  if(msr == "prop"){
    title = paste("Proportions of Reads from Undemultiplexed Files", sep = "\n")
  } else if(msr == "countSum"){
    title = paste("Number of Reads from Undemultiplexed Files", sep = "\n")
  }
  
  data %>%
    filter(measure == msr, cat2 == Category, species != "Olat") %>%
    ggplot(aes(x = species, y = mean)) +
    geom_bar(stat = "identity", width = 0.9, fill = "lightgrey") +
    geom_errorbar(aes(x = species, ymin = mean-sd, ymax = mean + sd),
                  width = 0.4, colour = "black", alpha = 0.9, size = 1.3) +
    theme_classic() +
    labs(title = title) +
    facet_grid(collectionTime ~ Trimming) 
}          

pdf(file = "Output/Cat2_Means - Read_Proportions_And_Counts.pdf", width = 14, height = 7)
poolBarPlot_cat2_means(Category = "NoInd", msr= "prop")
poolBarPlot_cat2_means(Category = "NoInd", msr = "countSum")
poolBarPlot_cat2_means(Category = "MotifWithIndels", msr= "prop")
poolBarPlot_cat2_means(Category = "MotifWithIndels", msr = "countSum")
poolBarPlot_cat2_means(Category = "NoMotif", msr= "prop")
poolBarPlot_cat2_means(Category = "NoMotif", msr = "countSum")
poolBarPlot_cat2_means(Category= "Filtered", msr= "prop")
poolBarPlot_cat2_means(Category= "Filtered", msr = "countSum")
dev.off()

#### Beta regressions and run lrtest, joint_tests, & summary (proportion data) ####
Filtered.data <- tbl_libID %>% 
  filter(cat == "Filtered", Trimming == "Trimmed", measure == "prop") %>%
  group_by(species, collectionTime, libID, cat) %>%
  dplyr::summarise(value = sum(value))
Filtered.model <- betareg(formula = value ~ species + collectionTime, data = Filtered.data)
sink(sprintf('%s/BetaRegression_Filtered_Proportions.txt', OUTDIR))
lrtest(Filtered.model)
writeLines(c("\n","####################################","\n","Joint Tests","\n"))
joint_tests(Filtered.model)
writeLines(c("\n","####################################","\n","Joint Tests by Species","\n"))
joint_tests(Filtered.model, by = "species")
writeLines(c("\n","####################################","\n","Summary","\n"))
summary(Filtered.model)
sink()

Motif.data <- tbl_libID %>%
  filter(cat == "Motif", Trimming == "Trimmed", measure == "prop") %>%
  group_by(species, collectionTime, libID, cat) %>%
  dplyr::summarise(value = sum(value))
Motif.model <- betareg(formula = value ~ species + collectionTime, data = Motif.data)
sink(sprintf('%s/BetaRegression_Motif_Proportions.txt', OUTDIR))
lrtest(Motif.model)
writeLines(c("\n","####################################","\n","Joint Tests","\n"))
joint_tests(Motif.model)
writeLines(c("\n","####################################","\n","Joint Tests by Species","\n"))
joint_tests(Motif.model, by = "species")
writeLines(c("\n","####################################","\n","Summary","\n"))
summary(Motif.model)
sink()

NoMotif.data <- tbl_libID %>%
  filter(cat == "NoMotif", Trimming == "Trimmed", measure == "prop") %>%
  group_by(species, collectionTime, libID, cat) %>%
  dplyr::summarise(value = sum(value))
NoMotif.model <- betareg(formula = value ~ species + collectionTime, data = NoMotif.data)
sink(sprintf('%s/BetaRegression_NoMotif_Proportions.txt', OUTDIR))
lrtest(NoMotif.model)
writeLines(c("\n","####################################","\n","Joint Tests","\n"))
joint_tests(NoMotif.model)
writeLines(c("\n","####################################","\n","Joint Tests by Species","\n"))
joint_tests(NoMotif.model, by = "species")
writeLines(c("\n","####################################","\n","Summary","\n"))
summary(NoMotif.model)
sink()

#### Proportion Analysis Matching Error Analysis ####
# massage data into submission
propAnalysisFiltered.data <- tbl_libID %>%
  filter(Trimming == "Trimmed", measure == "countSum") %>%
  group_by(species, collectionTime, libID, cat) %>%
  dplyr::summarise(counts = sum(value)) %>%
  pivot_wider(names_from=cat, values_from = counts) %>%
  mutate(NotFiltered = Motif + NoMotif)

# statistical model testing for diffs in proportions
propAnalysisFiltered.model1 <- mixed(cbind(Filtered, NotFiltered) ~ species * collectionTime 
                                    + (1 | libID),
                                    data = propAnalysisFiltered.data, 
                                    family = binomial, 
                                    method = 'LRT')
propAnalysisFiltered.model1

propAnalysisFiltered.model2 <- mixed(cbind(Filtered, NotFiltered) ~ species * collectionTime 
                            + (1 | libID),
                            data = propAnalysisFiltered.data, 
                            family = binomial, 
                            method = 'LRT',
                            control = glmerControl(optimizer = "bobyqa"))
propAnalysisFiltered.model2

propAnalysisFiltered.model3 <- mixed(cbind(Filtered, NotFiltered) ~ species * collectionTime 
                                    + (1 | libID),
                                    data = propAnalysisFiltered.data, 
                                    family = binomial, 
                                    method = 'LRT',
                                    control = glmerControl(optimizer = "Nelder_Mead"))
propAnalysisFiltered.model3

propAnalysisFiltered.model4 <- mixed(cbind(Filtered, NotFiltered) ~ species * collectionTime + (1 | libID),
                                    data = propAnalysisFiltered.data, 
                                    family = binomial,
                                    method = 'PB')
propAnalysisFiltered.model4

# estimated means from model & test pairwaise comps
propAnalysisFiltered.emm <- emmeans(propAnalysisFiltered.model1, ~ species * collectionTime)
propAnalysisFiltered_pairwise_postHoc <- contrast(propAnalysisFiltered.emm, 
                                  method = 'pairwise', 
                                  simple = 'each', 
                                  combine = FALSE, 
                                  adjust = 'BH')
propAnalysisFiltered_pairwise_postHoc

# write statistical result table to output
sink('./Aduo_ErrorCalc_Analysis_Output.txt', append = TRUE)
writeLines(c("\n","####################################","\n","\n","> summary(Aduo_error_model)"))
summary(Aduo_error_model_2)
writeLines(c("\n","####################################","\n","\n","> anova(Aduo_error_model)"))
anova(Aduo_error_model_2)
writeLines(c("\n","####################################","\n","\n","> Aduo_pairwise_postHoc"))
Aduo_pairwise_postHoc
sink()


propAnalysisMotif.model1 <- mixed(cbind(Motif, NoMotif) ~ species * collectionTime 
                                     + (1 | libID),
                                     data = propAnalysisFiltered.data, 
                                     family = binomial, 
                                     method = 'LRT')
propAnalysisMotif.model1

propAnalysisMotif.model2 <- mixed(cbind(Motif, NoMotif) ~ species * collectionTime 
                                 + (1 | libID),
                                 data = propAnalysisFiltered.data, 
                                 family = binomial, 
                                 method = 'LRT',
                                 control = glmerControl(optimizer = "bobyqa",
                                                        optCtrl = list(maxfun = 2e8)))
propAnalysisMotif.model2

cntAnalysisMotif.model1 <- mixed(Motif ~ species * collectionTime 
                                  + (1 | libID),
                                  data = propAnalysisFiltered.data, 
                                  family = poisson, 
                                  method = 'LRT',
                                  control = glmerControl(optimizer = "bobyqa"))
cntAnalysisMotif.model1

cntAnalysisMotif.model2 <- mixed(Motif ~ species * collectionTime 
                                  + (1 | libID),
                                  data = propAnalysisFiltered.data, 
                                  family = poisson, 
                                  method = 'LRT',
                                  control = glmerControl(optimizer = "Nelder_Mead"))
cntAnalysisMotif.model1


#### Format data for Negative Binomial Regression (count data) ####
Filtered.data <- tbl_libID %>%
  dplyr::mutate(species = factor(species,c("Aduo","Olat","Sspi")),
                collectionTime = factor(collectionTime,c("Historical","Contemporary")),
                Trimming = factor(Trimming,c("Trimmed","Untrimmed"))) %>%
  filter(cat == "Filtered", Trimming == "Trimmed", species != "Olat", measure == "countSum") %>%
  dplyr::select(value,species,collectionTime)

Filtered.model.nb = glm.nb(value ~ species * collectionTime,
                           data=Filtered.data,
                           control = glm.control(maxit=10000))
sink(sprintf('%s/NegativeBinomial_Filtered_Counts.txt', OUTDIR))
writeLines(c("\n","####################################","\n","Summary","\n"))
summary(Filtered.model.nb)
writeLines(c("\n","####################################","\n","LR Anova","\n"))
Anova(Filtered.model.nb,
      type="II",
      test="LR")
sink()

Motif.data <- tbl_libID %>%
  dplyr::mutate(species = factor(species,c("Aduo","Olat","Sspi")),
                collectionTime = factor(collectionTime,c("Historical","Contemporary")),
                Trimming = factor(Trimming,c("Trimmed","Untrimmed"))) %>%
  filter(cat == "Motif", Trimming == "Trimmed", species != "Olat", measure == "countSum")

Motif.model.nb = glm.nb(value ~ species * collectionTime,
                        data=Motif.data,
                        control = glm.control(maxit=10000))
sink(sprintf('%s/NegativeBinomial_Motif_Counts.txt', OUTDIR))
writeLines(c("\n","####################################","\n","Summary","\n"))
summary(Motif.model.nb)
writeLines(c("\n","####################################","\n","LR Anova","\n"))
Anova(Motif.model.nb,
      type="II",
      test="LR")
sink()

NoMotif.data <- tbl_libID %>%
  dplyr::mutate(species = factor(species,c("Aduo","Olat","Sspi")),
                collectionTime = factor(collectionTime,c("Historical","Contemporary")),
                Trimming = factor(Trimming,c("Trimmed","Untrimmed"))) %>%
  filter(cat == "NoMotif", Trimming == "Trimmed", species != "Olat", measure == "countSum") 

NoMotif.model.nb = glm.nb(value ~ species * collectionTime,
                          data=NoMotif.data,
                          control = glm.control(maxit=10000))
sink(sprintf('%s/NegativeBinomial_NoMotif_Counts.txt', OUTDIR))
writeLines(c("\n","####################################","\n","Summary","\n"))
summary(NoMotif.model.nb)
writeLines(c("\n","####################################","\n","LR Anova","\n"))
Anova(NoMotif.model.nb,
      type="III",
      test="LR")
sink()


#### Number of Reads bar plot for number of reads by cat, collectionTime, measure, and species. ####
MeanReadBarPlot <- function(data=tbl_libID_cat_means){
  title = paste("Number of Reads", sep ="\n")
  data %>%
    replace_na(list(partially_logged_Lower_StDev = 0)) %>%
    filter(species != "Olat", measure == "countSum", Trimming == "Trimmed") %>%
    ggplot(aes(x=species, y=mean)) +
    geom_col(color="lightgrey") +
    scale_fill_brewer(palette="Set1") +
    theme_classic() +
    geom_errorbar(aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.3) +
    labs(title = title) +
    facet_grid(collectionTime ~ cat) 
}
MeanReadBarPlot()
ggsave('Output/Mean_Read_Counts.png', height = 6, width = 8)


#### Read Proportions bar plot for number of reads by cat, collectionTime, measure, and species. ####
MeanPropBarPlot <- function(data=tbl_libID_cat_means){
  title = paste("Proportion of Reads", sep ="\n")
  data %>%
    replace_na(list(partially_logged_Lower_StDev = 0)) %>%
    filter(species != "Olat", measure == "prop", Trimming == "Trimmed") %>%
    ggplot(aes(x=species, y=mean)) +
    geom_col(color="lightgrey") +
    scale_fill_brewer(palette="Set1") +
    theme_classic() +
    geom_errorbar(aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.3) +
    labs(title = title) +
    facet_grid(collectionTime ~ cat) 
}
MeanPropBarPlot()
ggsave('Output/Mean_Read_Proportions.png', height = 6, width = 8)


#### 12 panel error bar plot of reads and proportions by cat, collectionTime, measure, and species. #### 
MeanReadCountAndPropBarplot <- function(data=tbl_libID_cat_means){
  # title = paste("Read Counts and Proportions", sep ="\n")
  data %>%
    replace_na(list(partially_logged_Lower_StDev = 0)) %>%
    filter(species != "Olat", Trimming == "Trimmed") %>%
    mutate(time_measure = paste(collectionTime,measure,sep = "_"),
           partially_logged_mean = case_when(measure == "countSum" ~ mean,
                                             measure == "prop" ~ mean)) %>%
    mutate(time_measure = factor(time_measure, c("Contemporary_countSum", "Historical_countSum" ,"Contemporary_prop", "Historical_prop"))) %>%
    ggplot(aes(x=species, y=mean)) +
    geom_bar(stat="identity", width=0.9, fill="lightgrey") +
    geom_errorbar(aes(x=species, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.3) +
    theme_classic() +
    theme(text = element_text(size = 20), axis.text.y = element_text(size = 12), axis.title.y = element_blank(),axis.title.x = element_blank(),strip.text.y = element_blank(),strip.text.x = element_blank()) +
    # labs(title = title) +
    xlab("Species") +
    facet_grid(time_measure ~ cat, scales = "free") 
}          
MeanReadCountAndPropBarplot()
ggsave('Output/Mean_Read_Counts_and_Proportions.png', height = 6.5, width = 10)

