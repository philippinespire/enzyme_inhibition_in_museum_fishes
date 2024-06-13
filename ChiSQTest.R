setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(janitor)
library(broom)

# read_csv("./Aduo_Contemp_Chi_matrix.csv") -> Aduo_Contemp_Chi_matrix
# read_csv("./Sspi_Contemp_Chi_matrix.csv") -> Sspi_Contemp_Chi_matrix
# read_csv("./Aduo_Museum_Chi_matrix.csv") -> Aduo_Museum_Chi_matrix
# read_csv("./Sspi_Museum_Chi_matrix.csv") -> Sspi_Museum_Chi_matrix
# 
# raw_data <- bind_rows(Aduo_Contemp = Aduo_Contemp_Chi_matrix,
#                       Sspi_Contemp = Sspi_Contemp_Chi_matrix,
#                       Aduo_Museum = Aduo_Museum_Chi_matrix,
#                       Sspi_Museum = Sspi_Museum_Chi_matrix,
#                       .id = 'Species_Time') %>%
#   separate(Species_Time, into = c('species', 'time')) %>%
#   mutate(id = as.factor(row_number()))

raw_data <- read_csv("../Chi_Table_ceb.csv") %>%
  clean_names() %>%
  rename(time = sampling) %>%
  mutate(id = as.factor(row_number()))

summarized_data <- raw_data %>%
  group_by(species, time) %>%
  summarize(obs_err15 = sum(err15),
            obs_err16 = sum(err16),
            obs_total = sum(total)) %>%
  mutate(pos15_err = obs_err15 / obs_total,
         pos16_err = obs_err16 / obs_total) %>%
  select(-starts_with("obs"))

summarized_data_2 <- raw_data %>%
  group_by(species, time, barcode) %>%
  summarize(obs_err15 = sum(err15),
            obs_err16 = sum(err16),
            obs_total = sum(total)) %>%
  mutate(pos15_err = obs_err15 / obs_total,
         pos16_err = obs_err16 / obs_total) %>%
  select(-starts_with("obs"))

#### Method from Zar for testing if A) pooling is "ok" and B) testing goodness-of-fit given expected error rates ####

sequence_run_chisq <- raw_data %>%
  left_join(summarized_data) %>%
  # mutate(pos15_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.02207,
  #                              species == 'Sspi' & time == 'Contemporary' ~ 0.01021,
  #                              species == 'Aduo' & time == 'Museum' ~ 0.05768,
  #                              species == 'Aduo' & time == 'Contemporary' ~ 0.00924),
  #        
  #        pos16_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.04695,
  #                              species == 'Sspi' & time == 'Contemporary' ~ 0.02547,
  #                              species == 'Aduo' & time == 'Museum' ~ 0.11971,
  #                              species == 'Aduo' & time == 'Contemporary' ~ 0.02334)) %>%
  rowwise() %>%
  mutate(n = sum(c_across(starts_with('observed'))),
         expected_freq = list(c((1 - pos15_err) * (1 - pos16_err), 
                                pos15_err * (1 - pos16_err), 
                                (1 - pos15_err) * pos16_err,
                                pos15_err * pos16_err))) %>%
  # filter(n * min(expected_freq) > 5) %>% #remove sequencing runs with very few reads
  filter(n > 20) %>%  #remove sequencing runs with very few reads
  mutate(chisq = list(chisq.test(x = c_across(starts_with('observed')), p = expected_freq, correct = FALSE)),
         chisq = list(tidy(chisq))) %>%
  unnest(chisq) %>%
  select(-p.value, -method, -pos15_err, -pos16_err) %>%
  group_by(species, time) %>%
  summarise(across(where(is.numeric), sum),
            expected_freq = unique(expected_freq),
            .groups = 'drop') %>%
  rename(total_chisq = statistic,
         total_df = parameter) %>%
  rowwise %>%
  mutate(chisq = list(chisq.test(x = c_across(starts_with('observed')), p = expected_freq, correct = FALSE)),
         chisq = list(tidy(chisq))) %>%
  unnest(chisq) %>%
  select(-method) %>%
  rename(pool_chisq = statistic,
         pool_df = parameter,
         pool_p = p.value) %>%
  mutate(hetero_chi = total_chisq - pool_chisq,
         hetero_df = total_df - pool_df) %>%
  mutate(hetero_p = pchisq(hetero_chi, hetero_df, lower.tail = FALSE)) %>%
  select(species, time, starts_with('total'), starts_with('pool'), starts_with('hetero'))

sequence_run_chisq

## Cannot assume that different sequencing runs come from the same population as the "hetero_chi" is significant and can't reject that null hypothesis
## Therefore the pooled chisquared test shouldn't be used as you are artificially merging multiple populations and it is unclear if the significance is related to that or the tested factors.

#### Bayesian model 1 - this is the best model ####
## Multinomial regression 
library(brms)
library(tidybayes)
library(emmeans)

test_data <- raw_data %>%
  rowwise %>%
  mutate(n = sum(c_across(starts_with('observed')))) %>%
  ungroup %>%
  as.data.frame()

test_data$y <- with(test_data, cbind(observed_neither_err, observed_only15, observed_only16, observed_both_err))

# multinom_test <- brm(bf(y | trials(n) ~ species * time + (1 | fq_base * barcode) ),
multinom_test <- brm(bf(y | trials(n) ~ species * time + (1 | id) ),
                     data = test_data,
                     family = multinomial(link = 'logit', refcat = 'observed_neither_err'),
                     chains = 4,
                     iter = 5000,
                     warmup = 2500,
                     cores = 4,
                     control = list(max_treedepth = 20, adapt_delta=0.8),
                     file = 'multinom_model_cl')
# 
# multinom_test <- brm(bf(y | trials(n) ~ species * time + (1 | fq_base) + (1 | barcode) ),
#                      data = test_data,
#                      family = multinomial(link = 'logit', refcat = 'observed_neither_err'),
#                      chains = 4,
#                      iter = 5000,
#                      warmup = 2500,
#                      cores = 4,
#                      control = list(max_treedepth = 20, adapt_delta=0.9999),
#                      file = 'multinom_model_cl_999')

# multinom_test <- brm(bf(y | trials(n) ~ species * time * barcode + (1 | fq_base) ),
#                      data = test_data,
#                      family = multinomial(link = 'logit', refcat = 'observed_neither_err'),
#                      chains = 4,
#                      iter = 5000,
#                      warmup = 2500,
#                      cores = 4,
#                      control = list(max_treedepth = 20, adapt_delta=0.80),
#                      file = 'multinom_model_cl2_80')

summary(multinom_test)
plot(multinom_test) #confirm trace plots are all converged


raw_data %>%
  select(species, time) %>%
  distinct() %>%
  mutate(n = 10000) %>%
  add_fitted_draws(multinom_test, re_formula = NA) %>%
  point_interval(.width=0.99) %>%
  ungroup %>%
  rename(error_type = .category) %>%
  left_join(summarized_data) %>% 
  mutate(expected = case_when(error_type == 'observed_both_err' ~ pos15_err * pos16_err,
                              error_type == 'observed_neither_err' ~ (1 - pos15_err) * (1 - pos16_err),
                              error_type == 'observed_only15' ~ pos15_err * (1 - pos16_err),
                              error_type == 'observed_only16' ~ (1 - pos15_err) * pos16_err),
         expected = expected * n) %>%
  mutate(error_type = factor(error_type, levels = c('observed_neither_err', 'observed_only15', 'observed_only16', 'observed_both_err'))) %>%
  mutate(signficantly_different_than_expected = expected < .lower | expected > .upper,
         marks = if_else(signficantly_different_than_expected, '*', ''),
         upper_mark = 1.1 * max(.upper),
         direction = case_when(expected < .lower ~ 'more errors',
                               expected > .upper ~ 'fewer errors',
                               TRUE ~ 'no diff')) %>%
  
  
  
  ggplot(aes(x = error_type, y = .value, ymin = .lower, ymax = .upper)) +
  geom_col(aes(y = expected), fill = "orange3") +
  geom_linerange() +
  geom_point() +
  geom_text(aes(y = upper_mark, label = marks, colour = direction)) +
  scale_y_log10() +
  facet_grid(time ~ species) +
  scale_colour_manual(values = c('more errors' = 'red', 'no diff' = 'white', 'fewer errors' = 'green'), 
                      breaks = c('more errors', 'fewer errors')) +
  labs(x = NULL,
       y = 'Number of Errors per 10,000 Reads') +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        strip.text.x = element_blank(), 
        strip.text.y = element_blank(),
        legend.position = "none") 
ggsave(filename = './NumberofErrors15_16.png', width = 10, height = 6.5)



raw_data %>%
  select(species, time) %>%
  # select(species, time, barcode) %>%
  distinct() %>%
  mutate(n = 1) %>%
  add_fitted_draws(multinom_test, re_formula = NA) %>%
  point_interval() %>%
  
  ungroup %>%
  rename(error_type = .category) %>%
  left_join(summarized_data) %>%  #barcode not a fixed factor
  # left_join(summarized_data_2) %>%  #barcode a fixed factor
  # mutate(pos15_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.02207,
  #                              species == 'Sspi' & time == 'Contemporary' ~ 0.01021,
  #                              species == 'Aduo' & time == 'Museum' ~ 0.05768,
  #                              species == 'Aduo' & time == 'Contemporary' ~ 0.00924),
  #        
  #        pos16_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.04695,
  #                              species == 'Sspi' & time == 'Contemporary' ~ 0.02547,
  #                              species == 'Aduo' & time == 'Museum' ~ 0.11971,
  #                              species == 'Aduo' & time == 'Contemporary' ~ 0.02334)) %>%
  
  mutate(expected = case_when(error_type == 'observed_both_err' ~ pos15_err * pos16_err,
                              error_type == 'observed_neither_err' ~ (1 - pos15_err) * (1 - pos16_err),
                              error_type == 'observed_only15' ~ pos15_err * (1 - pos16_err),
                              error_type == 'observed_only16' ~ (1 - pos15_err) * pos16_err),
         expected = expected * n) %>%
  select(species, time, error_type, .value:.upper, expected)
  # select(species, time, barcode, error_type, .value:.upper, expected)


# These are all on the link function scale to get the estimates/contrasts on the response scale you need to inverse logit transform
emmeans(multinom_test, pairwise ~ species, dpar = 'muobservedonly15', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ time, dpar = 'muobservedonly15', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species * time, dpar = 'muobservedonly15', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species, dpar = 'muobservedonly15', level=0.99, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ time, dpar = 'muobservedonly15', level=0.99, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species * time, dpar = 'muobservedonly15', level=0.99, adjust="bonferroni") 


emmeans(multinom_test, pairwise ~ species, dpar = 'muobservedonly16', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ time, dpar = 'muobservedonly16', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species * time, dpar = 'muobservedonly16', level=0.95, adjust="bonferroni")
emmeans(multinom_test, pairwise ~ species, dpar = 'muobservedonly16', level=0.99, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ time, dpar = 'muobservedonly16', level=0.99, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species * time, dpar = 'muobservedonly16', level=0.99, adjust="bonferroni") 


emmeans(multinom_test, pairwise ~ species, dpar = 'muobservedbotherr', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ time, dpar = 'muobservedbotherr', level=0.95, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species * time, dpar = 'muobservedbotherr', level=0.95, adjust="bonferroni")
emmeans(multinom_test, pairwise ~ species, dpar = 'muobservedbotherr', level=0.99, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ time, dpar = 'muobservedbotherr', level=0.99, adjust="bonferroni") 
emmeans(multinom_test, pairwise ~ species * time, dpar = 'muobservedbotherr', level=0.99, adjust="bonferroni") 

# If positive, the first is bigger, if negative, the second is bigger

#### Bayesian model 2 - not as good a model ####
## Negative binomial regression.
# This isn't as good as the above since it doesn't account for the fact that each read has a chance of falling into one of these 4 categories. 
# basically the above model if there are 50 total reads then the sum of all the different error types must sum to 50. The model below doesn't force this fact


brm_mod <- raw_data %>%
  rowwise %>%
  mutate(n = sum(c_across(starts_with('observed')))) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('observed'),
               names_to = 'error_type',
               values_to = 'count') %>%
  brm(count ~ error_type * time * species + (1 + error_type | id) + offset(log(n)),
      data = ., 
      family = negbinomial(),
      chains = 4,
      iter = 5000,
      warmup = 2500,
      cores = 4, 
      control = list(max_treedepth = 20))


summary(brm_mod)
plot(brm_mod) #confirm trace plots are all converged

raw_data %>%
  pivot_longer(cols = starts_with('observed'),
               names_to = 'error_type',
               values_to = 'count') %>%
  select(species, time, error_type) %>%
  distinct %>%
  mutate(n = 10000) %>%
  add_fitted_draws(brm_mod, re_formula = NA) %>%
  point_interval() %>%
  ungroup %>%
  
  mutate(pos15_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.02207,
                               species == 'Sspi' & time == 'Contemporary' ~ 0.01021,
                               species == 'Aduo' & time == 'Museum' ~ 0.05768,
                               species == 'Aduo' & time == 'Contemporary' ~ 0.00924),
         
         pos16_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.04695,
                               species == 'Sspi' & time == 'Contemporary' ~ 0.02547,
                               species == 'Aduo' & time == 'Museum' ~ 0.11971,
                               species == 'Aduo' & time == 'Contemporary' ~ 0.02334)) %>%

  mutate(expected = case_when(error_type == 'observed_both_err' ~ pos15_err * pos16_err,
                              error_type == 'observed_neither_err' ~ (1 - pos15_err) * (1 - pos16_err),
                              error_type == 'observed_only15' ~ pos15_err * (1 - pos16_err),
                              error_type == 'observed_only16' ~ (1 - pos15_err) * pos16_err),
         expected = expected * n) %>%
  mutate(error_type = factor(error_type, levels = c('observed_neither_err', 'observed_only15', 'observed_only16', 'observed_both_err'))) %>%
  ggplot(aes(x = error_type, y = .value, ymin = .lower, ymax = .upper)) +
  geom_col(aes(y = expected)) +
  geom_linerange() +
  geom_point() +
  scale_y_log10() +
  facet_grid(time ~ species) +
  labs(x = 'Error Type',
       y = 'Number of Errors per 10,000 Reads') +
  theme_classic()

# If expected value is within credible interval then it is not significantly different from what you expect. If outside it is significantly different from the prediction based on an assumed independent error rate.

#### Can also perform post-hoc tests not possible with chi-squared test version. e.g. differences in error rates by species/time point
emmeans(brm_mod, pairwise ~ species * time, by = 'error_type', offset = log(10000), type = 'response') #if 1 is in the confidence interval then it is non-significant difference


#### reformate Bayesian test 2 as a frequentist test ####
library(lme4)
library(emmeans)
lmer_predict <- function(x, model, N){
  make_preditions <- function(.) {
    predict(., newdata = x, re.form = ~0)
  }
  
  x$fit <- predict(model, newdata = x, re.form = ~0, type = 'link')
  
  tmp <- bootMer(model, FUN = make_preditions, nsim = N, re.form = ~0)
  
  x$se <- apply(tmp$t, 2, sd, na.rm = TRUE)
  
  x$lwr <- x$fit - 1.96 * x$se
  x$upr <- x$fit + 1.96 * x$se
  
  x
}

# This takes a long time and probably won't converge
lmer_mod <- raw_data %>%
  rowwise %>%
  mutate(n = sum(c_across(starts_with('observed')))) %>%
  ungroup %>%
  pivot_longer(cols = starts_with('observed'),
               names_to = 'error_type',
               values_to = 'count') %>%
  glmer.nb(count ~ error_type * time * species + (1 | id) + offset(log(n)),
      data = .)

summary(lmer_mod)


# This takes a while - 10 should be increased to more like 100
lmer_predictions <- raw_data %>%
  pivot_longer(cols = starts_with('observed'),
               names_to = 'error_type',
               values_to = 'count') %>%
  select(species, time, error_type) %>%
  distinct %>%
  mutate(n = 10000) %>%
  lmer_predict(lmer_mod, 10)

lmer_predictions %>%
  mutate(across(c(fit, lwr, upr), ~exp(.))) %>%
  
  mutate(pos15_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.02207,
                               species == 'Sspi' & time == 'Contemporary' ~ 0.01021,
                               species == 'Aduo' & time == 'Museum' ~ 0.05768,
                               species == 'Aduo' & time == 'Contemporary' ~ 0.00924),
         
         pos16_err = case_when(species == 'Sspi' & time == 'Museum' ~ 0.04695,
                               species == 'Sspi' & time == 'Contemporary' ~ 0.02547,
                               species == 'Aduo' & time == 'Museum' ~ 0.11971,
                               species == 'Aduo' & time == 'Contemporary' ~ 0.02334)) %>%

  mutate(expected = case_when(error_type == 'observed_both_err' ~ pos15_err * pos16_err,
                              error_type == 'observed_neither_err' ~ (1 - pos15_err) * (1 - pos16_err),
                              error_type == 'observed_only15' ~ pos15_err * (1 - pos16_err),
                              error_type == 'observed_only16' ~ (1 - pos15_err) * pos16_err),
         expected = expected * n) %>%
  mutate(error_type = factor(error_type, levels = c('observed_neither_err', 'observed_only15', 'observed_only16', 'observed_both_err'))) %>%
  ggplot(aes(x = error_type, y = fit, ymin = lwr, ymax = upr)) +
  geom_col(aes(y = expected)) +
  geom_linerange() +
  geom_point() +
  scale_y_log10() +
  facet_grid(time ~ species) +
  labs(x = 'Error Type',
       y = 'Number of Errors per 10,000 Reads') +
  theme_classic()


emmeans(lmer_mod, pairwise ~ species * time, by = 'error_type', offset = log(10000), type = 'response')
