# This script was developed by Luc Bussière
# In collaboration with all the authors of the following MS

# 
# Exploiting pathogen defence trade-offs to manage risks of crop pests 
# evolving resistance to biocontrol. 
# 
# 
# Authors: Rosie Mangan1 *, Matthew C. Tinsley1, Ester Ferrari2, 
# Ricardo A. Polanczyk3 and Luc F. Bussière4
# 
# 1 Biological and Environmental Sciences, School of Natural Sciences,
# University of Stirling, Stirling, FK94LA, UK.
# 2 Crops, Environment and Land Use Programme, Teagasc: 
# The Agriculture and Food Development Authority, 
# Ashtown, Dublin, D15DY05, Ireland.
# 3 Júlio de Mesquita Filho State University of São Paulo, 
# Faculty of Agrarian and Veterinary Sciences of Jaboticabal, 
# Jaboticabal, SP Brazil.
# 4 Biological and Environmental Sciences & 
# Gothenburg Global Biodiversity Centre, 
# The University of Gothenburg, 405 30 Gothenburg, Sweden.

# Script originated May 1, 2023
# Last modified Aug 23, 2024


# The next line is to clear the memory, but it is annotated out 
# to prevent accidental erasure of the big models that take a lot of time to find
rm(list = ls())

# load libraries ####
library(lme4)
library(brms)
library(lubridate)
library(graph4lg)
library(bayestestR)
library(performance)
library(reshape2)
library(data.table)
library(coda)
library(ggpubr)
library(Hmisc)
library(png)
library(ggtext)
library(binom)
library(parallel)
library(ggridges)
library(pbkrtest)
# library(QGglmm)
library(broom.mixed)
library(RColorBrewer)
library(tidybayes)

# The script uses some functions from the rethinking package by McElreath, 
# which can be downloaded using the annotated code below.
# install.packages("rethinking", 
#                  repos=c(cran="https://cloud.r-project.org",
#                          rethinking="http://xcelab.net/R"))
# library(rethinking)
library(broom)
library(tidyverse)

# set the default ggplot theme for graphics
theme_set(theme_bw())



# preamble ####
## Import data ####
# read data in frame and examine it
df_qgen_pre<-read_csv2("data/Quan_Gene_APRIL_12_2023_RM_LFB.csv")

names(df_qgen_pre)
df_qgen_pre %>% 
  filter(is.na(inoculationdate))




## housekeeping ####

#remove NAs inoculation dates
df_qgen<-df_qgen_pre %>% 
  filter(is.na(inoculationdate)==FALSE)

summary(df_qgen) #Summary of data

# CHANGE VARIABLES from num to factor
df_qgen$replicate<-as.factor(df_qgen$replicate)
df_qgen$matingpair<-as.factor(df_qgen$matingpair)

# MAKE NEW VARIABLES
df_qgen<-df_qgen %>%
  mutate(mother = paste("D",damID,sep="_"),
         father = paste("S",sireID,sep="_"),
         id = paste("O",replicate,sep="_")) %>%
  mutate(animal=id) %>% 
  select(-sireID,-damID)

df_qgen$mother<-as.factor(df_qgen$mother)
df_qgen$father<-as.factor(df_qgen$father)
df_qgen$animal<-as.factor(df_qgen$animal)

df_qgen$id<-as.factor(df_qgen$id)
df_qgen$plant<-as.factor(df_qgen$plant)
df_qgen$plant<-relevel(df_qgen$plant,ref="Tom")

df_qgen$isolate<-as.factor(df_qgen$isolate)
df_qgen$isolate<-relevel(df_qgen$isolate,ref="CON")

# can choose not to recode treatment as factor to suppress warnings later
df_qgen$treatment<-as.factor(df_qgen$treatment)
df_qgen$isolate<-relevel(df_qgen$isolate,ref="CON")

df_qgen$Pupal_sex<-as.factor(df_qgen$Pupal_sex)

names(df_qgen)



# MORTdayTENin the data file incorrectly does not include larvae that died on d10, 
# and that all NAs not returned in ifelse statements, 
# so these lines fix it and generate new mortality vars on several days (not just on d10)

## recreate death vars ####
df_qgen <- df_qgen %>% 
  group_by(plant, isolate) %>% 
  mutate(MORTdayFIVE = if_else(
    is.na(Days_to_death),
    0,
    if_else(Days_to_death < 6 & Days_to_death > 0, 
            1, 
            0))) %>%
  mutate(MORTdaySIX = if_else(
    is.na(Days_to_death),
    0,
    if_else(Days_to_death < 7 & Days_to_death > 0, 
            1, 
            0))) %>%
  mutate(MORTdaySEVEN = if_else(
    is.na(Days_to_death),
    0,
    if_else(
      Days_to_death < 8 & Days_to_death > 0, 
      1, 
      0))) %>% 
  mutate(MORTdayEIGHT = if_else(
    is.na(Days_to_death),
    0,
    if_else(
      Days_to_death < 9 & Days_to_death > 0, 
      1, 
      0))) %>% 
  mutate(MORTdayNINE = if_else(
    is.na(Days_to_death),
    0,
    if_else(
      Days_to_death < 10 & Days_to_death > 0, 
      1, 
      0))) %>% 
  mutate(MORTdayTEN = if_else(
    is.na(Days_to_death),
    0,
    if_else(
      Days_to_death < 11 & Days_to_death > 0, 
      1, 
      0))) %>% 
  mutate(MORTday14 = if_else(
    is.na(Days_to_death),
    0,
    if_else(
      Days_to_death < 15 & Days_to_death > 0, 
      1, 
      0))) 



## trim to essential cols ####
# variables plus exclude oil and early deaths
df_qgen_sub <- df_qgen %>% 
  mutate(Exclude = if_else(
    Notes %in% c("oil death",
                 "Oil death",
                 "missing",
                 "squashed",
                 "Escaped"), 1, 0)) %>% 
  filter(Exclude == 0) %>% 
  filter(Days_to_death != 1 & Days_to_death != 2 | is.na(Days_to_death)) %>%
  select(animal,mother,father,plant,isolate,treatment,MORTday14)


df_qgen %>% 
  mutate(oildeath = if_else(
    Notes %in% c("oil death",
                 "Oil death"), 1, 0)) %>% 
  group_by(oildeath) %>% 
  summarise(N = n())

# check that this catches all?
df_qgen %>% 
  mutate(oildeath = if_else(
    Notes %in% c("oil death",
                 "Oil death"), 1, 0)) %>% 
  filter(oildeath == 0) %>% 
  filter(Days_to_death != 1 & Days_to_death != 2 | is.na(Days_to_death)) %>% 
group_by(oildeath) %>% 
  summarise(N = n())



df_qgen_wblock <- df_qgen %>% 
  mutate(Exclude = if_else(
    Notes %in% c("oil death",
                 "Oil death",
                 "missing",
                 "squashed",
                 "Escaped"), 1, 0)) %>% 
  filter(Exclude == 0) %>% 
  filter(Days_to_death != 1 & Days_to_death != 2 | is.na(Days_to_death)) %>%
  select(animal,mother,father,plant,isolate,treatment,MORTday14, block)


# convert df to data.frame
df_qgen_sub_df <- as.data.frame(df_qgen_sub)


# write csv for records
# df_qgen_sub_df %>% 
#   write_csv("outputs/ManganetalQG2023.csv")



## reorder levels in main df ####
names(df_qgen_sub_df)
levels(factor(df_qgen_sub_df$plant))
levels(factor(df_qgen_sub_df$isolate))
levels(factor(df_qgen_sub_df$treatment))
df_qgen_sub_df <- df_qgen_sub_df %>% 
  mutate(plant = fct_relevel(plant, "Soya", "Maz", "Tom"),
         isolate = fct_relevel(isolate, "CON", "BB", "MT"),
         treatment = fct_relevel(treatment, "Soya_CON", "Soya_BB",  "Soya_MT", 
                                 "Maz_CON", "Maz_BB", "Maz_MT",
                                 "Tom_CON", "Tom_BB","Tom_MT" ))
names(df_qgen_wblock)  
df_qgen_wblock <- df_qgen_wblock %>% 
  mutate(plant = fct_relevel(plant, "Soya", "Maz", "Tom"),
         isolate = fct_relevel(isolate, "CON", "BB", "MT"),
         treatment = fct_relevel(treatment, "Soya_CON", "Soya_BB",  "Soya_MT", 
                                 "Maz_CON", "Maz_BB", "Maz_MT",
                                 "Tom_CON", "Tom_BB","Tom_MT" ))

names(df_qgen_wblock)
df_sum_qgfams <- df_qgen_wblock %>% 
  group_by(mother,father, block) %>% 
  summarise(N = n()) %>% 
  arrange(block,father)


# model of plant and isolate ####

# fixed version of model to check overdispersion
glm_mortality_fixed <- glm(MORTday14 ~ 
                             plant * 
                             isolate,
                           data = df_qgen_sub_df,
                           family = "binomial"
)

summary(glm_mortality_fixed)
# no overdispersion in fixed model



# generalised mixed model
glm_mortality <- glmer(MORTday14 ~ 
                         plant * 
                         isolate +
                         (1|mother) +
                         (1|father),
                       data = df_qgen_sub_df,
                       family = "binomial"
)

summary(glm_mortality)

glm_mortality_noint <- update(glm_mortality,
                              ~. - plant:isolate)

check_model(glm_mortality)

# the following lines compute the probability associated with the interaction 
# using parametric bootstrapping, but are annotated out to save processing time
# modcomp_glmint <- PBmodcomp(glm_mortality,
#                             glm_mortality_noint, nsim = 500)
# modcomp_glmint

# with 500 sims, p val is near minimum possible with 500 sims, 0.003



tidy_glm_mort <- tidy(glm_mortality) %>% 
  mutate(across(where(is.numeric), round, 3))

# the next line saves the table as a csv
# write.csv(tidy_glm_mort, "outputs/glm_mort.csv")


# fig for mortalities by treatment ####

names(df_qgen_sub_df)



df_qgen_sumbytreat <- df_qgen_sub_df %>% 
  group_by(treatment) %>% 
  summarise(Famsize = n(),
            Dead = sum(MORTday14, 
                       na.rm = TRUE),
            Alive = Famsize - Dead,
            prop_mort = mean(MORTday14, 
                             na.rm = TRUE))




df_qgen_trt_cis <-  binom.confint(
  x = df_qgen_sumbytreat$Dead, 
  n = df_qgen_sumbytreat$Famsize,
  methods = "exact") %>% 
  rename(PropMort = mean,
         FamSize2 = n)  

df_qgen_trt_data <- bind_cols(df_qgen_sumbytreat,                               df_qgen_trt_cis) %>%
  select(-method, -x, -PropMort, -FamSize2)  %>% 
  mutate(treatment2 = treatment) %>% 
  separate(treatment2, 
           into = c("Host_plant",
                    "Pathogen_trt"))


levels(factor(df_qgen_trt_data$treatment))
df_qgen_trt_data$treatment <- factor(
  df_qgen_trt_data$treatment, 
  levels = c("Soya_CON","Soya_BB","Soya_MT",
             "Maz_CON","Maz_BB","Maz_MT",
             "Tom_CON","Tom_BB","Tom_MT"))

levels(factor(df_qgen_trt_data$Host_plant))
df_qgen_trt_data$Host_plant <- factor(
  df_qgen_trt_data$Host_plant, 
  levels = c("Soya","Maz","Tom"))

levels(factor(df_qgen_trt_data$Pathogen_trt))
df_qgen_trt_data$Pathogen_trt <- factor(
  df_qgen_trt_data$Pathogen_trt, 
  levels = c("CON","BB","MT"))


df_qgen_trt_data %>% 
  ggplot(aes(x = treatment, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower,
                      ymax = upper)) +
  labs(y = "Mortality",
       x = "Treatment") +
  scale_x_discrete(labels = c(
    "Soya Control","Soya Beauveria","Soya Metarhizium",
    "Maize Control","Maize Beauveria","Maize Metarhizium",
    "Tomato Control","Tomato Beauveria","Tomato Metarhizium"))


Host_names <- c(
  Soya = "Soybean",
  Maz = "Maize",
  Tom = "Tomato")



#### plot mortalities #######
df_qgen_trt_data %>% 
  ggplot(aes(x = Pathogen_trt, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower,
                      ymax = upper,
                      fill = Host_plant),
                  pch = 21,
                  size = 2) +
  labs(y = "Mortality",
       x = "Treatment") +
  facet_wrap(~ Host_plant,
             labeller = labeller(Host_plant = Host_names)) +
  #scale_x_discrete(labels=expression(Control,italic(Beauveria),italic(Metarhizium))) +
  scale_x_discrete(labels=expression("","", "")) +
  scale_fill_manual(values = c("lightgreen","gold","red")) +
  lims(y = c(0,1)) +
  theme(legend.position = "none",
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text = element_text(size = 24),
        strip.background = element_rect(fill="white" ))


# ggsave("outputs/MortalityPhenotypes.png",
#        width = 20,
#        height = 20,
#        units = "cm")

df_qgen_trt_data %>% 
  ggplot(aes(x = Pathogen_trt, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower,
                      ymax = upper,
                      fill = Host_plant),
                  pch = 21) +
  labs(y = "Mortality",
       x = "Treatment") +
  facet_wrap(. ~ Host_plant,
             labeller = labeller(Host_plant = Host_names)) +
  scale_x_discrete(labels=expression(Control,italic(Beauveria),italic(Metarhizium))) +
  scale_fill_manual(values = c("lightgreen","gold","red")) +
  lims(y = c(0,1)) +
  theme(legend.position = "none")



# brms modelling

# these models take ages to run, so am providing a data object that includes the long chain models for convenience
# can unannotate and run for a few iterations to verify that these work

# load data object containing models ####
# if you want to examine the Bayesian models quickly without running your computer for days, 
# load the following data object
# load("ManganBigModels.RData")

# 9 trait brms model ####

# this is annotated out because it takes long to run
## long chain 9 trait model? ####

# to ensure that results are reproducible:
# set.seed(1234)
# model.9traits.brm.matnocov.long <-
#   brm(MORTday14 ~ treatment +
#         (treatment-1| p | gr(father)) +
#         (1 | mother),
#       data = df_qgen_sub_df,
#       family = bernoulli(),
#       chains = 8, cores = 8,
#       iter = 8000,
#       control = list(adapt_delta = 0.96))


summary(model.9traits.brm.matnocov.long)

prior_summary(model.9traits.brm.matnocov.long)

pp_check(model.9traits.brm.matnocov.long)
mcmc_plot(model.9traits.brm.matnocov.long)


plot(model.9traits.brm.matnocov.long)



# Extract variances ####
# VCVarray9traits <- VarCorr(model.9traits.brm.matnocov.long)

VCVarray9traits <- VarCorr(model.9traits.brm.matnocov.long,
                           probs = c(0.11,0.89))


## reorder levels ####
VCVarray9traits

names(VCVarray9traits$father)

VCVarray9traits$father$cor

VCVarray9traits$father$cov


GCovs9traits <- melt(VCVarray9traits) %>% 
  filter(L2 == "cov",
         Var2 == "Estimate") %>% 
  rename(Cov = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Cov,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)

Q11Covs9traitswide <- melt(VCVarray9traits) %>% 
  filter(L2 == "cov",
         Var2 == "Q11") %>% 
  rename(Q11 = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Q11,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)

Q89Covs9traitswide <- melt(VCVarray9traits) %>% 
  filter(L2 == "cov",
         Var2 == "Q89") %>% 
  rename(Q89 = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Q89,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)

UpperHalfCovmat9traitswide <- GCovs9traits %>% 
  column_to_rownames(var = "Trait1") %>% 
  as.matrix(.)

UpperHalfCovmat9traitswide[lower.tri(UpperHalfCovmat9traitswide)] <- NA




## Create VCV matrix ####
VCVarray9traits$father$cor[,,"treatmentTom_MT"]

Gmat9traitswide <- melt(VCVarray9traits) %>% 
  filter(L2 == "cor",
         Var2 == "Estimate") %>% 
  rename(Cor = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Cor,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)


order_levels <- c("treatmentSoya_CON", 
                  "treatmentSoya_BB",
                  "treatmentSoya_MT",
                  "treatmentMaz_CON",
                  "treatmentMaz_BB",
                  "treatmentMaz_MT",
                  "treatmentTom_CON",
                  "treatmentTom_BB",
                  "treatmentTom_MT" )

# convert to a matrix with row and column names


Gmat9traitswide_mat <- Gmat9traitswide %>% 
  column_to_rownames(var = "Trait1") %>% 
  as.matrix(.)


Gvar9traitswide <- melt(VCVarray9traits) %>% 
  filter(L2 == "sd",
         Var2 == "Estimate") %>% 
  rename(Cor = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Cor,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)

Gmat9traitswide <- Gmat9traitswide %>% 
  mutate(Trait1 = factor(Trait1, 
                         levels = c("treatmentSoya_CON", 
                                    "treatmentSoya_BB",
                                    "treatmentSoya_MT",
                                    "treatmentMaz_CON",
                                    "treatmentMaz_BB",
                                    "treatmentMaz_MT",
                                    "treatmentTom_CON",
                                    "treatmentTom_BB",
                                    "treatmentTom_MT" )))



# convert to half of the matrix
LowerHalfGmat9traitswide <- Gmat9traitswide_mat



# remove half of the matrix
LowerHalfGmat9traitswide[upper.tri(LowerHalfGmat9traitswide)] <- NA

diag(LowerHalfGmat9traitswide) <- NA

FullTableGMat <- UpperHalfCovmat9traitswide
FullTableGMat[lower.tri(FullTableGMat)] <- LowerHalfGmat9traitswide[lower.tri(LowerHalfGmat9traitswide)]

# rounding
FullTableGMat <- round(FullTableGMat, digits = 3)
# write out file
# write.csv(FullTableGMat, "outputs/FullGTable.csv")






# work with posteriors ####
# Extract sds for each trait for father and mother

## create new labels for treatments ####
newcroplabs <- c("Maize", "Soybean", 
                 "Tomato") 
names(newcroplabs) <- c("Maz", "Soya",  
                        "Tom") 

newfunglabs <- c("Beauveria", "Control", 
                 "Metarhizium") 
names(newfunglabs) <- c("BB", "CON",  
                        "MT") 
## wrangle data ####
summary(model.9traits.brm.matnocov.long)
Posteriors9Traits_lc <- as_tibble(as_draws_df(model.9traits.brm.matnocov.long))
names(Posteriors9Traits_lc)
Posteriors9Traits_lc %>% 
  rownames_to_column()

names(Posteriors9Traits_lc)
Posterior9Traits_sds_lc <- Posteriors9Traits_lc %>% 
  rownames_to_column() %>%
  select(rowname, sd_father__treatmentSoya_CON:
           sd_mother__Intercept) %>%
  pivot_longer(sd_father__treatmentSoya_CON:
                 sd_father__treatmentTom_MT,
               names_to = "traits",
               values_to = "SD") %>% 
  mutate(treat = traits) %>% 
  separate(traits, into = c("junk", "junk2", "crop",
                            "fungus")) %>% 
  select(-junk, -junk2) %>% 
  separate(crop,
           into = c("junk3","crop"),
           sep = 9) %>% 
  select(-junk3) %>% 
  separate(treat, into = c("junk4", "treat"), sep = 20) %>% 
  select(-junk4) 




# in following line compute heritability within treatments considering residual variance set to 1, 
# and square of iteration-level maternal variance 
## extract sample-level maternal variance  ####
Posterior9Traits_sds_lc <- Posterior9Traits_sds_lc %>% 
  mutate(herit = SD*SD / ((SD * SD) + 1 + (sd_mother__Intercept * sd_mother__Intercept)))


# extract 89%HDI for heritability
Posterior9Traits_sds_lc %>% 
  group_by(treat) %>% 
  summarise(HDI_herit = hdi(herit, ci = 0.89))

# Posterior9Traits_sds_lc %>% 
#   select(treat,herit) %>% 
#   group_by(treat) %>% 
#   summarise(HDI_MAP = map_estimate(herit))

map_estimate(model.9traits.brm.matnocov.long, effects = "random")

model.9traitslongRC <- model.9traits.brm.matnocov.long %>% 
  recover_types()

model.9traitslongRC %>%
  summarise_draws()




## heritability ridge plot ####
Posterior9Traits_sds_lc %>% 
  ggplot(aes(x = herit, 
             y = fct_rev(treat),
             # group = treat,
             fill = crop)) +
  geom_density_ridges(alpha = 0.3) +
  # facet_wrap(crop~fungus, nrow = 3, 
  #            labeller = labeller(fungus = newfunglabs)) +
  lims(x = c(0,1)) +
  labs(x = "Heritability",
       y = "") +
  scale_fill_manual(labels = c("Soybean", "Maize", "Tomato"),
                    values = c("darkgreen","gold",  "maroon")) +
  scale_y_discrete(labels=expression(italic(Metarhizium),italic(Beauveria),Control,italic(Metarhizium),italic(Beauveria),Control,italic(Metarhizium),italic(Beauveria),Control)) +
  theme(legend.position = "top",
        legend.title=element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 16))



ggsave("outputs/Heritability_ridges.png",
       height = 20,
       width = 12,
       units = "cm")





# Extract sds for each trait for father and mother
Gmat9traitslong <- Gmat9traitswide %>% 
  pivot_longer(treatmentMaz_BB:treatmentTom_MT,
             names_to = "Trait2",
             values_to = "Rg") %>% 
  separate(Trait1,
           into = c("crop1","fungus1")) %>% 
  separate(crop1,
           into = c("junk1","crop1"),
           sep = 9) %>% 
  separate(Trait2,
           into = c("crop2","fungus2")) %>% 
  separate(crop2,
           into = c("junk2","crop2"),
           sep = 9) %>% 
  select(-junk1,-junk2) %>% 
  distinct(Rg,
           .keep_all = TRUE) %>% 
  filter(Rg != 1) %>% 
  mutate(cropdiff = if_else(crop1 == crop2, 0, 1),
         fungdiff = if_else(fungus1 == fungus2, 0, 1),
         diff_dims = cropdiff+fungdiff,
         funguscontrol = if_else(fungus1 == "CON" | fungus2 == "CON", "control", "2 fungi")) %>% 
  mutate(diff_desc = if_else(diff_dims == 2, "both different",
                             if_else(cropdiff == 1, "crop different", "infection treatment different")))
# maybe plot differences in estimates for indiv correlations in mcmcglmm vs brms



# plot correlations 9-trait model long chains ####

Posteriors9Traits_lc <- as_tibble(as_draws_df(model.9traits.brm.matnocov.long))
names(Posteriors9Traits_lc)

Posterior9Traits_long_lc <- Posteriors9Traits_lc %>% 
  select(cor_father__treatmentSoya_CON__treatmentSoya_BB:
           cor_father__treatmentTom_BB__treatmentTom_MT) %>% 
  pivot_longer(cor_father__treatmentSoya_CON__treatmentSoya_BB:
                 cor_father__treatmentTom_BB__treatmentTom_MT,
               names_to = "traits",
               values_to = "Rg") %>% 
  mutate(treat = traits) %>% 
  separate(traits, into = c("junk", "junk2", "crop1",
                            "fungus1", "crop2", "fungus2")) %>% 
  select(-junk, -junk2) %>% 
  separate(crop1,
           into = c("junk3","crop1"),
           sep = 9) %>% 
  separate(crop2,
           into = c("junk4","crop2"),
           sep = 9) %>% 
  select(-junk3,-junk4) %>% 
  mutate(cropdiff = if_else(crop1 == crop2, 0, 1),
         fungdiff = if_else(fungus1 == fungus2, 0, 1),
         diff_dims = cropdiff+fungdiff,
         funguscontrol = if_else(fungus1 == "CON" | fungus2 == "CON", "control", "2 fungi")) %>% 
  mutate(diff_desc = if_else(diff_dims == 2, "both different",
                             if_else(cropdiff == 1, "crop different", "infection treatment different"))) %>% 
  mutate(treat1 = paste(crop1,fungus1, sep = "_"),
         treat2 = paste(crop2,fungus2, sep = "_"))


# change order of levels
levels(factor(Posterior9Traits_long_lc$treat1))
names(Posterior9Traits_long_lc)
Posterior9Traits_long_lc <- Posterior9Traits_long_lc %>% 
  mutate(crop1 = fct_relevel(crop1, "Soya", "Maz", "Tom"),
         fungus1 = fct_relevel(fungus1, "CON", "BB", "MT"),
         treat1 = fct_relevel(treat1, "Soya_CON", "Soya_BB",  "Soya_MT", 
                              "Maz_CON", "Maz_BB", "Maz_MT",
                              "Tom_CON", "Tom_BB" ),
         crop2 = fct_relevel(crop2, "Soya", "Maz", "Tom"),
         fungus2 = fct_relevel(fungus2, "CON", "BB", "MT"),
         diff_desc = fct_relevel(diff_desc, "crop different", "infection treatment different", "both different"),
         treat2 = fct_relevel(treat2, "Soya_BB",  "Soya_MT", 
                              "Maz_CON", "Maz_BB", "Maz_MT",
                              "Tom_CON", "Tom_BB","Tom_MT"),
         treatshort = paste(treat1,treat2),
         treatshort = fct_relevel(treatshort, 
                                  "Soya_CON Soya_BB",
                                  "Soya_CON Soya_MT",
                                  "Soya_CON Maz_CON",
                                  "Soya_CON Maz_BB",
                                  "Soya_CON Maz_MT",
                                  "Soya_CON Tom_CON",
                                  "Soya_CON Tom_BB",
                                  "Soya_CON Tom_MT",
                                  "Soya_BB Soya_MT",
                                  "Soya_BB Maz_CON",
                                  "Soya_BB Maz_BB",
                                  "Soya_BB Maz_MT",
                                  "Soya_BB Tom_CON",
                                  "Soya_BB Tom_BB", 
                                  "Soya_BB Tom_MT", 
                                  "Soya_MT Maz_CON",
                                  "Soya_MT Maz_BB",
                                  "Soya_MT Maz_MT",
                                  "Soya_MT Tom_CON" ,
                                  "Soya_MT Tom_BB" ,
                                  "Soya_MT Tom_MT" ,
                                  "Maz_CON Maz_BB",
                                  "Maz_CON Maz_MT",
                                  "Maz_CON Tom_CON" ,
                                  "Maz_CON Tom_BB" , 
                                  "Maz_CON Tom_MT" ,
                                  "Maz_BB Maz_MT" , 
                                  "Maz_BB Tom_CON" ,
                                  "Maz_BB Tom_BB" ,
                                  "Maz_BB Tom_MT" , 
                                  "Maz_MT Tom_CON" , 
                                  "Maz_MT Tom_BB" ,  
                                  "Maz_MT Tom_MT"   ,
                                  "Tom_CON Tom_BB",
                                  "Tom_CON Tom_MT",
                                  "Tom_BB Tom_MT"))


levels(factor(Posterior9Traits_long_lc$treatshort))




## create new labels for treatments ####

newcroplabs <- c("Maize", "Soybean", 
                 "Tomato") 
names(newcroplabs) <- c("Maz", "Soya",  
                        "Tom") 

newfunglabs <- c("Beauveria", "Control", 
                 "Metarhizium") 
names(newfunglabs) <- c("BB", "CON",  
                        "MT") 



## Genetic correlation ridges ####

facet_names <- as_labeller(c(`1`= "1D difference", 
                             `2`= "2D difference"))

names(Posterior9Traits_long_lc)
Posterior9Traits_long_lc %>% 
  ggplot(aes(x = Rg, 
             y = fct_rev(treatshort),
             group = treatshort,
             fill = factor(diff_desc))) +
  geom_density_ridges2(alpha = 0.3) +
  facet_wrap(~ diff_dims,
             labeller = facet_names)+
  geom_vline(xintercept = 0, 
             linetype = 1,
             size = 1,
             alpha = 0.4)+
  labs(x = expression(Genetic~correlation~(italic(r)[g])),
       y = "") +
  theme(legend.position = "top",
        legend.title=element_blank(),
        axis.text.x = element_text(size = 24),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill="white" ))


ggsave("outputs/RGsRidgeFacetDiffDesc.png",
       height = 30,
       width = 25,
       units = "cm",
       bg = "white")





Posterior9Traits_long_lc %>% 
  ggplot(aes(x = Rg, 
             y = fct_rev(treatshort),
             group = treatshort,
             fill = factor(diff_desc))) +
  geom_density_ridges2(alpha = 0.3) +
  geom_vline(xintercept = 0, 
             linetype = 1,
             size = 1,
             alpha = 0.4)

ggsave("outputs/RGsRidgeSingleDiffDesc.png")











Posterior9Traits_long_lc %>% 
  ggplot(aes(x = Rg, 
             y = ..density..,
             group = treat,
             fill = factor(cropdiff))) +
  geom_density(adjust = 1/2,
               alpha = 0.3,
               bw = 0.02) +
  facet_grid(treat1~treat2) +
  geom_vline(xintercept = 0, 
             linetype = 1,
             size = 1,
             alpha = 0.3)

ggsave("outputs/RGsMatrixLayout.png")







# reaction norm plots ####

## calc sire level survival by treat comb ####
df_qgen_sub_sirecols <- df_qgen_sub_df %>%
  group_by(father, treatment) %>% 
  summarise(Famsize = n(),
            Dead = sum(MORTday14, 
                       na.rm = TRUE),
            Alive = Famsize - Dead,
            prop_mort = mean(MORTday14, 
                             na.rm = TRUE)) %>% 
  ungroup()

# compute CIs for mortality
df_qgen_sire_cis <- binom.confint(
  x = df_qgen_sub_sirecols$Dead, 
  n = df_qgen_sub_sirecols$Famsize,
  methods = "exact") %>% 
  rename(PropMort = mean,
         FamSize2 = n)

df_qgen_sire_data <- bind_cols(df_qgen_sub_sirecols,
                               df_qgen_sire_cis) %>%
  select(-method, -x, -PropMort, -FamSize2)


## compute spearman correlations ####
names(df_qgen_sire_data)
df_qgen_sire_morts_wide <-  df_qgen_sire_data %>% 
  select(father, prop_mort, treatment) %>% 
  pivot_wider(names_from = treatment,
              values_from = prop_mort)

df_qgen_sire_morts_wide <- as.matrix(df_qgen_sire_morts_wide)

rcorr(df_qgen_sire_morts_wide[,2:10],
      type = "spearman")

# create subset and arrange by soyaBB
df_qgen_sire_data_soyBB <- df_qgen_sire_data %>%
  filter(treatment == "Soya_BB") %>% 
  arrange(prop_mort) 

df_qgen_sire_data$father <- factor(df_qgen_sire_data$father, levels = df_qgen_sire_data_soyBB$father[order(df_qgen_sire_data_soyBB$prop_mort)])

levels(df_qgen_sire_data_soyBB$father)

# colour fathers to highlight a few
df_qgen_sire_data <- df_qgen_sire_data %>% 
  mutate(father_colour = if_else(
    father == "S_288", "2",
    if_else(father == "S_250", "3",
            if_else(father == "S_1335", "4",
                    if_else(father == "S_214", "5",
                            "1")))
  ))

cbbPalette <- c("gray40", "#E69F00","#0072B2", "#CC79A7", "#009E73", "#F0E442", "#56B4E9", "#D55E00" )


colnames(df_qgen_sire_data)

## individual panels ####
SoyaBB_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Soya_BB") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  labs(x="Sire Identity",
       y="Beauveria",
       title = "Soybean")+
  lims(y = c(0,1)) +
  guides(colour="none") +
  scale_colour_manual(values=cbbPalette) +
  theme(axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y=element_text(face="italic"),
        # axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 0,
                                  size = 24),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )

SoyaMT_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Soya_MT") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  labs(x="Sire Identity",
       y="Metarhizium")+
  lims(y = c(0,1)) +
  scale_colour_manual(values=cbbPalette) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y=element_text(face="italic"),
        # axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )

SoyaCON_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Soya_CON") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  labs(x="Sire Identity",
       y="Control")+
  lims(y = c(0,1)) +  
  scale_colour_manual(values=cbbPalette) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        # axis.title.y=element_blank(),
        axis.title.x=element_blank(),        
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )



MaizeBB_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Maz_BB") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  labs(x="Sire Identity",
       y="Mortality",
       title = "Maize")+
  scale_colour_manual(values=cbbPalette) +
  lims(y = c(0,1)) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 0,
                                  size = 24),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )

MaizeMT_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Maz_MT") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  labs(x="Sire Identity",
       y="")+
  scale_colour_manual(values=cbbPalette) +
  lims(y = c(0,1)) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )

MaizeCON_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Maz_CON") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  labs(x="Sire Identity",
       y="")+
  scale_colour_manual(values=cbbPalette) +
  lims(y = c(0,1)) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )



TomBB_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Tom_BB") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  scale_colour_manual(values=cbbPalette) +
  labs(x="Sire Identity",
       y="Mortality",
       title = "Tomato")+
  lims(y = c(0,1)) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 0,
                                  size = 24),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )

TomMT_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Tom_MT") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  scale_colour_manual(values=cbbPalette) +
  labs(x="Sire Identity",
       y="")+
  lims(y = c(0,1)) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )

TomCON_SireMort <- df_qgen_sire_data %>%
  filter(treatment == "Tom_CON") %>% 
  ggplot(aes(x = father, y = prop_mort)) +
  geom_pointrange(aes(ymin = lower, 
                      ymax = upper,
                      colour = father_colour)) +
  scale_colour_manual(values=cbbPalette) +
  labs(x="Sire Identity",
       y="")+
  lims(y = c(0,1)) +
  guides(colour="none") +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24)
  )



## ggarrange 9 panels ####
rxnnorms9panel <- ggarrange(SoyaBB_SireMort, MaizeBB_SireMort,
                            TomBB_SireMort, SoyaMT_SireMort,
                            MaizeMT_SireMort, TomMT_SireMort,
                            SoyaCON_SireMort, MaizeCON_SireMort,
                            TomCON_SireMort,
                            labels = c("1.00",
                                       "-0.04",
                                       "0.01",
                                       "0.68",
                                       "-0.04",
                                       "0.01",
                                       "0.54",
                                       "-0.13",
                                       "-0.02"),
                            label.x = c(0.3,0.8,0.18,
                                        0.3,0.8,0.18,
                                        0.3,0.8,0.16),
                            label.y = c(0.83,0.16,0.16,
                                        0.95,0.16,0.16,
                                        0.95,0.16,0.16),
                            # hjust = c(-6,-2,-2,-4,-4,-4,-4,-4,-4),
                            # vjust = 3,
                            widths = c(1,0.92,0.92),
                            heights = c(1,0.9,0.9),
                            nrow = 3, 
                            ncol = 3
)


# Annotate the figure by adding a common labels
annotate_figure(rxnnorms9panel,
                top = text_grob("Host plant", 
                                # color = "red",
                                size = 24),
                bottom = text_grob("Sire Identity",
                                   # color = "blue",
                                   # hjust = 1, 
                                   # x = 1, 
                                   # face = "italic",
                                   size = 24),
                # left = text_grob("Pathogen treatment", 
                #                  # color = "green",
                #                  rot = 90),
                # right = "I'm done, thanks :-)!",
                # fig.lab = "Figure 1", fig.lab.face = "bold"
)

ggsave("outputs/NewReactionNorms3x3.png",
       height = 20,
       width = 30,
       units = "cm",
       bg = "white")




## ggarrange no corr numbers 9 panels ####
rxnnorms9panel <- ggarrange(SoyaBB_SireMort, MaizeBB_SireMort,
                            TomBB_SireMort, SoyaMT_SireMort,
                            MaizeMT_SireMort, TomMT_SireMort,
                            SoyaCON_SireMort, MaizeCON_SireMort,
                            TomCON_SireMort,
                            # labels = c("1.00",
                            #            "-0.04",
                            #            "0.01",
                            #            "0.68",
                            #            "-0.04",
                            #            "0.01",
                            #            "0.54",
                            #            "-0.13",
                            #            "-0.02"),
                            # label.x = c(0.3,0.8,0.18,
                            #             0.3,0.8,0.18,
                            #             0.3,0.8,0.16),
                            # label.y = c(0.83,0.16,0.16,
                            #             0.95,0.16,0.16,
                            #             0.95,0.16,0.16),
                            # hjust = c(-6,-2,-2,-4,-4,-4,-4,-4,-4),
                            # vjust = 3,
                            widths = c(1,0.92,0.92),
                            heights = c(1,0.9,0.9),
                            nrow = 3, 
                            ncol = 3
)


# Annotate the figure by adding a common labels
annotate_figure(rxnnorms9panel,
                top = text_grob("Host plant", 
                                # color = "red",
                                size = 24),
                bottom = text_grob("Sire Identity",
                                   # color = "blue",
                                   # hjust = 1, 
                                   # x = 1, 
                                   # face = "italic",
                                   size = 24),
                # left = text_grob("Pathogen treatment", 
                #                  # color = "green",
                #                  rot = 90),
                # right = "I'm done, thanks :-)!",
                # fig.lab = "Figure 1", fig.lab.face = "bold"
)

ggsave("outputs/NoCorrs_NewReactionNorms3x3.png",
       height = 20,
       width = 30,
       units = "cm",
       bg = "white")


# get summary stats on genetic corrs ####

VCVarray9traits_werr <- VarCorr(model.9traits.brm.matnocov.long,
                           probs = c(0.11,0.89))


names(VCVarray9traits_werr)

Q11Cors9traitswide <- melt(VCVarray9traits_werr) %>% 
  filter(L2 == "cor",
         Var2 == "Q11") %>% 
  rename(Q11 = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Q11,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)

Q89Cors9traitswide <- melt(VCVarray9traits_werr) %>% 
  filter(L2 == "cor",
         Var2 == "Q89") %>% 
  rename(Q89 = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Q89,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)



Gcor9traitswide <- melt(VCVarray9traits_werr) %>% 
  filter(L2 == "cor",
         Var2 == "Estimate") %>% 
  rename(Cor = Var3,
         Trait1 = Var1) %>% 
  pivot_wider(.,
              names_from = Cor,
              values_from = value) %>% 
  select(-Var2,-L2,-L1)
# 
# order_levels <- c("treatmentSoya_CON", 
#                   "treatmentSoya_BB",
#                   "treatmentSoya_MT",
#                   "treatmentMaz_CON",
#                   "treatmentMaz_BB",
#                   "treatmentMaz_MT",
#                   "treatmentTom_CON",
#                   "treatmentTom_BB",
#                   "treatmentTom_MT" )


Gcors9traitslong <- Gcor9traitswide %>% 
  pivot_longer(treatmentSoya_CON:treatmentTom_MT,
               names_to = "Trait2",
               values_to = "Rg") %>% 
  separate(Trait1,
           into = c("crop1","fungus1")) %>% 
  separate(crop1,
           into = c("junk1","crop1"),
           sep = 9) %>% 
  separate(Trait2,
           into = c("crop2","fungus2")) %>% 
  separate(crop2,
           into = c("junk2","crop2"),
           sep = 9) %>% 
  select(-junk1,-junk2) %>% 
  distinct(Rg,
           .keep_all = TRUE) %>% 
  filter(Rg != 1) %>% 
  mutate(cropdiff = if_else(crop1 == crop2, 0, 1),
         fungdiff = if_else(fungus1 == fungus2, 0, 1),
         diff_dims = cropdiff+fungdiff,
         funguscontrol = if_else(fungus1 == "CON" | fungus2 == "CON", "control", "2 fungi")) %>% 
  mutate(diff_desc = if_else(diff_dims == 2, "both different",
                             if_else(cropdiff == 1, "crop different", "infection treatment different")))


Gcors9traitslong

# relevel categories
Gcors9traitslong$diff_desc <- factor(Gcors9traitslong$diff_desc, 
                                    levels = c("infection treatment different","crop different","both different"))
levels(Gcors9traitslong$diff_desc)

## extract Q11 ####
Q11Gcors9traitslong <- Q11Cors9traitswide %>% 
  pivot_longer(treatmentSoya_CON:treatmentTom_MT,
               names_to = "Trait2",
               values_to = "Rg") %>% 
  separate(Trait1,
           into = c("crop1","fungus1")) %>% 
  separate(crop1,
           into = c("junk1","crop1"),
           sep = 9) %>% 
  separate(Trait2,
           into = c("crop2","fungus2")) %>% 
  separate(crop2,
           into = c("junk2","crop2"),
           sep = 9) %>% 
  select(-junk1,-junk2) %>% 
  distinct(Rg,
           .keep_all = TRUE) %>% 
  filter(Rg != 1) %>% 
  mutate(cropdiff = if_else(crop1 == crop2, 0, 1),
         fungdiff = if_else(fungus1 == fungus2, 0, 1),
         diff_dims = cropdiff+fungdiff,
         funguscontrol = if_else(fungus1 == "CON" | fungus2 == "CON", "control", "2 fungi")) %>% 
  mutate(diff_desc = if_else(diff_dims == 2, "both different",
                             if_else(cropdiff == 1, "crop different", "infection treatment different")))




## extract Q89 ####
Q89Cors9traitsLong <- Q89Cors9traitswide %>% 
  pivot_longer(treatmentSoya_CON:treatmentTom_MT,
               names_to = "Trait2",
               values_to = "Q89") %>% 
  separate(Trait1,
           into = c("crop1","fungus1")) %>% 
  separate(crop1,
           into = c("junk1","crop1"),
           sep = 9) %>% 
  separate(Trait2,
           into = c("crop2","fungus2")) %>% 
  separate(crop2,
           into = c("junk2","crop2"),
           sep = 9) %>% 
  select(-junk1,-junk2) %>% 
  distinct(Q89,
           .keep_all = TRUE) %>% 
  filter(Q89 != 1)

## joins frames ####

GmatwErrors <- left_join(Gcors9traitslong, Q11Gcors9traitslong) %>% 
  left_join(., Q89Cors9traitsLong)


GmatwErrors <- GmatwErrors %>% 
  mutate(diff_desc_2 = if_else(fungdiff == 0,
                               if_else(fungus1 == "CON",
                                       "uninfected",
                                       "infected"),
                               if_else(fungus1 == "CON" |
                                         fungus2 == "CON",
                                       "control-infected",
                                       "isolate genus")),
         fung_diff_2 = if_else(fungdiff == 0,
                               if_else(fungus1 == "CON",
                                       "both control",
                                       "same isolate"),
                               if_else(fungus1 == "CON" |
                                         fungus2 == "CON",
                                       "control-infected",
                                       "different isolate"))) 

names(GmatwErrors)

# compute bootstrapped CIs ####
GmatwErrors %>% 
  group_by(diff_desc) %>% 
  summarise(meanRg = mean_cl_boot(Rg, conf.int = .89))

# make a plot of genetic correlations by diff
GmatwBootsPlot <- GmatwErrors %>% 
  ggplot(aes(x = diff_desc, y = Rg)) +
  geom_violin(alpha = 0.2,
              fill = "gray20",
              trim = FALSE) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black", 
               linewidth = 1.2, size = 1.5,
               alpha = 0.5)+
  geom_jitter(height = 0,
              width = 0.04,
              alpha = 0.6,
              aes(color = diff_desc_2)) 
