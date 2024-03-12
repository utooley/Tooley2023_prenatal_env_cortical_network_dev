# Setup -------------------------------------------------------------------
library(raincloudplots)
library(tidyverse)
library(stats)
library(parallel)
library(lm.beta)
library(summarytools)
library(psych)
library(GGally)
library(cifti)
library(ciftiTools)
library(ggseg)
library(ggsegGordon)
library(visreg)
library(mgcv)
library(broom)
library(gratia)
library(cowplot)
library(gtsummary)
source("~/Box/projects/in_progress/Tooley2023_prenatal_env_cortical_network_dev/gamm_models/gamm_functions.R")

ciftiTools.setOption('wb_path', '/Applications/workbench/')
library(rgl) #to use ciftiTools graphics
gordon_parcel_matching <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.csv")
gordon_networks <- str_split_i(gordon_parcel_matching$label, "_", 3);gordon_networks <- as.factor(gordon_networks)

#color schemes
scale_fill_gordon <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c(c("#F434F4", "#962996", "#FF2929","#28FE29", "#F6F629", "#A967FD","gray","#FEFED5", 
                          "#282828", "#29D6D6","#F79029", "#299B9B", "#2929D0")), levels(gordon_networks)), 
    ...
  )
}
scale_color_gordon <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c(c("#F434F4", "#962996", "#FF2929","#28FE29", "#F6F629", "#A967FD","gray","#FEFED5", 
                          "#282828", "#29D6D6","#F79029", "#299B9B", "#2929D0")), levels(gordon_networks)), 
    ...
  )
}
# Load demographic and exclusion data ---------------------------------------------------------------
exclusions <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/subjects_all_inclusion_exclusion_01_25_2022.csv")
demographics <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Tooley_RC_20230118_updated08JUN23.csv")
demo_data_all <- left_join(demographics, exclusions, by="modid");dim(demo_data_all)

#exclude pre-term babies, MRI_injury exclusion, NICU exclusions, and birthweight exclusion, as well as the IRB exclusion
demo_data <- demo_data_all %>% filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0);dim(demo_data)
#remame PMA
demo_data$PMA_scan <- demo_data$mri_test_pma_scan_dob
#make sex a factor
demo_data$child_sex <- factor(demo_data$child_sex, labels = c("Male", "Female"))

# Load race data ----------------------------------------------------------
#1, White | 2, Black or African American | 3, American Indian or Alaska Native | 4, Asian Indian | 5, Chinese | 6, Filipino | 7, Japanese | 8, Korean | 9, Vietnamese | 10, Other Asian | 11, Native Hawaiian |
#12, Guamanian or Chamorro | 13, Samoan | 14, Other Pacific Islander | 15, Other
race_data <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Tooley_02_14_2024_race_hispanic.csv")
race_data$child_race <- ifelse(grepl( ",",race_data$child_race), "Multiracial", race_data$child_race)

race_data <- race_data %>% 
  mutate(
    #Replace with words for race
    child_race_word=case_match(child_race,
                         "1" ~ "White",
                         "2" ~ "Black",
                         "4" ~ "Asian Indian",
                         "5" ~ "Chinese",
                         "8" ~ "Korean",
                         "14"~ "Other Pacific Islander",
                         "15" ~ "Other",
                         .default=child_race))
                         
race_data$child_hispanic <- factor(race_data$child_hispanic, labels=c("Not Hispanic or Latino","Hispanic or Latino", "Unspecified"))
demo_data <- left_join(demo_data, race_data, by= "modid")

#demo_data %>% select(modid, GAWEEKS, child_race_word)
# Fill in the birth income and education data -------------------------
edu_income_prenatal <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Prenatal_Edu_LogIN.csv") #add in the full data for parental education
edu_income_prenatal$modid <- edu_income_prenatal$MODID
demo_data <- left_join(demo_data, edu_income_prenatal, by= "modid") #join it to demo_data

demo_data$income_needs_demo_b_log <- demo_data %>% dplyr::select(contains("LOG_")) %>% rowMeans(na.rm=T)
demo_data$income_needs_demo_b_unlogged <- 10^(demo_data$income_needs_demo_b_log)

#birth data
#0 – Less than high school 1 – Completed high school 2 – College graduate 3 – Advanced degree
#Y1-Y3 data
#1, Less than 12th grade | 2, High school degree/GED | 3, Some college/vocational school | 4, College degree (4 years) | 5, Graduate degree
demo_data$edu_birth <- case_match(demo_data$EDU, 0~1, 1~2, 2~4, 3~5) #recode so the two scales above match
demo_data$demo_edu_b_filled_in <- ifelse(is.na(demo_data$demo_edu_b),demo_data$edu_birth,demo_data$demo_edu_b) #fill in demo_edu_birth for only the values that are missing

median_prenatal_disadvantage <- median(demo_data$disadv_prenatal)

# Load additional SES indicators ------------------------------------------
SEM_birth_indicators <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/SEM_Birth_Ins_HEI.csv")
colnames(SEM_birth_indicators) <- c("MODID","insurance_status_birth", "HEI_birth");SEM_birth_indicators$MODID
SEM_birth_indicators$modid <- str_c("MOD",SEM_birth_indicators$MODID)
demo_data <- left_join(demo_data, SEM_birth_indicators, by="modid")

SEM_Y3_indicators <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/SEMY3_FS_JA14_ForRepository_06142023.csv")
colnames(SEM_Y3_indicators)[1:9] <- str_c(colnames(SEM_Y3_indicators[,1:9]),"_y3")
SEM_Y3_indicators$modid <- str_c("MOD",SEM_Y3_indicators$MODID)
demo_data <- left_join(demo_data, SEM_Y3_indicators, by="modid")

ADI_all_timepoints <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/ADI_ChildAge_23JUN23.csv")
cols <- ADI_all_timepoints %>% select(contains("neighborhood") | contains("age")) %>%  colnames()
ADI_all_timepoints <- ADI_all_timepoints %>% mutate(across(all_of(cols), as.numeric))
intersect(colnames(demo_data), colnames(ADI_all_timepoints))
demo_data <- demo_data %>% select(-contains("neighborhood")) 
demo_data <- left_join(demo_data,ADI_all_timepoints, by="modid")

# Load network data -------------------------------------------------------
y0_network_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/birth/fc/n319_nreg333_birth_within_between_gordon_withmodulpartcoef.csv") %>% rename(.,modid=Var1, system_segreg=system_segreg_prop, mean_within_sys=mean_within_sys_prop, mean_between_sys=mean_between_sys_prop)
y0_network_data$modid<- str_remove(y0_network_data$modid, "_V1_a");y0_network_data$modid<- str_remove(y0_network_data$modid, "_V1_b");
y0_network_data$timepoint <- "y0"
colnames(y0_network_data)

y2_network_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/y2/n114_y2_within_between_gordon_withmodulpartcoef.csv")%>% rename(.,modid=Var1, system_segreg=system_segreg_prop, mean_within_sys=mean_within_sys_prop, mean_between_sys=mean_between_sys_prop)
y2_network_data$timepoint <- "y2"
colnames(y2_network_data)

y3_network_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/y3/n88_y3_within_between_gordon_withmodulpartcoef.csv")%>% rename(.,modid=Var1,  system_segreg=system_segreg_prop, mean_within_sys=mean_within_sys_prop, mean_between_sys=mean_between_sys_prop)
y3_network_data$timepoint <- "y3"
y3_network_data$modid <- str_remove(y3_network_data$modid, "_Y3_a"); y3_network_data$modid <- str_remove(y3_network_data$modid, "_Y3_b")
colnames(y3_network_data)

all_network_data <- rbind(y0_network_data, y2_network_data)
all_network_data <- rbind(all_network_data, y3_network_data)
all_network_data$part_coef_avg <- (all_network_data$part_coef_neg+all_network_data$part_coef_pos)/2

# Load MRI age data -----------------------------------------------------------
#age is separated by modality for toddler timepoints
age_y2_mri <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Child_Age_MRI_Y2_Updated.csv")
age_y3_mri <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Child_Age_MRI_Y3.csv")
demo_data <- left_join(demo_data, age_y2_mri, by = "modid")
demo_data <- left_join(demo_data, age_y3_mri, by = "modid")

# Reshape data to long format --------------------------------------------
colnames(demo_data)
demo <- demo_data %>% dplyr::select(c("modid","child_birthweight","child_sex", "PMA_scan"),contains("age"), contains("disadv"),contains("income_needs_demo"), contains("demo_edu"), contains("child_birthweight"), contains("GAWEEKS")) 
demo$child_age_y0_mri <- (demo$PMA_scan-38)/54 #change PMA to age in years at scan for later plotting
demo_long <- reshape(demo, direction="long", varying =list(c("child_age_y0_mri","child_age_y2_mri_fun","child_age_y3_mri_fun")), v.names = c("child_age_mri"), times=c("y0","y2","y3"), idvar = "modid", timevar = "timepoint") %>% select(.,modid, timepoint,child_age_mri,child_sex,disadv_prenatal, disadv_y1,disadv_y2,disadv_y3, income_needs_demo_b, income_needs_demo_b_unlogged, income_needs_demo_b_log, demo_edu_b,demo_edu_b_filled_in, child_birthweight, GAWEEKS)

#merge network data in with demo data
network_demo_data_long <- left_join(all_network_data, demo_long, by=c("modid", "timepoint"))
colSums(table(network_demo_data_long$child_age_mri,network_demo_data_long$timepoint)) #how many healthy FT at eatch timepoint?

# Load motion data --------------------------------------------------------
y0_motion <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N319_avg_FD_retained_02_28_2023.csv")
y0_motion$timepoint <- "y0"
y0_motion$modid<- str_remove(y0_motion$MODID, "_V1_a");y0_motion$modid<- str_remove(y0_motion$modid, "_V1_b");
y2_motion <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N114_avg_FD_retained_03_01_2023.csv")
y2_motion$timepoint <- "y2"
y2_motion$modid<- str_remove(y2_motion$MODID, "_Y2_a");y2_motion$modid<- str_remove(y2_motion$modid, "_Y2_b");
y3_motion <- read_table("~/Box/Tooley 01_18_2023 eLABE Requested Data/N88_Y3_avg_FD_retained_06_07_2023.csv")
y3_motion$timepoint <- "y3"
y3_motion$modid <- y3_motion$MODID
all_motion_data <- rbind(y0_motion, y2_motion)
all_motion_data <- rbind(all_motion_data, y3_motion)

network_demo_data_long <- left_join(network_demo_data_long,all_motion_data, by=c("modid", "timepoint"))
colSums(table(network_demo_data_long$child_age_mri,network_demo_data_long$timepoint))

# Load pre-censored motion data --------------------------------------------------------
y0_motion_precensored <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N319_avg_FD_retained_no_censoring_06_16_2023.csv")
y0_motion_precensored$timepoint <- "y0"
y0_motion_precensored$modid<- str_remove(y0_motion_precensored$MODID, "_V1_a");y0_motion_precensored$modid<- str_remove(y0_motion_precensored$modid, "_V1_b");
y2_motion_precensored <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N114_avg_FD_retained_no_censoring_06_16_2023.csv")
y2_motion_precensored$timepoint <- "y2"
y2_motion_precensored$modid<- str_remove(y2_motion_precensored$MODID, "_Y2_a");y2_motion_precensored$modid<- str_remove(y2_motion_precensored$modid, "_Y2_b");
y3_motion_precensored <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N88_Y3_avg_FD_retained_no_censoring_06_16_2023.csv")
y3_motion_precensored$timepoint <- "y3"
y3_motion_precensored$modid <- y3_motion_precensored$MODID

all_motion_precensored_data <- rbind(y0_motion_precensored, y2_motion_precensored)
all_motion_precensored_data <- rbind(all_motion_precensored_data, y3_motion_precensored)

network_demo_data_long <- left_join(network_demo_data_long,all_motion_precensored_data, by=c("modid", "timepoint"))

#filter out anyone with < 5 min of uncensored data
network_demo_data_long$retained_frames_min <- (network_demo_data_long$retained_frames*0.8)/60
network_demo_data_long$total_frames_min <- (network_demo_data_long$total_frames.x*0.8)/60
network_demo_data_long$percent_retained <- network_demo_data_long$retained_frames/network_demo_data_long$total_frames.x
network_demo_data_long <- network_demo_data_long %>% filter(., retained_frames_min>5)
#check subject numbers
colSums(table(network_demo_data_long$child_age_mri,network_demo_data_long$timepoint))
#check avg amt of data
network_demo_data_long %>% filter(!is.na(child_sex) & retained_frames_min >10) %>% dim()
amt_of_data <-  network_demo_data_long %>% filter(!is.na(child_sex)) %>% select(child_sex,avg_FD_of_retained_frames, retained_frames, total_frames.x, retained_frames_min,total_frames_min,percent_retained)
amt_of_data %>% tbl_summary(statistic = list(
  all_continuous() ~ "{mean} ({min} - {max})",
  all_categorical() ~ "{n} / {N} ({p}%)"))

# Make Supplemental Table 1: Demographics  ----------------
#filter original birth demo data by those in the network analyses
birth_data <- left_join(y0_network_data, demo_data, by="modid")
birth_data$PMA_scan <- as.numeric(birth_data$PMA_scan)
birth_data$demo_edu_b <- factor(birth_data$demo_edu_b, labels = c("Less than 12th grade", "High school degree/GED",
  "Some college/vocational school","College degree (4 years)" ,"Graduate degree"))
birth_demo <-  birth_data %>% select(PMA_scan, child_race_word, child_hispanic, child_sex,GAWEEKS, child_birthweight, ADI, income_needs_demo_b, demo_edu_b,insurance_status_birth,HEI_birth, disadv_prenatal)
birth_demo$child_sex<- factor(birth_demo$child_sex, labels = c("Male", "Female"))
table <- birth_demo %>% tbl_summary(statistic = list(
    all_continuous() ~ "{mean} ({min} - {max})",
    all_categorical() ~ "{n} / {N} ({p}%)"
  ),
  digits = all_continuous() ~ 1,
  type = PMA_scan ~ "continuous",
  label = list(PMA_scan ~ "Age at scan (months)",
               GAWEEKS ~ "Gestational age (weeks)",
               child_sex ~ "Child sex",
               child_race_word ~ "Child race",
               child_hispanic ~ "Child ethnicity",
      child_birthweight ~ "Birthweight (g)",
      ADI ~ "Area Deprivation Index",
      income_needs_demo_b ~ "Income to Needs Ratio",
      demo_edu_b ~ "Highest level of parent education completed",
      disadv_prenatal ~ "Socioeconomic disadvantage factor score",
      insurance_status_birth ~ "Insurance status (private)",
      HEI_birth ~ "Healthy Eating Index"),
  missing = "no")
table %>% modify_header(label = "**Variable**", stat_0 = '**N = 261**') 
table %>% show_header_names()
theme_gtsummary_compact()
table %>%
  as_gt() %>%
  gt::gtsave(filename = "~/Box/projects/in_progress/within_between_network_longitudinal/Demographic_table.docx")

#Supplemental Table 2: correlation table of SES variables
library(datscience)
birth_demo$demo_edu_b_num <- as.numeric(birth_demo$demo_edu_b)
table <- apa_corrTable(birth_demo, rmDiag = F, summarystats = F, method = "spearman",
                       table_caption = c("Bivariate correlations between SES variables"))
table
table %>% save_flextable(., "~/Downloads/table.docx")

# Figure 1: GAMMs for age -------------------------------------------------------------
source("~/Box/projects/in_progress/Tooley2023_prenatal_env_cortical_network_dev/gamm_models/gamm_functions.R")
#system segregation
gam_age_segreg <- gamm(system_segreg ~ child_sex + s(child_age_mri, k = 4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_segreg$lme);summary(gam_age_segreg$gam)
visualize_model(gam_age_segreg , smooth_var = "child_age_mri", int_var= NULL, group_var = "modid",plabels = "System segregation", derivative=T,resid_spaghetti=T)
a <- resid_plot(gam_age_segreg ,"child_age_mri","modid", NULL, add.intercept = T)
a + xlab("Child age at MRI (years)") + ylab("System segregation (partial residuals)")

#modularity
gam_age_modul <- gamm(modul ~ child_sex + s(child_age_mri, k = 4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_modul$lme);summary(gam_age_modul$gam)
visualize_model(gam_age_modul , smooth_var = "child_age_mri", int_var= NULL, group_var = "modid", plabels = "Modularity", derivative = T, resid_spaghetti = F)
a <- resid_plot(gam_age_modul ,"child_age_mri","modid", NULL, add.intercept = T)
a + xlab("Child age at MRI (years)") + ylab("Modularity (partial residuals)")

#clustering coefficient 
gam_age_clust_coef <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k = 4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_clust_coef$lme);summary(gam_age_clust_coef$gam)
visualize_model(gam_age_clust_coef , smooth_var = "child_age_mri", int_var= NULL, group_var = "modid",plabels = "Average Clustering Coefficient", derivative=T, resid_spaghetti = F)
a <- resid_plot(gam_age_clust_coef ,"child_age_mri","modid", NULL, add.intercept = T)
a + xlab("Child age at MRI (years)") + ylab("Clustering coefficient (partial residuals)")+ggtitle("Clustering coefficient")

p.adjust(c(summary(gam_age_segreg$gam)$s.table[1,4],summary(gam_age_modul$gam)$s.table[1,4],summary(gam_age_clust_coef$gam)$s.table[1,4], summary(gam_age_part_coef$gam)$s.table[1,4]), method = "fdr")

# Fig 2: GAMMs for age x SES (prenatal disadvantage)-------------------------------------------------------------
## System segregation ##
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)

#then check with pboot if the interaction is better than no interaction
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight+ti(child_age_mri,disadv_prenatal, k=4), random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_segreg_ti);pb
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight + s(child_age_mri,by=disadv_prenatal, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_segreg_by);pb

resid_plot_int(modobj = gam_age_ses_segreg_by,term = "child_age_mri", int.term = "disadv_prenatal",add.intercept = T)
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+theme(legend.position = "none")

## Modularity ##
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames +retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)

#then check with pboot if the interaction is better than no interaction
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight+ s(child_age_mri,by=disadv_prenatal, k=4), random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_modul_by);pb
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight+ s(child_age_mri,by=disadv_prenatal, k=4), random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_modul_by);pb
a <- resid_plot_int(gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")+theme(legend.position = "none")

## clustering coefficient ## 
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_ti <- gamm(avgclustco_both ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_ti$lme);summary(gam_age_ses_clust_coef_ti$gam)
BIC(gam_age_ses_clust_coef_by$lme , gam_age_ses_clust_coef_ti$lme)

#then check with pboot if the interaction is better than no interaction
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight + s(child_age_mri,by=disadv_prenatal, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_clust_coef_by);pb
a <- resid_plot_int(gam_age_ses_clust_coef_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

#adjust p's
p.adjust(c(summary(gam_age_ses_segreg_by$gam)$s.table[2,4],summary(gam_age_ses_modul_by$gam)$s.table[2,4],summary(gam_age_ses_clust_coef_by$gam)$s.table[2,4], summary(gam_age_ses_part_coef_by$gam)$s.table[2,4]), method = "fdr")

# Prenatal SES moderates trajectories of cortical network segregation -------------------------------
#Correlations between segregation measures
anova(lmer(avgclustco_both ~ modul +(1|modid), data=network_demo_data_long))
cor.test(network_demo_data_long$avgclustco_both, network_demo_data_long$modul)
cor.test(network_demo_data_long$avgclustco_both, network_demo_data_long$system_segreg)
cor.test(network_demo_data_long$modul, network_demo_data_long$system_segreg)

#which measure of segregation accounts for the age x SES effect
## Check local segregation inclusion in other models ## 
#system segreg
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight +avgclustco_both , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight ++avgclustco_both+ s(child_age_mri,by=disadv_prenatal, k=4), random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_segreg_by);pb
#modul
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames +retained_frames + avgweight +avgclustco_both, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight+ +avgclustco_both + s(child_age_mri,by=disadv_prenatal, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_modul_by);pb

## check meso-scale or global segregation inclusion in local segregation model ##
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight + modul, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight + modul+ s(child_age_mri,by=disadv_prenatal, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_clust_coef_by);pb

gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight + system_segreg, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight + system_segreg+ s(child_age_mri,by=disadv_prenatal, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_clust_coef_by);pb

# Load data for regional local segregation age x SES effects -----------
#load regional clustco at birth
y0_clustco_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/birth/fc/n319_nreg333_birth_clust_coef_avg_nodewise.csv")
y0_clustco_data <- rename(y0_clustco_data, modid=Var1)
y0_clustco_data$modid<- str_remove(y0_clustco_data$modid, "_V1_a");y0_clustco_data$modid<- str_remove(y0_clustco_data$modid, "_V1_b");
y0_clustco_data$timepoint <- "y0"
y2_clustco_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/y2/n114_y2_clust_coef_avg_nodewise.csv") %>% rename(.,modid=Var1)
y2_clustco_data$timepoint <- "y2"
y3_clustco_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/y3/n80_y3_clust_coef_avg_nodewise.csv") %>% rename(.,modid=Var1)
y3_clustco_data$modid <- str_remove(y3_clustco_data$modid,"_Y3_a");y3_clustco_data$modid <- str_remove(y3_clustco_data$modid,"_Y3_b")
y3_clustco_data$timepoint <- "y3"

all_clustco_data <- rbind(y0_clustco_data, y2_clustco_data)
all_clustco_data <- rbind(all_clustco_data, y3_clustco_data)
network_demo_data_long_clustco <- left_join(network_demo_data_long,all_clustco_data, by=c("modid", "timepoint"))

# Fig 3a: Regional GAMMs for local segregation ----------------------------------------
## Using a modified version of Val's function ##
avgclustco_both.gordon.elabe <- network_demo_data_long_clustco #rename data to be what the function is looking for
gam.age.ses.clustco.gordon <- matrix(data=NA, nrow=333, ncol=4) #matrix to save gam.fit output to

for(row in c(1:nrow(gordon_parcel_matching))){ #for each region
  region <- str_c("avgclustco_all_", row)
  GAM.RESULTS <- gamm.fit.smooth.int(measure = "avgclustco_both", atlas = "gordon", dataset = "elabe", 
                                     region = region, smooth_var = "child_age_mri", int_var = "disadv_prenatal",
                                     id_var="modid", random_slope=FALSE,
                                     covariates = "child_sex + avg_FD_of_retained_frames + retained_frames + avgweight", 
                                     knots = 4, set_fx = TRUE, stats_only = T) #run the gam.fit.smooth function
  gam.age.ses.clustco.gordon[row,] <- GAM.RESULTS}
gam.age.ses.clustco.gordon.fx.TRUE <- as.data.frame(gam.age.ses.clustco.gordon)
colnames(gam.age.ses.clustco.gordon.fx.TRUE) <- c("label", #region name
                                                  "GAM.age.ses.Fvalue", #GAM F-value for the age smooth term
                                                  "GAM.age.ses.pvalue", #GAM p-value for the age smooth term
                                                  "GAM.age.ses.partialR2") #partial Rsq from age and age-null models
head(gam.age.ses.clustco.gordon)
cols = c(2:4)    
gam.age.ses.clustco.gordon.fx.TRUE[,cols] = apply(gam.age.ses.clustco.gordon.fx.TRUE[,cols], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(gam.age.ses.clustco.gordon.fx.TRUE, "~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_gamm_statistics_gordon_fx_T_no_neg.csv", row.names = F, quote = F)
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_gamm_statistics_gordon_fx_T_no_neg.csv")
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct

## Write out the age x SES effects to a pscalar so can look at it in workbench ##
values <- ifelse(gordon_networks=="None",0,gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.Fvalue)#remove 'None' network values from plotting
new.347.values <- rep(0, 347)#dim is 347, because of subcortical, so pad it with 0s
new.347.values[1:333] <- values
output_file=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_Fvalue_no_none_no_neg.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_Fvalue_no_none_no_neg.ptseries.nii")

command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

# Fig 3b: Plot the age x SES effect on local segregation by system-------------------------
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_gamm_statistics_gordon_fx_T_no_neg.csv")
head(gam.age.ses.clustco.gordon.fx.TRUE)
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")
gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks
kruskal.test(GAM.age.ses.Fvalue~network, data = gam.age.ses.clustco.gordon.fx.TRUE)

#make a boxplot
g <- ggplot(data = dplyr::filter(gam.age.ses.clustco.gordon.fx.TRUE, network !="None"), aes(x = reorder(network, -GAM.age.ses.Fvalue, median), y= GAM.age.ses.Fvalue, fill=network)) +
  geom_boxplot(width = .6, alpha = 0.8) +
  geom_point(aes(x = reorder(network, -GAM.age.ses.Fvalue, median), y = GAM.age.ses.Fvalue, color = network), alpha=0.3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_fill_gordon() +
  scale_color_gordon() +
  theme_bw() +
  ylab("Age x SES effects on local segregation (F)")+
  xlab("Networks")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=14)) #fix the offset
g

# Fig 3c: Spin test for age effects on local segregation to S-A axis --------------
source("~/Box/tools/rotate_parcellation/R/rotate.parcellation.R")
source("~/Box/tools/rotate_parcellation/R/perm.sphere.p.R")
library(matrixStats) #otherwise get rowMins error
library(scales)
load("~/Box/projects/in_progress/struct_funct_neonates/data/rotations_gordon333_volume_MNI_10000x.Rdata") #load presaved rotations
sa_axis_gordon <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/SensorimotorAssociation_Axis_Gordon333.csv")
print(perm.sphere.p(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.Fvalue,sa_axis_gordon$SA.Axis.ranks, rotations, "spearman")) #test the spin test
cor.test(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.Fvalue,sa_axis_gordon$SA.Axis.ranks, method = "spearman")

mydata <- data.frame(sa_axis_gordon,gam.age.ses.clustco.gordon.fx.TRUE)
mydata$sig_or_not <- factor(ifelse(mydata$GAM.age.ses.pvalue.fdr<0.05, 1,1))
ggplot(mydata, aes(x = SA.Axis.ranks, y = GAM.age.ses.Fvalue, fill = GAM.age.ses.Fvalue,color = sig_or_not)) + 
  geom_point(shape = 21, size = 3) +
  scale_fill_gradient2(low="darkblue",high = "darkorange", guide="colorbar",aesthetics = c( "fill"), midpoint = mean(mydata$GAM.age.ses.Fvalue), name=NULL)+
  scale_color_manual(values = c("black")) +
  labs(x="\nS-A Axis Ranks", y="Local segregation age-by-SES effect\n") +
  geom_smooth(method = 'lm', se = TRUE, fill = scales::alpha(c("gray70"),.7), col = "black", size = .75) +
  theme_classic()  
  #theme(legend.position = "none") 

# Fig 4a: Bayley over time and SES ----------------
library(data.table)
Bayley_scores <- demo_data %>% select(c("modid","child_birthweight","child_sex", "PMA_scan"),contains("age"), contains("disadv"),contains("income_needs_demo"), (contains("bayley_cog") | contains("bayley_mot")| contains("bayley_lang")) & contains("comp")) 
Bayley_scores$exclusions <- rep(F,dim(Bayley_scores)[1])
Bayley_long <- reshape(Bayley_scores, direction="long", varying =list(c("child_age_y1_assessment" ,"child_age_y2_assessment","child_age_y3_assessment"),c("bayley_cog_comp_y1" , "bayley_cog_comp_y2","bayley_cog_comp_y3"),c("bayley_mot_comp_y1","bayley_mot_comp_y2","bayley_mot_comp_y3"), c("bayley_lang_comp_y1","bayley_lang_comp_y2","bayley_lang_comp_y3")), v.names = c("child_age_assessment","bayley_cog_comp", "bayley_mot_comp","bayley_lang_comp"), times=c("y1","y2","y3"), idvar = "modid")

#Cognition
gam_bayley_cog_comp_by <- gamm(bayley_cog_comp ~ child_sex + s(child_age_assessment, k = 4) + s(child_age_assessment, k =4, by= disadv_prenatal), random = list(modid =~ 1), data=Bayley_long, method = "REML")
#ti() is singular
summary(gam_bayley_cog_comp_by$lme);summary(gam_bayley_cog_comp_by$gam)
visualize_model(gam_bayley_cog_comp_by, smooth_var = "child_age_assessment", int_var = "disadv_prenatal", group_var = "modid",plabels = "Bayley Cognitive Composite")
#Language
gam_bayley_lang_comp_by <- gamm(bayley_lang_comp ~ child_sex + s(child_age_assessment, k = 4) + s(child_age_assessment, k =4, by= disadv_prenatal), random = list(modid =~ 1), data=Bayley_long, method = "REML")
summary(gam_bayley_lang_comp_by$lme);summary(gam_bayley_lang_comp_by$gam)
gam_bayley_lang_comp_ti <- gamm(bayley_lang_comp ~ child_sex + ti(child_age_assessment, k = 4) +ti(disadv_prenatal, k=4)+ ti(child_age_assessment, disadv_prenatal, k =4), random = list(modid =~ 1), data=Bayley_long, method = "REML")
summary(gam_bayley_lang_comp_ti$lme);summary(gam_bayley_lang_comp_ti$gam)
BIC(gam_bayley_lang_comp_by$lme, gam_bayley_lang_comp_ti$lme)
visualize_model(gam_bayley_lang_comp_ti, smooth_var = "child_age_assessment", int_var = "disadv_prenatal", group_var = "modid",plabels = "Bayley Language Composite") 
p.adjust(c(summary(gam_bayley_cog_comp_by$gam)$s.table[2,4],summary(gam_bayley_lang_comp_ti$gam)$s.table[3,4]))
         
# Fig 4b: Cognitive and language outcomes predicted by local segregation, Y2 ----------------
y2_data_only <- left_join(y2_network_data, y2_motion, by = c("modid", "timepoint"))
y2_data_only <- left_join(y2_data_only, y2_motion_precensored, by = c("modid", "timepoint"))
y2_data_only <- left_join(y2_data_only, demo_data, by = c("modid"))
y2_data_only$part_coef_avg <- (y2_data_only$part_coef_neg+y2_data_only$part_coef_pos)/2
y2_data_only <- y2_data_only %>% filter(., bayley_dataquality_y2 >=2)#filter out bayley exclusions
describe(y2_data_only$bayley_cog_comp_y2)

#language
lm_bayley_lang_clustco <- lm(avgclustco_both ~ child_sex + child_age_y2_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_lang_comp_y2, data=y2_data_only)
summary(lm_bayley_lang_clustco)
lm.beta(lm_bayley_lang_clustco)
visreg(lm_bayley_lang_clustco)
#cognition
lm_bayley_cog_clustco <- lm(avgclustco_both ~ child_sex + child_age_y2_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_cog_comp_y2, data=y2_data_only)
summary(lm_bayley_cog_clustco)
lm.beta(lm_bayley_cog_clustco)
visreg(lm_bayley_cog_clustco)

# lm_bayley_lang_clustco <- lm(bayley_lang_comp_y2 ~ child_sex + avgclustco_both, data=y2_data_only)
# summary(lm_bayley_lang_clustco)
# visreg(lm_bayley_lang_clustco)
# 
# lm_bayley_cog_clustco <- lm(bayley_cog_comp_y2 ~ child_sex + avgclustco_both, data=y2_data_only)
# summary(lm_bayley_cog_clustco)
# visreg(lm_bayley_cog_clustco)

#adjusted ps
p.adjust(c(summary(lm_bayley_lang_clustco)$coefficients[7,4], summary(lm_bayley_cog_clustco)$coefficients[7,4]))

##when you control for SES, is this accounting for the effect or not?
#language
lm_bayley_lang_clustco <- lm(bayley_lang_comp_y2 ~ child_sex + disadv_prenatal + avgclustco_both, data=y2_data_only)
summary(lm_bayley_lang_clustco)
lm.beta(lm_bayley_lang_clustco)
visreg(lm_bayley_lang_clustco) #clust co significant when including disadvantage in the model, disadvantage also sig

# Is it Bayley expressive or receptive?
## Expressive
lm_bayley_elang_clustco <- lm(avgclustco_both ~ child_sex + child_age_y2_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_elang_scale_y2, data=y2_data_only)
lm.beta(lm_bayley_elang_clustco)
summary(lm_bayley_elang_clustco)
visreg(lm_bayley_elang_clustco, "bayley_elang_scale_y2")
visreg(lm_bayley_elang_clustco, "bayley_elang_scale_y2",xlim=c(1,13))

##when you control for SES, is this accounting for the effect or not?
lm_bayley_elang_clustco <- lm(bayley_elang_scale_y2 ~ child_sex + avgclustco_both + disadv_prenatal, data=y2_data_only)
summary(lm_bayley_elang_clustco)
lm.beta(lm_bayley_elang_clustco)
visreg(lm_bayley_elang_clustco) #clust co no longer significant when including disadvantage in the model, disadvantage also marginal
visreg(lm_bayley_elang_clustco, "avgclustco_both")

## Receptive
lm_bayley_rlang_clustco <- lm(avgclustco_both ~ child_sex + child_age_y2_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_rlang_scale_y2, data=y2_data_only)
lm.beta(lm_bayley_rlang_clustco)
summary(lm_bayley_rlang_clustco)
visreg(lm_bayley_rlang_clustco)
a <- visreg(lm_bayley_rlang_clustco, "bayley_rlang_scale_y2");plot(a)
plot(a, xlim=c(1,13))
p.adjust(c(summary(lm_bayley_rlang_clustco)$coefficients[7,4], summary(lm_bayley_elang_clustco)$coefficients[7,4]))

##when you control for SES, is this accounting for the effect or not?
lm_bayley_rlang_clustco <- lm(bayley_rlang_scale_y2 ~ child_sex + avgclustco_both + disadv_prenatal, data=y2_data_only)
summary(lm_bayley_rlang_clustco)
lm.beta(lm_bayley_rlang_clustco)
visreg(lm_bayley_rlang_clustco, "avgclustco_both")

# Cognitive and language outcomes predicted by local segregation, Y3 --------
y3_data_only <- left_join(y3_network_data, y3_motion, by = c("modid", "timepoint"))
y3_data_only <- left_join(y3_data_only, y3_motion_precensored, by = c("modid", "timepoint"))
y3_data_only <- left_join(y3_data_only, demo_data, by = c("modid"))
y3_data_only <- y3_data_only %>% filter(.,modid != "MOD1039") #take out the person with < 5 min of data remaining post-censoring
y3_data_only$part_coef_avg <- (y3_data_only$part_coef_neg+y3_data_only$part_coef_pos)/2
y3_data_only <- y3_data_only %>% filter(., bayley_dataquality_y3 >=2)#filter out bayley exclusions
describe(y3_data_only$bayley_cog_scale_y3)

#language
lm_bayley_lang_clustco <- lm(avgclustco_both ~ child_sex + child_age_y3_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_lang_comp_y3, data=y3_data_only)
summary(lm_bayley_lang_clustco)
lm.beta(lm_bayley_lang_clustco)
visreg(lm_bayley_lang_clustco)

#cognition
lm_bayley_cog_clustco <- lm(avgclustco_both ~ child_sex + child_age_y3_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_cog_comp_y3, data=y3_data_only)
summary(lm_bayley_cog_clustco)
lm.beta(lm_bayley_cog_clustco)
visreg(lm_bayley_cog_clustco)
p.adjust(c(summary(lm_bayley_lang_clustco)$coefficients[7,4], summary(lm_bayley_cog_clustco)$coefficients[7,4]))

# lm_bayley_lang_clustco <- lm(bayley_lang_comp_y3 ~ child_sex + avgclustco_both, data=y3_data_only)
# summary(lm_bayley_lang_clustco)
# visreg(lm_bayley_lang_clustco)
# lm_bayley_cog_clustco <- lm(bayley_cog_comp_y3 ~ child_sex + avgclustco_both, data=y3_data_only)
# summary(lm_bayley_cog_clustco)
# visreg(lm_bayley_cog_clustco)

##when you control for SES, is this accounting for the effect or not?
lm_bayley_lang_clustco <- lm(bayley_lang_comp_y3 ~ child_sex + disadv_prenatal + avgclustco_both, data=y3_data_only)
summary(lm_bayley_lang_clustco)
lm.beta(lm_bayley_lang_clustco)
visreg(lm_bayley_lang_clustco) #clust co no longer significant when including disadvantage in the model, disadvantage also marginal

# Is it Bayley expressive or receptive?
## Expressive
lm_bayley_elang_clustco <- lm(avgclustco_both ~ child_sex + child_age_y3_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_elang_scale_y3, data=y3_data_only)
lm.beta(lm_bayley_elang_clustco)
summary(lm_bayley_elang_clustco)
visreg(lm_bayley_elang_clustco, "bayley_elang_scale_y3")
visreg(lm_bayley_elang_clustco, "bayley_elang_scale_y3",xlim=c(1,13))

##when you control for SES, is this accounting for the effect or not?
lm_bayley_elang_clustco <- lm(bayley_elang_scale_y3 ~ child_sex + avgclustco_both + disadv_prenatal, data=y3_data_only)
summary(lm_bayley_elang_clustco)
lm.beta(lm_bayley_elang_clustco)
visreg(lm_bayley_elang_clustco) #clust co no longer significant when including disadvantage in the model, disadvantage also marginal
visreg(lm_bayley_elang_clustco, "avgclustco_both")

## Receptive
lm_bayley_rlang_clustco <- lm(avgclustco_both ~ child_sex + child_age_y3_mri_fun + avgweight + avg_FD_of_retained_frames + retained_frames+ bayley_rlang_scale_y3, data=y3_data_only)
lm.beta(lm_bayley_rlang_clustco)
summary(lm_bayley_rlang_clustco)
visreg(lm_bayley_rlang_clustco)
a <- visreg(lm_bayley_rlang_clustco, "bayley_rlang_scale_y3");plot(a)
plot(a, xlim=c(1,13))
p.adjust(c(summary(lm_bayley_rlang_clustco)$coefficients[7,4], summary(lm_bayley_elang_clustco)$coefficients[7,4]))

##when you control for SES, is this accounting for the effect or not?
lm_bayley_rlang_clustco <- lm(bayley_rlang_scale_y3 ~ child_sex + avgclustco_both + disadv_prenatal, data=y3_data_only)
summary(lm_bayley_rlang_clustco)
lm.beta(lm_bayley_rlang_clustco)
visreg(lm_bayley_rlang_clustco, "avgclustco_both")

# Figure 5a: Age x SES with composite SES (maternal ed + INR) -----------------------------
#make SES composite from the demo variables in SEM
network_demo_data_long$ses_composite_b_sem <- network_demo_data_long %>% select(c("demo_edu_b_filled_in","income_needs_demo_b_unlogged")) %>% rowMeans(na.rm=T)
describe(network_demo_data_long$ses_composite_b_sem)

## System segregation ##
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=ses_composite_b_sem, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
#then check with pboot if the interaction is better than no interaction
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight + s(child_age_mri,by=ses_composite_b_sem, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_segreg_by);pb
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "ses_composite_b_sem", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+theme(legend.position = "none")

## Modularity ##
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=ses_composite_b_sem, k=4) + avg_FD_of_retained_frames +retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
#then check with pboot if the interaction is better than no interaction
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight+ s(child_age_mri,by=ses_composite_b_sem, k=4), random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_modul_by);pb
a <- resid_plot_int(gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "ses_composite_b_sem", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")+theme(legend.position = "none")

## clustering coefficient ## 
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=ses_composite_b_sem, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
#then check with pboot if the interaction is better than no interaction
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight +modul+ s(child_age_mri,by=ses_composite_b_sem, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_clust_coef_by);pb
a <- resid_plot_int(gam_age_ses_clust_coef_by, "child_age_mri", "modid", int.term = "ses_composite_b_sem", add.intercept=T,high_color="#FAA820", low_color="#4D8AC8")
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## participation coefficient ##
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=ses_composite_b_sem, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight + s(child_age_mri,by=ses_composite_b_sem, k=4) , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
pb <- pboot(gam_age_ses_part_coef_by);pb
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "ses_composite_b_sem", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Participation coefficient (residuals)")+theme(legend.position = "none")

#adjust p's
format(p.adjust(c(summary(gam_age_ses_segreg_by$gam)$s.table[2,4],summary(gam_age_ses_modul_by$gam)$s.table[2,4],summary(gam_age_ses_clust_coef_by$gam)$s.table[2,4], summary(gam_age_ses_part_coef_by$gam)$s.table[2,4]), method = "fdr"), scientific= F)

# Fig 5b: Run a GAMM for clustco age x SES with mat ed + INR ----------------------------------------
ciftiTools.setOption('wb_path', '/Applications/workbench/')
library(rgl) #to use ciftiTools graphics
gordon_parcel_matching <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.csv")
gordon_networks <- str_split_i(gordon_parcel_matching$label, "_", 3);gordon_networks <- as.factor(gordon_networks)

#load regional clustco at birth
y0_clustco_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/birth/fc/n319_nreg333_birth_clust_coef_avg_nodewise.csv")
y0_clustco_data <- rename(y0_clustco_data, modid=Var1)
y0_clustco_data$modid<- str_remove(y0_clustco_data$modid, "_V1_a");y0_clustco_data$modid<- str_remove(y0_clustco_data$modid, "_V1_b");
y0_clustco_data$timepoint <- "y0"
y2_clustco_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/y2/n114_y2_clust_coef_avg_nodewise.csv") %>% rename(.,modid=Var1)
y2_clustco_data$timepoint <- "y2"
y3_clustco_data <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/y3/n88_y3_clust_coef_avg_nodewise.csv") %>% rename(.,modid=Var1)
y3_clustco_data$modid <- str_remove(y3_clustco_data$modid,"_Y3_a");y3_clustco_data$modid <- str_remove(y3_clustco_data$modid,"_Y3_b")
y3_clustco_data$timepoint <- "y3"

all_clustco_data <- rbind(y0_clustco_data, y2_clustco_data)
all_clustco_data <- rbind(all_clustco_data, y3_clustco_data)

#join clustco data to existing network data
network_demo_data_long_clustco <- left_join(network_demo_data_long,all_clustco_data, by=c("modid", "timepoint"))
avgclustco_both.gordon.elabe <- network_demo_data_long_clustco #rename data to be what the function is looking for
gam.age.ses.clustco.gordon <- matrix(data=NA, nrow=333, ncol=4) #matrix to save gam.fit output to

for(row in c(1:nrow(gordon_parcel_matching))){ #for each region
  region <- str_c("avgclustco_all_", row)
  GAM.RESULTS <- gamm.fit.smooth.int(measure = "avgclustco_both", atlas = "gordon", dataset = "elabe", 
                                     region = region, smooth_var = "child_age_mri", int_var = "ses_composite_b_sem",
                                     id_var="modid", random_slope=FALSE,
                                     covariates = "child_sex + avg_FD_of_retained_frames + retained_frames + avgweight", 
                                     knots = 4, set_fx = TRUE, stats_only = T) #run the gam.fit.smooth function
  gam.age.ses.clustco.gordon[row,] <- GAM.RESULTS}
gam.age.ses.clustco.gordon.fx.TRUE <- as.data.frame(gam.age.ses.clustco.gordon)
colnames(gam.age.ses.clustco.gordon.fx.TRUE) <- c("label", #region name
                                                  "GAM.age.ses.Fvalue", #GAM F-value for the age smooth term
                                                  "GAM.age.ses.pvalue", #GAM p-value for the age smooth term
                                                  "GAM.age.ses.partialR2") #partial Rsq from age and age-null models
head(gam.age.ses.clustco.gordon)
cols = c(2:4)    
gam.age.ses.clustco.gordon.fx.TRUE[,cols] = apply(gam.age.ses.clustco.gordon.fx.TRUE[,cols], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(gam.age.ses.clustco.gordon.fx.TRUE, "~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_SES_mat_ed_INR_gamm_statistics_gordon_fx_T_no_neg.csv", row.names = F, quote = F)
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_SES_mat_ed_INR_gamm_statistics_gordon_fx_T_no_neg.csv")
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct

## Write out the age x SES effects to a pscalar so can look at it in workbench ##
values <- ifelse(gordon_networks=="None",0,gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr)#remove 'None' network values
new.347.values <- rep(0, 347)#dim is 347, because of subcortical, so pad it with 0s
new.347.values[1:333] <- values
output_file=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_mat_ed_inr_pval_no_none_no_neg.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_mat_ed_inr_pval_no_none_no_neg.ptseries.nii")
command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

# Fig 5c: Plot the age x SES effect on clustering coefficient by system -------------------------
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_SES_mat_ed_INR_gamm_statistics_gordon_fx_T_no_neg.csv")
head(gam.age.ses.clustco.gordon.fx.TRUE)
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")
gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks

data <- gam.age.ses.clustco.gordon.fx.TRUE
kruskal.test(GAM.age.ses.Fvalue~network, data = data)
## Plot 
g <- ggplot(data = filter(data, network !="None"), aes(x = reorder(network, -GAM.age.ses.Fvalue, median), y= GAM.age.ses.Fvalue, fill=network)) +
  geom_boxplot(width = .6, alpha = 0.8) +
  geom_point(aes(x = reorder(network, -GAM.age.ses.Fvalue, median), y = GAM.age.ses.Fvalue, color = network), alpha=0.3, width = 0.2) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_fill_gordon() +
  scale_color_gordon() +
  theme_bw() +
  ylab("Age x SES effects on local segregation (F)")+
  xlab("Networks")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=14)) #fix the offset
g

# Fig 6: Network integration: participation coefficient -------------------
#age effect
gam_age_part_coef <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k = 4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_part_coef$lme);summary(gam_age_part_coef$gam)
visualize_model(gam_age_part_coef, smooth_var = "child_age_mri", int_var= NULL, group_var = "modid",plabels = "Average Participation Coefficient", derivative = T)
a <- resid_plot(gam_age_part_coef ,"child_age_mri","modid", NULL, add.intercept = T)
a + xlab("Child age at MRI (years)") + ylab("Participation coefficient (partial residuals)")+ggtitle("Participation coefficient")

#age x SES effect
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme) #neither are significant
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Participation coefficient (residuals)")+theme(legend.position = "none")

