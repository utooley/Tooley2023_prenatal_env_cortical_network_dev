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

# Fill in the birth income and education data -------------------------
#add in the full data for parental education
edu_income_prenatal <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Prenatal_Edu_LogIN.csv")
edu_income_prenatal$modid <- edu_income_prenatal$MODID
#join it to demo_data
demo_data <- left_join(demo_data, edu_income_prenatal, by= "modid")

demo_data$income_needs_demo_b_log <- demo_data  %>% dplyr::select(contains("LOG_")) %>% rowMeans(, na.rm=T)
demo_data$income_needs_demo_b_unlogged <- 10^(demo_data$income_needs_demo_b_log)

#birth
#0 – Less than high school 1 – Completed high school 2 – College graduate 3 – Advanced degree
#Y1-Y3
#1, Less than 12th grade | 2, High school degree/GED | 3, Some college/vocational school | 4, College degree (4 years) | 5, Graduate degree
demo_data$edu_birth <- case_match(demo_data$EDU, 0~1, 1~2, 2~4, 3~5) #recode so the two scales above match
#fill in demo_edu_birth for only the values that are missing
demo_data$demo_edu_b_filled_in <- ifelse(is.na(demo_data$demo_edu_b),demo_data$edu_birth,demo_data$demo_edu_b)

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
demo <- demo_data %>% dplyr::select(c("modid","child_birthweight","child_sex", "PMA_scan"),contains("age"), contains("disadv"),contains("income_needs_demo"), contains("demo_edu")) 
demo$child_age_y0_mri <- (demo$PMA_scan-38)/54 #change PMA to age in years at scan for later plotting
demo_long <- reshape(demo, direction="long", varying =list(c("child_age_y0_mri","child_age_y2_mri_fun","child_age_y3_mri_fun")), v.names = c("child_age_mri"), times=c("y0","y2","y3"), idvar = "modid", timevar = "timepoint") %>% select(.,modid, timepoint,child_age_mri,child_sex,disadv_prenatal, disadv_y1,disadv_y2,disadv_y3, income_needs_demo_b, income_needs_demo_b_unlogged, income_needs_demo_b_log, demo_edu_b,demo_edu_b_filled_in)

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
