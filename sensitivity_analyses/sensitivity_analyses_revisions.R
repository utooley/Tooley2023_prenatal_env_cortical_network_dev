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

demo_data$income_needs_demo_b_log <- demo_data  %>% dplyr::select(contains("LOG_")) %>% rowMeans(na.rm=T)
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

#check subject numbers for each timepoint
#Y0
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y0") %in% timepoint) & !any(c("y2","y3") %in% timepoint)) %>%  distinct(modid)
#Y2
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y2") %in% timepoint) & !any(c("y0","y3") %in% timepoint)) %>%  distinct(modid)
#Y3
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y3") %in% timepoint) & !any(c("y0","y2") %in% timepoint)) %>%  distinct(modid)
#Y0-Y2
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y0", "y2") %in% timepoint) & !any(c("y3") %in% timepoint)) %>%  distinct(modid)
#Y0-Y3
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y0", "y3") %in% timepoint) & !any(c("y2") %in% timepoint)) %>%  distinct(modid)
#Y0-Y3
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y2", "y3") %in% timepoint) & !any(c("y0") %in% timepoint)) %>%  distinct(modid)
#Y0-Y2-Y3
network_demo_data_long %>%
  filter(!is.na(child_age_mri)) %>%
  group_by(modid) %>%
  filter(all(c("y0","y2", "y3") %in% timepoint)) %>%  distinct(modid)

# Supplemental Figure 1: 2 or more timepoints of data/participant -----------------------------
modids_with_2_timepoints_more <- network_demo_data_long%>% group_by(modid) %>% filter(n_distinct(timepoint)>1) %>%  filter(!is.na(child_age_mri))%>% pull(modid)
modids_with_2_timepoints_more <- unique(modids_with_2_timepoints_more); length(modids_with_2_timepoints_more);length(unique(network_demo_data_long$modid)) #n=132 vs. original n=346 unique ids
network_demo_data_long_2_timepoints <- network_demo_data_long %>% filter(., modid %in% modids_with_2_timepoints_more)

#Age x SES GAMMs
## System segregation ##
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Modularity ##
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames +retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Clustering coefficient ## 
# gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
# summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_ti <- gamm(avgclustco_both~ child_sex + ti(child_age_mri, k=4)+  ti(child_age_mri,disadv_prenatal, k=4)+ti(disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_clust_coef_ti$lme);summary(gam_age_ses_clust_coef_ti$gam)
BIC(gam_age_ses_clust_coef_by$lme , gam_age_ses_clust_coef_ti$lme)# cannot use by= with this few datapoints
a <- resid_plot_int(gam_age_ses_clust_coef_ti, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Participation coefficient
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg ~ child_sex + ti(child_age_mri, k=4)+  ti(child_age_mri,disadv_prenatal, k=4)+ti(disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long_2_timepoints, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme)

#adjust p's
p.adjust(c(summary(gam_age_ses_segreg_by$gam)$s.table[2,4],summary(gam_age_ses_modul_by$gam)$s.table[2,4],summary(gam_age_ses_clust_coef_ti$gam)$s.table[2,4], summary(gam_age_ses_part_coef_by$gam)$s.table[2,4]), method = "fdr")

## Regional age x SES GAMMs
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
network_demo_data_long_2_timepoints_clustco <- left_join(network_demo_data_long_2_timepoints,all_clustco_data, by=c("modid", "timepoint"))

##Regional GAMMs for local segregation ##
avgclustco_both.gordon.elabe <- network_demo_data_long_2_timepoints_clustco #rename data to be what the function is looking for
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
write.csv(gam.age.ses.clustco.gordon.fx.TRUE, "~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_2_timepoints_gamm_statistics_gordon_fx_T_no_neg.csv", row.names = F, quote = F)
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_2_timepoints_gamm_statistics_gordon_fx_T_no_neg.csv")
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct

## Write out the age x SES effects to a pscalar so can look at it in workbench ##
values <- ifelse(gordon_networks=="None",0,gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr)#remove 'None' network values from plotting
new.347.values <- rep(0, 347)#dim is 347, because of subcortical, so pad it with 0s
new.347.values[1:333] <- values
output_file=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_2_timepoints_pvalue_no_none_no_neg.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_2_timepoints_pvalue_no_none_no_neg.ptseries.nii")

command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

#make a boxplot
gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks
kruskal.test(GAM.age.ses.Fvalue~network, data = gam.age.ses.clustco.gordon.fx.TRUE)
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

# Supplemental Figure 2: Using precensoring motion -----------------------------
#Age x SES GAMMs
## System segregation ##
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_all_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_all_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+theme(legend.position = "none")

## Modularity ##
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_all_frames +retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_all_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")+theme(legend.position = "none")

## clustering coefficient ## 
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_all_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_ti <- gamm(avgclustco_both ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_all_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_ti$lme);summary(gam_age_ses_clust_coef_ti$gam)
BIC(gam_age_ses_clust_coef_by$lme , gam_age_ses_clust_coef_ti$lme)

## Participation coefficient
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_all_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg  ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_prenatal, k=4) + ti(child_age_mri,disadv_prenatal, k=4)+ avg_FD_of_all_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme)

#adjust p's
p.adjust(c(summary(gam_age_ses_segreg_by$gam)$s.table[2,4],summary(gam_age_ses_modul_by$gam)$s.table[2,4],summary(gam_age_ses_clust_coef_by$gam)$s.table[2,4], summary(gam_age_ses_part_coef_by$gam)$s.table[2,4]), method = "fdr")

## Regional age x SES GAMMs
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

##Regional GAMMs for local segregation ##
avgclustco_both.gordon.elabe <- network_demo_data_long_clustco #rename data to be what the function is looking for
gam.age.ses.clustco.gordon <- matrix(data=NA, nrow=333, ncol=4) #matrix to save gam.fit output to

for(row in c(1:nrow(gordon_parcel_matching))){ #for each region
  region <- str_c("avgclustco_all_", row)
  GAM.RESULTS <- gamm.fit.smooth.int(measure = "avgclustco_both", atlas = "gordon", dataset = "elabe", 
                                     region = region, smooth_var = "child_age_mri", int_var = "disadv_prenatal",
                                     id_var="modid", random_slope=FALSE,
                                     covariates = "child_sex + avg_FD_of_all_frames + retained_frames + avgweight", 
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
write.csv(gam.age.ses.clustco.gordon.fx.TRUE, "~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_precensoring_motion_gamm_statistics_gordon_fx_T_no_neg.csv", row.names = F, quote = F)
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_precensoring_motion_gamm_statistics_gordon_fx_T_no_neg.csv")
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct

## Write out the age x SES effects to a pscalar so can look at it in workbench ##
values <- ifelse(gordon_networks=="None",0,gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr)#remove 'None' network values from plotting
new.347.values <- rep(0, 347)#dim is 347, because of subcortical, so pad it with 0s
new.347.values[1:333] <- values
output_file=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_precensoring_motion_pvalue_no_none_no_neg.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_sage_ses_precensoring_motion_pvalue_no_none_no_neg.ptseries.nii")

command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

#make a boxplot
gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks
kruskal.test(GAM.age.ses.Fvalue~network, data = gam.age.ses.clustco.gordon.fx.TRUE)
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

# Supplemental Figure 3a-c: Change in SES ------------------------------------
#make long data that includes all SES metrics
demo <- demo_data %>% select(c("modid","child_birthweight","child_sex", "PMA_scan"),contains("age"), contains("disadv"),contains("income_needs_demo"), contains("neighborhood"),contains("demo_edu_"), contains("bayley")) 
demo$child_age_y0_mri <- (demo$PMA_scan-38)/54
demo$exclusions <- rep(F,dim(demo)[1])
demo$bayley_cog_scale_y0 <- NA;demo$bayley_lang_scale_y0 <- NA;demo$bayley_mot_scale_y0 <- NA #make NA Bayley data for birth timepoint so no complaints reshaping
demo_long <- reshape(demo, direction="long", varying =list(c("child_age_y0_mri","child_age_y1_assessment","child_age_y2_mri_fun","child_age_y3_mri_fun"), 
                                                                  c("income_needs_demo_b","income_needs_demo_y1", "income_needs_demo_y2", "income_needs_demo_y3"),
                                                                  c("demo_edu_b","demo_edu_y1", "demo_edu_y2", "demo_edu_y3"),
                                                                  c("neighborhood_natlcentile_b","neighborhood_natlcentile_y1","neighborhood_natlcentile_y2","neighborhood_natlcentile_y3"),
                                                                  c("disadv_prenatal","disadv_y1","disadv_y2","disadv_y3"),
                                                                  c( "bayley_cog_scale_y0","bayley_cog_scale_y1","bayley_cog_scale_y2", "bayley_cog_scale_y3"),
                                                                  c("bayley_mot_scale_y0","bayley_mot_scale_y1","bayley_mot_scale_y2", "bayley_mot_scale_y3"), 
                                                                  c("bayley_lang_scale_y0","bayley_lang_scale_y1","bayley_lang_scale_y2", "bayley_lang_scale_y3")),
                            v.names = c("child_age_mri", "income_needs_demo","demo_edu", "ADI","disadv","bayley_cog_scale", "bayley_mot_scale","bayley_lang_scale"), times=c("y0","y1","y2","y3"), idvar = "modid", timevar = "timepoint") %>% 
  select(.,modid, timepoint,child_age_mri, income_needs_demo,demo_edu,ADI,disadv, bayley_mot_scale, bayley_cog_scale, bayley_lang_scale, child_sex)
demo_network_data_long <- left_join(demo_long,all_network_data, by=c("modid", "timepoint"))

#gam for parental education over time
gam_edu_visits <- gamm(demo_edu ~ s(child_age_mri, k = 4), random = list(modid =~ 1, modid=~0+child_age_mri), data=demo_network_data_long, method = "REML")
summary(gam_edu_visits$lme);summary(gam_edu_visits$gam)
gam_edu_visits$gam$data <- demo_network_data_long
visualize_model(gam_edu_visits$gam, smooth_var = "child_age_mri", int_var = NULL, group_var = "modid", resid_spaghetti = F)+xlab("Child age (years)")+ylab("Parental Education")

#gam for INR over time
gam_inr_visits <- gamm(income_needs_demo ~ s(child_age_mri, k = 4), random = list(modid =~ 1, modid=~0+child_age_mri), data=demo_network_data_long, method = "REML")
summary(gam_inr_visits$lme);summary(gam_inr_visits$gam)
gam_inr_visits$gam$data <- demo_network_data_long
visualize_model(gam_inr_visits$gam, smooth_var = "child_age_mri", int_var = NULL, group_var = "modid", resid_spaghetti = F)+xlab("Child age (years)")+ylab("Income to Needs Ratio")

#gam for ADI over time
gam_ADI_visits <- gamm(ADI ~ s(child_age_mri, k = 4), random = list(modid =~ 1, modid=~0+child_age_mri), data=demo_network_data_long, method = "REML")
summary(gam_ADI_visits$lme);summary(gam_ADI_visits$gam)
gam_ADI_visits$gam$data <- demo_network_data_long
visualize_model(gam_ADI_visits$gam, smooth_var = "child_age_mri", int_var = NULL, group_var = "modid", resid_spaghetti = F)+xlab("Child age (years)")+ylab("ADI")

# Supplemental Figure 3d-f: Disadvantage at Y1, Y2, and Y3 ------------------------------------
## System segregation ## 
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y1, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y1, k=4) + ti(child_age_mri,disadv_y1, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)

gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y2, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y2, k=4) + ti(child_age_mri,disadv_y2, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme) #ti is better

gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y3, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y3, k=4) + ti(child_age_mri,disadv_y3, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)

## Modularity ## 
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y1, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y1, k=4) + ti(child_age_mri,disadv_y1, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)

gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y2, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y2, k=4) + ti(child_age_mri,disadv_y2, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme) #ti is better

gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y3, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y3, k=4) + ti(child_age_mri,disadv_y3, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)

## Clustering coefficient ## 
## Y1 disadv ## 
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y1, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_ti <- gamm(avgclustco_both ~ child_sex + ti(child_age_mri, k=4)+ ti(child_age_mri,disadv_y1, k=4) + ti(disadv_y1, k=4) +avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_ti$lme);summary(gam_age_ses_clust_coef_ti$gam)
BIC(gam_age_ses_clust_coef_by$lme , gam_age_ses_clust_coef_ti$lme) #ti is better
a <- resid_plot_int(gam_age_ses_clust_coef_ti, "child_age_mri", "modid", int.term = "disadv_y1", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Y2 disadv ## 
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y2, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_ti <- gamm(avgclustco_both ~ child_sex + ti(child_age_mri, k=4)+ ti(child_age_mri,disadv_y2, k=4) + ti(disadv_y2, k=4) +avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_ti$lme);summary(gam_age_ses_clust_coef_ti$gam)
BIC(gam_age_ses_clust_coef_by$lme , gam_age_ses_clust_coef_ti$lme)
a <- resid_plot_int(gam_age_ses_clust_coef_by, "child_age_mri", "modid", int.term = "disadv_y2", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Y3 disadv ## 
gam_age_ses_clust_coef_by <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y3, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_by$lme);summary(gam_age_ses_clust_coef_by$gam)
gam_age_ses_clust_coef_ti <- gamm(avgclustco_both ~ child_sex + ti(child_age_mri, k=4)+ ti(child_age_mri,disadv_y3, k=4) + ti(disadv_y3, k=4) +avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clust_coef_ti$lme);summary(gam_age_ses_clust_coef_ti$gam)
BIC(gam_age_ses_clust_coef_by$lme , gam_age_ses_clust_coef_ti$lme)
a <- resid_plot_int(gam_age_ses_clust_coef_by, "child_age_mri", "modid", int.term = "disadv_y3", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Participation coefficient## 
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y1, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y1, k=4) + ti(child_age_mri,disadv_y1, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme)
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "disadv_y1", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y2, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y2, k=4) + ti(child_age_mri,disadv_y2, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme) #ti is better
a <- resid_plot_int(gam_age_ses_part_coef_ti, "child_age_mri", "modid", int.term = "disadv_y2", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_y3, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg ~ child_sex + ti(child_age_mri, k=4)+ ti(disadv_y3, k=4) + ti(child_age_mri,disadv_y3, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme)
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "disadv_y3", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

##Regional GAMMs for disadv ##
years=c("y1","y2","y3")
for (year in years){
  variable_name=paste0("disadv_",year)
  avgclustco_both.gordon.elabe <- network_demo_data_long_clustco #rename data to be what the function is looking for
  gam.age.ses.clustco.gordon <- matrix(data=NA, nrow=333, ncol=4) #matrix to save gam.fit output to
  
  for(row in c(1:nrow(gordon_parcel_matching))){ #for each region
    region <- str_c("avgclustco_all_", row)
    GAM.RESULTS <- gamm.fit.smooth.int(measure = "avgclustco_both", atlas = "gordon", dataset = "elabe", 
                                       region = region, smooth_var = "child_age_mri", int_var = variable_name,
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
  write.csv(gam.age.ses.clustco.gordon.fx.TRUE, paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_",year,"_disadv_gamm_statistics_gordon_fx_T_no_neg.csv"), row.names = F, quote = F)
  gam.age.ses.clustco.gordon.fx.TRUE <- read.csv(paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_",year,"_disadv_gamm_statistics_gordon_fx_T_no_neg.csv"))
  gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct
  #make a boxplot
  gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks
  kruskal.test(GAM.age.ses.Fvalue~network, data = gam.age.ses.clustco.gordon.fx.TRUE)
  g <- ggplot(data = dplyr::filter(gam.age.ses.clustco.gordon.fx.TRUE, network !="None"), aes(x = reorder(network, -GAM.age.ses.Fvalue, median), y= GAM.age.ses.Fvalue, fill=network)) +
    geom_boxplot(width = .6, alpha = 0.8) +
    geom_point(aes(x = reorder(network, -GAM.age.ses.Fvalue, median), y = GAM.age.ses.Fvalue, color = network), alpha=0.3) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_fill_gordon() +
    scale_color_gordon() +
    theme_bw() + main(variable_name)+
    ylab("Age x SES effects on local segregation (F)")+
    xlab("Networks")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size=14)) #fix the offset
  print(g)
}

# Supplemental Fig 6: Using psychosocial stress instead of social disadvantage -----------------------------

SEM_psychosocial_stress_y0 <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Tooley_09_26_2023_disadvantage_psychoscocial_stress.csv") %>% select(contains("psych_prenatal") | contains("MODID_prenatal"))
colnames(SEM_psychosocial_stress_y0) <- c("psych_prenatal", "MODID")
SEM_psychosocial_stress_y0$modid <- str_c("MOD",SEM_psychosocial_stress_y0$MODID)
network_demo_data_long<- left_join(network_demo_data_long, SEM_psychosocial_stress_y0, by="modid")
colSums(table(network_demo_data_long$child_age_mri,network_demo_data_long$timepoint))

birth_data <- network_demo_data_long %>% filter(timepoint=="y0")
cor.test(birth_data$psych_prenatal, birth_data$disadv_prenatal)

#Age x SES GAMMs
## System segregation ##
## just the regular analysis, psychosocial stress is associated with integration
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + avg_FD_of_retained_frames  + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+ggtitle("System segregation")+theme(plot.title = element_text(hjust = 0.5))

## Controlling for prenatal disadvantage as a fixed covariate, stress x age interaction is significant, same directionality as socioeconomic disadvantage
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + avg_FD_of_retained_frames + disadv_prenatal + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ avg_FD_of_retained_frames + disadv_prenatal + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+ggtitle("System segregation")+theme(plot.title = element_text(hjust = 0.5))

## Controlling for the interaction of SES x age, stress x age interaction is significant, same directionality as socioeconomic disadvantage
gam_age_ses_segreg_by <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by$lme);summary(gam_age_ses_segreg_by$gam)
gam_age_ses_segreg_ti <- gamm(system_segreg ~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ ti(child_age_mri,disadv_prenatal, k=4)+ ti(disadv_prenatal, k=4)+avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_ti$lme);summary(gam_age_ses_segreg_ti$gam)
BIC(gam_age_ses_segreg_by$lme , gam_age_ses_segreg_ti$lme)
a <- resid_plot_int(gam_age_ses_segreg_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+ggtitle("System segregation")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")

## Modularity ##
## Controlling for prenatal disadvantage as a fixed covariate, stress x age interaction not significant
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + avg_FD_of_retained_frames + disadv_prenatal + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ avg_FD_of_retained_frames + disadv_prenatal + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")

## Controlling for the interaction of SES x age, stress x age interaction not significant
gam_age_ses_modul_by <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by$lme);summary(gam_age_ses_modul_by$gam)
gam_age_ses_modul_ti <- gamm(modul ~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ ti(child_age_mri,disadv_prenatal, k=4)+ ti(disadv_prenatal, k=4)+avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_modul_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")+theme(legend.position = "none")

## Clustering coefficient ## 
## Controlling for prenatal disadvantage as a fixed covariate, stress x age interaction not significant
gam_age_ses_clustco_by <- gamm(avgclustco_both~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clustco_by$lme);summary(gam_age_ses_clustco_by$gam)
gam_age_ses_modul_ti <- gamm(avgclustco_both~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_clustco_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_clustco_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")

## Controlling for the interaction of SES x age, stress x age interaction not significant
gam_age_ses_clustco_by <- gamm(avgclustco_both~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clustco_by$lme);summary(gam_age_ses_clustco_by$gam)
gam_age_ses_modul_ti <- gamm(avgclustco_both~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ ti(child_age_mri,disadv_prenatal, k=4)+ ti(disadv_prenatal, k=4)+avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_clustco_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_clustco_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Participation coefficient
#just the regular analysis, psychosocial stress is associated with integration
gam_age_ses_part_coef_by <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_part_coef_ti <- gamm(part_coef_avg ~ child_sex + ti(child_age_mri, k=4)+  ti(child_age_mri,psych_prenatal, k=4)+ti(psych_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_ti$lme);summary(gam_age_ses_part_coef_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_part_coef_ti$lme)
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Participation coefficient (residuals)")+ggtitle("Participation coefficient")+theme(plot.title = element_text(hjust = 0.5))

## Controlling for prenatal disadvantage as a fixed covariate, stress x age interaction not significant
gam_age_ses_part_coef_by <- gamm(part_coef_avg~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + avg_FD_of_retained_frames + disadv_prenatal + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_modul_ti <- gamm(part_coef_avg~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ avg_FD_of_retained_frames + disadv_prenatal + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")+ggtitle("Participation coefficient")+theme(plot.title = element_text(hjust = 0.5))

## Controlling for the interaction of SES x age, stress x age interaction not significant
gam_age_ses_part_coef_by <- gamm(part_coef_avg~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=psych_prenatal, k=4) + s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_part_coef_by$lme);summary(gam_age_ses_part_coef_by$gam)
gam_age_ses_modul_ti <- gamm(part_coef_avg~ child_sex + ti(child_age_mri, k=4)+ ti(psych_prenatal, k=4) + ti(child_age_mri,psych_prenatal, k=4)+ ti(child_age_mri,disadv_prenatal, k=4)+ ti(disadv_prenatal, k=4)+avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_ti$lme);summary(gam_age_ses_modul_ti$gam)
BIC(gam_age_ses_part_coef_by$lme , gam_age_ses_modul_ti$lme)
a <- resid_plot_int(gam_age_ses_part_coef_by, "child_age_mri", "modid", int.term = "psych_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Participation coefficient (residuals)")+theme(legend.position = "none")

#adjust p's
p.adjust(c(summary(gam_age_ses_segreg_by$gam)$s.table[2,4],summary(gam_age_ses_modul_by$gam)$s.table[2,4],summary(gam_age_ses_clustco_by$gam)$s.table[2,4], summary(gam_age_ses_part_coef_by$gam)$s.table[2,4]), method = "fdr")

# Supplementary Fig 4: Random effect of age included in age x SES ----------------------
#First have to run all models allowing for both an random intercept and a random slope (uncorrelated bc correlated does not converge)!
#according to this answer https://stackoverflow.com/questions/25872488/how-to-add-a-random-intercept-and-random-slope-term-to-a-gamm-model-in-r

#segregation
gam_age_ses_segreg_by_random_slope <- gamm(system_segreg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1, modid=~0+child_age_mri), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_segreg_by_random_slope$lme);summary(gam_age_ses_segreg_by_random_slope$gam)

#modularity
gam_age_ses_modul_by_random_slope <- gamm(modul ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1, modid=~0+child_age_mri), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_modul_by_random_slope$lme);summary(gam_age_ses_modul_by_random_slope$gam)

#clustering coefficient
gam_age_ses_clustco_by_random_slope <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1, modid=~0+child_age_mri), data=network_demo_data_long, method = "REML", control =lmeControl(sing.tol=1e-20, niterEM=0))
#had to increase the singular tolerance of the lmeControl parameters to get this to converge, found here: https://groups.google.com/g/davis-rug/c/sbcjeMalCgc
#actually, it might just be that parameters are on very different scales, what if we scale all parameters?
network_demo_data_long$avgweight_sc <- as.numeric(scale(network_demo_data_long$avgweight));network_demo_data_long$child_age_mri_sc <- as.numeric(scale(network_demo_data_long$child_age_mri));network_demo_data_long$avg_FD_of_retained_frames_sc <- as.numeric(scale(network_demo_data_long$avg_FD_of_retained_frames));
network_demo_data_long$disadv_prenatal_sc<- as.numeric(scale(network_demo_data_long$disadv_prenatal));network_demo_data_long$retained_frames_sc <- as.numeric(scale(network_demo_data_long$retained_frames));network_demo_data_long$avgclustco_both_sc <- as.numeric(scale(network_demo_data_long$avgclustco_both));
network_demo_data_long$system_segreg_sc <- as.numeric(scale(network_demo_data_long$system_segreg));network_demo_data_long$modul_sc <- as.numeric(scale(network_demo_data_long$modul));network_demo_data_long$part_coef_avg_sc <- as.numeric(scale(network_demo_data_long$part_coef_avg));
network_demo_data_long$disadv_y1_sc<- as.numeric(scale(network_demo_data_long$disadv_y1));
gam_age_ses_clustco_by_random_slope <- gamm(avgclustco_both ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri, by=disadv_prenatal_sc, k=4) + avg_FD_of_retained_frames_sc + retained_frames_sc + avgweight_sc, random = list(modid =~ 1, modid=~0+child_age_mri), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_clustco_by_random_slope$lme);summary(gam_age_ses_clustco_by_random_slope$gam)

a <- resid_plot_int(gam_age_ses_clustco_by_random_slope, "child_age_mri", "modid", int.term = "disadv_prenatal_sc", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

#participation coefficient
gam_age_ses_partcoef_by_random_slope <- gamm(part_coef_avg ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight, random = list(modid =~ 1, modid=~0+child_age_mri), data=network_demo_data_long, method = "REML")
summary(gam_age_ses_partcoef_by_random_slope$lme);summary(gam_age_ses_partcoef_by_random_slope$gam)

#adjust p's
p.adjust(c(summary(gam_age_ses_segreg_by_random_slope$gam)$s.table[2,4],summary(gam_age_ses_modul_by_random_slope$gam)$s.table[2,4],summary(gam_age_ses_clustco_by_random_slope$gam)$s.table[2,4], summary(gam_age_ses_partcoef_by_random_slope$gam)$s.table[2,4]), method = "fdr")

## Regional age x SES GAMMs
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

##Regional GAMMs for local segregation ##
avgclustco_both_sc.gordon.elabe <- network_demo_data_long_clustco #rename data to be what the function is looking for
gam.age.ses.clustco.gordon <- matrix(data=NA, nrow=333, ncol=4) #matrix to save gam.fit output to

for(row in c(1:nrow(gordon_parcel_matching))){ #for each region
  region <- str_c("avgclustco_all_", row) #could also used the scaled variables here
  GAM.RESULTS <- gamm.fit.smooth.int(measure = "avgclustco_both_sc", atlas = "gordon", dataset = "elabe", 
                                     region = region, smooth_var = "child_age_mri_sc", int_var = "disadv_prenatal_sc",
                                     id_var="modid", random_slope=TRUE,
                                     covariates = "child_sex + avg_FD_of_retained_frames_sc + retained_frames_sc + avgweight_sc", 
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
write.csv(gam.age.ses.clustco.gordon.fx.TRUE, "~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_random_age_effect_scaled_covariates_gamm_statistics_gordon_fx_T.csv", row.names = F, quote = F)
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_random_age_effect_scaled_covariates_gamm_statistics_gordon_fx_T.csv")
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct

## Write out the age x SES effects to a pscalar so can look at it in workbench ##
values <- ifelse(gordon_networks=="None",0,gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.Fvalue)#remove 'None' network values from plotting
new.347.values <- rep(0, 347)#dim is 347, because of subcortical, so pad it with 0s
new.347.values[1:333] <- values
output_file=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_age_ses_random_age_effect_scaled_covariates_fvalue_no_none_no_neg.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_age_ses_random_age_effect_scaled_covariates_fvalue_no_none_no_neg.ptseries.nii")

command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

#make a boxplot
gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks
kruskal.test(GAM.age.ses.Fvalue~network, data = gam.age.ses.clustco.gordon.fx.TRUE)
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


# Response to reviewers: Checking network outcomes for outliers -----------------------------
#Exclude observations that are more than 3 SD away from the mean for each timepoint
#should exclude at that timepoint or exclude the participant entirely?

#3SD range for y0, modul
mean(y0_network_data$modul);sd(y0_network_data$modul)
max_modul_y0=mean(y0_network_data$modul)+(3*sd(y0_network_data$modul))
min_modul_y0=mean(y0_network_data$modul)-(3*sd(y0_network_data$modul))

#excluding from y0 observations only
network_demo_data_long_excl_y0 <- network_demo_data_long %>% filter(!(timepoint=="y0" & modul>max_modul_y0)) %>% filter(!(timepoint=="y0" & modul<min_modul_y0))
dim(network_demo_data_long);dim(network_demo_data_long_excl_y0)

colnames=c("system_segreg", "modul","avgclustco_both")
#colnames=c("system_segreg", "modul","avgclustco_both","part_coef_avg")
#can rotate which one is last to get exact numbers for each timepoint from the code below.
for (col_name in colnames) {
  # Access the column using df[, col_name] or df[[col_name]]
  max_y0 <- mean(all_network_data[all_network_data$timepoint=="y0",col_name]) + (3*sd(all_network_data[all_network_data$timepoint=="y0",col_name])) #make max and min values for y0
  min_y0 <- mean(all_network_data[all_network_data$timepoint=="y0",col_name]) - (3*sd(all_network_data[all_network_data$timepoint=="y0",col_name]))
  network_demo_data_long_excl_y0 <- network_demo_data_long %>% filter(!(timepoint=="y0" & get(col_name)>max_y0)) %>% filter(!(timepoint=="y0" & get(col_name)<min_y0))#remove these from dataset
  max_y2 <- mean(all_network_data[all_network_data$timepoint=="y2",col_name], na.rm=T) + (3*sd(all_network_data[all_network_data$timepoint=="y2",col_name])) #make max and min values for y2
  min_y2 <- mean(all_network_data[all_network_data$timepoint=="y2",col_name]) - (3*sd(all_network_data[all_network_data$timepoint=="y2",col_name]))
  network_demo_data_long_excl_y2 <- network_demo_data_long_excl_y0 %>% filter(!(timepoint=="y2" & get(col_name)>max_y2)) %>% filter(!(timepoint=="y2" & get(col_name)<min_y2))
  max_y3 <- mean(all_network_data[all_network_data$timepoint=="y3",col_name], na.rm=T) + (3*sd(all_network_data[all_network_data$timepoint=="y3",col_name])) #make max and min values for y3
  min_y3 <- mean(all_network_data[all_network_data$timepoint=="y3",col_name]) - (3*sd(all_network_data[all_network_data$timepoint=="y3",col_name]))
  network_demo_data_long_excl_y3 <- network_demo_data_long_excl_y2 %>% filter(!(timepoint=="y3" & get(col_name)>max_y3)) %>% filter(!(timepoint=="y3" & get(col_name)<min_y3))
  name<-paste0("gamm_gam_age_ses_",col_name,"_by")
  #gam_age_ses_segreg_by <- gamm(col_name ~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight , random = list(modid =~ 1), data=network_demo_data_long_excl_y3, method = "REML")
  formula<-formula(paste0(col_name, '~ child_sex + s(child_age_mri, k=4)+ s(child_age_mri,by=disadv_prenatal, k=4) + avg_FD_of_retained_frames + retained_frames + avgweight')) 
  assign(name, gamm(formula,random = list(modid =~ 1), data=network_demo_data_long_excl_y3, method = "REML"))
  print(name)
  print(summary(get(name)$gam))
}
network_demo_data_long_excl_y0 %>% filter(timepoint=="y0") %>% select(col_name) %>% summary(range=range)
dim(network_demo_data_long);dim(network_demo_data_long_excl_y0);dim(network_demo_data_long_excl_y2);dim(network_demo_data_long_excl_y3)
colSums(table(network_demo_data_long_excl_y3$child_age_mri,network_demo_data_long_excl_y3$timepoint)) #how many healthy FT at each timepoint?
modids_excl_y3 <- unique(network_demo_data_long_excl_y3$modid); length(modids_excl_y3);length(unique(network_demo_data_long$modid))

#everything looks fine here

#different outlier test
library(EnvStats)
rosnerTest(y0_network_data$system_segreg, k=5);rosnerTest(y0_network_data$modul, k=5);rosnerTest(y0_network_data$avgclustco_both, k=5)
rosnerTest(y2_network_data$system_segreg, k=5);rosnerTest(y2_network_data$modul, k=5);rosnerTest(y2_network_data$avgclustco_both, k=5)
rosnerTest(y3_network_data$system_segreg, k=5);rosnerTest(y3_network_data$modul, k=5);rosnerTest(y3_network_data$avgclustco_both, k=5)

#Age x SES GAMMs
## System segregation ##
a <- resid_plot_int(gamm_gam_age_ses_system_segreg_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("System segregation (residuals)")+theme(legend.position = "none")

## Modularity ##
a <- resid_plot_int(gamm_gam_age_ses_modul_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Modularity (residuals)")+theme(legend.position = "none")

## Clustering coefficient ## 
a <- resid_plot_int(gamm_gam_age_ses_avgclustco_both_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Clustering coefficient (residuals)")+theme(legend.position = "none")

## Participation coefficient
a <- resid_plot_int(gamm_gam_age_ses_part_coef_avg_by, "child_age_mri", "modid", int.term = "disadv_prenatal", add.intercept=T)
a+xlab("Child age at MRI (years)")+ylab("Participation coefficient (residuals)")+theme(legend.position = "none")

#adjust p's
p.adjust(c(summary(gamm_gam_age_ses_system_segreg_by$gam)$s.table[2,4],summary(gamm_gam_age_ses_modul_by$gam)$s.table[2,4],summary(gamm_gam_age_ses_avgclustco_both_by$gam)$s.table[2,4], summary(gamm_gam_age_ses_part_coef_avg_by$gam)$s.table[2,4]), method = "fdr")

## Regional age x SES GAMMs
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
network_demo_data_long_clustco <- left_join(network_demo_data_long_excl_y3,all_clustco_data, by=c("modid", "timepoint"))

##Regional GAMMs for local segregation ##
avgclustco_both.gordon.elabe <- network_demo_data_long_clustco #rename data to be what the function is looking for
gam.age.ses.clustco.gordon <- matrix(data=NA, nrow=333, ncol=4) #matrix to save gam.fit output to

for(row in c(1:nrow(gordon_parcel_matching))){ #for each region
  region <- str_c("avgclustco_all_", row) #could also used the scaled variables here
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
write.csv(gam.age.ses.clustco.gordon.fx.TRUE, "~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_no_outliers_gamm_statistics_gordon_fx_T.csv", row.names = F, quote = F)
gam.age.ses.clustco.gordon.fx.TRUE <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/GAMM_results/n419_avgclustco_age_ses_no_outliers_gamm_statistics_gordon_fx_T.csv")
gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue.fdr <- p.adjust(gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.pvalue, method = "fdr")#FDR correct

## Write out the age x SES effects to a pscalar so can look at it in workbench ##
values <- ifelse(gordon_networks=="None",0,gam.age.ses.clustco.gordon.fx.TRUE$GAM.age.ses.Fvalue)#remove 'None' network values from plotting
new.347.values <- rep(0, 347)#dim is 347, because of subcortical, so pad it with 0s
new.347.values[1:333] <- values
output_file=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_age_ses_no_outliers_fvalue_no_none_no_neg.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/within_between_network_longitudinal/data/wb_viz/clustco/avgclustco_gamm_age_ses_no_outliers_fvalue_no_none_no_neg.ptseries.nii")

command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

#make a boxplot
gam.age.ses.clustco.gordon.fx.TRUE$network <- gordon_networks
kruskal.test(GAM.age.ses.Fvalue~network, data = gam.age.ses.clustco.gordon.fx.TRUE)
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

# Included versus excluded participants, MRI and behavioral ---------------
#Included MRI versus excluded MRI participants
unusable_birth_data <- filter(demo_data_all, !(modid %in% y0_network_data$modid))
modids_missing_at_year2 <- filter(demo_data_all, !(modid %in% y2_network_data$modid)) %>% select(modid);modids_missing_at_year2 <- left_join(modids_missing_at_year2, demo_data_all, by="modid")# %>%filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0); 
modids_missing_at_year2 <- left_join(modids_missing_at_year2, SEM_birth_indicators, by="modid");modids_missing_at_year2 <- left_join(modids_missing_at_year2, SEM_psychosocial_stress_y0, by="modid"); 
modids_missing_at_year2 <- left_join(modids_missing_at_year2, ADI_all_timepoints, by="modid")
y2_data <- left_join(y2_network_data,demo_data_all, by="modid");y2_data <- left_join(y2_data, SEM_birth_indicators, by="modid");y2_data <- left_join(y2_data, SEM_psychosocial_stress_y0, by="modid"); 
y2_data <- left_join(y2_data, ADI_all_timepoints, by="modid")#%>%filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0)
modids_missing_at_year3 <- filter(demo_data_all, !(modid %in% y3_network_data$modid)) %>% select(modid);modids_missing_at_year3 <- left_join(modids_missing_at_year3,  demo_data_all, by="modid")# %>%filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0); 
modids_missing_at_year3 <- left_join(modids_missing_at_year3, SEM_birth_indicators, by="modid");modids_missing_at_year3 <- left_join(modids_missing_at_year3, SEM_psychosocial_stress_y0, by="modid"); 
modids_missing_at_year3 <- left_join(modids_missing_at_year3, ADI_all_timepoints, by="modid")
y3_data <- left_join(y3_network_data,demo_data_all, by="modid");y3_data <- left_join(y3_data, SEM_birth_indicators, by="modid");y3_data <- left_join(y3_data, SEM_psychosocial_stress_y0, by="modid"); 
y3_data <- left_join(y3_data, ADI_all_timepoints, by="modid")#%>%filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0)

#Differences in ADI, disadvantage, sex, gestational age, birthweight in those who have usable MRI data at a timepoint versus those who don't
t.test(unusable_birth_data$child_birthweight, birth_data$child_birthweight)
t.test(unusable_birth_data$GAWEEKS, birth_data$GAWEEKS)
t.test(unusable_birth_data$income_needs_demo_b, birth_data$income_needs_demo_b)
t.test(unusable_birth_data$disadv_prenatal, birth_data$disadv_prenatal)
#t.test(unusable_birth_data$psych_prenatal, birth_data$psych_prenatal)
t.test(unusable_birth_data$ADI, birth_data$ADI)

#Years 2 and 3 usable MRI vs not
t.test(y2_data$child_age_y2, modids_missing_at_year2$child_age_y2)
t.test(y2_data$disadv_prenatal, modids_missing_at_year2$disadv_prenatal)
t.test(y3_data$disadv_prenatal, modids_missing_at_year3$disadv_prenatal)
t.test(y2_data$neighborhood_natlcentile_b.x, modids_missing_at_year2$ADI)
t.test(y3_data$neighborhood_natlcentile_b.x, modids_missing_at_year3$ADI)
chisq.test(cbind(table(y2_data$child_sex), table(modids_missing_at_year2$child_sex)))
chisq.test(cbind(table(y3_data$child_sex), table(modids_missing_at_year3$child_sex)))
t.test(y2_data$GAWEEKS, modids_missing_at_year2$GAWEEKS)
t.test(y3_data$GAWEEKS, modids_missing_at_year3$GAWEEKS)
t.test(y2_data$child_birthweight, modids_missing_at_year2$child_birthweight)
t.test(y3_data$child_birthweight, modids_missing_at_year3$child_birthweight)
t.test(y2_data$psych_prenatal, modids_missing_at_year2$psych_prenatal)
t.test(y3_data$psych_prenatal, modids_missing_at_year3$psych_prenatal)
t.test(y2_data$income_needs_demo_b, modids_missing_at_year2$income_needs_demo_b)
t.test(y2_data$income_needs_demo_y2, modids_missing_at_year2$income_needs_demo_y2)
t.test(y3_data$income_needs_demo_b, modids_missing_at_year3$income_needs_demo_b)
t.test(y3_data$income_needs_demo_y3, modids_missing_at_year3$income_needs_demo_y3)

#look at those who came in at all at a time point versus those who didn't
modids_missing_at_year2 <- filter(demo_data_all, is.na(child_age_y2_assessment)); #who doesn't have age assessed at y2
modids_missing_at_year2 <- left_join(modids_missing_at_year2, SEM_birth_indicators, by="modid");modids_missing_at_year2 <- left_join(modids_missing_at_year2, SEM_psychosocial_stress_y0, by="modid"); 
modids_missing_at_year2 <- left_join(modids_missing_at_year2, ADI_all_timepoints, by="modid")
y2_data <- filter(demo_data_all, !is.na(child_age_y2_assessment)); y2_data <- left_join(y2_data, SEM_birth_indicators, by="modid");y2_data <- left_join(y2_data, SEM_psychosocial_stress_y0, by="modid"); 
y2_data <- left_join(y2_data, ADI_all_timepoints, by="modid");#%>%filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0)
modids_missing_at_year3 <- filter(demo_data_all, is.na(child_age_y3_assessment)); #who doesn't have age assessed at y2
modids_missing_at_year3 <- left_join(modids_missing_at_year3, SEM_birth_indicators, by="modid");modids_missing_at_year3 <- left_join(modids_missing_at_year3, SEM_psychosocial_stress_y0, by="modid"); 
modids_missing_at_year3 <- left_join(modids_missing_at_year3, ADI_all_timepoints, by="modid")
y3_data <- filter(demo_data_all, !is.na(child_age_y3_assessment));y3_data <- left_join(y3_data, SEM_birth_indicators, by="modid");y3_data <- left_join(y3_data, SEM_psychosocial_stress_y0, by="modid"); 
y3_data <- left_join(y3_data, ADI_all_timepoints, by="modid")#%>%filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0)

t.test(y2_data$disadv_prenatal, modids_missing_at_year2$disadv_prenatal)
t.test(y3_data$disadv_prenatal, modids_missing_at_year3$disadv_prenatal)
t.test(y2_data$neighborhood_natlcentile_b.x, modids_missing_at_year2$ADI)
t.test(y3_data$neighborhood_natlcentile_b.x, modids_missing_at_year3$ADI)
chisq.test(cbind(table(y2_data$child_sex), table(modids_missing_at_year2$child_sex)))
chisq.test(cbind(table(y3_data$child_sex), table(modids_missing_at_year3$child_sex)))
t.test(y2_data$GAWEEKS, modids_missing_at_year2$GAWEEKS)
t.test(y3_data$GAWEEKS, modids_missing_at_year3$GAWEEKS)
t.test(y2_data$child_birthweight, modids_missing_at_year2$child_birthweight)
t.test(y3_data$child_birthweight, modids_missing_at_year3$child_birthweight)
t.test(y2_data$psych_prenatal, modids_missing_at_year2$psych_prenatal)
t.test(y3_data$psych_prenatal, modids_missing_at_year3$psych_prenatal)
t.test(y2_data$income_needs_demo_b, modids_missing_at_year2$income_needs_demo_b)
t.test(y2_data$income_needs_demo_y2, modids_missing_at_year2$income_needs_demo_y2)
t.test(y3_data$income_needs_demo_b, modids_missing_at_year3$income_needs_demo_b)
t.test(y3_data$income_needs_demo_y3, modids_missing_at_year3$income_needs_demo_y3)

combined_data <- rbind(
  data.frame(Group = "Included", Value = y2_data$disadv_prenatal),
  data.frame(Group = "Missing", Value = modids_missing_at_year2$disadv_prenatal)
)

# Create boxplot using ggplot with jittered data points
ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.5) +  # Make boxplot slightly transparent
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +  # Jittered data points
  labs(title = "Included versus excluded participants at Y2",
       x = "Group",
       y = "Prenatal disadvantage") +
  scale_fill_manual(values = c("blue", "red"))
