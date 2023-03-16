# 3-15-23
# Forked from gating_112122.R
# Purpose: add g0 flow data, make lineplots up to g240.
# Previously have been doing up to including g203, thinking it was messy after
# but haven't actually looked in a while
# I'm not convinced it SHOULDN'T be included

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(ggridges)
library(docstring)

setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

####################################
##### Find timepoint with lowest median normalized GFP for each of the 4 genotypes ####

norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv")

medianGFP = norm_medians %>% filter(generation <=100) %>% group_by(Description, generation) %>% summarize(median = median(Med_B2A_FSC)) %>% slice(which.min(median)) #grouped by Description and took median across each population within its genotype

#medianGFP %>% write_csv("min-median-norm-GFP_112222.csv")
medianGFP = read_csv("min-median-norm-GFP_112222.csv")

# Description         generation    median   timepoint

#  0 copy control               58  0.403     06
#  1 copy control               29  1.74      03
#  2 copy control               29  1.99      03
#  GAP1 ARS KO                  50  1.86      05
#  GAP1 LTR + ARS KO            50  1.93      05
#  GAP1 LTR KO                  37  1.74      04
#  GAP1 WT architecture         21  1.73      02
###########################################################

# 3-15-23 Redo CNV proportion lineplot with generation 0-240 and see what it tells us
# Consequently, redo the median across population lineplot. 
# Might have to redo quantify dynamics boxplots after
#Select FCS directories for generation 0-240
folders = list.dirs()[c(2, 10:11,13,15,17,20:23,25:39)] #select the FSC file folders in your directory

# Choose a name to be used for all output files including the gating template and associated flow data and graphs.

version_name = "03_one_ctrl_2Gates"
version_name = "03_two_ctrl" #draw a two copy gate, not 0 or 1 copy gates.
version_name = "03_one_ctrl_1Gate" # just draw a one copy gate, no other gates.
version_name = "06_zero_ctrl" # just to draw a zero copy gate, no other gates
version_name = "05_ALL_121522" #gates drawn, need to apply
version_name = "04_LTR_112222" #applied to folders
version_name = "02_WT_112222" #applied to folders
#version_name = "03_112122_allko" #Timepoint 3 only because it had the lowest median GFP for ALLKO pops. THIS IS NOT TRUE. TIMEPOINT 5 had the lowest. Use this instead: "05_ALL_120822
#version_name = "05_112122_ars" # ARS KO samples only. Timepoint 5 because it has the lowest median GFP. ARS KO only. Because the ancestors have strain-specific GFP that is higher than that of the WT and LTR. Therefore, a separate gating template is needed.
#version_name = "03_112122_cons" #Timepoint 3 only because it has the lowest median GFP. cons = conservative meaning higher border between 1 and 2 copy gates.
#version_name = "03_112122_liberal" #Timepoint 3 only because it has the lowest median GFP. Liberal meaning lower border between 1 and 2 copy gates, which will lead to high positive rate and lowest threshold for CNV detection.
#version_name = "01_02_04_v2_111522"  #Used gating template 01_02_04_v2. This version reflect reanalysis with ordered_exp_details in Step 2 which now accurately attaches metadata to gating set.

# DO NOT USE these versions below.
# "01_02_04_v2_111522" replaces "01_02_04_v2"
# _____to do_____ replaces ""newGates_01_02_04_ars_all"
# version_name = "01_02_03_v5_wt_ltr" #timepoints 1 2 3 using WT and LTR KO experimental samples only - no controls - see script gating_111122_WTLTR.R DO NOT USE
#version_name = "01_02_04_v4_wt_ltr" #timepoints folders 1,2,4 using WT and LTR and control samples only to draw gates (folder 01_02_04_wt_ltr_ctrls_only_FSCfiles) using 0,1,2 copy as guides.
#version_name = "01_02_04_v3_wt_ltr" #timepoints folders 1,2,4 [using WT and LTR and control samples only](<- i don't think this is true. I think all samples were used to draw gates) using 0,1 copy as guides.
#version_name = "newGates_01_02_04_ars_all" #timepoint tolders 1,2,4, using only ARS and LTR+ARS samples (folder 01_02_04_ars_all_only_FSCfiles) only to draw gates
#version_name = "01_02_04_v2" #timepponts folders 1,2,4 all samples (experimental and controls) to draw gates
# other versions: 01_02_04_v2

#STEP 1: Generate experiment details file from folder and FCS file names
# Experiment details file is a .csv file that contains the list of .fcs files in the directory and the associated metadata for each sample
#Author: Grace

make_exp_details = function(folder_name, samplesheet) {
  pref = folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021")
  generation = folder_name %>% str_extract("[g]\\d+") %>% str_remove("g")

  files = as_tibble(list.files(paste0(folder_name))) %>%
    separate(value, into = c("well", "samp"), sep = " ", remove = F) %>%
    mutate(well = str_extract(well, "([A-Z])([0-9]){1,2}$")) %>%
    mutate(samp = str_remove(samp, ".fcs")) %>%
    mutate(sample = case_when(str_detect(value, "Unstained") ~ "ctrl0",
                              str_detect(value, "DGY500") ~ "ctrl1",
                              str_detect(value, "DGY1315") ~ "ctrl2",
                              TRUE ~ samp)) %>% 
    select(value,sample) %>% 
    rename(name = value) %>% 
    filter(!is.na(sample))

  all = files %>%
    left_join(read_csv(paste0("./",samplesheet)), by = c("sample" = "Sample name")) %>%
    mutate(generation = as.numeric(generation))

  write_csv(all, file = paste0(folder_name,"/",pref,"_experiment_details.csv"))

}

# Make Experiment Details for g0 folder
make_exp_details(folders[1], samplesheet = "EE_GAP1_ArchMuts_2021.csv")


# Skip STEPS 2-7 because we already hhave gating templates 
# and applied them to each flow genotype dataset 
# Just need to read in the freq and count data, select the timepoints we want
# and plot! 

##### Apply existing gates to g0 flow data ####

# g0 is folders[1] 

my_markers<-c("GFP") #list your marker name(s)
channel<-c("B2-A") #list your channel(s)
names(my_markers)<-channel

analyze_all_exp = function(folder_name, my_markers, gating_template="cytek_gating.csv") {
  
  path <- folder_name #gets relative path name for folder to be analyzed
  
  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name
  
  exp_details_path <- paste0(path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory
  #exp_details_path <- list.files(path = paste0(path), pattern = "_experiment_details.csv", full.names = T)
  #1. read in files and make a gating set
  print(path)
  timepoint_gating_set <- cyto_setup(path=path, select="fcs", details=F, markers = F)
  
  #2. read in experiment details for that gating set
  experiment_details <- read_csv(exp_details_path, show_col_types = F) #import experiment-details.csv
  #Write For Loop: for column in exp_details_path, add that column to timepoint_gating_set's metadata
  
  ordered_exp_details = pData(timepoint_gating_set) %>% left_join(experiment_details) #rerrange rows of data frame merging is correct. ie. fcs name matches the metadata
  for(i in 1:length(names(ordered_exp_details))){
    flowWorkspace::pData(timepoint_gating_set)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
  }
  
  #3. specify markers for that gating set
  markernames(timepoint_gating_set)<-my_markers
  
  #4. transform data
  GFP_trans <- cyto_transformer_logicle(timepoint_gating_set,
                                        channels = c("B2-A"),
                                        widthBasis = -10
  )#returns it as a list
  FSC_SSC_trans <- cyto_transformer_log(timepoint_gating_set,
                                        channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")
  )
  combined_trans <- cyto_transformer_combine(GFP_trans,FSC_SSC_trans)
  transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
                                                     trans = combined_trans) #applies the the transformation and returns it as a gatingSet
  
  #5. apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= gating_template)
  #  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")
  
  #6. Get cell counts and frequencies inside each gate
  #Julie
  
  #get cell count from each gate
  gs_pop_get_stats(transformed_timepoint_gating_set, c("Single_cells", "zero_copy", "one_copy", "two_or_more_copy")) %>%
    rename(Gate = pop, name = sample, Count = count) %>%
    left_join(experiment_details) %>%
    write_csv(paste0(version_name,"_counts_",prefix,".csv"))
  
  #get frequency of cells inside each gate
  gs_pop_get_stats(transformed_timepoint_gating_set, c("Single_cells","zero_copy", "one_copy", "two_or_more_copy"), type = "percent") %>%
    rename(Gate = pop, name = sample, Frequency = percent) %>%
    left_join(experiment_details) %>%
    write_csv(paste0(version_name,"_freq_",prefix,".csv"))
  
  #get single cell fluorescence normalized over cell size
  timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices
  
  map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
    mutate(name = as.factor(name)) %>% #convert `name` to factor
    left_join(experiment_details %>% #join by name column to add metadata
                mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
    mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor
    write_csv(paste0(version_name,"_SingleCellDistributions_",prefix,".csv"))
  
}
# Since we have four gating template we will just apply it four times,
# changing the gating_template.csv each time
version_name = "g0_WT_gateApply"
analyze_all_exp(folders[1], my_markers, gating_template = "cytek_gating_02_WT_112222.csv")

version_name = "g0_LTR_gateApply"
analyze_all_exp(folders[1], my_markers, gating_template = "cytek_gating_04_LTR_112222.csv")

version_name = "g0_ARS_gateApply"
analyze_all_exp(folders[1], my_markers, gating_template = "cytek_gating_05_112122_ars.csv")

version_name = "g0_ALL_gateApply"
analyze_all_exp(folders[1], my_markers, gating_template ="cytek_gating_05_ALL_121522.csv")

#### Pull in all the g0_count csv files, subset, and merge 
count_ALL = read_csv("g0_ALL_gateApply_counts_00_EE_GAP1_ArchMuts_2021.csv") %>%
filter(Description == "GAP1 LTR + ARS KO")
count_ARS = read_csv("g0_ARS_gateApply_counts_00_EE_GAP1_ArchMuts_2021.csv")  %>%
  filter(Description == "GAP1 ARS KO")
count_WT = read_csv("g0_WT_gateApply_counts_00_EE_GAP1_ArchMuts_2021.csv")  %>%
  filter(Description == "GAP1 WT architecture")
count_LTR = read_csv("g0_LTR_gateApply_counts_00_EE_GAP1_ArchMuts_2021.csv") %>%
  filter(Description == "GAP1 LTR KO")

g0_count = bind_rows(count_WT, count_ARS, count_ALL, count_LTR)

#### Pull in all the g0_Freq csv files, subset, and merge 
freq_ALL = read_csv("g0_ALL_gateApply_freq_00_EE_GAP1_ArchMuts_2021.csv") %>%
  filter(Description == "GAP1 LTR + ARS KO") %>% 
  mutate(gating_template = "cytek_gating_05_ALL_121522.csv")
freq_ARS = read_csv("g0_ARS_gateApply_freq_00_EE_GAP1_ArchMuts_2021.csv")  %>%
  filter(Description == "GAP1 ARS KO") %>%
  mutate(gating_template = "cytek_gating_05_112122_ars.csv")
freq_WT = read_csv("g0_WT_gateApply_freq_00_EE_GAP1_ArchMuts_2021.csv")  %>%
  filter(Description == "GAP1 WT architecture") %>%
  mutate(gating_template = "cytek_gating_02_WT_112222.csv")
freq_LTR = read_csv("g0_LTR_gateApply_freq_00_EE_GAP1_ArchMuts_2021.csv") %>%
  filter(Description == "GAP1 LTR KO") %>%
  mutate(gating_template = "cytek_gating_04_LTR_112222.csv")
g0_freq = bind_rows(freq_ALL, freq_ARS, freq_WT, freq_LTR)
  
g0_freq_and_counts =
  g0_count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(g0_freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

full_join(freq_and_counts, g0_freq_and_counts) %>% write_csv("freq_and_counts_031523_MERGED_Experimentals_noControls.csv")

##### STEP 8: Plot cells in the CNV gate #####
#freq_and_counts = read_csv("freq_and_counts_121622_MERGED_Experimentals_noControls.csv")
freq_and_counts = read_csv("freq_and_counts_031523_MERGED_Experimentals_noControls.csv")
freq_and_counts %>% filter(generation == 0, Gate == "two_or_more_copy") %>% View()
 
# Plot the data to identify outliers or general weird stuff 

my_facet_names <- as_labeller(c("GAP1 WT architecture" = "Wildtype architecture",
                                "GAP1 LTR KO" = "LTR KO",
                                "GAP1 ARS KO" = "ARS KO",
                                "GAP1 LTR + ARS KO" = "LTR and ARS KO"))
wtGrays = c("#354f52","#666666","#6b705c","#414833","#999999")
allGolds = c("#ffba08", "#faa307", "#dda15e", "#7f5539", "#9c6644", "#fdc409", "#9c7e1e","#D9BB59")
arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631", "#f85c46", "#bb3521","#d9402a" )
ltrBlues = c( "#6699cc", "#005f73", "#0a9396", "#4292C6", "#2171B5", "#3799fb", '#66b3cc', "#3a0ca3")

quartz()
freq_and_counts %>%
#clean_freq_and_counts %>%
  filter(Count>70000) %>%
  filter(generation != 252) %>% 
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~factor(Description,
                     levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
 # theme_classic() +
  #scale_x_continuous(breaks=seq(0,250,50)) +
  scale_x_continuous(breaks=seq(0,260,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
  )

ggsave("propCNV_0to260_grid_031523.png",bg = "#FFFFFF", height = 6, width = 10)

######## Lineplot per population ########
pop_lineplot = freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  facet_wrap(~sample) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
 # theme_classic() +
  scale_x_continuous(breaks=seq(0,260,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
  )

pop_lineplot

ggsave("propCNV_pop_0to260_grid_031523_10x14.pdf", bg = "#FFFFFF", height = 10, width = 14)
ggsave("propCNV_pop_0to260_grid_031523_10x14.png", bg = "#FFFFFF", height = 10, width = 14)

##### NEXT STEP ####
# redo normalized_median_plots_wControls.R to find the lowest median normalized GFP 
# hopefully its g0 hahaha? 
# Generally I would not trust g0 given that g8<50 GFP goes up which we think might be
# some cellular remodeling reacting to ggrowing in the chemostat, 
# when do we turn on the pumps? before or after g0? 
# ~g50 is our currentl lowest  


####### Look at ARS_7 after g182 #######
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  filter(sample == "gap1_ars_7") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_point(size = 4) +
  geom_line(size = 2.5) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications")+
  scale_x_continuous(breaks=seq(0,260,50))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
axis.text.y = element_text(size = 30, color = "black"))

ggsave("ARS7_grid_031523.png")
# Decision is to omit generations >=182 for ARS_7 

#########  Look at ARS 5 at g252 ########
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  filter(sample == "gap1_ars_5") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_point(size = 4) +
  geom_line(size = 2.5) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications")+
  scale_x_continuous(breaks=seq(0,260,50))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"))

##### Look at LTR_4 at g252 ####
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  filter(sample == "gap1_ltr_4") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_point(size = 4) +
  geom_line(size = 2.5) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications")+
  scale_x_continuous(breaks=seq(0,260,50))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"))

##### Look at LTR_3 at g240 #####
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  filter(sample == "gap1_ltr_3") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_point(size = 4) +
  geom_line(size = 2.5) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications")+
  scale_x_continuous(breaks=seq(0,260,50))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"))

##### Look at ALL_5 at g252 #####
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  filter(sample == "gap1_all_5") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_point(size = 4) +
  geom_line(size = 2.5) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications")+
  scale_x_continuous(breaks=seq(0,260,50))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"))
##### Look at ALL_7 at g252 #####
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  filter(sample == "gap1_all_7") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_point(size = 4) +
  geom_line(size = 2.5) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications")+
  scale_x_continuous(breaks=seq(0,260,50))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"))

#####  chose these timepoints by eye ####
#### I cannot cherry-pick timepoints out for the paper ####
weird_tp = freq_and_counts %>%
  filter(sample == "gap1_4" & Gate == "two_or_more_copy" & generation == 66 |  #use slope rule
           sample == "gap1_2" & Gate == "two_or_more_copy" & generation == 108|  #use slope rule
           sample == "gap1_ltr_2" |  
           sample == "gap1_ltr_4" & Gate == "two_or_more_copy" & generation == 21 | # 17.3%
           sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 166|
           sample == "gap1_all_5" & Gate == "two_or_more_copy" & generation == 116|
           sample == "gap1_all_6" & Gate == "two_or_more_copy" & generation == 124|
           sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 79 |
           sample == "gap1_ars_7" & Gate == "two_or_more_copy" & generation >= 182 |
           sample == "gap1_ars_5" & Gate == "two_or_more_copy" & generation == 252 |
           sample == "gap1_ars_3" & Gate == "two_or_more_copy" & generation == 252 |
           sample == "gap1_ltr_4" & Gate == "two_or_more_copy" & generation == 252 |
           sample == "gap1_ltr_3" & Gate == "two_or_more_copy" & generation == 240 |
           sample == "gap1_all_5" & Gate == "two_or_more_copy" & generation == 252 |
           sample == "gap1_all_7" & Gate == "two_or_more_copy" & generation == 252 
           )
clean_freq_and_counts = freq_and_counts %>%
  filter(Count>70000,
        Gate %in% c("two_or_more_copy"), 
        Type == "Experimental") %>%
  anti_join(weird_tp)

# clean_freq_and_counts %>% write_csv("freq_and_counts_merged_CLEAN_121622.csv")
# clean_freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

# ggsave("propCNVpop_clean_121622.pdf", bg = "#FFFFFF", height = 8, width = 12)
# ggsave("propCNVpop_clean_121622.png", bg = "#FFFFFF", height = 8, width = 12)



########## I can't use the following methods to justify outliers! ####### 
##### I can just omit them for presentation purposes like in Friday Seminar but NOT for my paper 
# Look at early timepoints 
# most are NOT outliers using the IQR rule 
# only LTR_4 g21 is outlier using IQR rule
# so don't anti_join(weird_pre50)
weird_pre50 =
  freq_and_counts %>%
  filter(generation < 50,
         Description %in% c("GAP1 WT architecture","GAP1 LTR KO"),
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample) %>%
  select(-name, -`Outflow well`, -Media, -Strain, -Type, -Description) %>%
  filter(Frequency > 6) %>% View()

###### IQR to detect outliers within a timepoint 
# Talking to David- I cant use this method of detecting outliers because this assumes each POPULATION acts the same or comes from the same distribution. Again because evolution and adaptation is stochastic, we CAN'T assume this. So this is NOT the way forward. 
x = freq_and_counts %>%
  filter(generation == 21,
         #  Description == "GAP1 WT architecture",
         Description == "GAP1 LTR KO",
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample) %>%
  select(sample, Description, generation, Gate, Parent, Count, Frequency, gating_template) 

boxplot(x$Frequency,ylab = "Frequency of GAP1 CNVs")
outlier_val = boxplot.stats(x$Frequency)$out  #the outlier Frequency value 
x[x$Frequency == outlier_val , ] # LTR_4, g21, 17.3% freq is outlier 





