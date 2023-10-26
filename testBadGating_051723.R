#workflow.R
#Started: Sept 21, 2021
#Authors: Julie Chuong, Titir De, Grace Avecilla, David Gresham

##Install CytoExplorer package and requirements (can be skipped if already installed)
library(BiocManager)
BiocManager::install("cytolib")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("openCyto")

##Install CytoExploreR from GitHub:
library(devtools)
devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(ggridges)
library(docstring)

setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_WT")  #Julie's WD

# Read in the names of those directories (data directories and one gating directory) here:
folders = list.dirs()

# Choose a name to be used for all output files including the gating template and associated flow data and graphs.
version_name = "g153_WT_gate"

gating_dir = folders[16] #change this folder for your gating directory
exp_details_path = list.files(path = paste0(gating_dir), pattern = "_experiment_details.csv", full.names = T)

timepoint_gating_set <- cyto_setup(path = paste0(gating_dir), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

# use flowWorkspace::pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv

ordered_exp_details = pData(timepoint_gating_set) %>% left_join(experiment_details) #rerrange rows of data frame so merging is correct. ie. check that the .fcs name matches the sample name in attached metadata
for(i in 1:length(names(ordered_exp_details))){
  flowWorkspace::pData(timepoint_gating_set)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
}

cyto_details(timepoint_gating_set) %>% View() #check correct attachment of metadata

# Rename the experiment-markers.csv file. Need to do once.
#file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv")

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir & Julie

## 3.1 transform the data
# looks useful if I want to choose different transformation: https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html
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

## quickly check the transformation by plotting the data
#my_samples = which(experiment_details$Description %in% c("0 copy control","1 copy control","2 copy control","GAP1 LTR KO", "GAP1 WT architecture"))
quartz()
cyto_plot_explore(transformed_timepoint_gating_set,
                  channels_x = "FSC-A",
                  channels_y = "B2-A",
                  axes_limits = "data")

## 3.2 Gating using the entire timepoint dataset or apply an existing gating template

# note:if you already have a gating template and don't need to draw gates, then skip cyto_draw, use cyto_gatingTemplate_apply to apply the gating template.csv to your gating set
#cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")

#First we gate for the cells
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)


#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)

# Draw zero-copy, one-copy, CNV gates 
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)

# Apply gate to the folder again and get cell frequency as sanity check




#STEP 5:  Use function to perform analysis
#A function that will
#1 Read in all the files in a folder
#2 Read in experiment details files using pData
#3 Specify experiment markers
#4 Transform gating set
#5 Apply existing gating file using cyto_gatingTemplate_apply
#6.Output stats file as .csv
#Author: David & Julie

my_markers<-c("GFP") #list your marker name(s)
channel<-c("B2-A") #list your channel(s)
names(my_markers)<-channel

analyze_all_exp = function(folder_name, my_markers, gating_template="cytek_gating_g153_WT_gate.csv") {

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



}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

analyze_all_exp(gating_dir, my_markers = my_markers, gating_template = "cytek_gating_g153_WT_gate.csv")


count = read_csv("g153_WT_gate_counts_16_EE_GAP1_ArchMuts_2021.csv")
freq = read_csv("g153_WT_gate_freq_16_EE_GAP1_ArchMuts_2021.csv")

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

g153_gate_freqs = freq_and_counts %>% filter(generation == 153, 
                           Gate == "two_or_more_copy",
                           Description == "GAP1 WT architecture") %>% arrange(sample) %>% View() 
#Frequencies are 69,61,74,76,66

# Compare to freq and counts g153 in my current table
current = read_csv("~/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files/freq_and_counts_031523_MERGED_Experimentals_noControls.csv")
current = current %>% filter(generation == 153, 
                             Gate == "two_or_more_copy",
                             Description == "GAP1 WT architecture") %>% arrange(sample)
#Frequencies are 69, 57, 72, 75. 63 

#Conclusion my current gates are fine! 