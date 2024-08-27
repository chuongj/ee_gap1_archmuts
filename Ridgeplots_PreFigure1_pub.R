# This script is for gating for flow data using CytoExploreR
# and making ridgeplots   

##Install CytoExplorer package and requirements (can be skipped if already installed)
library(BiocManager)
install.packages("cytolib", "flowCore", "flowWorkspace", "openCyto")

##Install CytoExploreR from GitHub:
library(devtools)
devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(ggridges)
library(docstring)

setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD


#In addition to having directories (one to many) containing data FSC files, make a gating directory, which is **ONE** directory that contains ALL the FSC files you want to overlay for drawing gates. Read in the names of those directories (data directories and one gating directory) here:
folders = list.dirs()[c(9:36)] #select the FSC file folders in your directory
gating_dir = folders[17] #change this folder for your gating directory

# Choose a name to be used for all output files including the gating template and associated flow data and graphs.

version_name = "05_ALL_121522" 
version_name = "04_LTR_112222" 
version_name = "02_WT_112222" 
version_name = "05_112122_ars" 


#### STEP 1: Generate experiment details file from folder and FCS file names ####
# Experiment details file is a .csv file that contains the list of .fcs files in the directory and the associated metadata for each sample

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
    drop_na()

  all = files %>%
    left_join(read_csv(paste0("./",samplesheet)), by = c("sample" = "Sample name")) %>%
    mutate(generation = as.numeric(generation))

  write_csv(all, file = paste0(folder_name,"/",pref,"_experiment_details.csv"))
}

#needs to be run once
map(folders[1:length(folders)], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv")

#### STEP 2: Read in all files in a directory and rename the channels.####
# NOTE!!!:Skip Step 2-3 if you already have a gating template and want to apply it to data. Proceed to Step 4.

#A directory contains an FCS file for each population.
#Results in 1 timepoint gating set containing all .fcs files, associated experiment details, and marker details

# here we will load in my gating directory.
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
file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv")

#### STEP 3:  Perform gating on gating set ####
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.

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
quartz()
cyto_plot_explore(transformed_timepoint_gating_set,
                  channels_x = "FSC-A",
                  channels_y = "B2-A",
                  axes_limits = "data")

## 3.2 Gating using the entire timepoint dataset or apply an existing gating template

# note:if you already have a gating template and don't need to draw gates, then skip cyto_draw, use cyto_gatingTemplate_apply to apply the gating template.csv to your gating set
# We have used 4 gating templates, one for each strain
# cytek_gating_02_WT_112222.csv
# cytek_gating_04_LTR_112222.csv
# cytek_gating_05_112122_ars.csv
# cytek_gating_05_ALL_121522.csv

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
# Draw zero-, one-, two-or-more- copy gates
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
               overlay = c(DGY1, DGY500),
               point_col = c("gray", "green","purple") #parent color then overlay colors
)

# gating for ARS KO populations only or ALL KO populations
cyto_gate_draw(transformed_timepoint_gating_set,
parent = "Single_cells", #first color
alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
channels = c("FSC-A","B2-A"),
axes_limits = "data",
#select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
overlay = DGY1,
point_col = c("gray", "green","purple") #parent color then overlay colors
)

# gating for WT populations
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
               overlay = DGY1,
               point_col = c("gray", "green","purple") #parent color then overlay colors
                )

# gating for LTR populations
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
               overlay = DGY1,
               point_col = c("gray", "green","purple") #parent color then overlay colors
)
#### STEP 4:  Use function to perform analysis####
#A function that will
#1 Read in all the files in a folder
#2 Read in experiment details files using pData
#3 Specify experiment markers
#4 Transform gating set
#5 Apply existing gating file using cyto_gatingTemplate_apply
#6.Output stats file as .csv

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

  #6. Get cell counts and frequencies inside each gate

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

#### STEP 5:  Apply function from STEP 4 to all subdirectories ####
folders = list.dirs(samples_dir)[-1]

try(map(folders[c(1:length(folders))],analyze_all_exp, my_markers, gating_template = paste0("cytek_gating_",version_name,".csv")))

#### DATA WRANGLING ####

# Pull in all counts or freq or single cell distribution files from directory and combine into a single dataframe #

list.files(path = ".", pattern = paste0(version_name,"_counts_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_counts_all_timepoints.csv"))

list.files(path = ".", pattern = paste0(version_name,"_freq_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_freq_all_timepoints.csv"))

# Do on hpc because large files, do once for each strain.
# Don't even try to run this command on your laptop, it's a ~12GB file!!!

list.files(path = ".", pattern = paste0(version_name,"_SingleCellDistributions")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_SingleCellDistributions_all_timepoints.csv"))

# on HPC: For Loop - for each sample, subset it and write a sc_distributions_SampleName_allTimepoints.csv
for(pop in unique(sc_distr_alltimepoints$sample)) {
  print(pop)
  sc_distr_alltimepoints %>%
    filter(sample == pop) %>%
    write_csv(paste0("sc_distributions_",pop,"_all_timepoints.csv"))
}

#### Supplement Figure - Ridgeplots ####
pop_files = list.files(pattern = "sc_distributions_") %>% sort()
make_ridgeplots = function(file_name){
  pop_name = sub("sc_distributions_", "", sub("_all_timepoints.csv","", file_name))
  
  pop_data = read.csv(file_name, stringsAsFactors = T) %>%
    filter(generation < 140) %>% 
    mutate(generation = factor(generation, levels = unique(generation)))
  
  lowcell = count(pop_data, generation, sample) %>% filter(n < 70000)
  
  pop_data %>%
    anti_join(lowcell) %>%
    ggplot(aes(x = B2A_FSC, y = generation, fill = after_stat(x), height=..density..)) +
    geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
    xlab("Cell-size normalized fluorescence (a.u.)") +
    ylab("Generation") +
    ggtitle(paste0(pop_name)) +
    theme_classic() +
    #scale_x_continuous(limits=c(0.0,3), breaks = c(0, 1, 2, 3.0)) +
    scale_x_continuous(limits=c(1.25,3.0), breaks = c(1.25, 1.5, 2, 3.0)) +
    scale_y_discrete(expand = expansion(add = c(0.2, 2.0))) + #expands the graph space or else the top is cut off
    scale_fill_distiller(type = "seq", palette = "Greens", direction = 1, guide = "colourbar") + #makes it green
    theme(
      legend.position = 'none', #remove the legend
      axis.text.x = element_text(family="Arial", size = 14, color = "black"), #edit x-tick labels
      axis.text.y = element_text(family="Arial", size = 14, color = "black"),
      axis.title=element_text(size=16)
    )
  ggsave(paste0(pop_name,"_ridgeplot_scale2_082724.png"))
}

map(pop_files, make_ridgeplots) #map() applies this ridgeplot function to all 32 population.csv files

make_ridgeplots()


#### Find timepoint with lowest median normalized GFP for each of the 4 strains ####
# Calculate Median Normalized Fluorescence for each timepoint and population
# on hpc, do once 
sc_distr_alltimepoints <- read.csv(paste0(version_name,"_SingleCellDistributions_all_timepoints.csv", stringsAsFactors = T)) %>% mutate(generation = factor(generation, levels = unique(generation)))

sc_distr_alltimepoints %>%
  group_by(sample, generation) %>%
  mutate(Med_B2A_FSC = median(B2A_FSC)) %>%
  distinct(Med_B2A_FSC, .keep_all = T) %>%
  select(-FSC.A, -B2.A, -B2A_FSC) %>%
  write_csv("medians_normalized_fluor_alltimepoints.csv")

sc_distr_alltimepoints %>%
  group_by(sample, generation) %>%
  mutate(Med_B2A_FSC = median(B2A_FSC)) %>%
  distinct(Med_B2A_FSC, .keep_all = T) %>%
  select(-FSC.A, -B2.A, -B2A_FSC) %>%
  write_csv("medians_normalized_fluor_alltimepoints.csv")

cell_numbers = count %>%
  filter(Gate == "Single_cells")

norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv") %>%
  left_join(cell_numbers) %>%
  filter(Count > 70000) %>% #filter out low cells
  filter(!(generation == 231 & Type == "0_copy_ctrl")) %>%   #filter out bad controls
  filter(!(generation == 182 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 203 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 252 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 260 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 79 & Type == "2_copy_ctrl")) %>%
  filter(!(generation == 95 & Type == "2_copy_ctrl")) %>%
  filter(!(generation == 108 & Type == "2_copy_ctrl")) %>%
  filter(!(generation == 116 & Type == "2_copy_ctrl"))

norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv")

medianGFP = norm_medians %>% filter(generation <=100) %>% group_by(Description, generation) %>% summarize(median = median(Med_B2A_FSC)) %>% slice(which.min(median)) #grouped by Description and took median across each population within its genotype

#medianGFP %>% write_csv("min-median-norm-GFP_112222.csv")
#medianGFP = read_csv("min-median-norm-GFP_112222.csv")

# Description         generation    median   timepoint

#  0 copy control               58  0.403     06
#  1 copy control               29  1.74      03
#  2 copy control               29  1.99      03
#  GAP1 ARS KO                  50  1.86      05
#  GAP1 LTR + ARS KO            50  1.93      05
#  GAP1 LTR KO                  37  1.74      04
#  GAP1 WT architecture         21  1.73      02


# Now that we have made gates based on the lowest median normalized GFP timepoint, applied them, 
# and outputted freq and counts files for each gating template, 
# let's export data frames for each strain that will be input for SBI 

# freq_and_counts.csv to input
# WT       02_WT_112222
# LTR KO   04_LTR_112222
# ARS KO   05_112122_ars
# ALL KO   05_ALL_121522

#### EXPORT DATA FRAMES that will be inputs for nnSBI #### 
# Filter for only CNV gates, until generation 116, samples with >70,000 cells. 

# Wildtype architecture
freq = read_csv("02_WT_112222_freq_all_timepoints.csv")
count = read_csv("02_WT_112222_counts_all_timepoints.csv")

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

WT_df = freq_and_counts %>% 
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 WT architecture") %>% 
  filter(generation <= 116, Type == "Experimental", Gate == "two_or_more_copy") %>% 
  filter(Count>60000) %>%
  select(sample, generation, Gate, Count, Frequency, Description) %>%
  arrange(generation, sample) %>%
  select(!Gate)

ggplot(df, aes(generation, Frequency, color = sample))+
  geom_point()+
  geom_line()+
  theme_classic()


# LTR∆ - remove population LTR_2 for sbi as well as for other CNV dynamics analysis 
freq = read_csv("04_LTR_112222_freq_all_timepoints.csv")
count = read_csv("04_LTR_112222_counts_all_timepoints.csv")

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

LTR = 
  freq_and_counts %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 LTR KO") %>% 
  filter(!sample == "gap1_ltr_2")  %>% 
  filter(generation <= 116, Type == "Experimental", Gate == "two_or_more_copy") %>% 
  select(sample, generation, Gate, Count, Frequency, Description) %>%
  arrange(generation, sample) %>%
  select(!Gate)

# ARS KO
freq = read_csv(paste0("05_112122_ars_freq_all_timepoints.csv"))
count= read_csv(paste0("05_112122_ars_counts_all_timepoints.csv"))

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

ARS = freq_and_counts %>%
  filter(Count>70000) %>% 
  filter(Description == "GAP1 ARS KO") %>% #distinct(generation)
  filter(generation <= 116, Type == "Experimental", Gate == "two_or_more_copy") %>%
  select(sample, generation, Gate, Count, Frequency, Description) %>%
  arrange(generation, sample) %>% 
  select(!Gate)

# ALL KO
freq = read_csv(paste0("05_ALL_121522_freq_all_timepoints.csv")) 
count= read_csv(paste0("05_ALL_121522_counts_all_timepoints.csv")) 

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

ALL = freq_and_counts %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 LTR + ARS KO") %>%
  filter(generation <= 116, Type == "Experimental", Gate == "two_or_more_copy") %>%
  select(sample, generation, Gate, Count, Frequency, Description) %>%
  arrange(generation, sample) %>%
  select(!Gate)

SBI_input = bind_rows(WT_df,LTR, ARS, ALL) 
SBI_input %>% write_csv("SBI_input_011723.csv")

####### MERGE DATA with no filtering ####### 

# Wildtype architecture
freq = read_csv("02_WT_112222_freq_all_timepoints.csv")
count = read_csv("02_WT_112222_counts_all_timepoints.csv")
WT = count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency) %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 WT architecture")

# LTR∆
freq = read_csv("04_LTR_112222_freq_all_timepoints.csv")
count = read_csv("04_LTR_112222_counts_all_timepoints.csv")
LTR =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency) %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 LTR KO")

# ARS∆
freq = read_csv(paste0("05_112122_ars_freq_all_timepoints.csv"))
count= read_csv(paste0("05_112122_ars_counts_all_timepoints.csv"))
ARS =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency) %>%
  filter(Count>70000) %>% 
  filter(Description == "GAP1 ARS KO")

# All∆
freq = read_csv(paste0("05_ALL_121522_freq_all_timepoints.csv")) 
count= read_csv(paste0("05_ALL_121522_counts_all_timepoints.csv"))
ALLKO =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency) %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 LTR + ARS KO")

#Merge them to one dataframe
merged_samps = bind_rows(WT, LTR, ARS, ALLKO)

write_csv(merged_samps, "freq_and_counts_121622_MERGED_Experimentals_noControls.csv")


#### CLEAN UP THE DATA  #### 

# Remove the timepoints weirdly high and spikey. 
freq_and_counts = read_csv("freq_and_counts_121622_MERGED_Experimentals_noControls.csv")

# Plot the data to identify outliers 
freq_and_counts %>%
  filter(Count>70000) %>%
#         generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~factor(Description,
                     levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
  theme_classic() +
  #scale_x_continuous(breaks=seq(0,250,50)) +
  scale_x_continuous(breaks=seq(0,203,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
  )


# IQR to detect outliers within a timepoint 
x = freq_and_counts %>%
  filter(generation == 50,
         Description == "GAP1 WT architecture",
         Gate == "two_or_more_copy") %>%
        arrange(generation, sample) %>%
        select(sample, Description, generation, Gate, Parent, Count, Frequency, gating_template) 

boxplot(x$Frequency,
       ylab = "Frequency of GAP1 CNVs")
boxplot.stats(x$Frequency)$out  #the outlier Frequency value   LTR_4, g21, 17.3% freq is outlier 

#chose these timepoints by eye
weird_tp = freq_and_counts %>%
  filter(sample == "gap1_4" & Gate == "two_or_more_copy" & generation == 66 |  #use slope rule
           sample == "gap1_2" & Gate == "two_or_more_copy" & generation == 108|  #use slope rule
           sample == "gap1_ltr_2" |  
           sample == "gap1_ltr_4" & Gate == "two_or_more_copy" & generation == 21 | # 17.3%
           sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 166|
           sample == "gap1_all_5" & Gate == "two_or_more_copy" & generation == 116|
           sample == "gap1_all_6" & Gate == "two_or_more_copy" & generation == 124|
           sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 79 |
           sample == "gap1_ars_7" & Gate == "two_or_more_copy" & generation == 182
  )

clean_freq_and_counts = freq_and_counts %>%
  filter(Count>70000,
         generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(weird_tp)

# clean_freq_and_counts %>% write_csv("freq_and_counts_merged_CLEAN_121622.csv")
 clean_freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

#### Can go to Figure1.R now ####
 
#### Supplement - Ridgeplots ####
pop_files = list.files(pattern = "sc_distributions_") %>% sort()
make_ridgeplots = function(file_name){
  pop_name = sub("sc_distributions_", "", sub("_all_timepoints.csv","", file_name))
  
  pop_data = read.csv(file_name, stringsAsFactors = T) %>%
    filter(generation < 140) %>% 
    mutate(generation = factor(generation, levels = unique(generation)))
  
  lowcell = count(pop_data, generation, sample) %>% filter(n < 70000)
  
  pop_data %>%
    anti_join(lowcell) %>%
    ggplot(aes(x = B2A_FSC, y = generation, fill = after_stat(x), height=..density..)) +
    geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
    xlab("Cell-size normalized fluorescence (a.u.)") +
    ylab("Generation") +
    ggtitle(paste0(pop_name)) +
    theme_classic() +
    #scale_x_continuous(limits=c(0.0,3), breaks = c(0, 1, 2, 3.0)) +
    scale_x_continuous(limits=c(1.25,3.0), breaks = c(1.25, 1.5, 2, 3.0)) +
    scale_y_discrete(expand = expansion(add = c(0.2, 2.0))) + #expands the graph space or else the top is cut off
    scale_fill_distiller(type = "seq", palette = "Greens", direction = 1, guide = "colourbar") + #makes it green
    theme(
      legend.position = 'none', #remove the legend
      axis.text.x = element_text(family="Arial", size = 14, color = "black"), #edit x-tick labels
      axis.text.y = element_text(family="Arial", size = 14, color = "black"),
      axis.title=element_text(size=16)
    )
  ggsave(paste0(pop_name,"_ridgeplot_scale2_082724.png"))
}

map(pop_files, make_ridgeplots) #map() applies this ridgeplot function to all 32 population.csv files

make_ridgeplots()