# Chuong et al. 2024
# DNA replication errors are a major source of adaptive gene amplification
# Figure 1 and Supplementary Information 

# Load required packages
library(chromoMap)
library(CytoExploreR)
library(tidyverse)
library(plotly)
library(scales)

# Figure 1A ####
chr11_feats = read.delim("Chromosome_XI_features.txt")
unique(chr11_feats$Feature.Type)

#format annotation file in the format for chromomap  
my_feats = chr11_feats %>% 
  filter(Feature.Type %in% c("telomere", "ARS", "tRNA gene", "long terminal repeat", "centromere") |
           Feature == "GAP1") %>% 
  separate(Coordinates,c("start", "end"), "-") %>% 
  mutate(name = rep("XI",50 )) %>%
  select(Feature, name, start, end, Feature.Type)
colnames(my_feats) <- NULL

#### Make a chromosome object as needed for chromoMap #
chrXI = data.frame(name = "XI",
                   start = 1,
                   end = 666816,
                   cent = 440129
)
colnames(chrXI)<-NULL

chromoMap(list(chrXI),
          list(my_feats),
          chr_color = c("#ede0d4"),
          #chr_width = 15,
          #chr_length = 10,
          #labels = T,
          y_chr_scale = 17, #bring ruler closer to chromosome 
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#F99B1C", "#B767AA", "#60C9DF", "#ede0d4","white", "black"))) #tRNA, ars, ltr, centromere, GAP1 ORF, telomere

# zoom in
zoom = my_feats[39:44,]
chrXI_zoom = data.frame(name = "XI", 
                        start = 512000, 
                        end = 520000)
colnames(chrXI_zoom) <- NULL
chromoMap(list(chrXI_zoom),
          list(zoom), 
          chr_width = 17,
          chr_length = 3,
          labels = T,
          chr_color = c("#ede0d4"),
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#F99B1C", "#60C9DF","white","#B767AA")))

# Figure 1B #####
# Step 1: Find timepoint with lowest median normalized GFP for each of the 4 strains #####
norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv")

medianGFP = norm_medians %>% filter(generation <=100) %>% group_by(Description, generation) %>% summarize(median = median(Med_B2A_FSC)) %>% slice(which.min(median)) #grouped by strain and took median across each population within strain

# Data structure. Have directories (one to many) containing raw FSC files.
# Also make a gating directory, which is **ONE** directory that contains ALL the FSC files you want to overlay to drawing gates. 
# Read in the names of those directories (data directories and one gating directory) here:
folders = list.dirs()

# Choose a name to be used for all output files including the gating template and associated flow data and graphs.
version_name = "version_name"

# STEP 2: Generate experiment details file from folder and FCS file names ####
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
    filter(!is.na(sample))
  
  all = files %>%
    left_join(read_csv(paste0("./",samplesheet)), by = c("sample" = "Sample name")) %>%
    mutate(generation = as.numeric(generation))
  
  write_csv(all, file = paste0(folder_name,"/",pref,"_experiment_details.csv"))
  
}

#only needs to be run once to make experiment details file
map(folders[1:length(folders)], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv")

# NOTE!!!:S kip Step 2-4 if you want to use our gating template.csv file and apply it to raw flow data. Proceed to Step 5.

# STEP 3: Read in all files in a directory and rename the channels. ####
#A directory contains an FCS file for each population.
#Results in 1 timepoint gating set containing all .fcs files, associated experiment details, and marker details

# here, load in my gating directory.
# It's only data from timepoint 3 to load in.
# Timepoint 3 was chosen because it had the lowest median fluorescence
# these data will guide us on drawing gates.
gating_dir = folders[3] #change this folder for your gating directory
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

# STEP 4:  Perform gating on gating set ####
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.

## 4.1 transform the data
# to choose different transformation see: https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html
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
cyto_plot_explore(transformed_timepoint_gating_set,
                  channels_x = "FSC-A",
                  channels_y = "B2-A",
                  axes_limits = "data")

## 4.2 Gating using the entire timepoint dataset or apply an existing gating template

# note:if you already have a gating templates and don't need to draw gates, then skip cyto_draw, use cyto_gatingTemplate_apply to apply the gating template.csv to your gating set. Note that we have one gating template per strain to a total of 4. 

cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating.csv")

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

#Gating for CNVs using the 0,1 and 2 copy controls:
indexes_ctr0 <- which(experiment_details$Description %in% c("0 copy control"))
transformed_timepoint_gating_set[[indexes_ctr0]] #check
DGY1 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[indexes_ctr0]] #DGY1 c(30,61,92)

ind_ctr1 <-as.numeric(which(experiment_details$Description %in% c("1 copy control")))
transformed_timepoint_gating_set[[ind_ctr1]] #check
DGY500 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[ind_ctr1]] #DGY500

indexes_ctr2 <- as.numeric(which(experiment_details$Description %in% c("2 copy control")))
transformed_timepoint_gating_set[[indexes_ctr2]]
DGY1315 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[indexes_ctr2]] #DGY1315

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

# STEP 5  - Use function to perform analysis ####
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
  
  exp_details_path <- paste0(path,"/",prefix,"_experiment_details.csv") 
  
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

# STEP 6:  Apply function from STEP 5 to strain subdirectories ####
# Repeat for each strain 
samples_dir = file.path("../FCS_LTR") # LTR∆ folder, for example
folders = list.dirs(samples_dir)[-1]
sample_folders = list.dirs(samples_dir)[-1] 

try(map(folders[c(1:length(folders))],analyze_all_exp, my_markers, gating_template = paste0("cytek_gating_",version_name,".csv")))
try(map(sample_folders[1:length(sample_folders)],analyze_all_exp, my_markers, gating_template = paste0("cytek_gating_",version_name,".csv")))

# STEP 7: Pull in all counts or freq or single cell distribution files from directory and combine into a single dataframe 
# Repeat for each strain
list.files(path = ".", pattern = paste0(version_name,"_counts_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_counts_all_timepoints.csv"))

list.files(path = ".", pattern = paste0(version_name,"_freq_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_freq_all_timepoints.csv"))

#STEP 8: Merge and clean up data

# Now that we have made gates for based on the lowest median normalized GFP timepoint
# and outputted freq and counts files for each gating template, let's merge them for the geno's we want.
# as export data frames for each genotype that will be input for SBI 

# freq_and_counts.csv to input
# WT       02_WT_112222
# LTR KO   04_LTR_112222
# ARS KO   05_112122_ars
# ALL KO   05_ALL_121522

#### Export data frames that will be inputs for SBI #### 
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
  #filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(Description == "GAP1 WT architecture") %>% 
  filter(generation <= 116, Type == "Experimental", Gate == "two_or_more_copy") %>% 
  filter(Count>60000) %>%
  select(sample, generation, Gate, Count, Frequency, Description) %>%
  arrange(generation, sample) %>%
  select(!Gate)

# LTR KO - remove population LTR_2 for sbi as well as for other CNV dynamics analysis 
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

####### Merge all experimental for each 4 genotypes into one dataframe with no filtering ####### 

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

# LTR KO
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

# ARS KO
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

# All KO
freq = read_csv(paste0("05_ALL_121522_freq_all_timepoints.csv")) #Do not use T03, use T05 
count= read_csv(paste0("05_ALL_121522_counts_all_timepoints.csv")) #Do not use T03, use T05 
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

# Plot the data to identify outliers 
merged_samps %>%
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

#chose these timepoints by eye
weird_tp = merged_samps %>%
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

freq_and_counts = merged_samps %>%
  filter(Count>70000,
         generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(weird_tp)

freq_and_counts %>% write_csv("freq_and_counts_merged_CLEAN_121622.csv")

#### Figure 1B after gating #####
#load CNV frequency data 
freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

#prep data for plot
med_freq_counts = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 166) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency),
         mad = mad(proportion),
         IQR = IQR(proportion))

#plot
med_freq_counts%>%
  filter(generation <= 137) %>%
ggplot(aes(x = generation, group = Description)) +
  geom_line(aes(y = med/100, color = Description), linewidth = 3) + 
  geom_ribbon(aes(y = med/100, ymin = med/100 - mad, ymax = med/100 + mad, fill = Description),alpha=0.3)+
  scale_color_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
                    limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+#third, now you can change legend labels
  scale_x_continuous(limits=c(0,125))+
  xlab("Generation")+
  ylab("Median proportion of cells
with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

# Figure 1C ####

# Define frequency threshold for CNV appearance 
  #intuition: inflection generation before it goes up (vertically) 
freq_and_counts %>%
  filter(Description == "GAP1 WT architecture",
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample)

thresh = 10  # Defining it at 10% seems to be the best to capture this inflection point 

# CNV Appearance 
Tup_per_pop_10 = 
  freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thresh) %>%
  select(Type, Strain, Description, sample, generation, Gate, Frequency) %>%
  group_by(sample) %>%
  slice(which.min(generation))

Tup_per_pop_10$Description <- factor(Tup_per_pop_10$Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO", "GAP1 LTR + ARS KO")) #reorder boxplots

ggplot(Tup_per_pop_10, aes(Description, generation, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("gray", "#6699cc", "#e26d5c", "#DEBD52"))+ #change order of colors
  ylab("   Generation of first CNV appearance") +
  #scale_x_discrete(labels=c("Wildtype architecture","LTR removed","ARS removed","LTR and ARS removed"))+
  scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, max(Tup_per_pop_10$generation)))+
  theme_classic() +
  theme(legend.position = "none",
        #axis.text.x = element_text(size = 16, color = "black"), #edit x-tick labels
        axis.text.x = element_blank(), #remove x-tick labels 
        axis.ticks.x=element_blank(), #remove x-ticks 
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 18),
        text = element_text(size=16))+
  geom_jitter(size = 2, alpha = 0.8, 
              color = c(rep("black", 5),  #wildtype, 5, gray
                        rep("#D9BB59", 8),  #LTR and ARS gold
                        rep("#e26d5c", 7),  #ARS, 7, salmon 
                        rep("#6699cc", 7)  #LTR,7, #blue
              ))

shapiro.test(Tup_per_pop_10$generation) #not normal

# Instead of ANOVA, do Krusal-Wallis test (non-parametric) 
kruskal.test(generation~Description, data = Tup_per_pop_10)

# Instead of pairwise t-tests, do pairwise Wilcoxon Mann-Whitney with Bonferroni correction
pairwise.wilcox.test(Tup_per_pop_10$generation, Tup_per_pop_10$Description, p.adjust.method = "bonferroni")

# Figure 1D #####
# Calculate CNV selection # 
#Compute natural log proportion of each population with CNV relative to that without CNV
ln_table = freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency),
         prop_NoCNV = 100-prop_CNV,
         CNV_NoCNV = prop_CNV/prop_NoCNV,
         logECNV_NoCNV = log(CNV_NoCNV))

pop_list = unique(ln_table$sample) %>% sort()
gens = unique(ln_table$generation)

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(summary(x)[4]$coefficients[2], digits = 4),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(slope == b~~~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));
}

#function to apply for loop to each of 28 populations
sliding_fit = function(num_fitpoints, population){
  timepoints = nrow(subset(ln_table, sample %in% c(population)))
  rounds = timepoints - num_fitpoints +1
  m <- matrix(ncol = 7, nrow = rounds) #nrow = number of iterations. number of iterations depend on the number of generations and the number of fitpoints. max num of generations = 24. minimum num of fitpoints is 2. therefore nrow max is 23.
  colnames(m) <- c("start", "end", "gen_start", "gen_end", "slope", "rsquared", "sample")
  start = 1
  end = num_fitpoints
  for (i in 1:rounds){
    print(i)
    pop_data <- subset(ln_table, sample %in% c(population))
    fit_points_df <- subset(ln_table, sample %in% c(population) & generation >= gens[start] & generation <= gens[end])
    if (is.na(gens[end]) == TRUE ){
      break
    }
    fit <- lm(logECNV_NoCNV ~ generation, fit_points_df) #linear model, lm(y~x, by the data)
    print(summary(fit))
    
    ggplot(pop_data, aes(x=generation,y=(as.numeric(logECNV_NoCNV)), colour=sample)) +
      geom_point() +
      geom_smooth(data=fit_points_df, method=lm, show.legend=FALSE) +
      scale_y_continuous(expand = c(0, 0), 'ln(Prop. CNV/Prop. non-CNV)', limits = c(min(pop_data$logECNV_NoCNV)-1, max(pop_data$logECNV_NoCNV)+1)) +
      annotate("text", x = 200, y = min(pop_data$logECNV_NoCNV)-0.5, label = equation(fit), parse = TRUE) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5), "Generations", limits=c(0,260)) +
      theme_classic() +
      scale_color_manual(values = c('black')) +
      guides(colour = guide_legend(override.aes = list(size=2))) +
      theme(legend.position = c(.15,.95), plot.title = element_text(size=14, hjust = 0.5), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
    
    ggsave(paste0(population,"_Sup_g",gens[start],"-",gens[end],"_",num_fitpoints,"pts.png"), width = 8, height = 5)
    
    print(paste0("From timepoints ",start," to ", end, ", generations ", gens[start], " to ", gens[end],", slope was ", as.numeric(coef(fit)[2]) %>% round(4)," and rsquared was ", as.numeric(summary(fit)[8]) %>% round(2) ))  
    
    m[i,1] <- start
    m[i, 2]<- end
    m[i, 3] <- gens[start]
    m[i,4] <- gens[end]
    m[i, 5] <- as.numeric(coef(fit)[2]) %>% round(4)
    m[i,6] <- as.numeric(summary(fit)[8]) %>% round(2)
    m[i, 7] <- population
    
    start = start + 1
    end = end + 1
  }
  m = m %>% na.omit() #remove NAs
  write_csv(as.data.frame(m), paste0(population,"_fits","_",num_fitpoints,"pts.csv"))
  return(m)
}

# call the function for all pops in the pop_list using map()
map(.x = pop_list[1:27], ~sliding_fit(4, .x))

# Pull in the fits tables per population and merged into one
slopes = list.files(path = ".", pattern = paste0("_fits","_",4,"pts")) %>%
  read_csv() %>%
  write_csv(file = "Sup_fits_4_pts_all_pops.csv")

# CNV Selection Boxplot 
meta = Tup_per_pop_10 %>% select(Description, sample)

Sup = slopes %>%
  right_join(meta) %>%
  group_by(sample) %>%
  mutate(slope = max(slope)) %>%
  ungroup() %>% 
  select(sample, slope, Description) %>% 
  distinct()

ggplot(Sup, aes(Description, slope, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("gray", "#6699cc", "#e26d5c", "#DEBD52"))+ #change order of colors
  ylab("Percent of increase in 
  CNVs per generation") +
  #scale_x_discrete(labels=c("Wildtype architecture","LTR removed","ARS removed","LTR and ARS removed"))+
  scale_y_continuous(labels = scales::label_number(scale = 100))+
  theme_classic() +
  theme(#plot.margin = unit(c(.5, .5, .5, .5), "cm"),
    legend.position = "none",
    #axis.text.x = element_text(size = 16, color = "black"), #edit x-tick labels
    axis.text.x = element_blank(), #remove x-tick labels 
    axis.ticks.x=element_blank(), #remove x-ticks 
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, vjust=2),
    text = element_text(size=16))+
  geom_jitter(size = 3, alpha = 0.8, 
              color = c(rep("black", 5),  #wildtype, 5, gray
                        rep("#D9BB59", 8),  #LTR and ARS gold
                        rep("#e26d5c", 7),  #ARS, 7, salmon 
                        rep("#6699cc", 7)  #LTR,7, #blue
              ))

hist(Sup$slope)
shapiro.test(Sup$slope) #normal distribution
Sup_anova = aov(slope~Description, data = Sup)
summary(Sup_anova) # p = 0.00318
pairwise.t.test(Sup$slope, Sup$Description, p.adjust.method = "bonferroni")

# Figure 1E #####
# Boxplot - CNV equilibrium, the other inflection point where the line pleateaus #
gen_maint = slopes %>% 
  right_join(meta) %>%
  group_by(sample) %>%
  filter(gen_start > 50, slope < 0.005) %>%
  slice(which.min(gen_start)) %>% 
  select(-start, -end)

ggplot(gen_maint, aes(Description, gen_start, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("gray", "#6699cc", "#e26d5c", "#DEBD52"))+ #change order of colors
  ylab("Generation of 
  CNV maintainence") +
  #scale_x_discrete(labels=c("Wildtype architecture","LTR removed","ARS removed","LTR and ARS removed"))+
  #scale_y_continuous(labels = label_number(scale = 100))+
  theme_classic() +
  theme(#plot.margin = unit(c(.5, .5, .5, .5), "cm"),
    legend.position = "none",
    #axis.text.x = element_text(size = 16, color = "black"), #edit x-tick labels
    axis.text.x = element_blank(), #remove x-tick labels 
    axis.ticks.x=element_blank(), #remove x-ticks 
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 18, vjust=2),
    text = element_text(size=16))+
  geom_jitter(size = 3, alpha = 0.8, 
              color = c(rep("black", 4),  #wildtype, 5, gray
                        rep("#D9BB59", 8),  #LTR and ARS gold
                        rep("#e26d5c", 6),  #ARS, 7, salmon 
                        rep("#6699cc", 7)  #LTR,7, #blue
              ))

hist(gen_maint$gen_start)
shapiro.test(gen_maint$gen_start) #normal
summary(aov(gen_start~Description, gen_maint))
pairwise.t.test(gen_maint$gen_start, gen_maint$Description, p.adjust.method = "bonferroni")
