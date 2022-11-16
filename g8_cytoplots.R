#### Graph WT pops only across 7 timepoints ####

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(ggridges)
library(docstring)

setwd("/Volumes/GoogleDrive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

folders = list.dirs()[c(9:32)]

folder_name = folders[1]

my_markers<-c("GFP") #list your marker name(s)
channel<-c("B2-A") #list your channel(s)
names(my_markers)<-channel


make_cyto_plots = function(folder_name, my_markers) {

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
  experiment_details <- read_csv(exp_details_path, show_col_types = F) #import experiment-details.csv
  for(i in 1:length(names(experiment_details))){
    print(experiment_details[i])
    #flowWorkspace::pData(timepoint_gating_set)[names(experiment_details[i])]<-experiment_details[i]  #edit this For Loop so that it searches for a match in the name before merging!
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
# Gate for Cells
  cyto_gate_draw(transformed_timepoint_gating_set,
                 parent = "root",
                 alias = "Cells",
                 channels = c("FSC-A","SSC-A"),
                 axes_limits = "data",
                 gatingTemplate = "cytek_gating_WTonly.csv")
  )
# Gate for Single Cells
  cyto_gate_draw(transformed_timepoint_gating_set,
                 parent = "Cells",
                 alias = "Single_cells",
                 channels = c("FSC-A","FSC-H"),
                 axes_limits = "data",
                 gatingTemplate = "cytek_gating_WTonly.csv")
  )
#  cyto_details(transformed_timepoint_gating_set) %>% View()
  indexes = which(str_detect(cyto_details(transformed_timepoint_gating_set)$Description,"WT"))
 # indexes = which(str_detect(experiment_details$Description,"WT"))
  cyto_plot_explore(transformed_timepoint_gating_set[indexes],
                    channels_x = "FSC-A",
                    channels_y = "B2-A",
                    axes_limits = "data")
### WHY are there TWO distributions of the WT pops at generation 8?! Something is not right!
  #let's overlay each population and see what's going on.
  # you know what!!!! those might be doublets!!!!!!!!!
  # I need to gate for singlets!!!!!

pop_names = experiment_details[indexes, ] %>% select(sample) %>% pull()
assign(pop_names[1],cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[22]]) #"gap1_2"
assign(pop_names[2],cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[23]]) # gap1_4
assign(pop_names[3],cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[26]]) # gap1_1
assign(pop_names[4],cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[27]]) # gap1_3
assign(pop_names[5],cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[28]]) #gap1_5

  # Point colour palette
cyto_plot(transformed_timepoint_gating_set[[indexes]],
            parent = "Single_cells", #first color
            overlay = c(gap1_1, gap1_2, gap1_3, gap1_4, gap1_5),
            channels = c("FSC-A","B2-A"),
            axes_limits = "data",
            point_col = "gray", # base layer - override density gradient
            point_cols = c("orange", # upper layers
                            "yellow",
                            "green",
                            "blue",
                            "purple"),
           #legend = T,
           title = "Wildtype Genome Architecture - All 5 Populations - Generation 8")

  cyto_gate_draw(transformed_timepoint_gating_set,
                 parent = "Single_cells", #first color
                 alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
                 channels = c("FSC-A","B2-A"),
                 axes_limits = "data",
                 #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
                 gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
                 #               overlay = c(DGY1, DGY500, DGY1315),
                 #overlay = c(DGY1, DGY500, exp_1_copy, DGY1315),
                 #point_col = c("gray", "green", "black", "purple","blue") #parent color then overlay colors
  )


map(folders[], make_cyto_plots)
