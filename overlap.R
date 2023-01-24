# calculate median or average % overlap between controls 

library(CytoExploreR)
library(tidyverse)

setwd("/Volumes/GoogleDrive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_controls")



folders = list.dirs()[-1]

timepoint <- cyto_setup(path = ".", restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

exp_details_path = list.files(path = folders[2], pattern = "_experiment_details.csv", full.names = T)
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
ordered_exp_details = pData(timepoint) %>% left_join(experiment_details)
for(i in 1:length(names(ordered_exp_details))){
  flowWorkspace::pData(timepoint)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
}

cyto_details(timepoint) %>% View() #check correct attachment of metadata

##### Transform #####
GFP_trans <- cyto_transformer_logicle(timepoint,
                                      channels = c("B2-A"),
                                      widthBasis = -10
)
FSC_SSC_trans <- cyto_transformer_log(timepoint,
                                      channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")
)
combined_trans <- cyto_transformer_combine(GFP_trans,FSC_SSC_trans)
transformed_timepoint <- cyto_transform(timepoint,
                                                   trans = combined_trans)

## quickly check the transformation by plotting the data
#my_samples = which(experiment_details$Description %in% c("0 copy control","1 copy control","2 copy control","GAP1 LTR KO", "GAP1 WT architecture"))
quartz()
cyto_plot_explore(transformed_timepoint,
                  channels_x = "FSC-A",
                  channels_y = "B2-A",
                  axes_limits = "data")


##### Draw Gates Around 0,1,2 to See Overlap % #### 
version_name="tp2"

#First we gate for the cells
cyto_gate_draw(transformed_timepoint,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)


#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)

#Gating for CNVs using the 0,1 and 2 copy controls:
indexes_ctr0 <- which(experiment_details$Description %in% c("0 copy control"))
transformed_timepoint[[indexes_ctr0]] #check
DGY1 <- cyto_extract(transformed_timepoint, "Single_cells")[[indexes_ctr0]] #DGY1 c(30,61,92)

ind_ctr1 <-as.numeric(which(experiment_details$Description %in% c("1 copy control")))
transformed_timepoint[[ind_ctr1]] #check
DGY500 <- cyto_extract(transformed_timepoint, "Single_cells")[[ind_ctr1]] #DGY500

indexes_ctr2 <- as.numeric(which(experiment_details$Description %in% c("2 copy control")))
transformed_timepoint[[indexes_ctr2]]
DGY1315 <- cyto_extract(transformed_timepoint, "Single_cells")[[indexes_ctr2]] #DGY1315

cyto_gate_draw(transformed_timepoint,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
               overlay = c(DGY1, DGY500, DGY1315),
               #overlay = c(DGY1, DGY500),
               point_col = c("gray", "green","purple","black") #parent color then overlay colors
)

###### Found a 3 copy strain ~~ #####

setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/0B5kOS1TfN-RpfnMwZV9KQjFGb0Z4aFJIRkYxNDMyeUV1RlV6RjVwU1NMdEo4WDJrLVhLdkU/Gresham Lab_Share/Ministat Evolution 2021/Clones_g125_071522_JCMandIS/overlap2_3_010323")

##### Draw Gates Around 2 3  to See Overlap % #### 
version_name="2v3"

#First we gate for the cells
cyto_gate_draw(transformed_timepoint,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)


#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)

transformed_timepoint[[1]] #check
DGY500 <- cyto_extract(transformed_timepoint, "Single_cells")[[1]]
indexes_ctr2 <- as.numeric(which(experiment_details$Description %in% c("2 copy control")))

transformed_timepoint[[2]]
DGY1315 <- cyto_extract(transformed_timepoint, "Single_cells")[[2]]

transformed_timepoint[[3]]
clone <- cyto_extract(transformed_timepoint, "Single_cells")[[3]]

cyto_gate_draw(transformed_timepoint,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
               #overlay = c(DGY500, DGY1315, clone),
               #overlay = c(DGY500, DGY1315),
               overlay = c(DGY500),
               point_col = c("gray", "green","purple","black") #parent color then overlay colors
)

### overlap
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/0B5kOS1TfN-RpfnMwZV9KQjFGb0Z4aFJIRkYxNDMyeUV1RlV6RjVwU1NMdEo4WDJrLVhLdkU/Gresham Lab_Share/Ministat Evolution 2021/EE_GAP1_ArchMuts_2021_Clones_g125_071522_JCMandIS/overlap")
