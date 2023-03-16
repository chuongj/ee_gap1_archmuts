# Purpose: Look at the raw data for the Wildtype evolution experiments
# I think this would be best displayed as the 2D plots (FSC vs GFP)  
# for each time point in a time series. If David and I could look at 
# them together we can make sure we understand the sources of
# variation in the data and think about how to explain the difference
# in results with prior experiments.

# Started Jan 10, 2023
# Modified Jan 26, 2023

# Load required packages
library(CytoExploreR)
library(tidyverse)

setwd("~/Google Drive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_WT")
folders = list.dirs()[-1] #select the FSC file folders in your directory

input_dir = folders[1] 
#make_plots = function(input_dir){
print(input_dir)
exp_details_path = list.files(path = paste0(input_dir), pattern = "_experiment_details.csv", full.names = T)

# tp = input_dir %>% substr(3,4)
# tp_name = paste("t", tp, sep ="")

gs_01 <- cyto_setup(path = paste(folders[1]), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close
# assign(tp_name, timepoint_gating_set)
# get(tp_name)
gs_02 <- cyto_setup(path = paste(folders[2]), restrict=TRUE, select="fcs", details=F)

# use flowWorkspace::pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path, show_col_types = F) #import experiment-details.csv

#ordered_exp_details = pData(timepoint_gating_set) %>% left_join(experiment_details)
ordered_exp_details = pData(gs_01) %>% left_join(experiment_details)
#ordered_exp_details = pData(timepoint_gating_set) %>% left_join(experiment_details, by = character()) #rerrange rows of data frame so merging is correct. ie. check that the .fcs name matches the sample name in attached metadata
# for(i in 1:length(names(ordered_exp_details))){
#   flowWorkspace::pData(timepoint_gating_set)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
# }
for(i in 1:length(names(ordered_exp_details))){
  flowWorkspace::pData(gs_01)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
}

cyto_details(gs_01) %>% View() #check correct attachment of metadata

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

## Plot the transformed data
# cyto_plot_save(paste0(input_dir,"_ungated_raw_plot.png"))
# cyto_plot_explore(transformed_timepoint_gating_set,
#                   channels_x = "FSC-A",
#                   channels_y = "B2-A",
#                   axes_limits = "data",
#                   title = paste0(input_dir))
#}

#map(folders[1:4], make_plots)
## 3.2 Gating using the entire timepoint dataset or apply an existing gating template

# apply gating template used for Wildtype populations. Gate for singlets. 
# Is it the same Cells gate and Single_cells gate used in the previous
# analysis to generate the prop_CNV lineplots? YES. 
cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_02_WT_112222_toSinglets.csv")

# Extract each population so we can color them distinctly
# pop1 <- which(experiment_details$sample %in% c("gap1_1"))
# transformed_timepoint_gating_set[[pop1]] #check
# pop1 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[pop1]]

# ind2 <-as.numeric(which(experiment_details$sample %in% c("gap1_2")))
# transformed_timepoint_gating_set[[ind2]] #check
# pop2 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[ind2]]

ind3 <-as.numeric(which(experiment_details$sample %in% c("gap1_3")))
transformed_timepoint_gating_set[[ind3]] #check
pop3 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[ind3]]

ind4 <-as.numeric(which(experiment_details$sample %in% c("gap1_4")))
transformed_timepoint_gating_set[[ind4]] #check
pop4 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[ind4]]

ind5 <-as.numeric(which(experiment_details$sample %in% c("gap1_5")))
transformed_timepoint_gating_set[[ind5]] #check
pop5 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[ind5]]

#quasi-facet plots with 2 overlays, in order, very bad data viz.
#not real facet plots because they don't truly subset the populations
#too big to fit in the plots pane, so use quartz()
# quartz()
# cyto_plot(transformed_timepoint_gating_set, 
#           parent = "Single_cells", 
#           axes_limits = "data",
#           channels = c("FSC-A", "GFP"),
#         #  overlay = c(pop1, pop2, pop3, pop4, pop5),
#          #point_col = c("red", "purple"),
#      #    legend = TRUE,
#       #   point_col = c("gray","purple","green", "magenta", "brown") 
#           )
cyto_plot_save(paste0(input_dir,"_color_pop.png"))
cyto_plot_explore(transformed_timepoint_gating_set, 
          parent = "Single_cells", 
          axes_limits = "data",
          channels_x = "FSC-A",
          channels_y = "GFP", 
      #   channels = c("FSC-A", "GFP"),
       #   overlay = c(pop1, pop2, pop3, pop4, pop5),
      #overlay = c(pop1, pop3, pop4, pop5),
      overlay = c(pop3, pop4, pop5),
          #  point_col = c("red", "purple"),
          title = input_dir,
          point_col = c("gray","purple","green", "magenta", "brown") #parent color then overlay colors
)
# }
#map(folders[22:length(folders)], make_plots)

