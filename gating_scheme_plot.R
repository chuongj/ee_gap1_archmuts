#workflow.R
#Started: Sept 21, 2021
#Authors: Julie Chuong, Titir De, Grace Avecilla, David Gresham

##Install CytoExplorer package and requirements (can be skipped if already installed)
#library(BiocManager)
#install.packages("cytolib", "flowCore", "flowWorkspace", "openCyto")

##Install CytoExploreR from GitHub:
#library(devtools)
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(ggridges)
library(docstring)

setwd("~/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

# Objective: Make a gating scheme plot by plotting a flow plot
# GFP vs. FSC plot with zero-, one-, two-copy control, colored grey, green, darker green respectively. like in my reporter cartoon. 
#hex codes: grey,  #99cc67 ,  #346734

folders = list.dirs() #select the FSC file folders in your directory
folders = folders[9]

# Gate for single cells only for this plot
# gatingtemplate is "cytek_gating_WTSinglesOnly.csv"

exp_details_path = list.files(path = folders, pattern = "_experiment_details.csv", full.names = T)

timepoint_gating_set <- cyto_setup(path = paste0(folders), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

# use flowWorkspace::pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv

ordered_exp_details = pData(timepoint_gating_set) %>% left_join(experiment_details) #rerrange rows of data frame so merging is correct. ie. check that the .fcs name matches the sample name in attached metadata
for(i in 1:length(names(ordered_exp_details))){
  flowWorkspace::pData(timepoint_gating_set)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
}

cyto_details(timepoint_gating_set) %>% View() #check correct attachment of metadata

#STEP 3:  Perform gating on gating set

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
# quartz()
# cyto_plot_explore(transformed_timepoint_gating_set,
#                   channels_x = "FSC-A",
#                   channels_y = "B2-A",
#                   axes_limits = "data")

## 3.2 Apply an existing gating template to gate for single cells

cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_WTSinglesOnly.csv")

#Extract the 0,1 and 2 copy controls so we can color them:
indexes_ctr0 <- which(experiment_details$Description %in% c("0 copy control"))
transformed_timepoint_gating_set[[indexes_ctr0]] #check
DGY1 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[indexes_ctr0]] #DGY1 c(30,61,92)

ind_ctr1 <-as.numeric(which(experiment_details$Description %in% c("1 copy control")))
transformed_timepoint_gating_set[[ind_ctr1]] #check
DGY500 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[ind_ctr1]] #DGY500

indexes_ctr2 <- as.numeric(which(experiment_details$Description %in% c("2 copy control")))
transformed_timepoint_gating_set[[indexes_ctr2]]
DGY1315 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[[indexes_ctr2]] #DGY1315

# We will get the plot be pretending to draw zero-, one- two-copy gates. 
# In fact we will not. We will execute the command, the plot will appear, 
# we will not draw plots and kill the command. 
cyto_gate_draw(transformed_timepoint_gating_set,
parent = "Single_cells", #first color
alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
channels = c("FSC-A","B2-A"),
axes_limits = "data",
gatingTemplate = "cytek_gating_fake.csv",
overlay = c(DGY1, DGY500,  DGY1315),
point_col = c("#99cc67","grey",  "#99cc67","#346734") #colors
)  #hex codes: grey,  #99cc67 ,  #346734

# from the Quartz, pop up window, click R Studio > Save 
# it will save the plot as a PDF, which is what we want for the poster. 

quartz.save(file = "0124plot.png", type = "png", dpi = 300, bg = "transparent")
quartz.save(file = "0124plot2.png", type = "png", dpi = 300, bg = "#FFFFFF")
