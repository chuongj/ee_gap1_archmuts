#CytoExploreR.R
#Author David Gresham
#Started 12/02/2020
#Modified 09/01/2021

#Julie's analysis of EE_GAP1_ArchMuts2021 flow cytometry FCS data files from the Cytek Aurora 
# where we measured GFP fluorescnece of GAP1 CNV's of clonal populations of
# various GAP1 architecture mutants (LTR delete, ARS delete, double delete) were
# experimentally evolved in glutamine-limited media for ~250 generations. 

#This code uses CytoExploreR to perform an analysis of our CNV reporter data.

#load library requirements
library(ggplot2)
library(ggforce)
library(tidyverse)
library(usethis)
library(ggridges)
library(readxl)
library(devtools)
library(shiny)

##Install CytoExplorer package and requirements (can be skipped if already installed)
#library(BiocManager)
#BiocManager::install("ncdfFlow","flowViz","cytolib","ggcyto","flowCore", "flowWorkspace", "openCyto","flowStats")
#BiocManager::install("openCyto")

##Install CytoExploreR from GitHub: skip if already installed
#devtools::install_github("DillonHammill/CytoExploreRData")
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(flowCore)
library(flowWorkspace)
library(CytoExploreR) 
#library(CytoExploreRData)

#Setup experiment using cyto_setup()
#This function will read .fcs files to a cutoset which are then added to a GatingSet
# Must click Save and Close in the Viewer to generate the .csv files
#in generating the cytoset a experiment-details.csv file is created. Additional columns can be added by right clicking on a mac and inserting columns
#Secondly Experiment-Markers.csv file is created
#accuri_gating_set <- cyto_setup(path="~/Projects/Data/flowdata/Accuri", select="fcs") #David's path to data
# We want to exclude all color channels except GFP from being read in in the beginning so we don't have so much data to process -- use argument restrict=TRUE to remove unused/unmarked channels.
# why do all three gates look the same when we plot them? It's as off no restriction occured. 
timept01_gating_set <- cyto_setup(path="~/nyudrive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files/01_EE_GAP1_ArchMuts_2021_061621_g8_GA", restrict=TRUE, select="fcs")
timept03_gs <- cyto_setup(path = "~/nyudrive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files/03_EE_GAP1_ArchMuts_2021_062321_g29_IS", restrict = TRUE, select="fcs")

#To interactively edit the experiment details
#cyto_details_edit(timept01_gating_set)

#To tqke a look at the data
cyto_plot_explore(timept01_gating_set,
                  channels_x = "FSC-A", #forward scatter on x-axis
                  channels_y = c("GFP") # GFP on y-axis
) #all cells from all strains 
cyto_plot_explore(timept03_gs,
                  channels_x = "FSC-A",
                  channels_y = "GFP")

#quartz.save(type = "png", file = "first_look_plot.png")

#cyto_plot_profile(accuri_gating_set,
#                  parent = "root",
#                  channels = c("FSC-A","GFP")
#                  )
cyto_plot_profile(timept01_gating_set,
                  parent = "root",
                  channels = c("FSC-A","GFP")
) #two expression profiles ridgeplots with % of mode of y-axis, FSC-A and GFP on x-axis. we chose the channels to display. still unsure what %of mode means like why are there so many colors. what do the colors mean? 

#quartz.save(type = "png", file = "init_expression_plot.png")

#cyto_plot_gating_tree(accuri_gating_set)
#cyto_plot_gating_tree(timept01_gating_set) #error


#Transform the data using a logicle transformation.
#The first step applies the transformation using cyto_transformer_logicle()
#By default, this function does not apply the transformation to FSC and SSC, so we do that in a second step.
#we then need to combine the transformation using cyto_transformer_combined()
#finally we apply the transformation to all the data using cyto_transform()
#accuri_transformed <- cyto_transformer_logicle(accuri_gating_set)

#this can rewritten as a function and then return only transformed_time01 since we don't need first 3  (intermediate) objects
time01_transformed <- cyto_transformer_logicle(timept01_gating_set) #initial transformation- what parameters is it using to do this transformation??? biexponential. Note how the expression profile plot changes to many distributions to just one. I think the transformation caused this. 
#accuri_FSC_SSC_transformed <- cyto_transformer_logicle(accuri_gating_set,
#                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")) #note the expression profiles had one for each channel
time03_trans = cyto_transformer_logicle(timept03_gs)
time03_fsc_ssc = cyto_transformer_logicle(timept03_gs, channels = c("FSC-A","FSC-H","SSC-A","SSC-H"))
time03_combined_trans = cyto_transformer_combine(time03_trans, time03_fsc_ssc)
transformed_time03 = cyto_transform(timept03_gs, trans = time03_combined_trans)
rm(time03_trans)
rm(time03_fsc_ssc)
rm(time03_combined_trans)

time01_FSC_SSC_transformed <- cyto_transformer_logicle(timept01_gating_set,
                                                       channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H"))
#accuri_combined_transformed <- cyto_transformer_combine(accuri_transformed,
#                                  accuri_FSC_SSC_transformed)
time01_combined_transformed <- cyto_transformer_combine(time01_transformed,
                                                        time01_FSC_SSC_transformed)
#transformed_accuri <- cyto_transform(accuri_gating_set,
#                                    trans = accuri_combined_transformed)
transformed_time01 <- cyto_transform(timept01_gating_set,
                                     trans = time01_combined_transformed)
#this can be a function and then return only transformed_time01 since we don't need the intermediate objects
#rm(time01_combined_transformed)
#rm(time01_FSC_SSC_transformed)
#rm(time01_transformed)
#quartz.save(type = "png", file = "combined_and_transformed_plot.png")

#To tqke a look at the transformed data
#cyto_plot_explore(transformed_accuri[1:3],
#                  density_modal = TRUE,
#                  axes_limits = "data",
              #    display = 100000,
#                  point_col_alpha = 0.5,
#                   point_col = "black",
#                  channels_x = "FSC-A",
#                  channels_y = c("GFP")
#)
cyto_plot_explore(transformed_time01,#all samples
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

cyto_plot_explore(transformed_time01[length(transformed_time01)-1],#DGY1, zero copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                      display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot_explore(transformed_time01[1],#DGY500, one copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot_explore(transformed_time01[length(transformed_time01)],#DGY1315, two copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot_explore(transformed_time01[c(length(transformed_time01)-1, 1, length(transformed_time01))],#zero_copy control, one copy control, two copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  point_col = c("black","red","blue"),
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

cyto_plot_explore(transformed_time03, #all sample
                  density_modal = TRUE,
                  #axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot_explore(transformed_time03[length(transformed_time01)-1], #DGY1 zero copy cells
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

cyto_plot_explore(transformed_time03[1], #DGY500 one copy cells
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

cyto_plot_explore(transformed_time03[,length(transformed_time01)], #DGY1315 two copy cells
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

cyto_plot_explore(transformed_time03[c(length(transformed_time01)-1, 1, length(transformed_time01))],#zero_copy control, one copy control, two copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = c("black","red","blue"),
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

#quartz.save(type = "png", file = "transformed_time01_plot.png")

#cyto_plot_profile(transformed_accuri,
#                  parent = "root",
#                  channels = c("FSC-A","GFP"), #add as many channels as needed.
#                  legend = "fill"
#)

#####Gating cells.
#To gate cells we use the interactive function of Cytoexplorer
#The details of the gating are recorded in a .cvs file, which must be specified in the function
#The gating is done in a hierarchical manner, so that there is a parent and a child for each gate.

#Cytoexplorer will merge the entire set of fcs files and plot them when gating is performed.  As a result it is not possible to distinguish the different samples.
#However, there are two ways that individual samples can be used for the purpose of gating as described below.

#Below I describe three different ways to gate the data using the gating set.

#1.  Gating using the ENTIRE SET of experimental data.
#This is the default behavior when the gating set is called.

#The first gate defines cells based on forward scatter and side scatter
#For this gate we use all the cells
#cyto_gate_draw(transformed_accuri,  #entire gating set is plotted for gating purposes.  It is downsamples to 25,000 events for
#               parent = "root",
#               alias = "Cells",
#               channels = c("FSC-A","SSC-A"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#              )
#press esc when done drawing to exit
cyto_gate_draw(transformed_time01,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               display = 100000,
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
)

#the next gate defines the singlets based on forward scatter height and width
#note parent=cells and alias = "single_cells" so we are going down one in the hierachy
#it adds the gate to the existing .csv "Adding newly constructed gate(s) to time01_gating.csv" 
#cyto_gate_draw(transformed_accuri,
#               parent = "Cells",
#               alias = "Single_cells",
#               channels = c("FSC-A","FSC-H"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv"
#               )
cyto_gate_draw(transformed_time01,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               display = 100000,
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv"
               )
#Identify the negative control/no GFP sample and use it to gate non-fluorescent cells
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "Negative",
#               channels = c("FSC-A","GFP"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv"
#               )
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "Negative",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv"
               )

#Define the one copy GFP gate
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "One_copy",
#               channels = c("FSC-A","FL1-A"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#)
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "One_copy",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
)

#Define the two copy GFP gate
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "Two_copy",
#               channels = c("FSC-A","FL1-A"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#)
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "Two_copy",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
)

#Define the three plus copy GFP gate
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "multi_copy",
#               channels = c("FSC-A","FL1-A"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#)
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "multi_copy",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
)

#To visualize the effect of the gating on each sample
# will flash through each FCS file, draw the gates, and report % of cells inside the gate
#cyto_plot_gating_scheme(transformed_accuri)
cyto_plot_gating_scheme(transformed_time01)

#To visualize the gating tree
#cyto_plot_gating_tree(transformed_accuri)
cyto_plot_gating_tree(transformed_time01)

#To visualize the freq of cells in each gate for a single sample
#cyto_plot_gating_tree(transformed_accuri[[1]], stat="freq")
cyto_plot_gating_tree(transformed_time01[[1]], stat="freq") #DGY500 , 1 copy control
cyto_plot_gating_tree(transformed_time01[[length(transformed_time01)]], stat="freq") #DGY1315,2copy control, which is last position of the list 
cyto_plot_gating_tree(rev(transformed_time01)[[1]], stat="freq") #DGY1315, 2copy which is at last of list
cyto_plot_gating_tree(dplyr::last(transformed_time01), stat="freq") transformed_time01[[length(transformed_time01)]]
cyto_plot_gating_tree(transformed_time01[[length(transformed_time01)-1]], stat="freq") #DGY1,0copy
cyto_plot_gating_tree(transformed_time01[[length(transformed_time01)-1]], stat="count") #DGY1,0copy
cyto_plot_gating_tree(rev(transformed_time01)[[1]], stat="freq") #DGY1315, 2copy which is at last of list

#stats_freq1 <- cyto_stats_compute(transformed_accuri,
#                             parent = "Single_cells",
#                             alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                             stat="freq",
#                             save_as = "stats_freq1.csv")

stats_freq1 <- cyto_stats_compute(transformed_time01,
                                  parent = "Single_cells",
                                  alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_freq1.csv")

#stats_count1 <- cyto_stats_compute(transformed_accuri,
#                                 parent = "Single_cells",
#                                 alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                 stat="count",
#                                 save_as = "stats_count1.csv")
stats_count1 <- cyto_stats_compute(transformed_time01,
                                   parent = "Single_cells",
                                   alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                   stat="count",
                                   save_as = "stats_count1.csv")

View(stats_freq1)
View(stats_count1)

######################################################################
#2. Gating strategy uses INDIVIDUAL samples
#In this case we visualize individual samples, which are used to define the gates

#Remove all gates to start clean with no gates applied

#cyto_gate_remove(transformed_accuri,
#                 gatingTemplate = "Accuri_gating.csv",
#                 alias = "Cells")  #will remove Cells gate and all descendant gates
cyto_gate_remove(transformed_time01,
                 gatingTemplate = "time01_gating.csv",
                 alias = "Cells") #will remove Cells gate and all descendant gates
#First we gate using all the cells to define the cells and singlets as in the first strategy
#cyto_gate_draw(transformed_accuri,  #entire gating set is plotted for gating purposes.  It is downsamples to 25,000 events for
#               parent = "root",
#               alias = "Cells",
#               channels = c("FSC-A","SSC-A"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#)
cyto_gate_draw(transformed_time01,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               display = 100000,
               #type = NULL,
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv"
) 


#the next gate defines the singlets based on forward scatter height and width
#cyto_gate_draw(transformed_accuri,
#               parent = "Cells",
#               alias = "Single_cells",
#               channels = c("FSC-A","FSC-H"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv"
#)
cyto_gate_draw(transformed_time01,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv"
)
####Now we use individual samples to define the gates, which are applied to all samples in the gatingset

#Identify the negative control/no GFP sample and use it to gate non-fluorescent cells
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "Negative",
#               channels = c("FSC-A","FL1-A"),
#               axes_limits = "data",
#               select = list(Strain = "DGY1"),  #strain used to define no GFP signal
#               gatingTemplate = "Accuri_gating.csv"
#               )
#cyto_gate_edit(transformed_time01, parent="Single_cells", alias="zero_copy",gatingTemplate = "time01_gating.csv", select = list(name = "Experiment_042-Plate_001-Reference Group-B3 Unstained (Cells).fcs"))

#cyto_gate_remove(transformed_time01,parent="Single_cells", alias="zero_copy",gatingTemplate = #"time01_gating.csv")
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "zero_copy",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               select = list(name = "Experiment_042-Plate_001-Reference Group-B3 Unstained (Cells).fcs"),  #strain used to define no GFP signal
               gatingTemplate = "time01_gating.csv"
)
#cyto_plot_gating_scheme(transformed_accuri2[1:2], stat="freq")

##To overlay the negative cells when gating the single copy GFP, we need to extract the data from the relevant sample, which is the first sample [[1]] in this case
#negative <- cyto_extract(transformed_accuri, "Single_cells")[[1]]
zero_copy <- cyto_extract(transformed_time01, "Single_cells")[[length(transformed_time01)-1]] #second to last sample is DGY1 zero copy control

#Define the one copy GFP gate using the relevant control sample
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "One_copy",
#               channels = c("FSC-A","FL1-A"),
#               axes_limits = "data",
#               select = list(Strain = "DGY500"),  #strain used to define one copy of GFP
#               gatingTemplate = "Accuri_gating.csv",
#               overlay=negative  #will plot the negative cells as gray plots on the same plot
#                )
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "one_copy",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               select = list(name="Experiment_042-Plate_001-1 copy control-D3 DGY500.fcs"),  #strain used to define one copy of GFP
               gatingTemplate = "time01_gating.csv",
               overlay=zero_copy  #will plot the negative cells as gray plots on the same plot
)
#visualize the freq of cells inside gates for a selected sample
#cyto_plot_gating_scheme(transformed_accuri[2], stat="freq")
cyto_plot_gating_scheme(dplyr::last(transformed_time01),stat="freq") #last sample is DGY1315, 2copy control strain
cyto_plot_gating_scheme(transformed_time01[[1]],stat="freq") #first sample is DGY500,1copy control strain
cyto_plot_gating_scheme(rev(transformed_time01)[[2]],stat="freq") #penultimate sample is DGY1, zero copy control strain


##To overlay the one GFP copy cells when gating the two copy GFP we need to extract the relevant data
#one_copy <- cyto_extract(transformed_accuri, "Single_cells")[[2]]
one_copy <- cyto_extract(transformed_time01, "Single_cells")[[1]] #first sample is DGY500,1copy control strain

#Define the two copy GFP gate using the relevant control sample
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "Two_copy",
#               channels = c("FSC-A","FL1-A"),
#               select = list( = "DGY1315"),  #control strain used to define two copies of GFP
#             axes_limits = "data",
#              gatingTemplate = "Accuri_gating.csv",
#               overlay=one_copy  #will plot the one copy cells as gray points on the same plot #for reference
#)
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "two_copy",
               channels = c("FSC-A","GFP"),
               select = list(name= "DGY1315"),  #control strain used to define two copies of GFP
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
               overlay=one_copy  #will plot the one copy cells as gray points on the same plot for reference
)
#cyto_plot_gating_scheme(transformed_accuri[3], stat="freq")
cyto_plot_gating_scheme(rev(transformed_time01)[2], stat="freq")
##To overlay the two GFP copy cells when gating the more than two copies we need to extract the data
#two_copy <- cyto_extract(transformed_accuri, "Single_cells")[[3]]
two_copy <- cyto_extract(transformed_time01, "Single_cells")[[length(transformed_time01)]]
#cyto_extract(dplyr::last(transformed_time01), "Single_cells")

#Define the three copy GFP gate using the relevant control sample
#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = "multi_copy",
#              channels = c("FSC-A","FL1-A"),
             #  select = list(name = "DGY2158"),  #if available use a known control
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#               overlay=two_copy  #will plot the two copy cells as gray points on the same plot for reference
              )
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = "multi_copy",
               channels = c("FSC-A","GFP"),
               #  select = list(name = "DGY2158"),  #if available use a known control
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
               overlay=two_copy  #will plot the two copy cells as gray points on the same plot for reference
)

#cyto_plot_gating_scheme(transformed_accuri2, stat="freq")

#cyto_plot_gating_tree(transformed_accuri[[7]],
#                      stat="freq")
cyto_plot_gating_tree(transformed_time01[[1]],
                      stat="freq")
#cyto_plot_gating_scheme(transformed_accuri[7],
#                        back_gate = TRUE,
#                        gate_track = TRUE)
cyto_plot_gating_scheme(transformed_time01[1],
                        back_gate = TRUE,
                        gate_track = TRUE)
cyto_plot_gating_scheme(dplyr::last(transformed_time01),
                        back_gate = TRUE,
                        gate_track = TRUE)
#stats_freq2 <- cyto_stats_compute(transformed_accuri,
#                                  parent = "Single_cells",
#                                  alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                  stat="freq",
#                                  save_as = "stats_freq2.csv")
stats_freq2 <- cyto_stats_compute(transformed_time01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_freq2.csv")
#stats_count2 <- cyto_stats_compute(transformed_accuri,
#                                   parent = "Single_cells",
#                                   alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                   stat="count",
#                                   save_as = "stats_count2.csv")

View(stats_freq2)
View(stats_count2)

############################################################################ 
# 3. Drawing multiple gates simultaneously.
#It is possible to draw multiple gates on the same plot by defining the samples, extracting the data and plotting them as overlays

#Remove all gates to start clean
#cyto_gate_remove(transformed_accuri,
#                 gatingTemplate = "Accuri_gating.csv",
#                 alias = "Cells")
#cyto_gate_remove(transformed_time01,
#                 gatingTemplate = "time01_gating.csv",
#                 alias = "Cells")
cyto_gate_remove(transformed_time01,
                 gatingTemplate = "time01_gating.csv",
                 #parent = "Single_cells",
                 alias = "zero_copy")
cyto_gate_remove(transformed_time01,
                 gatingTemplate = "time01_gating.csv",
                 alias = "one_copy")
cyto_gate_remove(transformed_time01,
                 gatingTemplate = "time01_gating.csv",
                 alias = "two_copy")

######First we have to gate the cells and singlets using the entire dataset
#cyto_gate_draw(transformed_accuri,  #entire gating set is plotted for gating purposes.  It is downsamples to 25,000 events for
#               parent = "root",
#               alias = "Cells",
#               channels = c("FSC-A","SSC-A"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv",
#)
cyto_gate_draw(transformed_time01,  #entire gating set is plotted for gating purposes. display 100,000 cells. Draw gate for entire cell population.
               display = 100000,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv",
)

#the next gate defines the singlets based on forward scatter height and width
#cyto_gate_draw(transformed_accuri,
#               parent = "Cells",
#               alias = "Single_cells",
#               channels = c("FSC-A","FSC-H"),
#               axes_limits = "data",
#               gatingTemplate = "Accuri_gating.csv"
#)
cyto_gate_draw(transformed_time01,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "time01_gating.csv"
)

##To overlay the negative cells we need to extract the data from the relevant sample, which is the first sample in this case
#negative <- cyto_extract(transformed_accuri, "Single_cells")[[1]] #DGY1
zero_copy <-cyto_extract(transformed_time01, "Single_cells")[[length(transformed_time01)-1]] #penultimate is DGY1
##To overlay the one GFP copy cells  we need to extract the relevant data
#one_copy <- cyto_extract(transformed_accuri, "Single_cells")[[2]] #DGY500
one_copy <- cyto_extract(transformed_time01, "Single_cells")[[1]] #first sample is DGY500
##To overlay the two GFP copy cells  we need to extract the data
#two_copy <- cyto_extract(transformed_accuri, "Single_cells")[[3]] #DGY1315
two_copy <- cyto_extract(transformed_time01, "Single_cells")[[length(transformed_time01)]] #DGY1315

#cyto_gate_draw(transformed_accuri,
#               parent = "Single_cells",
#               alias = c("Negative", "One_copy", "Two_copy","multi_copy"), #defines gate names (4 total in this case)
#               channels = c("FSC-A","FL1-A"),
#               axes_limits = "data",
#               select = list(Strain = c("DGY1","DGY500","DGY1315","DGY2158")),  #control strains used to different copy numbers of GFP.  May not be necessary to do this.
#               gatingTemplate = "Accuri_gating_three.csv",
#               overlay = c(negative, one_copy, two_copy),  #the corresponding data for each control which has been extracted
#               point_col = c("black", "green", "red", "blue")
#               )
cyto_gate_draw(transformed_time01,
               parent = "Single_cells",
               alias = c("zero_copy", "one_copy", "two_copy","multi_copy"), #defines gate names (4 total in this case)
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               select = list(name = c("Experiment_043-Plate_001-Reference Group-B3 Unstained (Cells).fcs",
                                        "Experiment_043-Plate_001-1 copy control-D3 DGY500.fcs",
                                        "Experiment_043-Plate_001-Reference Group-F3 DGY1315 mCitrine (Cells).fcs"
                                  )),  #control strains used to different copy numbers of GFP.  May not be necessary to do this.
               gatingTemplate = "time01_gating.csv",
               overlay = c(zero_copy, one_copy, two_copy),#the corresponding data for each control which has been extracted
               #point_col = c("black", "green", "red")
)
#quartz.save(type = "png", file = "draw_3_gates_tp01.png")

#cyto_plot_gating_tree(transformed_accuri[[7]],
#                      stat="freq")

#cyto_plot_gating_scheme(transformed_accuri[7],
#                        back_gate = TRUE,
#                        gate_track = TRUE)
cyto_plot_gating_scheme(transformed_time01[length(transformed_time01)-1], #check % of DGY1 control cells that fit in the zero_copy gate
                        back_gate = TRUE,
                        gate_track = TRUE)
cyto_plot_gating_scheme(transformed_time01[1], #check % of DGY500 control cells that fit in the one_copy gate
                        back_gate = TRUE,
                        gate_track = TRUE)

cyto_plot_gating_scheme(transformed_time01[length(transformed_time01)], #check % of DGY1315 control cells that fit in the two_copy gate
                        back_gate = TRUE,
                        gate_track = TRUE)

#stats3 <- cyto_stats_compute(transformed_accuri,
#                             parent = "Single_cells",
#                             alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                             stat="freq")
stats3 <- cyto_stats_compute(transformed_time01,
                             parent = "Single_cells",
                             alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                             stat="freq")

#stats_freq3 <- cyto_stats_compute(transformed_accuri,
#                                  parent = "Single_cells",
#                                  alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                  stat="freq",
#                                  save_as = "stats_freq3.csv")
stats_freq3 <- cyto_stats_compute(transformed_time01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_freq3.csv")
#stats_count3 <- cyto_stats_compute(transformed_accuri,
#                                   parent = "Single_cells",
#                                   alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                   stat="count",
#                                   save_as = "stats_count3.csv")
stats_count3 <- cyto_stats_compute(transformed_time01,
                                   parent = "Single_cells",
                                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                   stat="count",
                                   save_as = "stats_count3.csv")
View(stats3)
View(stats_freq3)
View(stats_count3)

# we cannot read-in all FSC files all at once into one gating set. Cytoexplorer error File Size exceeded. Therefore an alternate strategy is --  draw gates for one timepoint, apply gates one timepoint at a time, export a stats table for each timepoint (named the timepoint and gen time), then finally merge all freq tables together. We want to see that the 95% of cells fall into the gate across all time points ( we can also write a code for this, add a column in the table that says if >95% write Pass, is <95% write Fail). 
# To apply the gates we drew for timepoint 1, can we simply read-in subsequent timepoint files in a loop transform them and then use cyto_stats_compute()  (parent = "root" perhaps instead of single_cells) to apply the gating and get frequency of cells inside the gate? 
# maybe the functions cyto_gate_copy() or cyto_gate_extract() may come in handy
# what's the difference between cyto_gate_extract() and cyto_extract? cyto_extract extracts a cytoframe and can be used to downstream functiond that need a cytoframe. cyto_gate_extract i think extract a gate from a gatingTemplate... maybe to be used downstream to apply the gate to another gating set hopefully that is the case because that's what I want to do.  
# I think what we can do is draw a set of gates based on ONE timepoint then save that as  Gating Template. Use that gating template to apply the gates to other timepoint(gating sets) and get freq stats using cyto_plot_gating_scheme() use argument gatingTemplate = time01_gates  and cyto_stats_compute.. Write a loop for this. 

cyto_plot_gating_scheme(transformed_time03, 
                        gatingTemplate = "time01_gating.csv",
                        back_gate = TRUE,
                        gate_track = TRUE)
#flashes through each individuals sample (unless you specify the one sample), draws gates from the template, and reports the % of cells in each gate. How to export these % numbers though? 

cyto_plot_gating_scheme(transformed_time03[1], #DGY500 sample
                        gatingTemplate = "time01_gating.csv",
                        back_gate = TRUE,
                        gate_track = TRUE)

cyto_stats_compute(transformed_time03,
                   parent = "Single_cells",
                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                   stat="freq",
                   save_as = "time03_freq.csv")
cyto_stats_compute(transformed_time01,
                   parent = "Single_cells",
                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                   stat="freq",
                   save_as = "time01_freq.csv")

#Simply try cyto_plot_gate() which Plot Gate Objects onto an Existing cyto_plot
cyto_plot_explore(transformed_time03, #all sample
                  density_modal = TRUE,
                  #axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP"))


#############Additional examples of using functions

#if you need to redraw a gate you need to use cyto_gate_edit()
cyto_gate_edit(transformed_accuri,
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               gatingTemplate = "Accuri_test_data.csv",
               type = "boundary")

###Drawing a tree of the gated samples
#cyto_plot_gating_tree(transformed_accuri,
#                      stat="freq")
cyto_plot_gating_tree(transformed_time01,
                      stat="freq")

#Draw the gating scheme for a single sample
#cyto_plot_gating_scheme(transformed_accuri[[3]],
#                        back_gate = TRUE,
#                        gate_track = TRUE)
#cyto_plot_gating_scheme(transformed_time01[[1]],
#                        back_gate = TRUE,
#                        gate_track = TRUE)
cyto_plot_gating_scheme(transformed_time01[[length(transformed_time01)-1]],
                        back_gate = TRUE,
                        gate_track = TRUE)

#cyto_stats_compute(transformed_accuri,
#                   alias = c("One copy GFP"),
#                   stat = "median",
#                   channels = "GFP")
cyto_stats_compute(transformed_time01[length(transformed_time01)-1],
                   alias = c("zero_copy"),
                   stat = "median",
                   channels = "GFP")


#############Visualizing data
#Function for visualizing the data

#cyto_plot_explore(transformed_accuri,
#                  channels_x = "FSC-A",
#                  channels_y = "GFP"
#)
cyto_plot_explore(transformed_time01[length(transformed_time01)-1],
                  channels_x = "FSC-A",
                  channels_y = "GFP",
                  display=100000,
                  ) #after you draw gates the the cyto_plot_explore plots have changed. Looks like most of points are missing (got cut out). One solution is to restart the R code or delete the transformed_time01 object and start from beginning. But that's weird. 
#quartz.save(type = "png", file = "cyto_plot_explore_after_gateDraw")
#quartz.save(type = "png", file = ".png")

#cyto_plot(transformed_accuri,
#          parent="Cells",
#          channels = "FSC-A",
     #    alias = "",
#          xlim = c(10000, 2000000),
#          density_fill = rep("blue",19),
#          density_stack = 0.7)

cyto_plot(transformed_time01, #cyto_plot() works just fine to plot data but it default plots the first sample in the gating set (not the whole set)
          parent="Cells",
          display=100000,
          channels = c("FSC-A","GFP"))

cyto_plot(transformed_time01[length(transformed_time01)], #plot last sample in gating set
          parent="Cells",
          display=100000,
          channels = c("FSC-A","GFP"))

cyto_plot(transformed_time01[length(transformed_time01), length(transformed_time01)-1], #cyto_plot() CANNOT plot more than one sample. ie) here I was trying to plot both the zero copy and two copy control samples. It will default to plot the first sample in the bracket only. 
          parent = "Cells",
          channels = c("FSC-A", "GFP"),
          #point_col = c("black", "red","blue")
          )

cyto_plot(transformed_time01, #cyto_plot() CANNOT plot more than one sample even if you use select argument to plot two samples at a time
          parent = "Cells",
          channels = c("FSC-A", "GFP"),
          select= c("Experiment_042-Plate_001-Reference Group-B3 Unstained (Cells).fcs", 
                    "Experiment_042-Plate_001-Reference Group-F3 DGY1315 mCitrine (Cells).fcs"),
          #point_col = c("black", "red","blue")
)
cyto_plot(transformed_time01[length(transformed_time01)], #try overlay argument to plot all three controls 
          parent = "Single_cells",
          display=100000,
          channels = c("FSC-A", "GFP"),
          overlay = c("zero_copy","two_copy"), #had to cyto_extract them first as on object in order to use them in the overlay argument
          point_col = c("black", "red", "blue") #first color (black) assigned to the sample, next two colors for the overlay colors respectively listed red for zero copy then blue for two copy 
)

#zero_copy <- cyto_extract(transformed_time01, "Single_cells")[[length(transformed_time01)-1]]
cyto_plot(transformed_accuri[1:19],
          parent="Cells",
          channels = "GFP",
        #  label_text = "Strain",
        legend = "line",
        xlim = c(10000, 2000000),
          density_fill = rep("green",19),
          density_stack = 0.7)


cyto_plot(transformed_accuri[2],
          parent="Cells",
          alias = "",
          channels = c("FSC-A","GFP"),
          xlim = c(100000, 3000000),
          ylim = c(10000, 3000000),
      )

