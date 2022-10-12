#CytoExploreR.R
#Author David Gresham
#Started 12/02/2020
#Modified 09/21/2021
#Author Julie Chuong

###################################

#Julie's analysis of EE_GAP1_ArchMuts2021 flow cytometry FCS data files from the Cytek Aurora 
# where we measured GFP fluorescence of GAP1 CNV's of clonal populations of
# various GAP1 architecture mutants (LTR delete, ARS delete, double delete) that were
# experimentally evolved in glutamine-limited media for ~250 generations. We sampled every __ generations and have 27 timepoints, up to 32 samples per timepoint.  

#This code aims to gate cells for zero copy, one copy, two copy, and two plus copy and make plots of the CNV's over time. 

#This code uses CytoExploreR to perform an analysis of our CNV reporter, where GAP1 is fused to GFP. 

# This R notebook contain's David'd #3 method to draw gates simultaneously on one plot with the help of overlaying controls samples. 

#load library requirements
library(ggplot2)
library(ggforce)
library(tidyverse)
library(usethis)
library(ggridges)
library(readxl)
library(devtools)
library(shiny)
library(stringr)

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

#9-30-21 method to apply gates. 
# there are 2 ways 
# 1) use openCyto::gt_gating.gatingTemplate or 
# 2) cyto_gate_extract to extract gates to apply then 

#gt_gating(x, y, ...)
# x	   a gatingTemplate object
# y	   a GatingSet object

# 9-21-21
# David, Grace, Titir, and Julie met and agreed to a workflow and split up the coding. We will use Github to work together.

# Data Structure
#     27 timepoints = 27 folders 
#     Each timepoint has up to 32 samples. Each sample has an fcs file. Sample names are not consistent or Unique across timepoints due to Cytek default naming in the worksheet. ie) Experiment-042-Plate_001

# Step 1 - Grace. Write unique experiment-details.csv file for each timepoint.Experiment-detail.csv contain metadata about the fcs file name, sample name, generation, well, media,etc. Save the .csv file in respective timepoint folder that contains its respective FSC files. Grace can make these by writing a function that combines the sample sheet from google doc and the FCS files metadata by well number. Use the same naming convention of experiment-details.csv for stats_freq.csv and 



#### Setup experiment using cyto_setup()
#This function will read .fcs files to a cutoset which are then added to a GatingSet
# Must click Save and Close in the Viewer to generate the .csv files
#in generating the cytoset a experiment-details.csv file is created. Additional columns can be added by right clicking on a mac and inserting columns
#Secondly Experiment-Markers.csv file is created
#accuri_gating_set <- cyto_setup(path="~/Projects/Data/flowdata/Accuri", select="fcs") #David's path to data

# We want to exclude all color channels except GFP from being read in in the beginning so we don't have so much data to process. Use argument restrict=TRUE to remove unused/unmarked channels.

# Step 2) Read in timepoint 1 FSC files from its respective folder. Edit the experimental-markers.csv so that B2-A channel is marked ie) renamed to "GFP" and unused channels are not marked. This experimental-markers.csv is universal to all FSC files/timepoints.Store in the parent directory ie) /FCS_files

timept01_gating_set <- cyto_setup(path="~/nyudrive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files/01_EE_GAP1_ArchMuts_2021_061621_g8_GA", restrict=TRUE, select="fcs", details=F) #point it to timepoint 1. Set details = F because we have unique experiment-details.csv files for each timepoint and they were generated in Step 1. 
# save_as = "01_g8-Experiment-Details.csv"
    #save_as argument did not work. It still made a .csv with default naming convention.
file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #hard code rename the experiment-markers.csv file from default to whatever you want. Since this is a universal file to be used across all timepoints I gave it the experiment name. 

#To edit the markers table, use cyto_markers_edit() 

## STEP 3 and 4 (see github)

## Big Arrow Step 5 and 6 David 
#STEP 5:  Use function to perform analysis
#A function that will
#1 Read in all the files in a folder (ie, cytosetup() to make a gating set)
#2 Read in experiment details files .. which lives in the same folder the the fcs files are in. 
#3 Specify experiment markers .. which lives in the parent directory
#4 Transform gating set  .. transform_logicle()
#5 Apply existing gating file ... cyto_extract() the zero, one, two, multi copy gates
#6.Output stats file as .csv... cyto_stats_compute() applies the extracted gates and writes the stats file
#Author: David

# 5.1 -- in parallel read in files in timepoint folders and make the gating sets 
map(folders[2:4], cyto_setup, restrict=TRUE, select="fcs", details=F, markers = "EE_GAP1_ArchMuts_2021_Experiment-Markers.csv") #markers editor will still pop up and you will have to click Save & Close

big_arrow_steps = function(folders) {
 gating_set =  cyto_setup(folders, restrict=TRUE, select="fcs", details=F, markers = F) 
    #trans_gating_set = cyto_transformer_logicle(gating_set, channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"))
    return(gating_set)
}  
cyto_transformer_logicle(gating_set, channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"))
big_arrow_steps(folders[2])  

map(folders[2:4], cyto_setup, restrict=TRUE, select="fcs", details=F, markers = F)

# 5.2 -- read in experiment details file ... we have a exp_details_path variable. We generate them in STEP1. But must read them in here ? if we can. If not, we can join it with the stats.csv in STEP 7.  
# 5.3 -- "./EE_GAP1_ArchMuts_2021_Experiment-Markers.csv"
#5.4  Transform gating set 

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Grace

# my attempt here to get some practice
listFiles = list.files(path="~/nyudrive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files/", full.names = TRUE)

#For Loop to Make Gating Sets and Experiment CSV files 
      # I don't like that it asks us to 'Save & Close' the experiment-markers.csv and experment-details.csv files in each iteration of the loop...or we only need the FIRST iteration's of markers-csv files. Someway to stop it in the next iterations? -- Yes, set arguments markers = F and details = F in cyto_setup()
 
#move to parent directory of the timepoints FCS files
listFiles = list.files(path="~/nyudrive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files/", full.names = TRUE)

#remove the timepoint 1 in the list of files
listFiles=listFiles[-1]

#remove more to work with a subset of the files
listFiles=listFiles[-c(6:length(listFiles))]

#making gating sets for each timepoint
#Use apply() or purrr::map() instead of For-Loop which will run in parallel. and will be faster than a For Loop. 
# cyto_setup() detects the first (existing) experimental-details.csv from timepoint 1 and keeps it.
# I don't like that it asks us to 'Save & Close' the experiment csv files in each iteration of the loop...or we only need the FIRST iteration's csv files. Someway to stop it in the next iterations?
for (file in (listFiles)){
  #print(file)
  tp_number = str_sub(file, 93, 94) #extract timepoint number
  #print(tp_number)
  assign(paste("timept",tp_number,"_gs", sep = ""), cyto_setup(path = file,restrict = TRUE, select = "fcs", markers = F, details = F)) #make a gating set per file
}

map(listFiles, make ))

make_gsets <- function(file)
{
  tp_number = str_sub(file, 93, 94) #extract timepoint number
  assign(paste("timept",tp_number,"_gs", sep = ""), cyto_setup(path = file,restrict = TRUE, select = "fcs", markers = F, details = F)) #use cyto_setup to make gating set and name gating set object using timepoint number
}

cyto_setup(path ="~/nyudrive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files/01_EE_GAP1_ArchMuts_2021_061621_g8_GA",restrict = TRUE, select = "fcs", markers = F, details = "./01_EE_GAP1_ArchMuts_2021_experiment_details.csv")
#even if you give it a path to the detail.csv file, it won't recognize it or read it in. It will  generate and pop-up a brand new experiment-details.csv containing one name column, using the names of the file you provided in path= argument. 

#for (file in listFiles){
#  make_gsets(file)
#}

testset = lapply(listFiles, make_gsets)



#STEP 7:  Combine stats_freq.csv files into a single dataframe
#Pull in all stats_freq files from directories and assemble into a single dataframe
#Author: Julie






#To interactively edit the experiment details
# cyto_details_edit(timept01_gs)



#####  FIRST LOOKS  #############################
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

gsets = ls(pattern = "timept+.*_gs") #put all the gating sets in a container called gsets

# For Loop/Apply instead of For Loop - Function to Transform Each Gating Sets
      #return the transformed gating set. we do not need the intermediates. 


######Transform the data using a logicle transformation.
#The first step applies the transformation using cyto_transformer_logicle()
#By default, this function does not apply the transformation to FSC and SSC, so we do that in a second step.
#we then need to combine the transformation using cyto_transformer_combined()
#finally we apply the transformation to all the data using cyto_transform()
#accuri_transformed <- cyto_transformer_logicle(accuri_gating_set)

#this can rewritten as a function and then return only transformed_time01 since we don't need first 3  (intermediate) objects

time03_trans = cyto_transformer_logicle(timept03_gs)
time03_fsc_ssc = cyto_transformer_logicle(timept03_gs, channels = c("FSC-A","FSC-H","SSC-A","SSC-H"))
time03_combined_trans = cyto_transformer_combine(time03_trans, time03_fsc_ssc)
transformed_time03 = cyto_transform(timept03_gs, trans = time03_combined_trans)
rm(time03_trans)
rm(time03_fsc_ssc)
rm(time03_combined_trans)

fsc_ssc_transform_gs = function(gating_set) {
  initial_trans = cyto_transformer_logicle(gating_set)
  fsc_ssc_trans = cyto_transformer_logicle(gating_set, channels = c("FSC-A","FSC-H","SSC-A","SSC-H"))
  combine_the_two = cyto_transformer_combine(initial_trans, fsc_ssc_trans)
  final_transformed = cyto_transform(gating_set, trans = combine_the_two)
  return(final_transformed)
}
fsc_ssc_transform_gs(timept01_gating_set)
time01_transformed <- cyto_transformer_logicle(timept01_gating_set) #initial transformation- what parameters is it using to do this transformation??? biexponential. Note how the expression profile plot changes to many distributions to just one. I think the transformation caused this. 
#accuri_FSC_SSC_transformed <- cyto_transformer_logicle(accuri_gating_set,
#                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")) #note the expression profiles had one for each channel
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

cyto_transformer_logicle(timept01_gating_set,
                         channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H","B2-A"))

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

cyto_plot_explore(transformed_g79[length(transformed_g79)-1],#DGY1, zero copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                      display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot_explore(transformed_g79[1],#DGY500, one copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot_explore(transformed_g79[length(transformed_g79)],#DGY1315, two copy control
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
#However, there are two ways that individual samples can be used for the purpose of gating as described below in method 2 and method 3. 

# Julie using method 3 only
############################################################################ 
# 3. Drawing multiple gates simultaneously.
#It is possible to draw multiple gates on the same plot by defining the samples, extracting the data and plotting them as overlays

#Remove all gates to start clean
#cyto_gate_remove(transformed_accuri,
#                 gatingTemplate = "Accuri_gating.csv",
#                 alias = "Cells")
cyto_gate_remove(transformed_time01,
                 gatingTemplate = "time01_gating.csv",
                 alias = "Cells")
#if it doesn't work use aliases further down the hierarchy like so:
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
##To overlay the one GFP copy cells  we need to extract the data from the relevant sample
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
cyto_plot_gating_scheme(transformed_timept01[length(transformed_timept01)-1], #check % of DGY1 control cells that fit in the zero_copy gate
                        #back_gate = TRUE,
                        #gate_track = TRUE
                        )
cyto_plot_gating_scheme(transformed_timept01[1], #check % most of DGY500 control cells that fit in the one_copy gate
                        #back_gate = TRUE,
                        #gate_track = TRUE
                        )

cyto_plot_gating_scheme(transformed_timept01[length(transformed_timept01)], #check % of DGY1315 control cells that fit in the two_copy gate
                        #back_gate = TRUE,
                        #gate_track = TRUE
                        )

cyto_plot_gating_scheme(transformed_timept01[length(transformed_timept01)],gatingTemplate = "cytek_gating_JC_v3.csv")
#stats3 <- cyto_stats_compute(transformed_accuri,
#                             parent = "Single_cells",
#                             alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                             stat="freq")
stats_time01 <- cyto_stats_compute(transformed_time01,
                             parent = "Single_cells",
                             alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                             stat="freq")

#stats_freq3 <- cyto_stats_compute(transformed_accuri,
#                                  parent = "Single_cells",
#                                  alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                  stat="freq",
#                                  save_as = "stats_freq3.csv")
freq_time01 <- cyto_stats_compute(transformed_time01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "time01_freq.csv")
#stats_count3 <- cyto_stats_compute(transformed_accuri,
#                                   parent = "Single_cells",
#                                   alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
#                                   stat="count",
#                                   save_as = "stats_count3.csv")
count_time01 <- cyto_stats_compute(transformed_time01,
                                   parent = "Single_cells",
                                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                   stat="count",
                                   save_as = "time01_count.csv")
View(stats_time01)
View(freq_time01)
View(count_time01)

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

cyto_plot_gating_scheme(transformed_time03[2], #GAP1_all_1 sample
                        gatingTemplate = "time01_gating.csv",
                        back_gate = TRUE,
                        gate_track = TRUE)

freq_time03 = cyto_stats_compute(transformed_time03,
                   parent = "Single_cells",
                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"), #it computes stats from the alias argument
                   stat="freq",
                   save_as = "time03_freq.csv") #these numbers don't match the %numbers from the cyto_plot_gating_scheme() they should be the same because 


cyto_stats_compute(transformed_time01,
                   parent = "Single_cells",
                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                   stat="freq",
                   save_as = "time01_freq.csv")

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

## Nov 1, 2021
# More troubleshooting, drawing gates and checking that the 80% of control cells fall into the drawn gate to constitute a good gate. However, some failing controls timepoint could fail because the controls could be contaminated or because the gate is bad. Plotting some controls to visually see if they are contaminated. 

#### results
# g179 2 copy control is contaminated that cells that fluoresece at zero copy. which explains why it failed
# g116 
# g252

library(CytoExploreR)
library(tidyverse)
setwd("/Volumes/GoogleDrive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_files") #Julie's WD

#folders of data
folders = list.dirs()[-1]

#markers
my_markers<-c("GFP") #list your marker name(s)
channel<-c("B2-A") #list your channel(s)
names(my_markers)<-channel

#function
analyze_all_exp = function(folder_name, my_markers, gating_template="cytek_gating.csv") {
  
  my_path <- folder_name #gets relative path name for folder to be analyzed
  
  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name
  
  my_expt_details_path <- paste0(my_path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory
  
  #1. read in files and make a gating set
  timepoint_gating_set <- cyto_setup(path=my_path, select="fcs", details=F, markers = F)
  
  #2. read in experiment details for that gating set
  my_experiment_details <- read_csv(my_expt_details_path) #import experiment-details.csv
  flowWorkspace::pData(timepoint_gating_set)$name<-my_experiment_details$name
  flowWorkspace::pData(timepoint_gating_set)$sample<-my_experiment_details$sample
  flowWorkspace::pData(timepoint_gating_set)$`Outflow Well`<-my_experiment_details$`Outflow well`
  flowWorkspace::pData(timepoint_gating_set)$Media<-my_experiment_details$Media
  flowWorkspace::pData(timepoint_gating_set)$Strain<-my_experiment_details$Strain
  flowWorkspace::pData(timepoint_gating_set)$Type<-my_experiment_details$Type
  flowWorkspace::pData(timepoint_gating_set)$Description<-my_experiment_details$Description
  flowWorkspace::pData(timepoint_gating_set)$generation<-my_experiment_details$generation
  
  #3. specify markers for that gating set
  markernames(timepoint_gating_set)<-my_markers
  
  #4. transform data
  timepoint_gating_set_transformed <- cyto_transformer_log(timepoint_gating_set,
                                                           channels =c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #transforms but returns the gating set as a list
  transformed_timepoint_gating_set<- cyto_transform(timepoint_gating_set,
                                                    trans = timepoint_gating_set_transformed)
  
  #5. apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= gating_template)
  
  #6. write stats: freq file for % of cells inside each gate, median FSC and GFP for each population, median FSC and GFP for each gated population
  #Titir
  stats_freq <- cyto_stats_compute(transformed_timepoint_gating_set,
                                   parent = c("Single_cells"),
                                   alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                   stat="freq",
                                   save_as = paste0("g252_stats_freq_",prefix,".csv") #writes to working directory
  )
#  stats_median_overall <- cyto_stats_compute(transformed_timepoint_gating_set,
#                                             parent = c("Single_cells"),
#                                             alias  = c("Single_cells"),
#                                             channels = c("FSC-A", "B2-A"),
#                                             stat="median",
#                                             save_as = paste0("v3_stats_median_overall_", prefix,".csv"))
  
  stats_cell_number <- cyto_stats_compute(transformed_timepoint_gating_set,
                                          parent = c("Single_cells"),
                                          alias  = c("Single_cells"),
                                          #channels = c("FSC-A", "B2-A"),
                                          stat="count",
                                          save_as = paste0("g252_stats_cell_number_", prefix,".csv"))
  
#  stats_median_gatewise <- cyto_stats_compute(transformed_timepoint_gating_set,
 #                                             parent = c("Single_cells"),
  #                                            alias  = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
#                                              channels = c("FSC-A", "B2-A"),
#                                              stat="median",
#                                              save_as = paste0("v3_stats_median_gatewise_", prefix,".csv"))
  transformed_timepoint_gating_set<<-return(transformed_timepoint_gating_set)
}



## Function to input and transform gating set
analyze_all_exp = function(folder_name, my_markers, gating_template="cytek_gating.csv") {
  
  my_path <- folder_name #gets relative path name for folder to be analyzed
  
  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name
  
  my_expt_details_path <- paste0(my_path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory
  
  #1. read in files and make a gating set
  timepoint_gating_set <- cyto_setup(path=my_path, select="fcs", details=F, markers = F)
  
  #2. read in experiment details for that gating set
  my_experiment_details <- read_csv(my_expt_details_path) #import experiment-details.csv
  flowWorkspace::pData(timepoint_gating_set)$name<-my_experiment_details$name
  flowWorkspace::pData(timepoint_gating_set)$sample<-my_experiment_details$sample
  flowWorkspace::pData(timepoint_gating_set)$`Outflow Well`<-my_experiment_details$`Outflow well`
  flowWorkspace::pData(timepoint_gating_set)$Media<-my_experiment_details$Media
  flowWorkspace::pData(timepoint_gating_set)$Strain<-my_experiment_details$Strain
  flowWorkspace::pData(timepoint_gating_set)$Type<-my_experiment_details$Type
  flowWorkspace::pData(timepoint_gating_set)$Description<-my_experiment_details$Description
  flowWorkspace::pData(timepoint_gating_set)$generation<-my_experiment_details$generation
  
  #3. specify markers for that gating set
  markernames(timepoint_gating_set)<-my_markers
  
  #4. transform data
  timepoint_gating_set_transformed <- cyto_transformer_log(timepoint_gating_set,
                                                           channels =c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #transforms but returns the gating set as a list
  transformed_timepoint_gating_set<- cyto_transform(timepoint_gating_set,
                                                    trans = timepoint_gating_set_transformed)
  
}

#apply function to g116 t g252
analyze_all_exp(folders[23], my_markers = my_markers, gating_template = "cytek_gating_JC_v2.csv")

folders[23]

# read in frequency csv, median csvs for all timepoints
freq = read_csv("g116_stats_freq_12_EE_GAP1_ArchMuts_2021.csv") %>% rename(Gate = Population)
cell_numbers = read_csv("g116_stats_cell_number_12_EE_GAP1_ArchMuts_2021.csv")

# add cell number column to freq table
freq_w_count = freq %>%
  #select(Description, Strain, generation, Gate, Frequency, name) %>%
left_join(cell_numbers) %>%
  select(-Marker)
View(freq_w_count)

# exclude any well/timepoint with less than 70,000 single cells
freq_w_count %>%
 dplyr::filter(Count>70000) %>%
  # check controls are in their proper gates
  dplyr::filter(str_detect(Description, "control")) %>%
  select(Description, Strain, generation, Gate, Frequency, name, Count) %>%
  mutate(flag = case_when(Strain == "DGY1" & Gate == "zero_copy" & Frequency >= 87 ~ "pass",
                          Strain == "DGY1" & Gate == "zero_copy" & Frequency < 87 ~ "fail",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency >= 87 ~ "pass",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency < 87 ~ "fail",
                          Strain == "DGY1315" & Gate == "two_copy" & Frequency >= 87 ~ "pass",
                          Strain == "DGY1315" & Gate == "two_copy" & Frequency < 87 ~ "fail"
  ))%>%
  dplyr::filter(flag == "fail") %>%
  arrange(Description) %>%
  View()

#Plot the two copy control

cyto_plot_explore(transformed_timepoint_gating_set[length(transformed_timepoint_gating_set)],#DGY1315, two copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
quartz.save("g116_DGY1315_2copycontrol.png")

#plot 1 copy control for g252
cyto_plot_explore(transformed_timepoint_gating_set[1],#DGY500, one copy control
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)

cyto_plot_explore(transformed_timepoint_gating_set[length(transformed_timepoint_gating_set)],
                  density_modal = TRUE,
                  axes_limits = "data",
                  display = 100000,
                  point_col_alpha = 0.5,
                  #point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)
cyto_plot(transformed_timepoint_gating_set[length(transformed_timepoint_gating_set)], #plot last sample in gating set
          parent="Cells",
          display=100000,
          channels = c("FSC-A","GFP"))

                  