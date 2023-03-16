# Purpose: Look at the raw data for the Wildtype evolution experiments
# I think this would be best displayed as the 2D plots (FSC vs GFP)  
# for each time point in a time series. If David and I could look at 
# them together we can make sure we understand the sources of
# variation in the data and think about how to explain the difference
# in results with prior experiments.

# Started Jan 10, 2023
# Modified Feb 6, 2023

# Load required packages
library(CytoExploreR)
library(tidyverse)

setwd("~/Google Drive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_by_pop")
folders = list.dirs()[-1] #select the FSC file folders in your directory
folders = folders[6:13] #get the LTR+ARS folders
input_dir = folders[1] 
folder_name = input_dir
###### Make Experimental Details - Time Series Per Population Folder #### 
make_exp_details = function(folder_name, samplesheet) {
  files = as_tibble(list.files(paste0(folder_name))) %>%
    separate(value, into = c("well", "sample"), sep = " ", remove = F) %>%
    mutate(well = str_extract(well, "([A-Z])([0-9]){1,2}$")) %>%
    mutate(sample = str_remove(sample, ".fcs")) %>%
    mutate(generation = str_extract(value, ("[g]\\d+")),
           generation = generation %>% str_remove("g") %>% as.numeric()) %>%
    rename(name = value) %>% 
    select(name, sample, generation) %>% 
    filter(!is.na(sample))
  
  all = files %>%
    left_join(read_csv(paste0("./",samplesheet)), by = c("sample" = "Sample name")) %>% 
    arrange(generation)
  
  write_csv(all, file = paste0(folder_name,"/",folder_name,"_experiment_details.csv"))
}
## Only need to run once to make files
map(folders, make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021_gap1_all.csv")
############ end of making experiment details ##

#### Function to Make Flow Plot - One Population Over time ####
make_plots = function(pop_folder){
print(pop_folder)
exp_details_path = list.files(path = paste0(pop_folder), pattern = "_experiment_details.csv", full.names = T)

#### Read in flow data and attach metadata using CytoexploreR
pop_gs <- cyto_setup(path = paste(pop_folder), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

# use flowWorkspace::pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path, show_col_types = F) 
ordered_exp_details = pData(pop_gs) %>% left_join(experiment_details)
#ordered_exp_details = pData(timepoint_gating_set) %>% left_join(experiment_details, by = character()) 
for(i in 1:length(names(ordered_exp_details))){
  flowWorkspace::pData(pop_gs)[names(ordered_exp_details[i])]<-ordered_exp_details[i]
}
  cyto_details(pop_gs) %>% View() #check correct attachment of metadata

#### Transform Data
GFP_trans <- cyto_transformer_logicle(pop_gs,
                                      channels = c("B2-A"),
                                      widthBasis = -10
)
FSC_SSC_trans <- cyto_transformer_log(pop_gs,
                                      channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")
)
combined_trans <- cyto_transformer_combine(GFP_trans,FSC_SSC_trans)
transformed_pop_gs <- cyto_transform(pop_gs,trans = combined_trans)

#### Gate for singlets 
# It is the same Cells gate and Single_cells gate used in the previous
# analysis to generate the prop_CNV lineplots. 
cyto_gatingTemplate_apply(transformed_pop_gs, gatingTemplate= "cytek_gating_02_WT_112222_toSinglets.csv")

#### Make Graphs #### in 4 plots per page 
start = seq(1,nrow(experiment_details), by = 5)
end = seq(5,nrow(experiment_details), by = 5)

  for(i in 1:length(start)){
  cyto_plot_save(paste0(pop_folder,"_t",start[i],"-",end[i],"_FlowPlot.pdf"),width=6,height=9) #run   cyto_plot_save(), then plot
  cyto_plot(pop_gs[start[i]:end[i]],
            parent = "Single_cells",
            axes_limits = "data",
            ylim = c(10^3, 10^6),
            channels = c("FSC-A", "GFP"),
            axes_text_size = 1.5,
            axes_label_text_size = 1.5,
            )
  }

}

# Run function across folders
map(folders[2:8], make_plots)

