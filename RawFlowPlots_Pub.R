# Load required packages
library(CytoExploreR)
library(tidyverse)

setwd("~/Google Drive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/FCS_by_pop")
folders = list.dirs()[-1] #select the FSC file folders in your directory

# Run this script four times, once for each strain. 

# Data storage structure: Every population has a folder. Each population folder contains all its associated .fcs files, one for each timepoint, ie 00.fcs, 01.fcs, 02.fcs.  

#### PART 1 ####
###### Make Experimental Details - Time Series Per Population Folder #### 
# Takes in old samplesheet that is one row per population and 
# makes an experimental details sheet that is one timepoint per row, and 
# one experiment details sheet per population, and outputs it in the population folder.  
make_exp_details = function(folder_name, samplesheet) {
  files = as_tibble(list.files(paste0(folder_name))) %>%
    separate(value, into = c("well", "sample"), sep = " ", remove = F) %>% 
    mutate(well = str_extract(well, "([A-Z])([0-9]){1,2}$")) %>% 
    mutate(sample = str_remove(sample, ".fcs")) %>% 
    mutate(generation = str_extract(value, ("[g]\\d+")),
           generation = generation %>% str_remove("g") %>% as.numeric()) %>% 
    rename(name = value) %>% 
    select(name, sample, generation) %>% 
    drop_na()
    
  
  all = files %>%
    left_join(read_csv(paste0("./",samplesheet)), by = c("sample" = "Sample name")) %>% 
    arrange(generation)
  
  write_csv(all, file = paste0(folder_name,"/",folder_name,"_experiment_details.csv"))
}
## Only need to run once to make files, can skip to make_plots
map(folders[8:12], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021_gap1_WT.csv")
map(folders[13:20], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021_gap1_all.csv")
map(folders[21:27], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021_gap1_LTR.csv")
map(folders[1:7], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021_gap1_ARS.csv")
####### end of making experiment details ##

#### PART 2 #####
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
#### Choose one of these four gating templates to apply (depending on which strain's flow data) 
#cyto_gatingTemplate_apply(transformed_pop_gs, gatingTemplate= "cytek_gating_02_WT_112222_toSinglets.csv")
#cyto_gatingTemplate_apply(transformed_pop_gs, gatingTemplate= "cytek_gating_04_LTR_112222.csv")
#cyto_gatingTemplate_apply(transformed_pop_gs, gatingTemplate= "cytek_gating_05_112122_ars.csv")
# cyto_gatingTemplate_apply(transformed_pop_gs, gatingTemplate = "cytek_gating_05_ALL_121522.csv")

#### Make Graphs #### in 4 plots per page 
start = seq(1,nrow(experiment_details), by = 5)
end = seq(5,nrow(experiment_details), by = 5)

  for(i in 1:length(start)){
  cyto_plot_save(paste0(pop_folder,"_t",start[i],"-",end[i],"_FlowPlot_082024.pdf"),width=6,height=9) #run cyto_plot_save(), then plot
  cyto_plot(pop_gs[start[i]:end[i]],
            parent = "Single_cells",
            axes_limits = "data",
            ylim = c(10^3, 10^6),
            channels = c("FSC-A", "GFP"),
            point_col_scale = rev(c(
              "#fde725",
              "#5ec962",
              "#21918c",
              "#3b528b",
              "#440154")),
            axes_text_size = 1.5,
            axes_label_text_size = 1.5,
            )
  }

}

# Run function across folders of the same strain. 
map(folders[8:12], make_plots) #Wildtype
map(folders[21:27], make_plots) #LTR∆
map(folders[1:7], make_plots) #ARS∆ 
map(folders[13:20], make_plots) #ALL∆
