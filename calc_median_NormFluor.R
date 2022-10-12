library(tidyverse)

version_name = "newGates_01_02_04_ars_all"

sc_distr_alltimepoints <- read.csv(paste0(version_name,"_SingleCellDistributions_all_timepoints.csv"), stringsAsFactors = T)%>% 
	mutate(generation = factor(generation, levels = unique(generation)))

#sc_distr_alltimepoints %>%
#mutate(generation = factor(generation, levels = unique(generation))) %>%
#filter(Description == "0 copy control") %>%
#write_csv(file = "sc_distributions_0copyControl_all_timepoints.csv")

#one <- sc_distr_alltimepoints %>% filter(Description == "1 copy control") %>%
  #write_csv(file = "sc_distributions_1copyControl_all_timepoints.csv") #do once

#two <- sc_distr_alltimepoints %>% filter(Description == "2 copy control") %>%
# write_csv(file = "sc_distributions_2copyControl_all_timepoints.csv") #do once

sc_distr_alltimepoints %>%
  group_by(sample, generation) %>%
  mutate(Med_B2A_FSC = median(B2A_FSC)) %>%
  distinct(Med_B2A_FSC, .keep_all = T) %>%
  select(-FSC.A, -B2.A, -B2A_FSC) %>%
  write_csv(paste0(version_name,"medians_normalized_fluor_alltimepoints.csv"))