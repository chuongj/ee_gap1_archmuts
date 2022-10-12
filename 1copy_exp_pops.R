# Plot freq of one-copy experimental populations in CNV gate

# 10/11/21
# Julie chuong

# try to see the false negative rate of one-copy experimental population appearing
# in the CNV gate even though I think they only harbor 1 copy of GFP/GAP1.
# such behavior appears in the WT and LTR KO populations
# aka  not entirely flat lines before the incline

library(tidyverse)
library(ggridges)
library(docstring)

setwd("/Volumes/GoogleDrive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

freq = read_csv("01_02_04_v2_fw_freq_all_timepoints.csv")
freq2 = read_csv("newGates_01_02_04_ars_all_freq_all_timepoints.csv")

count= read_csv("01_02_04_v2_fw_counts_all_timepoints.csv")
count2=read_csv("newGates_01_02_04_ars_all_counts_all_timepoints.csv")

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

#Table of low cell observations, convenient to have to anti_join() in further steps
lowcell = freq_and_counts %>%
  filter(Count <7000) %>%
  mutate(generation = factor(generation, levels = unique(generation))) %>% #View()
  select(-Count)

## check controls are in their proper gates
fails = freq_and_counts %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(str_detect(Description, "control")) %>%
  select(Description, Strain, generation, Gate, Frequency, name, Count) %>%
  mutate(flag = case_when(Strain == "DGY1" & Gate == "zero_copy" & Frequency >= 95 ~ "pass",
                          Strain == "DGY1" & Gate == "zero_copy" & Frequency < 95 ~ "fail",
                          Strain == "DGY1" & Gate == "one_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY1" & Gate == "two_or_more_copy" & Frequency >=11 ~ "fail",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY500" & Gate == "zero_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY500" & Gate == "two_or_more_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY1315" & Gate == "zero_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY1315" & Gate == "one_copy" & Frequency >= 11 ~ "fail"
  ))%>%
  dplyr::filter(flag == "fail") %>%
  arrange(Description)
View(fails)
#fails %>% write_csv("01_02_04_v2_83_fail.csv")
#fails %>% write_csv("01_02_04_v2_fail_calc_thres_stringent_.csv")
#fails %>% write_csv("01_02_04_v2_79_10_fail_.csv")
#fails %>% write_csv("01_02_04_v2_fw_79_11_fail.csv")

###### plot proportion of control cells in control gates over time ######
# in order to asses the gates
freq_and_counts %>%
  filter(Count>70000,
         str_detect(Description, "control"),
         generation <250) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count) %>%
  dplyr::filter(!(Description == "1 copy control" & generation == 182 |
                    Description == "2 copy control" & generation == 79 |
                    Description == "2 copy control" & generation == 95 |
                    Description == "2 copy control" & generation == 108 |
                    Description == "2 copy control" & generation == 116)) %>% #exclude these controls timepoints that look weird on ridgeplots
  anti_join(fails) %>% #exclude the contaminated controls timepoints (the failed timepoints)
  ggplot(aes(generation, Frequency, color = Gate)) +
  geom_line() +
  facet_wrap(~Description) +
  ylab("% of cells in gate") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=12))


######  Plot proportion of the population with a CNV over time ######
my_facet_names <- as_labeller(c("GAP1 WT architecture" = "Wildtype architecture",
                                "GAP1 LTR KO" = "LTR KO",
                                "GAP1 ARS KO" = "ARS KO",
                                "GAP1 LTR + ARS KO" = "LTR and ARS KO"))
#colors
wtGrays = c("gray","#666666","#CCCCCC","gray","#999999")
allGolds = c("#DEBD52","#DBB741","#D7B02F","#dbb844","#D9BB59","#fdc409","#9c7e1e","#D9BB59")
arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631", "#f85c46", "#bb3521","#d9402a" )
ltrBlues = c("#6699cc", '#66b3cc',"#6BAED6" ,"#4292C6", "#2171B5","#3799fb","#3972ab","#4799eb")

propCNV = freq_and_counts %>%
  filter(Count>70000,
         generation <= 250) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  dplyr::filter(!(Description == "1 copy control" & generation == 182 |
                    Description == "2 copy control" & generation == 79 |
                    Description == "2 copy control" & generation == 95 |
                    Description == "2 copy control" & generation == 108 |
                    Description == "2 copy control" & generation == 116)) %>% #exclude these controls timepoints that look weird on ridgeplots
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~factor(Description,
                     levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
  )
propCNV

################################################################################
# plot proportion of population in each gate over time for each of 28 experimental populations
prop_plot_list = list()
i=1
for(exp in unique(freq_and_counts$Description)) {
  prop_plot_list[[i]] = freq_and_counts %>%
    filter(Count>70000) %>%
    #filter(generation != 79, generation != 116,generation != 182,generation != 252) %>%
    filter(Description==exp) %>%
    ggplot(aes(generation, Frequency, color = Gate)) +
    geom_line(size =1.5) +
    facet_wrap(~sample) +
    ylab("% of cells in gate") +
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,250,50))+
    theme(text = element_text(size=10))
  i = i+1
}
names(prop_plot_list) = unique(freq_and_counts$Description)
prop_plot_list$`GAP1 WT architecture` # change index to view replicates for different genetic backgrounds
prop_plot_list$`GAP1 ARS KO`
prop_plot_list$`GAP1 LTR KO`
prop_plot_list$`GAP1 LTR + ARS KO`
