library(tidyverse)
library(ggridges)

#Generation 8 Ridgeplots
sc_distributions_g8 <- read.csv("01_02_04_v2_SingleCellDistributions_01_EE_GAP1_ArchMuts_2021.csv", stringsAsFactors = T)

sc_distributions_g8 = sc_distributions_g8 %>% mutate(
  generation = factor(generation, levels = unique(sc_distributions_g8$generation)),
  sample = factor(sample, levels = unique(rev(c("ctrl0",
                                                "ctrl1",
                                                "ctrl2",
                                                "gap1_1",
                                                "gap1_2",
                                                "gap1_3",
                                                "gap1_4",
                                                "gap1_5",
                                                "gap1_ltr_1",
                                                "gap1_ltr_2",
                                                "gap1_ltr_3",
                                                "gap1_ltr_4",
                                                "gap1_ltr_5",
                                                "gap1_ltr_6",
                                                "gap1_ltr_7",
                                                "gap1_ltr_8",
                                                "gap1_ars_1",
                                                "gap1_ars_3",
                                                "gap1_ars_4",
                                                "gap1_ars_5",
                                                "gap1_ars_6",
                                                "gap1_ars_7",
                                                "gap1_ars_8",
                                                "gap1_all_1",
                                                "gap1_all_2",
                                                "gap1_all_3",
                                                "gap1_all_4",
                                                "gap1_all_5",
                                                "gap1_all_6",
                                                "gap1_all_7",
                                                "gap1_all_8"
  )
  ))),
  name = factor(name, levels = unique(rev(c("Experiment_042-Plate_001-Reference Group-B3 Unstained (Cells).fcs",
                                            "Experiment_042-Plate_001-1 copy control-D3 DGY500.fcs",
                                            "Experiment_042-Plate_001-Reference Group-F3 DGY1315 mCitrine (Cells).fcs",
                                            "Experiment_042-Plate_001-Experimental-H3 gap1_1.fcs",
                                            "Experiment_042-Plate_001-Experimental-G4 gap1_2.fcs",
                                            "Experiment_042-Plate_001-Experimental-H5 gap1_3.fcs",
                                            "Experiment_042-Plate_001-Experimental-G6 gap1_4.fcs",
                                            "Experiment_042-Plate_001-Experimental-H7 gap1_5.fcs",
                                            "Experiment_042-Plate_001-Experimental-C4 gap1_ltr_1.fcs",
                                            "Experiment_042-Plate_001-Experimental-D5 gap1_ltr_2.fcs",
                                            "Experiment_042-Plate_001-Experimental-C6 gap1_ltr_3.fcs",
                                            "Experiment_042-Plate_001-Experimental-D7 gap1_ltr_4.fcs",
                                            "Experiment_042-Plate_001-Experimental-C8 gap1_ltr_5.fcs",
                                            "Experiment_042-Plate_001-Experimental-B9 gap1_ltr_6.fcs",
                                            "Experiment_042-Plate_001-Experimental-H9 gap1_ltr_7.fcs",
                                            "Experiment_042-Plate_001-Experimental-E10 gap1_ltr_8.fcs",
                                            "Experiment_042-Plate_001-Experimental-D9 gap1_ars_6.fcs",
                                            "Experiment_042-Plate_001-Experimental-E4 gap1_ars_1.fcs",
                                            "Experiment_042-Plate_001-Experimental-E6 gap1_ars_3.fcs",
                                            "Experiment_042-Plate_001-Experimental-F7 gap1_ars_4.fcs",
                                            "Experiment_042-Plate_001-Experimental-E8 gap1_ars_5.fcs",
                                            "Experiment_042-Plate_001-Experimental-D9 gap1_ars_6.fcs",
                                            "Experiment_042-Plate_001-Experimental-A10 gap1_ars_7.fcs",
                                            "Experiment_042-Plate_001-Experimental-G10 gap1_ars_8.fcs",
                                            "Experiment_042-Plate_001-Experimental-A4 gap1_all_1.fcs",
                                            "Experiment_042-Plate_001-Experimental-A6 gap1_all_3.fcs",
                                            "Experiment_042-Plate_001-Experimental-B5 gap1_all_2.fcs",
                                            "Experiment_042-Plate_001-Experimental-B7 gap1_all_4.fcs",
                                            "Experiment_042-Plate_001-Experimental-A8 gap1_all_5.fcs",
                                            "Experiment_042-Plate_001-Experimental-G8 gap1_all_6.fcs",
                                            "Experiment_042-Plate_001-Experimental-F9 gap1_all_7.fcs",
                                            "Experiment_042-Plate_001-Experimental-C10 gap1_all_8.fcs"
  )))))

ggplot(sc_distributions_g8, aes(x = B2A_FSC, y = sample, fill = Description)) +
  geom_density_ridges(scale=1.0, quantile_lines = F, quantiles = 2) +
  xlab("Normalized fluorescence (a.u.)") +
  ylab("Sample") +
  ggtitle("Generation 8") +
  theme_classic() +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) + #expands the graph space or else the top is cut off
  #scale_fill_discrete(breaks=c("0 copy control",
  #                             "1 copy control",
  #                             "2 copy control",
  #                             "GAP1 WT architecture",
  #                             "GAP1 LTR KO",
  #                             "GAP1 ARS KO",
  #                             "GAP1 LTR + ARS KO") #change order of legend items
  #                   )+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(4, "Greens")[-1],"#DEBD52","#54DE79","#DE54B9","#5474DE"),
                    breaks=c("0 copy control",
                             "1 copy control",
                             "2 copy control",
                             "GAP1 WT architecture",
                             "GAP1 LTR KO",
                             "GAP1 ARS KO",
                             "GAP1 LTR + ARS KO")) +
  #scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52")) +  #c(ARS,LTR+ARS,LTR,WT)+
  theme(
    legend.text = element_text(family="Arial", size = 12),#edit legend text font and size
    legend.title = element_blank(), #remove legend title
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
ggsave("generation8_ridgeplots_niceColors_scale1.png")

#Ridgeplot - combine populations/replicates
sc_distributions_g8 %>%
  mutate(Description = factor(Description, levels = unique(rev(c("0 copy control",
                                                                 "1 copy control",
                                                                 "2 copy control",
                                                                 "GAP1 WT architecture",
                                                                 "GAP1 LTR KO",
                                                                 "GAP1 ARS KO",
                                                                 "GAP1 LTR + ARS KO"))))) %>%
  ggplot(aes(x = B2A_FSC, y = Description, fill = Description)) +
  geom_density_ridges(scale=1.0) +
  #scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) + #expands the graph space or else the top is cut off
  xlab("Normalized fluorescence (a.u.)") +
  ylab("Genotype") +
  ggtitle("Generation 8") +
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(4, "Greens")[-1],"#DEBD52","#54DE79","#DE54B9","#5474DE"),
                    breaks=c("0 copy control",
                             "1 copy control",
                             "2 copy control",
                             "GAP1 WT architecture",
                             "GAP1 LTR KO",
                             "GAP1 ARS KO",
                             "GAP1 LTR + ARS KO")) +
  #scale_fill_discrete(breaks=c("0 copy control", #change order of legend items
  #                             "1 copy control",
  #                             "2 copy control",
  #                             "GAP1 WT architecture",
  #                             "GAP1 LTR KO",
  #                             "GAP1 ARS KO",
  #                             "GAP1 LTR + ARS KO")) +
  theme_classic()+
  theme(
    legend.text = element_text(family="Arial", size = 12),#edit legend text font and size
    legend.title = element_blank(), #remove legend title
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
ggsave("Generation8_ridgeplots_byGenotype_scale1.png")
