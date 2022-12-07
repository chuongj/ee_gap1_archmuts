# 12-5-22
# Scrap the idea of using a geom_smooth or any sort of spline for now since
# generalized additive model presented during lab meeting was not a good fit. Did not fit the raw data.
# current parameters... span = 1, formula = 'y ~ s(x, bs = "cs")'

# set up freq_and_counts data frame
freq_and_counts = read_csv("freq_and_counts_Merged_080622_all_timepoints.csv")

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

weird_early = freq_and_counts %>%
  filter(generation < 30,
         Type %in% c("Experimental", "1_copy_ctrl"),
         Description %in% c("1 copy control", "GAP1 WT architecture","GAP1 LTR KO"),
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample) %>%
  #select(-name, -`Outflow well`, -Media)
  filter(Frequency > 15)

#chose these timepoints by eye
weird_tp = freq_and_counts %>%
  filter(sample == "gap1_4" & Gate == "two_or_more_copy" & generation == 66 |
           sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 166|
           sample == "gap1_all_5" & Gate == "two_or_more_copy" & generation == 116|
           sample == "gap1_all_6" & Gate == "two_or_more_copy" & generation == 124|
           sample == "gap1_ltr_2"
  )
############################################################################
# Plot 1
# all 5 wildtype lines drawn
# lineplot of proportion of CNVs over time

wtGrays = c("gray","#666666","#CCCCCC","gray","#999999")

freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  anti_join(fails) %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  filter(Description == "GAP1 WT architecture") %>%
  ggplot(aes(generation, Frequency, color = sample)) +
 # stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_line(size = 3)+
  scale_color_manual(values=wtGrays,  #custom colors
                     limits=c("GAP1 WT architecture"),
                     labels=c("Wild type architecture"))+
  scale_fill_manual(values=wtGrays, #custom colors
                    limits=c("GAP1 WT architecture"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture"))+ #third, now you can change legend labels
  scale_x_continuous(breaks=seq(0,200,50))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25))+
  #scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  xlab("Generation")+
  ylab("Percent of cells with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

ggsave(paste0("WT-propCNV_120622_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("WT-propCNV_120622_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)


# Plot 2
# median of only the WT populations
quartz()
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  anti_join(fails) %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  filter(Description == "GAP1 WT architecture") %>%
  group_by(generation) %>%
  mutate(med = median(Frequency)) %>% #calculate the medians per genotype per timepoint
  #arrange(generation, sample) %>%
  #  select(sample, Description, generation, Frequency, med) %>% View()
  #stat_summary(aes(generation, med, color = Description) ) +
  ggplot(aes(generation, med, color = Description)) +
#  stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_line(size = 3)+
  scale_color_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
                    limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+#third, now you can change legend labels
  scale_x_continuous(breaks=seq(0,200,50))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25))+
  #scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  xlab("Generation")+
  ylab("Percent of cells with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

ggsave(paste0("WT-MEDIAN_noBars_propCNV_120622_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("WT-MEDIAN_noBars_propCNV_120622_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)

# Plot 3 - WT and LTR removed populations - median lineplots with CI bars
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  anti_join(fails) %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  filter(Description %in% c("GAP1 WT architecture" ,"GAP1 LTR KO")) %>%
  group_by(Description, generation) %>%
  mutate(med = median(Frequency)) %>% #calculate the medians per genotype per timepoint
  #arrange(generation, sample) %>%
  #  select(sample, Description, generation, Frequency, med) %>% View()
  #stat_summary(aes(generation, med, color = Description) ) +
  ggplot(aes(generation, med, color = Description)) +
#  stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_line(size = 3)+
  scale_color_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
                    limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+#third, now you can change legend labels
  scale_x_continuous(breaks=seq(0,200,50))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25))+
  #scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  xlab("Generation")+
  ylab("Percent of cells with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

ggsave(paste0("WTLTR_MEDIAN_noBars_propCNV_120622_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("WTLTR_MEDIAN_noBars_propCNV_120622_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)

# Plot 4 - WT, LTR removed,  ARS removed populations - median lineplots with CI bars
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  anti_join(fails) %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  filter(Description %in% c("GAP1 WT architecture" ,"GAP1 LTR KO", "GAP1 ARS KO")) %>%
  group_by(Description, generation) %>%
  mutate(med = median(Frequency)) %>% #calculate the medians per genotype per timepoint
  #arrange(generation, sample) %>%
  #  select(sample, Description, generation, Frequency, med) %>% View()
  #stat_summary(aes(generation, med, color = Description) ) +
  ggplot(aes(generation, med, color = Description)) +
  #stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_line(size = 3)+
  scale_color_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
                    limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+#third, now you can change legend labels
  scale_x_continuous(breaks=seq(0,200,50))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25))+
  #scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  xlab("Generation")+
  ylab("Percent of cells with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

ggsave(paste0("WTLTRARS_noBars_MEDIAN_propCNV_120622_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("WTLTRARS_noBars_MEDIAN_propCNV_120622_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)


# Plot 5
# median of all population for each genotype lineplot  with 95% confidence intervals
# median line plot code -
# https://www.datanovia.com/en/lessons/ggplot-error-bars/
# http://rstudio-pubs-static.s3.amazonaws.com/28101_41a7995107d94c8dbb07bbf7cd7e8291.html
library(elucidate)

# bootstrap the 95% confidence interval
median_cl_boot <- function(x, conf = 0.95) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t,
                                                                          uconf))
}

quartz()
freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  anti_join(fails) %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency)) %>% #calculate the medians per genotype per timepoint
  #arrange(generation, sample) %>%
#  select(sample, Description, generation, Frequency, med) %>% View()
  #stat_summary(aes(generation, med, color = Description) ) +
  ggplot(aes(generation, med, color = Description)) +
  #stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_line(size = 3)+
  scale_color_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
                    limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+#third, now you can change legend labels
  scale_x_continuous(breaks=seq(0,200,50))+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25))+
  #scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
  xlab("Generation")+
  ylab("Percent of cells with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

ggsave(paste0("medianpropCNV_noBars_120622_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianpropCNV_noBars_120622_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)
