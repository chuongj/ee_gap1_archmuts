# 4-24-23
# Sanity check of what data my lineplots used, before I present committee meeting Apr 28 2023. 

library(tidyverse)
setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")
# -----------------------------------------------------
# 12-16-22
# Scrap the idea of using a geom_smooth or any sort of spline for now since
# generalized additive model presented during lab meeting was not a good fit. Did not fit the raw data.
# current parameters... span = 1, formula = 'y ~ s(x, bs = "cs")'

# Plot median population proportion of CNVs over time on a single plot

# use freq_and_counts data frame from 12-16-22
# which used genotype-based gating, ie) one custom gate for each genotype, using lowest median GFP timepoint as the gating guide AND NOT gates based on controls 
# and we removed outliers (by eye and IQR rule) --> clean_freq_and_counts
clean_freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

######
# median of all population for each genotype lineplot  with 95% confidence intervals
# median line plot code -
# https://www.datanovia.com/en/lessons/ggplot-error-bars/
# http://rstudio-pubs-static.s3.amazonaws.com/28101_41a7995107d94c8dbb07bbf7cd7e8291.html

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
clean_freq_and_counts %>% 
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
#  generation <= 116) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency)) %>% #calculate the medians per genotype per timepoint
  # filter(Description == "GAP1 WT architecture") %>%
  filter(Description %in% c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")) %>%
  ggplot(aes(generation, med, color = Description)) +
  #stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_point(size = 2)+
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
  ylab("Median percent of cells
with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

ggsave(paste0("medianPropCNV_WT_042423_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT_042423_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT+LTR_042423_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT+LTR_042423_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT+LTR+ARS_042423_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT+LTR+ARS_042423_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_042423_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_042423_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_g116_042523_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_g116_042523_8x16.pdf"), bg = "#FFFFFF", height = 8, width = 16)


##############################
# Population line plots just for the committee slides
quartz()
clean_freq_and_counts %>% 
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>% 
  filter(sample == "gap1_1") %>% 
  filter(Description %in% c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")) %>%
  ggplot(aes(generation, Frequency, color = Description)) +
  #stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_point(size = 5)+
  geom_line(size = 6)+
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
  ylab("CNV Frequency") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 50),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 50, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 50, color = "black"))
