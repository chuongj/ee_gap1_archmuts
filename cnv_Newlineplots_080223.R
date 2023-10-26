# Forked from gating_112122.R
# making new CNV lineplots because the 4 wildtype lines are under-reporting GAP1 CNVs due to broken CNV reporter subpop

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(ggridges)
library(docstring)

setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

#load colors
#wtGrays = c("#354f52","#666666","#6b705c","#414833","#999999")
wtGrays = "#6b705c"
allGolds = c("#ffba08", "#faa307", "#dda15e", "#7f5539", "#9c6644", "#fdc409", "#9c7e1e","#D9BB59")
arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631", "#f85c46", "#bb3521","#d9402a")
ltrBlues = c( "#6699cc", "#005f73", "#0a9396", "#4292C6", "#2171B5", "#3799fb", '#66b3cc', "#3a0ca3")

#load names
my_facet_names <- as_labeller(c("GAP1 WT architecture" = "Wildtype architecture",
                                "GAP1 LTR KO" = "LTR removed",
                                "GAP1 ARS KO" = "ARS removed",
                                "GAP1 LTR + ARS KO" = "LTR and ARS removed"))

#load frequency data
freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")


# proportion of CNV over time
propCNV = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 203) %>%
  dplyr::filter(!(sample == "gap1_1"  |
                    sample == "gap1_2" |
                    sample == "gap1_4" |
                    sample == "gap1_5")
                ) %>%
#exclude these controls timepoints that look weird on ridgeplots
  ggplot(aes(generation, proportion, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~factor(Description,
                     levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
 # facet_wrap(~sample) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
  theme_classic() +
#  scale_x_continuous(breaks=seq(0,250,50)) +
  scale_x_continuous(breaks=seq(0,203,50)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 25, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 25, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=18)
  )
propCNV

ggsave("propCNV_080223_8x12.pdf", bg = "#FFFFFF", height = 8, width = 12)

##############
# Lineplots per populations
pop_lineplot = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 203) %>%
  dplyr::filter(!(sample == "gap1_1"  |
                    sample == "gap1_2" |
                    sample == "gap1_4" |
                    sample == "gap1_5")
  ) %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  # facet_wrap(~factor(Description,
  #                    levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  facet_wrap(~sample) +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
  theme_classic() +
  #scale_x_continuous(breaks=seq(0,250,50)) +
  scale_x_continuous(breaks=seq(0,203,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
  )

pop_lineplot

ggsave(paste0("propCNV_pop_080223_10x14.pdf"), bg = "#FFFFFF", height = 10, width = 14)

## Now let's add median lines overlaid

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
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 203) %>%
  dplyr::filter(!(sample == "gap1_1"  |
                    sample == "gap1_2" |
                    sample == "gap1_4" |
                    sample == "gap1_5")
  ) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency)) %>% #calculate the medians per genotype per timepoint
  dplyr::filter(Description %in% c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO")) %>%
#  dplyr::filter(Description %in% c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")) %>%
#  dplyr::filter(Description == c("GAP1 WT architecture")) %>% #choose which genotypes to plot
  ggplot(aes(generation, med, color = Description)) +
#  stat_summary(fun.data = median_cl_boot, aes(generation, Frequency, color = Description), geom = "errorbar")+
  geom_point(size = 2)+
  geom_line(size = 3)+
  scale_color_manual(values=c("#6b705c", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("#6b705c", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
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

ggsave(paste0("medianPropCNV_noBars_080223_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT_080223_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT+LTR_080223_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)
ggsave(paste0("medianPropCNV_WT+LTR+ARS_080223_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)



