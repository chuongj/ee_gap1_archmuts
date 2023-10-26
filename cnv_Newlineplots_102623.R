# Forked from gating_112122.R
# making new CNV lineplots because the 4 wildtype lines are under-reporting GAP1 CNVs due to broken CNV reporter subpop
# Make the order of lines showing up: Wildtype, ARS removed, LTRs removed, LTR and ARS removed. 

# Load required packages
library(tidyverse)
library(plotly)

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
  geom_line(linewidth = 2.5) +
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

#### Median CNV lineplots ####
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
  dplyr::filter(Description %in% c("GAP1 WT architecture", "GAP1 ARS KO")) %>%
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

ggsave(paste0("medianPropCNV_WT+ARS_102623_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)


#### 10/26/23 add the alpha-fade to lineplots because people always ask what the variance is between populations ####
# https://plotly.com/r/line-charts/ 
# Plotly resource Pieter showed me! Also called filled line charts 

med_freq_counts = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 203) %>%
  dplyr::filter(!(sample == "gap1_1"  |
                    sample == "gap1_2" |
                    sample == "gap1_4" |
                    sample == "gap1_5")
  ) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency))

ltrko = med_freq_counts %>% dplyr::filter(Description == "GAP1 LTR KO")

# get lows and highs (min and max per timepoint)
low = ltrko %>% group_by(generation) %>% summarize(mins = min(Frequency))
high = ltrko %>% group_by(generation) %>% summarize(maxs = max(Frequency))
medians = ltrko %>% group_by(generation) %>% summarise(median = max(med))
                                    
fig = plot_ly(high, x = ~generation, y = high$maxs, type = 'scatter', mode = 'lines', line = list(color = 'transparent'), showlegend = FALSE, name = "High")
fig
fig = fig %>% add_trace(y= ~low$mins, type = 'scatter', mode = 'lines',
                         fill = 'tonexty', fillcolor = 'rgba(203, 221, 239,0.6)', line = list(color = 'transparent'), showlegend = FALSE, name = "Low")   ##CBDDEF  #BBD3EA # color='#CBDDEF', #RGBA the 4th number is the opacity! 
fig
fig <- fig %>% add_trace(x = ~generation, y = medians$median, type = 'scatter', mode = 'lines',
                         line = list(color='#6699cc', width = 5),
                         name = 'LTR removed') 

fig <- fig %>% layout(paper_bgcolor='white', plot_bgcolor='white',
                      xaxis = list(title = "Generation",
                                 
                                   showgrid = TRUE,
                                   showline = TRUE,
                                   showticklabels = TRUE,
                                   tickcolor = 'rgb(127,127,127)',
                                   ticks = 'outside',
                                   zeroline = TRUE),
                      yaxis = list(title = "Median percent of cells with GAP1 CNV",
                                   
                                   showgrid = TRUE,
                                   showline = TRUE,
                                   showticklabels = TRUE,
                                   tickcolor = 'rgb(127,127,127)',
                                   ticks = 'outside',
                                   zeroline = FALSE))
fig

ggsave(paste0("medianCNV_fade_LTR_102623_8x16.png"), bg = "#FFFFFF", height = 8, width = 16)

