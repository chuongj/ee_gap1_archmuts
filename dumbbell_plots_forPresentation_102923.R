#####this script is for estimating breakpoints and then making dumbbell plots using read depth####
# This version subsets ~25 clones representative of all mechanisms we have inferred so far  
#G Avecilla (but code mostly stolen from S Lauer)
#Jan 2019
# Modified Julie Chuong 12-19-22

## mCitrine is located 514934 - 515650
## coordinates of GAP1 CDS for our samples are 518204 - 520012
## LTR11 coordinates are 513532:513846, LTR12 is ~ 520925

#Set working directory and samplesheet
setwd("/Users/juliechuong/Lab_docs/cloneseq1_temp")

#load libraries
library(tidyverse)
library(docstring)
library(ggalt)
library(gridExtra)

# set color pallete
ars_color = "#e26d5c"
#arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631","#f85c46", "#bb3521","#d9402a")
wt_color = "black"
#wtGrays = c("#354f52", "#666666", "#6b705c", "#414833" ,"#999999")
all_color = "#DEBD52"
#allGolds=c("#ffba08", "#faa307", "#dda15e", "#7f5539", "#9c6644", "#fdc409", "#9c7e1e", "#D9BB59")
ltr_color = "#6699cc"
#ltrBlues = c("#6699cc", "#005f73", "#0a9396", "#4292C6", "#2171B5", "#3799fb", "#66b3cc","#3a0ca3")

####### STEP 7 - Make Dumbbell plots
# Dumbbell plots show the genomic position and span of CNV with a horizontal line, like a handle of dumbbell weight and two dots at either end are the CNV breakpoints, like a ends of dumbbell weight, like this: ()=========()

##### example from githib ggalt #####
#library(hrbrthemes)
#df <- data.frame(trt=LETTERS[1:5], l=c(20, 40, 10, 30, 50), r=c(70, 50, 30, 60, 80))
#ggplot(df, aes(y=trt, x=l, xend=r)) +
#  geom_dumbbell(size=3, color="#e3e2e1",
#                colour_x = "#5b8124", colour_xend = "#bad744",
#                dot_guide=TRUE, dot_guide_size=0.25) +
#  labs(x=NULL, y=NULL, title="ggplot2 geom_dumbbell with dot guide") +
#  theme_ipsum_rc(grid="X") +
#  theme(panel.grid.major.x=element_line(size=0.05))
###

# import clone data
clonedata = read_csv("cloneseq0and1_Breaks_Copies.csv")

# subset clones of interest
clonedata = clonedata %>% filter(sample %in% c(2694, 2697, 2839, 2843, 2927, 2929, 
                                   2691, 2852, 2854, 2855, 2945, 2946, 2948,
                                   2682, 2684, 2685, 2687, 2689, 2849, 
                                   2676, 2677, 2678, 2681, 2680))

#### Dumbbell Plot - Grouped by Genotype  ####
db_geno=
  clonedata %>%
  filter(gap1_rounded>1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=8, dot_guide=TRUE, dot_guide_size=0.25)+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670000),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 36, color = "black"),
        axis.text.y= element_blank(),
        axis.title = element_text(size = 30, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text=element_text(size=30, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
  #facet_grid(~gap1_rounded)
db_geno

ggsave(filename = "selectClones_dumbbell_plot_102923.png", width = 8, height = 10, bg = "white")
ggsave(filename = "selectClones_dumbbell_plot_102923.pdf", width = 8, height = 10, bg = "white") 

# Make dumbbell plot VERSION 2 
# Include aneuploidy clones and zoom into 400-660kb region
plot2 = 
  clonedata %>%
  filter(gap1_rounded>1) %>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=8, dot_guide=TRUE, dot_guide_size=0.25)+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),limits=c(400000,670191),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        #axis.text.y= element_blank(),
        axis.title = element_text(size = 26, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text=element_text(size=30, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
plot2

ggsave(filename = "selectClones_ZOOM2_dumbbell_plot_102923.png", width = 8, height = 10, bg = "white")
ggsave(filename = "selectClones_ZOOM2_dumbbell_plot_102923.pdf", width = 8, height = 10, bg = "white") 

### Barplot grouped by genotype, then sort by copy number low to high ####
bar = clonedata %>%
  filter(gap1_rounded>1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorders the strain order after arrange()
  ggplot(aes(clone, gap1_rounded, fill = Description)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                    limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  coord_flip()+
  theme_classic() +
  ylab("GAP1 copy number") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 30, color = "black"),
              #axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(size = 26, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.text=element_text(size =18, color = "black"),
        legend.box.margin=margin(10,10,10,10), #move legend away from plot
        plot.margin = margin(25, 25, 25, 25)
        )+
  guides(fill = guide_legend(title = "Genotype")) #change legend title

bar
### combine dumbbell plot and barplot side-by-side
both = grid.arrange(plot2, bar, ncol = 2)
ggsave("selectClones_dumbbell_and_barplot_v4_102923.png", plot = both, width = 14, height = 9)
ggsave("selectClones_dumbbell_and_barplot_v4_102923.pdf", plot = both, width = 14, height = 9)
