##### Figure 4 - Dumbbell plots and RD-estimated GAP1 copy numbers for sequenced clones ####
#### Combining data from cloneseq 0,1,2 repaired. 

## mCitrine is located 514934 - 515650
## coordinates of GAP1 CDS for our samples are 518204 - 520012
## LTR11 coordinates are 513532:513846, LTR12 is ~ 520925

#Set working directory and samplesheet
setwd("/Users/juliechuong/Lab_docs/cloneseq2_newBarcodes")

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
clones0_1 = read_csv("../cloneseq1_temp/cloneseq0and1_Breaks_Copies.csv")
clones2 = read_csv("cloneseq2_RD_results_repaired_092523.csv")
all_clones = rbind(clones0_1, clones2)

#### Dumbbell Plot - Grouped by Genotype  ####
db_geno=
  all_clones %>%
  filter(gap1_rounded>1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=3, dot_guide=TRUE, dot_guide_size=0.25)+
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
        #axis.text.y= element_blank(),
        axis.title = element_text(size = 30, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text=element_text(size=30, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
  #facet_grid(~gap1_rounded)
db_geno

# ggsave("all_clones0and1_103123.png", width = 8, height = 10, bg = "white")
# ggsave("all_clones0and1_103123.pdf", width = 8, height = 10, bg = "white")
 # ggsave(filename = "selectClones_dumbbell_plot_102923.png", width = 8, height = 10, bg = "white")
# ggsave(filename = "selectClones_dumbbell_plot_102923.pdf", width = 8, height = 10, bg = "white") 
ggsave(filename = "dumbbell_plot_112923.png", width = 8, height = 10, bg = "white")


# Make dumbbell plot VERSION 2 
# Include aneuploidy clones and zoom into 400-660kb region
# Try this to zoom into plot without removing data 
# https://www.geeksforgeeks.org/zoom-into-ggplot2-plot-without-removing-data-in-r/?ref=ml_lbp 
plot2 = 
  all_clones %>%
  filter(gap1_rounded>1) %>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=1, dot_guide=F, dot_guide_size=0.25)+
  coord_cartesian(xlim=c(400000,670191))+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670191),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        #axis.text.y= element_blank(),
        axis.title = element_text(size = 26, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text=element_text(size=30, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))

plot2

# ggsave(filename = "selectClones_ZOOM2_dumbbell_plot_102923.png", width = 8, height = 10, bg = "white")
# ggsave(filename = "selectClones_ZOOM2_dumbbell_plot_102923.pdf", width = 8, height = 10, bg = "white") 
ggsave(filename = "dumbbell_plot_ZOOM_112923.png", width = 8, height = 10, bg = "white")
ggsave(filename = "dumbbell_plot_ZOOM_112923.pdf", plot = plot2, width = 8, height = 10, bg = "white") 
ggsave(filename = "dumbbell_plot_ZOOM_112923_long.png", width = 8, height = 20, bg = "white")
ggsave(filename = "dumbbell_plot_ZOOM_112923_long.pdf", plot = plot2, width = 8, height = 20, bg = "white") 
### Barplot grouped by genotype, then sort by copy number low to high ####
bar = all_clones %>%
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
# ggsave("selectClones_dumbbell_and_barplot_v4_102923.png", plot = both, width = 14, height = 9)
# ggsave("selectClones_dumbbell_and_barplot_v4_102923.pdf", plot = both, width = 14, height = 9)
ggsave("dumbbell_and_barplot_112923.png", plot = both, width = 16, height = 20)
ggsave("dumbbell_and_barplot_112923.pdf", plot = both, width = 16, height = 20)

######## Split it up by generation ####### 

#### Generation 79 ####
g79_plot = all_clones %>%
  filter(gap1_rounded>1) %>%
  filter(generation == 79) %>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=2, dot_guide=F, dot_guide_size=0.25)+
  coord_cartesian(xlim=c(400000,670191))+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670191),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        #axis.text.y= element_blank(),
        axis.title = element_text(size = 26, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text=element_text(size=30, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
g79_plot

ggsave(filename = "dumbbell_g79_ZOOM_112923.png", plot = g79_plot, width = 8, height = 10, bg = "white")
ggsave(filename = "dumbbell_g79_ZOOM_112923.pdf", plot = g79_plot, width = 8, height = 10, bg = "white") 
ggsave(filename = "dumbbell_g79_ZOOM_112923_long.png", plot = g79_plot, width = 8, height = 20, bg = "white")
ggsave(filename = "dumbbell_g79_ZOOM_112923_long.pdf", plot = g79_plot, width = 8, height = 20, bg = "white") 


g125_plot = all_clones %>%
  filter(gap1_rounded>1) %>%
  filter(generation == 125) %>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=2, dot_guide=F, dot_guide_size=0.25)+
  coord_cartesian(xlim=c(400000,670191))+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670191),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        #axis.text.y= element_blank(),
        axis.title = element_text(size = 26, color = "black"),
        legend.title = element_text(size = 30, color = "black"),
        legend.text=element_text(size=30, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
g125_plot

ggsave(filename = "dumbbell_g125_ZOOM_112923.png", plot = g125_plot, width = 8, height = 10, bg = "white")
ggsave(filename = "dumbbell_g125_ZOOM_112923.pdf", plot = g125_plot, width = 8, height = 10, bg = "white") 
ggsave(filename = "dumbbell_g125_ZOOM_112923_long.png", plot = g125_plot, width = 8, height = 20, bg = "white")
ggsave(filename = "dumbbell_g125_ZOOM_112923_long.pdf", plot = g125_plot, width = 8, height = 20, bg = "white") 

g125bar = all_clones %>%
  filter(gap1_rounded>1)%>%
  filter(generation == 125) %>%
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
g125bar

g125_both = grid.arrange(g125_plot, g125bar, ncol = 2)
ggsave("g125_dumbbell_and_barplot_112923.pdf", plot = g125_both, width = 14, height = 12)
ggsave("g125_dumbbell_and_barplot_112923.png", dpi = 300, plot = g125_both, width = 14, height = 12)

g79bar = all_clones %>%
  filter(gap1_rounded>1)%>%
  filter(generation == 79) %>%
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
g79bar

g79_both = grid.arrange(g79_plot, g79bar, ncol = 2)
ggsave("g79_dumbbell_and_barplot_112923.pdf", plot = g79_both, width = 14, height = 12)
ggsave("g79_dumbbell_and_barplot_112923.png", plot = g79_both, width = 14, height = 12)
