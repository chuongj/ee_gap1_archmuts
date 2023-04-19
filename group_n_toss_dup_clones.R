##### Merge Cloneseq 0, 1, 2 CNV breakpoint + GAP1 Copy Number data frames #####
##### Remove duplicate clones (possibly not from unique lineages) ######
# Duplicate Clones (NON-Unique Clones) defined as if BOTH their start and ends points are within a 1000bp range of each other 
clonedata = read_csv("~Lab_docs/cloneseq1_temp/cloneseq0and1_Breaks_Copies.csv")
breaks = read_csv("~/Lab_docs/cloneseq2_temp/cloneseq2_rd_Breakpoints_Copies.csv")

all_clones = bind_rows(clonedata,breaks)
head(all_clones)
tail(all_clones)
unique(all_clones$Description)
unique(all_clones$pop_name)
length(unique(all_clones$pop_name)) #28 which is what it should be :) 
unique(all_clones$generation) #79 and 125

### THIS IS TOO HARD. I CANT APPEND EASILY. I NEED TO MOVE TO PYTHON SO I CAN APPEND.GROUPING ASSIGNMENTS #####

#for loop 
# for every row, get the start number 
# lefbounds = start number +/- 1kb
# rightbonds = end number +/- 1kb
# logical ask whether everyone in group has start in that range
# if yes group together
# if no exit. 
# 4-9-23 has some bugs - doesn't exactly do what i want it to do. work more on it later.
df = ltr2_125
get_groups = function(df){
#found_group = list(rep(NA, nrow(df))) #empty data frame
df$match_starts = rep(NA, nrow(df))
df$match_ends = rep(NA, nrow(df))
df$group_ID = 1:nrow(df) #set initial group ID


for (i in 1:nrow(df)){
  #starts
  start_lef = df$start[i]-1000
  start_right = df$start[i]+1000
  index_matching_starts = which(df$start < start_right & df$start > start_lef)
  df[ ,"match_starts"] <- df$start < start_right & df$start > start_lef

  if (index_matching_starts == 0) {
    next  #if none match startpoints, skip to next iteration of loop
  }
  #   else if (index_matching_starts == i) {
  #   next  #only in common with itself - next iteration 
  # }
  
  #ends
  end_lef = df$end[i]-1000
  end_right = df$end[i]+1000
  df[ ,"matching_ends"]=df$end < end_right & df$end > end_lef
  #index_matching_ends = which(df$end < end_right & df$end > end_lef)
  
  #group together indexes that appear in BOTH index_matching_starts and index_matching_ends
  results<-index_matching_starts %in% index_matching_ends
  group_index<-index_matching_starts[results]
 # group_ID < - paste(group_index) # need better way to assign / append new group ID after finding matches
  found_group[[i]] <- group_index
  
}
return(found_group[[i]])
}

find_groups = function(ltr2_125){
  ltr2_125$group_ID = 1:nrow(ltr2_125) #set initial group ID
  
  
  #compare each strain to the first strain
  row1<-ltr2_125[1,]
  start_lef = row1$start -1000
  start_rite = row1$start +1000
  end_lef = row1$end - 1000
  end_rite = row1$end + 1000
  for (i in nrow(df)){

  }
}
# hardcode from found_group list 
ltr1_79$grouping<-c("1", "2,3,4","2,3,4","2,3,4","5,6","5,6", "7", "8", "9")

### THIS IS TOO HARD. I CANT APPEND EASILY. I NEED TO MOVE TO PYTHON SO I CAN APPEND.GROUPING ASSIGNMENTS
### DO IT LATER

#### FOR NOW - JUST FIND DUPLICATES BY EYE  ####
# Duplicates defined as if both ends are within a +/- 1000 bp distance of each other

### LTR1 g79 ####
ltr1_79 = all_clones %>% filter(pop_name == "gap1_ltr_1", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 
# Toss  2686, 2689, 2687 

#ltr1_125
all_clones %>% filter(pop_name == "gap1_ltr_1", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 

# LTR2 g79 ##
ltr2_79<-all_clones %>% filter(pop_name=="gap1_ltr_2", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 
# Toss out 2850, 2851

# LTR2 g125
ltr2_125<-all_clones %>% filter(pop_name=="gap1_ltr_2", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 
# Toss out 2926 

#LTR3 g79
ltr3_79<-all_clones%>%filter(pop_name=="gap1_ltr_3", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# Toss none

#LTR3 g125
ltr3_125<-all_clones%>%filter(pop_name=="gap1_ltr_3", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss none

#"gap1_ltr_4"
all_clones%>%filter(pop_name=="gap1_ltr_4", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 
#toss none
all_clones%>%filter(pop_name=="gap1_ltr_4", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

#"gap1_ltr_5"
all_clones%>%filter(pop_name=="gap1_ltr_5", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 
#toss 2967, 2969
all_clones%>%filter(pop_name=="gap1_ltr_5", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) 
#Toss 3051

#"gap1_ltr_6"
all_clones%>%filter(pop_name=="gap1_ltr_6", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end) #2972 and 2973 are a group; 2971 and 2970 are a group
#Toss 2973, 2971 
all_clones%>%filter(pop_name=="gap1_ltr_6", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss none

#"gap1_ltr_7"
all_clones%>%filter(pop_name=="gap1_ltr_7", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss 2976
all_clones%>%filter(pop_name=="gap1_ltr_7", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

#"gap1_ltr_8"
all_clones%>%filter(pop_name=="gap1_ltr_8", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_ltr_8", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

##### Wildtype ##### 
#WT1 g79
all_clones%>%filter(pop_name=="gap1_1", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss none

# WT1 g125
all_clones%>%filter(pop_name=="gap1_1", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss 3024, 3025 

#WT2 g79
all_clones%>%filter(pop_name=="gap1_2", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss 2842
all_clones%>%filter(pop_name=="gap1_2", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss none

# gap1_3
all_clones%>%filter(pop_name=="gap1_3", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss out 2951, 2952
all_clones%>%filter(pop_name=="gap1_3", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss out 3030, 3115 

# "gap1_4" 
all_clones%>%filter(pop_name=="gap1_4", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss out 2954, 2956 
all_clones%>%filter(pop_name=="gap1_4", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss out 3034

# "gap1_5" 
all_clones%>%filter(pop_name=="gap1_5", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss out 2958
all_clones%>%filter(pop_name=="gap1_5", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss out 3037

##### ARS #### 
# "gap1_ars_1"
all_clones%>%filter(pop_name=="gap1_ars_1", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_ars_1", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss out 3065, 3066, 3067

"gap1_ars_3" 
all_clones%>%filter(pop_name=="gap1_ars_3", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# Toss out 2986
all_clones%>%filter(pop_name=="gap1_ars_3", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss out 3071

"gap1_ars_4"
all_clones%>%filter(pop_name=="gap1_ars_4", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_ars_4", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss out 3073, 3074, 3075

# "gap1_ars_5" 
all_clones%>%filter(pop_name=="gap1_ars_5", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss 2853 
all_clones%>%filter(pop_name=="gap1_ars_5", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss none

#"gap1_ars_6"
all_clones%>%filter(pop_name=="gap1_ars_6", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss out 2993, 2995
all_clones%>%filter(pop_name=="gap1_ars_6", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss out 3078, 3079 

#"gap1_ars_7"
all_clones%>%filter(pop_name=="gap1_ars_7", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#Toss out 2997
all_clones%>%filter(pop_name=="gap1_ars_7", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

#"gap1_ars_8" 
all_clones%>%filter(pop_name=="gap1_ars_8", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_ars_8", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss 3085

##### ALL Removed #####
# "gap1_all_1" 
all_clones%>%filter(pop_name=="gap1_all_1", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# Toss out 2680
all_clones%>%filter(pop_name=="gap1_all_1", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss out 3091 

# "gap1_all_2"
all_clones%>%filter(pop_name=="gap1_all_2", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# Toss none
all_clones%>%filter(pop_name=="gap1_all_2", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# Toss 2944 

# "gap1_all_3" 
all_clones%>%filter(pop_name=="gap1_all_3", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss none
all_clones%>%filter(pop_name=="gap1_all_3", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss 2939, 2940

"gap1_all_4"
all_clones%>%filter(pop_name=="gap1_all_4", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_all_4", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
# toss out 3095

"gap1_all_5"
all_clones%>%filter(pop_name=="gap1_all_5", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_all_5", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

"gap1_all_6"
all_clones%>%filter(pop_name=="gap1_all_6", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_all_6", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

#"gap1_all_7"
all_clones%>%filter(pop_name=="gap1_all_7", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_all_7", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

#"gap1_all_8"
all_clones%>%filter(pop_name=="gap1_all_8", generation == 79) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none
all_clones%>%filter(pop_name=="gap1_all_8", generation == 125) %>% arrange(cnv_length) %>% select(sample, generation, pop_name, cnv_length, start, end)
#toss none

######### Remove the duplicate clones ###### 
duplicated_clones<-c(2686,2687,2689,2850,2851,2926,
                     2853,2680,2944,2939,2940,3025,3024,2967,2969,3051,
                     2973,2971,2976,2842,2951,2952,3030,3115,2954,2956,3034,2958,3037,3065,3066,3067,2986,
                     3071,3073,3074,3075,2993,2995,3078,3079,2997,3085,3091,3095)

new = all_clones %>% dplyr::filter(!sample %in% duplicated_clones)

#add back the aneuploidy clones
aneu = read_csv("/Users/juliechuong/Lab_docs/cloneseq1_temp/cloneseq1_rd_Breakpoints_Copies.csv")
aneu = aneu %>% 
  filter(sample %in% c(2927, 2926, 2925, 2851,2850, 2849))
new = rbind(new, aneu)

######### Plot amplified region #######
library(ggalt)
library(gridExtra)
quartz()
  new %>%  
  filter(generation == 79) %>% 
#  filter(generation == 125) %>%     
#  filter(gap1_rounded>1, start !=1)%>% #filter out aneuploidies 
  filter(gap1_rounded>1) %>% #keep aneuploidy
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), 
#          desc(gap1_rounded),
          desc(cnv_length)) %>% #custom reorder by genotype, and then arrange copy number low to high 
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  mutate(pop_name = factor(pop_name)) %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
#  geom_dumbbell(size=3, dot_guide=TRUE, dot_guide_size=0.25)+  
  geom_dumbbell(size=2, dot_guide=TRUE, dot_guide_size=0.15)+
  scale_color_manual(#values = c(arsSalmons[1], allGolds[1],ltrBlues[1],wtGrays[3]), #custom colors
                     values = c(wtGrays[3], allGolds[1], ltrBlues[1], arsSalmons[1]),
                     limits = c("Wildtype architecture", "LTR and ARS KO","LTR KO","ARS KO"   ),
                     labels = c("Wildtype architecture", "LTR removed","ARS removed","LTR and ARS removed"))+ #third, now you can change legend labels
  scale_x_discrete(expand=c(0,1),
                  # limits=c(400000,670000),
                   labels = function(l) {trans = l / 1000},
                   "Position on Chromosome XI (kb)")+
  scale_x_continuous(expand=c(0,1),
                    # limits=c(400000,670000),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone ID") +
  # facet_grid(gap1_rounded ~ ., space = "fixed")+
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        strip.text.x = element_text(size = 12), #facet font size
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.text=element_text(size=12, color = "black"),
        legend.position = "none"
  )+
  guides(color = guide_legend(title = "Genotype"))
  
ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041023_noLegend_noAneu.png", width = 8, height = 10, bg = "white")
ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041023_noLegend_noAneu.pdf", width = 8, height = 10, bg = "white")   
ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041023_CnvLength_noLegend_noAneu.png", width = 8, height = 10, bg = "white")
ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041023_CnvLength_noLegend_noAneu.pdf", width = 8, height = 10, bg = "white")  
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041023_CnvLength_noLegend_noAneu.png", width = 8, height = 10, bg = "white")
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041023_CnvLength_noLegend_noAneu.pdf", width = 8, height = 10, bg = "white")  

ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041123_CnvLength_noLegend_noAneu.png", width = 8, height = 9, bg = "white")
ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041123_CnvLength_noLegend_noAneu.pdf", width = 8, height = 9, bg = "white")  
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041123_CnvLength_noLegend_noAneu.png", width = 8, height = 9, bg = "white")
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041123_CnvLength_noLegend_noAneu.pdf", width = 8, height = 9, bg = "white")  

ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041723_CnvLength_noLegend_noAneu.png", width = 12, height = 7, bg = "white")
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041723_CnvLength_noLegend_noAneu.png", width = 12, height = 7, bg = "white")

ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041823_CnvLength_noLegend_ANEU_WIDE.png", width = 20, height = 7, bg = "white")
ggsave(filename = "cloneseq0_1_2_g79_dumbbell_plot_041823_CnvLength_noLegend_noAeu.png", width = 12, height = 7, bg = "white")
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041823_CnvLength_noLegend_ANEU_WIDE.png", width = 20, height = 7, bg = "white")
ggsave(filename = "cloneseq0_1_2_g125_dumbbell_plot_041823_CnvLength_noLegend_noAeu.png", width = 12, height = 7, bg = "white")

##### Barplot for GAP1 Copy Number
### Barplot grouped by genotype, then sort by copy number low to high ####
#portion_g79%>%
new %>%  
  filter(generation == 125) %>% 
  filter(gap1_rounded>1)%>%
  #arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorders the strain order after arrange()
  ggplot(aes(clone, gap1_rounded, fill = Description)) +
  geom_bar(stat = "identity",  width = 0.7 )+
  scale_fill_manual(values = c(wtGrays[3], ltrBlues[1], arsSalmons[1], allGolds[1]), #custom colors
                    limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  coord_flip()+
  theme_classic() +
  ylab("GAP1 copy number") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.text=element_text(size=12, color = "black"),
        legend.box.margin=margin(10,10,10,10), #move legend away from plot
        plot.margin = margin(25, 25, 25, 25)
  )+
  guides(fill = guide_legend(title = "Genotype")) #change legend title

ggsave(filename = "cloneseq0_1_2_g125_barplot_041123_noAneu.png", width = 6, height = 9, bg = "white")
ggsave(filename = "cloneseq0_1_2_g125_barplot_041123_noAneu.pdf", width = 6, height = 9, bg = "white")   

#### histogram of CNV Length by Genotype
ggplot(all_clones, aes(log10(cnv_length), fill = Description)) + 
  geom_histogram(bins = 10,
                 aes(y = stat(width*density))) +
  #facet_wrap(~Description) +
  # scale_y_continuous(labels = scales::percent) +
  theme_bw()
