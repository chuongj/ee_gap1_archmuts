##### Figure 4 - Dumbbell plots and RD-estimated GAP1 copy numbers for sequenced clones ####
#### Combining data from cloneseq 0,1,2 repaired. 
### 1-31-24 Cleaned up the CNV mechanism categories data 
### Now do final Figure 4! 
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
library(scales)
library(RColorBrewer)

# set color pallete
ars_color = "#e26d5c"
arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631","#f85c46", "#bb3521","#d9402a")
wt_color = "gray50"
wtGrays = c("#354f52", "#666666", "#6b705c", "#414833" ,"#999999")
all_color = "#ffba08" 
all_color = "#D9BB59" 
allGolds=c("#ffba08", "#faa307", "#dda15e", "#7f5539", "#9c6644", "#fdc409", "#9c7e1e", "#D9BB59")
ltr_color = "#6699cc"
ltrBlues = c("#6699cc", "#005f73", "#0a9396", "#4292C6", "#2171B5", "#3799fb", "#66b3cc","#3a0ca3")

# brewer.pal(8, "Dark2")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
# "#D8BFD8" "#78184a" "#953553" "#800080"

my_palet = c("#E7298A","#66A61E", "#7570B3","#D95F02","#3799fb", "#1B9E77")

my.facet.labs = c("Wildtype", "LTR∆", "ARS∆", "ALL∆")
names(my.facet.labs) = c("Wildtype architecture", "LTR KO", "ARS KO", "LTR and ARS KO")

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
all_clones = read_csv("all_clones_RD_mech_020224.csv")
all_clones = all_clones %>% mutate_at('generation', as.factor)
all_clones$Description <- factor(all_clones$Description, levels=c('Wildtype architecture','LTR KO','ARS KO', 'LTR and ARS KO')) #reorder

#### calculate cnv lengths, breakpoint distances from GAP1 start codon ####
# 518204 is the start of GAP1 CDS
all_clones = all_clones %>% mutate(cnv_length = abs(start-end),
                      distance_start = start-518204,
                      distance_end = end-518204)
all_clones = all_clones %>%
  mutate(distance_start = distance_start*-1)%>%
  mutate(CNV_Mechanism4 = as_factor(if_else(grepl("ODIRA no ARS|ODIRA one end", CNV_Mechanism3), "ODIRA", CNV_Mechanism3))) %>% relocate(CNV_Mechanism4, .after = CNV_Mechanism3) 

#all_clones %>% write_csv("all_clones_RD_mech_020224.csv")

######## Test differences of CNV length between groups between timepoints #######
hist(all_clones$cnv_length)
shapiro.test(all_clones$cnv_length) #W = 0.59563, p-value < 2.2e-16 NOT NORMAL

##### Supplementary Figure 4_?_  ####
## Grouped boxplot of CNV Length by Genotypes & by Generation - Colorful ####
all_clones$Description <- factor(all_clones$Description, levels=c('Wildtype architecture','LTR KO','ARS KO', 'LTR and ARS KO')) #reorder

all_clones %>% group_by(Description) %>% summarize(median = median(cnv_length), iqr = IQR(cnv_length), min = min(cnv_length), max = max(cnv_length)) 

log10_plot = all_clones %>% ggplot(aes(x=Description, y=cnv_length, fill = interaction(Description,generation))) +
  geom_boxplot()+ 
  geom_point(size = 0.5, 
             position=position_jitterdodge())+
  # position=position_jitterdodge()
  # geom_jitter(aes(colour = c(Description,generation), x = Description), 
  #             position = position_jitter(width = .05), alpha = 0.5)+
  #scale_y_continuous(trans='log10')+
  scale_y_continuous(trans='log10',
                     breaks = c(3000, 10000, 30000, 100000, 300000, 670000),
                     labels = c(3000, 10000, 30000, 100000, 300000, 670000)/1000
                       )+
  #ylab(expression('log'[10]*'CNV Length (kb)')) +
  ylab("CNV length (kb)") +
  xlab("Genotype")+
  theme_classic()+
  scale_fill_manual(values=c("gray","#6699cc","#e26d5c","#ffba08",
                             "gray40","#0a9396","#da4631","#faa307"))+
  theme(
    #axis.text.x = element_blank(), #remove x-tick labels 
    #axis.ticks.x=element_blank(), #remove x-ticks 
    axis.text.x = element_text(size =12, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, vjust=2),
    text = element_text(size=16)
  )
log10_plot
ggsave("boxplot_CNVlength_Log_colorfulLabel_020124.png",
       width = 10, height = 4, bg = "white")
ggsave("boxplot_CNVlength_Log_colorfulLabel_021224.pdf",
       width = 10, height = 4, bg = "white")


####### Two way ANOVA #######
summary(aov(cnv_length~Description*generation, data = all_clones))

#                          Df    Sum Sq   Mean Sq F value  Pr(>F)   
# Description              3 1.587e+11 5.290e+10   4.023 0.00819 **
# generation               1 5.036e+07 5.036e+07   0.004 0.95071   
# Description:generation   3 2.556e+10 8.521e+09   0.648 0.58498   
# Residuals              221 2.906e+12 1.315e+10 

aov.org <- aov(
  cnv_length ~ Description * generation, data = all_clones,
  contrasts = list(
    Description = 'contr.sum',
    generation = 'contr.sum'
  )
)
summary(aov.org, type = 'III')


### Ln() and log10() transforming the data for ANOVA
shapiro.test(log(all_clones$cnv_length)) #not Normal
shapiro.test(log10(all_clones$cnv_length)) #not Normal
summary(aov(log(cnv_length)~Description*generation, data = all_clones))
summary(aov(log10(cnv_length)~Description*generation, data = all_clones))

aov.log <- aov(
  log10(cnv_length) ~ Description * generation, data = all_clones,
  contrasts = list(
    Description = 'contr.sum',
    generation = 'contr.sum'
  )
)
summary(aov.log, type = 'III')

###### Rank-Transformed Two way ANOVA #### 
aov.rnk = aov(rank(cnv_length)~Description*generation, data = all_clones, 
    contrasts = list(
      Description = 'contr.sum',
      generation = 'contr.sum'
    ))
summary(aov.rnk, type = 'III')

## Compare residuals from all 3 ## 
res.org = aov.org$resid
res.log = aov.log$resid
res.rnk = aov.rnk$resid
qqnorm(
  res.org, pch = 20, main = "Original Data",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(res.org)

qqnorm(
  res.log, pch = 20, main = "Log-Transformed",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(res.log)

qqnorm(
  res.rnk, pch = 20, main = "Rank-Transformed",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(res.rnk)

plot(aov.org, 1, main = "Original Data")
plot(aov.log, 1, main = "Log-Transformed")
plot(aov.rnk, 1, main = "Rank-Transformed")

###### Violin Boxplot CNV Length (Ungrouped) By Description only ########
ggplot(all_clones, aes(Description, cnv_length, fill = Description))+
  #geom_boxplot()+
  geom_violin(draw_quantiles = 0.5)+
 # geom_boxplot(width=.3)+
  geom_point(size = 0.4, 
             position=position_jitterdodge(1))+
  # scale_y_continuous(breaks=seq(0,669000,200000),
  #                    labels=seq(0,669000,200000)/1000)+
  scale_y_continuous(trans='log10',
                     breaks = c(3000, 10000, 30000, 100000, 300000, 670000),
                     labels = c(3000, 10000, 30000, 100000, 300000, 670000)/1000
  )+
  ylab("CNV length (kb)") +
  xlab("Genotype")+
  scale_fill_manual(values=c('gray40', ltr_color, ars_color, all_color))+
  theme_classic()+
  theme(
    legend.position="none",
    axis.text.x = element_blank(), #remove x-tick labels 
    axis.ticks.x=element_blank(), #remove x-ticks 
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 20, vjust=2),
    text = element_text(size=20)
  )
# ggsave("boxplot_jit_CNVlength_Log_UNgroup_011324.png",
#        width = 6, height = 5, bg = "white")
# ggsave("boxplot_jit_CNVlength_Log_UNgroup_011324.pdf",
#        width = 6, height = 5, bg = "white")
# ggsave("Violin_jit_CNVlength_Log_UNgroup_020124.png",
#        width = 6, height = 5, bg = "white")
# ggsave("Violin_jit_CNVlength_Log_UNgroup_020124.pdf",
#        width = 6, height = 5, bg = "white")
# ggsave("Violin_jit_CNVlength_NoLegend_020224.png",
#        width = 6, height = 5, bg = "white")
# ggsave("Violin_jit_CNVlength_NoLegend_020224.pdf",
#        width = 6, height = 5, bg = "white")

####### One Way ANOVA - CNV length ~ Description #####
summary(aov(cnv_length~Description, data = all_clones))
#               Df    Sum Sq   Mean Sq F value  Pr(>F)   
# Description   3 1.587e+11 5.290e+10    4.06 0.00778 **
#   Residuals   225 2.932e+12 1.303e+10 

####### Kruskal-Wallis Test - CNV length ~ Description ######
kruskal.test(cnv_length~Description, data = all_clones)
# Kruskal-Wallis chi-squared = 18.824, df = 3, p-value = 0.0002973

####### Pairwise Wilcox Mann Whitney Test with Bonferonni Correction ####### 
pairwise.wilcox.test(all_clones$cnv_length, all_clones$Description, p.adjust.method = "bonferroni")
#               Wildtype architecture LTR KO ARS KO
# LTR KO         0.0061                -      -     
# ARS KO         1.0000                0.0089 -     
# LTR and ARS KO 0.0416                1.0000 0.0167       

# Nonparametric Two-way ANOVA
# Rank Transformation  https://www.cfholbert.com/blog/nonparametric_two_way_anova/

# Poisson assumption is that mean = variance
# If that's not true, use negative binomial 
# Try both tests anyway for fun
all_clones %>% group_by(Description) %>% summarize(mean = mean(cnv_length), var = var(cnv_length))

#poisson.test
# glm.nb

##### Scatterplot of CNV start vs CNV ends  ########
all_clones %>%
  #filter(!start == 1) %>% #exclude aneuploids
  ggplot(aes(x = distance_start, y=distance_end, color=Description)) +
  geom_point(size = 2)+
  scale_x_continuous(#trans = scales::pseudo_log_trans(1),
                     trans = 'log10',
                     breaks = c(0,3000,10000,30000,70000,518000),
                     labels = scales::comma_format(scale = 0.001))+
  scale_y_continuous(trans='log10',
                     breaks = c(0,1e3, 1e4, 1e5, 153000),
                     labels = scales::comma_format(scale = 0.001)) +
  scale_color_manual(values = c(wt_color, ltr_color, ars_color, all_color))+
  geom_hline(yintercept=1090,linetype=3)+  #dotted line for end of GAP1 CDS
  geom_vline(xintercept=0,linetype=3)+ #dotted line for Start Codon of GAP1 CDS
  theme_classic()+
  xlab("Upstream breakpoint (kb)") +
  ylab("Downstream breakpoint (kb)") +
  theme(axis.text.y = element_text(size =14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        #axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0.5),
        text = element_text(size=16)
  )
# ggsave(filename = "scatter_cnvStartEnd_020124.png",
#        width = 7,
#        height = 5,
#        bg = "white")
# ggsave(filename = "scatter_cnvStartEnd_020124.pdf",
#        width = 7,
#        height = 5,
#        bg = "white")


##### Scatterplot by Mechanism and Copy number, Facet by Description #####
db_clones = all_clones %>%
  filter(gap1_rounded>1) %>%
  filter(sample %in% c(2927,2929,3023, 3029, 3025, 2695, 3027, 2926,2931,3050,
                       2924,2690,3056,2948,2945,3069,2691,2852,2983,2946,3076,
                       3019,3022,3107,3099,3101,3093,2968)
  )  
all_clones %>%   
  filter(CNV_Mechanism4 %in% c("LTR NAHR", "ODIRA", "transposon-mediated", "aneuploid", "NAHR", "complex CNV"))%>% 
  mutate(CNV_Mechanism4 = fct_relevel(CNV_Mechanism4, "LTR NAHR", "NAHR", "ODIRA", "transposon-mediated", "aneuploid", "complex CNV")) %>% #custom order
  ggplot(aes(x = distance_start, y=distance_end, color=CNV_Mechanism4, shape = as.factor(gap1_rounded))) +
  geom_point(size = 4, alpha = 1)+
  scale_x_continuous(#trans = scales::pseudo_log_trans(1),
    trans = 'log10',
    breaks = c(0,3000,10000,30000,70000,518000),
    labels = scales::comma_format(scale = 0.001))+
  scale_y_continuous(trans='log10',
                     breaks = c(0,1e3, 1e4, 1e5, 153000),
                     labels = scales::comma_format(scale = 0.001)) +
  #scale_color_manual(values = c(wt_color, ltr_color, ars_color, all_color))+
  geom_text( data = db_clones,
    aes(label = sample),
    size = 3,
    nudge_x = -.25, nudge_y = 0,
    check_overlap = T )+
  # geom_label(
  #   data = all_clones %>% filter(sample %in% c(2927,2929,3023, 3029, 3025, 2695, 3027, 2926,2931,3050,
  #                                 2924,2690,3056,2948,2945,3069,2691,2852,2983,2946,3076,
  #                                 3019,3022,3107,3099,3101,3093,2968)
  #   ),
  #   aes(label=sample)
  # ) +
  scale_color_manual(values= rev(my_palet))+
  scale_shape_manual(values=c(19, 17, 15, 8))+
  geom_hline(yintercept=1090,linetype=3)+  #dotted line for end of GAP1 CDS
  geom_vline(xintercept=0,linetype=3)+ #dotted line for Start Codon of GAP1 CDS
  theme_minimal()+
  xlab("Upstream breakpoint (kb)") +
  ylab("Downstream breakpoint (kb)")+
  labs(shape="GAP1 Copy Number",color="CNV Mechanism")+#edit legend title
  theme(axis.text.y = element_text(size =14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        panel.spacing = unit(1.5, "lines"),
        text = element_text(size=16)
  )+
  facet_wrap(~Description, 
             labeller = labeller(Description = my.facet.labs))

# ggsave("scatterplot_Mechs_CopyNum_FacetDescription_020124.png", width = 10, height = 8, bg = "white")
# ggsave("scatterplot_Mechs_CopyNum_FacetDescription_020124.pdf", width = 10, height = 8, bg = "white")
# ggsave("scatterplot_Mechs_CopyNum_FacetDescription_Label_020524.png", width = 12, height = 8, bg = "white")
# ggsave("scatterplot_Mechs_CopyNum_FacetDescription_Label_020524.pdf", width = 12, height = 8, bg = "white")

f=2

#SUPPLEMENTARY FIGURE#### Barplot of Mechanisms by Generation ####### 
all_clones %>%
  mutate(distance_start = distance_start*-1)%>%
  mutate(CNV_Mechanism4 = as_factor(if_else(grepl("ODIRA no ARS|ODIRA one end", CNV_Mechanism3), "ODIRA", CNV_Mechanism3))) %>% relocate(CNV_Mechanism4, .after = CNV_Mechanism3) %>% 
  filter(CNV_Mechanism4 %in% c("LTR NAHR", "ODIRA", "transposon-mediated", "aneuploid", "NAHR", "complex CNV"))%>% 
#  mutate(CNV_Mechanism4 = fct_relevel(CNV_Mechanism4, "LTR NAHR", "NAHR", "ODIRA", "transposon-mediated", "aneuploid", "complex CNV")) %>% #custom order
ggplot(aes(x=generation, fill = CNV_Mechanism4)) +
  geom_bar(position="dodge")+
  # geom_bar()+
  # scale_fill_manual(values=c('gray40', ltr_color, ars_color, all_color))+
  xlab("Generation")+
  ylab("Count")+
  theme_classic()+
  scale_fill_manual(values= rev(my_palet))+
  theme(
    axis.text.x = element_text(size = 14, color = "black",angle = 60, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, vjust=2),
    text = element_text(size=16)
  )

## Horizontal Barplot - Condensed Bins for Talk #####
all_clones %>%
  count(CNV_Mechanism4, Description) %>%
  group_by(Description) %>%   # now required with changes to dplyr::count()
  mutate(prop = prop.table(n)) %>% 
  ungroup() %>% 
  complete(CNV_Mechanism4, Description) %>% 
  mutate(across('CNV_Mechanism4', str_replace, 'transposon-mediated', 'Transposon mediated')) %>% 
  mutate(across('CNV_Mechanism4', str_replace, 'LTR NAHR', 'LTR Recombination')) %>% 
  mutate(across('CNV_Mechanism4', str_replace, 'aneuploid|NAHR|complex CNV', 'Other')) %>%  
  filter(CNV_Mechanism4 %in% c("LTR Recombination", "ODIRA", "Transposon mediated", "Other"))%>%
  mutate(CNV_Mechanism4 = fct_relevel(CNV_Mechanism4, "Other", "Transposon mediated", "LTR Recombination", "ODIRA" )) %>% #custom order
  mutate(Description = fct_relevel(Description, "LTR and ARS KO", "ARS KO", "LTR KO", "Wildtype architecture")) %>%   
  ggplot(aes(y=CNV_Mechanism4, x = n, fill = Description)) +
  geom_col(width = 0.75, position = position_dodge())+
  scale_fill_manual(values=rev(c('gray40', ltr_color, ars_color, all_color)),
                    guide = guide_legend(reverse = TRUE))+
  scale_x_continuous(limits = c(0, 50),
                     expand = c(0, 0),
                     position = "top"
  )+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 5))+
  ylab("Inferred CNV mechanism")+
  xlab("Count")+
  theme_classic()+
  geom_hline(yintercept=c(0.5,1.5,2.5,3.5, 4.5, 5.5, 6.5),color="#A8BAC4",size = 0.1)+
  theme(
    panel.grid.major.x = element_line(color = "gray70", size = 0.3), #vertical grid lines
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20),
    text = element_text(size=20),
    legend.title=element_blank()
  )

# ggsave("HorizontalBar_cnvMechs_043024.png", width = 10, height = 6, bg = "white")
# ggsave("HorizontalBar_cnvMechs_043024.pdf", width = 10, height = 6, bg = "white")

#### Horizontal Barplot - CNV Mechanisms per Genotype #### 
# all_clones %>%
#   # mutate(CNV_Mechanism4 = as_factor(if_else(grepl("ODIRA no ARS|ODIRA one end", CNV_Mechanism3), "ODIRA", CNV_Mechanism3))) %>% relocate(CNV_Mechanism4, .after = CNV_Mechanism3) %>%
#   # group_by(Description) %>% # now required with changes to dplyr::count()
#   mutate(prop = prop.table()) %>% View()
#   ungroup() %>%
#   complete(CNV_Mechanism4, Description) %>% View()
#   mutate(across('CNV_Mechanism4', str_replace, 'transposon-mediated', 'Transposon mediated')) %>% 
#   filter(CNV_Mechanism4 %in% c("LTR NAHR", "ODIRA", "Transposon mediated", "aneuploid", "NAHR", "complex CNV"))%>%  
#   mutate(across(CNV_Mechanism4, str_replace, 'LTR NAHR', 'LTR Recombination')) %>%
#   mutate(across(CNV_Mechanism4, str_replace, 'aneuploid|NAHR|complex CNV', 'Other')) %>% 
#   count(CNV_Mechanism4, Description) %>% 
#   # mutate(CNV_Mechanism4 = fct_relevel(CNV_Mechanism4, "NAHR", "aneuploid","transposon mediated", "complex CNV", "LTR NAHR", "ODIRA" )) %>% #custom order
#   # mutate(Description = fct_relevel(Description, "LTR and ARS KO", "ARS KO", "LTR KO", "Wildtype architecture")) %>% 
#     mutate(CNV_Mechanism4 = fct_relevel(CNV_Mechanism4, "Other", "Transposon mediated", "LTR Recombination", "ODIRA" )) %>% #custom order
#     mutate(Description = fct_relevel(Description, "LTR and ARS KO", "ARS KO", "LTR KO", "Wildtype architecture")) %>%   
  ggplot(aes(y=CNV_Mechanism4, x = n, fill = Description)) +
  geom_col(width = 0.75, position = position_dodge())+
  scale_fill_manual(values=rev(c('gray40', ltr_color, ars_color, all_color)),
                    guide = guide_legend(reverse = TRUE))+
  scale_x_continuous(limits = c(0, 50),
                     expand = c(0, 0),
                     position = "top"
                     )+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 5))+
  ylab("Inferred CNV mechanism")+
  xlab("Count")+
  theme_classic()+
  geom_hline(yintercept=c(0.5,1.5,2.5,3.5, 4.5, 5.5, 6.5),color="#A8BAC4",size = 0.1)+
  theme(
    panel.grid.major.x = element_line(color = "gray70", size = 0.3), #vertical grid lines
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20),
    text = element_text(size=20),
    legend.title=element_blank()
  )
# ggsave("HorizontalBar_cnvMechs_Genotype_020124.png", width = 8, height = 6, bg = "white")
# ggsave("HorizontalBar_cnvMechs_Genotype_020124.pdf", width = 8, height = 6, bg = "white")
# ggsave("HorizontalBar_cnvMechs_Genotype_020524.png", width = 10, height = 6, bg = "white")
# ggsave("HorizontalBar_cnvMechs_Genotype_020524.pdf", width = 10, height = 6, bg = "white")

#### chi-sq test for association of genotype with CNV mechanism ####
# Null: Genotypes and CNv mechanisms are independent, no association
# Alt Hypothesis: CNV mechanism counts observed are dependent on genotype

# make contingency table
conTbl = as.data.frame(with(all_clones, table(Description, CNV_Mechanism4))) %>% 
  pivot_wider(id_cols = CNV_Mechanism4, names_from = Description, values_from = Freq) %>% 
  filter(!(CNV_Mechanism4 %in% c("skip cuz clonal","unresolved","seq failed"))) %>%
  as.data.frame()
rownames(conTbl) = conTbl$CNV_Mechanism4
conTbl = select(conTbl, 2:5)

cs_test = chisq.test(conTbl)
cs_test 
# Warning message:
# In chisq.test(conTbl) : Chi-squared approximation may be incorrect
# Pearson's Chi-squared test
# data:  conTbl
# X-squared = 158.09, df = 15, p-value < 2.2e-16
cs_test$observed
round(cs_test$expected,1)
round(cs_test$residuals,2)

library(corrplot)
corrplot(cs_test$residuals, is.corr = F)

contrib = 100*cs_test$residuals^2/cs_test$statistic
round(contrib,2)
corrplot(contrib, is.corr = F)

chisq.test(conTbl, simulate.p.value = TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  conTbl
# X-squared = 158.09, df = NA, p-value = 0.0004998

#### Fisher's Exact Test ####

fisher.test(conTbl, simulate.p.value = T)
#Fisher's Exact Test for Count Data with simulated p-value (based on 2000 replicates)
# data:  conTbl
# p-value = 0.0004998
# alternative hypothesis: two.sided

library(rstatix)

conTbl[c(3,5), 1:2]
pairwise_fisher_test(conTbl[c(3,5), 1:2], p.adjust.method = "bonferroni", detailed = T)
#### poisson test ####
# First, check that my count data follows poisson distribution? Goodness of fit test



f=1


# SUPPLEMENTARY BARPLOT, Types of ODIRA we found ####
all_clones %>% select(CNV_Mechanism3, gap1_rounded) %>% filter(CNV_Mechanism3 %in% c("ODIRA one end", "ODIRA", "ODIRA no ARS")) %>% arrange(gap1_rounded) %>% summarize(ODIRA_2 = sum(CNV_Mechanism3 == "ODIRA" & gap1_rounded == 2),
            ODIRA_3 = sum(CNV_Mechanism3 == "ODIRA" & gap1_rounded == 3),
            ODIRA_oneEnd_2 = sum(CNV_Mechanism3 == "ODIRA one end" & gap1_rounded == 2),
            ODIRA_oneEnd_3 = sum(CNV_Mechanism3 == "ODIRA one end" & gap1_rounded == 3),
            ODIRA_noARS_3 = sum(CNV_Mechanism3 == "ODIRA no ARS" & gap1_rounded==3)
            ) %>% pivot_longer(cols = ODIRA_2:ODIRA_noARS_3, names_to = "ODIRA_Type", values_to = "Count") %>%
  mutate(Proportion = Count/sum(Count)) %>% 
  #ggplot(aes(x = reorder(ODIRA_Type, -Proportion), y=Proportion)) +
  ggplot(aes(x = reorder(ODIRA_Type, -Count), y=Count)) +
  ylim(0,60)+
  geom_col(fill = "#414833")+
  xlab("ODIRA Type")+
  ylab("Proportion")+
  theme_classic()

ggsave(filename = "barplot_ODIRAtypes_Count.png", width = 6, height = 4, bg = "white")
  

# SUPPLEMENTARY TABLE, Count of Mechanism for each Genotype ####
library(janitor)
with(all_clones,table(Description, CNV_Mechanism4)) %>% as.data.frame() %>%
#  filter(!(CNV_Mechanism4 %in% c("seq failed", "skip cuz clonal", "unresolved"))) %>%
  pivot_wider(names_from = Description, CNV_Mechanism4, values_from = Freq) %>% 
  mutate(total = rowSums(.[, 2:5])) %>%
  adorn_totals("row")%>% View() #write_csv("count_Mechs_perGeno.csv")

with(all_clones,table(CNV_Mechanism4, Description)) %>% as.data.frame() %>% 
  filter(!(CNV_Mechanism4 %in% c("seq failed", "skip cuz clonal", "unresolved"))) %>%
  group_by(Description) %>% 
  mutate(proportion_geno = Freq/sum(Freq)) %>% ungroup() %>% 
  mutate(across('CNV_Mechanism4', str_replace, 'transposon-mediated', 'transposon mediated')) %>% 
  mutate(CNV_Mechanism4 = fct_relevel(CNV_Mechanism4, "NAHR", "aneuploid","transposon mediated", "complex CNV", "LTR NAHR", "ODIRA" )) %>% #custom order
  mutate(Description = fct_relevel(Description, "LTR and ARS KO", "ARS KO", "LTR KO", "Wildtype architecture")) %>% 
  ggplot(aes(y=CNV_Mechanism4, x = proportion_geno, fill = Description)) +
  geom_col(width = 0.75, position = position_dodge())+
  scale_fill_manual(values=rev(c('gray40', ltr_color, ars_color, all_color)),
                    guide = guide_legend(reverse = TRUE))+
  scale_x_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     position = "top"
  )+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 5))+
  ylab("Inferred CNV mechanism")+
  xlab("Proportion")+
  theme_classic()+
  geom_hline(yintercept=c(0.5,1.5,2.5,3.5, 4.5, 5.5, 6.5),color="#A8BAC4",size = 0.1)+
  theme(
    panel.grid.major.x = element_line(color = "gray70", size = 0.3), #vertical grid lines
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20),
    text = element_text(size=20),
    legend.title=element_blank()
  )

ggsave("HorizontalBarProp_cnvMechs_Genotype_020624.png", width = 10, height = 6, bg = "white")
ggsave("HorizontalBarProp_cnvMechs_Genotype_020624.pdf", width = 10, height = 6, bg = "white")
f=1
#### Dumbbell Plot - Grouped by Genotype  ####
all_clones %>%
  filter(gap1_rounded>1,
         !(CNV_Mechanism4 %in% c("seq failed", "skip cuz clonal", "unresolved") 
         ))%>% 
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=3, dot_guide=F, dot_guide_size=0.25)+
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

#ggsave(filename = "dumbbell_plot_020124.png", width = 8, height = 20, bg = "white")

### Dumbbell plot of select clones for Fig 4  ####
all_clones %>%
  filter(gap1_rounded>1) %>%
  filter(sample %in% c(2927,2929,3023, 3029, 3025, 2695, 3027, 2926,2931,3050,
                       2924,2690,3056,2948,2945,3069,2691,2852,2983,2946,3076,
                       3019,3022,3107,3099,3101,3093,2968)
         ) %>% 
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))),desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=2, dot_guide=F, dot_guide_size=0.1)+
  coord_cartesian(xlim=c(360000,670191))+
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
        axis.text.x = element_text(size = 25, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        #axis.text.y= element_blank(),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 25, color = "black"),
        legend.text=element_text(size=25, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))

ggsave(filename = "dumbbell_Fig4_020224.png", width = 9, height = 4, bg = "white")
ggsave(filename = "dumbbell_Fig4_020224.pdf", width = 9, height = 4, bg = "white")

### Make dumbbell plot VERSION 2 ####  
# Include aneuploidy clones and zoom into 400-660kb region
# Try this to zoom into plot without removing data 
# https://www.geeksforgeeks.org/zoom-into-ggplot2-plot-without-removing-data-in-r/?ref=ml_lbp 
plot2 = 
  all_clones %>%
  filter(gap1_rounded>1) %>%
  filter(!(CNV_Mechanism4 %in% c("skip cuz clonal","seq failed","unresolved"))) %>%  #Only plot clones with RESOLVED mechs
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

plot2

# ggsave(filename = "dumbbell_plot_ZOOM_020124_long.png", width = 8, height = 20, bg = "white")
# ggsave(filename = "dumbbell_plot_ZOOM_020124_long.pdf", plot = plot2, width = 8, height = 20, bg = "white") 
# ggsave(filename = "dumbbell_plot_ZOOM_Resolved_020124.png", width = 8, height = 20, bg = "white")
# ggsave(filename = "dumbbell_plot_ZOOM_Resolved_020124.pdf", plot = plot2, width = 8, height = 20, bg = "white")

### Barplot grouped by genotype, then sort by copy number low to high ####
bar = all_clones %>%
  filter(gap1_rounded>1)%>%
  filter(!(CNV_Mechanism4 %in% c("skip cuz clonal","seq failed","unresolved"))) %>%  #Only plot clones with RESOLVED mechs
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
both
# ggsave("selectClones_dumbbell_and_barplot_v4_102923.png", plot = both, width = 14, height = 9)
# ggsave("selectClones_dumbbell_and_barplot_v4_102923.pdf", plot = both, width = 14, height = 9)
ggsave("dumbbell_and_barplot_Resolved_020224.png", plot = both, width = 16, height = 20)
ggsave("dumbbell_and_barplot_Resolved_020224.pdf", plot = both, width = 16, height = 20)

######## Split it up by generation ####### 

#### Generation 79 
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
# ggsave("g79_dumbbell_and_barplot_112923.pdf", plot = g79_both, width = 14, height = 12)
# ggsave("g79_dumbbell_and_barplot_112923.png", plot = g79_both, width = 14, height = 12)
