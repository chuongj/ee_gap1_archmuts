# Chuong et al. 2024
# Figure 4 

setwd("")

#load libraries
library(tidyverse)
library(ggalt)
library(scales)

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
my_palet = c("#E7298A","#66A61E", "#7570B3","#D95F02","#3799fb", "#1B9E77")
my.facet.labs = c("Wildtype", "LTR∆", "ARS∆", "ALL∆")
names(my.facet.labs) = c("Wildtype architecture", "LTR KO", "ARS KO", "LTR and ARS KO")

# import clone data
all_clones = read_csv("clone_sequencing_analysis.csv")
all_clones = all_clones %>% mutate_at('generation', as.factor)
all_clones$Description <- factor(all_clones$Description, levels=c('Wildtype architecture','LTR KO','ARS KO', 'LTR and ARS KO')) #reorder

#### Figure 4B ####
# Violin Boxplot CNV Length By Strain #
ggplot(all_clones, aes(Description, cnv_length, fill = Description))+
  geom_violin(draw_quantiles = 0.5)+
  geom_point(size = 0.4, 
             position=position_jitterdodge(1))+
  # scale_y_continuous(breaks=seq(0,669000,200000),
  #                    labels=seq(0,669000,200000)/1000)+
  scale_y_continuous(trans='log10',
                     breaks = c(3000, 10000, 30000, 100000, 300000, 670000),
                     labels = c(3000, 10000, 30000, 100000, 300000, 670000)/1000
  )+
  ylab("CNV length (kb)") +
  xlab("Strain")+
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

## Kruskal-Wallis Test - CNV length ~ Description ###
kruskal.test(cnv_length~Description, data = all_clones)

## Pairwise Wilcox Mann Whitney Test with Bonferonni Correction ## 
pairwise.wilcox.test(all_clones$cnv_length, all_clones$Description, p.adjust.method = "bonferroni")


#### Figure 4C ####
# Horizontal Barplot - CNV Mechanisms per Genotype #  
all_clones %>%
  count(CNV_Mechanism, Description) %>%
  group_by(Description) %>%   # now required with changes to dplyr::count()
  mutate(prop = prop.table(n)) %>% 
  ungroup() %>% 
  complete(CNV_Mechanism, Description) %>% 
  mutate(across('CNV_Mechanism', str_replace, 'transposon-mediated', 'transposon mediated')) %>% 
  filter(CNV_Mechanism %in% c("LTR NAHR", "ODIRA", "transposon mediated", "aneuploid", "NAHR", "complex CNV"))%>% 
  mutate(CNV_Mechanism = fct_relevel(CNV_Mechanism, "NAHR", "aneuploid","transposon mediated", "complex CNV", "LTR NAHR", "ODIRA" )) %>% #custom order
  mutate(Description = fct_relevel(Description, "LTR and ARS KO", "ARS KO", "LTR KO", "Wildtype architecture")) %>% 
  ggplot(aes(y=CNV_Mechanism, x = n, fill = Description)) +
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


## Fisher test for association of genotype with CNV mechanism #
# Null: Genotypes and CNV mechanisms are independent, no association
# Alt Hypothesis: CNV mechanism counts observed are dependent on strain

# make contingency table
conTbl = as.data.frame(with(all_clones, table(Description, CNV_Mechanism))) %>% 
  pivot_wider(id_cols = CNV_Mechanism, names_from = Description, values_from = Freq) %>% 
  as.data.frame()
rownames(conTbl) = conTbl$CNV_Mechanism
conTbl = select(conTbl, 2:5)

chisq.test(conTbl, simulate.p.value = TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# data:  conTbl
# X-squared = 158.09, df = NA, p-value = 0.0004998

### Fisher's Exact Test ##
fisher.test(conTbl, simulate.p.value = T)

#We observe an increase of LTR NAHR in the ARS∆ clones (27/52, 52%) 
# relative to wildtype clones (11/37, 39%) 
pairwise.prop.test(x = c(11, 27), n = c(37, 52), alternative = "greater", p.adjust.method = "holm")
prop.test(x = c(11, 27), n = c(37, 52), alternative = "less") # p-value = 0.03083

#chisq_test(conTbl[3,c(4,1)])

#Is there is significant increase in ODIRA in from WT to LTR∆?
# WT  22/37
# LTR∆ ODIRA (42/52, 81%) 
prop.test(x = c(22, 42), n = c(37, 52), alternative = "less")

#Is there is significant decrease in ODIRA in from WT to ARS∆?
# WT 22/37, ARS∆ 11/42
prop.test(x = c(22, 11), n = c(37, 42), alternative = "greater")

#Is there is significant decrease in ODIRA in from WT to ALL∆?
# WT 22/37, ARS∆ 11/42
prop.test(x = c(22, 12), n = c(37, 46), alternative = "greater")

###  Figure 4D  ####
# Dumbbell plots of select clones # 
db_clones = all_clones %>%
  filter(gap1_rounded>1) %>%
  filter(sample %in% c(2927,2929,3023, 3029, 3025, 2695, 3027, 2926,2931,3050,
                       2924,2690,3056,2948,2945,3069,2691,2852,2983,2946,3076,
                       3019,3022,3107,3099,3101,3093,2968)
  ) 

db_clones %>% 
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))),desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = cnv_start, xend=cnv_end, y = clone, color = Description))+
  geom_dumbbell(size=2, dot_guide=F, dot_guide_size=0.1)+
  coord_cartesian(xlim=c(360000,670191))+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype","LTR∆", "ARS∆","ALL∆") #third, now you can change legend labels
  )+

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
        #legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Strain"))

#### Figure 4E ##### 
### Scatterplot by Mechanism and Copy number, Facet by Description ##
all_clones %>%   
  filter(CNV_Mechanism %in% c("LTR NAHR", "ODIRA", "transposon-mediated", "aneuploid", "NAHR", "complex CNV"))%>% 
  mutate(CNV_Mechanism = fct_relevel(CNV_Mechanism, "LTR NAHR", "NAHR", "ODIRA", "transposon-mediated", "aneuploid", "complex CNV")) %>% #custom order
  ggplot(aes(x = distance_start, y=distance_end, color=CNV_Mechanism, shape = as.factor(gap1_rounded))) +
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

# Supplementary Figure S4A ####
## Grouped Boxplot of CNV Length by Strain & Generation 
sup_S4A_plot = all_clones %>% ggplot(aes(x=Description, y=cnv_length, fill = interaction(Description,generation))) +
  geom_boxplot()+ 
  geom_point(size = 0.5, 
             position=position_jitterdodge())+
  scale_y_continuous(trans='log10',
                     breaks = c(3000, 10000, 30000, 100000, 300000, 670000),
                     labels = c(3000, 10000, 30000, 100000, 300000, 670000)/1000
                       )+
  ylab("CNV length (kb)") +
  xlab("Strain")+
  theme_classic()+
  scale_fill_manual(values=c("gray","#6699cc","#e26d5c","#ffba08",
                             "gray40","#0a9396","#da4631","#faa307"))+
  theme(
    axis.text.x = element_text(size =12, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, vjust=2),
    text = element_text(size=16)
  )
sup_S4A_plot

# Two way ANOVA ##
summary(aov(log10(cnv_length)~Description*generation, data = all_clones))

# Supplementary Figure S4B - Types of ODIRA  ####
all_clones %>% select(ODIRAtype, gap1_rounded) %>% filter(ODIRAtype %in% c("ODIRA one end", "ODIRA", "ODIRA no ARS")) %>% arrange(gap1_rounded) %>% summarize(ODIRA_2 = sum(ODIRAtype == "ODIRA" & gap1_rounded == 2),
            ODIRA_3 = sum(ODIRAtype == "ODIRA" & gap1_rounded == 3),
            ODIRA_oneEnd_2 = sum(ODIRAtype == "ODIRA one end" & gap1_rounded == 2),
            ODIRA_oneEnd_3 = sum(ODIRAtype == "ODIRA one end" & gap1_rounded == 3),
            ODIRA_noARS_3 = sum(ODIRAtype == "ODIRA no ARS" & gap1_rounded==3)
            ) %>% pivot_longer(cols = ODIRA_2:ODIRA_noARS_3, names_to = "ODIRA_Type", values_to = "Count") %>% 
  mutate(Proportion = Count/sum(Count)) %>% 
  ggplot(aes(x = reorder(ODIRA_Type, -Proportion), y=Proportion)) +
  geom_col(fill = "#414833")+
  xlab("ODIRA Type")+
  ylab("Proportion")+
  theme_classic()
  
# Supplementary Table S4C - Count of Mechanism separated by strain ####
with(all_clones,table(Description, CNV_Mechanism)) %>% as.data.frame() %>%
  pivot_wider(names_from = Description, CNV_Mechanism, values_from = Freq) %>% 
  mutate(total = rowSums(.[, 2:5])) %>%
  adorn_totals("row")
