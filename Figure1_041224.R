# Chuong et al. 2024
# DNA replication errors are a major source of adaptive gene amplification
# Figure 1 and Supplementary Information 

# Load required packages
library(chromoMap)
library(tidyverse)
library(plotly)
library(scales)

# Figure 1A ####
chr11_feats = read.delim("Chromosome_XI_features.txt")
unique(chr11_feats$Feature.Type)

#format annotation file in the format for chromomap  
my_feats = chr11_feats %>% 
  filter(Feature.Type %in% c("telomere", "ARS", "tRNA gene", "long terminal repeat", "centromere") |
           Feature == "GAP1") %>% 
  separate(Coordinates,c("start", "end"), "-") %>% 
  mutate(name = rep("XI",50 )) %>%
  select(Feature, name, start, end, Feature.Type)
colnames(my_feats) <- NULL

#### Make a chromosome object as needed for chromoMap #
chrXI = data.frame(name = "XI",
                   start = 1,
                   end = 666816,
                   cent = 440129
)
colnames(chrXI)<-NULL

chromoMap(list(chrXI),
          list(my_feats),
          chr_color = c("#ede0d4"),
          #chr_width = 15,
          #chr_length = 10,
          #labels = T,
          y_chr_scale = 17, #bring ruler closer to chromosome 
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#F99B1C", "#B767AA", "#60C9DF", "#ede0d4","white", "black"))) #tRNA, ars, ltr, centromere, GAP1 ORF, telomere

# zoom in
zoom = my_feats[39:44,]
chrXI_zoom = data.frame(name = "XI", 
                        start = 512000, 
                        end = 520000)
colnames(chrXI_zoom) <- NULL
chromoMap(list(chrXI_zoom),
          list(zoom), 
          chr_width = 17,
          chr_length = 3,
          labels = T,
          chr_color = c("#ede0d4"),
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#F99B1C", "#60C9DF","white","#B767AA")))

# Figure 1B #####
#load CNV frequency data 
freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

#prep data for plot
med_freq_counts = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 166) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency),
         mad = mad(proportion),
         IQR = IQR(proportion))

#plot
med_freq_counts%>%
  filter(generation <= 137) %>%
ggplot(aes(x = generation, group = Description)) +
  geom_line(aes(y = med/100, color = Description), linewidth = 3) + 
  geom_ribbon(aes(y = med/100, ymin = med/100 - mad, ymax = med/100 + mad, fill = Description),alpha=0.3)+
  scale_color_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"),  #custom colors
                     limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"),
                     labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+
  scale_fill_manual(values=c("gray6", "#6699cc", "#e26d5c", "#DEBD52"), #custom colors
                    limits=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels=c("Wild type architecture", "LTR removed", "ARS removed", "LTR and ARS removed"))+#third, now you can change legend labels
  scale_x_continuous(limits=c(0,125))+
  xlab("Generation")+
  ylab("Median proportion of cells
with GAP1 CNV") +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 35),
        text = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=25), #change legend text font size
        axis.text.x = element_text(size = 40, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 44, color = "black"))

# Figure 1C ####

# Define frequency threshold for CNV appearance 
  #intuition: inflection generation before it goes up (vertically) 
freq_and_counts %>%
  filter(Description == "GAP1 WT architecture",
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample)

thresh = 10  # Defining it at 10% seems to be the best to capture this inflection point 

# CNV Appearance 
Tup_per_pop_10 = 
  freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thresh) %>%
  select(Type, Strain, Description, sample, generation, Gate, Frequency) %>%
  group_by(sample) %>%
  slice(which.min(generation))

Tup_per_pop_10$Description <- factor(Tup_per_pop_10$Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO", "GAP1 LTR + ARS KO")) #reorder boxplots

ggplot(Tup_per_pop_10, aes(Description, generation, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("gray", "#6699cc", "#e26d5c", "#DEBD52"))+ #change order of colors
  ylab("   Generation of first CNV appearance") +
  #scale_x_discrete(labels=c("Wildtype architecture","LTR removed","ARS removed","LTR and ARS removed"))+
  scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, max(Tup_per_pop_10$generation)))+
  theme_classic() +
  theme(legend.position = "none",
        #axis.text.x = element_text(size = 16, color = "black"), #edit x-tick labels
        axis.text.x = element_blank(), #remove x-tick labels 
        axis.ticks.x=element_blank(), #remove x-ticks 
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 18),
        text = element_text(size=16))+
  geom_jitter(size = 2, alpha = 0.8, 
              color = c(rep("black", 5),  #wildtype, 5, gray
                        rep("#D9BB59", 8),  #LTR and ARS gold
                        rep("#e26d5c", 7),  #ARS, 7, salmon 
                        rep("#6699cc", 7)  #LTR,7, #blue
              ))

shapiro.test(Tup_per_pop_10$generation) #not normal

# Instead of ANOVA, do Krusal-Wallis test (non-parametric) 
kruskal.test(generation~Description, data = Tup_per_pop_10)

# Instead of pairwise t-tests, do pairwise Wilcoxon Mann-Whitney with Bonferroni correction
pairwise.wilcox.test(Tup_per_pop_10$generation, Tup_per_pop_10$Description, p.adjust.method = "bonferroni")

# Figure 1D #####
# Calculate CNV selection # 
#Compute natural log proportion of each population with CNV relative to that without CNV
ln_table = freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency),
         prop_NoCNV = 100-prop_CNV,
         CNV_NoCNV = prop_CNV/prop_NoCNV,
         logECNV_NoCNV = log(CNV_NoCNV))

pop_list = unique(ln_table$sample) %>% sort()
gens = unique(ln_table$generation)

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(summary(x)[4]$coefficients[2], digits = 4),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(slope == b~~~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));
}

#function to apply for loop to each of 28 populations
sliding_fit = function(num_fitpoints, population){
  timepoints = nrow(subset(ln_table, sample %in% c(population)))
  rounds = timepoints - num_fitpoints +1
  m <- matrix(ncol = 7, nrow = rounds) #nrow = number of iterations. number of iterations depend on the number of generations and the number of fitpoints. max num of generations = 24. minimum num of fitpoints is 2. therefore nrow max is 23.
  colnames(m) <- c("start", "end", "gen_start", "gen_end", "slope", "rsquared", "sample")
  start = 1
  end = num_fitpoints
  for (i in 1:rounds){
    print(i)
    pop_data <- subset(ln_table, sample %in% c(population))
    fit_points_df <- subset(ln_table, sample %in% c(population) & generation >= gens[start] & generation <= gens[end])
    if (is.na(gens[end]) == TRUE ){
      break
    }
    fit <- lm(logECNV_NoCNV ~ generation, fit_points_df) #linear model, lm(y~x, by the data)
    print(summary(fit))
    
    ggplot(pop_data, aes(x=generation,y=(as.numeric(logECNV_NoCNV)), colour=sample)) +
      geom_point() +
      geom_smooth(data=fit_points_df, method=lm, show.legend=FALSE) +
      scale_y_continuous(expand = c(0, 0), 'ln(Prop. CNV/Prop. non-CNV)', limits = c(min(pop_data$logECNV_NoCNV)-1, max(pop_data$logECNV_NoCNV)+1)) +
      annotate("text", x = 200, y = min(pop_data$logECNV_NoCNV)-0.5, label = equation(fit), parse = TRUE) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5), "Generations", limits=c(0,260)) +
      theme_classic() +
      scale_color_manual(values = c('black')) +
      guides(colour = guide_legend(override.aes = list(size=2))) +
      theme(legend.position = c(.15,.95), plot.title = element_text(size=14, hjust = 0.5), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
    
    ggsave(paste0(population,"_Sup_g",gens[start],"-",gens[end],"_",num_fitpoints,"pts.png"), width = 8, height = 5)
    
    print(paste0("From timepoints ",start," to ", end, ", generations ", gens[start], " to ", gens[end],", slope was ", as.numeric(coef(fit)[2]) %>% round(4)," and rsquared was ", as.numeric(summary(fit)[8]) %>% round(2) ))  
    
    m[i,1] <- start
    m[i, 2]<- end
    m[i, 3] <- gens[start]
    m[i,4] <- gens[end]
    m[i, 5] <- as.numeric(coef(fit)[2]) %>% round(4)
    m[i,6] <- as.numeric(summary(fit)[8]) %>% round(2)
    m[i, 7] <- population
    
    start = start + 1
    end = end + 1
  }
  m = m %>% na.omit() #remove NAs
  write_csv(as.data.frame(m), paste0(population,"_fits","_",num_fitpoints,"pts.csv"))
  return(m)
}

# call the function for all pops in the pop_list using map()
map(.x = pop_list[1:27], ~sliding_fit(4, .x))

# Pull in the fits tables per population and merged into one
slopes = list.files(path = ".", pattern = paste0("_fits","_",4,"pts")) %>%
  read_csv() %>%
  write_csv(file = "Sup_fits_4_pts_all_pops.csv")

# CNV Selection Boxplot 
meta = Tup_per_pop_10 %>% select(Description, sample)

Sup = slopes %>%
  right_join(meta) %>%
  group_by(sample) %>%
  mutate(slope = max(slope)) %>%
  ungroup() %>% 
  select(sample, slope, Description) %>% 
  distinct()

ggplot(Sup, aes(Description, slope, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("gray", "#6699cc", "#e26d5c", "#DEBD52"))+ #change order of colors
  ylab("Percent of increase in 
  CNVs per generation") +
  #scale_x_discrete(labels=c("Wildtype architecture","LTR removed","ARS removed","LTR and ARS removed"))+
  scale_y_continuous(labels = scales::label_number(scale = 100))+
  theme_classic() +
  theme(#plot.margin = unit(c(.5, .5, .5, .5), "cm"),
    legend.position = "none",
    #axis.text.x = element_text(size = 16, color = "black"), #edit x-tick labels
    axis.text.x = element_blank(), #remove x-tick labels 
    axis.ticks.x=element_blank(), #remove x-ticks 
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, vjust=2),
    text = element_text(size=16))+
  geom_jitter(size = 3, alpha = 0.8, 
              color = c(rep("black", 5),  #wildtype, 5, gray
                        rep("#D9BB59", 8),  #LTR and ARS gold
                        rep("#e26d5c", 7),  #ARS, 7, salmon 
                        rep("#6699cc", 7)  #LTR,7, #blue
              ))

hist(Sup$slope)
shapiro.test(Sup$slope) #normal distribution
Sup_anova = aov(slope~Description, data = Sup)
summary(Sup_anova) # p = 0.00318
pairwise.t.test(Sup$slope, Sup$Description, p.adjust.method = "bonferroni")

# Figure 1E #####
# Boxplot - CNV equilibrium, the other inflection point where the line pleateaus #
gen_maint = slopes %>% 
  right_join(meta) %>%
  group_by(sample) %>%
  filter(gen_start > 50, slope < 0.005) %>%
  slice(which.min(gen_start)) %>% 
  select(-start, -end)

ggplot(gen_maint, aes(Description, gen_start, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("gray", "#6699cc", "#e26d5c", "#DEBD52"))+ #change order of colors
  ylab("Generation of 
  CNV maintainence") +
  #scale_x_discrete(labels=c("Wildtype architecture","LTR removed","ARS removed","LTR and ARS removed"))+
  #scale_y_continuous(labels = label_number(scale = 100))+
  theme_classic() +
  theme(#plot.margin = unit(c(.5, .5, .5, .5), "cm"),
    legend.position = "none",
    #axis.text.x = element_text(size = 16, color = "black"), #edit x-tick labels
    axis.text.x = element_blank(), #remove x-tick labels 
    axis.ticks.x=element_blank(), #remove x-ticks 
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 18, vjust=2),
    text = element_text(size=16))+
  geom_jitter(size = 3, alpha = 0.8, 
              color = c(rep("black", 4),  #wildtype, 5, gray
                        rep("#D9BB59", 8),  #LTR and ARS gold
                        rep("#e26d5c", 6),  #ARS, 7, salmon 
                        rep("#6699cc", 7)  #LTR,7, #blue
              ))

hist(gen_maint$gen_start)
shapiro.test(gen_maint$gen_start) #normal
summary(aov(gen_start~Description, gen_maint))
pairwise.t.test(gen_maint$gen_start, gen_maint$Description, p.adjust.method = "bonferroni")
