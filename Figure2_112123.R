# Forked from gating_112122.R
# Figure 2 R script, Chuong et al. 2023 
# 10/26/23making new CNV lineplots because the 4 wildtype lines are under-reporting GAP1 CNVs due to broken CNV reporter subpop
# Make the order of lines showing up: Wildtype, ARS removed, LTRs removed, LTR and ARS removed. 

# 11/21/23 - Adding back in the 4 wildtype lines! Because we include them in the SBI. 
# We changed the model to include a broken reporter parameter, p0. Since it's included in SBI,
# it's best to keep them in my own population dynamics and just explain it in the text. 
# Also, having more than 1 wildtype population lets me do pairwise t-tests. 

# Load required packages
library(tidyverse)
library(plotly)
library(scales)

setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

#load frequency data
freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

#### 10/26/23 add the alpha-fade to lineplots because people always ask what the variance is between populations ####
# https://plotly.com/r/line-charts/ 
# Plotly resource Pieter showed me! Also called filled line charts 
med_freq_counts = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 203) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency))

#subset by genotype for easier graphing
wt = med_freq_counts %>% dplyr::filter(Description == "GAP1 WT architecture")
arsko = med_freq_counts %>% dplyr::filter(Description == "GAP1 ARS KO")
ltrko = med_freq_counts %>% dplyr::filter(Description == "GAP1 LTR KO")
allko = med_freq_counts %>% dplyr::filter(Description == "GAP1 LTR + ARS KO")

# get lows and highs (min and max per timepoint) per genotype
wtLow = wt %>% group_by(generation) %>% summarize(mins = min(Frequency)/100)
wtHigh = wt %>% group_by(generation) %>% summarize(maxs = max(Frequency)/100)
wtMeds =  wt %>% group_by(generation) %>% summarise(median = max(med)/100)
low = ltrko %>% group_by(generation) %>% summarize(mins = min(Frequency)/100)
high = ltrko %>% group_by(generation) %>% summarize(maxs = max(Frequency)/100)
medians = ltrko %>% group_by(generation) %>% summarise(median = max(med)/100)
lowars = arsko %>% group_by(generation) %>% summarize(mins = min(Frequency)/100)
highars = arsko %>% group_by(generation) %>% summarize(maxs = max(Frequency)/100)
medars = arsko %>% group_by(generation) %>% summarise(median = max(med)/100)
alllow = allko %>% group_by(generation) %>% summarise(mins = min(Frequency)/100)
allhigh = allko %>% group_by(generation) %>% summarise(maxs = max(Frequency)/100)
allMedian = allko %>% group_by(generation) %>% summarise(median = max(med)/100)

#graph median line for Wildtype architecture and faded range for repeat evolutions (n = 5)
fig = plot_ly(wtHigh, x = wtHigh$generation, y = wtHigh$maxs, type = 'scatter', mode = 'lines', line = list(color = 'transparent'), showlegend = FALSE, name = "WT maxs") %>%
  add_trace(y=~wtLow$mins, type = 'scatter', mode = 'lines', fill = 'tonexty', fillcolor = 'lightgray',line = list(color = 'transparent'), showlegend = FALSE, name = "WT mins") %>%
  #tonexty is fills to next y (shades under y)
  add_trace(x=wtMeds$generation, y=wtMeds$median, type = 'scatter', mode = 'lines', line = list(color = 'black', width = 10), showlegend = FALSE, name = "Wildtype architecture")%>%
  layout(paper_bgcolor='white', plot_bgcolor='white',
                      xaxis = list(
                        title = list(text="Generation",
                                     font=list(size = 50, family="Arial", color = "black")),
                        dtick = 50, 
                        tick0 = 0, 
                        zerolinecolor = "black",
                        showgrid = TRUE,
                        showline = TRUE,
                        showticklabels = TRUE,
                        tickfont = list(size = 55, family = "Arial", color = "black"),
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = TRUE),
                      yaxis = list(
                        title = list(text="Proportion of cells with GAP1 CNV",
                                     font=list(size=40, family="Arial", color = "black")),
                        dtick = 0.25,
                        tick0 = 0,
                        zerolinecolor = "black",
                        showgrid = TRUE,
                        showline = TRUE,
                        showticklabels = TRUE,
                        tickfont = list(size = 60, family="Arial",color = "black"),
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE))
fig
#save Wildtype only plot
save_image(fig, file = "propCNV_WT_fade_112123_1200x800.png", scale = 1, width = 1200, height = 800)
save_image(fig, file = "propCNV_WT_fade_112123_1200x800.pdf", scale = 1, width = 1200, height = 800)

# add ARS KO median line and fade range
fig = fig %>% add_trace(highars, x = highars$generation, y = highars$maxs, type = 'scatter', mode = 'lines', line = list(color = 'transparent'), showlegend = FALSE, name = "ARS maxs")
fig = fig %>% add_trace(y=~lowars$mins, type = 'scatter', mode = 'lines', fill = 'tonexty', fillcolor = 'rgba(226, 109, 92, 0.2)',line = list(color = 'transparent'), showlegend = FALSE, name = "ARS mins") #tonexty is fills to next y (shades under y) 
fig = fig %>% add_trace(x=~medars$generation, y = ~medars$median, type = 'scatter', mode = 'lines', line = list(color='rgba(226, 109, 92, 1)',width=10), showlegend = FALSE,name = 'ARS removed')
fig

#save WT and ARS plot
# save_image(fig, file = "propCNV_WT+ARS_fade_102723_1200x800.png", scale = 1, width = 1200, height = 800)
# save_image(fig, file = "propCNV_WT+ARS_fade_102723_1200x800.pdf", scale = 1, width = 1200, height = 800)

# add LTR KO median line and fade range
fig = fig %>% add_trace(high, x = high$generation, y = high$maxs, type = 'scatter', mode = 'lines', line = list(color = 'transparent'), showlegend = FALSE, name = "LTRKO max")
fig = fig %>% add_trace(y= ~low$mins, type = 'scatter', mode = 'lines',
                         fill = 'tonexty', fillcolor = 'rgba(203, 221, 239,0.5)', line = list(color = 'transparent'), showlegend = FALSE, name = "LTRKO Low") #BBD3EA #color='#CBDDEF', #RGBA the 4th number is the opacity! 
fig <- fig %>% add_trace(x = medians$generation, y = medians$median, type = 'scatter', mode = 'lines',
                         line = list(color='#6699cc', width = 10),
                         showlegend = FALSE,
                         name = 'LTR removed')
fig

#save WT+ARS+LTR plot
# save_image(fig, file = "propCNV_WT+ARS+LTR_fade_102723_1200x800.png", scale = 1, width = 1200, height = 800)
# save_image(fig, file = "propCNV_WT+ARS+LTR_fade_102723_1200x800.pdf", scale = 1, width = 1200, height = 800)


# add ALL KO median line and fade region
fig = fig %>% add_trace(y = ~allhigh$maxs, type = 'scatter', mode = 'lines', line = list(color = 'transparent', showlegend = FALSE, name = "ALLko max"))
fig = fig %>% add_trace(y = ~alllow$mins, type = 'scatter', mode = 'lines', fill = "tonexty", fillcolor = 'rgba(222, 189, 82, 0.3)', line = list(color='transparent'), showlegend = FALSE, name = "All KO mins")
fig = fig %>% add_trace(y = ~allMedian$median, type = 'scatter', mode = 'lines', line = list(color="#DEBD52", width = 10), showlegend = FALSE, name ="LTR and ARS removed")
fig
fig <- fig %>% layout(paper_bgcolor='white', plot_bgcolor='white',
                      xaxis = list(
                        title = list(text="Generation",
                                    font=list(size = 50, family="Arial", color = "black")),
                        dtick = 50, 
                        tick0 = 0, 
                        zerolinecolor = "black",
                        showgrid = TRUE,
                        showline = TRUE,
                        showticklabels = TRUE,
                        tickfont = list(size = 55, family = "Arial", color = "black"),
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = TRUE),
                      yaxis = list(
                        title = list(text="Proportion of cells with GAP1 CNV",
                                     font=list(size=40, family="Arial", color = "black")),
                        nticks= 5,
                        dtick = 0.25,
                        tick0 = 0,
                        zerolinecolor = "black",
                        showgrid = TRUE,
                        showline = TRUE,
                        showticklabels = TRUE,
                        tickfont = list(size = 60, family="Arial",color = "black"),
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE))
fig

#save as static image
#reticulate::use_miniconda('r-reticulate')

#Plot with wildtype only
save_image()


#Plot with all four genotypes

# save_image(fig, file = "propCNV_full_fade_112123_1200x800.png", scale = 1, width = 1200, height = 800)
# save_image(fig, file = "propCNV_full_fade_112123_1200x800.pdf", scale = 1, width = 1200, height = 800)

######## Boxplots for populations dynamics #########
clean_freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

# Define frequency threshold for CNV appearance 

#intuition: try to find this inflection generation before it goes up (vertically) 
clean_freq_and_counts %>%
  filter(Description == "GAP1 WT architecture",
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample) %>% View()

# Second definition: using the one and two copy controls to define a false positive rate. 
# this second definition is NOT compatible for our analysis, since our controls are not good controls for CNVs for our genotypes. The controls are probabaly ONLY good controls for the Wildtype architecture populations since they share the same strain background. 

#Determine Tup for each population, the generation when CNVs first appear above the threshold frequency.
thresh = 10  # Defining it at 10% seems to be the best to capture this inflection point and is consistent with the story displayed in the median lineplot
Tup_per_pop_10 = 
  clean_freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thresh) %>%
  #filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= 15) %>%
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
#ggsave("Tup_boxplot_121622.png", bg = "#FFFFFF", height = 6, width = 10)
#ggsave("Tup_boxplot_121622.pdf", bg = "#FFFFFF", height = 6, width = 10)
# ggsave("Tup_boxplot_121622_3x4_smallFont.png", bg = "#FFFFFF", height = 3, width = 4)
# ggsave("Tup_boxplot_112123.png" , bg = '#FFFFFF', height = 2.5, width = 2.3)
# ggsave("Tup_boxplot_112123_2.pdf" , bg = '#FFFFFF', height = 2.5, width = 2.3)

# ANOVA to test for significance
# One Way ANOVA because there is only 1 independent variable, genotype.
# Anova Tut: https://www.scribbr.com/statistics/anova-in-r/
Tup_anova = aov(generation~Description, data = Tup_per_pop_10)
summary(Tup_anova)
#               Df Sum Sq Mean Sq F value   Pr(>F)    
#  Description  3   2122   707.4    8.84 0.000442 ***
#  Residuals   23   1840    80.0 
# Conclusion: There IS a significant difference in the means. Genotype has has significant effect on Tup.

hist(Tup_per_pop_10$generation) 
shapiro.test(Tup_per_pop_10$generation) #not normal
# Shapiro-Wilk normality test
# data:  Tup_per_pop_10$generation
# W = 0.90148, p-value = 0.01446
# Null distribution is normality 
# Because p-value is less than 0.05, the distribution is NOT normal

##### Instead of ANOVA, do Krusal-Wallis test (non-parametric) #####
kruskal.test(generation~Description, data = Tup_per_pop_10)
# Kruskal-Wallis rank sum test
# data:  generation by Description
# Kruskal-Wallis chi-squared = 15.577, df = 3, p-value = 0.001384
# Conclusion - Genotype had significant effect on generation of CNV appearance, but we don't know WHICH genotype has the effect

#### Tup is not a normal distribution. Instead of pairwise t-tests, do pairwise Wilcoxon Mann-Whitney with Bonferroni correction ###### 
# Try pairwise t-tests anyway for fun
pairwise.t.test(Tup_per_pop_10$generation, Tup_per_pop_10$Description, p.adjust.method = "bonferroni")
# Pairwise comparisons using t tests with pooled SD 
# data:  Tup_per_pop_10$generation and Tup_per_pop_10$Description 
#                   GAP1 WT architecture GAP1 LTR KO GAP1 ARS KO
# GAP1 LTR KO       1.00000              -           -          
# GAP1 ARS KO       0.01008              0.00036     -          
# GAP1 LTR + ARS KO 1.00000              0.09781     0.12857 

pairwise.wilcox.test(Tup_per_pop_10$generation, Tup_per_pop_10$Description, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# P value adjustment method: bonferroni 
# data:  Tup_per_pop_10$generation and Tup_per_pop_10$Description 
#                 GAP1 WT architecture GAP1 LTR KO GAP1 ARS KO
# GAP1 LTR KO       1.0000               -           -          
# GAP1 ARS KO       0.2765               0.0059      -          
# GAP1 LTR + ARS KO 1.0000               0.0124      0.2548     




# Calculate Sup
# "Sup is the rate of increase in CNV abundance during the initial expansion of the CNV subpopulation" Lauer et al. 2018
# S1 Text. Calculation of CNV dynamics parameters.
# Graphic representation of linear fit (and corresponding R2 values)
# during initial population expansion of CNV alleles.
# Slope of the linear fit corresponds to the dynamics parameter Sup shown in
# Table 1 and was calculated for the original evolution experiment and the
# barcode experiment. Data and code used to generate these figures can be
# accessed in OSF: https://osf.io/fxhze/. CNV, copy number variant.

#Calculate natural log proportion of each population with CNV relative to that without CNV

## well  maybe I Need to draw new gates, since I have to calculate prop of CNV relative to noCNV
# and I basically believe that my propCNv for WT and LTR are false positives.
# need gates that make the gen 0-50 smooth smooth! for WT and LTR KO
ln_table = clean_freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency),
         prop_NoCNV = 100-prop_CNV,
         CNV_NoCNV = prop_CNV/prop_NoCNV,
         logECNV_NoCNV = log(CNV_NoCNV)) #log() function is natural logarithm in R (even though  log() commonly thought as base10 )

# JULIE: A function to calculate Sup, Explained Variance, make graphs, ggsave graphs
# then, use map() to apply function to all 27 populations - I have 27
# as part of the function, include a variable that changes number of fit points (window width)
# I recycled some code from Lauer et al. 2018

pop_list = unique(ln_table$sample) %>% sort()
gens = unique(ln_table$generation)

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(summary(x)[4]$coefficients[2], digits = 4),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(slope == b~~~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));
}

#function to apply for loop to each of 28 populations, used some code from Lauer et al. 2018
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

# Test call the function and assign outputted matrix to a variable
result <- sliding_fit(4, pop_list[27]) #population 27
assign(paste0(pop_list[27],"_fits","_",4,"pts"), result) #Use assign() to rename and save to R environment

#call the function for all pops in the pop_list using map()
map(.x = pop_list[1:27], ~sliding_fit(4, .x)) #tilde inside functions (https://stackoverflow.com/questions/70665707/what-is-the-meaning-of-and-inside-the-function-map)

summary(fit)$coef[[4]] #fourth coefficient is the standard error of the linear model slope

# Pull in the fits tables per population and merged into one

# list.files(path = ".", pattern = paste0("_fits","_",4,"pts")) %>% 
#   read_csv() %>%
#   write_csv(file = "Sup_fits_4_pts_all_pops.csv")

slopes = read_csv("Sup_fits_4_pts_all_pops.csv")

meta = Tup_per_pop_10 %>% select(Description, sample)

Sup = slopes %>%
  right_join(meta) %>%
  group_by(sample) %>%
  mutate(slope = max(slope)) %>%
  ungroup() %>% 
  select(sample, slope, Description) %>% 
  distinct()

#write_csv(Sup, "Sup_121722.csv")

#Sup boxplot 
Sup

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

# ggsave("Sup_boxplot_112123.png" , bg = '#FFFFFF', height = 2.5, width = 2.3)
# ggsave("Sup_boxplot_112123.pdf" , bg = '#FFFFFF', height = 2.5, width = 2.3)

hist(Sup$slope)
shapiro.test(Sup$slope) #normal distribution
# Shapiro-Wilk normality test
# data:  Sup$slope
# W = 0.95825, p-value = 0.337   
# Conclusion Slopes follow a normal distribution
Sup_anova = aov(slope~Description, data = Sup)
summary(Sup_anova) # p = 0.00318

pairwise.t.test(Sup$slope, Sup$Description, p.adjust.method = "bonferroni")
#                   GAP1 WT architecture GAP1 LTR KO GAP1 ARS KO
# GAP1 LTR KO       1.0000               -           -          
# GAP1 ARS KO       1.0000               0.0577      -          
# GAP1 LTR + ARS KO 0.1635               0.0026      1.0000



#### Calculate Tdown ####

# No need to do this because only a few pops exhibit CNV loss, so it doesn't merit a whole figure.
# Mention it in the text

##### Calculate Sdown ####

# No need to do this because only a few pops exhibit CNV loss, so it doesn't merit a whole figure

#### Calculate, Tmaint, early generation of the CNV maintenance phase, the other inflection point where the line goes flat (horizontal) ####
# ie) what generation gives lowest slope values for a consecutive 4 timepoints .. 
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

# ggsave("Tmaint_boxplot_112123.png", bg = "#FFFFFF", height = 2.5, width = 2.3)
# ggsave("Tmaint_boxplot.112123_2.pdf", bg = "#FFFFFF", height = 2.5, width = 2.6)

hist(gen_maint$gen_start)
shapiro.test(gen_maint$gen_start) #normal

summary(aov(gen_start~Description, gen_maint))
# pval = 0.00833

pairwise.t.test(gen_maint$gen_start, gen_maint$Description, p.adjust.method = "bonferroni")
#                   GAP1 WT architecture GAP1 LTR KO GAP1 ARS KO
# GAP1 LTR KO       1.000                -           -          
# GAP1 ARS KO       0.050                0.083       -          
# GAP1 LTR + ARS KO 0.060                0.098       1.000 