####################################################
# STEP 9:  Quantify CNV dynamics (Lauer et al. 2018)
# Author: Julie
# 1) First, calculate Tup, the generation at which CNVs are initially detected, (Lang et al. 2011 and Lauer et al. 2018)
# To do that, calculate the false positive rate for CNV detection (threshold), which I will define as the median frequency of 1 copy control cells appearing in the two_copy_or_more gate and gate across generations 8-260 plus the interquartile range (IQR) if the distribution of one copy controls appearing in the CNV gate as NOT normal. If the distribution is normal, then use the mean plus one standard deviation like in Lauer et al. 2018. Like in Lauer et al. 2018, samples surpassing this threshold is considered to contain CNVs.

#CNV False Positive Rate is defined by the frequency of the 1 copy control strain appearing in the CNV gate which is called the Two_or_more copy gate.
CNV_false_pos_df = #freq %>%
  freq_and_counts %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "1_copy_ctrl") %>%
  filter(Gate %in% c("two_or_more_copy")) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count)

# draw a histogram, see if distribution is normal by eye
hist(CNV_false_pos_df$Frequency) #looks normal but could be left skewed
abline(v = median(CNV_false_pos_df$Frequency),col = "red",lwd = 1.5)
abline(v = mean(CNV_false_pos_df$Frequency), col = "blue", lwd = 1.5)

#Test for normality
shapiro.test(CNV_false_pos_df$Frequency) #null hypothesis is that the distribution is normal. if p <0.05, then it rejects the null hypothesis and so the distribution is NOT normal.
#W = 0.95715, p-value = 0.4886
# Null hypothesis cannot be rejected. Our distribution is normal. Use the mean + 1SD as the threshold value.
# mean = 4.29, sd = 2.32
thres_mean = mean(CNV_false_pos_df$Frequency) + sd(CNV_false_pos_df$Frequency) #6.615341

#Determine Tup for each population, the generation when CNVs first appear above the threshold frequency.
Tup_per_pop = #freq %>%
  freq_and_counts %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thres_mean) %>%
  select(Type, Strain, Description, sample, generation, Gate, Frequency) %>% #View()
  group_by(sample) %>%
  slice(which.min(generation))

#Tup_per_pop %>% write_csv(file = "01_02_04_v2_fw_Tup_per_pop.csv")
#Tup_per_pop = read_csv(file = "01_02_04_v2_fw_Tup_per_pop.csv")

ggplot(Tup_per_pop, aes(reorder(Description, -generation),generation, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("#e26d5c", "#DEBD52", "#6699cc", "gray"))+
  ylab("Generation of first CNV appearance") +
  scale_x_discrete(labels=c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))+
  scale_y_continuous(breaks=c(8,20, 30, 40, 50, 60, max(Tup_per_pop$generation)))+
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(family="Arial", size = 16, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 20, color = "black"),
        text = element_text(size=18))+
  geom_jitter(size = 2, alpha = 0.9, color = c(
    rep("black",5), #wildtype, 5, gray
    "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR and ARS gold
    "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", #ARS, 7, softer salmon repeats
    rep("black", 8) #LTR,8, #blue
  ))
#ggsave("01_02_04_v2_fw_Tup_boxplot_042022.png")
#ggsave("01_02_04_v2_fw_Tup_boxplot_051022.png")

# ANOVA to test for significance
# One Way ANOVA because there is only 1 independent variable, genotype.
# Anova Tut: https://www.scribbr.com/statistics/anova-in-r/
Tup_anova = aov(generation~Description, data = Tup_per_pop)
summary(Tup_anova)
#            Df  Sum Sq Mean Sq F value  Pr(>F)
#Description  3   4483  1494.3   6.772 0.00181 **
#Residuals   24   5296   220.7
# Conclusion: There IS a significant difference in the means. Genotype has has significant effect on Tup.

#Calculate Sup
# "Sup is the rate of increase in CNV abundance during the initial expansion of the CNV subpopulation" Lauer et al. 2018
# S1 Text. Calculation of CNV dynamics parameters.
# Graphic representation of linear fit (and corresponding R2 values)
# during initial population expansion of CNV alleles.
# Slope of the linear fit corresponds to the dynamics parameter Sup shown in
# Table 1 and was calculated for the original evolution experiment and the
# barcode experiment. Data and code used to generate these figures can be
# accessed in OSF: https://osf.io/fxhze/. CNV, copy number variant.

#Calculate natural log proportion of each population with CNV relative to that without CNV
ln_table = freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency),
         prop_NoCNV = 100-prop_CNV,
         CNV_NoCNV = prop_CNV/prop_NoCNV,
         logECNV_NoCNV = log(CNV_NoCNV)) #log() function is natural logarithm in R (even though  log() commonly thought as base10 )

# JULIE: A function to calculate Sup, Explained Variance, make graphs, ggsave graphs
# then, use map() to apply function to all 28 populations - I have 28
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
  rounds = timepoints- (timepoints/num_fitpoints) - 1
  m <- matrix(ncol = 6, nrow = rounds) #nrow = number of iterations. number of iterations depend on the number of generations and the number of fitpoints. max num of generations = 24. minimum num of fitpoints is 2. therefore nrow max is 23.
  colnames(m) <- c("start", "end", "gen_start", "gen_end", "slope", "rsquared")
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

    print(paste0("From timepoints ",start," to ", end, ", generations ", gens[start], " to ", gens[end],", slope was ", as.numeric(coef(fit)[2]) %>% round(4)," and rsquared was ", as.numeric(summary(fit)[8]) %>% round(2) ))  #populate a data frame? five columns: start, end, gen_start, gen_end,  Rsq.

    m[i,1] <- start
    m[i, 2]<- end
    m[i, 3] <- gens[start]
    m[i,4] <- gens[end]
    m[i, 5] <- as.numeric(coef(fit)[2]) %>% round(4)
    m[i,6] <- as.numeric(summary(fit)[8]) %>% round(2)

    start = start + 1
    end = end + 1
  }
  m = m %>% na.omit() #remove NAs
  write_csv(as.data.frame(m), paste0(population,"_fits","_",num_fitpoints,"pts.csv"))
  return(m)
}

#call the function and assign outputted matrix to a variable
result <- sliding_fit(4, pop_list[27])
assign(paste0(pop_list[27],"_fits","_",4,"pts"), result) #Use assign() to rename and save to R environment

#call the function for all pops in the pop_list using map()
map(.x = pop_list[1:28], ~sliding_fit(4, .x)) #tilde inside functions (https://stackoverflow.com/questions/70665707/what-is-the-meaning-of-and-inside-the-function-map)

summary(fit)$coef[[4]] #fourth coefficient is the standard error of the linear model slope


####### MY PALLETTE

#Gold Metallic (6)
#DEBD52
#DBB741
#D7B02F
#CAA426
#D9BB59
#B89523

#NEW PALLETE 4/20/22
#WT = gray  "gray","gray","gray","gray","gray",
#LTR and ARS = GOLD metalic "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR,8,gold
#ARS = SALMON  "#e26d5c"
#LTR = BABY BLUE "#6699cc"
