#Cytoexplorer analysis that I will not use for the publications.
# Median normalized fluoresence lineplots over time combining the replicate populations for each genotype.
# But want to keep because it's useful code to know how to combine replicate lineplots.
# In our LTEE experiments, we do not combine replicate populations because they are there own population with independent evolutionary trajectories.

# Combine the replicates/populations and plot median normalized fluorescence over time
# dashed gray controls lines on top of the experiment lineplot
clean_adj_norm_medians %>%
  filter(Type == "Experimental") %>%
  filter(!(Med_B2A_FSC<1.5 & Type == "Experimental")) %>% #filter out outliers (likely resulting from contamination)
  ggplot(mapping = aes(generation, Med_B2A_FSC, color = Description)) +
  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill = Description), alpha=0.3) + #experimentals loess regression with standard error cloud
  geom_line(mapping = aes(generation, Med_B2A_FSC, color = Type, linetype = Type),
            data = clean_adj_norm_medians %>% filter(Type != "Experimental")) + #controls lineplot
  scale_linetype_manual(values = c("dashed", "dashed", "dashed"))+ #dashed lines for controls
  scale_color_manual(values=c("gray", "gray", "gray","#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+ #line color
  scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+ #standard error cloud color
  facet_wrap(~Description) +
  theme_classic() +
  ggtitle("loess regression of combined populations") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Generation") +
  ylab("Median normalized fluorescence (a.u.)") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(legend.position = "none",
        text = element_text(size=12),
        strip.background = element_blank(), #remove box around facet title
        strip.text = element_text(size=12),
        axis.text.x = element_text(family="Arial", size = 12, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 12, color = "black")) #I did it!

ggsave("loes_regression_MedNormFluo_011222.png")


#Plot proportion of the populations with a CNV over time (collapse the replicates)
#freq %>%
fw_freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  #ggplot(aes(Frequency)) + geom_histogram(bins=50) #right skewed distribution
  ggplot(aes(generation, Frequency, color = Description)) +
  scale_color_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+
  scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+
  #geom_point(alpha = 0.5, size =1) +
  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=Description), alpha=0.3) +
  facet_wrap(~Description) +
  theme_minimal() +
  ggtitle("loess regression of the populations") +
  ylab("Proportion of the population with GAP1 CNV") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  scale_y_continuous(breaks=seq(0,100,25))+
  theme(text = element_text(size=12),
        legend.position = "none",
        axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 10, color = "black"))
