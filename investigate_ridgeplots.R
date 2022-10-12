###### RANDOM INVESTIGATIONS ############

# Note: These variables in this Rscript come from workflow.R'''

#### Investigate the zig zaggy timepoints by graphing the ridgeplots for those populations to see what the distribution is like. (Zigzaggy lines of any one population can seen in the plot_list plots)
#my idea is maybe the distribution shape can tel whether it's CNV dynamics or contamination.
sc_gap1_4 = read.csv(file = "sc_distributions_gap1_4_all_timepoints.csv", stringsAsFactors = T) %>%
  mutate(generation = factor(generation, levels = unique(generation)))

sc_gap1_4 %>%
  ggplot(aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
  xlab("Normalized fluorescence (a.u.)") +
  ylab("Generation") +
  ggtitle("GAP1 WT 4") +
  theme_classic() +
  scale_x_continuous(limits=c(0.0,2.5), breaks = c(0, 1, 2, 2.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 2.5))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
ggsave("gap1_4_ridgeplot_scale2.png")

sc_gap1_4 %>%
  group_by(generation) %>%
  mutate(norm_median = median(B2A_FSC)) %>%
  distinct() %>%
  ggplot(aes(generation, norm_median, group = 1)) +
  geom_line()+
  geom_point()+
  ggtitle("gap1_4")+
  ylab("normalized B2A/FSC then median")
theme_classic()
ggsave("gap1_4_normalized_median_lineplot.png")

sc_gap1_4 %>%
  group_by(generation) %>%
  mutate(med_B2A = median(B2.A)) %>% View()
ggplot(aes(generation, med_B2A, group = 1)) +
  geom_line() +
  geom_point()+
  ggtitle("gap1_4") +
  ylab("median B2-A fluorescence (a.u)") +
  theme_classic()
ggsave("gap1_4_raw-B2A_lineplot.png")

sc_gap1_4 %>%
  group_by(generation) %>%
  mutate(med_FSC = median(FSC.A)) %>%
  ggplot(aes(generation, med_FSC, group = 1)) +
  geom_line() +
  geom_point()+
  ggtitle("gap1_4") +
  ylab("median FSC fluorescence (a.u)") +
  theme_classic()
ggsave("gap1_4_raw-FSC_lineplot.png")

clean_adj_norm_medians %>% filter(sample == "gap1_4") %>%
  ggplot(aes(generation, Med_B2A_FSC, group = 1))+
  geom_line()+
  geom_point()+
  theme_classic()
ggsave("gap1_4_median_then_norm_lineplot.png")

clean_adj_norm_medians %>%
  filter(sample == "gap1_4", generation > 120) %>%
  select(sample, Description, generation, `FSC-A`, GFP, Med_B2A_FSC, Count) %>% View()

#experimental populations that dip down, investigate.
dips = clean_adj_norm_medians %>%
  filter(Med_B2A_FSC < 1.5 & Type == "Experimental") %>%
  write_csv("Med_B2A_FSC_dips.csv")

#investigate sample gap1_all_3. Graph ridgeplots.
gap1_all_3 <- read.csv(file = "sc_distributions_gap1_all_3_all_timepoints.csv", stringsAsFactors = T) %>%
  mutate(generation = factor(generation, levels = unique(generation)))
gap1_all_3 %>%
  ggplot(aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
  xlab("Normalized fluorescence (a.u.)") +
  ylab("Generation") +
  ggtitle("GAP1 LTR + ARS 3") +
  theme_classic() +
  #scale_x_continuous(limits=c(0.0,5), breaks = c(0, 1, 5, 0.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 2.5))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
ggsave("gap1_all_3_ridgeplots.png")
