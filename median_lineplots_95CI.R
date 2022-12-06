# 12-5-22
# Scrap the idea of using a geom_smooth or any sort of spline for now since
# generalized additive model presented during lab meeting was not a good fit. Did not fit the raw data.
# current parameters... span = 1, formula = 'y ~ s(x, bs = "cs")'

# set up freq_and_counts data frame
freq_and_counts = read_csv("freq_and_counts_Merged_080622_all_timepoints.csv")

fails = freq_and_counts %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(str_detect(Description, "control")) %>%
  select(Description, Strain, generation, Gate, Frequency, name, Count) %>%
  mutate(flag = case_when(Strain == "DGY1" & Gate == "zero_copy" & Frequency >= 95 ~ "pass",
                          Strain == "DGY1" & Gate == "zero_copy" & Frequency < 95 ~ "fail",
                          Strain == "DGY1" & Gate == "one_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY1" & Gate == "two_or_more_copy" & Frequency >=11 ~ "fail",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY500" & Gate == "zero_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY500" & Gate == "two_or_more_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY1315" & Gate == "zero_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY1315" & Gate == "one_copy" & Frequency >= 11 ~ "fail"
  ))%>%
  dplyr::filter(flag == "fail") %>%
  arrange(Description)
View(fails)

weird_early = freq_and_counts %>%
  filter(generation < 30,
         Type %in% c("Experimental", "1_copy_ctrl"),
         Description %in% c("1 copy control", "GAP1 WT architecture","GAP1 LTR KO"),
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample) %>%
  #select(-name, -`Outflow well`, -Media)
  filter(Frequency > 15)

#chose these timepoints by eye
weird_tp = freq_and_counts %>%
  filter(sample == "gap1_4" & Gate == "two_or_more_copy" & generation == 66 |
           sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 166|
           sample == "gap1_all_5" & Gate == "two_or_more_copy" & generation == 116|
           sample == "gap1_all_6" & Gate == "two_or_more_copy" & generation == 124|
           sample == "gap1_ltr_2"
  )

# Plot 1
# all 5 wildtype lines drawn
# lineplot of proportion of CNVs over time



# Plot 2
# median of all population for each genotype lineplot  with 95% confidence intervals



###########
# Current GAM is not a good fit with its parameter span = 1, formula = 'y ~ s(x, bs = "cs")'
# as discussed when I presented the GAM during lab meeting
# Instead try to just plot the median of the populations with 95% confidence intervals bars


# median line plot code -
# https://www.datanovia.com/en/lessons/ggplot-error-bars/

freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental",
         generation <= 203) %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  ggplot(aes(generation, Frequency, color = Description)) +
  geom_line(aes(linetype = supp, group = supp))+
  geom_point()+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = supp),
    width = 0.2
  )
ggplot(sampleDataFrame, aes(x=percentage, y=temp, colour=userId)) +
  geom_line() +
  geom_line(aes(y=means), size=2, color="black")
