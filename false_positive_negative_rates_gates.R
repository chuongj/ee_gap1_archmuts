# False Positive and False Negative Rates of CNVs
# example code that I might not really use in the analysis
# it's more like a sanity check

#determine upper threshold for zero copy gate, aka False Negative Rate of Detecting One Copy
##Everyone should be start as One Copy (atleast) except the Zero Copy Control
# one copy control falling into the zero copy or two+ copy gates
freq %>%
  filter(Count>70000) %>%
  #filter(str_detect(Description, "control"), Gate == "zero_copy") %>%
  filter(Description == "1 copy control" & Gate == "zero_copy" | Description == "1 copy control" & Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>% #View()
  #ggplot(aes(Frequency)) +
  #geom_histogram(bins = 50)
  #  geom_boxplot()
  summarize(IQR(Frequency)+median(Frequency)) #One Copy FN Rate is  6.86% = (median+IQR)

freq %>%
  filter(Count>70000) %>%
  filter(Description == "1 copy control", Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  #ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(median(Frequency) + 1.5*IQR(Frequency)) #8.58
  summarize(median(Frequency) + IQR(Frequency))

#Determine CNV False Negative Rate. 2 copy control falling into 1copy or 0copy Gates: median + IQR
freq %>%
  filter(Count>70000) %>%
  #filter(str_detect(Description, "control"), Gate == "zero_copy") %>%
  filter(Description == "2 copy control", Gate != "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  #ggplot(aes(Frequency)) + geom_boxplot()
  summarize(IQR(Frequency))
# CNV FN Rate is 1.33% = median + IQR

#Determine lower threshold for zero copy control gate
freq %>%
  filter(Count>70000) %>%
  filter(Description == "0 copy control", Gate == "zero_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>% #View()
  #ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(min(Frequency))
  summarize(median(Frequency)-1.5*IQR(Frequency))

#Determine lower threshold for one copy control gate
freq %>%
  filter(Count>70000) %>%
  filter(Description == "1 copy control", Gate == "one_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>% #View()
  #ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(min(Frequency))
  summarize(median(Frequency)-1.5*IQR(Frequency)) #87.4, the left end of the lower whisker
#summarize(1.5*IQR(Frequency))

#Determine lower threshold for two copy control gate
freq %>%
  filter(Count>70000) %>%
  filter(Description == "2 copy control", Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  ggplot(aes(Frequency)) + geom_boxplot()
#summarize(min(Frequency))
summarize(median(Frequency)-1.5*IQR(Frequency)) #92.3, the left end of the lower whisker
#summarize(1.5*IQR(Frequency))
