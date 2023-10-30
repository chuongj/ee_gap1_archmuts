# Forked from gating_112122.R
# making new CNV lineplots because the 4 wildtype lines are under-reporting GAP1 CNVs due to broken CNV reporter subpop
# Make the order of lines showing up: Wildtype, ARS removed, LTRs removed, LTR and ARS removed. 

# Load required packages
library(tidyverse)
library(plotly)

setwd("/Users/juliechuong/Library/CloudStorage/GoogleDrive-jc10007@nyu.edu/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

#load frequency data
freq_and_counts = read_csv("freq_and_counts_merged_CLEAN_121622.csv")

#### 10/26/23 add the alpha-fade to lineplots because people always ask what the variance is between populations ####
# https://plotly.com/r/line-charts/ 
# Plotly resource Pieter showed me! Also called filled line charts 
med_freq_counts = freq_and_counts %>%
  mutate(proportion = Frequency/100) %>%
  dplyr::filter(generation <= 203) %>%
  dplyr::filter(!(sample == "gap1_1"  |
                    sample == "gap1_2" |
                    sample == "gap1_4" |
                    sample == "gap1_5")
  ) %>%
  mutate(Description = factor(Description, levels=c("GAP1 WT architecture", "GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")))%>%
  group_by(generation, Description) %>%
  mutate(med = median(Frequency))

#subset by genotype for easier graphing
wt = med_freq_counts %>% dplyr::filter(Description == "GAP1 WT architecture")
arsko = med_freq_counts %>% dplyr::filter(Description == "GAP1 ARS KO")
ltrko = med_freq_counts %>% dplyr::filter(Description == "GAP1 LTR KO")
allko = med_freq_counts %>% dplyr::filter(Description == "GAP1 LTR + ARS KO")

# get lows and highs (min and max per timepoint) per genotype
low = ltrko %>% group_by(generation) %>% summarize(mins = min(Frequency)/100)
high = ltrko %>% group_by(generation) %>% summarize(maxs = max(Frequency)/100)
medians = ltrko %>% group_by(generation) %>% summarise(median = max(med)/100)
lowars = arsko %>% group_by(generation) %>% summarize(mins = min(Frequency)/100)
highars = arsko %>% group_by(generation) %>% summarize(maxs = max(Frequency)/100)
medars = arsko %>% group_by(generation) %>% summarise(median = max(med)/100)
alllow = allko %>% group_by(generation) %>% summarise(mins = min(Frequency)/100)
allhigh = allko %>% group_by(generation) %>% summarise(maxs = max(Frequency)/100)
allMedian = allko %>% group_by(generation) %>% summarise(median = max(med)/100)

#graph single line for Wildtype architecture (n = 1)
fig = plot_ly(x=wt$generation, y=wt$Frequency/100, type = 'scatter', mode = 'lines', line = list(color = 'black', width = 10), showlegend = FALSE, name = "Wildtype architecture") %>% layout(paper_bgcolor='white', plot_bgcolor='white',
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
# save_image(fig, file = "propCNV_WT_fade_102723_1200x800.png", scale = 1, width = 1200, height = 800)
# save_image(fig, file = "propCNV_WT_fade_102723_1200x800.pdf", scale = 1, width = 1200, height = 800)

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
save_image(fig, file = "propCNV_WT+ARS+LTR_fade_102723_1200x800.png", scale = 1, width = 1200, height = 800)
save_image(fig, file = "propCNV_WT+ARS+LTR_fade_102723_1200x800.pdf", scale = 1, width = 1200, height = 800)


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
# save_image(fig, file = "propCNV_full_fade_102723_1200x800.png", scale = 1, width = 1200, height = 800)
# save_image(fig, file = "propCNV_full_fade_102723_1200x800.pdf", scale = 1, width = 1200, height = 800)
 