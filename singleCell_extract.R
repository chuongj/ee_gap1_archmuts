sc = read_csv("06_zero_ctrl_SingleCellDistributions_27_EE_GAP1_ArchMuts_2021.csv")
sc2 = read_csv("06_zero_ctrl_SingleCellDistributions_26_EE_GAP1_ArchMuts_2021.csv")
sc_distr_alltimepoints = rbind(sc, sc2)
df = sc_distr_alltimepoints %>%
  sc %>% group_by(sample, generation) %>%
  mutate(Med_FSC = median(`FSC-A`)) %>%
  distinct(Med_FSC, .keep_all = T) %>%
  select(-`FSC-A`, -`B2-A`) %>%
  write_csv(paste0(version_name,"_median_FSC_alltimepoints.csv"))

df %>% ggplot(mapping = aes(generation, Med_FSC, color = sample)) +
  geom_line()

### add to script
## calculate Med_B2A_FSC and Med_FSC and Med_B2A at the same time

geno = args[2] #GAP1 LTR KO
#  "1 copy control"       "GAP1 LTR + ARS KO"    "GAP1 ARS KO"          "GAP1 LTR KO"
#  "GAP1 WT architecture" "0 copy control"
sc_distr_alltimepoints %>%
sm %>%
  filter(Description == "GAP1 LTR KO") %>% #genotype
  group_by(sample, generation) %>%
  mutate(Med_FSC = median(`FSC-A`),
         Med_B2A = median(`B2-A`),
         Med_B2A_FSC = median(B2A_FSC),
         ) %>%
  distinct(Med_FSC, .keep_all = T) %>%
  distinct(Med_B2A, .keep_all = T) %>%
  distinct(Med_B2A_FSC, .keep_all = T) %>%
  select(-`FSC-A`, -`B2-A`, -B2A_FSC) %>% View()
#  write_csv(paste0(version_name,"_medians_alltimepoints.csv"))



### for each sample, subset it and write a sc_distributions_SampleName_allTimepoints.csv
sc = sc %>%
  filter(Description == "GAP1 LTR KO")
for(pop in unique(sc$sample)) {
 print(pop)
 sc %>%
 filter(sample == pop) %>%
 write_csv(paste0(version_name,"sc_distributions_",pop,"_all_timepoints.csv"))
}

###
pop_files = list.files(pattern = "sc_distributions_") %>% sort()
make_ridgeplots = function(file_name){
  pop_name = sub("sc_distributions_", "", sub("_all_timepoints.csv","", file_name))

  pop_data = read.csv(file_name, stringsAsFactors = T) %>%
    mutate(generation = factor(generation, levels = unique(generation)))

  lowcell = count(pop_data, generation, sample) %>% filter(n < 70000)

  pop_data %>%
    anti_join(lowcell) %>%
    ggplot(aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
    geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
    xlab("Normalized fluorescence (a.u.)") +
    ylab("Generation") +
    ggtitle(paste0(pop_name)) +
    theme_classic() +
    scale_x_continuous(limits=c(0.0,3), breaks = c(0, 1, 2, 3.0)) +
    scale_y_discrete(expand = expansion(add = c(0.2, 2.5))) + #expands the graph space or else the top is cut off
    scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
    theme(
      legend.position = 'none', #remove the legend
      axis.text.x = element_text(family="Arial", size = 25, color = "black"), #edit x-tick labels
      axis.text.y = element_text(family="Arial", size = 25, color = "black")
    )
  ggsave(paste0(pop_name,"_ridgeplot_scale2.pdf"))
}

