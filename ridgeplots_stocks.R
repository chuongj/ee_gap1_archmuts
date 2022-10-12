# Ridgeplots of Stock Strains

#Draw histograms of the stock strains to see if CNVs appeared at gen8 or not
stock_folder_path = "/Volumes/GoogleDrive/My Drive/greshamlab/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/mCitrine CNV reporter test"
all_files = list.files(paste0(stock_folder_path,"/FCS_files"))
#files_pos <- match(sample_sheet$name, all_files) #find positions of files that match the sample_sheet names
#all_files[files_pos] #subset of files to read in
stocks_gs <- cyto_setup(path = paste0(stock_folder_path,"/FCS_files"), restrict = T, details = F) #must specify a folder, all fcs files will be read

sample_sheet = read_csv(file = paste0(stock_folder_path,"/sample_sheet.csv"))
sample_sheet = sample_sheet %>%
  mutate(name = paste0(name, ".fcs"))

#annotate experiment details using sample sheet
for(i in 1:length(names(sample_sheet))){
  flowWorkspace::pData(stocks_gs)[names(sample_sheet[i])]<-sample_sheet[i]
}
cyto_details(stocks_gs) %>% View()

#transform stock strain gating set
GFP_trans <- cyto_transformer_logicle(stocks_gs,
                                      channels = c("B2-A", "B3-A"),
                                      widthBasis = -10
)#returns it as a list
FSC_SSC_trans <- cyto_transformer_log(stocks_gs,
                                      channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")
)
combined_trans <- cyto_transformer_combine(GFP_trans,FSC_SSC_trans)
transformed_timepoint_gating_set <- cyto_transform(stocks_gs,
                                                   trans = combined_trans)

#apply gating template to stock strain gating set
cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")

#get raw single cell flow data of stock strains
timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = T, channels = c("FSC-A", "B2-A", "B3-A")) #raw flow data of each single cell as a list of matrices
sc_distributions <- map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
  mutate(name = as.factor(name), #convert `name` to factor
         B2A_FSC = `B2-A`/`FSC-A`,
         B3A_FSC = `B3-A`/`FSC-A`) %>% #compute normalized GFP over forward scatter
  left_join(sample_sheet) #join by name column to add other metadata to raw data


##Plot ridgeplots of stock strains
#B2-A channel ridgeplot of stock strains
sc_distributions %>%
  mutate(Strain = fct_reorder(Strain, desc(B2A_FSC))) %>% #reorder samples by y values low to high
  ggplot(aes(x = B2A_FSC, y = Strain, fill = `mCitrine copy number`)) +
  geom_density_ridges(scale=1.5, quantile_lines = F, quantiles = 2) +
  xlab("B2-A mCitrine fluorescence/forward scatter") +
  ylab("Strain") +
  ggtitle("Stock Strains Ridgeplots") +
  theme_minimal() +
  scale_x_continuous(limits=c(0,3))+
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) +
  scale_fill_discrete(breaks=c("zero copy",
                               "one copy",
                               "two copy"
  )) + #change order of legend items
  theme(axis.text.x = element_text(family="Arial", size = 10, color = "black"),
        axis.text.y = element_text(family="Arial", size = 10, color = "black"))
ggsave("stock_strains_B2A_ridgeplot_noMedians.png")

#B3-A channel ridgeplot of stock strains
sc_distributions %>%
  mutate(Strain = fct_reorder(Strain, desc(B3A_FSC))) %>% #reorder samples by y values low to high
  ggplot(aes(x = B3A_FSC, y = Strain, fill = `mCitrine copy number`)) +
  geom_density_ridges(scale=1.5, quantile_lines = F, quantiles = 2) +
  xlab("B3-A mCitrine fluorescence/forward scatter") +
  ylab("Strain") +
  ggtitle("Stock Strains Ridgeplots") +
  theme_minimal() +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) +
  scale_fill_discrete(breaks=c("zero copy","one copy","two copy")) + #change order of legend items
  theme(axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 10, color = "black"))
ggsave("stock_strains_ridgeplot_B3A_noMedians.png")
