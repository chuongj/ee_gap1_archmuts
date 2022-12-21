library(tidyverse)
library(chromoMap)
setwd("~/Lab_docs/ee_gap1_archmuts")

# From SGD, I downloaded ChrXI features

# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html

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
        

#### Make a chromosome file as needed for chromoMap ####
my_feats %>% View()
chrXI = data.frame(name = "XI",
                   start = 1,
                   end = 666816,
                   cent = 440129
                   )
colnames(chrXI)<-NULL

####
#unique(my_feats[,5]) 

chromoMap(list(chrXI),
          list(my_feats),
          chr_color = c("#ede0d4"),
          #chr_width = 15,
          #chr_length = 10,
          #labels = T,
          y_chr_scale = 17, #bring ruler closer to chromosome 
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("#FFC20A", "#e26d5c", "#5887FF", "#ede0d4","#999933", "black"))) # tRNA, ars, ltr, centromere, GAP1 ORF, telomere
 #unique(my_feats[,5]) 
# FFC20A = tRNA, yellow
#"#999933" = GAP1 ORF olive
#55C1FF
#ltr blue #6699cc

### zoom in! 
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
          data_colors = list(c("#FFC20A", "#5887FF","#999933","#e26d5c")))

1.#FFC20A"
2. "#5887FF"
3. "#e26d5c",
4 "#e26d5c",
