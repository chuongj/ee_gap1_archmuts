#####this script is for estimating breakpoints and then making dumbbell plots using read depth####
#G Avecilla (but code mostly stolen from S Lauer)
#Jan 2019
# Modified Julie Chuong 12-19-22

## mCitrine is located 514934 - 515650
## coordinates of GAP1 CDS for our samples are 518204 - 520012
## LTR11 coordinates are 513532:513846, LTR12 is ~ 520925

##Grace: For this set of samples, I used 3 SDs and >300 bp cutoff..
##However, for some samples, the CNV is clearly bigger than those bounds
##So I made a 2nd set of boundaries using relaxed criteria (>100bp) to refine the CNV
##

#Set working directory and samplesheet
setwd("/Users/juliechuong/Lab_docs/cloneseq1_temp")
samplesheet = "samplesheet_51_clones.csv"

#load libraries
library(tidyverse)
library(docstring)
library(ggalt)
library(gridExtra)

# set color pallete
ars_color = "#e26d5c"
#arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631","#f85c46", "#bb3521","#d9402a")
wt_color = "black"
#wtGrays = c("#354f52", "#666666", "#6b705c", "#414833" ,"#999999")
all_color = "#DEBD52"
#allGolds=c("#ffba08", "#faa307", "#dda15e", "#7f5539", "#9c6644", "#fdc409", "#9c7e1e", "#D9BB59")
ltr_color = "#6699cc"
#ltrBlues = c("#6699cc", "#005f73", "#0a9396", "#4292C6", "#2171B5", "#3799fb", "#66b3cc","#3a0ca3")

#### Function 1 #### 
#' Format read depth data. 
#' Subsets out chromosome 11, adds a column with mean read depth across chromosome 11, standard dev, and max RD
#'
#' @param a the read depth data list containing chromosomes, coordinate, read depths for every strain sequenced
#' @param name name of the item in a list, usually the strain/sample name
#'
#' @return list of tibbles containing only chromosome 11 data. Each row is a nucleotide position on chromosome 11. Columns are read depth, mean read depth, sd, and max
#' @export
#'
#' @examples
get_chrXI = function(x, name) {
  x=x %>% dplyr::filter(X1 == 'NC_001143.9')
  names(x) <- c("chromosome", "coordinate", "depth")
  x %>% mutate(mean = mean(depth), stdev = sd(depth), max = max(depth)) %>%
    mutate(sample = name) %>%
    group_by(sample) %>%
    mutate(linez = runmed(depth/mean,25001)) #for smooth lines later
}

#### Function 2 ### 
#' Get candidate CNV breakpoints on chromosome 11
#'
#' @param x
#' @param cnv_min_length minimum length of contiguous sequence you want to define as a CNV. 300 is liberal. 100 is even more relaxed.
#'
#' @return a tibble with start and end positions of CNVs breakpoints at
#' chromosome 11. The tibble gives CANDIDATE CNV breakpoints. Every row is a
#' single CNV of length cnv_min_length identified with our function. The start and end are the as a nucleotide positions of the upper and lower bounds of a CNV.
#' Manually plot x-intercept lines using the start and end values on a relative read depth plot.
#' Visually inspect and use other start end values to further define CNV breakpoints.
#'
#' @examples From Lauer et al. 2018 Methods: To manually estimate CNVs boundaries, we used a read depth–based approach.
#' For each sample sequenced, we used samtools to determine the read depth for each nucleotide in the genome. We liberally
#' defined CNVs by identifying greater or equal 300 base pairs of contiguous sequence when read depth was  greater or equal to 3 times the standard deviation
#' across Chromosome XI for GAP1 or Chromosome VIII for DUR3. These boundaries were further refined by visual inspection
#' of contiguous sequence greater or equal to 100 base pairs with read depth greater or equal to 3 times the standard deviation.
#' These analyses were only performed on sequenced clones because population samples are likely to have multiple CNVs
#' and breakpoints, thereby confounding read depth–based approaches. We compared manually estimated breakpoints to those
#' identified by the algorithms (S5 Table) and defined a set of “high-confidence breakpoints.”

get_breaks = function(x, cnv_min_length=300, SDtimes = 3) {
  hi_depth = x$depth >= sd(x$depth) * SDtimes #Boolean: is depth at each coordinate greater or equal to 3 sd's?
  test <- rle(hi_depth) #rle: runs in our case can be thought of as contiguous nucleotides that have read depth >= 3SD.
#' Contiguous DNA or runs that DON'T have same read depth are CNV breakpoints.
#' The output of rle() can be interpreted as "We saw the value Y for a run length of X. (equal to X consecutive times)
#' eg. We saw the value "FALSE" for run length of 196. FALSE in this case corresponds to when read depth is NOT equal to or
#' greater than 3SD of read depth.
  hi_depth_contig = which(test$values == TRUE & test$lengths >= cnv_min_length)
          #' groups runs of nucleotides into what I call contigs in which read depth is >= 3*SD.
          #' Each item is a contig (is a run).
          #' Each number reports the run length in bp, aka the length of the contig.
          #' Remember that values are T/F for contiguous.
          #' and that lengths are the number of nucleotide positions it was contigious for.S o the run.
  runs.lengths.cumsum = cumsum(test$lengths)
  ends = runs.lengths.cumsum[hi_depth_contig]
          #' the length of the run is simply the right end of the breakpoint.
          #' This commands maps the hi_depth_contig back to the chromosome position/index.
  #max(runs.lengths.cumsum) #sanity check - should be the length of the chromosome in bp
  newindex = ifelse(hi_depth_contig>1, hi_depth_contig-1, 0)
          #' Substract 1 from the length of contig, so that we can have the index of the previous contig positions, which is needed to determine start
          #' Doesn't make sense to want end position -1
  starts = runs.lengths.cumsum[newindex] + 1
          #'map it to previous contig's position on chromosome 11 and add back 1 bp run length
  if (0 %in% newindex) starts = c(1,starts)
  as_tibble(data.frame(starts, ends))
}

#### Function 3 ####
#' make_brkpt_plots
#' Makes chromosome 11 CNV breakpoints plot
#' @param name Strain name as a string to
#' @param left bp position of left CNV boundary from get_breaks()
#' @param right bp position of right CNV boundary from get_breaks()
#' @return outputs two plots. A relative read depth plot for inputted strain at Chromosome with red lines at CNV breakpoints. A second one zoomed in.

make_brkpt_plots = function(name, left, right){
  basePlot <-
    ggplot(allrd_chr11[sample(nrow(allrd_chr11), 5e3),], aes(x=coordinate,y=depth/mean)) + #downsample 5e4 points
    geom_point(size = .25) +
    theme_classic() +
    scale_color_manual(values = c("GAP1 CDS" = "dark green", "LTR" = ltr_color, "CNV breakpoint" = "red")) +
    #    geom_vline(xintercept = c(518204:520012), alpha = 0.01, color = 'dark green')+
    geom_vline(xintercept = mean(518204,520012), size = 3, alpha = 0.8, color = 'dark green')+
    #    geom_vline(xintercept = 520012, alpha = 0.8, color = 'dark green')+
    #  geom_vline(aes(xintercept = c(518204:520012), color='GAP1 CDS'), size = 1, alpha = 0.8,show.legend = TRUE) + #GAP1
    # geom_vline(aes(xintercept = 520012, color='GAP1 CDS'), size = 1, alpha = 0.8,show.legend = TRUE) + #GAP1
    #  geom_vline(aes(xintercept = 513532, color='LTR'), size = 1, alpha = 0.8, show.legend = TRUE) + #blue line to mark LTR
    #  geom_vline(aes(xintercept = 513846, color='LTR'), size=1, alpha = 0.8,show.legend = TRUE) + #blue line to mark LTR
    #    geom_vline(xintercept = c(513532:513846), alpha = 0.05, color = ltr_color)+
    geom_vline(xintercept = 513532, alpha = 0.8, size = 1, color = ltr_color)+
    #  geom_vline(aes(xintercept = 520594, color='LTR'), size=1, alpha = 0.8,show.legend = TRUE) + #blue line to mark LTR
    #geom_vline(aes(xintercept = 520925, color='LTR'), size=1, alpha = 0.8,show.legend = TRUE) + #blue line to mark LTR
    geom_vline(xintercept = 520925, alpha = 0.8, size = 1, color = ltr_color)+
    #    geom_vline(aes(xintercept = 520158, color = 'ARS'),size = 0.5, show.legend = TRUE)+
    #    geom_vline(aes(xintercept = 520406, color = 'ARS'),size = 0.5, show.legend = TRUE)+
    geom_vline(aes(xintercept = left, color = 'CNV breakpoint'), size = 1, alpha = 0.8, show.legend = TRUE) +
    geom_vline(aes(xintercept = right, color = 'CNV breakpoint'), size = 1, alpha = 0.8, show.legend = TRUE)+
    scale_y_continuous(expand = c(0, 0.5), 'Relative Read Depth', limits=c(-0.5, 6), breaks = c(0, 1, 2, 3, 4, 5, 6)) +
    scale_x_continuous(expand = c(0, 0),
                       limits=c(0,max(allrd_chr11$coordinate)),
                       labels = function(l){trans = l / 1000},
                       breaks = scales::pretty_breaks(n = 5), "Position on Chromosome XI (kb)") +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          plot.margin = margin(25, 25, 25, 25))+
    #    annotate(geom="text", x=-Inf,y=-Inf,hjust=-0.25,vjust=-1, label=paste0(left," & ",right), color="black")
    annotate("text", x = Inf, y = -0.5, label = paste0(left," & ",right), color="black") +
    coord_cartesian(clip = "off")
  
  genotype = allrd_chr11 %>% filter(sample == paste0(name)) %>% distinct(Description) %>% pull() %>% as.character() 
  
  pop = allrd_chr11 %>% filter(sample == paste0(name)) %>% distinct(pop_name) %>% pull() 
  
  generation = allrd_chr11 %>% filter(sample == paste0(name)) %>% distinct(generation) %>% pull() 
  
  #### whole chromosome plot (optional)
  # plot = basePlot %+% subset(allrd_chr11[sample(nrow(allrd_chr11), 5e4),], sample %in% c(name)) +
  #   ggtitle(paste0(name,"      ", genotype)) +
  #   annotate(geom="text",x=-Inf,y=-Inf,hjust=-6.5,vjust=-2, label=paste0(left," & ",right),color="black")
  
  #Save Whole Chromosome Plot
  #ggsave(plot = plot, filename = paste0(name,"_chr11_breakplot.pdf"), width = 8, height = 5, path = "./breakpoint_plots")
  
  # zoomed in plot at the left and right boundaries
  In = basePlot %+% subset(allrd_chr11[sample(nrow(allrd_chr11), 5e5),], sample %in% c(name)) +
    ggtitle(paste0("Clone ",name,"      ", genotype,"  ",pop,"    Generation ", generation)) +
    annotate(geom="text",x=-Inf,y=-Inf,hjust=-6.5,vjust=-2, label=paste0(left," & ",right),color="black") +
    scale_x_continuous(expand = c(0, 0),
                       limits=c(left-0.3e5,right+0.3e5),
                       labels = function(l){trans = l / 1000},
                       breaks = scales::pretty_breaks(n = 5), "Position on Chromosome XI (kb)")
  
  #output zoom-in plot only
  ifelse(!dir.exists("breakpoint_plots"), dir.create("breakpoint_plots"), FALSE)
  ggsave(plot = In, filename = paste0(name,"_zoomin_chr11_breakplot.pdf"), bg='#ffffff', width = 8, height = 5, path = "./breakpoint_plots")
  return(In)
#  return(basePlot)
  
  
}

######## RUN FUNCTIONS ##########


#### Get GAP1 copy numbers
copies = read_csv("cnv_cloneseq1_rd_summary.csv") %>% #outfile from rd_cnv_estimation.R script.
  rename(strain = sample) %>%
  mutate(gap1_rounded = round(gap1_copies)) %>%
  arrange(factor(Description, levels = c("GAP1 WT architecture","GAP1 LTR KO","GAP1 ARS KO","GAP1 LTR + ARS KO")), gap1_rounded) %>%
  relocate(gap1_rounded, .after = gap1_copies)

write_csv(copies, "cnv_cloneseq1_rd_summary_rounded.csv")

# Identify samples that have less than 2 copies of GAP1 
noCNV = copies %>% filter(gap1_rounded < 2) %>% select(strain) %>% pull()

#### STEP 1 ####
# Read in RD.txt files that was outputted from the `samtools depth -a ${bam} > ${NAME}_RD.txt`
# columns are chromosome number, coordinate, and read depth
files = list.files(pattern="*_RD.txt", path = './RD_files') #get vector of files
files = files[-which(str_detect(files, noCNV))] #remove noCNV samples
files = files[(-22:-66)] #exclude Titir's samples
files = files[22:47] #work with subset at a time because not enough RAM
data_list = lapply(paste0('./RD_files/', files), read_tsv, col_names = F, show_col_types = FALSE) #read in all _RD.txt files and puts them all in a list.
data_names = str_sub(files, 0, -8) #extract name string from files names
names(data_list) = data_names #assign name string to each items in our data list

#### STEP 2 ####
# Make list of read depth files for every strain resulting in
# chromosome 11 read depth data for every strain
rd_chr11 = map2(data_list, data_names, get_chrXI)

allrd_chr11 = do.call(rbind, rd_chr11) #row bind all items in our list one gigantic tibble

#merge with samplesheet so we have population and genotype metadata
metadata = read_csv(paste0("./",samplesheet), show_col_types = F)
allrd_chr11 = left_join(allrd_chr11, metadata, by = "sample") %>%
  relocate(c(sample, generation, pop_name, Description), .before = chromosome)

## STEP 3 ### 
#get all the candidate starts and ends
brkpts_300 = map(rd_chr11, get_breaks)
brkpts_300_4 = map(rd_chr11, get_breaks, SDtimes=4) #higher depth cutoff with 4 standard deviations, more stringent

# havent run these yet but might later
# brkpts_100 = map(rd_chr11, get_breaks, cnv_min_length=100)
# brkpts_500 = map(rd_chr11, get_breaks, cnv_min_length=500)
# brkpts_1000 = map(rd_chr11, get_breaks, cnv_min_length=1000)
# brkpts_5000 = map(rd_chr11, get_breaks, cnv_min_length = 5000)


### STEP 4 - FACET PLOT Relative Read Depth for ALL CLONES ######

### Add a black smooth/trend line
fPlot <- allrd_chr11[sample(nrow(allrd_chr11),5e4),] %>% #downsample to 5e6
  ggplot(aes(x=coordinate,y=depth/mean, color = Description)) +
  geom_point(size = .25) +
  geom_line(aes(x= coordinate, y=linez), color = "dimgray")+ #custom smooth line in dark gray
  ggtitle("Read Depth of Chromosome XI for 22 Isolated Clones")+
  theme_classic() +
  scale_color_manual(values = c("#DEBD52","#6699cc","#e26d5c", "gray"), #custom colors  # third, first, second
                     limits = c("GAP1 LTR + ARS KO","GAP1 LTR KO","GAP1 ARS KO", "GAP1 WT architecture"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("LTR and ARS KO", "LTR KO", "ARS KO", "Wildtype architecture") #third, now you can change legend labels
                     )+
  scale_y_continuous(expand = c(0, 0.5), 'Relative Read Depth', limits=c(-0.5, max(allrd_chr11$max/allrd_chr11$mean)), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)) +
  scale_x_continuous(expand = c(0, 5),
                     limits=c(0, 6.7e5),
                     labels = function(l) {trans = l / 1000},
                     breaks = scales::pretty_breaks(n = 5), "Position on Chromosome XI (kb)") +
  facet_wrap(~strain, ncol = 6) +
  theme(#strip.background = element_blank(),strip.text.x = element_blank(), # remove facet titles for the Poster
        axis.text.y = element_text(size = 14, color = "black"), #make y axis numbers bigger
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16,angle=90)
        )+
  guides(color = guide_legend(override.aes = list(size=5),
                              title = "Genotype")
         )
fPlot

#save plots
#ggsave(fPlot,"facetplots_4.pdf", bg='#ffffff', path = "./breakpoint_plots")
#ggsave(fPlot,"facetplots_5.png", bg='#ffffff', path = "./breakpoint_plots")
#quartz.save("./breakpoint_plots/facetplots_5.png", type = "png", device = dev.cur(), bg='#ffffff')
#quartz.save("./breakpoint_plots/facetplots_5.pdf", type = "pdf", device = dev.cur(), bg='#ffffff')

##### STEP 5 # Plot breakpoints and visually identify and refine breakpoints. ######

# Manually choose the left and right values outputted from get_breaks() aka brkpts_300 and plot it

# add breakpoints coordinates to dataframe

startend = data.frame(sample = character(), start=numeric() , end =numeric())

brkpts_300$`2836`[1,1]
brkpts_300$`2836`[37, 2]
make_brkpt_plots("2836",504098,652500)
startend[1,] = c("2836",504098,652500)

brkpts_300$`2837` %>% View()
make_brkpt_plots("2837",511958,535586)
startend[2,] = c("2837", 511958,535586)

brkpts_300$`2838` %>% View()
make_brkpt_plots("2838",511966,635322)
startend[3,] = c("2838", 511966,635322)

brkpts_300$`2839` %>% View()
make_brkpt_plots("2839",513658,552077)
startend[4,]=c("2839",513658,552077)

brkpts_300$`2840` %>% View()
make_brkpt_plots("2840",508699,527376)
startend[5,]=c("2840",508699,527376)

brkpts_300$`2841`%>% View()
make_brkpt_plots("2841",511000,586033)
startend[6,]=c("2841",511000,586033)

brkpts_300$`2842`%>% View()
make_brkpt_plots("2842",511453,586170)
startend[7,]=c("2842",511453,586170)

brkpts_300$`2843`%>% View()
make_brkpt_plots("2843",511972,538360)
startend[8,]=c("2843",511972,538360)

brkpts_300$`2844`%>% View()
brkpts_300_4$`2844`%>% View()
make_brkpt_plots("2844",513104,544644) 
startend[9,]=c("2844",513104,544644) 

brkpts_300$`2845`%>% View()
make_brkpt_plots("2845",443000,567931) 
startend[10,]=c("2845",443000,567931) 

brkpts_300$`2846`%>% View()
make_brkpt_plots("2846",512861,539737) 
startend[11,]=c("2846",512861,539737)

brkpts_300$`2847`%>% View()
make_brkpt_plots("2847",512457,521751)
startend[12,]=c("2847",512457,521751)

brkpts_300$`2848`%>% View()
brkpts_300_4$`2848`%>% View()
make_brkpt_plots("2848",491720,532607)
startend[13,]=c("2848",491720,532607)

brkpts_300$`2849`%>% View() #aneuploid, gain of chr11, see genome_rd_plot
startend[14,] = c("2849",1,670191)

brkpts_300$`2850`%>% View() #aneuploid, gain of chr11,see genome_rd_plot
startend[15,]=c("2850",1,670191)

brkpts_300$`2851`%>% View() #aneuploid, gain of chr11, see genome_rd_plot
startend[16,]=c("2851",1,670191)

brkpts_300$`2852`%>% View()
brkpts_300_4$`2852` %>% View()
make_brkpt_plots("2852",448362,631285)
startend[17,]=c("2852",448362,631285)

brkpts_300$`2853`%>% View()
make_brkpt_plots("2853",445401,569905)
startend[18,]=c("2853",445401,569905)

brkpts_300$`2854`%>% View()
make_brkpt_plots("2854",448401,521845) 
startend[19,]=c("2854",448401,521845) 

brkpts_300$`2855`%>% View()
make_brkpt_plots("2855",445398,569997) 
startend[20,]=c("2855",445398,569997) 

brkpts_300$`2856`%>% View()
make_brkpt_plots("2856",447070,569641)
startend[21,]=c("2856",447070,569641)

brkpts_300$`2923`%>% View()
make_brkpt_plots("2923",512093,574510) 
startend[22,]=c("2923",512093,574510) 

brkpts_300$`2924`%>% View()
make_brkpt_plots("2924",512392,577434) 
startend[23,]=c("2924",512392,577434) 

brkpts_300$`2925`%>% View() #aneuploidy, see genome plots
startend[24,]=c("2925",1,670191)

brkpts_300$`2926`%>% View() #aneuploidy, see genome plots
startend[25,]=c("2926",1,670191)

brkpts_300$`2927`%>% View() #aneuploidy, see genome plots
startend[26,]=c("2927",1,670191)

brkpts_300$`2928`%>% View()
make_brkpt_plots("2928",512395,544558)
startend[27,]=c("2928",512395,544558)

brkpts_300$`2929`%>% View()
make_brkpt_plots("2929",513053,521595)#LTR-bound
startend[28,]=c("2929",513053,521595)

brkpts_300$`2930`%>% View()
make_brkpt_plots("2930",511546,521576) 
startend[29,]=c("2930",511546,521576) 

brkpts_300$`2931`%>% View()
make_brkpt_plots("2931",508569,527471) 
startend[30,]=c("2931",508569,527471) 

brkpts_300$`2932`%>% View()
make_brkpt_plots("2932",505192,521926) 
startend[31,]=c("2932",505192,521926) 

brkpts_300$`2933`%>% View()
make_brkpt_plots("2933",508093,521497)
startend[32,]=c("2933",508093,521497)

brkpts_300$`2934`%>% View()
make_brkpt_plots("2934",510943,563017)
startend[33,]=c("2934",510943,563017)

brkpts_300$`2935`%>% View()
make_brkpt_plots("2935",456147,522246)
startend[34,]=c("2935",456147,522246)

brkpts_300$`2936`%>% View()
make_brkpt_plots("2936",460480,522246)
startend[35,]=c("2936",460480,522246)

brkpts_300$`2937`%>% View()
make_brkpt_plots("2937",459785,521993)
startend[36,]=c("2937",459785,521993)

brkpts_300$`2938`%>% View()
make_brkpt_plots("2938",489885,669931) 
startend[37,]=c("2938",489885,669931)

brkpts_300$`2939`%>% View()
make_brkpt_plots("2939",489864,669931) 
startend[38,]=c("2939",489864,669931) 

brkpts_300$`2940`%>% View()
make_brkpt_plots("2940",489864,669931)
startend[39,]=c("2940",489864,669931)

brkpts_300$`2941`%>% View()
make_brkpt_plots("2941",489873,669931) 
startend[40,]=c("2941",489873,669931) 

brkpts_300$`2942`%>% View()
make_brkpt_plots("2942",489899,522086) 
startend[41,]=c("2942",489899,522086) 

brkpts_300$`2943`%>% View()
make_brkpt_plots("2943",458141,522246) 
startend[42,]=c("2943",458141,522246) 

brkpts_300$`2944`%>% View()
make_brkpt_plots("2944",490211,522351)
startend[43,]=c("2944",490211,522351)

brkpts_300$`2945`%>% View()
make_brkpt_plots("2945",490211,522351) 
startend[44,]=c("2945",490211,522351) 

brkpts_300$`2946`%>% View()
make_brkpt_plots("2946", 457020,669931)
startend[45,]=c("2946", 457020,669931)

brkpts_300$`2947`%>% View()
make_brkpt_plots("2947",446040,	572227)
startend[46,]=c("2947",446040,	572227)

brkpts_300$`2948`%>% View()
make_brkpt_plots("2948",453977,525289)
startend[47,]=c("2948",453977,525289)

#### STEP 6 - Output a table with breakpoints with estimated copy numbers ####
startend_c = startend %>%
  mutate(cnv_length = as.numeric(end)-as.numeric(start),
         start = as.numeric(start), 
         end = as.numeric(end))

# merge `startend` with `copies` 
breaks = copies %>% rename(sample = strain) %>% 
  right_join(startend_c) %>% write_csv("cloneseq1_rd_Breakpoints_Copies.csv")

breaks = read_csv("cloneseq1_rd_Breakpoints_Copies.csv")
####### STEP 7 - Make Dumbbell plots
# Dumbbell plots show the genomic position and span of CNV with a horizontal line, like a handle of dumbbell weight and two dots at either end are the CNV breakpoints, like a ends of dumbbell weight, like this: ()=========()

#devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

##### example from githib ggalt #####
#library(hrbrthemes)
#df <- data.frame(trt=LETTERS[1:5], l=c(20, 40, 10, 30, 50), r=c(70, 50, 30, 60, 80))
#ggplot(df, aes(y=trt, x=l, xend=r)) +
#  geom_dumbbell(size=3, color="#e3e2e1",
#                colour_x = "#5b8124", colour_xend = "#bad744",
#                dot_guide=TRUE, dot_guide_size=0.25) +
#  labs(x=NULL, y=NULL, title="ggplot2 geom_dumbbell with dot guide") +
#  theme_ipsum_rc(grid="X") +
#  theme(panel.grid.major.x=element_line(size=0.05))
###

#### Dumbbell Plot - Grouped by Genotype  ####
db_geno=
  breaks %>%
  filter(gap1_rounded>1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=3, dot_guide=TRUE, dot_guide_size=0.25)+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670000),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text=element_text(size=14, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
  #facet_grid(~gap1_rounded)

ggsave(filename = "cloneseq1_dumbbell_plot_122022_aneu.png", width = 8, height = 10, bg = "white")
ggsave(filename = "cloneseq1_dumbbell_plot_122022_aneu.pdf", width = 8, height = 10, bg = "white") 

# Make dumbbell plot version 2 
# exclude aneuplodies clones and zoom into 400-660kb region
db_noAneu = 
  breaks %>%
  filter(gap1_rounded>1, start != 1) %>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  ggplot(aes(x = start, xend=end, y = clone, color = Description))+
  geom_dumbbell(size=3, dot_guide=TRUE, dot_guide_size=0.25)+
  scale_color_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                     limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  # scale_x_discrete(expand=c(0,1),
  #                  labels = function(l) {trans = l / 1000},
  #                  "Position on Chromosome XI (kb)")+ 
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670000),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_wrap(~generation) #facet by generation? 
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text=element_text(size=14, color = "black"),
        legend.box.margin=margin(20,20,20,20), #move legend away from plot
        legend.position = "none",
        plot.margin = margin(25, 25, 25, 25))+
  guides(color = guide_legend(title = "Genotype"))
db_noAneu
ggsave(filename = "cloneseq1_dumbbell_plot_122022_noAneu.png", width = 8, height = 10, bg = "white")
ggsave(filename = "cloneseq1_dumbbell_plot_122022_noAneu.pdf", width = 8, height = 10, bg = "white") 

### Barplot grouped by genotype, then sort by copy number low to high ####
bar = breaks%>%
  filter(gap1_rounded>1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% #reorder genotype custom order, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorders the strain order after arrange()
  ggplot(aes(clone, gap1_rounded, fill = Description)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c(wt_color,ltr_color,ars_color,all_color), #custom colors
                    limits = c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                    labels = c("Wildtype architecture","LTR removed", "ARS removed","LTR and ARS removed") #third, now you can change legend labels
  )+
  coord_flip()+
  theme_classic() +
  ylab("GAP1 copy number") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.text.x = element_text(size = 12, color = "black"),
              #axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.text=element_text(size=12, color = "black"),
        legend.box.margin=margin(10,10,10,10), #move legend away from plot
        plot.margin = margin(25, 25, 25, 25)
        )+
  guides(fill = guide_legend(title = "Genotype")) #change legend title

bar
### combine dumbbell plot and barplot side-by-side
library(gridExtra)
both = grid.arrange(db_geno, bar, ncol = 2)
ggsave("dumbbell_and_barplot_122022.png", plot = both, width = 14, height = 9)

#add table to side the graph 
copy_tbl<-tableGrob(breaks%>%
                    filter(gap1_rounded>1)%>%
                    arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), desc(gap1_rounded), desc(cnv_length)) %>% select(gap1_rounded), 
                    theme_classic()
  )
quartz()
grid.arrange(db_geno,copy_tbl, ncol=2)

### STEP 7 - 
# Make dumbbell plot version 3, color by population 

#db_pop=
  breaks %>%
  filter(gap1_rounded>1, start !=1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), 
          factor(pop_name, levels = c("gap1_2","gap1_ltr_3","gap1_ltr_2","gap1_ars_5","gap1_all_3","gap1_all_2")),
          desc(gap1_rounded), 
          desc(cnv_length)) %>% #custom reorder by genotype, then population, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
   # mutate(pop_name = factor(pop_name)) %>% 
  ggplot(aes(x = start, xend=end, y = clone, color = pop_name))+
  geom_dumbbell(size=3, dot_guide=TRUE, dot_guide_size=0.25)+
  scale_color_manual(values = c("gray", "blue", "darkblue", "pink", "yellow", "orange"), #custom colors
                     limits = c("gap1_2","gap1_ltr_2","gap1_ltr_3","gap1_ars_5","gap1_all_2","gap1_all_3"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture, Population 2","LTR removed, Population 2", "LTR removed, Population 3", "ARS removed, Population 5","LTR and ARS removed, Population 2", "LTR and ARS removed, Population 3") #third, now you can change legend labels
  )+
  scale_x_discrete(expand=c(0,1),
                   labels = function(l) {trans = l / 1000},
                   "Position on Chromosome XI (kb)")+
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670000),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  #facet_grid(generation ~ ., scale = "free")+
 # facet_grid(gap1_rounded ~ ., space = "fixed")+
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.text=element_text(size=12, color = "black"),
        #legend.position = "none"
        )+
  guides(color = guide_legend(title = "Population"))
  
ggsave(filename = "cloneseq1_dumbbell_plot_POP_Gen.png", width = 8, height = 10, bg = "white")  

ggsave(filename = "cloneseq1_dumbbell_plot_POP.png", width = 8, height = 10, bg = "white")  


#######################
# Import cloneseq0 aka Dec2021 g79 clones 
g79 = read.csv("../WGS_clones_g79/WGS_g79_clones_rd_estimated_breakpoints_GeneAnnotations.csv") %>% rename(sample = strain, pop_name = population) %>%
  select(-start_genes.....1kb., -end_genes......1kb.) 

g79 = g79 %>% mutate(Description = recode(Description, 'GAP1 WT architecture' = "Wildtype architecture", 'GAP1 LTR KO' = "LTR KO", 'GAP1 ARS KO' = "ARS KO", 'GAP1 LTR + ARS KO' = "LTR and ARS KO"),
         sample = as.numeric(sub("DGY","",g79$sample))
         )

clonedata = full_join(breaks, g79)
#write_csv(clonedata, "cloneseq0and1_Breaks_Copies.csv")

plot = clonedata %>%
  filter(gap1_rounded>1, start !=1)%>%
  arrange(factor(Description, levels = rev(c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))), 
          factor(pop_name, levels = c("gap1_2","gap1_1",
                                      "gap1_ltr_3","gap1_ltr_2","gap1_ltr_1", "gap1_ars_5","gap1_ars_1","gap1_all_3","gap1_all_2", "gap1_all_1")),
          desc(gap1_rounded), 
          desc(cnv_length)) %>% #custom reorder by genotype, then population, and then arrange copy number low to high
  mutate(clone = factor(sample, levels = unique(sample))) %>% #reorder strain order %>%
  # mutate(pop_name = factor(pop_name)) %>% 
  ggplot(aes(x = start, xend=end, y = clone, color = pop_name))+
  geom_dumbbell(size=3, dot_guide=TRUE, dot_guide_size=0.25)+
  scale_color_manual(values = c(wtGrays[4:5],ltrBlues[1:3], arsSalmons[1:2], allGolds[1:3]), #custom colors
                     limits = c("gap1_1","gap1_2","gap1_ltr_1","gap1_ltr_2","gap1_ltr_3","gap1_ars_1","gap1_ars_5","gap1_all_1","gap1_all_2","gap1_all_3"), #second, change order of legend items, by listing in the order you want em. using the real names in the aes(color =  ) argument
                     labels = c("Wildtype architecture, Population 1","Wildtype architecture, Population 2","LTR removed, Population 1","LTR removed, Population 2","LTR removed, Population 3", "ARS removed, Population 1","ARS removed, Population 5","LTR and ARS removed, Population 1","LTR and ARS removed, Population 2", "LTR and ARS removed, Population 3") #third, now you can change legend labels
  )+
  scale_x_discrete(expand=c(0,1),
                   labels = function(l) {trans = l / 1000},
                   "Position on Chromosome XI (kb)")+
  scale_x_continuous(expand=c(0,1),#limits=c(400000,670000),
                     labels = function(l) {trans = l / 1000},"Position on Chromosome XI (kb)")+ #breaks = scales::pretty_breaks(n=4)
  #scale_y_discrete(limits=rev)+
  labs(x = "Position on Chromosome XI", y= "Clone") +
  # facet_grid(gap1_rounded ~ ., space = "fixed")+
  theme_bw()+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        strip.text.x = element_text(size = 12), #facet font size
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.text=element_text(size=12, color = "black"),
        legend.position = "none"
  )+
  guides(color = guide_legend(title = "Genotype, Population"))

plot

ggsave("cloneseq0and1_dumbbell_plot_POP_fullWidth.png", bg = "white", width = 8, height = 10)
ggsave("cloneseq0and1_dumbbell_plot_POP_fullWidth.pdf", bg = "white", width = 8, height = 10)

#Facet dumbbell plot by generation 
plot + ggforce::facet_col(vars(generation), scales = 'free', space = 'free')
#plot + ggforce::facet_col(~generation, scales = 'free', space = 'free')
ggsave("cloneseq0and1_dumbbell_POP_fullWidth_FacetGen.png", bg = "white", width = 8, height = 10)
ggsave("cloneseq0and1_dumbbell_POP_fullWidth_FacetGen.pdf", bg = "white", width = 8, height = 10)
#YAY!! we did it :) great job boo 


