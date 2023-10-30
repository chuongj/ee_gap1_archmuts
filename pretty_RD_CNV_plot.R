# title: "Read depth based estimation of gene copy numbers"
# author: G. Avecilla. Modified Julie Chuong. 
# Date: 10/04/22

# Purpose: Estimate GAP1 copy numbers from read depth of WGS data of isolated clones from the EE_GAP1_ArchMuts_2021.

#Upstream of this, I aligned reads by running the mapcheck.sh bash/slurm script
#which outputted read depth files, one for each sequenced strain `{NAME}_RD.txt`. 
# Column headers are not there but are chromosome number, coordinate, and read depth.

# This R script will take the output from samtools depth as the input
#`samtools depth -as ${bam} > ${ID}_RD.txt`

#On the hpc, I transferred the {name}_RD.txt files to my laptop using scp, secure copy
#scp jc10007@greene.hpc.nyu.edu:/scratch/cgsb/gresham/julie/WGS_clones/cloneseq1/try2/d_*/*_RD.txt ./RD_files

#Other inputs: samplesheet.csv about the clones that were sequenced. 
## START HERE for R analysis

#Set up working directory and load packages
setwd("~/Lab_docs/cloneseq1_temp")
library(tidyverse)
#library(docstring) #like docstring in python for R

#Note, the ref sequence used for alignment and mapping was
#/scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna which contains the GFP KanMX GAP1 CNV reporter (JC manually searched for the sequences to verify)   

# Function to read in read depth file
get_rd_data = function(x) {
#' This function reads in a file, names columns, removes the mitochondrial chromosome, and changes chromosome names from RefSeq chromosomes to numeric
#' @param x {name}_RD.txt file outputted from samtools
#' @return a tibble with three columns - chromosome number, coordinate, and read depth
#' @example DGY2071_read_data = get_rd_data("DGY2017_RD.txt")
  print(paste0("reading in", x))
  rd <- read_tsv(x, col_names = FALSE, show_col_types = FALSE)
  names(rd) <- c("chromo", "coordinate", "depth")
  rd <- rd[rd$chromo !="NC_001224.1",] #mitochondrial
  #change chromosome names
  rd <- rd %>% mutate(chromosome = case_when(chromo == 'NC_001133.9' ~ 1,
                                             chromo == 'NC_001134.8' ~ 2,
                                             chromo == 'NC_001135.5' ~ 3,
                                             chromo == 'NC_001136.10' ~ 4,
                                             chromo == 'NC_001137.3' ~ 5,
                                             chromo == 'NC_001138.5' ~ 6,
                                             chromo == 'NC_001139.9' ~ 7,
                                             chromo == 'NC_001140.6' ~ 8,
                                             chromo == 'NC_001141.2' ~ 9,
                                             chromo == 'NC_001142.9' ~ 10,
                                             chromo == 'NC_001143.9' ~ 11,
                                             chromo == 'NC_001144.5' ~ 12,
                                             chromo == 'NC_001145.3' ~ 13,
                                             chromo == 'NC_001146.8' ~ 14,
                                             chromo == 'NC_001147.6' ~ 15,
                                             chromo == 'NC_001148.4' ~ 16)) %>%
    select(chromosome, coordinate, depth)
  return(rd)
}



# The next set of functions calculate the relative copy number at GAP1 and other specified genes.  
# Genes are specified by giving the coordinates from the reference sequence.
#' The next set of functions calculate the relative copy number at GAP1 and other specified genes
#' Genes are specified by giving the coordinates from the reference sequence
#' @param x the tibble returned from get_rd_data()
#'
#' @return numeric; copy number at specified gene relative to mean whole genome coverage
#' @export
#'
#' @examples get_gap1_reporter_copies(DGY2071_read_data)
get_gap1_reporter_copies = function(x) {
  #coordinates are adjusted to take into account the PrACT1-GFP-KanMX
  mean((x$depth[x$chromosome=="11" & x$coordinate>=513946 & x$coordinate<=520246])/mean(x$depth))
}
#this function gets relative GFP copy number
get_GFP_copies = function(x){
  #coordinates of GFP (717bp) located upstream of GAP1
  mean((x$depth[x$chromosome=="11" & x$coordinate>=514934 & x$coordinate<=515650])/mean(x$depth))

}

#this function gets the relative copy number at GAP1
get_gap1_copies = function(x) {
  #coordinates are GAP1 CDS location our reference sequence GCF_000146045.2_R64_genomic_GAP1.fna
  mean((x$depth[x$chromosome=="11" & x$coordinate>=518438 & x$coordinate<=520246])/mean(x$depth))
}

#This function makes plots of read depth normalized to chromosome 11.
#This lets us see GAP1 CNVs. 
make_rd_plot_chr11 = function(df, output, dir) {
#' this function makes plots read depth normalized to chromosome 11 read depth
#'
#' @param df is a dataframe of read data,
#' @param output name of outputted plot file
#' @param dir is the directory in which to save it
#'
#' @return save a .tiff plot read depth normalized to chromosome 11
#'
#' @examples
  #check to see if the directory exists and if it doesn't, make it
  if(dir.exists(dir) == FALSE) {
    dir.create(dir)
  }
  #get chromosome 11 only
  chr11 <- df[which(df$chromosome == 11),]
  avg11 <- mean(chr11$depth)
  linez <- runmed(chr11$depth/avg11,25001)
  newdata <- cbind(chr11, linez)
  newdata <- data.frame(newdata)
  #downsample for ease of plotting
  chr11plot <- newdata[sample(nrow(newdata), 20000),]
  #plot
  plot <- ggplot(chr11plot, aes(x=chr11plot[,2],y=chr11plot[,3]/avg11)) +
    geom_point(size = .25) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), 'Relative Read Depth', limits=c(-0.5, 6)) +
    # scale_x_continuous(expand = c(0, 0), limits=c(5000,655000), breaks = scales::pretty_breaks(n = 5), "Position on Chromosome 11 (kb)") +
    scale_x_continuous(expand = c(0, 0),
                       limits=c(0,655000),
                       labels = function(l){trans = l / 1000},
                       breaks = scales::pretty_breaks(n = 5), "Position on Chromosome XI (kb)") +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          plot.margin = margin(25, 25, 25, 25))+
    
    
    geom_line(aes(x= chr11plot[,2], y=chr11plot[,4]), color = "red")
  ggsave(filename = paste0(output, "_rd_chr11_plot.png"), plot = plot, path = dir)
}

####### RUN FUNCTIONS #######

# Get RD file for 2924 cloneseq1 only to make pretty graph
setwd("~/Lab_docs/cloneseq1_temp")
files = list.files(path = './RD_files', pattern = '2924_RD.txt')
files

#perform get_rd_data() 
clones_rd_data = get_rd_data(paste0('./RD_files/',files))

#make plot
output = "2924_RD_pretty"
my_dir = "~/Lab_docs/cloneseq1_temp"
make_rd_plot_chr11(df = clones_rd_data, output, my_dir)