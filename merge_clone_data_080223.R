setwd("/Users/juliechuong/Lab_docs")

dat = read_csv("/Users/juliechuong/Lab_docs/workingcopy_080123_isolated_clones_seq_database_EE_GAP1_ArchMuts_2021 .csv")

seq0and1 = read_csv("/Users/juliechuong/Lab_docs/cloneseq1_temp/cloneseq0and1_Breaks_Copies.csv")
seq2 = read_csv("/Users/juliechuong/Lab_docs/cloneseq2_temp/cloneseq2_rd_Breakpoints_Copies.csv")
seq2_other = read_csv("~/Lab_docs/cloneseq2_temp/cnv_cloneseq2_rd_summary_rounded.csv")


colnames(seq2)
colnames(seq0and1)
colnames(dat)

#change colname of in dat from DGY Strain Number to strain so we can merge via this common column
dat =rename(dat, sample = `DGY Strain Number`)
dat = rename(dat, generation = Generation)
dat = dat %>% dplyr::filter(`Pool ID` %in% c("cloneseq0", "cloneseq1", "cloneseq_2"))
seq2_other = rename(seq2_other, sample = strain)

nrow(dat)
nrow(seq2)  #160
nrow(seq0and1) #69
nrow(seq2) + nrow(seq0and1) == nrow(dat) #FALSE

dat %>% dplyr::filter(`Pool ID`=="cloneseq_2") %>% nrow() #168
dat %>% dplyr::filter(`Pool ID`=="cloneseq0") %>% nrow() #22
dat %>% dplyr::filter(`Pool ID`=="cloneseq1") %>% nrow() #48

dat %>% dplyr::filter(!(`Pool ID`=="cloneseq_2"| 
                          `Pool ID`=="cloneseq0"|
                          `Pool ID`=="cloneseq1"
                          )) %>% View()

dat_sub = dat %>% select(sample, Timepoint, `Ancestor Strain`, `Pool ID`, `Sequencing Run(s)`, `Platform - Read Length - Cycle`, `Filepath to SR Seq Data`, `MultiQC Report`)

seq2and2 = left_join(seq2_other, seq2)
seq012 = rbind(seq0and1, seq2and2)
dat_0_1_2 = left_join(dat_sub, seq012)
dat_0_1_2 %>% write_csv(file = "clones_data_080123.csv") 

# 2835 has one copy GAP1 


