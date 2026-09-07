#originally created on: October 18th, 2017
#by: Nick Sard

#ABOUT: This script was written to evaluate GT data associated with a VCF

#loading necessary libraries
library(tidyverse)
library(gridExtra)
library(grid)

#setting working directory and reading in data
#setwd("C:/Users/sard/Google Drive/R/Data analysis/2017/Sea Lamprey/Round 1 VCF/")
#setwd("G:/My Drive/R/Data analysis/2017/Sea Lamprey/Round 1 VCF")
setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/post_RAD analysis")
#listing files
list.files("Input/")

######
### note you need to delete the # before CHROM, and the number of lines to be skipped before the header my change ###
######
#reading in file
df <- read.table(file = "Input/batch_0_ns_mod.vcf",header = T,sep = "\t",strip.white = T,skip = 9)
head(df[,1:10])

#getting nuclotide position of SNPs
df2 <- df  %>%
  mutate(ID1 = ID) %>% 
  separate(col = "ID1",into = c("misc","nt_pos"),sep = "_") %>%
  select(ID,nt_pos)
head(df2)

#now a graphic
p1 <- ggplot(df2, aes(x=as.numeric(nt_pos))) +
  geom_bar()+
  labs(x="Nucleotide position on the read",y="Number of loci")+
  theme_bw()
p1

#making into long form and extracting information in each cell
geno <- df %>%
  select(-FORMAT,-FILTER,-QUAL) %>%
  gather(key = "id",value = "all",-CHROM:-INFO) %>%
  mutate(loc = substr(x = id,start = 3,stop = 5)) %>%
  separate(col = all,into = c("gt","dp","allele_dp","lh"),sep = ":") %>%
  separate(col = allele_dp, into = c("allele_a","allele_b"),sep = ",") %>%
  separate(col = lh,into = c("lh1","lh2","lh3"),sep = ",")
head(geno)

#need to change the locations
unique(geno$loc)
length(unique(geno$loc))

#change genotype call format to something more amenable to subsetting and that can be used as column names below
geno$gt <- gsub(pattern = "\\./\\.",replacement = NA,x = geno$gt)
geno$gt <- gsub(pattern = "0/0",replacement = "AA",x = geno$gt)
geno$gt <- gsub(pattern = "0/1",replacement = "AB",x = geno$gt)
geno$gt <- gsub(pattern = "1/1",replacement = "BB",x = geno$gt)
head(geno)

#replacing . in NA for allele depths, so that I can make them numeric below
geno$allele_a <- gsub(pattern = "\\.",replacement = NA,x = geno$allele_a)
geno$allele_b <- gsub(pattern = "\\.",replacement = NA,x = geno$allele_b)
head(geno)

#identifying which individuals are heterozygous
geno <- geno %>% mutate(het = ifelse(test = gt == "AB",1,0))

#identifying individuals that were gentoyped
geno$test <- ifelse(is.na(geno$gt),"No","Yes")
head(geno)

#now for a graphic to visual how complete genotype data are for each location
p2 <- ggplot(geno,aes(x=id,y=ID,fill=gt))+
  facet_wrap(~loc,nrow = 1,scales = "free_x")+
  geom_tile()+
  labs(x= "Samples",y="Loci",fill = "Genotyped")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size=8),
        legend.position = "bottom")
p2

#creating a table for each location, summarizing the number of individuals
#that were or were not genotyped at more than %50 of the loci
loc_count <- geno %>%
  group_by(loc,id,test) %>%
  count() %>%
  spread(key = test,value = n,fill = 0) %>%
  mutate(total = No + Yes) %>%
  mutate(prop = No/total) %>%
  mutate(fail = ifelse(prop >= .95,1,0)) %>%
  group_by(loc) %>%
  summarize(count = n(),
            n_failp95 = sum(fail))
loc_count
p2.1 <- tableGrob(loc_count,rows = NULL)

#saving them to file
rem_ind <- geno %>%
  group_by(loc,id,test) %>%
  count() %>%
  spread(key = test,value = n,fill = 0) %>%
  mutate(total = No + Yes) %>%
  mutate(prop = round(No/total,2)) %>%
  mutate(fail = ifelse(prop >= .95,1,0)) %>%
  filter(fail == 1) %>%
  select(-fail)
rem_ind
p2.2 <- tableGrob(rem_ind,rows = NULL)
#write.table(x = rem_ind,file = "Output/check_these_individuals.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#distribution of median coverage per locus
#making things numeric first
head(geno)
geno$dp <- as.numeric(geno$dp)
geno$allele_a <- as.numeric(geno$allele_a)
geno$allele_b <- as.numeric(geno$allele_b)

#calculating median coverage per locus
geno1 <- geno %>%
  filter(!is.na(gt)) %>%
  group_by(ID) %>%
  summarize(med_cov = median(dp))
head(geno1)

#calculating the center of that distribution
med_med_coverage <- median(geno1$med_cov)

#figure 
p3 <- ggplot(geno1, aes(x=med_cov)) +
  geom_histogram()+
  annotate(geom = "text",label = med_med_coverage,x = med_med_coverage,y = 50,color="white")+
  labs(x="Median coverage per locus",y="Number of loci")+
  theme_bw()
p3

#now looking at loci genotype to calculate global major and minor allele frequencies and observed heterozygoties
head(geno)
geno2 <- geno %>%
  filter(!is.na(gt)) %>%
  select(ID,gt) %>%
  group_by(ID,gt) %>%
  count() %>%
  spread(key = gt,value = n,fill = 0) %>%
  mutate(total_gt = AA + AB + BB) %>%
  mutate(Hobs = round(AB/ total_gt, 3)) %>%
  mutate(p = round( ((2*AA) + AB)/ (2*total_gt),3)) %>%
  mutate(q = round(1-p,3)) %>%
  mutate(ID2 = ID) %>%
  separate(col = ID2,into = c("misc","read_pos"),sep = "_") %>%
  select(-misc)
geno2$read_pos <- as.numeric(geno2$read_pos)
head(geno2)
mean(geno2$Hobs)
sd(geno2$Hobs)
mean(geno2$q)
sd(geno2$q)
max(geno2$read_pos)
min(geno2$read_pos)
#plotting minor allele freq vs read pos
p3.1 <- ggplot(geno2, aes(x=read_pos,y=q))+
  geom_point(alpha=.5,size=3)+
  geom_smooth()+
  theme_bw()+
  labs(x="Nucleotide position on the read", y= "Minor allele frequency")
p3.1
p3.2 <- ggplot(geno2, aes(x=read_pos,y=Hobs))+
  geom_point(alpha=.5,size=3)+
  geom_smooth()+
  theme_bw()+
  labs(x="Nucleotide position on the read", y= "Observed heterozygosity")
p3.2

#looking to count unusual loci
geno2$check <- "good"
geno2$check[geno2$AA != 0 & geno2$AB == 0 & geno2$BB ==0] <- "AA_only"
geno2$check[geno2$AA == 0 & geno2$AB != 0 & geno2$BB ==0] <- "AB_only"
geno2$check[geno2$AA == 0 & geno2$AB == 0 & geno2$BB !=0] <- "BB_only"
geno2$check[geno2$AA != 0 & geno2$AB == 0 & geno2$BB !=0] <- "No_hets"

#getting a summary table
g2_sum <- geno2 %>%
  group_by(check) %>%
  count()
g2_sum
p4 <- tableGrob(g2_sum,rows = NULL)

#saving them to file
no_good <- geno2 %>%
  filter(check != "good") %>%
  arrange(check) %>%
  as.data.frame()
no_good
write.table(x = no_good,file = "Output/check_these_loci.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#getting distrbution of observed heterozygosity per locus

#calculating the center of that distribution
med_med_coverage <- round(median(geno2 %>% filter(check == "good") %>% pull(Hobs)),2)
head(geno2)
p5 <- ggplot(geno2 %>% filter(check == "good"), aes(x=Hobs))+
  geom_histogram()+
  annotate(geom = "text",label = med_med_coverage,x = med_med_coverage,y = 50,color="white")+
  labs(x="Observed heterozygosity",y="Number of loci")+
  theme_bw()
p5

#looking at the number of populations where a locus is genotyped in more than half the individuals
head(geno)
p6 <- geno %>%
  mutate(gt2 = ifelse(test = is.na(gt),0,1)) %>%
  group_by(loc,ID) %>%
  summarise(counts = n(),
            gt3 = sum(gt2),
            prop  = round(gt3/counts,2)) %>%
  mutate(test  = ifelse(prop >= .5,1,0)) %>%
  select(ID,loc,test) %>%
  group_by(ID) %>%
  summarize(total = sum(test)) %>%
  ggplot(aes(x=total)) +
  geom_histogram(binwidth = 1)+
  theme_bw()+
  labs(x="Number of populations where\n>=50% of ind. are genotyped",y="Number of loci")
p6

#looking at the number of SNPs genotyped at >=50% within a population
p7 <- geno %>%
  mutate(gt2 = ifelse(test = is.na(gt),0,1)) %>%
  group_by(loc,ID) %>%
  summarise(counts = n(),
            gt3 = sum(gt2),
            prop  = round(gt3/counts,2))%>%
  mutate(test  = ifelse(prop >= .5,1,0)) %>%
  group_by(loc) %>%
  summarise(nloc = sum(test)) %>% 
  ggplot(aes(x=loc,y=nloc)) +
  geom_bar(stat="identity",color="black",fill="grey")+
  theme_bw()+
  labs(x="Populations",y="Number of loci where\n>=50% of ind. are genotyped")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1))
p7

#now creating a plot to identify putative paralogs using McKinney et al. (2017)
#method, with Ho and D cut-offs choosen based on their simulated data
head(geno)

geno3 <- geno %>%
  filter(gt == "AB") %>%
  select(ID,INFO,id:het) %>%
  group_by(ID) %>%
  summarize(count = n(),
            a1 = sum(allele_a),
            a2 = sum(allele_b),
            tot_reads = sum(a1,a2),
            sig = sqrt(tot_reads*.5*.5),
            tr_2 = tot_reads/2,
            z = (tr_2 - a1)/sig)
head(geno3)

#now merging geno3 to with geno2 to make graph
head(geno2)
head(geno3)
g2 <- merge(x = geno2,y = geno3,by = "ID")
head(g2)

#the row numbers differed between gsum and gsum2, I'm guessing gsum had some loci with no hets
geno2 %>% filter(Hobs == 0) %>% nrow()
nrow(geno2)-nrow(geno3)
#yup

#based on McKinney et al's simulated data
#singleton loci are those that have
#H <= 0.60 & Z <= abs(7)

#mkaing a column to putative singletons and non-singletons
head(g2)
g2$singleton <- ifelse(test = (g2$Hobs <= 0.6 & abs(g2$z) <=  7),yes = "Yes",no = "No")
table(g2$singleton)

#now the graphic
p8 <- ggplot(g2, aes(x=Hobs,y=z,color=singleton)) +
  geom_point(size=3)+
  theme_bw()+
  labs(x="Proportion of heterozygotes", y="Read-ratio deviation",color="Putative\nsingleton\nlocus")
p8

#saving g2 to file for an analysis of Fis on singletons and not
write.table(x = g2,file = "Output/locus_type.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#saving them all to PDFs
pdf("Output/VCF_mussles_plots.pdf", onefile = TRUE)
grid.arrange(p1,p3.1,p3.2,p4,nrow=2)
grid.arrange(p3,p5,p6,p7,nrow=2)
grid.arrange(p8)
#grid.arrange(p2.1,p2.2,ncol=2)
dev.off()

pdf("Output/VCF_mussels_gt_graphic.pdf", onefile = TRUE)
grid.arrange(p2)
dev.off()
#fin!