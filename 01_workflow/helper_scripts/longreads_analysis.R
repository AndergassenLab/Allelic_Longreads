#######################################
#   Long read allelic ratio analysis (& short read corr)
#   Author: Lison Lemoine
#   Created: 27.04.2023
#   Notes:  
#######################################
library(gtools)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(biomaRt)
library(ggrepel)
library(VennDiagram)

# var -----------------------------------------------------------------
path_to_folder <- "~/longread_ASE"
name <- "brain"

input <- paste0(path_to_folder, "/output/")
output <- =paste0(path_to_folder, "/plots/")

short_reads <- paste0(path_to_folder, "/input/locus_table.txt")
refseqID <- paste0(path_to_folder, "/input/Refseq_IDs_unique_max_locus.txt")

long_read1=paste0(input, name, "_h1.count")
long_read2=paste0(input, name, "_h2.count")

pat <- "#283D90" # "#40509E"
mat <- "#BE222F" # "#CE3730"
bae <- "#cccccc" 

# read LR -------------------------------------------------------
count_mat<-read.delim(long_read1, header=F)
count_pat<-read.delim(long_read2, header=F)
colnames(count_mat)<-c("name", "reads_h1")
colnames(count_pat)<-c("name", "reads_h2")

long_counts<-full_join(count_mat, count_pat, by=join_by(name))
long_counts<-long_counts[long_counts$name!="__ambiguous"&long_counts$name!="__no_feature",]
rm(count_mat, count_pat)

# calculate ratio -----------------------------------------------------
long_counts$ratio_lr <- round(apply(long_counts[,c("reads_h1","reads_h2")], 1,
                                    function(x) {return(x["reads_h1"]/(x["reads_h1"]+x["reads_h2"]))}),2)
long_counts<-na.omit(long_counts) #this removes 0 reads
long_counts$lr_total <- apply(long_counts[,c("reads_h1","reads_h2")], 1,
                              function(x) {return((x["reads_h1"]+x["reads_h2"]))})

# filter lr -----------------------------------------------------------
long_counts <- long_counts[long_counts$lr_total>=10,]

gen_names <- read.delim(refseqID, header=T)

long_counts <- left_join(long_counts, gen_names, by=join_by(name))%>%
dplyr::select(-c(strand, start, end))
remove(gen_names)

# short reads ---------------------------------------------------------
short_counts<-read.delim(short_reads)
colnames(short_counts)<-c("chr", "start", "end", "name", "A1_reads", 
                          "A2_reads", "total_reads", "allelic_score", 
                          "allelic_ratio")
short_counts<-na.omit(short_counts[short_counts$total_reads>=20,])

# merge ---------------------------------------------------------------
lrsr<-left_join(short_counts,long_counts, by=join_by(name), copy=T) %>%
mutate(chr=chr.x)%>%
dplyr::select(-c(chr.x, chr.y))
lrsr<-na.omit(lrsr)

# autosome mega table -------------------------------------------------
autosomes<-rbind(long_counts%>%
                      filter(chr!="chrX")%>%
                      mutate(read="long", total=lr_total, ratio=ratio_lr)%>%
                      dplyr::select(c(name, chr, total, ratio, read)),
                    short_counts%>%
                      mutate(read="short", total=total_reads, ratio=allelic_ratio)%>%
                      filter(chr!="chrX"&chr!="chrY")%>%
                      dplyr::select(c(name, chr, total, ratio, read)) 
)

# sex chr filter ------------------------------------------------------
lrsr<-lrsr%>%filter(chr!="chrX")%>%filter(chr!="chrY")
rm(long_counts, short_counts)

# allelic score -------------------------------------------------------
calc.allelic.score <- function(x) {
  if (all(x==0) || any(is.na(x)))
    return(0)
  if (sum(x) < 10)
    return(0)
  allelic.score <- log10(pbinom(min(matrix(x,ncol=2)), x[1]+x[2], 0.5, 
                                lower.tail=TRUE, log = FALSE)) 
  ifelse(x[1] > x[2], -allelic.score,allelic.score) 
}

lrsr$score_lr <- round(apply(lrsr[,c("reads_h1","reads_h2")], 1,
                             calc.allelic.score),3)

# Loci: Replacing +/-infinitive with the highest/lowest finite value
max_value <- range(abs(lrsr$score_lr),finite=TRUE)[2]
lrsr$score_lr[lrsr$score_lr == "Inf"] <- max_value
lrsr$score_lr[lrsr$score_lr == "-Inf"] <- -max_value

lrsr<-na.omit(lrsr)

rm(max_value, calc.allelic.score)

# table clean up ------------------------------------------------------
## thinner table -----
lrsr<-subset(lrsr, select = -c(start, end, A1_reads, A2_reads))

## reorganize -----
lrsr<-lrsr%>%relocate(name)%>%relocate(chr)%>%relocate(lr_total, .after=reads_h2)%>%
  relocate(score_lr, .after=ratio_lr)%>%
  relocate(total_reads, .after=score_lr)%>%
  relocate(allelic_ratio, .after=total_reads)%>%
  relocate(allelic_score, .after=allelic_ratio)

colnames(lrsr)<-c("chr", "name", "reads_h1", "reads_h2", "t_reads_lr", "ratio_lr", "score_lr", "t_reads_sr", "ratio_sr", "score_sr")


