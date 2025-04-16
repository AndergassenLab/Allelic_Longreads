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
path_to_folder="/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/ge92cem2/longread_ASE"
name="brain"

input=paste0(path_to_folder, "/output/")
output=paste0(path_to_folder, "/plots/")

# I included these in the input folder
short_reads=paste0(path_to_folder, "/input/locus_table.txt")
refseqID=paste0(path_to_folder, "/input/Refseq_IDs_unique_max_locus.txt")

long_read1=paste0(input, name, "_h1.count")
long_read2=paste0(input, name, "_h2.count")

pat <- "#283D90" # "#40509E"
mat <- "#BE222F" # "#CE3730"
bae <- "#cccccc" 

# read LR -------------------------------------------------------
br_count_1<-read.delim(long_read1, header=F)
br_count_2<-read.delim(long_read2, header=F)
colnames(br_count_1)<-c("name", "reads_h1")
colnames(br_count_2)<-c("name", "reads_h2")

br_count<-full_join(br_count_1, br_count_2, by=join_by(name))
br_count<-br_count[br_count$name!="__ambiguous"&br_count$name!="__no_feature",]
rm(br_count_1, br_count_2)

# calculate ratio -----------------------------------------------------
br_count$ratio_lr <- round(apply(br_count[,c("reads_h1","reads_h2")],1,function(x) {return(x["reads_h1"]/(x["reads_h1"]+x["reads_h2"]))}),2)
br_count<-na.omit(br_count) #this removes 0 reads
br_count$lr_total <- apply(br_count[,c("reads_h1","reads_h2")],1,function(x) {return((x["reads_h1"]+x["reads_h2"]))})

# filter lr -----------------------------------------------------------
br_count<-br_count[br_count$lr_total>=10,]

gen_names<-read.delim(refseqID, header=T)

br_count<-left_join(br_count, gen_names, by=join_by(name))%>%dplyr::select(-c(strand, start, end))
remove(gen_names)

# short reads ---------------------------------------------------------
br_tb<-read.delim(short_reads)
colnames(br_tb)<-c("chr", "start", "end", "name", "A1_reads", "A2_reads", "total_reads", "allelic_score", "allelic_ratio")
br_tb<-na.omit(br_tb[br_tb$total_reads>=20,])

# merge ---------------------------------------------------------------
br_lrsr<-left_join(br_tb,br_count, by=join_by(name), copy=T)%>%mutate(chr=chr.x)%>%dplyr::select(-c(chr.x, chr.y))
br_lrsr<-na.omit(br_lrsr)

# autosome mega table -------------------------------------------------
br_autosomes<-rbind(br_count%>%
                      filter(chr!="chrX")%>%
                      mutate(read="long", total=lr_total, ratio=ratio_lr)%>%
                      dplyr::select(c(name, chr, total, ratio, read)),
                    br_tb%>%
                      mutate(read="short", total=total_reads, ratio=allelic_ratio)%>%
                      filter(chr!="chrX"&chr!="chrY")%>%
                      dplyr::select(c(name, chr, total, ratio, read)) 
)

# sex chr filter ------------------------------------------------------
br_lrsr<-br_lrsr%>%filter(chr!="chrX")%>%filter(chr!="chrY")
rm(br_count, br_tb)

# allelic score -------------------------------------------------------
calc.allelic.score <- function(x) {
  if (all(x==0) || any(is.na(x)))
    return(0)
  if (sum(x) < 10)
    return(0)
  allelic.score <- log10(pbinom(min(matrix(x,ncol=2)), x[1]+x[2], 0.5, lower.tail=TRUE, log = FALSE)) 
  ifelse(x[1] > x[2], -allelic.score,allelic.score) 
}

br_lrsr$score_lr <- round(apply(br_lrsr[,c("reads_h1","reads_h2")],1,calc.allelic.score),3)

# Loci: Replacing +/-infinitive with the highest/lowest finite value
max_value <- range(abs(br_lrsr$score_lr),finite=TRUE)[2]
br_lrsr$score_lr[br_lrsr$score_lr == "Inf"] <- max_value
br_lrsr$score_lr[br_lrsr$score_lr == "-Inf"] <- -max_value

br_lrsr<-na.omit(br_lrsr)

rm(max_value, calc.allelic.score)

# table clean up ------------------------------------------------------
#thinner table
br_lrsr<-subset(br_lrsr, select = -c(start, end, A1_reads, A2_reads))

#reorganize
br_lrsr<-br_lrsr%>%relocate(name)%>%relocate(chr)%>%relocate(lr_total, .after=reads_h2)%>%
  relocate(score_lr, .after=ratio_lr)%>%
  relocate(total_reads, .after=score_lr)%>%
  relocate(allelic_ratio, .after=total_reads)%>%
  relocate(allelic_score, .after=allelic_ratio)

colnames(br_lrsr)<-c("chr", "name", "reads_h1", "reads_h2", "t_reads_lr", "ratio_lr", "score_lr", "t_reads_sr", "ratio_sr", "score_sr")

#normalize & filter
br_lrsr_norm<-mutate(br_lrsr, ratio_sr_norm=abs(ratio_sr-0.5), ratio_lr_norm=abs(ratio_lr-0.5))#%>%filter(t_reads_lr>=20)


# imprinted genes -----------------------------------------------------

imp_genes<-data.frame(name=c("Grb10", "Peg3", "Meg3", "Zdbf2", "Airn",
                             "Impact", "Rian", "Mirg", "Ago2", "Peg13",  
                             "Zrsr1", "Ndn", "A330076H08Rik", "Igf2", "Kcnk9", 
                             "Kcnq1ot1", "Snurf", "Snrpn", "Peg10", "Sgce", 
                             "B230209E15Rik", "Snhg14", "Zim1", "Nap1l5", "Usp29", 
                             "Plagl1", "Rasgrf1", "Ipw", "A230057D06Rik"), 
                      type=c("PAT", "PAT", "MAT", "PAT", "PAT",
                             "PAT", "MAT", "MAT", "MAT", "PAT",
                             "PAT", "PAT", "PAT", "MAT", "MAT", 
                             "PAT", "PAT", "PAT", "PAT", "PAT", 
                             "PAT", "PAT", "MAT", "PAT", "PAT", 
                             "PAT", "PAT", "PAT", "PAT"))

br_lrsr_norm<-left_join(br_lrsr_norm, imp_genes, by=join_by(name))
br_lrsr_norm<-na.replace(br_lrsr_norm, "BAE")

# plots ---------------------------------------------------------------
make_graph_0_1<-function(df, titl, method){
  exact=NULL
  if (method=="spearman"){
    exact=F
  }
  corr.norm<-cor.test(df$ratio_lr, df$ratio_sr, method = method, exact=exact)
  
  p<-ggplot(df, aes(x=ratio_sr, y=ratio_lr, color=type, shape=type)) +
    geom_point(position=position_jitter(w=0.005))+
    scale_shape_manual(values=c(20, 24, 25))+
    scale_color_manual(values=c(bae, mat, pat))+
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip = 'off')+
    ggtitle(titl)+
    annotate("text", x=0.5, y=0.05, label=paste0("R = ", round(corr.norm$estimate, 3)), fontface =2)+
    annotate("text", x=0.3, y=df[df$name%in%c("Gnas", "Ube2g2", "Malat1"),]$ratio_lr, label=df[df$name%in%c("Gnas", "Ube2g2", "Malat1"),]$name, fontface =2)+
    annotate("segment", x=0.3,
             xend=df[df$name%in%c("Gnas", "Ube2g2", "Malat1"),]$ratio_sr,
             y=df[df$name%in%c("Gnas", "Ube2g2", "Malat1"),]$ratio_sr, yend = df[df$name%in%c("Gnas", "Ube2g2", "Malat1"),]$ratio_lr)+
    theme_classic()+
    theme(aspect.ratio=1)
  return(p)
}
make_graph_imp<-function(df, titl, method){
  exact=NULL
  if (method=="spearman"){
    exact=F
  }
  
  corr.norm<-cor.test(df$ratio_lr, df$ratio_sr, method = method, exact=exact)
  
  p<-ggplot(df, aes(x=ratio_sr, y=ratio_lr, label=name, color=type, shape=type)) +
    geom_point(position=position_jitter(h=0.005,w=0.005))+
    scale_shape_manual(values=c(24, 25))+
    scale_color_manual(values=c(mat, pat))+
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip = 'off')+
    annotate("text", x=0.45, y=0, label=paste0("R = ", abs(round(corr.norm$estimate, 3))), fontface =2)+
    ggtitle(titl)+
    geom_text_repel(
      force_pull   = 0, # do not pull toward data points
      nudge_y      = -0.1,
      direction    = "x",
      hjust        = 3.5,
      segment.size = 0.1,
      size = 2.5) +
    theme_classic()+
    theme(aspect.ratio=1)
}

br_imp<-filter(br_lrsr_norm, name%in%imp_genes$name)

p_all<-make_graph_0_1(br_lrsr_norm, paste0("BA_BxC all, n=", nrow(br_lrsr_norm)), "pearson")

p_imp<-make_graph_imp(br_imp, paste0("BA_BxC imprinted, n=", nrow(br_imp)), "pearson")

#pdf(file=paste0(output,"scatter_imp.pdf"), width=4, height=3)
#grid.arrange(p_all, p_imp, nrow=1)
#dev.off()  

pdf(file=paste0(output,"scatter_all.pdf"), width=4, height=3)
p_all
dev.off()  
pdf(file=paste0(output,"scatter_imp.pdf"), width=4, height=3)
p_imp
dev.off() 

# pies ----------------------------------------------------------------
ratio.type<- function(gen.sum){
  gen.sum$type<-"BAE"
  
  for (i in 1:nrow(gen.sum)) {
    if (gen.sum[i,]$ratio>=0.7){
      gen.sum[i,]$type<-"MAT"
    }
    else if (gen.sum[i,]$ratio<=0.3){
      gen.sum[i,]$type<-"PAT"
    }
  }
  return(gen.sum)
}
br_autosomes<-ratio.type(br_autosomes)

#write_csv(br_autosomes, "./Documents/table.csv")

br_over<-br_autosomes%>%filter(name%in%br_lrsr$name)%>%pivot_wider(names_from = read, values_from = c(total, ratio, type))

stat.all<-data.frame(count=c(nrow(br_over%>%filter(type_long!="BAE")),
                             nrow(br_over%>%filter(type_long=="BAE")),
                             nrow(br_over%>%filter(type_short!="BAE")),
                             nrow(br_over%>%filter(type_short=="BAE"))), 
                     ratio_type=rep(c("ASE", "BAE"), 2), 
                     read=rep(c("Long read", "Short read"), each=2))


pdf(file=paste0(output,"pie_shortread_ratios.pdf"), width=4, height=3)

ggplot(subset(stat.all, read=="Short read"), aes(x="", y=count, fill=ratio_type))+
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values=c("orange", bae))+
  coord_polar("y", start=0)+
  #facet_grid(cols = vars(read)) +
  theme_classic()

dev.off()

pdf(file=paste0(output,"pie_longread_ratios.pdf"), width=4, height=3)

ggplot(subset(stat.all, read=="Long read"), aes(x="", y=count, fill=ratio_type))+
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values=c("orange", bae))+
  coord_polar("y", start=0)+
  #facet_grid(cols = vars(read)) +
  theme_classic()

dev.off()

# venn diagrams -------------------------------------------------------

v<-venn.diagram(
  x = list(filter(br_autosomes, read=="long")$name, filter(br_autosomes, read=="short")$name),
  category.names = c("Long reads" , "Short reads"),
  filename = NULL
)

pdf(file=paste0(output,"venn_all_overlap.pdf"), width=4, height=4)

grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), v))

dev.off()

br_over_ase<-subset(br_over, type_long!="BAE"|type_short!="BAE")
v<-venn.diagram(
  x = list(filter(br_over_ase, type_long!="BAE")$name, filter(br_over_ase, type_short!="BAE")$name),
  category.names = c("Long reads" , "Short reads"),
  filename = NULL
)

pdf(file=paste0(output,"venn_ASE_overlap.pdf"), width=4, height=4)

grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), v))

dev.off()

br_over_bae<-subset(br_over, type_long=="BAE"|type_short=="BAE")
v<-venn.diagram(
  x = list(filter(br_over_bae, type_long=="BAE")$name, filter(br_over_bae, type_short=="BAE")$name),
  category.names = c("Long reads" , "Short reads"),
  filename = NULL
)

pdf(file=paste0(output,"venn_BAE_overlap.pdf"), width=4, height=4)

grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), v))

dev.off()

#draw.pairwise.venn(area1 = length(filter(br_autosomes, read=="long")), area2 = length(filter(br_autosomes, read=="short")), cross.area = , category = c("Long reads" , "Short reads"))

# violin --------------------------------------------------------------
pdf(file=paste0(output,"violin_autosomes.pdf"), width=4, height=3)

ggplot(br_autosomes, aes(x=read, y=ratio))+
  geom_violin(trim = F)+
  scale_fill_manual(values=c("#565656", "#ADBFC5"))+
  geom_boxplot(width=0.1, coef = 6, aes(x=read, y=ratio), fill="white")+
  geom_jitter(data= subset(br_autosomes, ratio%in%boxplot.stats((br_autosomes%>%
                                                                   filter(!(name%in%imp_genes$name))%>%
                                                                   filter(read=="long"))$ratio)$out),
              color=bae, shape=20, size=0.5, alpha=0.25, width=0.05)+
  geom_jitter(data= subset(br_autosomes, ratio%in%boxplot.stats((br_autosomes%>%
                                                                   filter(!(name%in%imp_genes$name))%>%
                                                                   filter(read=="short"))$ratio)$out),
              color=bae, shape=20, size=0.5, alpha=0.25, width=0.05)+
  ylim(0, 1)+
  theme(aspect.ratio=1)+
  theme_classic()

dev.off()

