#**********************************************************************
#* Long read allelic ratio analysis (& short read corr)
#* Project: allelic_longreads
#* Author: Lison Lemoine
#* Date: 04.2023
#**********************************************************************
# functions -----------------------------------------------------------
load.packages <- function(packages){
  has.package <- packages %in% rownames(installed.packages())
  if(any(!has.package)){
    message("Installing packages...")
    setRepositories(ind = c(1:5))
    install.packages(packages[!has.package], 
                     repos = "https://cran.uni-muenster.de/", 
                     dependencies = T)
  } 
  for(r in packages){
    suppressMessages(require(r,character.only = TRUE))
  }
}

calc.allelic.score <- function(x) {
  if (all(x==0) || any(is.na(x)))
    return(0)
  if (sum(x) < 10)
    return(0)
  allelic.score <- log10(pbinom(min(matrix(x,ncol=2)), x[1]+x[2], 0.5, 
                                lower.tail=TRUE, log = FALSE)) 
  ifelse(x[1] > x[2], -allelic.score,allelic.score) 
}

# init ----------------------------------------------------------------
load.packages(c("gtools", "dplyr", "tidyverse", "openxlsx"))

# var -----------------------------------------------------------------
path_to_folder <- "~/longread_ASE"
name <- "brain"

input <- paste0(path_to_folder, "/output/")
output <- paste0(path_to_folder, "/plots/")

short_reads <- paste0(path_to_folder, "/input/locus_table.txt")
refseqID <- paste0(path_to_folder, "/input/Refseq_IDs_unique_max_locus.txt")

long_read1=paste0(input, name, "_h1.count")
long_read2=paste0(input, name, "_h2.count")

# read LR -------------------------------------------------------
count_mat <- read.delim(long_read1, header=F)
count_pat <- read.delim(long_read2, header=F)
colnames(count_mat) <- c("name", "reads_h1")
colnames(count_pat) <- c("name", "reads_h2")

long_counts <- full_join(count_mat, count_pat, by=join_by(name))
long_counts <- long_counts[long_counts$name!="__ambiguous"&long_counts$name!="__no_feature",]
rm(count_mat, count_pat)

# calculate ratio -----------------------------------------------------
long_counts$ratio_lr <- round(apply(long_counts[,c("reads_h1","reads_h2")], 1,
                                    function(x) {return(x["reads_h1"]/(x["reads_h1"]+x["reads_h2"]))}),2)
long_counts <- na.omit(long_counts) #this removes 0 reads
long_counts$lr_total <- apply(long_counts[,c("reads_h1","reads_h2")], 1,
                              function(x) {return((x["reads_h1"]+x["reads_h2"]))})

# filter lr -----------------------------------------------------------
long_counts <- long_counts[long_counts$lr_total>=10,]

gen_names <- read.delim(refseqID, header=T)

long_counts <- left_join(long_counts, gen_names, by=join_by(name)) %>%
  dplyr::select(-c(strand, start, end))
remove(gen_names)

# short reads ---------------------------------------------------------
short_counts <- read.delim(short_reads)
colnames(short_counts) <- c("chr", "start", "end", "name", "A1_reads", 
                          "A2_reads", "total_reads", "allelic_score", 
                          "allelic_ratio")
short_counts <- na.omit(short_counts[short_counts$total_reads>=20,])

# merge ---------------------------------------------------------------
lrsr <- left_join(short_counts,long_counts, by=join_by(name), copy=T) %>%
  mutate(chr=chr.x)%>%
  dplyr::select(-c(chr.x, chr.y))
lrsr <- na.omit(lrsr)

# all overlapping autosomes table -------------------------------------
autosomes <- rbind(long_counts %>%
                     filter(chr!="chrX") %>%
                     mutate(read = "long", total = lr_total, 
                            ratio = ratio_lr) %>%
                     dplyr::select(c(name, chr, total, ratio, read)),
                   short_counts %>%
                     mutate(read = "short", total = total_reads, 
                            ratio=allelic_ratio) %>%
                     filter(chr!="chrX"&chr!="chrY") %>%
                     dplyr::select(c(name, chr, total, ratio, read)) 
)

# sex chr filter ------------------------------------------------------
lrsr <- lrsr %>% filter(chr!="chrX") %>% filter(chr!="chrY")
rm(long_counts, short_counts)

# allelic score -------------------------------------------------------
lrsr$score_lr <- round(apply(lrsr[,c("reads_h1","reads_h2")], 1,
                             calc.allelic.score),3)

# Loci: Replacing +/-infinitive with the highest/lowest finite value
max_value <- range(abs(lrsr$score_lr),finite=TRUE)[2]
lrsr$score_lr[lrsr$score_lr == "Inf"] <- max_value
lrsr$score_lr[lrsr$score_lr == "-Inf"] <- -max_value

lrsr <- na.omit(lrsr)

rm(max_value, calc.allelic.score)

# table clean up ------------------------------------------------------
## thinner table
lrsr <- subset(lrsr, select = -c(start, end, A1_reads, A2_reads))

## reorganize
lrsr <- lrsr %>% relocate(name) %>% relocate(chr) %>%
  relocate(lr_total, .after=reads_h2) %>%
  relocate(score_lr, .after=ratio_lr) %>%
  relocate(total_reads, .after=score_lr) %>%
  relocate(allelic_ratio, .after=total_reads) %>%
  relocate(allelic_score, .after=allelic_ratio)

colnames(lrsr) <- c("chr", "name", "reads_h1", "reads_h2", "t_reads_lr", 
                  "ratio_lr", "score_lr", "t_reads_sr", "ratio_sr", "score_sr")

# to excel sheet ------------------------------------------------------
ase_overlap <- autosomes %>% filter(name%in%br_lrsr$name) %>%
  pivot_wider(names_from = read, values_from = c(read_BL6, read_CAST, 
                                                 total, ratio, type)) %>%
  filter(type_long!="BAE"&type_short!="BAE")
ase_overlap$consistent <- ase_overlap$type_short==ase_overlap$type_long

br_lrsr <- dplyr::select(br_lrsr, c(name, chr, start, end, reads_h1, reads_h2, 
                                  t_reads_lr, ratio_lr, A1_reads, A2_reads, 
                                  t_reads_sr, ratio_sr))
colnames(br_lrsr) <- c("name", "chr", "start", "end", "BL6 long reads", 
                     "CAST long reads", "total long reads", 
                     "long reads allelic ratio", 
                     "BL6 short reads", "CAST short reads",
                     "total short reads", "short reads allelic ratio")

long_counts <- dplyr::select(long_counts, c(name, chr, start, end, reads_h1, reads_h2, 
                                    lr_total, ratio_lr))
colnames(long_counts) <- c("name", "chr", "start", "end", "BL6 reads", "CAST reads", 
                      "total reads", "allelic ratio")

short_counts <- dplyr::select(short_counts, c(name, chr, start, end, A1_reads, 
                              A2_reads, total_reads, allelic_ratio))
colnames(short_counts) <- c("name", "chr", "start", "end", "BL6 reads", "CAST reads", 
                   "total reads", "allelic ratio")

ase_overlap <- dplyr::select(ase_overlap, c(name, chr, start, end, read_BL6_long, 
                                          read_CAST_long, total_long, 
                                          ratio_long, type_long, 
                                          read_BL6_short, read_CAST_short, 
                                          total_short, ratio_short, type_short, 
                                          consistent))
colnames(ase_overlap) <- c("name", "chr", "start", "end", "BL6 long reads", 
                         "CAST long reads", "total long reads", 
                         "long reads allelic ratio", 
                         "allelic ratio category long reads", 
                         "BL6 short reads", "CAST short reads", 
                         "total short reads", 
                         "short reads allelic ratio", 
                         "allelic ratio category short reads",
                         "is consistent between short and long reads")

dataset_names <- list('a Long reads' = long_counts, 
                      'b Short reads' = short_counts, 
                      'c Overlap' = lrsr, 'd ASE Overlap' = ase_overlap)
write.xlsx(dataset_names, file = 'Supplementary_table1.xlsx')


# R session -----------------------------------------------------------
devtools::session_info()

## ─ Session info ──────────────────────────────────────────────────────────────────────
## setting  value
## version  R version 4.4.1 (2024-06-14)
## os       macOS 15.1
## system   aarch64, darwin20
## ui       RStudio
## language (EN)
## collate  en_US.UTF-8
## ctype    en_US.UTF-8
## tz       Europe/Berlin
## date     2024-11-15
## rstudio  2024.09.0+375 Cranberry Hibiscus (desktop)
## pandoc   NA
## 
## ─ Packages ──────────────────────────────────────────────────────────────────────────
## package          * version date (UTC) lib source
## abind              1.4-5   2016-07-21 [1] CRAN (R 4.4.0)
## backports          1.5.0   2024-05-23 [1] CRAN (R 4.4.0)
## bit                4.0.5   2022-11-15 [1] CRAN (R 4.4.0)
## bit64              4.0.5   2020-08-30 [1] CRAN (R 4.4.0)
## blob               1.2.4   2023-03-17 [1] CRAN (R 4.4.0)
## broom              1.0.6   2024-05-17 [1] CRAN (R 4.4.0)
## cachem             1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
## car                3.1-2   2023-03-30 [1] CRAN (R 4.4.0)
## carData            3.0-5   2022-01-06 [1] CRAN (R 4.4.0)
## cli                3.6.3   2024-06-21 [1] CRAN (R 4.4.0)
## colorspace         2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
## crayon             1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
## curl               5.2.1   2024-03-01 [1] CRAN (R 4.4.0)
## DBI                1.2.3   2024-06-02 [1] CRAN (R 4.4.0)
## dbplyr             2.5.0   2024-03-19 [1] CRAN (R 4.4.0)
## devtools           2.4.5   2022-10-11 [1] CRAN (R 4.4.0)
## digest             0.6.37  2024-08-19 [1] CRAN (R 4.4.1)
## dplyr            * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
## ellipsis           0.3.2   2021-04-29 [1] CRAN (R 4.4.0)
## fansi              1.0.6   2023-12-08 [1] CRAN (R 4.4.0)
## fastmap            1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
## filelock           1.0.3   2023-12-11 [1] CRAN (R 4.4.0)
## forcats          * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
## fs                 1.6.4   2024-04-25 [1] CRAN (R 4.4.0)
## generics           0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
## ggplot2            3.5.1   2024-04-23 [1] CRAN (R 4.4.0)
## ggpubr             0.6.0   2023-02-10 [1] CRAN (R 4.4.0)
## ggrepel            0.9.5   2024-01-10 [1] CRAN (R 4.4.0)
## ggsignif           0.6.4   2022-10-13 [1] CRAN (R 4.4.0)
## glue               1.7.0   2024-01-09 [1] CRAN (R 4.4.0)
## gridExtra          2.3     2017-09-09 [1] CRAN (R 4.4.0)
## gtable             0.3.5   2024-04-22 [1] CRAN (R 4.4.0)
## gtools           * 3.9.5   2023-11-20 [1] CRAN (R 4.4.0)
## hms                1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
## htmltools          0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
## htmlwidgets        1.6.4   2023-12-06 [1] CRAN (R 4.4.0)
## httpuv             1.6.15  2024-03-26 [1] CRAN (R 4.4.0)
## httr               1.4.7   2023-08-15 [1] CRAN (R 4.4.0)
## httr2              1.0.2   2024-07-16 [1] CRAN (R 4.4.0)
## jsonlite           1.8.8   2023-12-04 [1] CRAN (R 4.4.0)
## later              1.3.2   2023-12-06 [1] CRAN (R 4.4.0)
## lifecycle          1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
## lubridate        * 1.9.3   2023-09-27 [1] CRAN (R 4.4.0)
## magrittr           2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
## memoise            2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
## mime               0.12    2021-09-28 [1] CRAN (R 4.4.0)
## miniUI             0.1.1.1 2018-05-18 [1] CRAN (R 4.4.0)
## munsell            0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
## openxlsx         * 4.2.7   2024-08-30 [1] CRAN (R 4.4.1)
## pheatmap           1.0.12  2019-01-04 [1] CRAN (R 4.4.0)
## pillar             1.9.0   2023-03-22 [1] CRAN (R 4.4.0)
## pkgbuild           1.4.4   2024-03-17 [1] CRAN (R 4.4.0)
## pkgconfig          2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
## pkgload            1.4.0   2024-06-28 [1] CRAN (R 4.4.0)
## png                0.1-8   2022-11-29 [1] CRAN (R 4.4.0)
## prettyunits        1.2.0   2023-09-24 [1] CRAN (R 4.4.0)
## profvis            0.3.8   2023-05-02 [1] CRAN (R 4.4.0)
## progress           1.2.3   2023-12-06 [1] CRAN (R 4.4.0)
## promises           1.3.0   2024-04-05 [1] CRAN (R 4.4.0)
## purrr            * 1.0.2   2023-08-10 [1] CRAN (R 4.4.0)
## R6                 2.5.1   2021-08-19 [1] CRAN (R 4.4.0)
## rappdirs           0.3.3   2021-01-31 [1] CRAN (R 4.4.0)
## RColorBrewer       1.1-3   2022-04-03 [1] CRAN (R 4.4.0)
## Rcpp               1.0.13  2024-07-17 [1] CRAN (R 4.4.0)
## readr            * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
## remotes            2.5.0   2024-03-17 [1] CRAN (R 4.4.0)
## rlang              1.1.4   2024-06-04 [1] CRAN (R 4.4.0)
## RSQLite            2.3.7   2024-05-27 [1] CRAN (R 4.4.0)
## rstatix            0.7.2   2023-02-01 [1] CRAN (R 4.4.0)
## rstudioapi         0.16.0  2024-03-24 [1] CRAN (R 4.4.0)
## S4Vectors          0.42.1  2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
## scales             1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
## sessioninfo        1.2.2   2021-12-06 [1] CRAN (R 4.4.0)
## shiny              1.9.1   2024-08-01 [1] CRAN (R 4.4.0)
## stringi            1.8.4   2024-05-06 [1] CRAN (R 4.4.0)
## stringr          * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
## tibble           * 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
## tidyr            * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
## tidyselect         1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
## tidyverse        * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
## timechange         0.3.0   2024-01-18 [1] CRAN (R 4.4.0)
## tzdb               0.4.0   2023-05-12 [1] CRAN (R 4.4.0)
## UCSC.utils         1.0.0   2024-05-06 [1] Bioconductor 3.19 (R 4.4.0)
## urlchecker         1.0.1   2021-11-30 [1] CRAN (R 4.4.0)
## usethis            3.0.0   2024-07-29 [1] CRAN (R 4.4.0)
## utf8               1.2.4   2023-10-22 [1] CRAN (R 4.4.0)
## vctrs              0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
## withr              3.0.1   2024-07-31 [1] CRAN (R 4.4.0)
## xml2               1.3.6   2023-12-04 [1] CRAN (R 4.4.0)
## xtable             1.8-4   2019-04-21 [1] CRAN (R 4.4.0)
## XVector            0.44.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
## zlibbioc           1.50.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
## 
## [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
## 
## ─────────────────────────────────────────────────────────────────────────────────────

