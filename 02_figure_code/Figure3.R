#**********************************************************************
#* Figure 3 script
#* Project: allelic_longreads
#* Author: Lison Lemoine
#* Date: 04.2025
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
allelic_ratio.type <- function(gen.sum){
  gen.sum$type <- "BAE"
  
  for (i in 1:nrow(gen.sum)) {
    if (gen.sum[i,]$allelic_ratio >= 0.7){
      gen.sum[i,]$type <- "MAT"
    }
    else if (gen.sum[i,]$allelic_ratio <= 0.3){
      gen.sum[i,]$type <- "PAT"
    }
  }
  return(gen.sum)
}
create_venn_diagram <- function(n_long, n_short, n_overlap, 
                                percentage, adjust_overlap) {
  percent_overlap <- round(n_overlap / n_long * 100)
  
  scale_factor <- 0.3
  rl <- scale_factor * sqrt(n_long) / max(sqrt(n_long), sqrt(n_short))
  rs <- scale_factor * sqrt(n_short) / max(sqrt(n_long), sqrt(n_short))  
  
  d <- adjust_overlap(percent_overlap) * (rl + rs)
  
  if (n_long > n_short){
    short_xpos <- 0.4 + d 
    long_xpos <- 0.4
  } else {
    short_xpos <- 0.4
    long_xpos <- 0.4 + d 
  }
  
  if (percentage){
    overlap_lab <- bquote(atop(.(percent_overlap)*"%",
                               italic(n) == .(comma(n_overlap))))
  } else {
    overlap_lab <- bquote(atop("Overlap", (italic(n) == .(comma(n_overlap)))))
  }
  
  plot <- grobTree(
    rectGrob(gp = gpar(fill = "white", col = NA)),
    
    # Draw full colors first
    circleGrob(x = short_xpos, y = 0.5, r = rs, 
               gp = gpar(fill = SHORT_COL, col = "black", lwd = 1)) , 
    circleGrob(x = long_xpos, y = 0.5, r = rl, 
               gp = gpar(fill = LONG_COL, col = "black", lwd = 1)) ,
    
    # Then add the circles again but transparent
    circleGrob(x = short_xpos, y = 0.5, r = rs, 
               gp = gpar(fill = alpha(SHORT_COL, 0.7), col = "black", lwd = 1)),
    circleGrob(x = long_xpos, y = 0.5, r = rl, 
               gp = gpar(fill = alpha(LONG_COL, 0.7), col = "black", lwd = 1)),   
    
    textGrob(
      label = overlap_lab,
      x = long_xpos / 1.4, y = 0.5, 
      gp = gpar(col = "white", cex = 1.6, fontface = "bold"),
    ),
    textGrob(
      label = bquote("Long reads" ~ (italic(n) == .(comma(n_long)))),
      x = 0.75, y = 0.1,
      gp = gpar(col = LONG_COL, cex = 1.4)
    ),
    textGrob(
      label = bquote("Short reads" ~ (italic(n) == .(comma(n_short)))),
      x = 0.3, y = 0.9,
      gp = gpar(col = SHORT_COL, cex = 1.4)
    )
  )
  return(plot)
}

# init ---------------------------------
load.packages(c("ggplot2", "tidyverse"))

# config --------------------------------------------------------------
input <- "/Users/lisonlemoine/Documents/GitHub/Allelic_Longreads/02_figure_code/"

PAT <- "#283D90" # "#40509E"
MAT <- "#BE222F" # "#CE3730"
BAE <- "#cccccc" 
LONG_COL <- "#565656"
SHORT_COL <- "#ADBFC5"
ALT_BAE <- "#1D6433"
ASE <- "#F68840"

# read tables ----------------------------------------------------------
isoforms<-data.frame(isoform=rep(c("Nesp", "Ex1A", "Gsα"), 2), 
                     num_reads=c(6, 0, 109, #MAT counts manually counted
                                 0, 26, 98), #PAT counts manually counted
                     type=rep(c("MAT", "PAT"), each=3))

isoforms$isoform <- factor(isoforms$isoform, levels = c( "Gsα", "Nesp", "Ex1A"))  


# Figure 3c -----------------------------------------------------------
figure3c <- ggplot(isoforms, aes(x=isoform, y=num_reads, fill=type)) +
  geom_bar( stat="identity", position=position_dodge()) +  
  scale_fill_manual(values=c(mat, pat)) +
  geom_text(aes(label = num_reads), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 4)+
  labs(x="Gnas isoforms", y = "Number of reads")+
  theme_classic()

# print figures -------------------------------------------------------
pdf(file="Figur3.pdf", width=10, height=5)
grid.draw(figure2a_right)
grid.text("C", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

figure3c

dev.off()


# R session -----------------------------------------------------------
devtools::session_info()

# ─ Session info ───────────────────────────────────────────────────
# setting  value
# version  R version 4.4.3 (2025-02-28)
# os       macOS Sequoia 15.3.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/Berlin
# date     2025-04-17
# rstudio  2024.12.1+563 Kousa Dogwood (desktop)
# pandoc   NA
# quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ───────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# abind         1.4-8   2024-09-12 [1] CRAN (R 4.4.1)
# backports     1.5.0   2024-05-23 [1] CRAN (R 4.4.0)
# broom         1.0.8   2025-03-28 [1] CRAN (R 4.4.1)
# cachem        1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
# car           3.1-3   2024-09-27 [1] CRAN (R 4.4.1)
# carData       3.0-5   2022-01-06 [1] CRAN (R 4.4.0)
# cli           3.6.4   2025-02-13 [1] CRAN (R 4.4.1)
# colorspace    2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
# crayon        1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
# devtools      2.4.5   2022-10-11 [1] CRAN (R 4.4.0)
# digest        0.6.37  2024-08-19 [1] CRAN (R 4.4.1)
# dplyr       * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
# ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.4.0)
# farver        2.1.2   2024-05-13 [1] CRAN (R 4.4.0)
# fastmap       1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
# Formula       1.2-5   2023-02-24 [1] CRAN (R 4.4.0)
# fs            1.6.6   2025-04-12 [1] CRAN (R 4.4.1)
# generics      0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
# ggplot2     * 3.5.2   2025-04-09 [1] CRAN (R 4.4.1)
# ggpubr        0.6.0   2023-02-10 [1] CRAN (R 4.4.0)
# ggsignif      0.6.4   2022-10-13 [1] CRAN (R 4.4.0)
# glue          1.8.0   2024-09-30 [1] CRAN (R 4.4.1)
# gtable        0.3.6   2024-10-25 [1] CRAN (R 4.4.1)
# htmltools     0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
# htmlwidgets   1.6.4   2023-12-06 [1] CRAN (R 4.4.0)
# httpuv        1.6.15  2024-03-26 [1] CRAN (R 4.4.0)
# labeling      0.4.3   2023-08-29 [1] CRAN (R 4.4.0)
# later         1.4.2   2025-04-08 [1] CRAN (R 4.4.1)
# lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
# memoise       2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
# mime          0.13    2025-03-17 [1] CRAN (R 4.4.1)
# miniUI        0.1.1.1 2018-05-18 [1] CRAN (R 4.4.0)
# munsell       0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
# pillar        1.10.2  2025-04-05 [1] CRAN (R 4.4.1)
# pkgbuild      1.4.7   2025-03-24 [1] CRAN (R 4.4.1)
# pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
# pkgload       1.4.0   2024-06-28 [1] CRAN (R 4.4.0)
# profvis       0.4.0   2024-09-20 [1] CRAN (R 4.4.1)
# promises      1.3.2   2024-11-28 [1] CRAN (R 4.4.1)
# purrr         1.0.4   2025-02-05 [1] CRAN (R 4.4.1)
# R6            2.6.1   2025-02-15 [1] CRAN (R 4.4.1)
# Rcpp          1.0.14  2025-01-12 [1] CRAN (R 4.4.1)
# remotes       2.5.0   2024-03-17 [1] CRAN (R 4.4.0)
# rlang         1.1.6   2025-04-11 [1] CRAN (R 4.4.1)
# rstatix       0.7.2   2023-02-01 [1] CRAN (R 4.4.0)
# rstudioapi    0.17.1  2024-10-22 [1] CRAN (R 4.4.1)
# scales        1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
# sessioninfo   1.2.3   2025-02-05 [1] CRAN (R 4.4.1)
# shiny         1.10.0  2024-12-14 [1] CRAN (R 4.4.1)
# tibble        3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
# tidyr       * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
# tidyselect    1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
# urlchecker    1.0.1   2021-11-30 [1] CRAN (R 4.4.0)
# usethis       3.1.0   2024-11-26 [1] CRAN (R 4.4.1)
# utf8          1.2.4   2023-10-22 [1] CRAN (R 4.4.0)
# vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
# withr         3.0.2   2024-10-28 [1] CRAN (R 4.4.1)
# xtable        1.8-4   2019-04-21 [1] CRAN (R 4.4.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
# * ── Packages attached to the search path.
# 
# ──────────────────────────────────────────────────────────────────
