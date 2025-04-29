#**********************************************************************
#* Figure 2 script
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

# init ----------------------------------------------------------------
load.packages(c("gtools", "ggpubr", "RColorBrewer", "dplyr", 
                "tidyverse", "gridExtra", "ggrepel", "readxl", 
                "scales", "VennDiagram"))

# config --------------------------------------------------------------
input <- "~/Documents/GitHub/Allelic_Longreads/02_figure_code/"

PAT <- "#283D90" # "#40509E"
MAT <- "#BE222F" # "#CE3730"
BAE <- "#cccccc" 
LONG_COL <- "#565656"
SHORT_COL <- "#ADBFC5"
ALT_BAE <- "#1D6433"
ASE <- "#F68840"

# read tables ----------------------------------------------------------
lr_ratios <- read_excel(paste0(input,"Supplementary Table S1.xlsx"), 
                        sheet = "a Long reads")
colnames(lr_ratios) <- c("name", "chr", "start", "end", "read_BL6", 
                         "read_CAST", "total_reads", "allelic_ratio")

sr_ratios <- read_excel(paste0(input,"Supplementary Table S1.xlsx"), 
                        sheet = "b Short reads")

colnames(sr_ratios) <- c("name", "chr", "start", "end", "read_BL6", 
                         "read_CAST", "total_reads", "allelic_ratio")

# merge ---------------------------------------------------------------
all_ratios <- left_join(sr_ratios, lr_ratios, 
                        by = c("name", "chr", "start", "end"), 
                        suffix = c("_sr", "_lr"))
all_ratios <- na.omit(all_ratios)

all_ratios <- all_ratios %>% filter(chr != "chrX") %>% filter(chr != "chrY")

# all overlapping autosomes table -------------------------------------
autosomes <- rbind(lr_ratios  %>% 
                     filter(chr != "chrX"&chr != "chrY")  %>% 
                     mutate(read = "Long reads"),
                   sr_ratios  %>% 
                     mutate(read = "Short reads")  %>% 
                     filter(chr != "chrX"&chr != "chrY")
)

# imprinted genes table -----------------------------------------------
imp_genes <- data.frame(name = c("Grb10", "Peg3", "Meg3", "Zdbf2", 
                                 "Airn", "Impact", "Rian", "Mirg", 
                                 "Ago2", "Peg13", "Zrsr1", "Ndn", 
                                 "A330076H08Rik", "Igf2", "Kcnk9", 
                                 "Kcnq1ot1", "Snurf", "Snrpn", "Peg10", 
                                 "Sgce", "B230209E15Rik", "Snhg14", 
                                 "Zim1", "Nap1l5", "Usp29", "Plagl1", 
                                 "Rasgrf1", "Ipw", "A230057D06Rik"), 
                        type = c("PAT", "PAT", "MAT", "PAT", "PAT",
                                 "PAT", "MAT", "MAT", "MAT", "PAT",
                                 "PAT", "PAT", "PAT", "MAT", "MAT", 
                                 "PAT", "PAT", "PAT", "PAT", "PAT", 
                                 "PAT", "PAT", "MAT", "PAT", "PAT", 
                                 "PAT", "PAT", "PAT", "PAT"))

all_ratios <- left_join(all_ratios, imp_genes, by = join_by(name))
all_ratios <- na.replace(all_ratios, "BAE")

imp_ratios <- filter(all_ratios, name%in%imp_genes$name)

# Figure 2a (Right) Venn ----------------------------------------------
figure2a_right <- create_venn_diagram(
  length(filter(autosomes, read == "Long reads")$name), 
  length(filter(autosomes, read == "Short reads")$name), 
  length(intersect(filter(autosomes, read == "Long reads")$name,
                   filter(autosomes, read == "Short reads")$name)), 
  TRUE, function(x){0.21}
)

figure2a_right_raw <- venn.diagram(
  x = list(filter(autosomes, read=="Long reads")$name, 
           filter(autosomes, read=="Short reads")$name),
  category.names = c("Long reads" , "Short reads"),
  filename = NULL
)

# Figure 2b Violin ----------------------------------------------------
figure2b <- ggplot(
  autosomes, 
  aes(x = read, y = allelic_ratio, fill = read)) +
  geom_violin(trim = FALSE, linewidth = 0.3) +
  scale_fill_manual(values = c(LONG_COL, SHORT_COL)) +
  
  geom_boxplot(
    width = 0.1, 
    coef = 6, 
    fill = "white", 
    size = 0.3
  ) +
  
  geom_jitter(
    data = subset(
      autosomes, 
      allelic_ratio %in% boxplot.stats(
        (autosomes %>% 
           filter(!(name %in% imp_genes$name)) %>% 
           filter(read == "Long reads"))$allelic_ratio
      )$out
    ),
    color = BAE, 
    shape = 20, 
    size = 0.5, 
    alpha = 0.25, 
    width = 0.05
  ) +
  
  geom_jitter(
    data = subset(
      autosomes, 
      allelic_ratio %in% boxplot.stats(
        (autosomes %>% 
           filter(!(name %in% imp_genes$name)) %>% 
           filter(read == "Short Reads"))$allelic_ratio
      )$out
    ),
    color = BAE, 
    shape = 20, 
    size = 0.5, 
    alpha = 0.25, 
    width = 0.05
  ) +
  
  ylim(0, 1) +
  
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  
  xlab("") +
  ylab("Allelic Ratio")

# Overlap -------------------------------------------------------------
autosomes <- allelic_ratio.type(autosomes)

overlap <- autosomes %>% 
  filter(name %in% all_ratios$name) %>% 
  mutate(read = ifelse(read=="Long reads", "long", "short")) %>% 
  pivot_wider(
    id_cols = name,
    names_from = read, 
    values_from = c(total_reads, allelic_ratio, type)
  )

# Figure 2c Pies  -----------------------------------------------------
stat.all <- data.frame(count = c(nrow(overlap %>% 
                                        filter(type_long != "BAE")),
                                 nrow(overlap %>% 
                                        filter(type_long == "BAE")),
                                 nrow(overlap %>% 
                                        filter(type_short != "BAE")),
                                 nrow(overlap %>% 
                                        filter(type_short == "BAE"))), 
                       allelic_ratio_type = rep(c("0.7 ≤ Allelic-specific ≤ 0.3", 
                                                  "0.7 > Biallelic > 0.3"), 2), 
                       read = rep(c("Long read", "Short read"), each = 2))

stat.all <- stat.all %>%
  group_by(read) %>%
  mutate(
    percent = count / sum(count) * 100,
    label = paste0(round(percent, 1), "%")
  )

figure2c <- ggplot(
  stat.all, 
  aes(x = "", y = count, fill = allelic_ratio_type)) +
  geom_bar(width = 1, stat = "identity") + 
  scale_fill_manual(values = c(ALT_BAE, ASE))+
  coord_polar("y", start = 0) +
  facet_wrap(~ read, ncol = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 4) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(2, "lines"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Figure 2d Venn ------------------------------------------------------
over_ase <- subset(overlap, type_long != "BAE"|type_short != "BAE")
over_bae <- subset(overlap, type_long == "BAE"|type_short == "BAE")

d_left <- create_venn_diagram(
  length(filter(over_bae, type_long == "BAE")$name), 
  length(filter(over_bae, type_short == "BAE")$name), 
  length(intersect(filter(over_bae, type_long == "BAE")$name, 
                   filter(over_bae, type_short == "BAE")$name)), 
  FALSE, function(x){log10(1 + 9 * (1 - x / 100))}
)

d_right <- create_venn_diagram(
  length(filter(over_ase, type_long != "BAE")$name), 
  length(filter(over_ase, type_short != "BAE")$name), 
  length(intersect(filter(over_ase, type_long != "BAE")$name, 
                   filter(over_ase, type_short != "BAE")$name)), 
  FALSE, function(x){(1 - x / 100)}
)

figure2d_left <- arrangeGrob(d_left, top = 
                             textGrob("Biallelic", 
                                      gp = gpar(fontsize = 16, 
                                                fontface = "bold")))
figure2d_right <- arrangeGrob(d_right, top = 
                              textGrob("Allele-specific", 
                                       gp = gpar(fontsize = 16, 
                                                 fontface = "bold")))

d_left_raw <- venn.diagram(
  x = list(filter(over_bae, type_long == "BAE")$name, 
           filter(over_bae, type_short == "BAE")$name),
  category.names = c("Long reads" , "Short reads"),
  filename = NULL
)

d_right_raw <- venn.diagram(
  x = list(filter(over_ase, type_long != "BAE")$name, 
           filter(over_ase, type_short != "BAE")$name),
  category.names = c("Long reads" , "Short reads"),
  filename = NULL
)

figure2d_left_raw <- arrangeGrob(
  grobs = list(d_left_raw[[1]]), 
  top = textGrob("Biallelic", 
                 gp = gpar(fontsize = 16, 
                           fontface = "bold"))
)

figure2d_right_raw <- arrangeGrob(
  grobs = list(d_right_raw[[1]]), 
  top = textGrob("Biallelic", 
                 gp = gpar(fontsize = 16, 
                           fontface = "bold"))
)

# Figure 2e Scatter ---------------------------------------------------
corr.norm_all <- cor.test(all_ratios$allelic_ratio_lr, 
                          all_ratios$allelic_ratio_sr, 
                          method = "pearson")
corr.norm_imp <- cor.test(imp_ratios$allelic_ratio_lr, 
                          imp_ratios$allelic_ratio_sr, 
                          method = "pearson")

goi <- c("Gnas", "Ube2g2", "Malat1")

figure2e_base_layer <- ggplot(
  all_ratios %>% filter(!(name %in% imp_genes$name)), 
  aes(x = allelic_ratio_sr, y = allelic_ratio_lr, color = type, shape = type)
) +
  geom_point(position = position_jitter(w = 0.005), alpha = 0.2) +
  
  scale_shape_manual(values = c(20, 24, 25)) +
  scale_color_manual(values = c(BAE, MAT, PAT)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = 'off') +
  
  annotate("text", 
           x = 0.2, 
           y = all_ratios[all_ratios$name %in% goi, ]$allelic_ratio_sr, 
           label = all_ratios[all_ratios$name %in% goi, ]$name, 
           fontface = 1) +
  
  annotate("segment", 
           x = 0.3,
           xend = all_ratios[all_ratios$name %in% goi, ]$allelic_ratio_sr,
           y = all_ratios[all_ratios$name %in% goi, ]$allelic_ratio_sr, 
           yend = all_ratios[all_ratios$name %in% goi, ]$allelic_ratio_lr) 

figure2e <- figure2e_base_layer +
  geom_point(
    data = imp_ratios, 
    aes(x = allelic_ratio_sr, y = allelic_ratio_lr, color = type, shape = type),
    position = position_jitter(h = 0.005, w = 0.005)
  ) +
  
  geom_text_repel(
    data = imp_ratios,
    aes(x = allelic_ratio_sr, y = allelic_ratio_lr, label = name, color = type),
    force_pull = 0,
    nudge_y = 0.1,
    direction = "x",
    hjust = 4,
    segment.size = 0.1,
    size = 2.5
  ) + 
  
  ggtitle(
    paste0(
      "All genes (r = ", round(corr.norm_all$estimate, 3), "), ",
      "Imprinted genes (r = ", round(corr.norm_imp$estimate, 3), ")"
    )
  ) +
  xlab("Allelic ratio short reads") +
  ylab("Allelic ratio long reads") +
  
  scale_x_continuous(labels = function(x) sprintf("%.2g", x)) +
  scale_y_continuous(labels = function(y) sprintf("%.2g", y)) +
  
  theme_classic() +
  theme(axis.text=element_text(size=8),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background=element_rect(
          fill = NA,
          color = NA
        ),
        axis.title = element_text(size=12),
        aspect.ratio = 1)

# print figures -------------------------------------------------------
pdf(file="Figure2.pdf", width=10, height=5)
grid.draw(figure2a_right)
grid.text("A", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

figure2b
grid.text("B", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

figure2c
grid.text("C", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))
                     
grid.arrange(figure2d_left, figure2d_right, ncol = 2)
grid.text("D", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

figure2e
grid.text("E", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

dev.off()
                     
pdf(file="Figure2_raw.pdf", width=10, height=5)
grid.draw(figure2a_right_raw)
grid.text("A raw", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))

grid.arrange(figure2d_left_raw, figure2d_right_raw, ncol = 2)
grid.text("D raw", x = 0.02, y = 0.98, just = c("left", "top"),
          gp = gpar(fontsize = 12, fontface = "bold"))
             
dev.off()
                     
                

# R session -----------------------------------------------------------
devtools::session_info()

# ─ Session info ───────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.3 (2025-02-28)
# os       macOS Sequoia 15.3.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/Berlin
# date     2025-04-16
# rstudio  2024.12.1+563 Kousa Dogwood (desktop)
# pandoc   NA
# quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ───────────────────────────────────────────────────────────
# package        * version date (UTC) lib source
# abind            1.4-8   2024-09-12 [1] CRAN (R 4.4.1)
# backports        1.5.0   2024-05-23 [1] CRAN (R 4.4.0)
# broom            1.0.8   2025-03-28 [1] CRAN (R 4.4.1)
# cachem           1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
# car              3.1-3   2024-09-27 [1] CRAN (R 4.4.1)
# carData          3.0-5   2022-01-06 [1] CRAN (R 4.4.0)
# cellranger       1.1.0   2016-07-27 [1] CRAN (R 4.4.0)
# cli              3.6.4   2025-02-13 [1] CRAN (R 4.4.1)
# colorspace       2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
# devtools         2.4.5   2022-10-11 [1] CRAN (R 4.4.0)
# digest           0.6.37  2024-08-19 [1] CRAN (R 4.4.1)
# dplyr          * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
# ellipsis         0.3.2   2021-04-29 [1] CRAN (R 4.4.0)
# fastmap          1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
# forcats        * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
# formatR          1.14    2023-01-17 [1] CRAN (R 4.4.0)
# Formula          1.2-5   2023-02-24 [1] CRAN (R 4.4.0)
# fs               1.6.6   2025-04-12 [1] CRAN (R 4.4.1)
# futile.logger  * 1.4.3   2016-07-10 [1] CRAN (R 4.4.0)
# futile.options   1.0.1   2018-04-20 [1] CRAN (R 4.4.0)
# generics         0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
# ggplot2        * 3.5.2   2025-04-09 [1] CRAN (R 4.4.1)
# ggpubr         * 0.6.0   2023-02-10 [1] CRAN (R 4.4.0)
# ggrepel        * 0.9.6   2024-09-07 [1] CRAN (R 4.4.1)
# ggsignif         0.6.4   2022-10-13 [1] CRAN (R 4.4.0)
# glue             1.8.0   2024-09-30 [1] CRAN (R 4.4.1)
# gridExtra      * 2.3     2017-09-09 [1] CRAN (R 4.4.0)
# gtable           0.3.6   2024-10-25 [1] CRAN (R 4.4.1)
# gtools         * 3.9.5   2023-11-20 [1] CRAN (R 4.4.0)
# hms              1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
# htmltools        0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
# htmlwidgets      1.6.4   2023-12-06 [1] CRAN (R 4.4.0)
# httpuv           1.6.15  2024-03-26 [1] CRAN (R 4.4.0)
# lambda.r         1.2.4   2019-09-18 [1] CRAN (R 4.4.0)
# later            1.4.2   2025-04-08 [1] CRAN (R 4.4.1)
# lifecycle        1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
# lubridate      * 1.9.4   2024-12-08 [1] CRAN (R 4.4.1)
# magrittr         2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
# memoise          2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
# mime             0.13    2025-03-17 [1] CRAN (R 4.4.1)
# miniUI           0.1.1.1 2018-05-18 [1] CRAN (R 4.4.0)
# munsell          0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
# pheatmap       * 1.0.12  2019-01-04 [1] CRAN (R 4.4.0)
# pillar           1.10.2  2025-04-05 [1] CRAN (R 4.4.1)
# pkgbuild         1.4.7   2025-03-24 [1] CRAN (R 4.4.1)
# pkgconfig        2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
# pkgload          1.4.0   2024-06-28 [1] CRAN (R 4.4.0)
# profvis          0.4.0   2024-09-20 [1] CRAN (R 4.4.1)
# promises         1.3.2   2024-11-28 [1] CRAN (R 4.4.1)
# purrr          * 1.0.4   2025-02-05 [1] CRAN (R 4.4.1)
# R6               2.6.1   2025-02-15 [1] CRAN (R 4.4.1)
# RColorBrewer   * 1.1-3   2022-04-03 [1] CRAN (R 4.4.0)
# Rcpp             1.0.14  2025-01-12 [1] CRAN (R 4.4.1)
# readr          * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
# readxl         * 1.4.5   2025-03-07 [1] CRAN (R 4.4.1)
# remotes          2.5.0   2024-03-17 [1] CRAN (R 4.4.0)
# rlang            1.1.6   2025-04-11 [1] CRAN (R 4.4.1)
# rstatix          0.7.2   2023-02-01 [1] CRAN (R 4.4.0)
# rstudioapi       0.17.1  2024-10-22 [1] CRAN (R 4.4.1)
# scales           1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
# sessioninfo      1.2.3   2025-02-05 [1] CRAN (R 4.4.1)
# shiny            1.10.0  2024-12-14 [1] CRAN (R 4.4.1)
# stringi          1.8.7   2025-03-27 [1] CRAN (R 4.4.1)
# stringr        * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
# tibble         * 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
# tidyr          * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
# tidyselect       1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
# tidyverse      * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
# timechange       0.3.0   2024-01-18 [1] CRAN (R 4.4.0)
# tzdb             0.5.0   2025-03-15 [1] CRAN (R 4.4.1)
# urlchecker       1.0.1   2021-11-30 [1] CRAN (R 4.4.0)
# usethis          3.1.0   2024-11-26 [1] CRAN (R 4.4.1)
# vctrs            0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
# VennDiagram    * 1.7.3   2022-04-12 [1] CRAN (R 4.4.0)
# withr            3.0.2   2024-10-28 [1] CRAN (R 4.4.1)
# xtable           1.8-4   2019-04-21 [1] CRAN (R 4.4.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
# * ── Packages attached to the search path.
# 
# ──────────────────────────────────────────────────────────────────────
