#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

fileloc <- args[1]

library(tidyverse)

# Load data
df <- read_tsv(fileloc)

#Get separate gene columns.
df <- df %>%
  mutate(gene = sub(":.*", "", Target_ID))
  
# Helper function to print counts
count_info <- function(data, label) {
  cat(label)
  cat(sprintf("\n  Rows: %d \n", nrow(data)))
  cat(sprintf("  Unique miRNA: %d \n", n_distinct(data$miRNA)))
  cat(sprintf("  Unique mRNA: %d \n", n_distinct(data$gene)))
  cat("\n")
  cat("\n")
}

# Default thresholds
f <- df %>% filter(Max_Score >= 140, Max_Energy_kcal_mol <= -20)
count_info(f, "Default threshold: Score > = 140 AND Energy < = -20")

# a) minimum score 150
f_a <- df %>% filter(Max_Score >= 150)
count_info(f_a, "a) Score > = 150")

# b) minimum score 160
f_b <- df %>% filter(Max_Score >= 160)
count_info(f_b, "b) Score > = 160")

# c) minimum score 150 AND maximum energy -25
f_c <- df %>% filter(Max_Score >= 150, Max_Energy_kcal_mol <= -25)
count_info(f_c, "c) Score > = 150 AND Energy < = -25")

# d) minimum score 160 AND maximum energy -25
f_d <- df %>% filter(Max_Score >= 160, Max_Energy_kcal_mol <= -25)
count_info(f_d, "d) Score > = 160 AND Energy < = -25")

# e) maximum energy -25
f_e <- df %>% filter(Max_Energy_kcal_mol <= -25)
count_info(f_e, "e) Energy < = -25")

# f) minimum score 150 AND maximum energy -30
f_f <- df %>% filter(Max_Score >= 150, Max_Energy_kcal_mol <= -30)
count_info(f_f, "f) Score > = 150 AND Energy < = -30")

# g) minimum score 160 AND maximum energy -30
f_g <- df %>% filter(Max_Score >= 160, Max_Energy_kcal_mol <= -30)
count_info(f_g, "g) Score > = 160 AND Energy < = -30")

# h) maximum energy -30
f_h <- df %>% filter(Max_Energy_kcal_mol <= -30)
count_info(f_h, "h) Energy < = -30")

  

# Scatter plot with trendline
# Extract directory of input file
outdir <- dirname(fileloc)
outfile <- file.path(outdir, "Score_vs_Energy.jpeg")

jpeg(outfile, width = 1280, height = 720)
df %>%
  ggplot(aes(Max_Score, Max_Energy_kcal_mol)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "blue", linewidth = 1.2) +
  theme_minimal()
dev.off()

