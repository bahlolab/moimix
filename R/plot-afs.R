# plot-afs.R
# Author: Stuart Lee
# Date: 16-02-2015
# Description: This script produces 152 plots of the B-allele frequencies
# by chromosome for each isolate in the Pfalciparum using the varscan 2
# calls (although it could be adpated to the other callers).

# get libraries
source("../lib/setup.R")

var.frq <- read.table(paste0(varscan.dir, "variant_frq.txt"),
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      comment.char = "",
                      na.strings  = c(".", "NA"))

stripPercent <- function(frq.df) {
  # strip percent sign from frequency data frame
  # Args:
  #   frq data.frame containing variant frequencies accross isolates
  # Returns:
  # A data frame containing numeric percentages for variant frequnecies
  rmPercent <- function(sample) {
    as.numeric(sub(pattern = "%", replacement = "", sample))
  }
  
  mutate_each(frq.df, funs(rmPercent), starts_with("PN"))
}

var.frq <- var.frq %>% stripPercent

chrom.names <- list("Pf3D7_01_v3" = "01",
                    "Pf3D7_02_v3" = "02",
                    "Pf3D7_03_v3" = "03",
                    "Pf3D7_04_v3" = "04",
                    "Pf3D7_05_v3" = "05",
                    "Pf3D7_06_v3" = "06",
                    "Pf3D7_07_v3" = "07",
                    "Pf3D7_08_v3" = "08",
                    "Pf3D7_09_v3" = "09",
                    "Pf3D7_10_v3" = "10",
                    "Pf3D7_11_v3" = "11",
                    "Pf3D7_12_v3" = "12",
                    "Pf3D7_13_v3" = "13",
                    "Pf3D7_14_v3" = "14")

chrom_labeller = function(variable, value) {
  return(chrom.names[value])
}

dir.create("../figures/afs_plots/")

plotAFS <- function(frq.df, sample.id) {
  # plot the B-allele frequency spectrum
  # Args: 
  #   frq.df: data.frame contain SNV positions and frequencies
  #   sample.id: isolate identifier
  # Returns:
  #   NULL but saves a png to the figures/afs_plot/ directory
  # find proportion of NC on each chromosome
  sample.df <- frq.df[, c("X.CHROM", "POS", sample.id)]
  #na.counts <- aggregate(sample.df[, sample.id],
  #                       by = list(sample.df[, "X.CHROM"]),
  #                       FUN = function(x) sum(is.na(x)) / length(x))
  #print(na.counts)
  
  p <- ggplot(data = sample.df, 
              aes_string(x = "POS", y = sample.id)) + 
    geom_jitter(alpha = I(1/20), size = 0.5) +
    facet_grid(X.CHROM ~ ., labeller = chrom_labeller) + 
    theme_bw() + 
    labs(y = "B-allele frequency",
         x = "Position (bp)") +
    theme(strip.text = element_text(size = 6),
          axis.text.y = element_text(size = 4),
          axis.text.x = element_text(size = 4),
          axis.title.x = element_text(size = 4),
          axis.title.y = element_text(size = 4))
  
  ggsave(plot = p, 
         filename = paste0("../figures/afs_plots/", sample.id, "-bAFS.png"),
         width = 16,
         height = 12,
         units = "cm")
}

# loop over each isolate and draw the plot
for (sample in colnames(var.frq)[108:ncol(var.frq)]) {
  # check that not all values are missing
  print(paste(sample, "is now plotting"))
  if (sum(is.na(var.frq[, sample])) / nrow(var.frq) >= 0.9) {
    print(paste(sample, "has only missing allele frequencies"))
    next
  }
  
  suppressWarnings(plotAFS(var.frq, sample))
}
