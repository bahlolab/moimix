# clean_metadata.R
library(readr)
library(dplyr)

mgen_all <- read_tsv("raw_data/pf3k_release_5_metadata.txt")

# extract just mixtures plus lab strains
mgen_mixtures_only <- mgen_all %>% 
    filter(contact_name == "Jason Wendler" | is.na(study_title)) %>% 
    select(sample, acc, bases, bases_mapped, bases_duplicated,
           mean_coverage, mean_fragment_size, sd_fragment_size, `%callable` )

# read in mixture metadata
mgen_mix <- read_tsv("raw_data/pf3k_release_5_mixtures_metadata.txt")

mgen_mixtures_all <- left_join(mgen_mixtures_only, mgen_mix,
                               by = c("sample", "acc"))

# fix headers sample labeling
mgen_mixtures_all$`7G8`[mgen_mixtures_all$sample == '7G8'] <- 100
mgen_mixtures_all$`3D7`[mgen_mixtures_all$sample == '7G8'] <- 0
mgen_mixtures_all$`Dd2`[mgen_mixtures_all$sample == '7G8'] <- 0
mgen_mixtures_all$`HB3`[mgen_mixtures_all$sample == '7G8'] <- 0


write_rds(mgen_mixtures_all, "processed_data/mixtures_all.rds")
