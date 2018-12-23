# set wd and load packages
library(vegan)
library(tidyverse)
library(data.table)
library(stringr)

#### Metadata ####
# metadata from Yassour et al. 2016. Natural history of the infant gut microbiome and impact of antibiotic treatment on bacterial strain diversity and stability. Science translational medicine 8:343ra81. [www.sciencetranslationalmedicine.org/cgi/content/full/8/343/343ra81/DC1 as supplementary table S1]

#note that the supplementary file is an excel file with three worksheets; below, each worksheet was saved as three individual csv files
d1_meta1 <- fread("data\\diab_yassour_meta_antibiotics.csv")
d1_meta2 <- fread("data\\diab_yassour_meta_feeding.csv")[, c(1:5)]
d1_meta3 <- fread("data\\diab_yassour_meta_general.csv")[, c(1:6)]

# metadata from Kostic et al. 2016. The Dynamiccs of the human gut microbiome in development and in progression towards type 1 diabetes 17:260-273. [https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(15)00021-9]
d2_meta <- fread("data\\diab_t1d_meta.csv")

#merge d1 metadata files
d1_meta <- d1_meta1 %>%
  transmute(subject = Subject, days = `Duration (days)`) %>%
  mutate(days = as.numeric(ifelse(days == 'Not known', 7, days))) %>%
  group_by(subject) %>%
  summarise(antibiotic_days = sum(days)) %>%
  full_join(d1_meta2 %>% rename(subject = Subject), by = c('subject')) %>%
  full_join(d1_meta3, by = 'subject') %>%
  transmute(
    subject, 
    antibiotic_days = ifelse(is.na(antibiotic_days), 0, antibiotic_days),
    bf_end = `Final age of breastfeeding (months)` * 30.5,
    country = Country,
    delivery = `Delivery type`,
    treatment_group = ifelse(antibiotic_days > 0, 'Antibiotics', 'Control'))

#clean up, convert days to months by dividing by 30.25
d2_meta <- d2_meta %>%
  group_by(Subject_ID) %>%
  mutate(n_abx = ifelse(sum(AbxAtCollection != 'no_abx') > 0, sum(AbxAtCollection != 'no_abx') * 7,
                 ifelse(TRUE %in% AbxExposureAbsolute, 1, 0))) %>%
  ungroup(Subject_ID) %>%
  transmute(
    sampleID = G_id, 
    subject = Subject_ID,
    antibiotic_days = n_abx,
    delivery = Delivery_Route,
    t = round(Age_at_Collection/30.25, 2),
    country = Country,
    treatment_group = ifelse(Case_Control == 'control', 'Control', 'T1D'))

#merge metadata
meta <- rbind(
  d1_meta %>% gather(var, val, -subject),
  d2_meta %>%
    select(-sampleID, -t) %>%
    distinct() %>%
    gather(var, val, -subject))

#clean metadata
meta <- meta %>%
  spread(var, val) %>%
  mutate(
    study = ifelse(subject %in% d1_meta$subject, 'Yassour2016', NA),
    study = ifelse(subject %in% d2_meta$subject, 'Kostic2015', study),
    delivery = ifelse(is.na(delivery), 'unknown', delivery),
    delivery = ifelse(delivery %in% c('vaginal','Vaginal delivery'), 'vaginal', 'caesaren'),
    antibiotic_days = as.numeric(antibiotic_days),
    bf_end = as.numeric(bf_end))

#### OTU data ####
# otu file generated using succession_seqs_to_otus.sh and usearch pipeline
otus <- read.table('data\\seqs_filtered_uniques_OTU_table.txt', header = TRUE, sep = '\t', row.names = 1, comment.char = '')
otus <- otus[, colSums(otus) > 5000]
set.seed(7)
otus <- as.data.frame(t(rrarefy(t(otus), 5000)))
otus <- otus %>%
  mutate(otu = row.names(otus)) %>%
  gather(sample, abun, -otu) %>%
  select(sample, otu, abun) %>%
  filter(grepl("time", sample) | sample %in% d2_meta$sampleID) %>%
  filter(abun > 0) %>%
  mutate(
    subject = sub("time.*", "", sample),
    subject = ifelse(subject %in% d1_meta$subject, subject, d2_meta$subject[match(subject, d2_meta$sampleID)]),
    t = ifelse(grepl("time", sample),
               sub("dot", ".", fixed = TRUE, sub(".*time", "", sample)),
               d2_meta$t[match(sample, d2_meta$sampleID)]),
    t = as.numeric(t)) %>%
  transmute(sampleID = paste(subject, t, sep = '_'), subject, t, otu, abun) %>%
  arrange(subject, t)

#remove subjects with fewer than 10 samples and/or less than 900 days sampling
otus <- otus %>%
  group_by(subject) %>%
  filter(length(unique(t)) >= 10) %>%
  filter(max(t) - min(t) > 30) %>%
  ungroup()

#### Taxonomic data ####
# SILVA sintax (taxonomy) file generated using succession_seqs_to_otus.sh and usearch pipeline
#Note: I also tried to use LTP to assign taxonomy. But it has WAY more 'unclassified'. Since we're just interested in getting taxon traits by mining the literature -- i.e., not on taxonomy/phylogeny per se -- we get way more hits if we use the whole SILVA database. Later, we map the OTUs onto the LTP phylogeny based on 16S similarity.
#I clean up some of the messier taxonomic designations
tax_succ_SILVA <- read.table('data\\seqs_filtered_uniques_otus_taxonomy_SILVA.sintax', fill = TRUE) %>%
  select(otu = V1, tax = V4) %>%
  separate(tax, c('Kingdom','Phylum','Class','Order','Family','Genus','Species'), sep = ',', fill = 'right') %>%
  mutate_all(funs(gsub('\\w:', '', .))) %>%
  mutate(
    Family = ifelse(Family %in% c('Clostridiaceae_1','Clostridium_sp._CA306'), 'Clostridiaceae', Family),
    Family = ifelse(grepl("_|[0-9]+", Family), NA, Family),
    Genus = trimws(Genus),
    Genus = ifelse(Genus == 'Incertae_Sedis', NA, Genus),
    Genus = gsub('\\[|\\]', '', Genus),
    Genus = gsub('-.*', '', Genus),
    Genus = gsub('Candidatus_', '', Genus),
    Genus = gsub('_subsp.*', '', Genus),
    Genus = gsub('_group', '', Genus),
    Genus = gsub('_[0-9]+', '', Genus),
    Species = ifelse(grepl('Eubacterium_|Sclerotinia_|Ruminococcus_', Genus), 
                     gsub('Eubacterium_|Sclerotinia_|Ruminococcus_', '', Genus), Species),
    Genus = gsub('_.*| .*', '', Genus),
    Genus = ifelse(Genus == 'Family'| grepl('.*eae', Genus), NA, Genus),
    Species = trimws(Species),
    Species = ifelse(grepl('[0-9]+', Species), NA, Species),
    Species = ifelse(grepl('sp.', Species, fixed = TRUE), NA, Species),
    Species = ifelse(grepl('_', Species), sub('.*_', '', Species), Species),
    Species = gsub('\\[|\\]', '', Species),
    Species = ifelse(Species %in% c('KT','LYH','unidentified'), NA, Species))

#load raw taxonomy data from LTP version 132, which is used for searching trait databases, 
#downloaded at [https://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_132/LTPs132_SSU.csv]
tax_LTP <- read.table('data\\LTPs132_SSU.csv', sep = '\t')
tax_LTP <- tax_LTP[!duplicated(tax_LTP$V1), ]
tax_LTP <- transmute(tax_LTP, otu = as.character(V1), 
                     Genus = unlist(lapply(strsplit(as.character(V5), ' ', fixed = TRUE), function(x) x[1])), 
                     Species = word(gsub("subsp.*", "", V5), 2))
tax_LTP$Genus[tax_LTP$otu == 'Y18189'] <- 'Clostridium'
tax_LTP$Species[tax_LTP$otu == 'Y18189'] <- 'polyendosporum'

#remove a c-section infant that also received heavy antibiotics
meta <- meta[meta$subject != 'E004628', ]
otus <- filter(otus, subject != 'E004628')

#remove metadata with no associated subjects or subjects with < 10 samples
meta <- filter(meta, subject %in% otus$subject)

#remove taxonomic data with no associated subjects or subjects with < 10 samples
tax_succ_SILVA <- filter(tax_succ_SILVA, otu %in% otus$otu)

#replace taxonomic NA with 'unclassified'
tax_succ_SILVA[is.na(tax_succ_SILVA)] <- 'unclassified'

#rename T1D treatment group as Control, because we are only looking at antibiotics as a predictor. T1D is one of many possible confounding factors, but we can't address all of them in this very-general study
meta$treatment_group[meta$treatment_group == 'T1D'] <- 'Control'

#fix one missing country entry
meta$country[meta$subject == 'E028794'] <- 'Finland'
meta$treatment_group[meta$subject == 'E028794'] <- 'Control'

#fix treatment group to exclude infants that received less than 50 days ab
meta <- mutate(meta,
  treatment_group = ifelse(delivery == 'caesaren', 'C-section',
    ifelse(!is.na(antibiotic_days) & antibiotic_days > 49, 'Antibiotics',
    ifelse(!is.na(antibiotic_days) & antibiotic_days == 0, 'Control', 'None'))))

#?& subject != 'E022137'? At some point I removed this infant because it's finnish, unlike all the others, but now... i'm leaving it in.


#### Save data ####
saveRDS(otus, file = 'data\\otus.RDS')
saveRDS(meta, file = 'data\\subject_metadata.RDS')
saveRDS(tax_succ_SILVA, 'data\\succession_tax_SILVA_all.RDS')
saveRDS(tax_LTP, 'data\\tax_LTP.RDS')

