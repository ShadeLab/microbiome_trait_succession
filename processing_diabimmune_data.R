# set wd and load packages
wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\'
setwd(wd)
library(tidyverse)
library(data.table)
library(stringr)
library(vegan)

#Read data
d1_otus <- fread('data\\Yassour_otu_table_filtered.txt', skip = 1)
d2_otus <- fread("data\\diab_t1d_otus.txt", skip = 1)
d3_otus <- readRDS("data\\diab_karelia_otus.RDS")

#read metadata (tmp1:3 are for d1_meta)
d1_meta1 <- fread("data\\diab_yassour_meta_antibiotics.csv")
d1_meta2 <- fread("data\\diab_yassour_meta_feeding.csv")[, c(1:5)]
d1_meta3 <- fread("data\\diab_yassour_meta_general.csv")[, c(1:6)]
d2_meta <- fread("data\\diab_t1d_meta.csv")
d3_meta <- readRDS("data\\diab_karelia_meta.RDS")

#list of taxonomic tiers used throughout
taxvars <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# for each data source, clean otus, split taxonomy, replace 0 with NA
d1_otus <- d1_otus %>%
  rename(tax = `Consensus Lineage`, otu = `#OTU ID`) %>%
  mutate(tax = gsub("\\s|.__", "", tax)) %>%
  separate(tax, sep = ';', taxvars) %>%
  mutate(otu = paste0('otu', otu)) %>%
  mutate_all(funs(ifelse(. == 0, NA, .)))

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

d2_otus <- d2_otus %>%
  mutate(
    otu = paste0('otu', `#OTU ID`),
    tax = gsub("\\s|.__", "", ConsensusLineage)) %>%
  select(-`#OTU ID`, -ConsensusLineage) %>%
  separate(tax, sep = ';', taxvars) %>%
  mutate_at(.vars = vars(Domain, Phylum, Class, Order, Family, Genus, Species),
            .funs = funs(ifelse(. == "", "unclassified", .))) %>%
  mutate_all(funs(ifelse(. == 0, NA, .)))

d2_meta <- d2_meta %>%
  transmute(
    sampleID = G_id, 
    subject = Subject_ID,
    delivery = Delivery_Route,
    t = Age_at_Collection,
    country = Country,
    treatment_group = ifelse(Case_Control == 'control', 'Control', 'T1D'))

d3_otus <- d3_otus %>%
  filter(str_count(sample, '\\|') == 7) %>%
  mutate(sample = gsub("\\s|.__", "", sample)) %>%
  separate(sample, sep = '\\|', c(taxvars,"otu")) %>%
  mutate(otu = paste0('otu', otu)) %>%
  mutate_all(funs(ifelse(. == 0, NA, .)))

d3_meta <- d3_meta %>%
  transmute(
    subject = subjectID,
    sampleID = SampleID,
    t = age_at_collection,
    delivery,
    country = c('Finland','Estonia','Russia')[match(country, c('FIN','EST','RUS'))],
    bf_end = bf_length * 30.5,
    formula_end = Any_baby_formula,
    solid_food = Other_than_BF,
    antibiotic_days = num_abx_treatments * 7) %>%
  group_by(subject) %>%
  mutate(
    formula_end = ifelse(sum(formula_end) > 0, max(t[formula_end]), NA),
    solid_food = ifelse(sum(solid_food) > 0, min(t[solid_food]), NA),
    antibiotic_days = max(antibiotic_days),
    treatment_group = 'Control') %>%
  as.data.frame()

#create tax
tax <- rbind(
  d1_otus[, c(taxvars, 'otu')],	
  d2_otus[, c(taxvars, 'otu')],	
  d3_otus[, c(taxvars, 'otu')]) %>%
  distinct(otu, .keep_all = TRUE)

#replace brackets and fill spaces
tax <- sapply(tax, function(x) gsub('\\[|\\]', '', as.character(x)))
tax <- as.data.frame(apply(tax, 2, function(x) ifelse(x == '', 'unclassified', x)), stringsAsFactors = FALSE)

#fix nomenclature issue
tax$Genus[tax$Genus == "Pseudoramibacter_Eubacterium"] <- "Pseudoramibacter"

otus <- rbind(
  d1_otus[, !colnames(d1_otus) %in% taxvars] %>% 
    gather(sampleID, abun, -otu, na.rm = T) %>%
    separate(sampleID, sep = '_', c("subject", "t"), remove = FALSE) %>%
    mutate(t = as.numeric(t) * 30.5) %>%
    transmute(sampleID, subject, t, otu, abun),
  d2_otus[, !colnames(d2_otus) %in% taxvars] %>% 
    gather(sampleID, abun, -otu, na.rm = T) %>%
    left_join(d2_meta[, c('sampleID', 'subject', 't')], 
              by = c('sampleID')),
  d3_otus[, !colnames(d3_otus) %in% taxvars] %>% 
    gather(sampleID, abun, -otu, na.rm = T) %>%
    left_join(d3_meta[, c('sampleID', 'subject', 't')],
              by = c('sampleID'))
)

#filter and subsample at 12000
otus_wide <- otus %>%
  group_by(sampleID) %>%
  filter(sum(abun) >= 12000) %>%
  spread(otu, abun, fill = 0) %>%
  ungroup()

otus <- data.frame(otus_wide[, c(1:3)], rrarefy(otus_wide[, -c(1:3)], sample = 12000)) %>%
  gather(otu, abun, -sampleID, -subject, -t)

otus_wide <- spread(otus, otu, abun, fill = 0)

meta <- rbind(
  d1_meta %>% gather(var, val, -subject),
  d2_meta %>%
    select(-sampleID, -t) %>%
    distinct() %>%
    gather(var, val, -subject),
  d3_meta %>%
    select(-sampleID, -t) %>%
    distinct() %>%
    gather(var, val, -subject))

meta <- meta %>%
  spread(var, val) %>%
  mutate(
    study = ifelse(subject %in% d1_meta$subject, 'Yassour2016', NA),
    study = ifelse(subject %in% d2_meta$subject, 'Kostic2015', study),
    study = ifelse(subject %in% d3_meta$subject, 'Kostic2016', study),
    delivery = ifelse(is.na(delivery), 'unknown', delivery),
    delivery = ifelse(delivery %in% c('vaginal','Vaginal delivery'), 'vaginal', 'caesaren'),
    antibiotic_days = as.numeric(antibiotic_days),
    bf_end = as.numeric(bf_end),
    formula_end = as.numeric(formula_end),
    food_start = as.numeric(solid_food)) %>%
  filter(subject %in% otus$subject) %>%
  group_by(subject) %>%
  mutate(milk_end = ifelse(!is.na(bf_end)|!is.na(formula_end), max(bf_end, formula_end, na.rm = T), NA)) %>%
  ungroup()

#add column for number of time points per sample
tmp <- otus_wide %>%
  group_by(subject) %>%
  summarise(n = length(t), t_min = min(t), t_max = max(t), t_length = t_max - t_min)

meta <- left_join(meta, tmp[, c('subject','n','t_min', 't_max', 't_length')], by = 'subject')

saveRDS(otus_wide, file = 'data\\otus_wide.RDS')
saveRDS(meta, file = 'data\\meta.RDS')
saveRDS(tax, file = 'data\\tax.RDS')

print("saved otus_wide.RDS, meta.RDS, and tax.RDS to msu/data/")

rm(list = ls())

