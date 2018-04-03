# set wd and load packages
wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\'
setwd(wd)
library(tidyverse)
library(data.table)
library(stringr)
library(vegan)

#load full picrust/greengenes taxonomy
ggtax <- fread('bigdata_unsynced\\gg_13_5_taxonomy.txt', stringsAsFactors = FALSE)
ggtax <- ggtax %>%
  rename(otu = V1, tax = V2) %>%
  mutate(tax = gsub("\\s|.__|\\[|\\]|Other", "", tax)) %>%
  mutate(otu = paste0('otu',otu)) %>%
  separate(tax, sep = ';', c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), fill = 'right')
ggtax <- as.data.frame(apply(ggtax, 2, function(x) ifelse(x == '', 'unclassified', x)))
ggtax$otu <- as.character(ggtax$otu)

#Read data
d1_otus <- fread('data\\Yassour_otu_table_filtered.txt', skip = 1)
d2_otus <- fread("data\\diab_t1d_otus.txt", skip = 1)

#read metadata
d1_meta1 <- fread("data\\diab_yassour_meta_antibiotics.csv")
d1_meta2 <- fread("data\\diab_yassour_meta_feeding.csv")[, c(1:5)]
d1_meta3 <- fread("data\\diab_yassour_meta_general.csv")[, c(1:6)]
d2_meta <- fread("data\\diab_t1d_meta.csv")

# for each data source, clean otus, split taxonomy, replace 0 with NA
d1_otus <- d1_otus %>%
  mutate(otu = paste0('otu', `#OTU ID`)) %>%
  select(-`Consensus Lineage`, -`#OTU ID`)

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
  mutate(otu = paste0('otu', `#OTU ID`)) %>%
  select(-`ConsensusLineage`, -`#OTU ID`)
  
d2_meta <- d2_meta %>%
  transmute(
    sampleID = G_id, 
    subject = Subject_ID,
    delivery = Delivery_Route,
    t = Age_at_Collection,
    country = Country,
    treatment_group = ifelse(Case_Control == 'control', 'Control', 'T1D'))

#combine otus from the two studies
otus <- rbind(
  d1_otus %>% 
    gather(sampleID, abun, -otu, na.rm = T) %>%
    separate(sampleID, sep = '_', c("subject", "t"), remove = FALSE) %>%
    mutate(t = as.numeric(t) * 30.5) %>%
    transmute(sampleID, subject, t, otu, abun),
  d2_otus %>% 
    gather(sampleID, abun, -otu, na.rm = T) %>%
    left_join(d2_meta[, c('sampleID', 'subject', 't')], 
              by = c('sampleID')))

#filter and subsample at 12000
otus_wide <- otus %>%
  group_by(sampleID) %>%
  filter(sum(abun) >= 12000) %>%
  spread(otu, abun, fill = 0) %>%
  ungroup()

set.seed(7)
otus <- data.frame(otus_wide[, c(1:3)], rrarefy(otus_wide[, -c(1:3)], sample = 12000)) %>%
  gather(otu, abun, -sampleID, -subject, -t)

otus_wide <- spread(otus, otu, abun, fill = 0)

#ensure there are not empty columns
min(apply(otus_wide[,c(4:ncol(otus_wide))], 2, sum)) > 0

#create taxonomy file (for knowing which trait data we need to infer...)
tax <- data.frame(otu = names(otus_wide)[c(4:ncol(otus_wide))], stringsAsFactors = FALSE) %>%
  select(otu) %>%
  left_join(ggtax, by = 'otu')

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
    bf_end = as.numeric(bf_end)) %>%
  filter(subject %in% otus$subject) %>%
  group_by(subject) %>%
  ungroup()

#add column for number of time points per sample
tmp <- otus_wide %>%
  group_by(subject) %>%
  summarise(n = length(t), t_min = min(t), t_max = max(t), t_length = t_max - t_min)

meta <- left_join(meta, tmp[, c('subject','n','t_min', 't_max', 't_length')], by = 'subject')

#now, remove samples with fewer than 10 time points and at least 900 days of sampling
if (TRUE) {
  meta <- filter(meta, n > 10 & t_length > 900)
  
  otus_wide <- otus_wide %>%
    filter(subject %in% meta$subject) %>%
    gather(otu, abun, -sampleID, -subject, -t) %>%
    group_by(otu) %>%
    filter(sum(abun) > 0) %>%
    ungroup() %>%
    spread(otu, abun, fill = 0)
  
  tax <- filter(tax, otu %in% names(otus_wide)[c(4:ncol(otus_wide))])

}

#creating biom table for PICRUSt/bugbase
if (FALSE) {
  tmp <- otus_wide[, c(4:ncol(otus_wide))]
  x <- data.frame(otu = gsub("otu", "", colnames(tmp)), t(tmp))
  colnames(x) <- c('#OTU ID',otus_wide$sampleID)
  write.table(x, 'data\\otu_table.txt', sep = '\t', row.names = FALSE, quote = FALSE)
  system('biom convert -i C:/Users/John/Documents/msu/microbiome_trait_succession/data/otu_table.txt -o C:/Users/John/Documents/msu/microbiome_trait_succession/data/otu_table.biom --table-type="OTU table" --to-json')
  # then, upload the resulting biom file to https://bugbase.cs.umn.edu/upload.html and download results
}

if (TRUE) {
  x <- fread('data\\BugBase_Analyses\\normalized_otus\\16s_normalized_otus.txt')
  x <- as.data.frame(x)
  tmp <- paste0('otu', x$V1)
  x$V1 <- NULL
  x <- as.data.frame(t(x))
  colnames(x) <- tmp
  x <- mutate(x, sampleID = row.names(x))
  otus_wide_normalized <- otus_wide[, c('sampleID','subject','t')] %>%
    left_join(x, by = 'sampleID')
}

#save data
saveRDS(otus_wide, file = 'data\\otus_wide.RDS')
saveRDS(otus_wide_normalized, file = 'data\\otus_wide_normalized.RDS')
saveRDS(meta, file = 'data\\meta.RDS')
saveRDS(tax, file = 'data\\tax.RDS')
saveRDS(ggtax, file = 'data\\ggtax.RDS')

print("saved otu datafiles to msu/data/")

#cleanup
rm(list = ls())

