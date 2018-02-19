#Microbial data scubbing and organization... work in progress
#John Guittar
#last edit: Feb 16 2018

wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\'
setwd(wd)
source('custom_functions.R')
library(tidyverse)
library(data.table)

#load taxonomy
tax <- readRDS('data\\tax.RDS')

#load full picrust/greengenes taxonomy
pitax <- fread('C:\\Users\\John\\Documents\\msu\\analysis\\picrust-1.1.2\\picrust\\data\\gg_13_8_99.gg.tax', stringsAsFactors = FALSE)
pitax <- pitax %>%
  rename(otu = V1, tax = V2) %>%
  mutate(tax = gsub("\\s|.__|\\[|\\]|Other", "", tax)) %>%
  mutate(otu = paste0('otu',otu)) %>%
  separate(tax, sep = ';', c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","drop"), fill = 'right') %>%
  select(-drop)
pitax <- as.data.frame(apply(pitax, 2, function(x) ifelse(x == '', NA, x)))

####first: Edits from by Barberan et al 2016
####Then: my own edits

#read table
ijsem<-read.delim("data\\IJSEM_pheno_db_v1.0.txt", sep="\t", header=T, check.names=F, fill=T,
                  na.strings=c("NA", "", "Not indicated", " Not indicated","not indicated", "Not Indicated", "n/a", "N/A", "Na", "Not given", "not given","Not given for yeasts", "not indicated, available in the online version", "Not indicated for yeasts", "Not Stated", "Not described for yeasts", "Not determined", "Not determined for yeasts"))

#simplify column names
colnames(ijsem)<-c("Habitat", "Year", "DOI", "rRNA16S", "GC", "Oxygen",
                   "Length", "Width", "Motility", "Spore", "MetabAssays", "Genus", "Species", "Strain", "pH_optimum", "pH_range", "Temp_optimum", "Temp_range", "Salt_optimum", "Salt_range", "Pigment", "Shape", "Aggregation", "FirstPage", "CultureCollection", "CarbonSubstrate", "Genome", "Gram", "Subhabitat", "Biolog")

#clean Habitat column
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater (river, lake, pond)"]<-"freshwater"
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater sediment (river, lake, pond"]<-"freshwater sediment"

#clean Oxygen column
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="aerobic"]<-"obligate aerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="anerobic"]<-"obligate anerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="microerophile"]<-"microaerophile"

#clean pH_optimum column
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$pH_optimum<-as.character(ijsem$pH_optimum)
ijsem$pH_optimum<-sapply(ijsem$pH_optimum, simplify=T, function(x){log10(10^mean(swan(unlist(strsplit(x, split="-", fixed=T)))))})

#remove pH values <0 and >10
ijsem$pH_optimum[ijsem$pH_optimum<0 | ijsem$pH_optimum>10]<-NA

#clean Temp_optimum column
ijsem$Temp_optimum<-as.character(ijsem$Temp_optimum)

#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$Temp_optimum<-sapply(ijsem$Temp_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))))})

#clean Salt_optimum column
ijsem$Salt_optimum<-as.character(ijsem$Salt_optimum)
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs

ijsem$Salt_optimum<-sapply(ijsem$Salt_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))))})
#there are some formatting issues that should be solved


##############Below are my additional edits

#Assign oxygen score
ijsem$Oxygen_score <- c(5,4,3,2,1)[match(ijsem$Oxygen, 
                                         c('obligate aerobe','microaerophile','facultative anaerobe','facultative aerobe','anaerobic'))]

#trun Aggregation into a binary
ijsem$Aggregation <- c(0,1,1)[match(ijsem$Aggregation, c('none','chain','clump'))]

#clean Length and Width
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
x <- as.character(ijsem$Length)
x <- gsub(' in diameter| \\(diameter\\)|mm|Indicated|diameter: | ', '', x)
x <- gsub(',', '.', x, fixed = TRUE)
x[grepl(">", x)] <- NA
x[x == ''] <- NA
x[x == c('0.61.6')] <- '0.6-1.6'
x[x == c('0.4 0.9 ')] <- '0.4-0.9'
x <- gsub('дус *', '-', x)

x <- ifelse(grepl("-[a-zA-Z]+", x) & !is.na(x), 
            as.character(format(as.Date(x, "%d-%b"), "%d.%m")), 
            x)
x <- ifelse(grepl("[a-zA-Z]+-", x) & !is.na(x), 
            as.character(format(as.Date(paste(1, x), "%d %b-%y"), "%m.%y")), 
            x)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", x) & !is.na(x)
x[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[2]]))))

#if it has a backsclash, assume it was a decimal point
filt <- grepl("\\/", x) & !is.na(x)
x[filt] <- paste(sapply(strsplit(x[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(x[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Length <- swan(x)

# remove the inordinately large values...
ijsem$Length[ijsem$Length > 250] <- NA
ijsem$Length[ijsem$Length == 0] <- NA

######################

x <- as.character(ijsem$Width)
x <- gsub(' |\\(|\\)|not.+|indiameter|Filament.+', '', x)
x <- gsub(',', '.', x, fixed = TRUE)
x[grepl("<", x)] <- NA
x[x == ''] <- NA
x <- gsub('дус *', '-', x)
x <- ifelse(grepl("-[a-zA-Z]+", x) & !is.na(x), 
            as.character(format(as.Date(x, "%d-%b"), "%d.%m")), 
            x)
x <- ifelse(grepl("[a-zA-Z]+-", x) & !is.na(x), 
            as.character(format(as.Date(paste(1, x), "%d %b-%y"), "%m.%y")), 
            x)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", x) & !is.na(x)
x[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[2]]))))

#if it has a backsclash, assume it was a decimal point
filt <- grepl("\\/", x) & !is.na(x)
x[filt] <- paste(sapply(strsplit(x[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(x[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Width <- as.numeric(x)
ijsem$Width[ijsem$Width > 250] <- NA
ijsem$Width[ijsem$Width == 0] <- NA

#############################################

#Turn motility into a binary
ijsem$Motility <- c(0,1,1,1,1)[match(ijsem$Motility, c('non-motile', 'flagella',
                                                       'motile, but unspecified structure', 'gliding', 'axial filament'))]

#turn spore into binary
ijsem$Spore <- c(0,1)[match(ijsem$Spore, c('no','yes'))]

#turn gram status into a binary
ijsem$Gram <- c(0, 0.5, 1)[match(ijsem$Gram, c('negative','variable','positive'))]

# fix pH
# NOTE: THE WAY THAT THEY CALCULATED PH WHEN THERE IS A RANGE - BY TAKING THE MEAN -- IS INCORRECT
# MUST UNLOG PH BEFORE CALCULATING MEAN, AND THEN RE-LOG.
# I HAVEN'T DONT THAT YET HERE... 

x <- as.character(ijsem$pH_range)
x <- gsub(' |\\(|\\)|at 6\\.0 but not 5\\.5|>|<|not inidcated', '', x)
x <- gsub('MinimumpHof5|t|4\\.0\\+|5\\.0\\+;notat4\\.5|14\\.08', '', x)
x[x %in% c('13.05','4.3','7.1','7.2','5.5','7.5','30.04','5.1')] <- NA
x[x == '5.7and6.8'] <- '5.7-6.8'
x[x == '4.0and8.5'] <- '4-8.5'
x[x == '5,6,10,12,'] <- '5-12'
x[x == '5and9.5'] <- '5-9.5'
x[x == '5.59.5'] <- '5.5-9.5'

# remove anything with an inequality
x[grepl("<|>", x)] <- NA
x[x == ''] <- NA
x <- gsub('дус *', '-', x)
x <- ifelse(grepl("-[a-zA-Z]+", x) & !is.na(x), 
            as.character(format(as.Date(x, "%d-%b"), "%d.%m")), 
            x)

filt <- grepl("\\/", x) & !is.na(x)
x[filt] <- paste(sapply(strsplit(x[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(x[filt], '\\/'), function(x) x[[2]]), sep = '.')

x[!grepl("-", x) & !is.na(x) & !grepl("\\.", x)] <- NA
filt <- !grepl("-", x) & !is.na(x) & grepl("\\.", x)
x[filt] <- gsub("\\.", "-", x[filt])
lows <- swan(sapply(strsplit(x, '-'), function(x) x[[1]]))
highs <- swan(sapply(strsplit(x, '-'), function(x) if (length(x) > 1) x[[2]] else NA))

ijsem$pH_lows <- lows
ijsem$pH_highs <- highs

####GC
x <- as.character(ijsem$GC)
x <- gsub('дус *', '-', x)
x <- gsub(',', '\\.', x)
x[x %in% c('DQ987877','marine organism','BAOL01000001','40.50%','GU323338')] <- NA
x <- gsub('O±.+|\\+/-.+|%|o', '', x)
filt <- grepl("-", x) & !is.na(x)
x[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[2]]))))

ijsem$GC <- swan(x)

# at the moment, i only want these traits...
ijsem <- select(ijsem, Genus, Species, GC_content = GC, 
                Oxygen_tolerance = Oxygen_score, Length, Width, 
                Motility = Motility, Spore, 
                pH_optimum, Temp_optimum, Salt_optimum, Aggregation_score = Aggregation, 
                Gram_positive = Gram) %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val)) %>%
  filter(!Genus %in% c("10.1099/ijs.0.006320-0","10.1099/ijs.0.02424-0")) %>%
  mutate_if(is.factor, as.character)

###########################################################

spo <- read.csv('data/Browne2016_sporulationTable.csv', header = T, skip = 1)

# fix long colnames
colnames(spo) <- substr(colnames(spo), 1, 2)

spo <- spo %>%
  transmute(
    Genus = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[1]])),
    Species = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[2]])),
    trait = 'Spore_score',
    val = si)

#############

# load 16s gene copy number data, derived from PICRUSt
gcn <- read.table('C:\\Users\\John\\Documents\\msu\\analysis\\picrust-1.1.2\\picrust\\data\\16S_13_5_precalculated.tab', stringsAsFactors = FALSE)
gcn <- gcn %>%
  transmute(
    otu = paste0('otu', V1), 
    trait = 'Copies_16S', 
    val = V2)

#add full 16S copy number from picrust to picrust taxonomy
#I will remove the subset of 16S data later
pitax <- pitax %>%
  mutate(otu = as.character(otu)) %>%
  left_join(gcn, by = 'otu') %>% 
  rename(Copies_16S = val) %>% 
  select(-trait)

gcn <- gcn %>%
  inner_join(tax[, c('otu','Genus','Species')], by = 'otu')

####################################### Exploring BacDrive

#library(RCurl)
#library(rjson)

#needed <- sort(as.character(unique(filter(tax, is.na(Gram_positive))$Genus)))
needed <- sort(unique(tax$Genus))
needed <- needed[!needed %in% c("unclassified","SMB53","02d06","rc4-4","WAL_1855D","1-68","cc_115")]

if (FALSE) {
  
  bacdat <- bac_search(needed) #re-calculate bacdata
  
} else {
  
  load('data\\bacdrive_genera.Rdata')#load bacdata
  
}

#taxonomic info
tax_bacdat <- sapply(bacdat, function(x)
  sapply(x, function(y)
    sapply(y$taxonomy_name$strains_tax_PNU, function(sp) sp$species_epithet)))

#oxygen tolerance
o2 <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$oxygen_tolerance, function(val) val$oxygen_tol)))
o2 <- bac_extract(o2, 'Oxygen_tolerance')
o2$val <- c(5,4,3,1,1)[match(o2$val, 
                             c('aerobe','microaerophile','facultative anaerobe','obligate anaerobe','anaerobe'))]

#gram stain
gs <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$gram_stain)))
gs <- bac_extract(gs, 'Gram_positive')
gs$val <- ifelse(gs$val == 'positive', 1, 0)

#spore forming
spofo <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$spore_formation, function(val) val$ability)))
spofo <- bac_extract(spofo, "Spore")
spofo$val <- ifelse(spofo$val == 'TRUE', 1, 0)

#motility
motility <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$motility)))
motility <- bac_extract(motility, 'Motility')
motility$val <- ifelse(motility$val == TRUE, 1, 0)

#cell length
cell_length <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$cell_len)))
cell_length <- bac_extract(cell_length, 'Length')
cell_length$val <- ifelse(grepl('<|>', cell_length$val), NA,
                          ifelse(grepl('-', cell_length$val), sapply(strsplit(cell_length$val, '-'), function(y) 
                            mean(swan(as.character(unlist(y))))), swan(cell_length$val)))
cell_length <- cell_length[!is.na(cell_length$val), ]

#cell width
cell_width <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$cell_width)))
cell_width <- bac_extract(cell_width, 'Width')
cell_width$val <- ifelse(grepl('-', cell_width$val), 
                         sapply(strsplit(cell_width$val, '-'), function(y) 
                           mean(swan(as.character(unlist(y))))), 
                         swan(as.numeric(cell_width$val)))

#aggregation (\any data at all? Am I doing it right?)
aggreg <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$multicellular_morphology, function(val) val$ability)))
aggreg <- bac_extract(aggreg, 'Aggregation_score')
aggreg$val <- ifelse(aggreg$val == TRUE, 1, 0)

bd <- rbind(o2, gs, spofo, motility, cell_length, cell_length, aggreg)
bd <- mutate_if(bd, is.factor, as.character)

###############################

## Genome size, GC content, and Gene Number from NCBI...
genos <- fread('data\\NCBI_prokaryotes.txt') %>% 
  filter(Status == 'Complete Genome') %>%
  transmute(
    sp = gsub("'|\\[|\\]", "", `Organism/Name`),
    Genus = sapply(sapply(sp, strsplit, ' '), function(x) x[[1]]),
    Species = sapply(sapply(sp, strsplit, ' '), function(x) x[[2]]),
    Genome_Mb = as.numeric(ifelse(`Size (Mb)` == '-', NA, `Size (Mb)`)),
    GC_content = `GC%`, 
    Gene_number = as.numeric(ifelse(Genes == '-', NA, Genes))) %>%
  gather(trait, val, Genome_Mb, GC_content, Gene_number) %>%
  select(-sp) %>%
  filter(!is.na(val))

####################
# 16s gene copy numbers from rrnDB

rrnDB <- read.csv('data//rrnDB-5.3.csv') %>%
  mutate(
    sp = gsub("'|\\[|\\]", "", NCBI.scientific.name),
    Genus = sapply(sapply(sp, strsplit, ' '), function(x) x[[1]]),
    Species = sapply(sapply(sp, strsplit, ' '), function(x) ifelse(length(x) > 1, x[[2]], NA)),
    trait = 'Copies_16S') %>%
  transmute(Genus, Species, trait, val = X16S.gene.count) %>%
  filter(!is.na(val))

################ IGA???
iga <- read.csv('data\\Palm2014_IGA.csv', stringsAsFactors = FALSE)
iga <- t(iga)
colnames(iga) <- iga[1, ]
iga <- as.data.frame(iga[c(2:nrow(iga)), ])
iga[, c(4:ncol(iga))] <- sapply(iga[, c(4:ncol(iga))], function(x) as.numeric(as.character(x)))

iga <- iga %>%
  gather(tax, val, -var, -status, -subj) %>%
  spread(var, val) %>%
  mutate(tax = gsub("\\s|.__|\\[|\\]|Other", "", tax)) %>%
  separate(tax, sep = ';', c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  filter(!(Genus == 'unclassified' | Genus == '')) %>%
  mutate(Species = ifelse(Species == '', 'unclassified', Species)) %>%
  group_by(Genus, Species, trait = 'IgA') %>%
  summarise(val = log(mean(ici) + 1))

#################
#b-vitamins?
bvit <- read.csv('data\\Bvitamins_Magnusdottir2015.csv', stringsAsFactors = FALSE) %>%
  mutate(Genus = sapply(sapply(tax, strsplit, ' '), function(x) x[[1]]),
         Species = sapply(sapply(tax, strsplit, ' '), function(x) x[[2]])) %>%
  mutate(Species = ifelse(Species == 'sp.', 'unclassified', Species)) %>%
  select(-NCBI.Taxonomy.ID, -tax, -Body.Site) %>%
  gather(Bvit, val, -Genus, -Species) %>%
  group_by(Genus, Species, trait = 'B_vitamins') %>%
  summarise(val = length(unique(Bvit[val > 0])))

# what about GOLD - JGI?


# put it all together
x <- bind_rows(
  mutate(ijsem, source = 'IJSEM'),
  mutate(spo, source = 'Browne2016'),
  mutate(bd, source = 'BacDrive'),
  transmute(gcn, Genus, Species, trait, val, source = 'PICRUSt'),
  mutate(genos, source = 'NCBI'),
  mutate(rrnDB, source = 'rrnDB'),
  mutate(iga, source = 'Palm2014'),
  mutate(bvit, source = 'Mag2015')
)

uncleaned <- nrow(x)
### I looked at all the points plotted using 
#ggplot(x, aes(x = 1, y = val)) + geom_jitter() + facet_wrap(~trait, scales = 'free')
#ggplot(x, aes(x = val)) + geom_histogram(bins = 100) + facet_wrap(~trait, scales = 'free') + scale_y_sqrt()
# identifying and removing outliers
x <- x %>%
  filter(!(trait == 'Copies_16S' & val > 20)) %>%
  filter(!(trait == 'GC_content' & (val > 80 | val < 20))) %>%
  filter(!(trait == 'Gene_number' & val > 11000)) %>%
  filter(!(trait == 'Genome_Mb' & val > 14)) %>%
  filter(!(trait == 'Length' & val > 30)) %>%
  filter(!(trait == 'pH_optimum' & val < 2.5)) %>%
  filter(!(trait == 'Salt_optimum' & val > 25)) %>%
  filter(!(trait == 'Temp_optimum' & val > 80)) %>%
  filter(!(trait == 'Width' & val > 8))

print(paste(uncleaned - nrow(x), "of", nrow(x), "outliers were removed"))

#for plotting source-wise comparisons
x1 <- x

x <- x %>%
  group_by(Genus, Species, trait) %>%
  summarise(val = mean(val)) %>%
  mutate(val = ifelse(trait %in% c('Length','Width'), log(val), val))

x_all <- spread(x, trait, val) %>%
  filter(Genus != 'unclassified')

x_sparse <- x %>%
  spread(trait, val) %>%
  filter(Species != 'unclassified') %>%
  filter(Genus %in% tax$Genus & Species %in% tax$Species)

if (FALSE) {
  ## plotting coverages 
  x1 <- x1 %>%
    filter(!Genus %in% c('02d06','1-68')) %>%
    group_by(Genus, Species, trait, source) %>%
    summarise(val = mean(val)) %>%
    spread(source, val)
  
  x1 %>%
    filter(!is.na(IJSEM) & !is.na(NCBI)) %>%
    ggplot(aes(x = IJSEM, y = NCBI)) +
    geom_point() +
    facet_wrap(~trait, scales = 'free')
  
  x1 %>%
    filter(!is.na(IJSEM) & !is.na(BacDrive)) %>%
    ggplot(aes(x = IJSEM, y = BacDrive)) +
    stat_smooth() +
    geom_point() +
    facet_wrap(~trait, scales = 'free')
  
  x1 %>%
    filter(!is.na(PICRUSt) & !is.na(rrnDB)) %>%
    ggplot(aes(x = PICRUSt, y = rrnDB)) +
    geom_point() +
    geom_abline(slope = 1) +
    facet_wrap(~trait, scales = 'free')
  
  sims <- x %>%
    group_by(Genus, trait) %>%
    filter(length(unique(Species)) > 2 | (trait == 'IgA' & length(unique(Species)) > 1)) %>%
    summarise(
      n = length(Genus),
      mean = mean(val),
      var = var(val)) %>%
    group_by(trait, n) %>%
    mutate(var_random = mean(replicate(100, var(sample(x$val[x$trait == trait[1]], size = n)))))
  
  sims <- filter(sims, trait != 'Spore')
  
  sims %>%
    group_by(trait) %>%
    filter(Genus %in% tax$Genus | trait == 'Spore_score') %>%
    mutate(
      rank = as.numeric(as.factor(jitter(mean))),
      sd = sqrt(var)) %>%
    ggplot(aes(x = mean, y = rank, color = trait)) +
    geom_point() +
    geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) +
    facet_wrap(~trait, scales = 'free') +
    theme(legend.position = 'none')
  
  sims %>% 
    group_by(trait) %>%
    arrange(mean) %>%
    mutate(max = max(var, var_random)) %>%
    ggplot(aes(x = var, y = var_random, color = trait)) +
    geom_point(aes(size = n)) +
    geom_point(aes(x = 0, y = 0), alpha = 0) +
    geom_point(aes(x = max, y = max), alpha = 0) +
    geom_abline(slope = 1, lty = 2) +
    facet_wrap(~trait, scales = 'free')
  
  sims %>%
    ungroup() %>%
    filter(Genus %in% tax$Genus) %>%
    mutate(trait = ifelse(trait == 'Spore_score', 'Sporulation', trait)) %>%
    group_by(trait) %>%
    #filter(!(trait == 'Length' & var_random - var < -2)) %>%
    #filter(!(trait == 'Salt_optimum' & var_random - var < -5)) %>%
    #filter(!(trait == 'pH_optimum' & var_random - var < -0.5)) %>%
    #filter(!(trait == 'Width' & var_random - var < -0.5)) %>%
    mutate(
      diff = var_random - var,
      max = max(abs(diff))) %>%
    ggplot(aes(x = diff, fill = trait)) +
    geom_density() +
    geom_point(aes(x = -max, y = 0), alpha = 0) +
    geom_point(aes(x = max, y = 0), alpha = 0) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_wrap(~trait, scales = 'free') +
    theme(legend.position = 'none',
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = 'Null variance - Observed variance',
         y = 'Density of genera')
  
  sims %>%
    group_by(trait) %>%
    mutate(
      diff = var_random - var,
      max = max(abs(diff)),
      in_study = factor(as.character(Genus %in% tax$Genus), levels = c(TRUE, FALSE))) %>%
    #filter(in_study == TRUE) %>%
    ggplot(aes(x = diff, fill = trait)) +
    geom_histogram(aes(y =..density.., color = trait), bins = 30, alpha = 0) + 
    geom_density(alpha = 0.7, lty = 2, color = NA) +
    geom_point(aes(x = -max, y = 0), alpha = 0) +
    geom_point(aes(x = max, y = 0), alpha = 0) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_wrap(~trait, scales = 'free') +
    labs(x = 'Null variance - Observed variance', y = 'Density of genera') +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank())
  
  tax_genus <- x %>%
         filter(Genus %in% tax$Genus) %>%
         filter(Genus != 'unclassified') %>%
         group_by(Genus, trait) %>%
         summarise(val_genus = mean(val)) %>%
         filter(!is.na(val_genus)) %>%
         ungroup()
    
    tax_traits <- x %>%
      group_by(Genus, Species, trait) %>%
      summarise(val = mean(val)) %>%
      ungroup() %>%
      spread(trait, val) %>%
      right_join(tax[, c('Genus','Species', 'otu')], by = c('Genus','Species')) %>%
      gather(trait, val_species, -Genus, -Species, -otu) %>%
      mutate(
        val_genus = tax_genus$val_genus[match(paste(Genus, trait), paste(tax_genus$Genus, tax_genus$trait))],
        data_level = ifelse(!is.na(val_species), 'Species', ifelse(!is.na(val_genus), 'Genus', 'None')),
        val = ifelse(data_level == 'Species', val_species, ifelse(data_level == 'Genus', val_genus, NA))) %>%
      select(-val_species, -val_genus) %>%
      mutate_if(is.character, as.factor)

    tax_traits %>%
      filter(trait != 'Spore') %>%
      filter(data_level != 'None') %>%
      group_by(Genus, trait) %>%
      summarise(
        Measured = -length(unique(Species[data_level == 'Species'])),
        Inferred = length(unique(Species[data_level == 'Genus']))) %>%
      gather(Source, val, Measured, Inferred) %>%
      ungroup() %>%
      mutate(Genus = as.numeric(as.factor(Genus))) %>%
      ggplot(aes(x = Genus, y = val, fill = Source)) +
        geom_bar(stat = 'identity') +
        facet_wrap(~trait, ncol = 4)

  }
###################################################################################

sims <- x %>%
  group_by(Genus, trait) %>%
  filter(length(unique(Species)) > 2 | (trait == 'IgA' & length(unique(Species)) > 1)) %>%
  summarise(
    n = length(Genus),
    mean = mean(val),
    var = var(val)) %>%
  group_by(trait, n) %>%
  mutate(var_random = mean(replicate(100, var(sample(x$val[x$trait == trait[1]], size = n)))))

# calculate our best idea of trait based on taxonomic info
# species when possible, otherwise genus
x <- x %>%
  ungroup() %>%
  filter(Genus %in% tax$Genus) %>%
  mutate(
    Species = ifelse(Species == 'unclassified', '', Species),
    data_level = ifelse(paste(Genus, Species) %in% paste(tax$Genus, tax$Species), 'Species','Genus'),
    Species = ifelse(data_level == 'Genus', 'unclassified', Species)) %>%
  group_by(trait, Genus) %>%
  left_join(sims[, c('Genus','var','var_random')], by = "Genus") %>%
  mutate(val = ifelse(data_level == 'Genus' & var < var_random, mean(val), val)) %>%
  group_by(Genus, Species, trait) %>%
  summarise(val = mean(val)) %>%
  spread(trait, val)

tax <- tax %>%
  left_join(x, by = c('Genus','Species'))

# Add a few last minute manual entries
tax$Aggregation_score[tax$Genus == 'Bifidobacterium'] <- 1 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3346390/

#Merge spore scores (if no score and ijsem says 0, then 0; 
# otherwise the median of spore score when we know for ijsem spore == 1

tax <- tax %>% mutate(
  Sporulation = ifelse(is.na(Spore_score), 
    ifelse(Spore < 0.5, 0, median(tax$Spore_score[tax$Spore > 0], na.rm = T)), Spore_score))

traits <- select(tax, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species, -Spore, -Spore_score)

x_all <- left_join(pitax, x_all[, names(x_all) != 'Copies_16S'], by = c("Genus", "Species"))
write.csv(x_all, file = 'data\\traits_all.csv', row.names = FALSE)
write.csv(x_sparse, file = 'data\\traits_sparse.csv', row.names = FALSE)
write.csv(traits, file = 'data\\traits.csv', row.names = FALSE)
print("saved traits.csv, traits_sparse.csv, and traits_all.csv to msu/data/")

rm(list = ls())

#########PROTRAIT DATA STUFF
#prot <- fread('data\\ProTraits_binaryIntegratedPr0.90.txt')
#prot[prot == '?'] <- NA
#prot[, c(2:ncol(prot))] <- apply(prot[,c(2:ncol(prot))], 2, as.numeric)
#prot <- as.data.frame(prot)
#prot <- gather(prot, substrate, use, -Organism_name, -Tax_ID)
#prot <- filter(prot, !is.na(use))
#OTU ID doesn't work: table(tax$otu %in% paste0('otu',prot$Tax_ID))
# Genus species is a little better:
#table(paste(tax$Genus, tax$Species) %in% prot$Organism_name)