#Succession gut microbiome traits -- data gathering and cleaning
#John Guittar

##Overview
#1. load raw taxonomic data from usearch output
#2. load taxonomy from LTP (Living Tree Project)
#3. remove any ambiguities/special characters
#4. gather and combine trait data from various sources
#5. match trait data to as many taxa as possible, either from usearch or LTP
#6. save raw trait data (no phylogenetic inference yet) 

######################################################################
######################################################################
######################################################################

#set working directory, load packages
wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\data'
setwd(wd)
source('C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\succession_custom_functions.R')
loadpax(c('data.table','tidyverse'))

#load taxonomy data, which will be used to filter otu unnecessary taxa
tax_succ <- readRDS('succession_tax_SILVA.RDS')
tax_LTP <- readRDS('tax_LTP.RDS')

### Begin amassing trait data...
####First: Edits drawn directly from Barberan et al 2016 script: 
#https://figshare.com/articles/International_Journal_of_Systematic_and_Evolutionary_Microbiology_IJSEM_phenotypic_database/4272392

#read ijsem table
ijsem<-read.delim('IJSEM_pheno_db_v1.0.txt', sep="\t", header=T, check.names=F, fill=T,
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
ijsem$pH_optimum<-sapply(ijsem$pH_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))), na.rm = T)})

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


############## Now my additional edits to the IJSEM data

#Assign oxygen score
ijsem$Oxygen_score <- c(5,4,3,2,1)[match(ijsem$Oxygen, 
                                         c('obligate aerobe','microaerophile','facultative anaerobe','facultative aerobe','anaerobic'))]

#trun Aggregation into a binary
ijsem$Aggregation <- c(0,1,1)[match(ijsem$Aggregation, c('none','chain','clump'))]

#clean Length and Width
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
j <- as.character(ijsem$Length)
j <- gsub(' in diameter| \\(diameter\\)|mm|Indicated|diameter: | ', '', j)
j <- gsub(',', '.', j, fixed = TRUE)
j[grepl(">", j)] <- NA
j[j == ''] <- NA
j[j == c('0.61.6')] <- '0.6-1.6'
j[j == c('0.4 0.9 ')] <- '0.4-0.9'
j <- gsub('äóñ *', '-', j)

j <- ifelse(grepl("-[a-zA-Z]+", j) & !is.na(j), 
            as.character(format(as.Date(j, "%d-%b"), "%d.%m")), 
            j)
j <- ifelse(grepl("[a-zA-Z]+-", j) & !is.na(j), 
            as.character(format(as.Date(paste(1, j), "%d %b-%y"), "%m.%y")), 
            j)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", j) & !is.na(j)
j[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[2]]))))

#if it has a backslash, assume it was a decimal point
filt <- grepl("\\/", j) & !is.na(j)
j[filt] <- paste(sapply(strsplit(j[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(j[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Length <- swan(j)


#now width

j <- as.character(ijsem$Width)
j <- gsub(' |\\(|\\)|not.+|indiameter|Filament.+', '', j)
j <- gsub(',', '.', j, fixed = TRUE)
j[grepl("<", j)] <- NA
j[j == ''] <- NA
j <- gsub('äóñ *', '-', j)
j <- ifelse(grepl("-[a-zA-Z]+", j) & !is.na(j), 
            as.character(format(as.Date(j, "%d-%b"), "%d.%m")), 
            j)
j <- ifelse(grepl("[a-zA-Z]+-", j) & !is.na(j), 
            as.character(format(as.Date(paste(1, j), "%d %b-%y"), "%m.%y")), 
            j)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", j) & !is.na(j)
j[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[2]]))))

#if it has a backsclash, assume it was a decimal point
filt <- grepl("\\/", j) & !is.na(j)
j[filt] <- paste(sapply(strsplit(j[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(j[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Width <- as.numeric(j)

#Turn motility into a binary
ijsem$Motility <- c(0,1,1,1,1)[match(ijsem$Motility, c('non-motile', 'flagella',
                                                       'motile, but unspecified structure', 'gliding', 'axial filament'))]

#turn spore-forming into binary
ijsem$Spore <- c(0,1)[match(ijsem$Spore, c('no','yes'))]

#turn gram status into a binary
ijsem$Gram <- c(0, 0.5, 1)[match(ijsem$Gram, c('negative','variable','positive'))]

# fix pH
j <- as.character(ijsem$pH_range)
j <- gsub(' |\\(|\\)|at 6\\.0 but not 5\\.5|>|<|not inidcated', '', j)
j <- gsub('MinimumpHof5|t|4\\.0\\+|5\\.0\\+;notat4\\.5|14\\.08', '', j)
j[j %in% c('13.05','4.3','7.1','7.2','5.5','7.5','30.04','5.1')] <- NA
j[j == '5.7and6.8'] <- '5.7-6.8'
j[j == '4.0and8.5'] <- '4-8.5'
j[j == '5,6,10,12,'] <- '5-12'
j[j == '5and9.5'] <- '5-9.5'
j[j == '5.59.5'] <- '5.5-9.5'

# remove anything with an inequality, or other special characters
j[grepl("<|>", j)] <- NA
j[j == ''] <- NA
j <- gsub('äóñ *', '-', j)
j <- ifelse(grepl("-[a-zA-Z]+", j) & !is.na(j), 
            as.character(format(as.Date(j, "%d-%b"), "%d.%m")), 
            j)

filt <- grepl("\\/", j) & !is.na(j)
j[filt] <- paste(sapply(strsplit(j[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(j[filt], '\\/'), function(x) x[[2]]), sep = '.')

j[!grepl("-", j) & !is.na(j) & !grepl("\\.", j)] <- NA
filt <- !grepl("-", j) & !is.na(j) & grepl("\\.", j)
j[filt] <- gsub("\\.", "-", j[filt])
lows <- swan(sapply(strsplit(j, '-'), function(x) x[[1]]))
highs <- swan(sapply(strsplit(j, '-'), function(x) if (length(x) > 1) x[[2]] else NA))

ijsem$pH_lows <- lows
ijsem$pH_highs <- highs

#if there isn't a pH optimum, but there are ph high and lows, we take the midpoint (same as done in Barberan 2016)
ijsem$pH_optimum <- ifelse(is.na(ijsem$pH_optimum), (lows + highs) / 2, ijsem$pH_optimum)

####GC
j <- as.character(ijsem$GC)
j <- gsub('äóñ *', '-', j)
j <- gsub(',', '\\.', j)
j[j %in% c('DQ987877','marine organism','BAOL01000001','40.50%','GU323338')] <- NA
j <- gsub('O±.+|\\+/-.+|%|o', '', j)
filt <- grepl("-", j) & !is.na(j)
j[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[2]]))))

ijsem$GC <- swan(j)

# select traits for analylsis...
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
# now, sporulation data from Browne et al 2016
# https://www.nature.com/articles/nature17645#supplementary-information 
# [data edited in excel for easier processing)

spo <- read.csv('Browne2016_sporulationTable.csv', header = T, skip = 1)

# fix long colnames
colnames(spo) <- substr(colnames(spo), 1, 2)

spo <- spo %>%
  transmute(
    Genus = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[1]])),
    Species = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[2]])),
    trait = 'Spore_score',
    val = si)

####################################### Exploring bacdive, https://bacdive.dsmz.de/

if (FALSE) {
  # currently, this takes a really long time (~20 hrs). I copy paste it into a new R window, which is why I am keeping the wd and packages
  
  #setup
  wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\data'
  setwd(wd)
  source('C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\succession_custom_functions.R')
  loadpax(c('data.table','tidyverse','RCurl','rjson','stringr'))
  
  #load tax
  tax_LTP <- readRDS('LTP_tax.RDS')
  tax_succ <- readRDS('succession_tax_SILVA.RDS')
  tax <- tax_succ %>%
    select(otu, Genus, Species) %>%
    bind_rows(tax_LTP) %>%
    mutate(source = ifelse(grepl('OTU', otu), 'DIAB', 'LTP'))
  
  #make list of needed bacdat queries
  #these include a Genus/Species combo when available, and a list of all genera (which may or may not be used in the final version)
  needed <- tax %>%
    filter(!is.na(Genus)) %>%
    mutate(needed = ifelse(!is.na(Genus) & !is.na(Species), paste0(Genus, '/', Species), NA)) %>%
    filter(!is.na(needed))
  needed <- sort(unique(c(needed$Genus, needed$needed)))
  
  #read previous backdat queries to know when to start
  bacdat <- readRDS('bacdive_Sept10_2018.RDS')
  needed <- needed[!needed %in% names(bacdat)]
  
  ##sometimes there is an error and a particular entry must be skipped
  #needed <- needed[2:length(needed)]
  
  #loop through each needed entry. Future large downloads may want to save each entry, because saving/resaving to hard-disk is unnecessarily resource intensive (but necessary to avoid loss in case of (semi-frequent) errors)
  for (i in needed) {
    bacdat[[i]] <- bac_search(i)
    saveRDS(bacdat, 'bacdive_Sept10_2018.RDS')
  }
  
} 

#load bacdive data
bacdat <- readRDS('bacdive_Sept10_2018.RDS')

#to determine what to filter for, refer to this list of fields:
bacdat_fields <- unlist(bacdat[[1]][[1]])
bacdat_fields <- bacdat_fields[!grepl('molecular_biology.sequence', bacdat_fields)]

#extract taxonomic data and desired traits from downloaded bacdive database entries
bacdat <- lapply(bacdat, function(x)
  lapply(x, function(y) 
    lapply(y, function(z) {
      z <- unlist(z)
      data.frame(
        Genus = z['taxonomy_name.strains.genus'],
        Species = z['taxonomy_name.strains.species_epithet'],
        Oxygen_tolerance = z['morphology_physiology.oxygen_tolerance.oxygen_tol'],
        Gram_positive = z['morphology_physiology.cell_morphology.gram_stain'],
        Motility = z['morphology_physiology.cell_morphology.motility'],
        Length = z['morphology_physiology.cell_morphology.cell_len'],
        Width = z['morphology_physiology.cell_morphology.cell_width'],
        Aggregation_score = z['morphology_physiology.multicellular_morphology.ability'],
        Spore = z['morphology_physiology.spore_formation.ability'],
        row.names = NULL, stringsAsFactors = FALSE)})))

bacdat <- bind_rows(lapply(bacdat, function(x) bind_rows(lapply(x, bind_rows))))

refs <- list(
  Oxygen_tolerance = c('obligate aerobe'      = 5,
                       'aerobe'               = 5,
                       'facultative aerobe'   = 4,
                       'aerotolerant'         = 4,
                       'microaerophile'       = 3,
                       'microaerotolerant'    = 3,
                       'facultative anaerobe' = 2,
                       'anaerobe'             = 1,
                       'obligate anaerobe'    = 1),
  Gram_positive = c('positive' = 1,
                    'variable' = 0.5,
                    'negative' = 0),
  Spore = c('TRUE'  = 1,
            'FALSE' = 0),
  Motility = c('TRUE'  = 1,
               'FALSE' = 0),
  Aggregation_score = c('TRUE'  = 1,
                        'FALSE' = 0)
)

bacdat <- bacdat %>%
  mutate_all(funs(ifelse(. == '', NA, .))) %>%
  mutate(
    Oxygen_tolerance = refs$Oxygen_tolerance[match(Oxygen_tolerance, names(refs$Oxygen_tolerance))],
    Gram_positive = refs$Gram_positive[match(Gram_positive, names(refs$Gram_positive))],
    Spore = refs$Spore[match(Spore, names(refs$Spore))],
    Motility = refs$Motility[match(Motility, names(refs$Motility))],
    Aggregation_score = refs$Aggregation_score[match(Aggregation_score, names(refs$Aggregation_score))],
    Length = ifelse(grepl('<|>', Length), NA, Length),
    Length = ifelse(grepl('-', Length), 
                    rowMeans(cbind(as.numeric(sub('-.*', '', Length)),
                                   as.numeric(sub('.*-', '', Length)))), Length),
    Width = ifelse(grepl('<|>', Width), NA, Width),
    Width = ifelse(grepl('-', Width), rowMeans(cbind(as.numeric(sub('-.*', '', Width)),
                                                     as.numeric(sub('.*-', '', Width)))), Width))

#convert to narrow table
bacdat <- bacdat %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val)) %>%
  mutate(val = as.numeric(val))

###############################

## get Genome size, GC content, and Gene Number from NCBI ftp site... (Which is down now??)
genos <- fread('NCBI_prokaryotes.txt') %>% 
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

# also found this -- not sure how much overlap there is between them, but am merging them to be safe
#https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/

genos2 <- fread('genome_size_NCBI.csv') %>%
  transmute(sp = `#Organism Name`, val = `Size(Mb)`) %>%
  mutate(sp = gsub(' bacterium.*', '', sp),
         sp = gsub("'.*", "", sub("'?", "", sp))) %>%
  filter(!grepl('\\d', sp)) %>%
  mutate(sp = ifelse(str_count(sp, "\\S+") > 1, word(sp, 1, 2), sp)) %>%
  separate(sp, c("Genus", "Species"), " ", fill = 'right') %>%
  transmute(Genus, Species, trait = 'Genome_Mb', val)

genos <- bind_rows(genos, genos2)

####################
# extract 16s gene copy numbers from rrnDB
#https://rrndb.umms.med.umich.edu/static/download/
#i modified it a bit in excel before processing

rrnDB <- read.csv('rrnDB-5.3.csv') %>%
  mutate(
    sp = gsub("'|\\[|\\]", "", NCBI.scientific.name),
    Genus = sapply(sapply(sp, strsplit, ' '), function(x) x[[1]]),
    Species = sapply(sapply(sp, strsplit, ' '), function(x) ifelse(length(x) > 1, x[[2]], NA)),
    trait = 'Copies_16S') %>%
  transmute(Genus, Species, trait, val = X16S.gene.count) %>%
  filter(!is.na(val)) %>%
  mutate(Species = ifelse(is.na(Species), 'unclassified', Species))

################ extract IgA data from Palm et al 2014
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4174347/
#i modified it a bit in excel before processing

iga <- read.csv('Palm2014_IGA.csv', stringsAsFactors = FALSE)
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
  summarise(val = log(mean(ici) + 1)) %>%
  as.data.frame()

#################
#extract number of b-vitamins in genomes from Magnusdottir 2015
#https://www.frontiersin.org/articles/10.3389/fgene.2015.00148/full
bvit <- read.csv('Bvitamins_Magnusdottir2015.csv', stringsAsFactors = FALSE) %>%
  mutate(Genus = sapply(sapply(tax, strsplit, ' '), function(x) x[[1]]),
         Species = sapply(sapply(tax, strsplit, ' '), function(x) x[[2]])) %>%
  mutate(Species = ifelse(Species == 'sp.', 'unclassified', Species)) %>%
  select(-NCBI.Taxonomy.ID, -tax, -Body.Site) %>%
  gather(Bvit, val, -Genus, -Species) %>%
  group_by(Genus, Species, trait = 'B_vitamins') %>%
  summarise(val = length(unique(Bvit[val > 0]))) %>%
  as.data.frame()

################################

#load traits from JGI
# currently I drop otu-level taxonomy, but I could probably figure it out with
# https://gold.jgi.doe.gov/organisms?id=Go0384442
# (ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/ could be useful later)
jgi <- fread("GOLD-328_Organism_Metadata_02232018.txt", quote = "")

#remove duplicate oxygen column
jgi <- jgi[, -9]

#convert null to NA
jgi[jgi == '(null)'] <- NA

#clean columns
jgi <- jgi %>%
  transmute(
    Species = SPECIES, 
    Width = CELL_DIAMETER, 
    Length = CELL_LENGTH, 
    Gram_positive = GRAM_STAIN,
    Motility = MOTILITY,
    Oxygen_tolerance = OXYGEN_REQUIREMENT, 
    pH_optimum = PH,
    Spore = SPORULATION)

#Fix taxonomy
jgi <- separate(jgi, Species, c("Genus","Species"), extra = 'merge', fill = "right") %>%
  mutate(Species = gsub("sp. ", "", Species))

# deal with width - take averages of ranges when necessary.
jgi <- jgi %>%
  mutate(
    Width = ifelse(Width == '803nm', 0.803, Width),
    Width = gsub(" |[[:alpha:]]|¿|¼|\\?|\\&#956;|`", "", Width)
  ) %>%
  separate(Width, c("Width0","Width1"), fill = 'right', sep = '-') %>%
  mutate(Width = ifelse(is.na(Width0), Width0, 
                        ifelse(is.na(Width1), as.numeric(Width0),
                               (as.numeric(Width1) + as.numeric(Width0)) / 2))) %>%
  mutate(Width = as.numeric(Width)) %>%
  select(-Width0, -Width1)

# same with length
jgi <- jgi %>%
  mutate(
    Length = ifelse(Length == '803nm', 0.803, Length),
    Length = gsub(" |[[:alpha:]]|¿|¼|\\?|\\&#956;|`|\\&#8197;", "", Length),
    Length = gsub(".5.0", "5.0", Length, fixed = T),
    Length = gsub("1.55.0", "1.55", Length, fixed = T),
    Length = gsub("0.81.2", "0.8-1.2", Length, fixed = T)
    
  ) %>%
  separate(Length, c("Length0","Length1"), fill = 'right', sep = '-') %>%
  mutate(Length = ifelse(is.na(Length0), Length0, 
                         ifelse(is.na(Length1), as.numeric(Length0),
                                (as.numeric(Length1) + as.numeric(Length0)) / 2))) %>%
  mutate(Length = as.numeric(Length)) %>%
  select(-Length0, -Length1)

# gram positive, motility, oxygen tolerance, sporulation
jgi <- jgi %>%
  mutate(
    Gram_positive = c(1,0)[match(Gram_positive, c('Gram+','Gram-'))],
    Motility = c(1,1,0,0)[match(Motility, c('Motile','Chemotactic','Non-motile','Nonmotile'))],
    Oxygen_tolerance = c(5,5,4,3,2,1,1)[match(jgi$Oxygen_tolerance, 
                                              c('Obligate aerobe','Aerobe','Microaerophilic','Facultative','Facultative anaerobe','Anaerobe','Obligate anaerobe'))],
    Spore = c(1,1,0)[match(Spore, c('Non-sporulating','Nonsporulating','Sporulating'))])

#ph optimum
jgi <- jgi %>%
  mutate(
    pH_optimum = gsub(" |~", "", pH_optimum),
    pH_optimum = ifelse(pH_optimum %in% c('acido-sensible','Notknown'), NA, pH_optimum)) %>%
  separate(pH_optimum, c('pH0','pH1'), fill = 'right', sep = '-') %>%
  mutate(pH_optimum = ifelse(is.na(pH0), pH0, 
                             ifelse(is.na(pH1), as.numeric(pH0),
                                    (as.numeric(pH1) + as.numeric(pH0)) / 2))) %>%
  mutate(pH_optimum = as.numeric(pH_optimum)) %>%
  select(-pH0, -pH1)

jgi <- jgi %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val))

###################################################

# put it all together
x <- bind_rows(
  mutate(ijsem, source = 'IJSEM'),
  mutate(spo, source = 'Browne2016'),
  mutate(bacdat, source = 'BacDive'),
  mutate(genos, source = 'NCBI'),
  mutate(rrnDB, source = 'rrnDB'),
  mutate(iga, source = 'Palm2014'),
  mutate(bvit, source = 'Mag2015'),
  mutate(jgi, source = 'JGI')
)

#manual entry
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3346390/
x$val[x$Genus == 'Bifidobacterium' & x$Species != 'unclassified' & x$trait == 'Aggregation_score'] <- 1 


x_unclean <- nrow(x)
### I looked at all the points plotted using 
#ggplot(x, aes(x = val)) + geom_density() + facet_wrap(~trait, scales = 'free')
# identifying and removing outliers
# remove inordinately large values...

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

print(paste(x_unclean - nrow(x), "outliers of", nrow(x), "data points were removed"))

#Genus/Species fixes
#any species with a special character, including 'sp.', becomes 'unclassified'
x <- x %>% 
  mutate(
    Genus = trimws(Genus),
    Genus = gsub('\\[|\\]', '', Genus),
    Genus = ifelse(Genus == 'Sul¥å\u0081tobacter', 'Sulfitobacter', Genus),
    Genus = ifelse(grepl('Type-', Genus), NA, Genus),
    Genus = gsub('-.*', '', Genus),
    Species = ifelse(grepl('_| ', Genus), gsub('.*_', '', Genus), Species),
    Genus = gsub('_.*| .*', '', Genus),
    Species = trimws(Species),
    Species = ifelse(grepl('(?!")[[:punct:]]', Species, perl = TRUE), 'unclassified', Species))

#for plotting comparisons by source
x_by_source <- x %>%
  group_by(Genus, Species, source, trait) %>%
  summarise(val = ifelse(trait[1] %in% c('Length','Width'), log(mean(val)), mean(val))) %>%
  spread(trait, val)

study_names <- c(Mag2015 = "Magnúsdóttir et al. 2015", 
                 Palm2014 = "Palm et al. 2014", 
                 Browne2016 = "Browne et al. 2016")

trait_sources <- x_by_source %>% 
  ungroup() %>%
  mutate(source = ifelse(source %in% names(study_names), study_names[match(source, names(study_names))], source)) %>%
  gather(trait, val, -Genus, -Species, -source) %>% 
  filter(!is.na(val)) %>% 
  mutate(trait = ifelse(trait %in% c('Spore','Spore_score'), 'Sporulation', trait)) %>%
  group_by(Trait = trait) %>% 
  summarise(Sources = paste(sort(unique(source)), collapse = '; '))

#calculate means among species. Log length/width data. trim whitespace
x <- x %>%
  group_by(Genus, Species, trait) %>%
  summarise(val = ifelse(trait[1] %in% c('Length','Width'), log(mean(val)), mean(val))) %>%
  spread(trait, val)

# Merge spore scores (if no score and ijsem says 0, then 0; 
# otherwise the median of spore score when we know for ijsem spore == 1
x <- x %>% 
  mutate(
    Sporulation = ifelse(is.na(Spore_score), 
                         ifelse(Spore > 0, median(Spore_score[Spore > 0], na.rm = T), 0), Spore_score)) %>%
  ungroup() %>%
  select(-Spore, -Spore_score)

#remove unclassified and NAs
x <- x %>%
  filter(!(is.na(Genus) | Genus == 'unclassified')) %>%
  filter(!(is.na(Species) | Species == 'unclassified'))

#remove taxa not present in the Living tree Project or the Silva-derived taxonomy file from our usearch pipeline
tmp <- bind_rows(select(tax_succ, Genus, Species), select(tax_LTP, Genus, Species)) %>%
  filter(!(Genus == 'unclassified' | Species == 'unclassified'))
x <- filter(x, paste(Genus, Species) %in% paste(tmp$Genus, tmp$Species))

###save RDS
saveRDS(x, file = 'traits_sparse.RDS')
saveRDS(trait_sources, file = 'trait_sources.RDS')
