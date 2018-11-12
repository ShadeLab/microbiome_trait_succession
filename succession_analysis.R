#### Setup ####
#setwd
wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\data\\'

#load custom functions
source('C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\succession_custom_functions.R')

#load packages
loadpax(pkg = c('Hmisc','gtable','dplyr','tidyr','knitr','data.table','vegan','GGally','grid','gridExtra','betapart','scales','stringr', 'phyloseq','ape','kableExtra','poilog','ggtree','ggplot2','scales'))

#load data
otus <- readRDS(file = paste0(wd, 'otus.RDS'))
meta <- readRDS(file = paste0(wd, 'subject_metadata.RDS'))
tax <- readRDS(file = paste0(wd, 'succession_tax_SILVA.RDS'))
tax_LTP <- readRDS(file = paste0(wd, 'tax_LTP.RDS'))
trait_deltas <- readRDS(file = paste0(wd, 'trait_deltas.RDS'))
trait_delta_models <- readRDS(file = paste0(wd, 'trait_delta_models.RDS'))
traits <- readRDS(file = paste0(wd, 'traits.RDS'))
traits_sparse <- readRDS(file = paste0(wd, 'traits_sparse.RDS'))
trait_predictions <- readRDS(file = paste0(wd, 'traits_all_predictions.RDS'))
trait_sources <- readRDS(file = paste0(wd, 'trait_sources.RDS'))
tree <- read.tree(file = paste0(wd, 'LTP_succession.tree'))

#create versions of otu abundances with and without c-section infants
#First, assign C-section data, skinny table
otus_cs <- otus %>% left_join(meta[, c('subject','delivery')], by = 'subject')

# Without C-section data, skinny table
otus <- otus_cs %>%
  #filter(delivery != 'caesaren') %>%
  select(-delivery)

#wide versions of otus
otus_wide_cs <- otus_cs %>% spread(otu, abun, fill = 0)
otus_wide <- otus %>% spread(otu, abun, fill = 0)

#wide version of traits
traits_wide <- traits %>% spread(trait, val)

#trait renaming dictionary
trait_names <- c(
  "Aggregation_score" = "Aggregation score",
  "B_vitamins" = "B vitamins",
  "Copies_16S" = "16S gene copies",
  "GC_content" = "GC content",
  "Gene_number"= "Gene number",
  "Genome_Mb" = "Genome size",
  "Gram_positive" = "Gram-positive",
  "IgA"        = "IgA binding affinity",
  "Length"     = "Length",
  "Motility"   = "Motility",
  "Oxygen_tolerance" = "Oxygen tolerance",
  "pH_optimum" = "pH optimum",
  "Salt_optimum" = "Salt optimum",
  "Sporulation" = "Sporulation score",
  "Temp_optimum" = "Temperature optimum",
  "Width"      = "Width"
)

#### Trait sources ####
tabTraitSources <- function(){}

trait_units <- data.frame(stringsAsFactors = FALSE,
                          "Aggregation_score" = "0 (never) to 1 (observed aggregation)",
                          "B_vitamins" = "No. B-vitamin pathways in genome",
                          "Copies_16S" = "No. in 16S rRNA gene copies in genome",
                          "GC_content" = "Percent (\\%) guanine and cytosine in genome",
                          "Gene_number" = "No. genes in genome",
                          "Gram_positive" = "0 (Gram-negative) to 1 (Gram-positive)",
                          "IgA" = "log ([IgA+]/[IgA-] + 1)",  
                          "Length" = "log ($\\mu$m)",
                          "Motility" = "0 (never motile) to 1 (always motile)",  
                          "Oxygen_tolerance" = "0 (obligate anaerobe) to 5 (obligate aerobe)",
                          "pH_optimum" = "pH",    
                          "Salt_optimum" = "g/l",  
                          "Sporulation" = "0 (never sporulates) to 1 (sporulates easily)",  
                          "Temp_optimum" = "$^{\\circ}$C",  
                          "Width" = "log ($\\mu$m)"
)

tabTraitSources <- trait_units %>%
  gather(Trait, Units) %>%
  left_join(trait_sources, by = 'Trait') %>%
  mutate(Trait = trait_names[match(Trait, names(trait_names))]) %>%
  mutate(Sources = gsub("ú", "\\'{u}", fixed = TRUE, Sources),
         Sources = gsub("ó", "\\'{o}", fixed = TRUE, Sources))

figTraitCorr <- function(){}

j <- expand.grid(unique(traits$trait), unique(traits$trait)) %>%
  group_by(Var1, Var2) %>%
  do(ct = cor.test(traits_wide[[as.character(.$Var1)]], traits_wide[[as.character(.$Var2)]])) %>%
  mutate(Var1 = factor(trait_names[match(Var1, names(trait_names))], levels = trait_names),
         Var2 = factor(trait_names[match(Var2, names(trait_names))], levels = trait_names),
         cor = ct$estimate[[1]],
         pval = ct$p.value,
         symb = ifelse(pval < 0.001, '***',
                       ifelse(pval < 0.01, '**',
                              ifelse(pval < 0.05, '*', ''))))

j <- filter(j, as.numeric(Var1) > as.numeric(Var2))


figTraitCorr <- j %>%
  ggplot(aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = 'black') +
  geom_text(aes(label = symb)) +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', name = 'Pears. Corr.') +
  labs(x = '', y = '') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

#### Coverage ####
figSampling <- function(){}

figSampling <- otus_cs %>% 
  distinct(subject, t, delivery) %>% 
  left_join(meta, by = c("subject", "delivery")) %>% 
  group_by(subject) %>%
  mutate(t_min = min(t), t_max = max(t), len = t_max - t_min) %>%
  ungroup() %>%
  arrange(desc(len)) %>%
  mutate(
    delivery = ifelse(delivery == 'caesaren', 'Caesarean delivery','Vaginal delivery'),
    Subject = as.numeric(factor(subject, levels = unique(subject))),
    t_min = min(t)) %>%
  ggplot(aes(y = Subject, fill = delivery, color = delivery)) + 
    geom_segment(aes(x = t_min, xend = t_max, yend = Subject)) + 
    geom_point(aes(x = t), shape = 21) + 
    scale_fill_manual(values = c('black','grey60'), name = '') +
    scale_color_manual(values = c('black','grey60'), name = '') +
    theme_bw() +
    theme(
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank()) +
    labs(x = 'Months after birth', y = 'Subject')

figTraitCoverage <- function(){}

ts <- traits_sparse %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val))

trait_coverage <- otus_cs %>%
  group_by(otu) %>%
  summarise(abun = sum(abun)) %>%
  left_join(tax[, c('otu', 'Genus', 'Species')], by = 'otu') %>%
  left_join(traits_sparse[, c('Genus','Species','Aggregation_score','pH_optimum','Salt_optimum','IgA')], by = c('Genus','Species')) %>%
  left_join(traits_wide, by = 'otu') %>%
  gather(trait, val, -otu, -abun, -Genus, -Species) %>%
  mutate(Source = ifelse(paste(trait, Genus, Species) %in% paste(ts$trait, ts$Genus, ts$Species), 'Observed data',
                         ifelse(!is.na(val), 'Inferred data', 'Missing data'))) %>%
  group_by(trait, Source) %>%
  summarise(abun = sum(abun)) %>%
  group_by(trait) %>%
  mutate(abun = abun / sum(abun),
         Source = factor(Source, levels = c('Missing data','Inferred data','Observed data'))) %>%
  ungroup() %>%
  mutate(
    trait = trait_names[match(trait, names(trait_names))],
    trait = factor(trait, rev(sort(unique(trait)))))

figTraitCoverage <- ggplot(trait_coverage, aes(x = trait, y = abun, fill = Source)) + 
  geom_bar(color = 'black', stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c('lightgrey','skyblue','blue'), name = '') +
  labs(x = "", y = 'Proportion of total abundance') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom') +
  guides(fill = guide_legend(reverse = TRUE))

#......................................................
figTreeTraits <- function(){}

#load tax
tax_tmp <- bind_rows(tax_LTP, tax[, c('otu','Genus','Species')]) %>%
  mutate(spp = paste(Genus, Species)) %>%
  select(-Genus, -Species)

# format traits
tmp <- traits_sparse %>% 
  mutate(spp = paste(Genus, Species)) %>%
  select(-Genus, -Species) %>%
  gather(trait, val, -spp)

tmp2 <- list(
  `OTU from this study` = tree$tip.label[grepl('OTU', tree$tip.label)],
  `OTU from the Living Tree Project` = tree$tip.label[!grepl('OTU', tree$tip.label)])
tree <- groupOTU(tree, tmp2)

j <- data.frame(otu = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(tax_tmp[, c('otu','spp')], by = 'otu') %>%
  left_join(tmp, by = 'spp') %>%
  mutate(val = ifelse(!is.na(val), jitter(as.numeric(as.factor(trait))), NA)) %>%
  filter(!is.na(val)) %>%
  mutate(trait = trait_names[match(trait, names(trait_names))]) %>%
  arrange(trait) %>%
  mutate(trait = paste0(as.numeric(as.factor(trait)), ': ', trait),
         trait = factor(trait, levels = unique(trait))) %>%
  rename(id = otu)

#ggtree(tree, size = 0.05, layout = 'circular') +
#  geom_tiplab2(label = '_______', 
#               size = 3, aes(color = group))

p <- ggtree(tree, size = 0.05) +
  geom_tiplab(label = '____________________________________________________', vjust = -.29, offset = 0.01, 
              size = 1, aes(color = group))

figTreeTraits <- facet_plot(p, panel = "Trait observations", data = j, geom = geom_point, 
           aes(x = val, fill = trait), size = 0.3, shape = 21, 
           color = 'transparent', stroke = 0) + 
  theme_tree2() + 
  scale_color_manual(values = c('black','red')) +
  scale_x_continuous(breaks = c(1:16)) +
  guides(
    fill = guide_legend(override.aes = list(size = 4, alpha = 1, color = 'black'), 
                        order = 0, ncol = 2), 
    color = guide_legend(ncol = 1, order = 1, override.aes = list(label = '-', size = 7), label.vjust = .35)) +
  theme(legend.position = 'bottom',
        legend.title = element_text(color = NA))

figTreeTraits

#### HSP comparisons ####
#Pearson correlations among exploring differences among different methods of hidden state character predictions
#note that I remove traits that were not amenable to hidden state prediction, and remove actual observations (i.e., when dist == 0)
tabHSP <- function(){}

tabHSP <- trait_predictions %>%
  filter(!trait %in% c("Aggregation_score","IgA","pH_optimum","Salt_optimum")) %>%
  filter(dist > 0) %>%
  group_by(trait) %>%
  do(data.frame(
    trait = .$trait[1], 
    method1 = c('IC','SBT','SCP'), 
    cor(.[, c('IC','SBT','SCP')]), stringsAsFactors = FALSE)) %>%
  gather(method2, cor, IC, SBT, SCP) %>%
  filter(method1 != method2) %>%
  arrange(trait, method1) %>%
  group_by(trait, method1) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  transmute(Trait = trait_names[match(trait, names(trait_names))],
            xy = paste(method1, 'x', method2), 
            cor) %>%
  spread(xy, cor)

#### Phylogenetic signal ####

#create a dataframe of null expectations of delta for each trait
#To remove phylogenetic effects, null delta is calculated as the mean pairwise trait-based distance for all randoly selected pairs of otus with greater than 10% difference in their 16S V4 region -- except for Salt optimum and ph optimum, which have no observable phylogenetic effect
#.................................................................................
tabThresholds <- function(){}

#data frame of maximum distances used to infer traits, for each trait
tdms <- mutate(trait_delta_models, trait = trait_names[match(trait, names(trait_names))])

max_dists <- tdms %>%
  distinct(trait, max_dist, .keep_all = FALSE) %>%
  mutate(max_dist = format(round(max_dist, 3), nsmall = 3)) %>%
  rename(Trait = trait, `Max. distance` = max_dist) %>%
  as.data.frame ()

tabThresholds <- max_dists

figThresholds <- function(){}

nulls <- trait_deltas %>%
  mutate(trait = as.vector(trait_names[match(trait, names(trait_names))])) %>%
  left_join(max_dists, by = c('trait' = 'Trait')) %>%
  group_by(trait) %>%
  filter(dist > `Max. distance`) %>%
  summarise(null = mean(delta))

#calculate means and sd of bins. Add a row for null deltas. Clean up trait names.
delta_bins <- trait_deltas %>%
  mutate(trait = trait_names[match(trait, names(trait_names))]) %>%  
  filter(dist_bin <= 0.2) %>%
  mutate(dist_bin = dist %/% 0.005 * 0.005) %>% 
  group_by(trait, dist_bin) %>%
  summarise(
    mean = mean(delta),
    sd = sd(delta),
    null = nulls$null[match(trait[1], nulls$trait)]) %>%
  ungroup()

#plot with standard deviations
figThresholds <- delta_bins %>%
  filter(!(trait == 'Sporulation' & mean > 0.2)) %>%
  filter(!(trait == 'Length' & mean > 1.5)) %>%
  filter(!(trait == 'Width' & mean > 1.2)) %>%
 ggplot(aes(x = dist_bin, y = mean)) +
    #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, color = 'grey', data = filter(delta_bins, !is.na(sd))) +
    geom_point() +
    geom_line(aes(y = delta, color = type), data = filter(tdms, max_dist > 0.03), lwd = 1.5) +
    geom_hline(aes(yintercept = null), linetype = 3) +
    geom_vline(aes(xintercept = max_dist), data = tdms, lty = 2) +
    facet_wrap(~trait, scales = 'free') +
    scale_color_discrete(name = '') +
    theme_bw() +
    theme(
      legend.position = 'bottom',
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank()) +
    labs(x = 'Phylogenetic distance between tips', y = 'Mean difference in traits between OTUs')


#.................................................................................
figTraitAbunCoverage <- function(){}

#direct observations
obs_dir <- traits_sparse %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val))

#determine if nearest_measured_otu is close enough to phylogenetically infer trait values
#append overall OTU abundances
tmp <- otus %>%
  left_join(tax[, c('Genus','Species','otu')], by = 'otu') %>%
  left_join(traits_wide, by = 'otu') %>%
  gather(trait, val, -sampleID, -subject, -t, -otu, -abun, -Genus, -Species) %>%
  mutate(Type = ifelse(paste(Genus, Species, trait) %in% paste(obs_dir$Genus, obs_dir$Species, obs_dir$trait), 'Observed', ifelse(!is.na(val), 'Phylogenetically inferred', 'Insufficient data'))) %>%
  mutate(Type = factor(Type, levels = c( 'Insufficient data','Phylogenetically inferred','Observed'))) %>%
  group_by(trait, Type) %>%
  summarise(rich = length(unique(otu)), abun = sum(abun)) %>%
  group_by(trait) %>%
  mutate(relabun = abun / sum(abun))

#summary barplot (abundances of OTUS by trait/data type)
figTraitAbunCoverage <- ggplot(tmp, aes(x = trait, y = relabun, fill = Type)) + 
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c('lightgrey','#6666ff','blue')) +
  labs(x = "Percent Relative abundance across all samples", y = '') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE)) +
  th

#### Taxonomic patterns ####

figTaxOverTime <- function(){}

# Bump up any t<2 to 2 to avoid early noise due to low-sampling.
# Determine centers of colonization periods of OTUs by taking the abundance-weighted mean occurence. 
# Focus only on taxa more than .1% abundant, 
#       only taxa in at least 5 individuals (~10% of analysis population), 
#       only taxa that appear for 3 or more consecutive months.

#after filtering, calculate start, mid, and end times of each otu
mods <- otus_wide %>%
  gather(otu, abun, -sampleID, -subject, -t) %>%
  group_by(otu) %>%
  do(mod = summary(lm(abun ~ t, .))) %>%
  mutate(
    pval = mod$coefficients[[8]],
    tval = mod$coefficients[[6]],
    group = ifelse(pval < 0.05 & tval < 0, 'Early successional', 'Mid-successional / No trend'),
    group = ifelse(pval < 0.05 & tval > 0, 'Late successional', group),
    group = factor(group, levels = c('Early successional','Mid-successional / No trend','Late successional'))) %>%
  select(otu, group)

j <- otus %>%
  left_join(mods, by = 'otu') %>%
  group_by(group, t = round(t)) %>%
  filter(t >= 2 & t <= 36) %>%
  ungroup()

p1 <- j %>%
  group_by(t) %>%
  mutate(abun = abun / sum(abun)) %>%
  ggplot(aes(x = t, y = abun, fill = group, color = group)) +
  geom_bar(stat = 'identity', width = 1) +
  facet_wrap(~group, ncol = 1) +
  labs(x = 'Months after birth', y = 'Relative abundance') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none')

p2 <- j %>%
  left_join(tax, by = 'otu') %>%
  group_by(group, Family) %>%
  summarise(abun = sum(abun)) %>%
  group_by(group) %>%
  mutate(rank = ifelse(Family == 'unclassified', 0, -as.numeric(as.factor(abun))),
         rank = as.numeric(as.factor(rank)),
         Family = ifelse(rank > 5, 'All others', Family)) %>%
  group_by(group, Family, rank) %>%
  summarise(abun = sum(abun)) %>%
  ungroup() %>%
  mutate(Family = as.character(Family),
         Family = ifelse(group == 'Early successional', paste0(' ', Family), 
                         ifelse(group == 'Mid-successional / No trend', paste0('  ', Family), Family))) %>%
  arrange(rank) %>%
  mutate(Family = factor(Family, unique(rev(Family)))) %>%
  mutate(abun = abun / sum(abun)) %>%
  group_by(group, Family) %>%
  summarise(abun = sum(abun)) %>%
  ggplot(aes(x = Family, y = abun, fill = group)) +
    geom_bar(stat = 'identity', color = 'black') +
    geom_text(aes(label = Family), hjust = 'left', nudge_y = 0.005) +
    facet_wrap(~group, ncol = 1, scales = 'free_y') +
    labs(x = "", y = "Relative abundance across all samples") +
    scale_y_continuous(limits = c(0, 0.325)) +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = 'none')


figTaxOverTime <- grid.arrange(grobs = list(p1, p2), widths = c(1, 1))

#......................................................................
figCWMsSuccGroup <- function(){}

labs <- c(Early = 'Early successional',
          Other = 'Mid-successional / No trend', 
          Late = 'Late successional')

tmp <- traits %>%
  mutate(trait = trait_names[match(trait, names(trait_names))]) %>%
  left_join(distinct(j, otu, group), by = 'otu') %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(names(labs)[match(group, labs)], names(labs)))

statz1 <- tmp %>%
  filter(group != 'Late') %>%
  group_by(trait) %>%
  do(mod = t.test(val ~ group, data = .)) %>%
  mutate(pval = mod$p.value, group = '1/2')

statz2 <- tmp %>%
  filter(group != 'Early') %>%
  group_by(trait) %>%
  do(mod = t.test(val ~ group, data = .)) %>%
  mutate(pval = mod$p.value, group = '2/3')
  
statz3 <- tmp %>%
  filter(group != 'Other') %>%
  group_by(trait) %>%
  do(mod = t.test(val ~ group, data = .)) %>%
  mutate(pval = mod$p.value, group = '1/3')

statz <- bind_rows(statz1, statz2, statz3) %>%
  filter(pval < 0.05) %>%
  ungroup() %>%
  group_by(trait) %>%
  summarise(label = 'NA',
            label = ifelse('1/2' %in% group & '1/3' %in% group & '2/3' %in% group, 'a_b_c', label),
            label = ifelse('1/2' %in% group & '1/3' %in% group & !'2/3' %in% group, 'a_b_b', label), 
            label = ifelse('1/2' %in% group & !'1/3' %in% group & '2/3' %in% group, 'a_b_a', label),
            label = ifelse('1/2' %in% group & !'1/3' %in% group & !'2/3' %in% group, 'a_b_ab', label),
            label = ifelse(!'1/2' %in% group & '1/3' %in% group & '2/3' %in% group, 'a_a_b', label),
            label = ifelse(!'1/2' %in% group & '1/3' %in% group & !'2/3' %in% group, 'a_ab_b', label),
            label = ifelse(!'1/2' %in% group & !'1/3' %in% group & '2/3' %in% group, 'ab_a_b', label),
            label = ifelse(!'1/2' %in% group & !'1/3' %in% group & !'2/3' %in% group, 'a_a_a', label)) %>%
  separate(label, c('Early','Other','Late'), sep = '_') %>%
  gather(group, label, -trait)

tmp$label = statz$label[match(paste(tmp$trait, tmp$group), paste(statz$trait, statz$group))]

yvals <- c(
  `16S gene copies` = 5.35,
  `GC content` = 51.5,
  Genes = 3600,
  `Gram-positive` = 0.68,
  Length = 1.15,
  Motility = 0.4,
  `Oxygen tolerance` = 2.5,
  `Sporulation score` = 0.25,
  `Temperature optimum` = 40
)

tmp <- tmp %>% mutate(yval = yvals[match(trait, names(yvals))])

figCWMsSuccGroup <- ggplot(tmp, aes(x = group, y = val)) +
  facet_wrap(~trait, scales = 'free') +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_text(aes(label = label, y = yval), data = filter(tmp, !is.na(label))) +
  th +
  labs(x = 'OTU successional group', y = 'Trait value')

figCWMsSuccGroup

#### Traits over time ####
figCWMs <- function(){}

j <- otus %>%
  left_join(traits_wide, by = 'otu') %>%
  gather(trait, val, -sampleID, -subject, -t, -otu, -abun) %>%
  filter(!is.na(val)) %>%
  mutate(trait = ifelse(trait %in% names(trait_names), trait_names[match(trait, names(trait_names))], trait)) %>%
  left_join(meta, by = 'subject') %>%
  group_by(subject, t, trait) %>%
  summarise(cwm = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(t = round(t)) %>%
  filter(t < 37 & t > 1)

figCWMs <- ggplot(j, aes(x = t, y = cwm)) +
  stat_summary(fun.y = mean, color = 'red', geom = "point") + 
  stat_summary(fun.data = mean_cl_boot, color = 'red', geom = "linerange") +
  stat_smooth(color = 'black', se = FALSE) +
  facet_wrap(~trait, scale = 'free_y', ncol = 3) +
  labs(x = "Months after birth", y = "Abundance-weighted community mean") +
  th

figCWMsTreat <- function(){}

#prep data for delivery mode
tmp1 <- otus_cs %>% 
  filter(subject %in% meta$subject[meta$treatment_group == 'C-section']) %>%
  mutate(treatment = 'C-section', delivery = NULL)

#prep data for antibiotics
tmp2 <- otus_cs %>% 
  filter(subject %in% meta$subject[meta$treatment_group == 'Antibiotics']) %>%
  mutate(treatment = 'Antibiotics', delivery = NULL)
  
#prep control group
tmp3 <- otus %>% 
  filter(subject %in% meta$subject[meta$treatment_group == 'Control']) %>%
  mutate(treatment = 'Control')

#append controls
tmp1 <- bind_rows(tmp1, tmp3) %>% mutate(group = 'Delivery mode')
tmp2 <- bind_rows(tmp2, tmp3) %>% mutate(group = 'Antibiotics')

#put it together, join traits, calculate CWMs
cwmsx <- bind_rows(tmp1, tmp2) %>%
  left_join(traits, by = 'otu') %>%
  filter(!is.na(trait)) %>%
  group_by(group, treatment, sampleID, subject, t, trait) %>%
  summarise(val = weighted.mean(val, w = abun)) %>%
  ungroup() %>%
  filter(round(t) > 1 & round(t) < 36) %>%
  mutate(trait = trait_names[match(trait, names(trait_names))], 
         treatment = factor(treatment, levels = c('Control','Antibiotics','C-section')),
         tbin = ceiling(t/6) * 6 - 3)

#perform ttests, adjust pval for multiple comparisons
statz <- cwmsx %>%
  group_by(trait, group, tbin) %>%
  do(mod = t.test(val ~ treatment, data = .)) %>%
  mutate(
    pval = mod$p.value, p.adjust(mod$p.value, method = 'holm', n = 22))

#add pvals to cwm data
#add custom tbin for staggered plotting
tmp <- cwmsx %>%
  filter(!(group == 'Antibiotics' & treatment == 'Control')) %>%
  mutate(pval = statz$pval[match(paste(trait, group, tbin), 
                                 paste(statz$trait, statz$group, statz$tbin))],
         sig = ifelse(treatment != 'Control', pval < 0.05, FALSE),
         tbin0 = tbin,
         tbin = ifelse(treatment == 'Control', tbin,
                       ifelse(treatment == 'Antibiotics', tbin - 0.5, tbin + 0.5)))

#calculate mean CWMs within each tbin
tmp2 <- tmp %>%
  group_by(trait, tbin, group, treatment) %>%
  summarise(val = mean(val))

#find otpimal locations for plotting significance symbols above lines
tmp3 <- tmp %>%
  group_by(trait, tbin0, group, treatment, sig, pval) %>%
  summarise(val = mean(val)) %>%
  group_by(trait) %>%
  mutate(yval = max(val) + 0.5 * (max(val) - min(val)),
         val = ifelse(treatment == 'Antibiotics',
                      max(val) + 0.2 * (max(val) - min(val)),
                      max(val) + 0.4 * (max(val) - min(val)))) %>%
  filter(sig) %>%
  distinct(trait, tbin0, treatment, group, yval, val, pval) %>%
  mutate(symb = ifelse(pval < 0.001, '***',
                       ifelse(pval < 0.01, '**',
                              ifelse(pval < 0.05, '*', '')))) %>%
  ungroup() %>%
  rename(tbin = tbin0)

figCWMsTreat <- ggplot(tmp, aes(x = tbin, y = val, color = treatment)) +
  stat_summary(fun.data = "mean_cl_boot", size = .25) +
  geom_line(data = tmp2) +
  scale_shape_discrete(guide = FALSE) +
  geom_text(aes(label = symb), data = tmp3, show.legend = FALSE, size = 6) +
  geom_point(aes(y = yval), data = tmp3, alpha = 0) +
  facet_wrap(~trait, scales = 'free', ncol = 3) +
  scale_color_manual(values = c('black',"#00A480","#FF6200"), name = '') +
  scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
  expand_limits(x = 2) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  labs(x = 'Months after birth', y = 'CWM trait value')


figCWMsVariance <- function(){}

tmp <- otus %>%
  left_join(traits, by = 'otu') %>%
  filter(!is.na(val)) %>%
  left_join(meta, by = 'subject') %>%
  group_by(trait) %>%
  mutate(val = scale(val)) %>%
  group_by(subject, t, trait) %>%
  summarise(var = var(rep(val, times = abun))) %>%
  ungroup() %>%
  mutate(t = round(t)) %>%
  filter(t < 37 & t >= 2) %>%
  mutate(trait = trait_names[match(trait, names(trait_names))])

figCWMsVariance <- ggplot(tmp, aes(x = t, y = var, color = trait)) +
  stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange") +
  stat_smooth(se = FALSE) +
  facet_wrap(~trait, ncol = 3) +
  labs(x = "Months after birth", y = "Within-sample communtiy trait variance") +
  theme_bw() +
  theme(
    legend.position = 'none',
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

figSporulatorsTime <- function(){}

#number of sporulating taxa over time
figSporulatorsTime <- otus %>%
  left_join(traits, by = 'otu') %>%
  filter(trait == 'Sporulation') %>%
  mutate(val = ifelse(val > 0.1, 'Sporulation score > 0.1', 'Sporulation score < 0.1')) %>%
  group_by(subject, t, val) %>%
  summarise(n = length(unique(otu))) %>%
  ggplot(aes(x = t, y = n, color = factor(val))) +
  geom_point(alpha = 0.5) +
  stat_smooth() +
  scale_color_manual(name = '', values = c('red','black')) +
  labs(x = 'Months after birth', y = 'Number of OTUs per infant') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom')

#### Intrababy dissimilarity ####
figIntrababyDiss <- function(){}

#use meltdist() to calculate dissimilarities within subjects for each method
j1 <- otus_cs %>%
  group_by(subject) %>%
  do(meltdist(x = unique(.$sampleID), y = otus_wide_cs, method = 'bray')) %>%
  mutate(method = 'OTU-based dissimilarity',
         data = 'obs')

#weighted UniFrac
#j2 <- otus_cs %>%
#  group_by(subject) %>%
#  do(meltdist(x = unique(.$sampleID), y = otus_wide_cs, tree = tree, method = 'unifrac')) %>%
#  mutate(method = 'Weighted UniFrac distance',
#         data = 'obs')

#second method: we use CWMs for columns, and calculate scaled euclidean distance

j2 <- otus_cs %>%
  left_join(traits, by = 'otu') %>%
  filter(!is.na(trait)) %>%
  group_by(sampleID, subject, trait, t) %>%
  summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
  group_by(trait) %>%
  mutate(val = scale(val)) %>%
  spread(trait, val) %>%
  group_by(subject) %>%
  do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
  mutate(method = 'Trait-based dissimilarity',
         data = 'obs')

#now, let's create dummy versions of unifrac and trait-based euclidean with scrambled traits
#set.seed(7)
#j4 <- j2[0, ]

#for (i in 1:10) {
#  
#  random_tree <- tree
#  random_tree$tip.label <- sample(random_tree$tip.label)
#  
#  tmp <- otus_cs %>%
#    group_by(subject) %>%
#    do(meltdist(x = unique(.$sampleID), 
#                y = otus_wide_cs, 
#                tree = random_tree, 
#                method = 'unifrac')) %>%
#    mutate(method = 'Weighted UniFrac distance',
#           data = 'ran')
#  
#  j4 <- bind_rows(j4, tmp)
#  print(paste('done with', i))
#  flush.console()
#
#} 

#j4 <- j4 %>% 
#  group_by(subject, sample1, sample2, t1, t2, method, data) %>%
#  summarise(dist = mean(dist))

#now, let's create a dummy version with scrambled traits
set.seed(7)
j3 <- j2[0,]

for (i in 1:10) {
  
  traits_randomized <- traits_wide %>% 
    gather(trait, val, -otu) %>%
    group_by(trait) %>% 
    mutate(val = sample(val)) %>%
    filter(!is.na(val))
  
  tmp <- otus_cs %>%
    left_join(traits_randomized, by = 'otu') %>%
    filter(!is.na(trait)) %>%
    group_by(sampleID, subject, trait, t) %>%
    summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
    group_by(trait) %>%
    mutate(val = scale(val)) %>%
    spread(trait, val) %>%
    group_by(subject) %>%
    do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
    ungroup() %>%
    mutate(method = 'Trait-based dissimilarity',
           data = 'ran')
  
  j3 <- bind_rows(j3, tmp)
  print(paste('done with', i))
  flush.console()

} 

j3 <- j3 %>% 
  group_by(subject, sample1, sample2, t1, t2, method, data) %>%
  summarise(dist = mean(dist))

#combine
intradist <- bind_rows(j1, j2, j3)

#prune to dissimilarities between subsequent samples, as a measure of rate of temporal turnover
j4 <- intradist %>%
  group_by(method, subject, t1) %>%
  filter(t2 > t1) %>%
  filter(t2 == min(t2)) %>%
  mutate(group = 'Dissimilarity to next sample')

#prune to dissimilarity to final sample
j5 <- intradist %>%
  group_by(method, subject, t1) %>%
  filter(t2 == max(t2) & t1 != max(t2)) %>%
  mutate(group = 'Dissimilarity to final sample')

tb_bump <- j5 %>%
  ungroup() %>%
  filter(method == 'Trait-based dissimilarity') %>%
  mutate(t1 = t1 %/% 6 * 6) %>%
  filter(t1 == min(t1)) %>%
  summarise(diff = mean(dist[data == 'obs']) - 
                   mean(dist[data == 'ran'])) %>%
  pull(diff)

j5 <- j5 %>% 
  mutate(
    dist = ifelse(group == 'Dissimilarity to final sample' &
                  method == 'Trait-based dissimilarity' &
                  data == 'ran', 
                  dist + tb_bump, dist))

#combine and modify for plotting
j <- bind_rows(j4, j5) %>%
  ungroup() %>%
  mutate(tbin = (t1 %/% 6) * 6 + 3) %>%
  filter(tbin < 36) %>%
  mutate(
    group = factor(group, levels = unique(group)),
    method = factor(method, levels = c("OTU-based dissimilarity", 
      "Weighted UniFrac distance", "Trait-based dissimilarity")),
    data = factor(c('Observed','Null model')[match(data, c('obs','ran'))],
                  c('Observed','Null model')))

statz <- j %>%
  filter(method == "Trait-based dissimilarity") %>%
  group_by(method, group, subject, tbin, data) %>%
  summarise(dist = mean(dist)) %>%
  group_by(method, group, tbin) %>%
  do(mod = t.test(dist ~ data, .)) %>%
  mutate(pval = mod$p.value,
         sig = pval < 0.1,
         symb = ifelse(pval < 0.001, '***',
                 ifelse(pval < 0.01, '**',
                   ifelse(pval < 0.05, '*',
                     ifelse(pval < 0.1, '·', '')))))

j <- j %>% left_join(statz, by = c("method", "group", "tbin"))

#do I want to add a zero? nah...
#j <- bind_rows(mutate(j, group = as.character(group)), 
#               distinct(j, data, method) %>% 
#                 mutate(group = 'Dissimilarity to final sample', tbin = 39, dist = 0)) %>%
#  mutate(group = factor(group, levels = c('Dissimilarity to next sample', 'Dissimilarity to final sample')))

tmp <- j %>% 
  group_by(tbin, method, group, data, sig, symb) %>% 
  summarise(dist = mean(dist)) %>%
  group_by(tbin, method, group) %>%
  mutate(symb = ifelse(sig & dist == max(dist), symb, NA),
         sig = ifelse(is.na(symb), NA, sig))

figIntrababyDiss <- ggplot(j, aes(x = tbin, y = dist, lty = data)) +
    stat_summary(fun.y = mean, geom = "point", aes(shape = data)) + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange", data = filter(j, tbin < 39)) +
    geom_text(aes(label = symb, y = dist + 0.5), data = filter(tmp, sig), size = 6) +
    geom_line(data = tmp) +
    facet_grid(method ~ group, scales = 'free') +
    labs(x = 'Months after birth', y = '') +
    scale_linetype_discrete(name = '') +
    scale_shape_manual(values = c(16,1), name = '') +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      legend.position = 'bottom')

figIntrababyDiss

#### Interbaby dissimilarity ####
figInterbabyDiss <- function(){}

#use meltdist() to calculate dissimilarities within subjects for each method
#filter out low density tbins
tmp <- otus_cs %>%
  mutate(tbin = round(t)) %>%
  filter(tbin < 36)

j1 <- tmp %>%
  group_by(tbin) %>%
  do(meltdist(x = unique(.$sampleID), y = otus_wide_cs, method = 'bray')) %>%
  mutate(method = 'OTU-based dissimilarity',
         data = 'obs')

#second method: we use CWMs for columns, and calculate scaled euclidean distance
j2 <- tmp %>%
  left_join(traits, by = 'otu') %>%
  filter(!is.na(trait)) %>%
  group_by(sampleID, subject, trait, tbin, t) %>%
  summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
  group_by(trait) %>%
  mutate(val = scale(val)) %>%
  spread(trait, val) %>%
  group_by(tbin) %>%
  do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
  mutate(method = 'Trait-based dissimilarity',
         data = 'obs')

set.seed(7)
j3 <- j2[0, ]

for (i in 1:10) {
  
  traits_randomized <- traits_wide %>% 
  gather(trait, val, -otu) %>%
  group_by(trait) %>% 
  mutate(val = sample(val)) %>%
  filter(!is.na(val))

  tmp1 <- tmp %>%
    left_join(traits_randomized, by = 'otu') %>%
    filter(!is.na(trait)) %>%
    group_by(sampleID, subject, trait, tbin, t) %>%
    summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
    group_by(trait) %>%
    mutate(val = scale(val)) %>%
    spread(trait, val) %>%
    group_by(tbin) %>%
    do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
    mutate(method = 'Trait-based dissimilarity',
           data = 'ran')
  
  j3 <- bind_rows(j3, tmp1)
  print(paste('done with', i))
  flush.console()
}

j3 <- j3 %>% 
  group_by(tbin, sample1, sample2, t1, t2, method, data) %>%
  summarise(dist = mean(dist))


#weighted UniFrac
#j3 <- tmp %>%
#  group_by(tbin) %>%
#  do(meltdist(x = unique(.$sampleID), y = otus_wide_cs, tree = tree, method = 'unifrac')) %>%
#  mutate(method = 'Weighted UniFrac distance')

#combine
interdist <- bind_rows(j1, j2, j3) %>%
  mutate(subject1 = gsub("_.*", "", sample1),
         subject2 = gsub("_.*", "", sample2)) %>%
  filter(subject1 != subject2) %>%
  ungroup() %>%
  mutate(method = factor(method, 
                         levels = c("OTU-based dissimilarity", 
                                    "Weighted UniFrac distance", 
                                    "Trait-based dissimilarity")))

j <- interdist %>%
  group_by(method, data, tbin) %>%
  summarise(dist = mean(dist)) %>%
  mutate(tbin = (tbin %/% 6) * 6 + 3) %>%
  ungroup() %>%
  mutate(data = factor(c('Observed','Null model')[match(data, c('obs','ran'))], c('Observed','Null model')))

statz <- j %>%
  filter(method == "Trait-based dissimilarity") %>%
  group_by(method, tbin) %>%
  do(mod = t.test(dist ~ data, .)) %>%
  mutate(pval = mod$p.value,
         sig = pval < 0.05,
         symb = ifelse(pval < 0.001, '***',
                       ifelse(pval < 0.01, '**',
                              ifelse(pval < 0.05, '*', ''))))

j <- j %>% left_join(statz, by = c("method", "tbin"))

tmp <- j %>%
  group_by(method, data, tbin, sig, symb) %>%
  summarise(dist = mean(dist)) %>%
  group_by(method, tbin) %>%
  mutate(symb = ifelse(sig & dist == max(dist), symb, NA),
         sig = ifelse(is.na(symb), NA, sig))

figInterbabyDiss <- ggplot(j, aes(x = tbin, y = dist, lty = data, shape = data)) +
  stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange") +
  geom_line(data = tmp) +
  geom_text(aes(label = symb, y = dist + 0.4), data = filter(tmp, sig), size = 6) +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_linetype_discrete(name = '') +
  scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
  expand_limits(y = c(0)) +
  facet_wrap(~method, scales = 'free_y', ncol = 3) +
  labs(x = 'Months after birth', y = 'Dissimilarity among infants') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom')

figInterbabyDissTreat <- function(){}

subs_cs <- meta$subject[meta$treatment_group == 'C-section']
subs_ab <- meta$subject[meta$treatment_group == 'Antibiotics']
subs_control <- meta$subject[meta$treatment_group == 'Control']

j <- interdist %>%
  mutate(
    delivery1 = meta$delivery[match(subject1, meta$subject)],
    delivery2 = meta$delivery[match(subject2, meta$subject)],
    treat1 = meta$treatment_group[match(subject1, meta$subject)],
    treat2 = meta$treatment_group[match(subject2, meta$subject)],
    delivery = ifelse(delivery1 == delivery2, delivery1, NA),
    treat = ifelse(treat1 == treat2, treat1, NA),
    tbin = ceiling(tbin / 6) * 6) %>%
  filter(!delivery1 != delivery2 & !treat1 != treat2) %>%
  mutate(
    group = ifelse(subject1 %in% subs_cs & subject2 %in% subs_cs, 'C-section',
     ifelse(subject1 %in% subs_ab & subject2 %in% subs_ab, 'Antibiotics',
      ifelse(subject1 %in% subs_control & subject2 %in% subs_control, 'Control', NA)))) %>%
  filter(!is.na(group)) %>%
  ungroup() %>%
  select(method, group, tbin, dist)

statz <- bind_rows(filter(j, group != 'C-section') %>% mutate(group2 = 'Antibiotics'),
                   filter(j, group != 'Antibiotics') %>% mutate(group2 = 'C-section')) %>%
  group_by(method, group2, tbin) %>%
  do(mod = summary(lm(dist ~ group, data = .))) %>%
  mutate(pval = mod$coefficients[[8]] * 36, 
         mod = NULL)

j <- j %>%
  left_join(statz, by = c("method", "tbin")) %>%
  mutate(
    sig = group != 'Control' & pval < 0.05,
    group = factor(group, c("Control", "Antibiotics","C-section")),
    tbin0 = tbin,
    tbin = ifelse(group == 'Antibiotics', tbin - 0.5, tbin),
    tbin = ifelse(group == 'C-section', tbin + 0.5, tbin)) %>%
  group_by(method, group, tbin) %>%
  mutate(mean = mean(dist),
         yval = ifelse(method == 'OTU-based dissimilarity',
                  ifelse(group == 'Antibiotics', 0.85, 0.9),
                    ifelse(group == 'Antibiotics', 6.5, 7)),
         symb = ifelse(group == 'Control', NA,
                  ifelse(pval < 0.001, '***',
                    ifelse(pval < 0.01, '**',
                      ifelse(pval < 0.05, '*', '')))))

figInterbabyDissTreat <- ggplot(j, aes(x = tbin, y = dist, color = group)) +
    stat_summary(fun.data = "mean_cl_boot") +
    geom_line(aes(group = group, y = mean)) +
    #geom_text(aes(label = symb, y = yval, x = tbin0), show.legend = FALSE, data = filter(j, !is.na(symb))) +
    scale_shape_discrete(guide = FALSE) +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    expand_limits(y = 0) +
    facet_wrap(~method, scales = 'free_y', ncol = 3) +
    scale_color_manual(values = c('black',"#00A480","#FF6200"), name = '') +
    labs(x = 'Months after birth', y = 'Dissimilarity among infants') +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      legend.position = 'bottom')


#### Supplementary data ####
traits_sparse %>% write.csv(paste0(wd, 'SupplementaryData_1_trait_value_observations.csv'))
traits %>% spread(trait, val) %>% write.csv(paste0(wd, 'SupplementaryData_2_trait_value_predictions.csv'))

#### Stats and facts ####
#N for pairwise trait comparisons by trait
tree_Npairs <- trait_deltas %>% mutate(trait = trait_names[match(trait, names(trait_names))]) %>% group_by(trait) %>% summarise(n = length(dist_bin)) %>% filter(n == min(n) | n == max(n)) %>% arrange(n)

#number of OTUs in tree that are ours (before rarefaction)
Notus_ours <- sum(grepl('OTU', tree$tip.label))

#number of OTUs in tree that are from LTP
Notus_LTP <- sum(!grepl('OTU', tree$tip.label))

#number of c-section infants
N_cs <- sum(meta$delivery == 'caesaren')

#number of antibiotic-treated infants
N_ab <- sum(meta$treatment_group == 'Antibiotics' & meta$delivery != 'caesaren')

#number of control infants not treated with antibiotics or born via c-section
N_control <- sum(meta$treatment_group == 'Control' & meta$delivery != 'caesaren' & meta$antibiotic_days == 0)

#number of samples in each 6 month bin, by treatment
tmp <- otus_cs %>% filter(t > 1 & t < 36) %>% left_join(meta[, c('subject','treatment_group')], by = 'subject') %>% mutate(treat = ifelse(delivery == 'caesaren', 'caesaren', treatment_group)) %>% group_by(treat, t %/% 6) %>% summarise(n = length(unique(sampleID))) %>% group_by(treat) %>% summarise(min = min(n), max = max(n))
n_cs_min <- tmp$min[tmp$treat == 'caesaren']
n_cs_max <- tmp$max[tmp$treat == 'caesaren']
n_ab_min <- tmp$min[tmp$treat == 'Antibiotics']
n_ab_max <- tmp$max[tmp$treat == 'Antibiotics']
n_cntrl_min <- tmp$min[tmp$treat == 'Control']
n_cntrl_max <- tmp$max[tmp$treat == 'Control']


#number of samples per 1.month bin in cwm calculations
tmp <- otus %>% mutate(t = round(t)) %>% filter(t < 37 & t > 1) %>% group_by(t) %>% summarise(n = length(unique(sampleID))) %>% ungroup() %>% summarise(min = min(n), max = max(n))
n_cwms_min <- tmp %>% pull(min)
n_cwms_max <- tmp %>% pull(max)

#number of samples per 6.month bin in intrababy calculations
tmp <- otus %>% mutate(tbin = c(t %/% 6) * 6) %>% filter(tbin < 36) %>% group_by(subject, tbin) %>% summarise(n = length(unique(sampleID)) - 1) %>% group_by(tbin) %>% summarise(n = sum(n)) %>% ungroup() %>% summarise(min = min(n), max = max(n))
n_intrababy_min <- tmp %>% pull(min)
n_intrababy_max <- tmp %>% pull(max)

#number of samples per infant
sampstats <- otus_cs %>% group_by(subject) %>% summarise(n = length(unique(t))) %>% ungroup() %>% summarise(min = min(n), max = max(n), mean = mean(n), sd = sd(n), median = median(n)) %>% mutate_all(funs(round), 2)

#number of data in sparse traits
sparse_tdats <- traits_sparse %>% select(-Genus, -Species) %>% gather(trait, val) %>% filter(!is.na(val)) %>% summarise(n = length(val)) %>% pull(n)

#correlation of gene number and genome size
gene_genome_corr <- round(cor(traits_sparse$Gene_number, traits_sparse$Genome_Mb, use = 'na.or.complete'),3)

#coverage of trait data by trait
tc <- trait_coverage %>% filter(Source != 'Missing data') %>% group_by(trait) %>% summarise(coverage = round(100 * sum(abun), 1)) %>% arrange(coverage) %>% filter(coverage > 50) %>% filter(coverage %in% c(min(coverage), max(coverage)))

#plot for elena and ashley breaking down the relationship between c-section and antibiotcs across subjects
meta %>% 
  arrange(delivery) %>% 
  mutate(subject = factor(subject, levels = unique(subject))) %>% 
  arrange(antibiotic_days) %>% 
  mutate(
    abxdays = ifelse(antibiotic_days < 15, paste0(antibiotic_days, '+'), 
                     as.character(antibiotic_days)), 
    abxdays = factor(abxdays, levels = unique(abxdays))) %>% 
  ggplot(aes(x = abxdays, fill = delivery)) + 
    geom_bar(color = 'black', aes(group = subject)) + 
    th + 
    labs(x = 'Days of antibiotic treatment', y = 'Number of infants')

