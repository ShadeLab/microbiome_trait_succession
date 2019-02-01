#### Setup ####

#setwd
setwd('~\\msu\\microbiome_trait_succession\\')

#load custom functions
source('succession_custom_functions.R')

#reprocess data?
if(FALSE) {
  
  #process metadata
  source('succession_processing_metadata.R')
  
  #process trait data (to re-download, refer to script)
  source('succession_processing_trait_data.R')
  
  #generate/infer trait predictions for unknown trait data
  source('succession_trait_predictions.R')
  
}

#load packages
loadpax(pkg = c('Hmisc','gtable','dplyr','tidyr','knitr','data.table','vegan',
                'GGally','grid','gridExtra','betapart','scales','stringr', 'phyloseq',
                'ape','kableExtra','poilog','ggtree','ggplot2','cowplot'))

#load data
otus <- readRDS(file = 'data\\otus.RDS')
meta <- readRDS(file = 'data\\subject_metadata.RDS')
tax_succ <- readRDS(file = 'data\\succession_tax_SILVA.RDS')
tax_LTP <- readRDS(file = 'data\\tax_LTP.RDS')
trait_deltas <- readRDS(file = 'data\\trait_deltas.RDS')
trait_delta_models <- readRDS(file = 'data\\trait_delta_models.RDS')
traits <- readRDS(file = 'data\\traits.RDS')
traits_sparse <- readRDS(file = 'data\\traits_sparse.RDS')
trait_predictions <- readRDS(file = 'data\\traits_all_predictions.RDS')
trait_sources <- readRDS(file = 'data\\trait_sources.RDS')
tree <- read.tree(file = 'data\\LTP_succession.tree')

#wide version of traits
traits_wide <- traits %>% spread(trait, val)

#trait renaming dictionary
trait_names <- c(
  "Aggregation_score" = "Aggregation score",
  "B_vitamins"        = "B vitamins",
  "Copies_16S"        = "16S gene copies",
  "GC_content"        = "GC content",
  "Gene_number"       = "Gene number",
  "Genome_Mb"         = "Genome size",
  "Gram_positive"     = "Gram-positive",
  "IgA"               = "IgA binding affinity",
  "Length"            = "Length",
  "Motility"          = "Motility",
  "Oxygen_tolerance"  = "Oxygen tolerance",
  "pH_optimum"        = "pH optimum",
  "Salt_optimum"      = "Salt optimum",
  "Sporulation"       = "Sporulation score",
  "Temp_optimum"      = "Temperature optimum",
  "Width"              = "Width"
)

trait_units <- data.frame(stringsAsFactors = FALSE,
  "Aggregation_score" = "0 (never) to 1 (observed aggregation)",
  "B_vitamins"        = "No. B-vitamin pathways in genome",
  "Copies_16S"        = "No. in 16S rRNA gene copies in genome",
  "GC_content"        = "Percent (\\%) guanine and cytosine in genome",
  "Gene_number"       = "No. genes in genome",
  "Genome_Mb"         = "Genome size in megabases",
  "Gram_positive"     = "0 (Gram-negative) to 1 (Gram-positive)",
  "IgA"               = "log ([IgA+]/[IgA-] + 1)",  
  "Length"            = "log ($\\mu$m)",
  "Motility"          = "0 (never motile) to 1 (always motile)",  
  "Oxygen_tolerance"  = "0 (obligate anaerobe) to 5 (obligate aerobe)",
  "pH_optimum"        = "pH",    
  "Salt_optimum"      = "g per l",  
  "Sporulation"       = "0 (never sporulates) to 1 (sporulates easily)",  
  "Temp_optimum"      = "$^{\\circ}$C",  
  "Width"             = "log ($\\mu$m)"
)

trait_names_units <- c(
  "Aggregation_score" = "Aggregation score",
  "B_vitamins"        = "B vitamins (#)",
  "Copies_16S"        = "16S gene copies (#)",
  "GC_content"        = "GC content (%)",
  "Gene_number"       = "Gene number (#)",
  "Genome_Mb"         = "Genome size (Mb)",
  "Gram_positive"     = "Gram-positive (%)",
  "IgA"               = "IgA binding score",
  "Length"            = "Length",
  "Motility"          = "Motility (%)",
  "Oxygen_tolerance"  = "O2 tolerance score", 
  "pH_optimum"        = "pH optimum",
  "Salt_optimum"      = "Salt optimum",
  "Sporulation"       = "Sporulation score",
  "Temp_optimum"      = "Temp. optimum",
  "Width"             = "Width"
)

#### Trait sources ####

##Table of trait sources
tabTraitSources <- function(){}

tabTraitSources <- trait_units %>%
  gather(Trait, Units) %>%
  left_join(trait_sources, by = 'Trait') %>%
  mutate(Trait = trait_names[match(Trait, names(trait_names))]) %>%
  mutate(Sources = gsub("ú", "\\'{u}", fixed = TRUE, Sources),
         Sources = gsub("ó", "\\'{o}", fixed = TRUE, Sources))

## Figure of trait correlations
figTraitCorr <- function(){}

j <- expand.grid(unique(traits$trait), unique(traits$trait)) %>%
  group_by(Var1, Var2) %>%
  do(ct = cor.test(traits_wide[[as.character(.$Var1)]], traits_wide[[as.character(.$Var2)]])) %>%
  mutate(Var1 = factor(trait_names[match(Var1, names(trait_names))], levels = trait_names),
         Var2 = factor(trait_names[match(Var2, names(trait_names))], levels = trait_names),
         cor = ct$estimate[[1]],
         pval = ct$p.value,
         pval = p.adjust(pval, method = 'BH', n = 55),
         df = ct$parameter,
         symb = ifelse(pval < 0.001, '***',
                       ifelse(pval < 0.01, '**',
                              ifelse(pval < 0.05, '*', ''))))

#convert to numeric to eliminate redundancies
j <- filter(j, as.numeric(Var1) > as.numeric(Var2))
j$Var1 <- factor(j$Var1, levels = rev(levels(j$Var1)))

figTraitCorr <- j %>%
  ggplot(aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = 'black') +
  geom_text(aes(label = symb)) +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', name = 'Pears. Corr.') +
  labs(x = '', y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    legend.position = c(1, 1), 
    legend.justification = c(1, 1))

#### Coverage ####

##A figure showing sampling timing by subject; very reminiscent of Yassour et al. 2016
figSampling <- function(){}

figSampling <- otus %>%
  left_join(meta[, c('subject','delivery')], by = 'subject') %>%
  distinct(subject, t, delivery) %>% 
  group_by(subject) %>%
  mutate(t_min = min(t), t_max = max(t), len = t_max - t_min) %>%
  ungroup() %>%
  arrange(desc(len)) %>%
  mutate(
    delivery = ifelse(delivery == 'caesaren', 'C-section delivery','Vaginal delivery'),
    subject = as.factor(subject),
    subject_num = as.numeric(subject)) %>%
  ggplot(aes(y = subject, fill = delivery, color = delivery)) + 
    geom_point(aes(x = t, y = subject), shape = 21) + 
    geom_segment(aes(x = t_min, xend = t_max, y = subject_num, yend = subject_num)) + 
    scale_fill_manual(values = c('black','grey60'), name = '') +
    scale_color_manual(values = c('black','grey60'), name = '') +
    theme_classic() + 
    theme(legend.position = 'bottom',
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    labs(x = 'Months after birth', y = 'Subject')

## Figure showing data coverage (in terms of the proportion of total abundance) by trait
figTraitCoverage <- function(){}

#create data with 'directly observed' traits and their associated Genus species labels
ts <- traits_sparse %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val))

#for each Genus-Species combination... and its abundance...
#do we have directly observed trait data?
#do we have inferred trait data? If not, it is missing data.
trait_coverage <- otus %>%
  group_by(otu) %>%
  summarise(abun = sum(abun)) %>%
  left_join(tax_succ[, c('otu', 'Genus', 'Species')], by = 'otu') %>%
  left_join(traits_sparse[, c('Genus','Species','Aggregation_score','pH_optimum','Salt_optimum','IgA')], by = c('Genus','Species')) %>%
  left_join(traits_wide, by = 'otu') %>%
  gather(trait, val, -otu, -abun, -Genus, -Species) %>%
  mutate(Source = ifelse(paste(trait, Genus, Species) %in% paste(ts$trait, ts$Genus, ts$Species),
                         'Observed data',
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

#plot
figTraitCoverage <- ggplot(trait_coverage, aes(x = trait, y = abun, fill = Source)) + 
  geom_bar(color = 'black', stat = 'identity') +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c('lightgrey','skyblue','blue'), name = '') +
  labs(x = "", y = 'Proportion of total abundance') +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom',
    plot.margin = unit(c(12, 12, 12, 12), "pt")) +
  guides(fill = guide_legend(reverse = TRUE))

## Tree with all taxa from LTP and in this study, with associated trait data
figTreeTraits <- function(){}

#load tax, create Genus species binomial
tax_tmp <- bind_rows(tax_LTP, tax_succ[, c('otu','Genus','Species')]) %>%
  mutate(spp = paste(Genus, Species)) %>%
  select(-Genus, -Species)

# gather sparse traits and Genus species binomial
tmp <- traits_sparse %>% 
  mutate(spp = paste(Genus, Species)) %>%
  select(-Genus, -Species) %>%
  gather(trait, val, -spp)

# mamke two lists of (1) OTU names in this study and (2) OTU names from LTP
tmp2 <- list(
  `OTU from this study` = tree$tip.label[grepl('OTU', tree$tip.label)],
  `OTU from the Living Tree Project` = tree$tip.label[!grepl('OTU', tree$tip.label)])

# add these descriptors to the tree (for ggtree to be able to understand them)
tree <- groupOTU(tree, tmp2)

#create a vector of numbers that roughly correspond to taxa with trait data, to map onto tree
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

#create ggtree
#note i originally modified the ggtree internal function add_panel (used in facet_plot) to have the panel label 'a', rather than 'Tree'. Such a useful hack! But now i just rename the tree panel data frame.
p <- ggtree(tree, size = 0.05) +
  geom_tiplab(label = paste0(rep('_', 52), collapse = ''), vjust = -.29, offset = 0.01, 
              size = 1, aes(color = group))

#rename tree panel to facet label 'a'
p$data[[".panel"]] <- "a"

#add facet plot to ggtree
figTreeTraits <- facet_plot(p, panel = "b", data = j, geom = geom_point, 
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
        legend.title = element_text(color = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 'bold'))

figTreeTraits

#### HSP comparisons ####
## Testing the differences among different methods of hidden state character predictions
tabHSP <- function(){}

#note that by using trait predictions, we only consider those traits that were amenable at all to hidden state prediction.
#And instances where dist == 0 were also removed because these are "observations" not inferences
tabHSP <- trait_predictions %>%
  filter(trait %in% traits$trait) %>%
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

##A table with the maximum allowable distances used to infer triat values, for each trait
#these were calculated in the succession_trait_predictions.R script
tabThresholds <- function(){}

max_dists <- trait_delta_models %>%
  distinct(trait, max_dist, .keep_all = FALSE) %>%
  mutate(max_dist = format(round(max_dist, 3), nsmall = 3))

tabThresholds <- max_dists %>%
  mutate(trait = trait_names[match(trait, names(trait_names))]) %>%
  rename(Trait = trait, `Max. distance` = max_dist) %>%
  as.data.frame()

## A figure illustrating the process of how we determined the maximum thresholds (above)
figThresholds <- function(){}

#First, calculate null expectations of trait-based differences among taxa with no phylogenetic signal 
# (i.e., taxa that more distantly related than the maximum allowable phylogenetic distance (see above)
nulls <- trait_deltas %>%
  left_join(max_dists) %>%
  group_by(trait) %>%
  filter(dist > max_dist) %>%
  summarise(null = mean(delta))

#load all trait deltas (pairwise trait-based differences among taxa)
#calculate mean trait-based differences within each bin for each trait
#only consider bins below 20% different in 16S V4 rRNA region -- beyond that things get sparse and don't bin well
delta_bins <- trait_deltas %>%
  filter(dist_bin <= 0.2) %>%
  mutate(dist_bin = dist %/% 0.005 * 0.005) %>% 
  group_by(trait, dist_bin) %>%
  summarise(
    mean = mean(delta),
    null = nulls$null[match(trait[1], nulls$trait)]) %>%
  ungroup()

#plot with standard deviations
#remove a few outliers that **importantly** fall well beyond the maximum allowable differences, so aren't significant anyway
j <- delta_bins %>%
  filter(!(trait == 'Sporulation' & mean > 0.2)) %>%
  filter(!(trait == 'Length' & mean > 1.5)) %>%
  filter(!(trait == 'Width' & mean > 1.2)) %>%
  left_join(distinct(trait_delta_models[, c('trait','type')]), by = 'trait') %>%
  mutate(type_lev = as.numeric(factor(type)),
         color = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", NA)[type_lev],
         trait = trait_names_units[match(trait, names(trait_names_units))])

tdms <- mutate(trait_delta_models, trait = trait_names_units[match(trait, names(trait_names_units))])

myplots <- list()
ts <- sort(unique(j$trait))
for (i in ts) {
  
  mycol <- j %>% filter(trait == i) %>% summarise(color = color[1]) %>% pull(color)

  myplots[[match(i, ts)]] <- j %>%
    filter(trait %in% i) %>%
    ggplot(aes(x = dist_bin, y = mean)) +
      geom_point() +
      geom_line(aes(y = delta), color = mycol, data = filter(tdms, max_dist > 0.03 & trait == i), lwd = 1.5) +
      geom_hline(aes(yintercept = null), linetype = 3) +
      geom_vline(aes(xintercept = max_dist), data = filter(tdms, trait == i), lty = 2) +
      theme_classic() +
      scale_x_continuous(breaks = c(0,0.08,0.16)) +
      labs(x = '', y = bquote(paste(Delta, .(i)))) +
      theme(
        strip.background = element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(0, 5.5, 0, 0), "pt"),
        axis.title.x=element_blank()) +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
  
}

myleg <- j %>%
  ggplot(aes(x = dist_bin, y = mean)) +
    geom_line(aes(y = delta, color = type), data = filter(tdms, max_dist > 0.03), lwd = 1.5) +
  facet_wrap(~trait, ncol = 3) +
    scale_color_manual(name = '', values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  theme(legend.direction = "horizontal",
        legend.justification="center" ,
        legend.box.just = "bottom") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

myleg <- get_legend(myleg)

myplots[1:13] <- lapply(myplots[1:13], function(x) {
  x + theme(axis.text.x = element_blank())
})

myplots[[9]] <- myplots[[9]] + labs(y = expression(Delta ~ Length ~ "(log["*mu*"m])"))
myplots[[13]] <- myplots[[13]] + labs(y = expression(Delta ~ Salt ~ optimum ~ (g ~ l^{-1})))
myplots[[15]] <- myplots[[15]] + labs(y = expression(Delta~'Temp.'~'optimum'~'('*degree*'C)'))
myplots[[16]] <- myplots[[16]] + labs(y = expression(Delta~Width~'(log['*mu*'m])'))

ps <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]], myplots[[9]], myplots[[10]], myplots[[11]], myplots[[12]], myplots[[13]], myplots[[14]], myplots[[15]], myplots[[16]], ncol = 4, align = 'hv')

xtitle <- ggdraw() + draw_label("Phylogenetic distance between tips")

figThresholds <- plot_grid(ps, xtitle, myleg, ncol = 1, rel_heights = c(1,0.06, 0.06))


#### Taxonomic patterns ####
## figure that (1) groups taxa in to early/mid/late successional specialists,
# and then (2) shows their abundance patterns over time
figTaxOverTime <- function(){}

#generate zeroes by spread(..., fill = 0) and then regathering
#filter out zeroes only when taxa *never* appear in subjects
#determine overall trends over succession using lm() and pvals and tvals
mods <- otus %>%
  spread(otu, abun, fill = 0) %>%
  gather(otu, abun, -sampleID, -subject, -t) %>%
  group_by(otu, subject) %>%
  filter(sum(abun) > 0) %>%
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
    filter(group == 'Early successional') %>%
    ggplot(aes(x = t, y = abun)) +
      geom_bar(stat = 'identity', width = 1, fill = '#66c2a5', color = '#66c2a5') +
      expand_limits(y = 0.7) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = '', y = 'Relative abundance', tag = 'a') +
      ggtitle('Early successional') +
      theme_classic() +
      theme(legend.position = 'none',
            axis.text.x = element_blank(),
            plot.margin = unit(c(0, 0, 0, 5.5), "pt"))
    
p3 <- j %>%
  group_by(t) %>%
  mutate(abun = abun / sum(abun)) %>%
  filter(group == 'Mid-successional / No trend') %>%
  ggplot(aes(x = t, y = abun)) +
  geom_bar(stat = 'identity', width = 1, fill = '#fc8d62', color = '#fc8d62') +
  expand_limits(y = 0.7) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Relative abundance', tag = 'c') +
  ggtitle('Mid-successional / No trend') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 5.5), "pt"))


p5 <- j %>%
  group_by(t) %>%
  mutate(abun = abun / sum(abun)) %>%
  filter(group == 'Late successional') %>%
  ggplot(aes(x = t, y = abun)) +
  geom_bar(stat = 'identity', width = 1, fill = '#8da0cb', color = '#8da0cb') +
  expand_limits(y = 0.7) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Months after birth', y = 'Relative abundance', tag = 'e') +
  ggtitle('Late successional') +
  theme_classic() +
  theme(legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 5.5), "pt"))

j <- j %>%
  left_join(tax_succ, by = 'otu') %>%
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
  summarise(abun = sum(abun))

p2 <- j %>%
  filter(group == 'Early successional') %>%
  ggplot(aes(x = Family, y = abun)) +
    geom_bar(stat = 'identity', color = 'black', fill = '#66c2a5') +
    geom_text(aes(label = Family), hjust = 'left', nudge_y = 0.005) +
    labs(x = "", y = "", tag = 'b   ') +
    ggtitle("Early successional") +
    scale_y_continuous(limits = c(0, 0.325)) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.text=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = 'none',
      axis.title.y=element_blank(),
      plot.margin = unit(c(0, 5.5, 0, 0), "pt"))


p4 <- j %>%
  filter(group == 'Mid-successional / No trend') %>%
  ggplot(aes(x = Family, y = abun)) +
  geom_bar(stat = 'identity', color = 'black', fill = '#fc8d62') +
  geom_text(aes(label = Family), hjust = 'left', nudge_y = 0.005) +
  labs(x = "", y = "", tag = 'd   ') +
  ggtitle('Mid-successional / No trend') +
  scale_y_continuous(limits = c(0, 0.325)) +
  coord_flip() +
  theme_classic() + 
  theme(
    axis.text=element_blank(),
    axis.ticks.y=element_blank(),
    legend.position = 'none',
    axis.title.y=element_blank(),
    plot.margin = unit(c(0, 5.5, 0, 0), "pt"))


p6 <- j %>%
  filter(group == 'Early successional') %>%
  ggplot(aes(x = Family, y = abun)) +
  geom_bar(stat = 'identity', color = 'black', fill = '#8da0cb') +
  geom_text(aes(label = Family), hjust = 'left', nudge_y = 0.005) +
  labs(x = "", y = "Relative abundance across all samples", tag = 'f   ') +
  ggtitle("Late successional") +
  scale_y_continuous(limits = c(0, 0.325)) +
  coord_flip() +
  theme_classic() + 
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    legend.position = 'none',
    axis.title.y=element_blank(),
    plot.margin = unit(c(0, 5.5, 0, 0), "pt"))

figTaxOverTime <- plot_grid(p1, p2, p3, p4, p5, p6, align = 'h', ncol = 2)

#### Traits over time ####

##A figure of trait CWMs over time
figCWMs <- function(){}

#note: We calculate CWM for each sample, but plot averages of all samples within months
j <- otus %>%
  left_join(traits_wide, by = 'otu') %>%
  gather(trait, val, -sampleID, -subject, -t, -otu, -abun) %>%
  filter(!is.na(val)) %>%
  mutate(trait = trait_names_units[match(trait, names(trait_names_units))]) %>%
  left_join(meta, by = 'subject') %>%
  group_by(subject, t, trait) %>%
  summarise(cwm = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(t = round(t)) %>%
  filter(t < 37 & t > 1)

myplots <- list()
ts <- sort(unique(j$trait))
for (i in ts) {
  
  myplots[[match(i, ts)]] <- j %>%
    filter(trait %in% i) %>%
    ggplot(aes(x = t, y = cwm)) +
      stat_summary(fun.y = mean, color = 'red', geom = "point") + 
      stat_summary(fun.data = mean_cl_boot, color = 'red', geom = "linerange") +
      stat_smooth(color = 'black', se = FALSE) +
      scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
      labs(x = "", y = i, tag = paste(letters[match(i, ts)], ' ')) +
      theme_classic() +
      theme(
        strip.background = element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.title.x=element_blank()) +
      annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
      annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

  
}

myplots[1:9] <- lapply(myplots[1:9], function(x) {
  x + theme(axis.text.x = element_blank())
})

myplots[[7]] <- myplots[[7]] + labs(y = expression('Length'~'(log['*mu*'m])'))
myplots[[11]] <- myplots[[11]] + labs(y = expression('Temp.'~'optimum'~'('*degree*'C)'))
myplots[[12]] <- myplots[[12]] + labs(y = expression('Width'~'(log['*mu*'m])'))

ps <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]], myplots[[9]], myplots[[10]], myplots[[11]], myplots[[12]], ncol = 3, align = 'hv')

xtitle <- ggdraw() + draw_label("Months after birth", size = 12)

figCWMs <- plot_grid(ps, xtitle, ncol = 1, rel_heights = c(1,0.06))

## figure showing CWMs separated by treatment
figCWMsTreat <- function(){}

#prep data for delivery mode, antibiotics, control groups
tmp <- otus %>% 
  left_join(meta[, c('subject','treatment_group')], by = "subject") %>%
  rename(treatment = treatment_group)
tmp1 <- filter(tmp, treatment == 'C-section')
tmp2 <- filter(tmp, treatment == 'Antibiotics')
tmp3 <- filter(tmp, treatment == 'Control')

#append controls to two treatments for pairwise statistical testing purposes
tmp1 <- bind_rows(tmp1, tmp3) %>% mutate(group = 'Delivery mode')
tmp2 <- bind_rows(tmp2, tmp3) %>% mutate(group = 'Antibiotics')

#join traits, calculate CWMs for each sample
#note: i bin samples within half-years (0-6M, 6-12M, 12-18M, etc.) because 
# (1) performing t-tests within each month didn't have sufficient sample sizes and (2) it didn't seem meaningful.
cwmsx <- bind_rows(tmp1, tmp2) %>%
  left_join(traits, by = 'otu') %>%
  filter(!is.na(trait)) %>%
  group_by(group, treatment, sampleID, subject, t, trait) %>%
  summarise(val = weighted.mean(val, w = abun)) %>%
  ungroup() %>%
  filter(round(t) > 1 & round(t) < 36) %>%
  mutate(trait = trait_names_units[match(trait, names(trait_names_units))], 
         treatment = factor(treatment, levels = c('Control','Antibiotics','C-section')),
         tbin = ceiling(t/6) * 6 - 3)

#perform ttests with 6 mo periods (tbin), 
#option: adjust pval for false positives duemultiple comparisons
statz <- cwmsx %>%
  group_by(trait, group, tbin) %>%
  do(mod = t.test(val ~ treatment, data = .)) %>%
  mutate(pval = mod$p.value,
         pval = p.adjust(pval, method = 'BH', n = 12))

#add pvals to cwm data
#add custom tbin for staggered plotting
tmp <- cwmsx %>%
  filter(!(group == 'Antibiotics' & treatment == 'Control')) %>%
  mutate(pval = statz$pval[match(paste(trait, group, tbin), 
                                 paste(statz$trait, statz$group, statz$tbin))],
         sig = ifelse(treatment != 'Control', pval < 0.1, FALSE),
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

myplots <- list()
ts <- sort(unique(tmp$trait))
for (i in ts) {
  
  j1 <- tmp %>% filter(trait %in% i)
  j2 <- tmp2 %>% filter(trait %in% i)
  j3 <- tmp3 %>% filter(trait %in% i)
  
  myplots[[match(i,ts)]] <- ggplot(j1, aes(x = tbin, y = val, color = treatment)) +
    stat_summary(fun.data = "mean_cl_boot", size = .25) +
    geom_line(data = j2) +
    scale_shape_discrete(guide = FALSE) +
    geom_text(aes(label = symb), data = j3, show.legend = FALSE, size = 6) +
    geom_point(aes(y = yval), data = j3, alpha = 0) +
    scale_color_manual(values = c('black',"#00A480","#FF6200"), name = '') +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    expand_limits(x = 2) +
    labs(x = '', y = i,
         tag = paste(letters[match(i, ts)], ' ')) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.title.x = element_blank(),
          legend.position = 'bottom',
          plot.margin = unit(c(0, 0, 0, 0), "pt")) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

}

myleg <- get_legend(myplots[[1]])

myplots <- lapply(myplots, function(x) x + theme(legend.position = 'none'))

myplots[1:9] <- lapply(myplots[1:9], function(x) {
  x + theme(axis.text.x = element_blank(), axis.title.x=element_blank())
})


myplots[[7]] <- myplots[[7]] + labs(y = expression('Length'~'(log['*mu*'m])'))
myplots[[11]] <- myplots[[11]] + labs(y = expression('Temp.'~'optimum'~'('*degree*'C)'))
myplots[[12]] <- myplots[[12]] + labs(y = expression('Width'~'(log['*mu*'m])'))

ps <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]], myplots[[9]], myplots[[10]], myplots[[11]], myplots[[12]], ncol = 3, align = 'hv')

xtitle <- ggdraw() + draw_label("Months after birth", size = 12)

figCWMsTreat <- plot_grid(ps, xtitle, myleg, ncol = 1, rel_heights = c(1,0.06, 0.06))


## Figure showing the variance of trait values of *cells* (not species)
figCWMsVariance <- function(){}

#scale each trait, calculate variance of trait values of *cells* not species
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
  mutate(trait = trait_names[match(trait, names(trait_names))],
         trait = ifelse(trait == 'Temperature optimum', 'Temp. optimum', trait))

myplots <- list()
ts <- sort(unique(tmp$trait))
for (i in ts) {
  
  myplots[[i]] <- ggplot(tmp[tmp$trait == i, ], aes(x = t, y = var)) +
    stat_summary(fun.y = mean, geom = "point", color = 'darkcyan') + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange", color = 'darkcyan') +
    stat_summary(fun.data = mean_cl_boot, geom = "linerange", color = 'darkcyan', data = tmp,
                 aes(group = trait), alpha = 0) +
    stat_smooth(se = FALSE, color = 'black') +
    scale_x_continuous(breaks= seq(10,30,by=10)) +
    labs(x = '', y = '', tag = letters[match(i, ts)]) +
    ggtitle(i) +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      strip.background = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  
}

myplots[1:9] <- lapply(myplots[1:9], function(x) {
  x + theme(axis.text.x = element_blank())
})

ps <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]], myplots[[9]], myplots[[10]], myplots[[11]], myplots[[12]], ncol = 3, align = 'hv')

xtitle <- ggdraw() + draw_label("Months after birth")
ytitle <- ggdraw() + draw_label(" Abundance-weighted community variance", angle = 90)

figCWMsVariance <- plot_grid(ps, xtitle, ncol = 1, rel_heights = c(1,0.06))
figCWMsVariance <- plot_grid(ytitle, figCWMsVariance, ncol = 2, rel_widths = c(0.06,1))


## Figure showing trait-based differences among early/mid/late succession specialists
figCWMsSuccGroup <- function(){}

labs <- c(Early = 'Early successional',
          Other = 'Mid-successional / No trend', 
          Late = 'Late successional')

#join trait values to successional specialities
j <- otus %>%
  left_join(mods, by = 'otu') %>%
  group_by(group, t = round(t)) %>%
  filter(t >= 2 & t <= 36) %>%
  ungroup()

tmp <- traits %>%
  left_join(distinct(j, otu, group), by = 'otu') %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(names(labs)[match(group, labs)], names(labs)))

#Three sets of stats for the three possible combinations of early/mid/late comparisons
statz1 <- tmp %>%
  filter(group %in% c('Early', 'Other')) %>%
  group_by(trait) %>%
  do(mod = t.test(val ~ group, data = .)) %>%
  mutate(pval = mod$p.value, group = '1/2')

statz2 <- tmp %>%
  filter(group %in% c('Other', 'Late')) %>%
  group_by(trait) %>%
  do(mod = t.test(val ~ group, data = .)) %>%
  mutate(pval = mod$p.value, group = '2/3')

statz3 <- tmp %>%
  filter(group %in% c('Early', 'Late')) %>%
  group_by(trait) %>%
  do(mod = t.test(val ~ group, data = .)) %>%
  mutate(pval = mod$p.value, group = '1/3')

#combine
statz <- bind_rows(statz1, statz2, statz3) %>%
  mutate(pval = p.adjust(pval, method = 'BH', n = 12)) %>%
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

#create vector of significance letters
tmp$label = statz$label[match(paste(tmp$trait, tmp$group), paste(statz$trait, statz$group))]

#create vector for proper y-value locations for significance letters
yvals <- c(
  B_vitamins = 4.25,
  Copies_16S = 5.55,
  GC_content = 51,
  Gene_number = 3550,
  Genome_Mb = 3.6,
  Gram_positive = 0.65,
  Length = 1.05,
  Motility = 0.375,
  Oxygen_tolerance = 2.4,
  Sporulation = 0.25,
  Temp_optimum = 39.1,
  Width = -0.31
)

#add yvals to dataframe
tmp <- tmp %>% mutate(yval = yvals[match(trait, names(yvals))])

myplots <- list()
ts <- sort(unique(tmp$trait))
for (i in ts) {
  
  ii <- trait_names_units[match(i, names(trait_names_units))]
  j1 <- tmp %>% filter(trait %in% i)
  myplots[[match(i, ts)]] <- ggplot(j1, aes(x = group, y = val)) +
    stat_summary(fun.data = "mean_cl_boot") +
    geom_text(aes(label = label, y = yval), size = 3, data = filter(j1, !is.na(label) & !is.na(yval))) +
    labs(x = '', y = ii) + 
    theme_classic() +
    theme(
      strip.background = element_blank(),
      legend.position = 'none',
      #plot.margin = unit(c(0, 0, 0, 0), "pt"),
      axis.title.x=element_blank()) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)

}

myplots[1:9] <- lapply(myplots[1:9], function(x) {
  x + theme(axis.text.x = element_blank())
})

myplots[[7]] <- myplots[[7]] + labs(y = expression('Length'~'(log['*mu*'m])'))
myplots[[11]] <- myplots[[11]] + labs(y = expression('Temp.'~'optimum'~'('*degree*'C)'))
myplots[[12]] <- myplots[[12]] + labs(y = expression('Width'~'(log['*mu*'m])'))

ps <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]], myplots[[9]], myplots[[10]], myplots[[11]], myplots[[12]], ncol = 3, align = 'hv')

xtitle <- ggdraw() + draw_label("OTU successional group")

figCWMsSuccGroup <- plot_grid(ps, xtitle, ncol = 1, rel_heights = c(1,0.06))

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
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom')

#### Dissimilarity over time ####
figIntrababyDiss <- function(){}

#use meltdist() to calculate dissimilarities among samples within subjects for each method
#first based on otu abundances
j1 <- otus %>%
  spread(otu, abun, fill = 0) %>%
  group_by(subject) %>%
  do(meltdist(x = .$sampleID, y = ., method = 'bray')) %>%
  mutate(method = 'OTU-based dissimilarity', data = 'obs')

#next functional dissimilarities; we calculate community weighted mean for each trait, 
#then calculate the euclidean distance among samples across (scaled) traits
j2 <- otus %>%
  left_join(traits, by = 'otu') %>%
  filter(!is.na(trait)) %>%
  group_by(sampleID, subject, trait, t) %>%
  summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
  group_by(trait) %>%
  mutate(val = scale(val)) %>%
  spread(trait, val) %>%
  group_by(subject) %>%
  do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
  mutate(method = 'Trait-based dissimilarity', data = 'obs')

#calculate null model predictions? With 1000 reps it takes ~30 min on my machine.
if (FALSE) {
  
  #empty list to store reps
  j3 <- list()
  
  #number of reps
  reps <- 1000
  for (rep in 1:reps) {
    
    #sample() shuffles values in a vector
    tmp <- otus %>%
      left_join(traits, by = 'otu') %>%
      group_by(trait) %>%
      mutate(val = sample(val)) %>%
      filter(!is.na(val)) %>%
      group_by(sampleID, subject, trait, t) %>%
      summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
      group_by(trait) %>%
      mutate(val = scale(val)) %>%
      spread(trait, val) %>%
      group_by(subject) %>%
      do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
      ungroup() %>%
      mutate(method = 'Trait-based dissimilarity', data = 'ran')
    
    j3[[rep]] <- tmp
    print(paste(rep, 'of', reps, 'done'))
    flush.console()
  
  if (rep == reps) {
    
    #save null simulation data
    bind_rows(j3) %>%
      group_by(subject, sample1, sample2, t1, t2, method, data) %>%
      summarise(dist = mean(dist)) %>%
      saveRDS(file = 'data\\null_intrababy_predictions.RDS')

}}}
  
j3 <- readRDS('data\\null_intrababy_predictions.RDS')

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

#determine difference, if any, between simulated and observed initial trait dissimilarity
#later this will be used to adjust simulated and observed data start together
tb_bump <- j5 %>%
  ungroup() %>%
  filter(method == 'Trait-based dissimilarity') %>%
  mutate(t1 = t1 %/% 6 * 6) %>%
  filter(t1 == min(t1)) %>%
  summarise(diff = mean(dist[data == 'obs']) - 
                   mean(dist[data == 'ran'])) %>%
  pull(diff)

#apply tb_bump
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
    method = factor(method),
    data = factor(c('Observed','Null model')[match(data, c('obs','ran'))],
                  c('Observed','Null model')))

#perform t.tests between simulated nulls and observed trait dissimilarity
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
                     ifelse(pval < 0.1, '', '')))))

#add significance values to data
j <- j %>% left_join(statz, by = c("method", "group", "tbin"))

#determine positioning for significance markers
tmp <- j %>% 
  group_by(tbin, method, group, data, sig, symb) %>% 
  summarise(dist = mean(dist)) %>%
  group_by(tbin, method, group) %>%
  mutate(symb = ifelse(sig & dist == max(dist), symb, NA),
         sig = ifelse(is.na(symb), NA, sig))


p1 <- j %>%
  filter(method == 'OTU-based dissimilarity' & group == 'Dissimilarity to next sample') %>%
  ggplot(aes(x = tbin, y = dist, lty = data)) +
    stat_summary(fun.y = mean, geom = "point", aes(shape = data)) + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange", data = filter(j, method == 'OTU-based dissimilarity' & group == 'Dissimilarity to next sample' & tbin < 39)) +
    geom_line(data = filter(tmp, method == 'OTU-based dissimilarity' & group == 'Dissimilarity to next sample')) +
    labs(x = 'Months after birth', y = 'OTU-based dissimilarity\nto next sample', tag = 'a  ') +
    scale_linetype_discrete(name = '') +
    scale_shape_manual(values = c(16,1), name = '') +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    scale_y_continuous(labels=scaleFUN) +
    expand_limits(y = 1) +
    theme_classic() +
  theme(legend.position = 'none')

p2 <- j %>%
  filter(method == 'OTU-based dissimilarity' & group == 'Dissimilarity to final sample') %>%
  ggplot(aes(x = tbin, y = dist, lty = data)) +
    stat_summary(fun.y = mean, geom = "point", aes(shape = data)) + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange", data = filter(j, method == 'OTU-based dissimilarity' & group == 'Dissimilarity to final sample' & tbin < 39)) +
    geom_line(data = filter(tmp, method == 'OTU-based dissimilarity' & group == 'Dissimilarity to final sample')) +
    labs(x = 'Months after birth', y = 'OTU-based dissimilarity\nto final sample', tag = 'b  ') +
    scale_linetype_discrete(name = '') +
    scale_shape_manual(values = c(16,1), name = '') +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
  scale_y_continuous(labels=scaleFUN) +
  expand_limits(y = 1) +
  theme_classic() +
  theme(legend.position = 'none')
  
p3 <- j %>%
  filter(method == 'Trait-based dissimilarity' & group == 'Dissimilarity to next sample') %>%
  ggplot(aes(x = tbin, y = dist, lty = data)) +
  stat_summary(fun.y = mean, geom = "point", aes(shape = data)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", data = filter(j, method == 'Trait-based dissimilarity' & group == 'Dissimilarity to next sample' & tbin < 39)) +
  geom_text(aes(label = symb, y = dist + 0.5), data = filter(tmp, method == 'Trait-based dissimilarity' & group == 'Dissimilarity to next sample'), size = 6) +
  geom_line(data = filter(tmp, method == 'Trait-based dissimilarity' & group == 'Dissimilarity to next sample')) +
  labs(x = 'Months after birth', y = 'Trait-based dissimilarity\nto next sample', tag = 'c  ') +
  scale_linetype_discrete(name = '') +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
  scale_y_continuous(labels=scaleFUN) +
  expand_limits(y = 6) +
  theme_classic() +
  theme(legend.position = 'none')

p4 <- j %>%
  filter(method == 'Trait-based dissimilarity' & group == 'Dissimilarity to final sample') %>%
  ggplot(aes(x = tbin, y = dist, lty = data)) +
  stat_summary(fun.y = mean, geom = "point", aes(shape = data)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", data = filter(j, method == 'Trait-based dissimilarity' & group == 'Dissimilarity to final sample' & tbin < 39)) +
  geom_text(aes(label = symb, y = dist + 0.5), data = filter(tmp, method == 'Trait-based dissimilarity' & group == 'Dissimilarity to final sample'), size = 6) +
  geom_line(data = filter(tmp, method == 'Trait-based dissimilarity' & group == 'Dissimilarity to final sample')) +
  labs(x = 'Months after birth', y = 'Trait-based dissimilarity\nto final sample', tag = 'd  ') +
  scale_linetype_discrete(name = '') +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
  scale_y_continuous(labels=scaleFUN) +
  theme_classic() +
  theme(legend.position = 'bottom', legend.justification="center")

myleg <- get_legend(p4)

p4 <- p4 + theme(legend.position = 'none')

ps <- plot_grid(p1, p2, p3, p4, ncol = 2)
figIntrababyDiss <- plot_grid(ps, myleg, ncol = 1, rel_heights = c(1, .07))

##Figure of interbaby dissimilarity
figInterbabyDiss <- function(){}

#use meltdist() to calculate dissimilarities within subjects for each method
#filter out low density tbins
j1 <- otus %>%
  spread(otu, abun, fill = 0) %>%
  mutate(tbin = round(t)) %>%
  filter(tbin < 36) %>%
  group_by(tbin) %>%
  do(meltdist(x = .$sampleID, y = ., method = 'bray')) %>%
  mutate(method = 'OTU-based dissimilarity', data = 'obs')

#second method: we use CWMs for columns, and calculate scaled euclidean distance
j2 <- otus %>%
  left_join(traits, by = 'otu') %>%
  mutate(tbin = round(t)) %>%
  filter(tbin < 36) %>%
  filter(!is.na(trait)) %>%
  group_by(sampleID, subject, trait, tbin, t) %>%
  summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
  group_by(trait) %>%
  mutate(val = scale(val)) %>%
  spread(trait, val) %>%
  group_by(tbin) %>%
  do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
  mutate(method = 'Trait-based dissimilarity', data = 'obs')

#recalculate null model predictions? With 1000 reps it takes ~1hr on my machine.
if (FALSE) {
  
  #empty list to store reps
  j3 <- list()
  
  #number of reps
  reps <- 1000
  for (rep in 1:reps) {
  
    #sample() shuffles values in a vector
    tmp <- otus %>%
      left_join(traits, by = 'otu') %>%
      group_by(trait) %>%
      mutate(val = sample(val)) %>%
      filter(!is.na(val)) %>%
      mutate(tbin = round(t)) %>%
      group_by(sampleID, subject, trait, tbin, t) %>%
      summarise(val = weighted.mean(val, w = abun, na.rm = TRUE)) %>%
      group_by(trait) %>%
      mutate(val = scale(val)) %>%
      spread(trait, val) %>%
      group_by(tbin) %>%
      do(meltdist(x = .$sampleID, y = ., method = 'euclidean')) %>%
      mutate(method = 'Trait-based dissimilarity', data = 'ran')
    
    j3[[rep]] <- tmp
    print(paste(rep, 'of', reps, 'done'))
    flush.console()
    
    if (rep == reps) {
      
      #save null simulation data
      bind_rows(j3) %>%
        group_by(tbin, sample1, sample2, t1, t2, method, data) %>%
        summarise(dist = mean(dist)) %>%
        saveRDS(file = 'data\\null_interbaby_predictions.RDS')
      
}}}

j3 <- readRDS('data\\null_interbaby_predictions.RDS')

#determine difference, if any, between simulated and observed initial trait dissimilarity
#later this will be used to adjust simulated and observed data to start together
tb_bump <- bind_rows(j2, j3) %>%
  ungroup() %>%
  mutate(t1 = t1 %/% 6 * 6 + 3) %>%
  filter(t1 == min(t1)) %>%
  summarise(diff = mean(dist[data == 'obs']) - 
              mean(dist[data == 'ran'])) %>%
  pull(diff)

#apply tb_bump
j3$dist <- j3$dist + tb_bump

#combine
interdist <- bind_rows(j1, j2, j3) %>%
  mutate(subject1 = gsub("_.*", "", sample1),
         subject2 = gsub("_.*", "", sample2)) %>%
  filter(subject1 != subject2) %>%
  ungroup()

#calculate means in 6 mo bins
j <- interdist %>%
  group_by(method, data, tbin) %>%
  summarise(dist = mean(dist), n = length(data)) %>%
  mutate(tbin = (tbin %/% 6) * 6 + 3) %>%
  ungroup() %>%
  mutate(data = factor(c('Observed','Null model')[match(data, c('obs','ran'))], c('Observed','Null model')))

#perform statistical tests, add symbols
statz <- j %>%
  filter(method == "Trait-based dissimilarity") %>%
  filter(tbin < 39) %>%
  group_by(method, tbin) %>%
  do(mod = t.test(dist ~ data, .)) %>%
  mutate(pval = mod$p.value,
         sig = pval < 0.05,
         symb = ifelse(pval < 0.001, '***',
                       ifelse(pval < 0.01, '**',
                              ifelse(pval < 0.05, '*', ''))))

j <- j %>% left_join(statz, by = c("method", "tbin")) %>% filter(tbin < 39)

#remove non-significant symbols
tmp <- j %>%
  group_by(method, data, tbin, sig, symb) %>%
  summarise(dist = mean(dist)) %>%
  group_by(method, tbin) %>%
  mutate(symb = ifelse(sig & dist == max(dist), symb, NA),
         sig = ifelse(is.na(symb), NA, sig))

#plot
p1 <- j %>%
  filter(method == 'OTU-based dissimilarity') %>%
  ggplot(aes(x = tbin, y = dist, lty = data, shape = data)) +
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange") +
    geom_line(data = filter(tmp, method == 'OTU-based dissimilarity')) +
    geom_text(aes(label = symb, y = dist + 0.4), data = filter(tmp, sig & method == 'OTU-based dissimilarity'), size = 6) +
    scale_shape_manual(values = c(16,1), name = '') +
    scale_linetype_discrete(name = '') +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    expand_limits(y = c(0)) +
    labs(x = 'Months after birth', y = 'OTU-based dissimilarity\namong infants', tag = 'a') +
    theme_classic() +
    theme(legend.position = 'none')

p2 <- j %>%
  filter(method == 'Trait-based dissimilarity') %>%
  ggplot(aes(x = tbin, y = dist, lty = data, shape = data)) +
  stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange") +
  geom_line(data = filter(tmp, method == 'Trait-based dissimilarity')) +
  geom_text(aes(label = symb, y = dist + 0.4), data = filter(tmp, sig & method == 'Trait-based dissimilarity'), size = 6) +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_linetype_discrete(name = '') +
  scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
  expand_limits(y = c(0)) +
  labs(x = 'Months after birth', y = 'Trait-based dissimilarity\namong infants', tag = 'b') +
  theme_classic() +
  theme(legend.position = 'bottom')

myleg <- get_legend(p2)

p2 <- p2 + theme(legend.position = 'none')

ps <- plot_grid(p1, p2, ncol = 2)

figInterbabyDiss <- plot_grid(ps, myleg, ncol = 1, rel_heights = c(1, .1))


##Figure showing interbaby distance by treatment
figInterbabyDissTreat <- function(){}

#create lists of subjects for each treatment group
subs_cs <- meta %>% filter(treatment_group == 'C-section') %>% pull(subject)
subs_ab <- meta %>% filter(treatment_group == 'Antibiotics') %>% pull(subject)
subs_control <- meta %>% filter(treatment_group == 'Control') %>% pull(subject)

#add metadata to pairwise subjects
#restrict to within-group comparisons
j <- interdist %>%
  mutate(
    delivery1 = meta$delivery[match(subject1, meta$subject)],
    delivery2 = meta$delivery[match(subject2, meta$subject)],
    treat1 = meta$treatment_group[match(subject1, meta$subject)],
    treat2 = meta$treatment_group[match(subject2, meta$subject)],
    delivery = ifelse(delivery1 == delivery2, delivery1, NA),
    treat = ifelse(treat1 == treat2, treat1, NA),
    tbin = floor(tbin / 6) * 6 + 3) %>%
  mutate(
    group = ifelse(subject1 %in% subs_cs & subject2 %in% subs_cs, 'C-section',
     ifelse(subject1 %in% subs_ab & subject2 %in% subs_ab, 'Antibiotics',
      ifelse(subject1 %in% subs_control & subject2 %in% subs_control, 'Control', NA)))) %>%
  filter(!is.na(group)) %>%
  ungroup() %>%
  select(method, group, tbin, dist)

#perform t-tests between each treatment and the control
statz <- bind_rows(j %>% filter(group %in% c('Antibiotics', 'Control')) %>% mutate(group2 = 'Antibiotics'),
                   j %>% filter(group %in% c('C-section', 'Control')) %>% mutate(group2 = 'C-section')) %>%
  group_by(method, group2, tbin) %>%
  do(mod = summary(lm(dist ~ group, data = .))) %>%
  mutate(pval = mod$coefficients[[8]])

#add stats to data, calculate locations for significance markers
j <- j %>%
  filter(tbin < 38) %>%
  left_join(statz, by = c("method", "tbin")) %>%
  mutate(
    sig = group != 'Control' & pval < 0.05,
    group = factor(group, c("Control", "Antibiotics","C-section")),
    tbin0 = tbin) %>%
  group_by(method, group, tbin) %>%
  mutate(mean = mean(dist),
         yval = ifelse(method == 'OTU-based dissimilarity',
                  ifelse(group == 'Antibiotics', 0.85, 0.9),
                    ifelse(group == 'Antibiotics', 6.5, 7)),
         symb = ifelse(group == 'Control', NA,
                  ifelse(pval < 0.001, '***',
                    ifelse(pval < 0.01, '**',
                      ifelse(pval < 0.05, '*', '')))))

# plot (almost every single comparison is significant due to the high number of pairwise comparisons, 
# so... i omit significance symbols)
p1 <- j %>%
  filter(method == 'OTU-based dissimilarity') %>%
  ggplot(aes(x = tbin, y = dist, color = group)) +
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange") +
    geom_line(aes(group = group, y = mean)) +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    scale_color_manual(values = c('black',"#00A480","#FF6200"), name = '') +
    expand_limits(y = c(0)) +
    labs(x = 'Months after birth', y = 'OTU-based dissimilarity among infants', tag = 'a  ') +
    theme_classic() +
    theme(legend.position = 'none')

p2 <- j %>%
  filter(method == 'Trait-based dissimilarity') %>%
  ggplot(aes(x = tbin, y = dist, color = group)) +
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_cl_boot, geom = "linerange") +
    geom_line(aes(group = group, y = mean)) +
    scale_x_continuous(breaks = c(3,9,15,21,27,33)) +
    scale_color_manual(values = c('black',"#00A480","#FF6200"), name = '') +
    expand_limits(y = c(0)) +
    labs(x = 'Months after birth', y = 'Trait-based dissimilarity among infants', tag = 'b  ') +
    theme_classic() +
    theme(legend.position = 'bottom')


myleg <- get_legend(p2)

p2 <- p2 + theme(legend.position = 'none')

ps <- plot_grid(p1, p2, ncol = 2)

figInterbabyDissTreat <- plot_grid(ps, myleg, ncol = 1, rel_heights = c(1, .06))


#### Supplementary data ####
write.csv(traits_sparse, 'images\\SupplementaryData_1_trait_value_observations.csv')
write.csv(traits_wide, 'images\\SupplementaryData_2_trait_value_predictions.csv')

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
tmp <- otus %>% 
  filter(t > 1 & t < 36) %>% 
  left_join(meta, by = 'subject') %>% 
  mutate(treat = ifelse(delivery == 'caesaren', 'caesaren', treatment_group)) %>% 
  group_by(treat, t %/% 6) %>% summarise(n = length(unique(sampleID))) %>% 
  group_by(treat) %>% 
  summarise(min = min(n), max = max(n))

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
sampstats <- otus %>% group_by(subject) %>% summarise(n = length(unique(t))) %>% ungroup() %>% summarise(min = min(n), max = max(n), mean = mean(n), sd = sd(n), median = median(n)) %>% mutate_all(funs(round), 2)

#number of data in sparse traits
sparse_tdats <- traits_sparse %>% select(-Genus, -Species) %>% gather(trait, val) %>% filter(!is.na(val)) %>% summarise(n = length(val)) %>% pull(n)

#correlation of gene number and genome size
gene_genome_corr <- round(cor(traits_sparse$Gene_number, traits_sparse$Genome_Mb, use = 'na.or.complete'),3)

#coverage of trait data by trait
tc <- trait_coverage %>% filter(Source != 'Missing data') %>% group_by(trait) %>% summarise(coverage = round(100 * sum(abun), 1)) %>% arrange(coverage) %>% filter(coverage > 50) %>% filter(coverage %in% c(min(coverage), max(coverage)))

#plot breaking down the relationship between c-section and antibiotics across subjects (not included in manuscript currently)
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

