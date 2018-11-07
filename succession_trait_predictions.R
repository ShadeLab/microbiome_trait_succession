#Trait-based Succession -- inferring unknown traits using phylogeny
#John Guittar

#load packages and set working directory
wd <- 'C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\data'
setwd(wd)
source('C:\\Users\\John\\Documents\\msu\\microbiome_trait_succession\\succession_custom_functions.R')
loadpax(c('tidyverse','ape','castor'))

#load tree with LTP and mouse tips, calculated with usearch
tre <- read_tree(write.tree(read.tree('LTP_succession.tree')))

#load rarefied OTU abundances
otus <- readRDS(file = 'otus.RDS')

#load traits, remove genome size because it is essentially the same as gene number
traits <- as.data.frame(readRDS('traits_sparse.RDS')) %>% 
  mutate(binomial = paste(Genus, Species)) %>%
  select(-Genome_Mb)
trait_names <- names(traits)[!names(traits) %in% c('Genus','Species','binomial')]

#load taxonomy files (used to link otus and traits)
#files were created in succession_processing_trait_data.R
tax_LTP <- readRDS('LTP_tax.RDS')
tax_succ <- readRDS('succession_tax_SILVA.RDS') %>%
  select(otu, Genus, Species)
tax <- bind_rows(tax_LTP, tax_succ) %>%
  mutate(binomial = paste(Genus, Species))

#join traits and OTUs using curated latin binomials
tax_traits <- tax %>% left_join(select(traits, -Genus, -Species), by = c("binomial"))

##Estimating unknown trait values, and evaluating confidence in those estimations
#I made this an option because it can take a really long time (days), depending on the number of queries
#loop through each trait...
for(i in trait_names) {
  
  #identify trait values of each tip
  tip_states <- tax_traits[[i]][match(tre$tip.label, tax_traits$otu)]
  
  #estimate unknown trait values
  IC <- hsp_independent_contrasts(tre, tip_states)
  SBT <- hsp_subtree_averaging(tre, tip_states)
  SCP <- hsp_squared_change_parsimony(tre, tip_states)
  
  #save temporary data for this trait 
  dat_trait <- data.frame(
    otu = tre$tip.label, 
    trait = i,
    IC = IC$states[1:length(tip_states)],
    SBT = SBT$states[1:length(tip_states)], 
    SCP = SCP$states[1:length(tip_states)], 
    stringsAsFactors = FALSE)
  
  #identify nearest otu with directly observed trait data
  nt <- find_nearest_tips(tre, target_tips = tre$tip.label[!is.na(tip_states)])
  dat_trait$nearest_measured_otu <- tre$tip.label[nt$nearest_tip_per_tip]
  dat_trait$dist <- nt$nearest_distance_per_tip
  
  ##manually calculate confidence of estimates
  #create tree with only measured tips
  tre_tmp <- drop.tip(tre, tre$tip.label[is.na(tip_states)])
  
  #create list of all possible combinations of measured otus, and their distances
  obs_trait <- as.data.frame(t(combn(tre_tmp$tip.label, 2)), stringsAsFactors = FALSE) %>%
    rename(otu1 = V1, otu2 = V2)
  
  j <- get_pairwise_distances(tre_tmp, obs_trait$otu1, obs_trait$otu2, 
                              as_edge_counts=FALSE, check_input=FALSE)
  
  #attach trait data, calculate difference in trait values
  obs_trait <- data.frame(obs_trait, dist = j) %>%
    select(otu1, otu2, dist) %>%
    mutate(
      trait = i,
      dist_bin = dist %/% 0.005 * 0.005,
      delta = abs(tax_traits[[i]][match(otu1, tax_traits$otu)] -
                    tax_traits[[i]][match(otu2, tax_traits$otu)]))
  
  #remove pairwise distances larger than dist <= 0.2, because we already know those are 'random')
  # drop data beyond 1000 randomly selected pairwise differences for each 0.001 increment in dist
  obs_trait <- obs_trait %>%
    group_by(dist_bin) %>%
    filter(seq_along(dist) %in% sample(seq_along(dist), min(length(dist), 1e4))) %>%
    ungroup()
  
  #progress bar
  print(i); flush.console()
  
  if (i == trait_names[1]) {
    dat <- dat_trait 
    obs <- obs_trait
  } else {
    dat <- bind_rows(dat, dat_trait)
    obs <- bind_rows(obs, obs_trait)
  }
}

# before fitting models, I first focus on just the space between 0 and 0.2 phylogenetic distance
# I revert to a condensed version of Gram-positive data, in which I take the average delta for each 0.001 increment of dist (i.e., each 'dist_bin') and then remove an outlier. I do this because the logistic model (obviously the right model) wasn't converging.
mods <- obs %>% filter(dist_bin <= 0.2 & dist > 0)
mods <- mods %>%
  filter(trait != 'Gram_positive') %>%
  bind_rows(mods %>% 
              filter(trait == 'Gram_positive') %>%
              mutate(dist_bin = dist %/% 0.005 * 0.005) %>%
              group_by(trait, dist_bin) %>% 
              mutate(delta = mean(delta)) %>% 
              distinct(dist, delta, .keep_all = TRUE) %>%
              filter(delta < 0.6))

#fit various models to delta~dist for each trait
#try_default ensures that do() finishes, even if the model does not converge or generates errors
mods <- mods %>%
  group_by(trait) %>%
  do(
    null = lm(delta ~ 1, data = .),
    linear = lm(delta ~ dist, data = .),
    linear_asym = plyr::try_default(nls(delta~SSasympOrig(dist, Asym, lrc), data = .),
                                    default = NA, quiet = TRUE),
    log = lm(delta~log(dist), data = .),
    logistic = plyr::try_default(nls(delta ~ SSlogis(dist, Asym, xmid, scal), data = .), default = NA, quiet = TRUE))

#model names for plotting later
mod_names <- c(null = 'Null', linear = 'Linear regression', linear_asym = 'Asymptotic regression',
               log = 'Logarithmic regression', logistic = 'Logistic regression')

#remove failed models 
mods <- mods %>%
  gather(type, mod, -trait) %>%
  rowwise() %>%
  filter(class(mod) != 'logical')

#calculate AIC and coefficients to see if the model was so bad that it predicted a negative slope
#if coefficient is negative, only keep null model
#normalize AIC values to lowest value for each set of models
mods <- mods %>%
  mutate(
    AIC = AIC(mod),
    est = ifelse(class(mod) == 'nls', summary(mod)$coef[[3]], summary(mod)$coef[[2]])) %>%
  ungroup() %>%
  mutate(type = mod_names[match(type, names(mod_names))]) %>%
  group_by(trait) %>%
  filter(if (min(est) < 0) type == 'Null' else est > 0) %>%
  mutate(AIC = AIC - min(AIC)) %>%
  arrange(trait, AIC)

#create dataframe/table with AIC scores for all models
mod_aic <- mods %>%
  select(Trait = trait, type, AIC) %>%
  spread(type, AIC) %>%
  select(Trait, Null, `Linear regression`, `Asymptotic regression`, `Logarithmic regression`, `Logistic regression`)

#create list of best mods
best_mods <- mods %>%
  filter(AIC == 0) %>% 
  mutate(best = paste(trait, type)) %>% 
  pull(sort(best))

#create a dataframe of null expectations of delta for each trait
#To remove phylogenetic effects, null delta is calculated as the mean pairwise trait-based distance for all randoly selected pairs of otus with greater than 10% difference in their 16S V4 region -- except for Salt optimum and ph optimum, which have no observable phylogenetic effect
null_thresh <- obs %>%
  filter(dist_bin > 0.1) %>%
  group_by(trait, dist_bin) %>%
  mutate(bin_mean = mean(delta)) %>%
  group_by(trait) %>%
  summarise(
    min_bin_mean = min(bin_mean), 
    null = mean(delta),
    threshold = null - (null - min_bin_mean) * 0.1)

#predict fits for models
#select only the best fitting model
#drop any foolish predicted deltas less than 0
#append null expectations for each trait
#create a delta threshold at the point when the model reaches the delta null expectation
#identitify the maximum phylogenetic distance that falls below that threshold (I included a zero to avoid -Inf results)
preds <- mods %>%
  rowwise() %>%
  do(data.frame(stringsAsFactors = FALSE,
                trait = .$trait,
                type = .$type,
                dist = seq(0, 0.2, length.out = 1000),
                delta = predict(.$mod, data.frame(dist = seq(0, 0.2, length.out = 1000))))) %>%
  filter(paste(trait, type) %in% best_mods) %>%
  filter(delta >= 0) %>%
  ungroup() %>%
  group_by(trait) %>%
  mutate(
    dist_bin = dist,
    threshold = null_thresh$threshold[match(trait, null_thresh$trait)],
    max_dist = max(dist[delta < threshold], 0),
    max_dist = ifelse(max_dist == max(dist) | max_dist < 0.03, 0, max_dist)) %>%
  ungroup()

#create version of trait data for which values that cannot be safely inferred are removed
traits <- dat %>%
  left_join(distinct(preds[, c('trait','max_dist')]), by = 'trait') %>%
  filter(otu %in% otus$otu) %>%
  filter(dist < max_dist) %>%
  transmute(otu, trait, val = SCP)

#write files
saveRDS(preds, 'trait_delta_models.RDS')
saveRDS(dat, 'traits_all_predictions.RDS')
saveRDS(traits, 'traits.RDS')
saveRDS(obs, 'trait_deltas.RDS')
